#!/usr/bin/env python

"""Train MultiVI from a notebook-prepared MuData and upload each HVG run."""

from __future__ import annotations

import gc
import json
import os
import re
import shutil
from pathlib import Path

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import scvi
import s3fs
import torch

ROOT = Path(__file__).resolve().parents[1]
INPUT_LOCAL = ROOT / "data/GSE245998/processed/multivi_ready.h5mu"
REQUIRE_CUDA = os.environ.get("MULTIVI_REQUIRE_CUDA", "1") == "1"
HVG_VALUES = [
    int(value.strip())
    for value in os.environ.get("MULTIVI_HVG_VALUES", "1500,2000").split(",")
    if value.strip()
]
ATAC_IMPUTE_THRESHOLD = float(os.environ.get("MULTIVI_ATAC_THRESHOLD", "0.05"))
IMPUTE_OUTPUT_CHUNK_SIZE = int(
    os.environ.get("MULTIVI_IMPUTE_CHUNK_SIZE", "2000")
)
MAX_EPOCHS = int(os.environ.get("MULTIVI_MAX_EPOCHS", "500"))
TRAIN_BATCH_SIZE = int(os.environ.get("MULTIVI_BATCH_SIZE", "256"))
N_LATENT = int(os.environ.get("MULTIVI_N_LATENT", "30"))

if not HVG_VALUES or any(value <= 0 for value in HVG_VALUES):
    raise ValueError(
        "MULTIVI_HVG_VALUES must be a comma-separated list of positive integers."
    )
if len(HVG_VALUES) != len(set(HVG_VALUES)):
    raise ValueError("MULTIVI_HVG_VALUES contains duplicate tiers.")
if not 0 < ATAC_IMPUTE_THRESHOLD <= 1:
    raise ValueError("MULTIVI_ATAC_THRESHOLD must be in (0, 1].")
if N_LATENT <= 0:
    raise ValueError("MULTIVI_N_LATENT must be a positive integer.")

INPUT_S3 = os.environ.get("MULTIVI_INPUT_S3")
if not INPUT_S3:
    raise RuntimeError(
        "Set MULTIVI_INPUT_S3 to the notebook-prepared integration-ready H5MU."
    )
INPUT_S3 = INPUT_S3.removeprefix("s3://")

OUTPUT_S3 = os.environ.get(
    "MULTIVI_OUTPUT_S3",
    (
        "stan-sequencing-data/processed-active/"
        "gse245998-tbxt-multiome-integration/runs"
    ),
).removeprefix("s3://")


def run_name(n_hvg: int, include_minn: bool) -> str:
    name = f"multivi_heemskerk_48_96h_42h_{n_hvg}hvg"
    if include_minn:
        name += "_minn_0_24h"
    if N_LATENT != 30:
        name += f"_nlatent{N_LATENT}"
    return name


def stage_input(fs: s3fs.S3FileSystem) -> None:
    INPUT_LOCAL.parent.mkdir(parents=True, exist_ok=True)
    expected = fs.info(INPUT_S3)["size"]
    if INPUT_LOCAL.exists() and INPUT_LOCAL.stat().st_size == expected:
        print(f"Already staged: {INPUT_LOCAL}", flush=True)
        return
    print(f"Downloading s3://{INPUT_S3} -> {INPUT_LOCAL}", flush=True)
    fs.get(INPUT_S3, str(INPUT_LOCAL))
    if INPUT_LOCAL.stat().st_size != expected:
        raise IOError(f"Incomplete download: {INPUT_LOCAL}")


def stream_file(
    fs: s3fs.S3FileSystem,
    output_local: Path,
    local: Path,
    current_run_name: str,
) -> None:
    """Upload one artifact to S3 and delete it locally to keep disk flat."""
    relative = local.relative_to(output_local).as_posix()
    remote = f"{OUTPUT_S3}/{current_run_name}/{relative}"
    print(f"Uploading {relative} -> s3://{remote}", flush=True)
    fs.put_file(str(local), remote)
    local.unlink()


def stream_dir(
    fs: s3fs.S3FileSystem,
    output_local: Path,
    local_dir: Path,
    current_run_name: str,
) -> None:
    """Upload every file under a directory to S3, then remove the directory."""
    for local in sorted(local_dir.rglob("*")):
        if not local.is_file():
            continue
        relative = local.relative_to(output_local).as_posix()
        remote = f"{OUTPUT_S3}/{current_run_name}/{relative}"
        print(f"Uploading {relative} -> s3://{remote}", flush=True)
        fs.put_file(str(local), remote)
    shutil.rmtree(local_dir, ignore_errors=True)


def prepare_hvg_tier(
    prepared: md.MuData, n_hvg: int
) -> tuple[md.MuData, list[str], pd.DataFrame]:
    rna_full = prepared.mod["rna"]
    atac_full = prepared.mod["atac"]
    required_var = {
        "gene_symbol",
        "highly_variable_rank",
        "macosko_cell_cycle",
    }
    missing_var = required_var - set(rna_full.var.columns)
    if missing_var:
        raise KeyError(
            "Prepared RNA is missing feature annotations: "
            + ", ".join(sorted(missing_var))
        )
    if rna_full.n_vars < n_hvg:
        raise ValueError(
            f"Prepared RNA has only {rna_full.n_vars} genes; cannot select {n_hvg}."
        )

    selected_genes = (
        rna_full.var["highly_variable_rank"].nsmallest(n_hvg).index.tolist()
    )
    if len(selected_genes) != n_hvg:
        raise ValueError(f"Expected {n_hvg} ranked HVGs; got {len(selected_genes)}.")
    cc_selected_genes = [
        gene
        for gene in selected_genes
        if bool(rna_full.var.at[gene, "macosko_cell_cycle"])
    ]
    cc_covariate_keys = [f"cc_expr__{gene}" for gene in cc_selected_genes]

    # Match the established micropattern scVI behavior: cell-cycle genes are
    # first selected from the HVG tier, then their log1p(CP10K) expression is
    # inserted into obs as individual continuous covariates.
    totals = np.asarray(rna_full.X.sum(axis=1)).ravel()
    if np.any(totals <= 0):
        raise ValueError("Cells with zero RNA counts cannot be normalized.")
    cc_counts = sp.csr_matrix(
        rna_full[:, cc_selected_genes].X,
        dtype=np.float32,
    )
    cc_normalized = cc_counts.multiply((1e4 / totals)[:, None]).tocsr()
    cc_normalized.data = np.log1p(cc_normalized.data)

    obs = rna_full.obs.copy()
    obs[cc_covariate_keys] = cc_normalized.toarray()
    rna = rna_full[:, selected_genes].copy()
    rna.obs = obs.copy()
    rna.X = sp.csr_matrix(rna.X)
    rna.layers["counts"] = rna.X.copy()
    atac = atac_full.copy()
    atac.obs = obs.copy()

    mdata = md.MuData({"rna": rna, "atac": atac})
    mdata.update()
    mdata.obs = obs.copy()
    if not mdata.mod["rna"].obs_names.equals(mdata.mod["atac"].obs_names):
        raise ValueError("Prepared RNA and ATAC observation axes are not aligned.")

    feature_selection = rna_full.var[
        [
            "gene_symbol",
            "highly_variable",
            "highly_variable_rank",
            "top1500",
            "top2000",
            "macosko_cell_cycle",
            "macosko_phases",
        ]
    ].copy()
    feature_selection["selected_for_run"] = feature_selection.index.isin(
        selected_genes
    )
    feature_selection["cell_cycle_covariate"] = feature_selection.index.isin(
        cc_selected_genes
    )
    print(
        f"top {n_hvg:,}: {len(cc_covariate_keys):,} Macosko "
        "cell-cycle continuous covariates",
        flush=True,
    )
    return mdata, cc_covariate_keys, feature_selection


def train_tier(
    fs: s3fs.S3FileSystem,
    prepared: md.MuData,
    n_hvg: int,
    output_local: Path,
    current_run_name: str,
) -> None:
    mdata, cc_covariate_keys, feature_selection = prepare_hvg_tier(
        prepared, n_hvg
    )

    scvi.model.MULTIVI.setup_mudata(
        mdata,
        batch_key="modality_batch",
        categorical_covariate_keys=["sample_id"],
        continuous_covariate_keys=cc_covariate_keys,
        modalities={
            "rna_layer": "rna",
            "atac_layer": "atac",
            "batch_key": "rna",
            "categorical_covariate_keys": "rna",
            "continuous_covariate_keys": "rna",
        },
    )
    model = scvi.model.MULTIVI(
        mdata,
        n_genes=mdata.mod["rna"].n_vars,
        n_regions=mdata.mod["atac"].n_vars,
        n_hidden=128,
        n_latent=N_LATENT,
        n_layers_encoder=2,
        n_layers_decoder=2,
        gene_likelihood="nb",
    )
    model.view_anndata_setup()
    model.train(
        accelerator="gpu" if torch.cuda.is_available() else "cpu",
        devices=1,
        batch_size=TRAIN_BATCH_SIZE,
        max_epochs=MAX_EPOCHS,
        early_stopping=True,
        check_val_every_n_epoch=5,
    )

    # Each artifact is streamed to S3 and deleted locally as soon as it is
    # written, so peak disk stays near a single file (the EC2 root volume is
    # small and accumulating the full run overflows it).
    output_local.mkdir(parents=True, exist_ok=True)

    model.save(output_local / "model", overwrite=True)
    stream_dir(fs, output_local, output_local / "model", current_run_name)

    feature_path = output_local / "rna_feature_selection.csv.gz"
    feature_selection.to_csv(feature_path, compression="gzip")
    stream_file(fs, output_local, feature_path, current_run_name)

    preparation = dict(prepared.uns.get("preparation", {}))
    run_config = {
        "scvi_tools_version": scvi.__version__,
        "torch_version": torch.__version__,
        "prepared_input_s3": f"s3://{INPUT_S3}",
        "n_cells": mdata.n_obs,
        "n_genes": mdata.mod["rna"].n_vars,
        "n_regions": mdata.mod["atac"].n_vars,
        "n_latent": N_LATENT,
        "n_hvg": n_hvg,
        "hvg_selection_source": preparation.get(
            "hvg_selection_source", "heemskerk_48_96h"
        ),
        "hvg_flavor": preparation.get("hvg_flavor", "seurat_v3"),
        "hvg_pool_size": preparation.get("hvg_pool_size"),
        "cell_cycle_source": preparation.get("cell_cycle_source"),
        "cell_cycle_covariate_count": len(cc_covariate_keys),
        "cell_cycle_covariates": cc_covariate_keys,
        "atac_imputation_threshold": ATAC_IMPUTE_THRESHOLD,
        "max_epochs": MAX_EPOCHS,
        "batch_size": TRAIN_BATCH_SIZE,
        "include_minn": bool(preparation.get("include_minn", False)),
    }
    config_path = output_local / "run_config.json"
    config_path.write_text(json.dumps(run_config, indent=2))
    stream_file(fs, output_local, config_path, current_run_name)

    latent = model.get_latent_representation()
    mdata.obsm["X_multivi"] = latent
    latent_adata = ad.AnnData(X=latent, obs=mdata.obs.copy())
    latent_adata.obsm["X_multivi"] = latent
    sc.pp.neighbors(latent_adata, use_rep="X_multivi", n_neighbors=15)
    sc.tl.umap(latent_adata, min_dist=0.2)
    mdata.obsm["X_umap"] = latent_adata.obsm["X_umap"]
    meta_path = output_local / "cell_metadata_latent_umap.csv.gz"
    latent_adata.obs.assign(
        multivi_1=latent[:, 0],
        multivi_2=latent[:, 1],
        umap_1=latent_adata.obsm["X_umap"][:, 0],
        umap_2=latent_adata.obsm["X_umap"][:, 1],
    ).to_csv(meta_path, compression="gzip")
    stream_file(fs, output_local, meta_path, current_run_name)

    for metric, history in model.history.items():
        safe_metric = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(metric))
        history_path = output_local / f"history_{safe_metric}.csv"
        history.to_csv(history_path)
        stream_file(fs, output_local, history_path, current_run_name)

    integrated_path = output_local / "integrated_multivi.h5mu"
    mdata.write_h5mu(integrated_path)
    stream_file(fs, output_local, integrated_path, current_run_name)

    imputed_dir = output_local / "imputed_atac_heemskerk"
    imputed_dir.mkdir(exist_ok=True)
    heem_mask = mdata.obs["source"].astype(str).eq("heemskerk").to_numpy()
    heem_obs = mdata.obs.loc[heem_mask]
    peak_means = []
    for (dataset, timepoint), group in heem_obs.groupby(
        ["dataset", "timepoint_hours"],
        observed=True,
        sort=True,
    ):
        group_indices = mdata.obs_names.get_indexer(group.index)
        group_sum = np.zeros(mdata.mod["atac"].n_vars, dtype=np.float64)
        for part, start in enumerate(
            range(0, len(group_indices), IMPUTE_OUTPUT_CHUNK_SIZE)
        ):
            indices = group_indices[start : start + IMPUTE_OUTPUT_CHUNK_SIZE]
            accessibility = model.get_normalized_accessibility(
                indices=indices,
                threshold=ATAC_IMPUTE_THRESHOLD,
                batch_size=TRAIN_BATCH_SIZE,
                return_numpy=True,
            )
            accessibility = sp.csr_matrix(accessibility)
            group_sum += np.asarray(accessibility.sum(axis=0)).ravel()
            imputed = ad.AnnData(
                X=accessibility,
                obs=mdata.obs.iloc[indices].copy(),
                var=mdata.mod["atac"].var.copy(),
            )
            imputed.uns[
                "multivi_accessibility_threshold"
            ] = ATAC_IMPUTE_THRESHOLD
            imputed.uns[
                "values"
            ] = "MultiVI accessibility probabilities below threshold set to zero"
            name = f"{dataset}_{int(timepoint)}h_part{part:03d}.h5ad"
            chunk_path = imputed_dir / name
            imputed.write_h5ad(chunk_path, compression="gzip")
            stream_file(fs, output_local, chunk_path, current_run_name)
        peak_means.append(
            pd.Series(
                group_sum / len(group_indices),
                index=mdata.mod["atac"].var_names,
                name=f"{dataset}_{int(timepoint)}h",
            )
        )

    peak_path = (
        output_local
        / "heemskerk_mean_imputed_accessibility_by_timepoint.csv.gz"
    )
    pd.concat(peak_means, axis=1).to_csv(peak_path, compression="gzip")
    stream_file(fs, output_local, peak_path, current_run_name)

    # Everything is on S3 now; drop the (now-empty) local run tree.
    shutil.rmtree(output_local, ignore_errors=True)


def main() -> None:
    if REQUIRE_CUDA and not torch.cuda.is_available():
        raise RuntimeError(
            "CUDA is unavailable in the EC2 container. MultiVI training was not started."
        )

    scvi.settings.seed = 0
    torch.set_float32_matmul_precision("high")
    fs = s3fs.S3FileSystem()
    stage_input(fs)
    prepared = md.read_h5mu(INPUT_LOCAL)
    if set(prepared.mod) != {"rna", "atac"}:
        raise ValueError(
            f"Expected prepared RNA and ATAC modalities; found {list(prepared.mod)}"
        )
    # The MuData is fully in memory now; free the staged input (several GB) so
    # the small EC2 root volume has room for the streamed per-tier outputs.
    INPUT_LOCAL.unlink(missing_ok=True)
    preparation = dict(prepared.uns.get("preparation", {}))
    include_minn = bool(preparation.get("include_minn", False))

    print(f"HVG sweep: {HVG_VALUES}", flush=True)
    for n_hvg in HVG_VALUES:
        current_run_name = run_name(n_hvg, include_minn)
        output_local = (
            ROOT / "data/GSE245998/processed" / current_run_name
        )
        print(f"Starting {current_run_name}", flush=True)
        train_tier(fs, prepared, n_hvg, output_local, current_run_name)
        print(f"Completed {current_run_name}", flush=True)
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()


if __name__ == "__main__":
    main()
