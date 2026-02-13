#!/usr/bin/env python

import os
from datetime import datetime, timezone

import pandas as pd
import s3fs
import scanpy as sc

import micropattern as mic

run_tag = os.environ.get(
    "RUN_TAG", datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
)

S3_MINN = "stan-sequencing-data/processed-active/Minn-micropattern-0-24h"
S3_HEEMSKERK_D2_D10 = "stan-sequencing-data/processed-active/Heemskerk-micropattern-2-10d"
S3_HEEMSKERK_42H = "stan-sequencing-data/processed-active/Heemskerk-micropattern-42h"

s3 = s3fs.S3FileSystem()

print("Loading Minn adata from S3...")
with s3.open(f"{S3_MINN}/{mic.io.MINN_ADATA_LIST['qc']}", "rb") as f:
    adata_minn = sc.read_h5ad(f)

print("Loading Heemskerk D2-D10 adata from S3...")
with s3.open(f"{S3_HEEMSKERK_D2_D10}/{mic.io.HEEMSKERK_ADATA_D2_D10['adata']}", "rb") as f:
    adata_heemskerk_d2_d10 = sc.read_h5ad(f)

print("Loading Heemskerk D2-D10 metadata from S3...")
with s3.open(f"{S3_HEEMSKERK_D2_D10}/{mic.io.HEEMSKERK_ADATA_D2_D10['meta']}", "rb") as f:
    meta_heemskerk_d2_d10 = pd.read_csv(f, index_col=0)
adata_heemskerk_d2_d10.obs = meta_heemskerk_d2_d10.loc[adata_heemskerk_d2_d10.obs.index]

print("Loading Heemskerk 42h reps from S3...")
adata_heemskerk_42h_list = []
for s, f in mic.io.HEEMSKERK_ADATA_42H.items():
    print(f"  Loading {f} ({s})")
    with s3.open(f"{S3_HEEMSKERK_42H}/{f}", "rb") as fh:
        adata_tmp = sc.read_h5ad(fh)
    adata_tmp.obs["sample_labels"] = adata_tmp.obs["sample_labels"].str.cat(
        ["_42h"] * adata_tmp.shape[0]
    )
    adata_heemskerk_42h_list.append(adata_tmp)
adata_heemskerk_42h = sc.concat(adata_heemskerk_42h_list)

print("Loading cell cycle genes from S3...")
with s3.open(f"{S3_HEEMSKERK_D2_D10}/{mic.io.MISC_DATA['cc_genes']}", "rb") as f:
    cc = pd.read_csv(f, sep="\t")
cc = cc.iloc[:, :5].to_dict(orient="list")
cc_genes = {cat: [g for g in gl if not pd.isna(g)] for cat, gl in cc.items()}

adata_heemskerk_d2_d10.obs["source"] = "heemskerk"
adata_heemskerk_42h.obs["source"] = "heemskerk"
adata_minn.obs["source"] = "minn"
adata = sc.concat([adata_heemskerk_d2_d10, adata_heemskerk_42h, adata_minn])
adata.obs = pd.concat([
    adata_heemskerk_d2_d10.obs,
    adata_heemskerk_42h.obs,
    adata_minn.obs,
])
print(f"Combined dataset: {adata.shape[0]} cells x {adata.shape[1]} genes")

for n_hvg in range(500, 1001, 100):
    print(f"\n{'=' * 60}")
    print(f"Running with {n_hvg} HVGs")
    print(f"{'=' * 60}")

    adata_copy = adata.copy()
    mic.utils.normalize_and_select_hvg(adata_copy, n_top_genes=n_hvg)

    cc_genes_all = [g for gl in cc_genes.values() for g in gl]
    cc_genes_all = [
        g
        for g in adata_copy.var_names
        if g.split()[0] in cc_genes_all and g not in adata_copy.obs.columns
    ]
    adata_copy.obs = pd.concat(
        [adata_copy.obs, adata_copy[:, cc_genes_all].to_df()], axis=1
    )

    mic.integration.integration_multi(
        adata_copy,
        base_run_name=f"minn_heemskerk_micropattern_0-10d_{n_hvg}hvg_{run_tag}",
        categorical_covariates_keys=None,
        continuous_covariates_keys=cc_genes_all,
        lr_values=[1e-5],
        n_latent_values=[30, 50, 70],
        n_hidden_values=[128, 256],
        n_layers_values=[2],
        dropout_values=[0.01, 0.1],
        likelihood_values=["nb", "zinb"],
        compute_umap=True,
        n_neighbors=15,
        umap_min_dist=0.01,
    )
