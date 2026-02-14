#!/usr/bin/env python

import gc
import itertools
import os
from datetime import datetime, timezone

import pandas as pd
import s3fs
import scanpy as sc

import micropattern as mic

run_tag = os.environ.get(
    "RUN_TAG", datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
)

S3_MERGED_HVG = "stan-sequencing-data/processed-active/Minn-Heemskerk-micropattern"
S3_HEEMSKERK_D2_D10 = (
    "stan-sequencing-data/processed-active/Heemskerk-micropattern-2-10d"
)

s3 = s3fs.S3FileSystem()

print("=" * 80)
print("Loading pre-computed HVG adata from S3...")
print("=" * 80)
with s3.open(
    f"{S3_MERGED_HVG}/adata_minn_heemskerk_micropattern_0-10d_hvg.h5ad", "rb"
) as f:
    adata = sc.read_h5ad(f)
print(f"Loaded dataset: {adata.shape[0]} cells x {adata.shape[1]} genes")
print(f"Available HVG tiers: {[col for col in adata.var.columns if col.startswith('top')]}")

print("\nLoading cell cycle genes from S3...")
with s3.open(f"{S3_HEEMSKERK_D2_D10}/{mic.io.MISC_DATA['cc_genes']}", "rb") as f:
    cc = pd.read_csv(f, sep="\t")
cc = cc.iloc[:, :5].to_dict(orient="list")
cc_genes = {cat: [g for g in gl if not pd.isna(g)] for cat, gl in cc.items()}

hvg_values = [750, 1000]
lr_values = [1e-5]
n_latent_values = [30, 60]
n_hidden_values = [128]
n_layers_values = [2]
dropout_values = [0.1]
likelihood_values = ["nb", "zinb"]

total_param_combinations = len(
    list(
        itertools.product(
            lr_values,
            n_latent_values,
            n_hidden_values,
            n_layers_values,
            dropout_values,
            likelihood_values,
        )
    )
)
total_combinations = len(hvg_values) * total_param_combinations
print(f"\n{'=' * 80}")
print("PARAMETER SWEEP CONFIGURATION")
print("=" * 80)
print(f"HVG values: {hvg_values}")
print(f"Learning rates: {lr_values}")
print(f"Latent dimensions: {n_latent_values}")
print(f"Hidden dimensions: {n_hidden_values}")
print(f"Layers: {n_layers_values}")
print(f"Dropout: {dropout_values}")
print(f"Likelihoods: {likelihood_values}")
print(f"\nTotal HVG tiers: {len(hvg_values)}")
print(f"Param combinations per HVG tier: {total_param_combinations}")
print(f"TOTAL COMBINATIONS: {total_combinations}")
print(f"Estimated runtime: ~{total_combinations * 0.75:.1f}-{total_combinations * 1:.1f} hours")
print(f"{'=' * 80}\n")

global_combination_idx = 0

for hvg_idx, n_hvg in enumerate(hvg_values, 1):
    hvg_col = f"top{n_hvg}"

    print(f"\n{'=' * 80}")
    print(f"HVG TIER {hvg_idx}/{len(hvg_values)}: {n_hvg} genes")
    print(f"{'=' * 80}")

    if hvg_col not in adata.var.columns:
        print(
            f"WARNING: Column '{hvg_col}' not found in adata.var, skipping {n_hvg} HVGs"
        )
        print(f"Available columns: {list(adata.var.columns)}")
        global_combination_idx += total_param_combinations
        print(f"Skipped {total_param_combinations} combinations\n")
        continue

    print(f"Subsetting to {n_hvg} highly variable genes using '{hvg_col}' column...")
    adata_hvg = adata[:, adata.var[hvg_col]].copy()
    print(f"Subsetted dataset: {adata_hvg.shape[0]} cells x {adata_hvg.shape[1]} genes")

    cc_genes_all = [g for gl in cc_genes.values() for g in gl]
    cc_genes_all = [
        g
        for g in adata_hvg.var_names
        if g.split()[0] in cc_genes_all and g not in adata_hvg.obs.columns
    ]
    print(f"Found {len(cc_genes_all)} cell cycle genes in HVG set")

    adata_hvg.obs = pd.concat(
        [adata_hvg.obs, adata_hvg[:, cc_genes_all].to_df()], axis=1
    )

    print(f"\nRunning {total_param_combinations} parameter combinations for {n_hvg} HVGs")
    print(f"Progress: {global_combination_idx}/{total_combinations} total combinations completed so far\n")

    mic.integration.integration_multi(
        adata_hvg,
        base_run_name=f"minn_heemskerk_micropattern_0-10d_{n_hvg}hvg_{run_tag}",
        categorical_covariates_keys=None,
        continuous_covariates_keys=cc_genes_all,
        lr_values=lr_values,
        n_latent_values=n_latent_values,
        n_hidden_values=n_hidden_values,
        n_layers_values=n_layers_values,
        dropout_values=dropout_values,
        likelihood_values=likelihood_values,
        compute_umap=True,
        n_neighbors=15,
        umap_min_dist=0.01,
    )

    del adata_hvg
    gc.collect()

    global_combination_idx += total_param_combinations
    print(f"\nCompleted HVG tier {hvg_idx}/{len(hvg_values)}")
    print(f"Overall progress: {global_combination_idx}/{total_combinations} combinations\n")

print(f"\n{'=' * 80}")
print("ALL PARAMETER COMBINATIONS COMPLETED!")
print(f"Total runs: {global_combination_idx}/{total_combinations}")
print(f"Run tag: {run_tag}")
print(f"{'=' * 80}")
