#!/usr/bin/env python

from pathlib import Path
import scanpy as sc
import pandas as pd
import micropattern as mic

cwd = Path(".").resolve()

mic.io.fetch_minn(cwd, type="adata", stage="qc")
mic.io.fetch_heemskerk(cwd)

minn_data_path = cwd / "minn_data"
heemskerk_data_path = cwd / "heemskerk_data"

adata_minn = sc.read_h5ad(minn_data_path / mic.io.MINN_ADATA_LIST["qc"])
adata_heemskerk = sc.read_h5ad(heemskerk_data_path / mic.io.HEEMSKERK_ADATA["adata"])
meta_heemskerk = pd.read_csv(
    heemskerk_data_path / mic.io.HEEMSKERK_ADATA["meta"], index_col=0
)
adata_heemskerk.obs = meta_heemskerk.loc[adata_heemskerk.obs.index]

adata_heemskerk.obs["source"] = "heemskerk"
adata_minn.obs["source"] = "minn"
adata = sc.concat([adata_heemskerk, adata_minn])
adata.obs = pd.concat([adata_heemskerk.obs, adata_minn.obs])

cc_genes = mic.io.read_cc_list(cwd)

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
        base_run_name=f"minn_heemskerk_micropattern_0-10d_{n_hvg}hvg",
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
