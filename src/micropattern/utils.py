from typing import Any

import numpy as np
import scanpy as sc
from anndata import AnnData
from loguru import logger

__all__ = ["normalize_and_select_hvg"]


def normalize_and_select_hvg(
    adata: AnnData,
    do_norm: bool = True,
    compute_hvg: bool = True,
    target_sum: float = 1e4,
    n_top_genes: int = 2000,
    batch_key: str | None = "sample_labels",
    subset: bool = True,
    verbose: bool = True,
):
    """
    normalize_and_select_hvg(adata, do_norm, compute_hvg, target_sum, n_top_genes, batch_key)

    Perform preprocessing on the input anndata.

    Parameters
    ----------
    adata: anndata
        Input anndata object to be processed.
    do_norm: bool
        Whether to perform library size normalization and log1p transformation.
    compute_hvg: bool
        Whether to compute highly variable genes (HVG). (seurat_v3)
    target_sum: int
        Target sum for library size normalization.
    n_top_genes: int
        Number of top highly variable genes to select.
    batch_key: str | None
        Key for batch information in adata.obs.
    """
    done_hvg = "hvg" in adata.uns
    done_norm = "log1p" in adata.uns
    counts_available = "counts" in adata.layers
    if do_norm:
        if done_norm:
            if verbose:
                logger.info("adata already normalized (log1p), skipping normalization")
        else:
            X = adata.X
            values = X.data if hasattr(X, "data") else X  # type: ignore[union-attr]
            if not np.all(np.equal(np.mod(np.asarray(values), 1), 0)):
                raise ValueError(
                    "adata.X contains non integer values, check if it is library normalized"
                )
            if not counts_available and not done_norm:
                if verbose:
                    logger.info("Copying X to counts layer before normalization")
                adata.layers["counts"] = adata.X.copy()  # type: ignore[union-attr]
                counts_available = True
            if verbose:
                logger.info(
                    f"Normalizing to target sum {target_sum} and applying log1p transformation"
                )
            sc.pp.normalize_total(adata, target_sum=target_sum)
            sc.pp.log1p(adata)

    if compute_hvg:
        if done_hvg:
            if verbose:
                logger.info("adata already has HVG selection, skipping HVG computation")
        else:
            if not counts_available:
                raise ValueError("counts layer is not available, cannot compute HVG")
            if batch_key not in adata.obs:
                if verbose:
                    logger.warning(
                        f"batch_key '{batch_key}' not found in adata.obs, computing HVG on all cells"
                    )
                batch_key = None

            if not adata.raw:
                if verbose:
                    logger.info("Saving X to raw before computing HVG")
                adata.raw = adata

            if verbose:
                logger.info(
                    f"Computing top {n_top_genes} highly variable genes (seurat_v3)"
                )
            sc.pp.highly_variable_genes(
                adata,
                n_top_genes=n_top_genes,
                subset=subset,
                flavor="seurat_v3",
                layer="counts",
                batch_key=batch_key,
            )
