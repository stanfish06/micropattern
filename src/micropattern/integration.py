from anndata import AnnData
import scvi
import scanpy as sc
from datetime import datetime
import tarfile
from typing import Literal
import tempfile
from pathlib import Path
import s3fs
from fsspec.callbacks import TqdmCallback
import itertools
from lightning.pytorch.callbacks import Callback

__all__ = ["integration", "reference_query_mapping", "integration_multi"]


class _EpochLogger(Callback):
    def on_train_epoch_end(self, trainer, pl_module):
        m = trainer.callback_metrics
        train = m.get("train_loss_epoch", float("nan"))
        val = m.get("validation_loss", float("nan"))
        print(
            f"  Epoch {trainer.current_epoch}: "
            f"train_loss={train:.4f}, val_loss={val:.4f}"
        )


def _upload_model_to_s3(
    model_dir: str,
    run_name: str,
    date_str: str,
    s3_path: str = "stan-sequencing-data/processed-active/models",
):
    temp_dir = tempfile.mkdtemp()
    tarball_name = f"{date_str}_SCVI_{run_name}.tar.gz"
    tarball_path = Path(temp_dir) / tarball_name

    print(f"Creating tarball {tarball_name}")
    with tarfile.open(tarball_path, "w:gz") as tar:
        tar.add(model_dir, arcname=Path(model_dir).name)

    s3 = s3fs.S3FileSystem()
    remote_path = f"{s3_path}/{tarball_name}"
    print(f"Uploading to s3://{remote_path}")
    s3.put(
        lpath=str(tarball_path),
        rpath=remote_path,
        callback=TqdmCallback(),
    )
    print(f"Upload complete: s3://{remote_path}")


def integration(
    adata: AnnData,
    run_name: str = "integration",
    batch_key: str = "sample_labels",
    layer: str | None = "counts",
    categorical_covariates_keys: list | None = None,
    continuous_covariates_keys: list | None = None,
    lr: float = 1e-5,
    n_latent: int = 30,
    n_hidden: int = 128,
    n_layers: int = 2,
    dropout: float = 0.1,
    likelihood: Literal["nb", "zinb"] = "nb",
    compute_umap: bool = False,
    n_neighbors: int = 15,
    umap_min_dist: float = 0.01,
):
    scvi.model.SCVI.setup_anndata(
        adata,
        layer=layer,
        batch_key=batch_key,
        categorical_covariate_keys=categorical_covariates_keys,
        continuous_covariate_keys=continuous_covariates_keys,
    )
    # some parameters i found that work well
    model = scvi.model.SCVI(
        adata,
        n_latent=n_latent,
        n_hidden=n_hidden,
        n_layers=n_layers,
        gene_likelihood=likelihood,
        deeply_inject_covariates=False,
        use_batch_norm="both",
        dropout_rate=dropout,
    )
    model.train(
        max_epochs=400,
        early_stopping=True,
        early_stopping_patience=50,
        check_val_every_n_epoch=1,
        plan_kwargs={"lr": lr},
        enable_progress_bar=False,
        callbacks=[_EpochLogger()],
    )

    date_str = datetime.now().strftime("%Y%m%d")
    adata.obsm["X_scvi"] = model.get_latent_representation(give_mean=True)
    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    model.adata.obsm["X_scvi_qzm"] = qzm
    model.adata.obsm["X_scvi_qzv"] = qzv
    adata.obsm["X_scvi_qzm"] = qzm
    adata.obsm["X_scvi_qzv"] = qzv
    model.minify_adata(use_latent_qzm_key="X_scvi_qzm", use_latent_qzv_key="X_scvi_qzv")

    if compute_umap:
        sc.pp.neighbors(
            adata, n_neighbors=n_neighbors, use_rep="X_scvi", random_state=0
        )
        sc.tl.umap(adata, min_dist=umap_min_dist)
        model.adata.obs["umap_x"] = adata.obsm["X_umap"][:, 0]
        model.adata.obs["umap_y"] = adata.obsm["X_umap"][:, 1]

    model_dir = f"{date_str}_SCVI_{run_name}"
    model.save(model_dir, save_anndata=True, overwrite=True)

    # save to s3
    _upload_model_to_s3(model_dir, run_name, date_str)


def reference_query_mapping():
    pass


def integration_multi(
    adata: AnnData,
    base_run_name: str = "integration",
    batch_key: str = "sample_labels",
    layer: str | None = "counts",
    categorical_covariates_keys: list | None = None,
    continuous_covariates_keys: list | None = None,
    lr_values: list[float] | None = None,
    n_latent_values: list[int] | None = None,
    n_hidden_values: list[int] | None = None,
    n_layers_values: list[int] | None = None,
    dropout_values: list[float] | None = None,
    likelihood_values: list[Literal["nb", "zinb"]] | None = None,
    compute_umap: bool = True,
    n_neighbors: int = 15,
    umap_min_dist: float = 0.01,
):
    if lr_values is None:
        lr_values = [1e-5]
    if n_latent_values is None:
        n_latent_values = [30]
    if n_hidden_values is None:
        n_hidden_values = [128]
    if n_layers_values is None:
        n_layers_values = [2]
    if dropout_values is None:
        dropout_values = [0.1]
    if likelihood_values is None:
        likelihood_values = ["nb"]

    param_combinations = list(
        itertools.product(
            lr_values,
            n_latent_values,
            n_hidden_values,
            n_layers_values,
            dropout_values,
            likelihood_values,
        )
    )

    print(f"Exploring {len(param_combinations)} parameter combinations")

    for i, (lr, n_latent, n_hidden, n_layers, dropout, likelihood) in enumerate(
        param_combinations, 1
    ):
        run_name = f"{base_run_name}_lr{lr}_lv{n_latent}_hd{n_hidden}_ly{n_layers}_dr{dropout}_{likelihood}"
        print(f"\n[{i}/{len(param_combinations)}] Running: {run_name}")
        print(
            f"  Parameters: lr={lr}, n_latent={n_latent}, n_hidden={n_hidden}, "
            f"n_layers={n_layers}, dropout={dropout}, likelihood={likelihood}"
        )

        try:
            integration(
                adata=adata,
                run_name=run_name,
                batch_key=batch_key,
                layer=layer,
                categorical_covariates_keys=categorical_covariates_keys,
                continuous_covariates_keys=continuous_covariates_keys,
                lr=lr,
                n_latent=n_latent,
                n_hidden=n_hidden,
                n_layers=n_layers,
                dropout=dropout,
                likelihood=likelihood,
                compute_umap=compute_umap,
                n_neighbors=n_neighbors,
                umap_min_dist=umap_min_dist,
            )
            print(f"  ✓ Completed: {run_name}")
        except Exception as e:
            print(f"  ✗ Failed: {run_name} - {str(e)}")

    print(f"\nCompleted all {len(param_combinations)} combinations")
