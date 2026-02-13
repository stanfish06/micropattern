from anndata import AnnData
import scvi
from datetime import datetime
import tarfile
import tempfile
from pathlib import Path
import s3fs
from fsspec.callbacks import TqdmCallback

__all__ = ["integration", "reference_query_mapping"]


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
    layer: str = "counts",
    continuous_covariates_keys: list | None = None,
):
    scvi.model.SCVI.setup_anndata(
        adata,
        layer=layer,
        batch_key=batch_key,
        continuous_covariate_keys=continuous_covariates_keys,
    )
    # some parameters i found that work well
    model = scvi.model.SCVI(
        adata,
        n_latent=30,
        n_hidden=128,
        n_layers=2,
        gene_likelihood="nb",
        deeply_inject_covariates=False,
        use_batch_norm="both",
    )
    model.train(
        max_epochs=400,
        early_stopping=True,
        early_stopping_patience=50,
        check_val_every_n_epoch=1,
        plan_kwargs={"lr": 1e-3},
    )

    date_str = datetime.now().strftime("%Y%m%d")
    adata.obsm["X_scvi"] = model.get_latent_representation(give_mean=True)
    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    model.adata.obsm["X_scvi_qzm"] = qzm
    model.adata.obsm["X_scvi_qzv"] = qzv
    adata.obsm["X_scvi_qzm"] = qzm
    adata.obsm["X_scvi_qzv"] = qzv
    model.minify_adata(use_latent_qzm_key="X_scvi_qzm", use_latent_qzv_key="X_scvi_qzv")

    model_dir = f"{date_str}_SCVI_{run_name}"
    model.save(model_dir, save_anndata=True, overwrite=True)

    # save to s3
    _upload_model_to_s3(model_dir, run_name, date_str)


def reference_query_mapping():
    pass
