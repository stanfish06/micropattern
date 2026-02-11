import s3fs
from pathlib import Path
from fsspec.callbacks import TqdmCallback
from anndata import AnnData
import tempfile
from datetime import datetime

__all__ = ["fetch_minn", "upload_adata", "upload_minn_adata", "MINN_FILE_LIST"]

MINN_FILE_LIST = [
    "GSM5176718_gastruloid_0h_1.barcodes.tsv.gz",
    "GSM5176718_gastruloid_0h_1.genes.tsv.gz",
    "GSM5176718_gastruloid_0h_1.matrix.tsv.gz",
    "GSM5176719_gastruloid_0h_2.barcodes.tsv.gz",
    "GSM5176719_gastruloid_0h_2.genes.tsv.gz",
    "GSM5176719_gastruloid_0h_2.matrix.tsv.gz",
    "GSM5176720_gastruloid_12h_1.barcodes.tsv.gz",
    "GSM5176720_gastruloid_12h_1.genes.tsv.gz",
    "GSM5176720_gastruloid_12h_1.matrix.tsv.gz",
    "GSM5176721_gastruloid_12h_2.barcodes.tsv.gz",
    "GSM5176721_gastruloid_12h_2.genes.tsv.gz",
    "GSM5176721_gastruloid_12h_2.matrix.tsv.gz",
    "GSM5176722_gastruloid_24h_1.barcodes.tsv.gz",
    "GSM5176722_gastruloid_24h_1.genes.tsv.gz",
    "GSM5176722_gastruloid_24h_1.matrix.tsv.gz",
    "GSM5176723_gastruloid_24h_2.barcodes.tsv.gz",
    "GSM5176723_gastruloid_24h_2.genes.tsv.gz",
    "GSM5176723_gastruloid_24h_2.matrix.tsv.gz",
    "GSE169074_gastruloid.all.time.metadata.csv.gz",
]


def fetch_minn(path: Path):
    s3 = s3fs.S3FileSystem()
    for f in MINN_FILE_LIST:
        print(f"Downloading {f}")
        s3.get(
            rpath=f"stan-sequencing-data/processed-active/Minn-micropattern-0-24h/{f}",
            lpath=str(path / "minn_data" / f),
            callback=TqdmCallback(),
        )


def upload_adata(
    adata: AnnData,
    *,
    s3_path: str,
    local_path: Path | None = None,
    filename: str,
):
    if local_path is None:
        temp_dir = tempfile.mkdtemp()
        local_path = Path(temp_dir)

    local_file = local_path / filename
    print(f"Saving AnnData to {local_file}")
    adata.write_h5ad(local_file)

    s3 = s3fs.S3FileSystem()
    remote_path = f"{s3_path}/{filename}"
    print(f"Uploading to s3://{remote_path}")
    s3.put(
        lpath=str(local_file),
        rpath=remote_path,
        callback=TqdmCallback(),
    )

    print(f"Upload complete: s3://{remote_path}")


def upload_minn_adata(
    adata: AnnData,
    *,
    processing_stage: str = "qc",
    local_path: Path | None = None,
):
    date_str = datetime.now().strftime("%Y%m%d")
    filename = f"minn_gastruloid_0-24h_{processing_stage}_{date_str}.h5ad"
    upload_adata(
        adata,
        s3_path="stan-sequencing-data/processed-active/Minn-micropattern-0-24h",
        local_path=local_path,
        filename=filename,
    )
