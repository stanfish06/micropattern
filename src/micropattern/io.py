import s3fs
from pathlib import Path
from fsspec.callbacks import TqdmCallback

__all__ = ["fetch_minn", "MINN_FILE_LIST"]

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
