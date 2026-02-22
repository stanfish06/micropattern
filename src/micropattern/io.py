import s3fs
import tarfile
from pathlib import Path
from fsspec.callbacks import TqdmCallback
from anndata import AnnData
import tempfile
from typing import Literal
import pandas as pd
from datetime import datetime

__all__ = [
    "fetch_minn",
    "fetch_minn_heemskerk_adata_hvg",
    "fetch_liang_amanda",
    "upload_adata",
    "upload_minn_adata",
    "MINN_FILE_LIST_0_24H",
    "MINN_FILE_LIST_44H",
    "MINN_ADATA_LIST",
    "HEEMSKERK_ADATA_D2_D10",
    "HEEMSKERK_ADATA_42H",
    "LIANG_AMANDA_MESO_DATA_LIST",
    "LIANG_AMANDA_MESO_ADATA",
    "SCVI_MODELS",
    "fetch_scvi_model",
]

MINN_FILE_LIST_0_24H = [
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

LIANG_AMANDA_MESO_DATA_LIST = [
    "novogene_demult/chicago/per_sample_outs",
    "novogene_demult/michigan/per_sample_outs",
]

LIANG_AMANDA_MESO_ADATA = {
    "qc": "adata_micropattern_meso_liang_amanda_qc_20260217.h5ad",
    "hvg": "adata_heemskerk_liang_amanda_micropattern_0-10d_hvg.h5ad",
}

SCVI_MODELS = {
    "heemskerk_liang_amanda_0-10d_1000hvg_lv30_nb": "20260217_SCVI_heemskerk_liang_amanda_micropattern_0-10d_1000hvg_20260217-123409_lr1e-05_lv30_hd128_ly2_dr0.1_nb.tar.gz",
    "heemskerk_liang_amanda_0-10d_1000hvg_lv30_zinb": "20260217_SCVI_heemskerk_liang_amanda_micropattern_0-10d_1000hvg_20260217-123409_lr1e-05_lv30_hd128_ly2_dr0.1_zinb.tar.gz",
    "heemskerk_liang_amanda_0-10d_1000hvg_lv60_nb": "20260217_SCVI_heemskerk_liang_amanda_micropattern_0-10d_1000hvg_20260217-123409_lr1e-05_lv60_hd128_ly2_dr0.1_nb.tar.gz",
    "heemskerk_liang_amanda_0-10d_1000hvg_lv60_zinb": "20260217_SCVI_heemskerk_liang_amanda_micropattern_0-10d_1000hvg_20260217-123409_lr1e-05_lv60_hd128_ly2_dr0.1_zinb.tar.gz",
    "heemskerk_liang_amanda_0-10d_750hvg_lv30_nb": "20260217_SCVI_heemskerk_liang_amanda_micropattern_0-10d_750hvg_20260217-123409_lr1e-05_lv30_hd128_ly2_dr0.1_nb.tar.gz",
    "heemskerk_liang_amanda_0-10d_750hvg_lv30_zinb": "20260217_SCVI_heemskerk_liang_amanda_micropattern_0-10d_750hvg_20260217-123409_lr1e-05_lv30_hd128_ly2_dr0.1_zinb.tar.gz",
    "heemskerk_liang_amanda_0-10d_750hvg_lv60_nb": "20260217_SCVI_heemskerk_liang_amanda_micropattern_0-10d_750hvg_20260217-123409_lr1e-05_lv60_hd128_ly2_dr0.1_nb.tar.gz",
    "heemskerk_liang_amanda_0-10d_750hvg_lv60_zinb": "20260217_SCVI_heemskerk_liang_amanda_micropattern_0-10d_750hvg_20260217-123409_lr1e-05_lv60_hd128_ly2_dr0.1_zinb.tar.gz",
}

MINN_FILE_LIST_44H = [
    "GSM4300502_gastruloid1.barcodes.tsv.gz",
    "GSM4300502_gastruloid1.genes.tsv.gz",
    "GSM4300502_gastruloid1.matrix.mtx.gz",
    "GSM4300503_gastruloid2.barcodes.tsv.gz",
    "GSM4300503_gastruloid2.genes.tsv.gz",
    "GSM4300503_gastruloid2.matrix.mtx.gz",
]

MINN_ADATA_LIST = {"qc": "minn_gastruloid_0-24h_qc_20260213.h5ad"}

HEEMSKERK_ADATA_D2_D10 = {
    "adata": "adata_timeseries_old_48-96h_new_D6-10_filtered_qc.h5ad",
    "meta": "MP_old_48-96h_new_D6-10_meta.csv",
}

HEEMSKERK_ADATA_42H = {
    "rep1": "adata_2020_force_9000.h5ad",
    "rep2": "adata_2021_BMP_contorl.h5ad",
}

MISC_DATA = {
    "minn_heemskerk_adata_hvg": "adata_minn_heemskerk_micropattern_0-10d_hvg.h5ad",
    "cc_genes": "Macosko_cell_cycle_genes.txt",
}


def fetch_minn_heemskerk_adata_hvg(path: Path):
    s3 = s3fs.S3FileSystem()
    s3.get(
        rpath=f"stan-sequencing-data/processed-active/Minn-Heemskerk-micropattern/{MISC_DATA['minn_heemskerk_adata_hvg']}",
        lpath=str(path / "misc_data" / MISC_DATA["minn_heemskerk_adata_hvg"]),
        callback=TqdmCallback(),
    )


def _batch_download(
    s3: s3fs.S3FileSystem,
    rpaths: list[str],
    lpaths: list[str],
) -> None:
    print(f"Downloading {len(rpaths)} files...")
    s3.get(rpath=rpaths, lpath=lpaths, callback=TqdmCallback())


def fetch_liang_amanda(
    path: Path,
    type: Literal["raw", "adata"] = "raw",
    stage: str = "qc",
) -> None:
    s3 = s3fs.S3FileSystem()
    base = "stan-sequencing-data/processed-active/micropattern-meso-liang-amanda"
    match type:
        case "raw":
            rpaths: list[str] = []
            lpaths: list[str] = []
            for loc in LIANG_AMANDA_MESO_DATA_LIST:
                for f in s3.find(f"{base}/{loc}"):
                    rpaths.append(f)
                    lpaths.append(
                        str(path / "liang_amanda_data" / f.removeprefix(f"{base}/"))
                    )
            _batch_download(s3, rpaths, lpaths)
        case "adata":
            if stage in LIANG_AMANDA_MESO_ADATA:
                f = LIANG_AMANDA_MESO_ADATA[stage]
                print(f"Downloading {f}")
                s3.get(
                    rpath=f"{base}/{f}",
                    lpath=str(path / "liang_amanda_data" / f),
                    callback=TqdmCallback(),
                )


def fetch_scvi_model(path: Path, name: str) -> Path:
    if name not in SCVI_MODELS:
        raise ValueError(
            f"Unknown model '{name}'. Available: {list(SCVI_MODELS.keys())}"
        )
    filename = SCVI_MODELS[name]
    s3 = s3fs.S3FileSystem()
    local_tar = path / filename
    print(f"Downloading {filename}")
    s3.get(
        rpath=f"stan-sequencing-data/processed-active/models/{filename}",
        lpath=str(local_tar),
        callback=TqdmCallback(),
    )
    extract_dir = path
    print(f"Extracting to {extract_dir}")
    with tarfile.open(local_tar) as tar:
        tar.extractall(path=extract_dir)
    local_tar.unlink()
    extract_dir = extract_dir / filename.removesuffix(".tar.gz")
    return extract_dir


def fetch_heemskerk(path: Path):
    s3 = s3fs.S3FileSystem()
    rpaths = [
        f"stan-sequencing-data/processed-active/Heemskerk-micropattern-2-10d/{dat}"
        for dat in HEEMSKERK_ADATA_D2_D10.values()
    ] + [
        f"stan-sequencing-data/processed-active/Heemskerk-micropattern-42h/{dat}"
        for dat in HEEMSKERK_ADATA_42H.values()
    ]
    lpaths = [
        str(path / "heemskerk_data" / dat) for dat in HEEMSKERK_ADATA_D2_D10.values()
    ] + [str(path / "heemskerk_data" / dat) for dat in HEEMSKERK_ADATA_42H.values()]
    _batch_download(s3, rpaths, lpaths)


def read_cc_list(path: Path) -> dict | None:
    s3 = s3fs.S3FileSystem()
    f = MISC_DATA["cc_genes"]
    print(f"Downloading {f}")
    s3.get(
        rpath=f"stan-sequencing-data/processed-active/Heemskerk-micropattern-2-10d/{f}",
        lpath=str(path / "heemskerk_data" / f),
        callback=TqdmCallback(),
    )
    cc = pd.read_csv(path / "heemskerk_data" / f, sep="\t")
    # there is a dummy 6th column, probably due to trailing tab
    cc = cc.iloc[:, :5].to_dict(orient="list")
    for cat, gl in cc.items():
        cc[cat] = [g for g in gl if not pd.isna(g)]
    return cc


def fetch_minn(
    path: Path,
    type: Literal["mtx", "adata"],
    stage: str = "qc",
):
    s3 = s3fs.S3FileSystem()
    match type:
        case "mtx":
            rpaths = [
                f"stan-sequencing-data/processed-active/Minn-micropattern-0-24h/{f}"
                for f in MINN_FILE_LIST_0_24H
            ] + [
                f"stan-sequencing-data/processed-active/Minn-micropattern-44h/{f}"
                for f in MINN_FILE_LIST_44H
            ]
            lpaths = [str(path / "minn_data" / f) for f in MINN_FILE_LIST_0_24H] + [
                str(path / "minn_data" / f) for f in MINN_FILE_LIST_44H
            ]
            _batch_download(s3, rpaths, lpaths)
        case "adata":
            if stage in MINN_ADATA_LIST:
                f = MINN_ADATA_LIST[stage]
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
