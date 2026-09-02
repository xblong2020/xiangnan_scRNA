from __future__ import annotations

import csv
import gzip
import re
import zipfile
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


ROOT = Path(__file__).resolve().parents[1]
PUBLIC = ROOT / "data" / "public"
EXTRACTED = ROOT / "data" / "processed" / "extracted"
OUT = ROOT / "data" / "processed" / "h5ad"
META = ROOT / "metadata" / "conversion"


def ensure_dirs() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    META.mkdir(parents=True, exist_ok=True)


def unique_index(values: list[str]) -> list[str]:
    seen: dict[str, int] = {}
    out: list[str] = []
    for value in values:
        value = str(value)
        n = seen.get(value, 0)
        seen[value] = n + 1
        out.append(value if n == 0 else f"{value}_{n + 1}")
    return out


def write_h5ad(adata: ad.AnnData, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.obs_names = unique_index(list(adata.obs_names))
    adata.var_names = unique_index(list(adata.var_names))
    for frame in (adata.obs, adata.var):
        for col in frame.columns:
            if pd.api.types.is_object_dtype(frame[col]) or isinstance(frame[col].dtype, pd.CategoricalDtype):
                frame[col] = frame[col].fillna("").astype(str)
    if path.exists():
        path.unlink()
    adata.write_h5ad(path, compression="gzip")
    print(f"WROTE {path} obs={adata.n_obs} var={adata.n_vars}", flush=True)


def valid_h5ad(path: Path) -> bool:
    if not path.exists():
        return False
    try:
        with h5py.File(path, "r") as handle:
            return "X" in handle and "obs" in handle and "var" in handle
    except OSError:
        return False


def read_tsv_gz(path: Path, header: bool = False) -> list[list[str]]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = list(reader)
    return rows[1:] if header else rows


def read_10x_mtx(matrix: Path, features: Path, barcodes: Path, dataset: str, sample: str) -> ad.AnnData:
    mat = mmread(str(matrix)).tocsr()
    rows = read_tsv_gz(features)
    barcodes_rows = read_tsv_gz(barcodes)
    gene_ids = [r[0] for r in rows]
    gene_names = [r[1] if len(r) > 1 else r[0] for r in rows]
    feature_type = [r[2] if len(r) > 2 else "Gene Expression" for r in rows]
    barcodes_list = [r[0] for r in barcodes_rows]

    # Convert genes x cells to cells x genes and remove empty raw droplets.
    cell_sum = np.asarray(mat.sum(axis=0)).ravel()
    keep = cell_sum > 0
    x = mat[:, keep].T.tocsr()
    obs_names = [f"{sample}_{bc}" for bc, ok in zip(barcodes_list, keep) if ok]
    obs = pd.DataFrame(
        {"dataset": dataset, "sample": sample, "n_counts_raw": cell_sum[keep].astype(np.int64)},
        index=obs_names,
    )
    var = pd.DataFrame({"gene_id": gene_ids, "feature_type": feature_type}, index=gene_names)
    return ad.AnnData(X=x, obs=obs, var=var)


def convert_gse174748() -> list[Path]:
    dataset = "GSE174748"
    indir = EXTRACTED / dataset
    outdir = OUT / dataset
    paths: list[Path] = []
    for matrix in sorted(indir.glob("*_matrix.mtx.gz")):
        prefix = matrix.name.replace("_matrix.mtx.gz", "")
        sample = re.sub(r"^GSM\d+_", "", prefix)
        features = indir / f"{prefix}_features.tsv.gz"
        barcodes = indir / f"{prefix}_barcodes.tsv.gz"
        out = outdir / f"{sample}.h5ad"
        if valid_h5ad(out):
            print(f"SKIP {out}", flush=True)
            paths.append(out)
            continue
        adata = read_10x_mtx(matrix, features, barcodes, dataset, sample)
        write_h5ad(adata, out)
        paths.append(out)
    return paths


def convert_gse151530() -> list[Path]:
    dataset = "GSE151530"
    indir = PUBLIC / "geo" / dataset
    out = OUT / dataset / f"{dataset}.h5ad"
    if valid_h5ad(out):
        print(f"SKIP {out}", flush=True)
        return [out]
    adata = read_10x_mtx(
        indir / "GSE151530_matrix.mtx.gz",
        indir / "GSE151530_genes.tsv.gz",
        indir / "GSE151530_barcodes.tsv.gz",
        dataset,
        dataset,
    )
    info = pd.read_csv(indir / "GSE151530_Info.txt.gz", sep="\t")
    first_col = info.columns[0]
    info.index = info[first_col].astype(str)
    adata.obs = adata.obs.join(info, how="left")
    write_h5ad(adata, out)
    return [out]


def read_10x_h5(path: Path, dataset: str, sample: str) -> ad.AnnData:
    with h5py.File(path, "r") as handle:
        group = handle["matrix"]
        shape = tuple(group["shape"][:])
        data = group["data"][:]
        indices = group["indices"][:]
        indptr = group["indptr"][:]
        barcodes = [x.decode("utf-8") for x in group["barcodes"][:]]
        features = group["features"]
        gene_ids = [x.decode("utf-8") for x in features["id"][:]]
        gene_names = [x.decode("utf-8") for x in features["name"][:]]
        feature_type = [x.decode("utf-8") for x in features["feature_type"][:]]
    mat = sparse.csc_matrix((data, indices, indptr), shape=shape)
    cell_sum = np.asarray(mat.sum(axis=0)).ravel()
    keep = cell_sum > 0
    x = mat[:, keep].T.tocsr()
    obs_names = [f"{sample}_{bc}" for bc, ok in zip(barcodes, keep) if ok]
    obs = pd.DataFrame(
        {"dataset": dataset, "sample": sample, "n_counts_raw": cell_sum[keep].astype(np.int64)},
        index=obs_names,
    )
    var = pd.DataFrame({"gene_id": gene_ids, "feature_type": feature_type}, index=gene_names)
    return ad.AnnData(X=x, obs=obs, var=var)


def convert_gse212046() -> list[Path]:
    dataset = "GSE212046"
    indir = EXTRACTED / dataset
    outdir = OUT / dataset
    paths: list[Path] = []
    for h5 in sorted(indir.glob("*.h5")):
        sample = h5.name.replace("_raw_feature_bc_matrix.h5", "")
        sample = re.sub(r"^GSM\d+_", "", sample)
        out = outdir / f"{sample}.h5ad"
        if valid_h5ad(out):
            print(f"SKIP {out}", flush=True)
            paths.append(out)
            continue
        adata = read_10x_h5(h5, dataset, sample)
        write_h5ad(adata, out)
        paths.append(out)
    return paths


def convert_csv_counts(path: Path, dataset: str, sample: str) -> ad.AnnData:
    frame = pd.read_csv(path, index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = [f"{sample}_{c}" for c in frame.columns.astype(str)]
    x = sparse.csr_matrix(frame.to_numpy(dtype=np.float32).T)
    obs = pd.DataFrame({"dataset": dataset, "sample": sample}, index=frame.columns)
    var = pd.DataFrame(index=frame.index)
    return ad.AnnData(X=x, obs=obs, var=var)


def convert_gse202379(max_files: int | None = None) -> list[Path]:
    dataset = "GSE202379"
    indir = EXTRACTED / dataset
    outdir = OUT / dataset
    paths: list[Path] = []
    files = sorted(indir.glob("*-raw_counts.csv.gz"))
    if max_files:
        files = files[:max_files]
    for path in files:
        sample = re.sub(r"^GSM\d+_", "", path.name).replace("-raw_counts.csv.gz", "")
        out = outdir / f"{sample}.h5ad"
        if valid_h5ad(out):
            print(f"SKIP {out}", flush=True)
            paths.append(out)
            continue
        adata = convert_csv_counts(path, dataset, sample)
        write_h5ad(adata, out)
        paths.append(out)
    return paths


def convert_gse149614() -> list[Path]:
    dataset = "GSE149614"
    indir = PUBLIC / "geo" / dataset
    out = OUT / dataset / f"{dataset}_raw_counts.h5ad"
    if valid_h5ad(out):
        print(f"SKIP {out}", flush=True)
        return [out]
    path = indir / "GSE149614_HCC.scRNAseq.S71915.count.txt.gz"
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    cells = header
    genes: list[str] = []
    blocks: list[sparse.csr_matrix] = []
    reader = pd.read_csv(path, sep="\t", header=None, skiprows=1, index_col=0, chunksize=1000)
    for chunk in reader:
        genes.extend(chunk.index.astype(str).tolist())
        arr = chunk.to_numpy(dtype=np.float32, copy=False)
        blocks.append(sparse.csr_matrix(arr))
        print(f"  {dataset}: read genes={len(genes)}", flush=True)
    gene_by_cell = sparse.vstack(blocks, format="csr")
    x = gene_by_cell.T.tocsr()
    obs = pd.DataFrame({"dataset": dataset}, index=cells)
    var = pd.DataFrame(index=genes)
    meta = pd.read_csv(indir / "GSE149614_HCC.metadata.updated.txt.gz", sep="\t")
    if meta.shape[1] and "Cell" in meta.columns:
        meta.index = meta["Cell"].astype(str)
        obs = obs.join(meta, how="left")
    elif meta.shape[1]:
        meta.index = meta.iloc[:, 0].astype(str)
        obs = obs.join(meta, how="left", rsuffix="_meta")
    adata = ad.AnnData(X=x, obs=obs, var=var)
    write_h5ad(adata, out)
    return [out]


def find_10x_dirs(indir: Path) -> list[Path]:
    dirs: list[Path] = []
    for matrix in indir.rglob("matrix.mtx.gz"):
        if any(part.lower().startswith("spatial") for part in matrix.parts):
            continue
        parent = matrix.parent
        if (parent / "barcodes.tsv.gz").exists() and (
            (parent / "features.tsv.gz").exists() or (parent / "genes.tsv.gz").exists()
        ):
            dirs.append(parent)
    return sorted(set(dirs))


def convert_gse185477() -> list[Path]:
    dataset = "GSE185477"
    indir = EXTRACTED / dataset
    outdir = OUT / dataset
    paths: list[Path] = []
    for folder in find_10x_dirs(indir):
        sample = folder.name
        if sample == "raw_feature_bc_matrix":
            sample = folder.parent.name
        out = outdir / f"{sample}.h5ad"
        if valid_h5ad(out):
            print(f"SKIP {out}", flush=True)
            paths.append(out)
            continue
        features = folder / "features.tsv.gz"
        if not features.exists():
            features = folder / "genes.tsv.gz"
        adata = read_10x_mtx(folder / "matrix.mtx.gz", features, folder / "barcodes.tsv.gz", dataset, sample)
        write_h5ad(adata, out)
        paths.append(out)
    return paths


def main() -> int:
    ensure_dirs()
    manifest: list[tuple[str, str]] = []
    for converter in [
        convert_gse151530,
        convert_gse174748,
        convert_gse212046,
        convert_gse185477,
        convert_gse202379,
        convert_gse149614,
    ]:
        for path in converter():
            manifest.append((path.parent.name, str(path)))
    manifest_path = META / "h5ad_conversion_manifest.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["dataset", "h5ad_path"])
        writer.writerows(manifest)
    print(f"WROTE {manifest_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
