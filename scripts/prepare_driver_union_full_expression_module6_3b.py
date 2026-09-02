from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 6.3b: prepare full-gene expression matrix for canonical SCENIC.")
    parser.add_argument(
        "--driver-cells",
        type=Path,
        default=ROOT / "metadata/driver/driver_module6_1_cells.tsv.gz",
    )
    parser.add_argument(
        "--cellrank-fate",
        type=Path,
        default=ROOT / "metadata/driver/driver_module6_2_cellrank_fate_probabilities.tsv.gz",
    )
    parser.add_argument("--processed-dir", type=Path, default=ROOT / "data/processed/driver")
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/driver")
    parser.add_argument("--output-name", default="driver_union_full_expression.module6_3b.h5ad")
    parser.add_argument("--min-cells", type=int, default=20)
    parser.add_argument("--min-mean", type=float, default=0.005)
    parser.add_argument("--max-genes", type=int, default=12000)
    parser.add_argument("--compression", default="gzip")
    return parser.parse_args()


def normalize_cell_id(value: object, study_sample: str) -> str:
    text = str(value)
    prefix = f"{study_sample}__"
    if text.startswith(prefix):
        return text
    return f"{study_sample}__{text}"


def x_to_counts(matrix) -> sparse.csr_matrix:
    if sparse.issparse(matrix):
        return matrix.tocsr()
    return sparse.csr_matrix(np.asarray(matrix))


def read_subset_from_source(source_h5ad: Path, cells: pd.DataFrame) -> ad.AnnData:
    backed = ad.read_h5ad(source_h5ad, backed="r")
    obs_names = pd.Index(backed.obs_names.astype(str))
    study_sample = str(cells["study_sample"].iloc[0])
    wanted = cells["cell_id"].astype(str).tolist()

    direct = pd.Index(wanted)
    if direct.isin(obs_names).sum() == len(wanted):
        source_ids = direct
    else:
        stripped = pd.Index([cell[len(f"{study_sample}__") :] if cell.startswith(f"{study_sample}__") else cell for cell in wanted])
        if stripped.isin(obs_names).sum() == len(wanted):
            source_ids = stripped
        else:
            normalized_lookup = pd.Series(obs_names.to_numpy(), index=[normalize_cell_id(x, study_sample) for x in obs_names])
            missing = [cell for cell in wanted if cell not in normalized_lookup.index]
            if missing:
                backed.file.close()
                raise ValueError(f"{len(missing)} cells from {study_sample} are missing in {source_h5ad}")
            source_ids = pd.Index(normalized_lookup.loc[wanted].to_numpy())

    sub = backed[source_ids, :].to_memory()
    backed.file.close()
    sub.obs_names = pd.Index(wanted)
    return sub


def choose_matrix(adata: ad.AnnData):
    if "counts" in adata.layers:
        return adata.layers["counts"]
    return adata.X


def filter_genes(
    matrix: sparse.csr_matrix,
    genes: pd.Index,
    tf_list: set[str],
    min_cells: int,
    min_mean: float,
    max_genes: int,
) -> pd.Index:
    detected = np.asarray((matrix > 0).sum(axis=0)).ravel()
    mean_expr = np.asarray(matrix.mean(axis=0)).ravel()
    keep = (detected >= min_cells) & (mean_expr >= min_mean)
    scores = detected * np.log1p(mean_expr)
    genes_arr = genes.astype(str).to_numpy()
    # Preserve expressed TFs for GRN inference while excluding catalog entries
    # that are completely absent from the selected driver-union matrix.
    tf_mask = np.array([gene in tf_list for gene in genes_arr]) & (detected > 0)
    keep = keep | tf_mask
    kept_idx = np.flatnonzero(keep)
    if kept_idx.size > max_genes:
        tf_idx = np.flatnonzero(tf_mask)
        non_tf = np.array([idx for idx in kept_idx if idx not in set(tf_idx)])
        slots = max(0, max_genes - len(tf_idx))
        top_non_tf = non_tf[np.argsort(scores[non_tf])[::-1]][:slots]
        kept_idx = np.unique(np.concatenate([tf_idx, top_non_tf]))
    return pd.Index(genes_arr[kept_idx])


def main() -> None:
    start = time.time()
    args = parse_args()
    args.processed_dir.mkdir(parents=True, exist_ok=True)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)

    cells = pd.read_csv(args.driver_cells, sep="\t")
    cells["cell_id"] = cells["cell_id"].astype(str)
    cells["source_h5ad"] = cells["source_h5ad"].astype(str)
    fate = pd.read_csv(args.cellrank_fate, sep="\t") if args.cellrank_fate.exists() else pd.DataFrame({"cell_id": []})
    if not fate.empty:
        fate["cell_id"] = fate["cell_id"].astype(str)

    parts = []
    source_rows = []
    for source, sub_cells in cells.groupby("source_h5ad", sort=True):
        source_path = Path(source)
        if not source_path.exists():
            raise FileNotFoundError(source_path)
        sub = read_subset_from_source(source_path, sub_cells)
        matrix = x_to_counts(choose_matrix(sub))
        part = ad.AnnData(X=matrix, obs=sub_cells.set_index("cell_id").loc[sub.obs_names].copy(), var=pd.DataFrame(index=sub.var_names.astype(str)))
        parts.append(part)
        source_rows.append(
            {
                "source_h5ad": str(source_path),
                "n_cells": int(part.n_obs),
                "n_genes_source": int(sub.n_vars),
                "matrix_source": "counts_layer" if "counts" in sub.layers else "X",
            }
        )

    combined = ad.concat(parts, join="outer", fill_value=0, merge="same")
    combined.X = x_to_counts(combined.X)
    combined.obs_names_make_unique()
    if not fate.empty:
        fate_aligned = fate.set_index("cell_id").reindex(combined.obs_names.astype(str))
        for col in fate_aligned.columns:
            combined.obs[col] = fate_aligned[col].to_numpy()
    for col in combined.obs.columns:
        series = combined.obs[col]
        if pd.api.types.is_bool_dtype(series):
            combined.obs[col] = series.fillna(False).astype(bool)
        elif pd.api.types.is_numeric_dtype(series):
            combined.obs[col] = pd.to_numeric(series, errors="coerce")
        else:
            combined.obs[col] = series.astype("object").where(series.notna(), "Unknown").astype(str)

    tf_path = args.metadata_dir / "scenic_resources/allTFs_hg38.txt"
    tf_list = set(tf_path.read_text(encoding="utf-8").splitlines()) if tf_path.exists() else set()
    keep_genes = filter_genes(
        combined.X,
        combined.var_names,
        tf_list=tf_list,
        min_cells=args.min_cells,
        min_mean=args.min_mean,
        max_genes=args.max_genes,
    )
    filtered = combined[:, keep_genes].copy()
    for col in filtered.obs.columns:
        series = filtered.obs[col]
        if pd.api.types.is_bool_dtype(series):
            filtered.obs[col] = series.fillna(False).astype(bool)
        elif pd.api.types.is_numeric_dtype(series):
            filtered.obs[col] = pd.to_numeric(series, errors="coerce")
        else:
            filtered.obs[col] = pd.Categorical(series.astype("object").where(series.notna(), "Unknown").astype(str))
    filtered.layers["counts"] = filtered.X.copy()
    counts = filtered.X.copy()
    totals = np.asarray(counts.sum(axis=1)).ravel()
    scale = np.divide(1e4, totals, out=np.zeros_like(totals, dtype=float), where=totals > 0)
    norm = sparse.diags(scale).dot(counts)
    norm.data = np.log1p(norm.data)
    filtered.X = norm.tocsr()
    filtered.uns["module6_3b_full_expression"] = {
        "module": "6.3b",
        "method": "driver union full-gene expression matrix assembled from source QC h5ad files",
        "n_cells": int(filtered.n_obs),
        "n_genes_before_filter": int(combined.n_vars),
        "n_genes_after_filter": int(filtered.n_vars),
        "min_cells": int(args.min_cells),
        "min_mean": float(args.min_mean),
        "max_genes": int(args.max_genes),
    }

    output_h5ad = args.processed_dir / args.output_name
    filtered.write_h5ad(output_h5ad, compression=args.compression)

    source_manifest = pd.DataFrame(source_rows)
    source_manifest_path = args.metadata_dir / "driver_module6_3b_full_expression_sources.tsv"
    source_manifest.to_csv(source_manifest_path, sep="\t", index=False)

    gene_stats = pd.DataFrame(
        {
            "gene": filtered.var_names.astype(str),
            "detected_cells": np.asarray((filtered.layers["counts"] > 0).sum(axis=0)).ravel().astype(int),
            "mean_counts": np.asarray(filtered.layers["counts"].mean(axis=0)).ravel(),
            "is_tf": [gene in tf_list for gene in filtered.var_names.astype(str)],
        }
    )
    gene_stats_path = args.metadata_dir / "driver_module6_3b_full_expression_gene_stats.tsv"
    gene_stats.to_csv(gene_stats_path, sep="\t", index=False)

    report = {
        "module": "6.3b",
        "method": "prepare full expression matrix for canonical SCENIC",
        "input_cells": str(args.driver_cells),
        "output_h5ad": str(output_h5ad),
        "n_cells": int(filtered.n_obs),
        "n_genes_before_filter": int(combined.n_vars),
        "n_genes_after_filter": int(filtered.n_vars),
        "n_tfs_after_filter": int(gene_stats["is_tf"].sum()),
        "source_manifest": str(source_manifest_path),
        "gene_stats": str(gene_stats_path),
        "elapsed_seconds": round(time.time() - start, 3),
    }
    report_path = args.metadata_dir / "driver_module6_3b_full_expression_report.json"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
