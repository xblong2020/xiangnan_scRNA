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
    parser = argparse.ArgumentParser(description="Create a preserved formal 6.3b input without all-zero genes.")
    parser.add_argument(
        "--input-h5ad",
        type=Path,
        default=ROOT / "data/processed/driver/driver_union_full_expression.module6_3b.h5ad",
    )
    parser.add_argument(
        "--output-h5ad",
        type=Path,
        default=ROOT / "data/processed/driver/driver_union_full_expression.module6_3b.formal.h5ad",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "data/processed/driver/scenic_module6_3b_formal",
    )
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/driver")
    return parser.parse_args()


def counts_matrix(adata: ad.AnnData):
    matrix = adata.layers["counts"] if "counts" in adata.layers else adata.X
    return matrix.tocsr() if sparse.issparse(matrix) else sparse.csr_matrix(np.asarray(matrix))


def main() -> None:
    start = time.time()
    args = parse_args()
    args.output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.input_h5ad)
    counts = counts_matrix(adata)
    detected = np.asarray((counts > 0).sum(axis=0)).ravel()
    keep = detected > 0
    removed = adata.var_names[~keep].astype(str).tolist()
    formal = adata[:, keep].copy()
    formal.layers["counts"] = counts[:, keep].copy()
    formal.uns["module6_3b_formal_input"] = {
        "module": "6.3b",
        "method": "preserved full-expression checkpoint cleaned of all-zero genes",
        "source_h5ad": str(args.input_h5ad),
        "n_cells": int(formal.n_obs),
        "n_genes_before": int(adata.n_vars),
        "n_genes_after": int(formal.n_vars),
        "removed_all_zero_genes": int(len(removed)),
        "removed_gene_examples": removed[:100],
    }
    formal.write_h5ad(args.output_h5ad, compression="gzip")

    genes = pd.Index(formal.var_names.astype(str))
    cells = pd.Index(formal.obs_names.astype(str))
    matrix = formal.layers["counts"]
    matrix_t = matrix.T.tocsr()
    loom_path = args.output_dir / "driver_union_full_expression_counts.loom"
    import loompy

    loompy.create(
        str(loom_path),
        layers=matrix_t,
        row_attrs={"Gene": genes.to_numpy()},
        col_attrs={"CellID": cells.to_numpy()},
    )
    tf_catalog = ROOT / "metadata/driver/scenic_resources/allTFs_hg38.txt"
    tf_set = {line.strip() for line in tf_catalog.read_text(encoding="utf-8").splitlines() if line.strip()}
    tf_in_matrix = [gene for gene in genes if gene in tf_set]
    tf_path = args.output_dir / "driver_union_tfs_in_matrix.txt"
    tf_path.write_text("\n".join(tf_in_matrix) + "\n", encoding="utf-8")

    stats_path = args.metadata_dir / "driver_module6_3b_formal_input_gene_stats.tsv"
    stats = pd.DataFrame(
        {
            "gene": genes,
            "detected_cells": np.asarray((matrix > 0).sum(axis=0)).ravel().astype(int),
            "mean_counts": np.asarray(matrix.mean(axis=0)).ravel(),
            "is_tf": [gene in tf_set for gene in genes],
        }
    )
    stats.to_csv(stats_path, sep="\t", index=False)
    report = {
        "module": "6.3b",
        "method": "preserved formal full-expression input cleanup and loom export",
        "source_h5ad": str(args.input_h5ad),
        "output_h5ad": str(args.output_h5ad),
        "loom": str(loom_path),
        "tf_file": str(tf_path),
        "n_cells": int(formal.n_obs),
        "n_genes_before": int(adata.n_vars),
        "n_genes_after": int(formal.n_vars),
        "n_removed_all_zero_genes": int(len(removed)),
        "n_tfs": int(len(tf_in_matrix)),
        "removed_gene_examples": removed[:100],
        "gene_stats": str(stats_path),
        "elapsed_seconds": round(time.time() - start, 3),
    }
    report_path = args.metadata_dir / "driver_module6_3b_formal_input_report.json"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
