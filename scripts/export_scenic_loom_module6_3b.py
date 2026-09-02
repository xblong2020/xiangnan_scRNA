from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import anndata as ad
import loompy
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 6.3b: export full-expression driver union to SCENIC loom.")
    parser.add_argument(
        "--input-h5ad",
        type=Path,
        default=ROOT / "data/processed/driver/driver_union_full_expression.module6_3b.h5ad",
    )
    parser.add_argument("--tf-list", type=Path, default=ROOT / "metadata/driver/scenic_resources/allTFs_hg38.txt")
    parser.add_argument("--processed-dir", type=Path, default=ROOT / "data/processed/driver/scenic_module6_3b")
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/driver")
    return parser.parse_args()


def main() -> None:
    start = time.time()
    args = parse_args()
    args.processed_dir.mkdir(parents=True, exist_ok=True)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.input_h5ad)
    genes = pd.Index(adata.var_names.astype(str))
    cells = pd.Index(adata.obs_names.astype(str))
    matrix = adata.layers["counts"] if "counts" in adata.layers else adata.X
    matrix = matrix.T.tocsr() if sparse.issparse(matrix) else sparse.csr_matrix(np.asarray(matrix).T)

    loom_path = args.processed_dir / "driver_union_full_expression_counts.loom"
    if loom_path.exists():
        loom_path.unlink()
    loompy.create(
        str(loom_path),
        layers=matrix,
        row_attrs={"Gene": genes.to_numpy()},
        col_attrs={"CellID": cells.to_numpy()},
    )

    tf_set = set(args.tf_list.read_text(encoding="utf-8").splitlines())
    tf_in_matrix = [gene for gene in genes if gene in tf_set]
    tf_path = args.processed_dir / "driver_union_tfs_in_matrix.txt"
    tf_path.write_text("\n".join(tf_in_matrix) + "\n", encoding="utf-8")

    report = {
        "module": "6.3b",
        "method": "export SCENIC loom input",
        "input_h5ad": str(args.input_h5ad),
        "loom": str(loom_path),
        "tf_file": str(tf_path),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_tfs": int(len(tf_in_matrix)),
        "elapsed_seconds": round(time.time() - start, 3),
    }
    report_path = args.metadata_dir / "driver_module6_3b_scenic_loom_export_report.json"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
