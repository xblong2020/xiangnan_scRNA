#!/usr/bin/env python3
"""Export the HNF4A identity-high normal-reference matrix for scTenifoldKnk."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import pandas as pd
from scipy import io, sparse


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def run(h5ad_path: Path, out_dir: Path, report_path: Path, target_tf: str) -> dict:
    adata = ad.read_h5ad(h5ad_path)
    if "celloracle_state" not in adata.obs:
        raise ValueError("celloracle_state is required")
    mask = adata.obs["celloracle_state"].astype(str).eq("normal_reference")
    subset = adata[mask.to_numpy()].copy()
    if target_tf not in subset.var_names:
        raise ValueError(f"{target_tf} is absent from the identity-high matrix")
    counts = subset.layers["counts"] if "counts" in subset.layers else subset.X
    counts = counts.tocsr() if sparse.issparse(counts) else sparse.csr_matrix(counts)
    genes_x_cells = counts.transpose().tocsr()
    out_dir.mkdir(parents=True, exist_ok=True)
    matrix_path = out_dir / "figure2e_hnf4a_normal_reference_counts_genes_x_cells.mtx"
    genes_path = out_dir / "figure2e_hnf4a_normal_reference_genes.tsv"
    cells_path = out_dir / "figure2e_hnf4a_normal_reference_cells.tsv"
    # A binary handle avoids scipy path-encoding issues on Windows projects
    # whose absolute path contains non-ASCII characters.
    with matrix_path.open("wb") as handle:
        io.mmwrite(handle, genes_x_cells)
    pd.DataFrame({"gene": subset.var_names.astype(str)}).to_csv(genes_path, sep="\t", index=False)
    pd.DataFrame({"cell_id": subset.obs_names.astype(str)}).to_csv(cells_path, sep="\t", index=False)
    if not matrix_path.exists() or matrix_path.stat().st_size == 0:
        raise RuntimeError(f"Matrix Market export failed: {matrix_path}")
    report = {
        "module": "Figure 2E HNF4A scTenifoldKnk input",
        "target_tf": target_tf,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "subset": "normal_reference",
        "selection": "celloracle_state == normal_reference",
        "matrix_orientation": "genes_x_cells",
        "n_cells": int(subset.n_obs),
        "n_genes": int(subset.n_vars),
        "input_h5ad": str(h5ad_path.resolve()),
        "input_layer": "counts" if "counts" in subset.layers else "X",
        "outputs": {
            "matrix": str(matrix_path.resolve()),
            "genes": str(genes_path.resolve()),
            "cells": str(cells_path.resolve()),
        },
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-tf", default="HNF4A")
    parser.add_argument("--h5ad", type=Path, default=PROJECT_ROOT / "data/processed/driver/celloracle_module6_6/celloracle_module6_6_input.h5ad")
    parser.add_argument("--out-dir", type=Path, default=PROJECT_ROOT / "data/processed/driver/figure2e_hnf4a_sctenifoldknk/normal_reference")
    parser.add_argument("--report", type=Path, default=PROJECT_ROOT / "metadata/driver/figure2e_hnf4a_sctenifoldknk/figure2e_hnf4a_normal_reference_input_report.json")
    args = parser.parse_args()
    print(json.dumps(run(args.h5ad, args.out_dir, args.report, args.target_tf), ensure_ascii=False))


if __name__ == "__main__":
    main()
