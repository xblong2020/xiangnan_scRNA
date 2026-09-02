from __future__ import annotations

import argparse
import json
import platform
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io, sparse


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_H5AD = PROJECT_ROOT / "data/processed/driver/celloracle_module6_6/celloracle_module6_6_input.h5ad"
DEFAULT_METADATA_DIR = PROJECT_ROOT / "metadata/driver"
DEFAULT_OUT_DIR = PROJECT_ROOT / "data/processed/driver/sctenifoldknk_module7_1"
DEFAULT_TF_LIST = DEFAULT_METADATA_DIR / "celloracle_input_tfs.module6_4.txt"

CELLoRACLE_TFS = [
    "ATF3",
    "CEBPB",
    "EGR1",
    "FOS",
    "HLF",
    "HNF4A",
    "IRF1",
    "JUN",
    "JUNB",
    "JUND",
    "MAFB",
    "MAFF",
    "MYC",
    "PPARA",
    "SOX4",
]


def read_tf_list(path: Path | None = None) -> list[str]:
    if path is None or not path.exists():
        return CELLoRACLE_TFS.copy()
    tfs = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    return tfs


def validate_tf_coverage(genes: Iterable[str], tfs: Iterable[str]) -> dict:
    gene_set = {str(gene) for gene in genes}
    tf_list = [str(tf) for tf in tfs]
    retained = [tf for tf in tf_list if tf in gene_set]
    missing = [tf for tf in tf_list if tf not in gene_set]
    summary = {
        "n_input_tfs": int(len(tf_list)),
        "n_retained_tfs": int(len(retained)),
        "retained_tfs": retained,
        "missing_tfs": missing,
    }
    if missing:
        raise ValueError(f"Missing required scTenifoldKnk knockout TFs: {', '.join(missing)}")
    return summary


def build_cell_subset_mask(obs: pd.DataFrame, subset_name: str) -> pd.Series:
    subset = subset_name.lower()
    if subset in {"driver_union_all", "all", "driver_union_all_cells"}:
        return pd.Series(True, index=obs.index)
    if subset in {"malignant_like", "malignant_or_malignant_like"}:
        if "celloracle_state" not in obs.columns:
            raise ValueError("celloracle_state is required for malignant_like subset")
        return obs["celloracle_state"].astype(str).eq("malignant_or_malignant_like")
    if subset in {"main_strict", "main_strict_cells"}:
        if "celloracle_main_strict" not in obs.columns:
            raise ValueError("celloracle_main_strict is required for main_strict subset")
        return obs["celloracle_main_strict"].astype(bool)
    raise ValueError(f"Unsupported Module 7.1 subset: {subset_name}")


def _as_csr(matrix) -> sparse.csr_matrix:
    if sparse.issparse(matrix):
        return matrix.tocsr()
    return sparse.csr_matrix(matrix)


def _write_matrix_market(matrix: sparse.spmatrix, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        io.mmwrite(handle, matrix)


def _write_index(values: Iterable[str], path: Path, column: str) -> None:
    pd.DataFrame({column: list(values)}).to_csv(path, sep="\t", index=False)


def build_export_report(
    subset_name: str,
    n_genes: int,
    n_cells: int,
    retained_tfs: list[str],
    output_files: dict,
    extra: dict | None = None,
) -> dict:
    report = {
        "module": "7.1",
        "method": "Export CellOracle prepared expression matrix for scTenifoldKnk",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "subset": subset_name,
        "matrix_orientation": "genes_x_cells",
        "shape": {"n_genes": int(n_genes), "n_cells": int(n_cells)},
        "n_retained_tfs": int(len(retained_tfs)),
        "retained_tfs": list(retained_tfs),
        "outputs": output_files,
        "python_runtime": {
            "version": platform.python_version(),
            "platform": platform.platform(),
        },
    }
    if extra:
        report.update(extra)
    return report


def export_sctenifoldknk_input(
    h5ad_path: Path,
    tf_list_path: Path | None,
    out_dir: Path,
    subset_name: str,
    layer: str = "counts",
) -> dict:
    adata = ad.read_h5ad(h5ad_path)
    tfs = read_tf_list(tf_list_path)
    validate = validate_tf_coverage(adata.var_names.astype(str), tfs)
    mask = build_cell_subset_mask(adata.obs, subset_name)
    if not bool(mask.any()):
        raise ValueError(f"Subset {subset_name!r} has no cells")

    subset = adata[mask.to_numpy(), :].copy()
    if layer in subset.layers:
        counts = _as_csr(subset.layers[layer])
    elif layer == "X":
        counts = _as_csr(subset.X)
    else:
        raise ValueError(f"Layer {layer!r} not found in h5ad")

    genes_x_cells = counts.transpose().tocsr()
    log_norm = genes_x_cells.copy().astype(np.float64)
    log_norm.data = np.log1p(log_norm.data)

    subset_dir = out_dir / subset_name
    files = {
        "counts_mtx": str(subset_dir / "sctenifoldknk_counts_genes_x_cells.mtx"),
        "log1p_mtx": str(subset_dir / "sctenifoldknk_log1p_genes_x_cells.mtx"),
        "genes": str(subset_dir / "sctenifoldknk_genes.tsv"),
        "cells": str(subset_dir / "sctenifoldknk_cells.tsv"),
        "cell_metadata": str(subset_dir / "sctenifoldknk_cell_metadata.tsv"),
        "tf_list": str(subset_dir / "sctenifoldknk_celloracle_tfs.txt"),
        "report": str(DEFAULT_METADATA_DIR / f"sctenifoldknk_module7_1_{subset_name}_export_report.json"),
    }
    _write_matrix_market(genes_x_cells, Path(files["counts_mtx"]))
    _write_matrix_market(log_norm, Path(files["log1p_mtx"]))
    _write_index(subset.var_names.astype(str), Path(files["genes"]), "gene")
    _write_index(subset.obs_names.astype(str), Path(files["cells"]), "cell_id")
    subset.obs.copy().assign(cell_id=subset.obs_names.astype(str)).to_csv(files["cell_metadata"], sep="\t", index=False)
    Path(files["tf_list"]).write_text("\n".join(tfs) + "\n", encoding="utf-8")

    report = build_export_report(
        subset_name=subset_name,
        n_genes=genes_x_cells.shape[0],
        n_cells=genes_x_cells.shape[1],
        retained_tfs=validate["retained_tfs"],
        output_files=files,
        extra={
            "input_h5ad": str(h5ad_path),
            "input_layer": layer,
            "n_missing_tfs": int(len(validate["missing_tfs"])),
        },
    )
    Path(files["report"]).parent.mkdir(parents=True, exist_ok=True)
    Path(files["report"]).write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 7.1 export CellOracle h5ad for scTenifoldKnk")
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--tf-list", type=Path, default=DEFAULT_TF_LIST)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--subset", default="driver_union_all", choices=["driver_union_all", "malignant_like", "main_strict"])
    parser.add_argument("--layer", default="counts")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    report = export_sctenifoldknk_input(
        h5ad_path=args.h5ad,
        tf_list_path=args.tf_list,
        out_dir=args.out_dir,
        subset_name=args.subset,
        layer=args.layer,
    )
    print(json.dumps({"report": report["outputs"]["report"], "shape": report["shape"]}, indent=2))


if __name__ == "__main__":
    main()
