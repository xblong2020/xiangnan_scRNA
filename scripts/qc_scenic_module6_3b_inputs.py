from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="QC formal full-expression SCENIC 6.3b inputs.")
    parser.add_argument(
        "--loom",
        type=Path,
        default=ROOT / "data/processed/driver/scenic_module6_3b/driver_union_full_expression_counts.loom",
    )
    parser.add_argument(
        "--h5ad",
        type=Path,
        default=ROOT / "data/processed/driver/driver_union_full_expression.module6_3b.h5ad",
    )
    parser.add_argument(
        "--tf-list",
        type=Path,
        default=ROOT / "data/processed/driver/scenic_module6_3b/driver_union_tfs_in_matrix.txt",
    )
    parser.add_argument(
        "--driver-cells",
        type=Path,
        default=ROOT / "metadata/driver/driver_module6_1_cells.tsv.gz",
    )
    parser.add_argument(
        "--fate",
        type=Path,
        default=ROOT / "metadata/driver/driver_module6_2_cellrank_fate_probabilities.tsv.gz",
    )
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/driver")
    parser.add_argument("--reports-dir", type=Path, default=ROOT / "reports")
    parser.add_argument("--expected-cells", type=int, default=9512)
    parser.add_argument("--hvg-reference", type=int, default=2000)
    return parser.parse_args()


def make_metric(name: str, value: object, expected: object, details: str) -> dict[str, object]:
    if callable(expected):
        passed = bool(expected(value))
        expected_text = getattr(expected, "__name__", "predicate")
    elif expected is None:
        passed = True
        expected_text = "informational"
    else:
        passed = bool(value == expected)
        expected_text = expected
    return {
        "metric": name,
        "value": value,
        "expected": expected_text,
        "status": "PASS" if passed else "FAIL",
        "details": details,
    }


def classify_gene_space(n_genes: int, hvg_reference: int = 2000) -> str:
    return "full_expression" if int(n_genes) > int(hvg_reference) * 1.5 else "hvg_like"


def summarize_metadata_columns(required: list[str], available: list[str]) -> dict[str, object]:
    available_set = set(available)
    present = [column for column in required if column in available_set]
    missing = [column for column in required if column not in available_set]
    return {"present": len(present), "missing": missing}


def decode_values(values: object) -> list[str]:
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


def audit_loom(path: Path) -> dict[str, object]:
    if not path.exists():
        raise FileNotFoundError(path)
    result: dict[str, object] = {"path": str(path), "chunk_size": 256}
    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        result["matrix_shape_genes_cells"] = [int(matrix.shape[0]), int(matrix.shape[1])]
        if "Gene" not in handle["row_attrs"]:
            raise KeyError("loom row_attrs/Gene is missing")
        if "CellID" not in handle["col_attrs"]:
            raise KeyError("loom col_attrs/CellID is missing")
        genes = decode_values(handle["row_attrs"]["Gene"][:])
        cells = decode_values(handle["col_attrs"]["CellID"][:])
        result["genes"] = genes
        result["cells"] = cells
        result["n_genes"] = len(genes)
        result["n_cells"] = len(cells)
        result["unique_genes"] = len(set(genes)) == len(genes)
        result["unique_cells"] = len(set(cells)) == len(cells)
        result["empty_gene_symbols"] = sum(not gene.strip() for gene in genes)
        result["ensembl_like_genes"] = sum(gene.upper().startswith("ENSG") for gene in genes)
        result["hgnc_like_genes"] = sum(not gene.upper().startswith("ENSG") for gene in genes)
        finite = True
        negative = False
        all_zero = 0
        min_value = np.inf
        max_value = -np.inf
        for start in range(0, matrix.shape[0], int(result["chunk_size"])):
            chunk = np.asarray(matrix[start : start + int(result["chunk_size"]), :])
            finite = finite and bool(np.isfinite(chunk).all())
            negative = negative or bool((chunk < 0).any())
            all_zero += int((chunk.sum(axis=1) == 0).sum())
            if chunk.size:
                min_value = min(min_value, float(np.nanmin(chunk)))
                max_value = max(max_value, float(np.nanmax(chunk)))
        result["finite_values"] = finite
        result["negative_values"] = negative
        result["all_zero_genes"] = all_zero
        result["min_value"] = float(min_value) if np.isfinite(min_value) else None
        result["max_value"] = float(max_value) if np.isfinite(max_value) else None
    return result


def relative(root: Path, path: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve())).replace("\\", "/")
    except ValueError:
        return str(path)


def build_report(root: Path, metrics: pd.DataFrame, audit: dict[str, object], details: dict[str, object]) -> str:
    failed = metrics.loc[metrics["status"].eq("FAIL")]
    lines = [
        "# Module 6.3b SCENIC input QC",
        "",
        "## Result",
        "",
        f"- Overall status: **{'PASS' if failed.empty else 'FAIL'}**",
        f"- Formal expression space: `{details['gene_space']}`",
        f"- Loom dimensions: `{audit['n_cells']} cells x {audit['n_genes']} genes`",
        f"- TFs in matrix: `{details['n_tfs']}`; TFs missing from loom: `{details['missing_tfs']}`",
        f"- CellRank fate overlap: `{details['fate_overlap_cells']}` cells",
        f"- Known datasets in CellRank-eligible metadata: `{details['n_datasets']}`; known samples: `{details['n_samples']}`; metadata-unknown cells: `{details['unknown_metadata_cells']}`.",
        "",
        "## Metric table",
        "",
        "| Metric | Value | Expected | Status | Details |",
        "|---|---:|---:|---|---|",
    ]
    for row in metrics.itertuples(index=False):
        lines.append(f"| {row.metric} | {row.value} | {row.expected} | {row.status} | {row.details} |")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- This input is the formal 6.3b full-expression matrix and is distinct from the old 6.3 2,000-HVG exploratory analysis.",
            "- Existing old 6.3 co-expression adjacency and 6.3c outputs are historical/exploratory references until a new GRNBoost2 adjacency is generated.",
            "- The next formal stage is GRNBoost2 on this validated loom input.",
            "",
            "## Input paths",
            "",
        ]
    )
    for label, path in details["paths"].items():
        lines.append(f"- `{label}`: `{relative(root, Path(path))}`")
    return "\n".join(lines) + "\n"


def main() -> None:
    start = time.time()
    args = parse_args()
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    args.reports_dir.mkdir(parents=True, exist_ok=True)

    audit = audit_loom(args.loom)
    adata = ad.read_h5ad(args.h5ad, backed="r")
    obs = adata.obs.copy()
    obs.index = adata.obs_names.astype(str)
    h5ad_genes = adata.var_names.astype(str).tolist()
    h5ad_cells = adata.obs_names.astype(str).tolist()

    tf_list = [line.strip() for line in args.tf_list.read_text(encoding="utf-8").splitlines() if line.strip()]
    driver_cells = pd.read_csv(args.driver_cells, sep="\t")
    driver_cell_ids = set(driver_cells["cell_id"].astype(str))
    fate = pd.read_csv(args.fate, sep="\t")
    fate_cell_ids = set(fate["cell_id"].astype(str))

    required_metadata = [
        "dataset",
        "sample_id",
        "driver_main_strict__pseudotime_phase",
        "driver_main_strict__eligible",
        "driver_primary_cnv_evidence_tier",
        "cellrank_fate_prob_cnv_supported_malignant",
    ]
    metadata_summary = summarize_metadata_columns(required_metadata, list(obs.columns))
    loom_genes = set(audit["genes"])
    loom_cells = set(audit["cells"])
    tf_set = set(tf_list)
    metrics = [
        make_metric("n_cells", audit["n_cells"], args.expected_cells, "loom column count"),
        make_metric("n_genes", audit["n_genes"], lambda value: int(value) > args.hvg_reference * 1.5, "full-expression gene count"),
        make_metric("gene_space", classify_gene_space(int(audit["n_genes"]), args.hvg_reference), "full_expression", "not a 2,000-HVG matrix"),
        make_metric("loom_cell_ids_unique", audit["unique_cells"], True, "loom col_attrs/CellID"),
        make_metric("loom_gene_symbols_unique", audit["unique_genes"], True, "loom row_attrs/Gene"),
        make_metric("empty_gene_symbols", audit["empty_gene_symbols"], 0, "empty gene symbols"),
        make_metric("all_zero_genes", audit["all_zero_genes"], 0, "all-zero genes in loom"),
        make_metric("finite_expression_values", audit["finite_values"], True, "chunked matrix scan"),
        make_metric("negative_expression_values", audit["negative_values"], False, "chunked matrix scan"),
        make_metric("tf_count", len(tf_set), lambda value: int(value) > 0, "driver_union_tfs_in_matrix.txt"),
        make_metric("tf_ids_unique", len(tf_set), len(tf_list), "TF list uniqueness"),
        make_metric("tf_missing_from_loom", len(tf_set - loom_genes), 0, "TF-to-expression gene overlap"),
        make_metric("loom_h5ad_cell_overlap", len(loom_cells & set(h5ad_cells)), audit["n_cells"], "loom-to-h5ad cell mapping"),
        make_metric("loom_driver_union_cell_overlap", len(loom_cells & driver_cell_ids), audit["n_cells"], "loom-to-module6.1 driver union mapping"),
        make_metric("h5ad_obs_names_unique", adata.obs_names.is_unique, True, "h5ad obs_names"),
        make_metric("h5ad_var_names_unique", adata.var_names.is_unique, True, "h5ad var_names"),
        make_metric("required_metadata_columns", metadata_summary["present"], len(required_metadata), "h5ad obs metadata"),
        make_metric("cellrank_fate_overlap_cells", len(fate_cell_ids & set(h5ad_cells)), 5000, "fate probability to formal h5ad"),
        make_metric("fate_cell_ids_unique", len(fate_cell_ids), len(fate), "CellRank fate file"),
    ]
    metrics_df = pd.DataFrame(metrics)
    metrics_path = args.metadata_dir / "scenic_module6_3b_input_qc.tsv"
    metrics_df.to_csv(metrics_path, sep="\t", index=False)
    details = {
        "gene_space": classify_gene_space(int(audit["n_genes"]), args.hvg_reference),
        "n_tfs": len(tf_set),
        "missing_tfs": len(tf_set - loom_genes),
        "fate_overlap_cells": len(fate_cell_ids & set(h5ad_cells)),
        "n_datasets": int(obs.loc[obs["dataset"].astype(str).ne("Unknown"), "dataset"].nunique()) if "dataset" in obs else 0,
        "n_samples": int(obs.loc[obs["sample_id"].astype(str).ne("Unknown"), "sample_id"].nunique()) if "sample_id" in obs else 0,
        "unknown_metadata_cells": int(obs["dataset"].astype(str).eq("Unknown").sum()) if "dataset" in obs else 0,
        "metadata_summary": metadata_summary,
        "paths": {
            "loom": str(args.loom),
            "h5ad": str(args.h5ad),
            "tf_list": str(args.tf_list),
            "driver_cells": str(args.driver_cells),
            "fate": str(args.fate),
        },
    }
    report_path = args.reports_dir / "module6_3b_input_qc_report.md"
    report_path.write_text(build_report(ROOT, metrics_df, audit, details), encoding="utf-8")
    audit_json = {
        "module": "6.3b",
        "status": "INPUT_QC_COMPLETE" if metrics_df["status"].eq("PASS").all() else "INPUT_QC_FAILED",
        "loom": {key: value for key, value in audit.items() if key not in {"genes", "cells"}},
        "metadata": details,
        "metrics": metrics_df.to_dict(orient="records"),
        "outputs": {
            "metrics": str(metrics_path),
            "report": str(report_path),
        },
        "elapsed_seconds": round(time.time() - start, 3),
    }
    (args.metadata_dir / "driver_module6_3b_input_qc_report.json").write_text(
        json.dumps(audit_json, indent=2), encoding="utf-8"
    )
    adata.file.close()
    if not metrics_df["status"].eq("PASS").all():
        raise SystemExit("SCENIC 6.3b input QC failed; inspect the report before continuing.")


if __name__ == "__main__":
    main()
