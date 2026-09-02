#!/usr/bin/env python3
"""Audit every input required for Figure 3 EGR1 before plotting or rerunning models."""

from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

try:
    from figure3_egr1_common import (
        PROJECT_ROOT,
        PSEUDOTIME_COLUMN,
        SEED,
        STATE_COLUMN,
        STATE_ORDER,
        STRICT_COLUMN,
        TARGET_TF,
        as_bool,
        fingerprint_paths,
        json_safe,
        write_json,
    )
except ModuleNotFoundError:
    from scripts.figure3_egr1_common import (
        PROJECT_ROOT,
        PSEUDOTIME_COLUMN,
        SEED,
        STATE_COLUMN,
        STATE_ORDER,
        STRICT_COLUMN,
        TARGET_TF,
        as_bool,
        fingerprint_paths,
        json_safe,
        write_json,
    )


DEFAULT_H5AD = PROJECT_ROOT / "data/processed/driver/celloracle_module6_6/celloracle_module6_6_input.h5ad"
DEFAULT_FITTED_ORACLE = (
    PROJECT_ROOT / "data/processed/driver/celloracle_module6_7/celloracle_module6_7_fitted.celloracle.oracle"
)
DEFAULT_OUT_DIR = PROJECT_ROOT / "metadata/driver/figure3_egr1_preflight"


def add_check(
    checks: list[dict],
    check_id: str,
    category: str,
    status: str,
    value,
    details: str,
    source: Path | str | None = None,
) -> None:
    if status not in {"pass", "warning", "fail", "info"}:
        raise ValueError(f"Invalid preflight status: {status}")
    checks.append(
        {
            "check_id": check_id,
            "category": category,
            "status": status,
            "value": json.dumps(json_safe(value), ensure_ascii=False) if isinstance(value, (dict, list)) else value,
            "details": details,
            "source": str(Path(source).resolve()) if source else "",
        }
    )


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def protected_paths(root: Path) -> list[Path]:
    paths: list[Path] = []
    for top in ["scripts", "data", "metadata", "figures", "reports"]:
        base = root / top
        if not base.exists():
            continue
        for path in base.rglob("*"):
            if path.is_file() and any(token in str(path).lower() for token in ("sox4", "hnf4a")):
                paths.append(path)
    return paths


def state_expression_table(adata: ad.AnnData, gene: str) -> pd.DataFrame:
    if gene not in adata.var_names:
        return pd.DataFrame(
            columns=["celloracle_state", "n_cells", "mean_expression", "detection_rate", "n_detected"]
        )
    gene_ix = int(adata.var_names.get_loc(gene))
    matrix = adata.layers["counts"] if "counts" in adata.layers else adata.X
    column = matrix[:, gene_ix]
    values = column.toarray().ravel() if sparse.issparse(column) else np.asarray(column).ravel()
    frame = pd.DataFrame(
        {
            "celloracle_state": adata.obs[STATE_COLUMN].astype(str).to_numpy(),
            "expression": values.astype(float),
        }
    )
    rows = []
    for state in STATE_ORDER:
        sub = frame.loc[frame["celloracle_state"].eq(state), "expression"]
        rows.append(
            {
                "celloracle_state": state,
                "n_cells": int(len(sub)),
                "mean_expression": float(sub.mean()) if len(sub) else None,
                "detection_rate": float(sub.gt(0).mean()) if len(sub) else None,
                "n_detected": int(sub.gt(0).sum()),
            }
        )
    return pd.DataFrame(rows)


def subset_network_inventory(root: Path) -> pd.DataFrame:
    rows = []
    input_root = root / "data/processed/driver/sctenifoldknk_module7_1"
    result_root = root / "data/processed/driver/sctenifoldknk_module7_2"
    for subset_dir in sorted(input_root.glob("*")):
        if not subset_dir.is_dir():
            continue
        cells_path = subset_dir / "sctenifoldknk_cells.tsv"
        genes_path = subset_dir / "sctenifoldknk_genes.tsv"
        metadata_path = subset_dir / "sctenifoldknk_cell_metadata.tsv"
        n_cells = sum(1 for _ in cells_path.open("r", encoding="utf-8")) - 1 if cells_path.exists() else 0
        genes = (
            set(pd.read_csv(genes_path, sep="\t").iloc[:, 0].astype(str))
            if genes_path.exists()
            else set()
        )
        egr1_result = result_root / subset_dir.name / f"sctenifoldknk_{subset_dir.name}_EGR1_perturbation_genes.tsv"
        report_path = root / f"metadata/driver/sctenifoldknk_module7_2_{subset_dir.name}_report.json"
        report = read_json(report_path) if report_path.exists() else {}
        rows.append(
            {
                "subset": subset_dir.name,
                "n_cells": int(max(n_cells, 0)),
                "n_genes": int(len(genes)),
                "egr1_in_network_genes": TARGET_TF in genes,
                "egr1_result_exists": egr1_result.exists() and egr1_result.stat().st_size > 0,
                "cell_metadata_exists": metadata_path.exists(),
                "nc_nNet": report.get("parameters", {}).get("nc_nNet"),
                "nc_nCells": report.get("parameters", {}).get("nc_nCells"),
                "seed": report.get("parameters", {}).get("seed"),
                "input_dir": str(subset_dir.resolve()),
                "egr1_result": str(egr1_result.resolve()),
            }
        )
    return pd.DataFrame(rows)


def candidate_inventory(root: Path) -> pd.DataFrame:
    candidates = ["JUN", "JUNB", "JUND", "FOS", "ATF3", "CEBPB", "EGR1"]
    quantitative = pd.read_csv(
        root / "metadata/driver/celloracle_module6_9b_quantitative_tf_scores.tsv", sep="\t"
    )
    stability = pd.read_csv(root / "metadata/driver/celloracle_module6_10_rank_stability.tsv", sep="\t")
    phase = pd.read_csv(root / "metadata/driver/celloracle_module6_10_phase_wide_summary.tsv", sep="\t")
    integrated = pd.read_csv(
        root / "metadata/driver/sctenifoldknk_module7_3_integrated_evidence_matrix.tsv", sep="\t"
    )
    selection = pd.read_csv(root / "metadata/driver/celloracle_tf_selection.module6_4.tsv", sep="\t")
    keep_selection = [
        "tf",
        "rank_by_total_score",
        "detection_rate_main",
        "mean_expression_main",
        "phase_early_mean_auc",
        "phase_middle_mean_auc",
        "phase_late_mean_auc",
        "loo_min_directional_r",
        "dataset_direction_consistency_fraction",
    ]
    keep_integrated = [
        "tf",
        "integrated_rank",
        "module7_integrated_rank",
        "scTenifoldKnk_rank",
        "n_significant_perturbed_genes",
        "top5_lodo_fraction",
        "proliferation_dependency" if "proliferation_dependency" in integrated.columns else None,
    ]
    keep_integrated = [column for column in keep_integrated if column]
    out = (
        pd.DataFrame({"tf": candidates})
        .merge(selection[keep_selection], on="tf", how="left")
        .merge(
            quantitative[
                [
                    "tf",
                    "quantitative_rank",
                    "quantitative_perturbation_score",
                    "state_specificity_ratio",
                    "proliferation_module_rescue_score",
                ]
            ],
            on="tf",
            how="left",
        )
        .merge(stability, on="tf", how="left")
        .merge(phase, on="tf", how="left")
        .merge(integrated[keep_integrated], on="tf", how="left")
    )
    phase_score_cols = ["phase_early_score", "phase_intermediate_score", "phase_late_score"]
    out["peak_phase"] = out[phase_score_cols].idxmax(axis=1).str.replace("phase_", "").str.replace("_score", "")
    out["peak_pseudotime_proxy"] = out["peak_phase"].map(
        {"early": 1 / 6, "intermediate": 0.5, "late": 5 / 6}
    )
    out["generic_stress_risk"] = out["tf"].isin(["JUN", "JUNB", "JUND", "FOS", "ATF3"]).astype(float)
    return out


def run(h5ad_path: Path, fitted_oracle: Path, out_dir: Path) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    checks: list[dict] = []
    review_risks: list[dict] = []

    protected = fingerprint_paths(protected_paths(PROJECT_ROOT))
    protected_path = out_dir / "protected_sox4_hnf4a_fingerprints_before.tsv"
    protected.to_csv(protected_path, sep="\t", index=False)
    add_check(
        checks,
        "protected_asset_manifest",
        "isolation",
        "pass" if len(protected) else "warning",
        int(len(protected)),
        "SHA-256 fingerprints recorded before any Figure 3 execution.",
        protected_path,
    )

    add_check(
        checks,
        "fitted_oracle_exists",
        "celloracle",
        "pass" if fitted_oracle.exists() and fitted_oracle.stat().st_size > 0 else "fail",
        fitted_oracle.stat().st_size if fitted_oracle.exists() else 0,
        "Saved fitted CellOracle Oracle object required for reproducible EGR1 KO.",
        fitted_oracle,
    )

    adata = ad.read_h5ad(h5ad_path)
    required_obs = [STRICT_COLUMN, PSEUDOTIME_COLUMN, STATE_COLUMN]
    missing_obs = [column for column in required_obs if column not in adata.obs.columns]
    add_check(
        checks,
        "required_h5ad_columns",
        "baseline",
        "pass" if not missing_obs and "X_celloracle_umap" in adata.obsm else "fail",
        {"obs_missing": missing_obs, "cell_id_source": "obs_names", "obsm_X_celloracle_umap": "X_celloracle_umap" in adata.obsm},
        "cell_id is represented by AnnData obs_names; required strict, pseudotime, state, and UMAP fields were audited.",
        h5ad_path,
    )
    pseudotime = pd.to_numeric(adata.obs[PSEUDOTIME_COLUMN], errors="coerce")
    strict_mask = as_bool(adata.obs[STRICT_COLUMN]) & pseudotime.notna()
    add_check(
        checks,
        "strict_5000_cells",
        "baseline",
        "pass" if int(strict_mask.sum()) == 5000 else "fail",
        int(strict_mask.sum()),
        "Figure 3B-F must use the exact strict-main 5,000-cell subset.",
        h5ad_path,
    )

    existing_baseline = PROJECT_ROOT / "metadata/driver/figure2b_sox4/figure2b_sox4_plot_cells.tsv.gz"
    if existing_baseline.exists():
        baseline = pd.read_csv(existing_baseline, sep="\t")
        strict_ids = adata.obs_names[strict_mask.to_numpy()].astype(str).tolist()
        identical_ids = baseline["cell_id"].astype(str).tolist() == strict_ids
        add_check(
            checks,
            "common_baseline_cell_order",
            "baseline",
            "pass" if identical_ids else "fail",
            {"n_existing": int(len(baseline)), "n_strict": int(len(strict_ids)), "identical_order": identical_ids},
            "Existing verified SOX4 baseline is eligible for exact reuse only if cell order is identical.",
            existing_baseline,
        )
    else:
        add_check(
            checks,
            "common_baseline_cell_order",
            "baseline",
            "fail",
            False,
            "Existing common baseline cell table is missing.",
            existing_baseline,
        )

    grn_links_path = PROJECT_ROOT / "metadata/driver/celloracle_module6_7_grn_links_filtered.tsv.gz"
    grn_links = pd.read_csv(grn_links_path, sep="\t", usecols=["celloracle_state", "source", "target"])
    egr1_links = grn_links.loc[grn_links["source"].astype(str).eq(TARGET_TF)]
    add_check(
        checks,
        "egr1_in_fitted_grn",
        "celloracle",
        "pass" if len(egr1_links) else "fail",
        {"n_links": int(len(egr1_links)), "n_states": int(egr1_links["celloracle_state"].nunique())},
        "EGR1 must have fitted outgoing regulatory links.",
        grn_links_path,
    )

    perturb_path = PROJECT_ROOT / "metadata/driver/celloracle_module6_8_cell_shift_summary.tsv.gz"
    perturb = pd.read_csv(perturb_path, sep="\t")
    egr1 = perturb.loc[perturb["tf"].astype(str).eq(TARGET_TF)].copy()
    delta_columns = ["delta_embedding_1", "delta_embedding_2"]
    complete_delta = (
        len(egr1) == adata.n_obs
        and not egr1["cell_id"].astype(str).duplicated().any()
        and np.isfinite(egr1[delta_columns].to_numpy(dtype=float)).all()
    )
    add_check(
        checks,
        "egr1_cell_level_delta_embedding",
        "celloracle",
        "pass" if complete_delta else "fail",
        {"n_rows": int(len(egr1)), "n_unique_cells": int(egr1["cell_id"].nunique()), "finite": bool(complete_delta)},
        "True saved EGR1 CellOracle displacement, never expression-derived displacement.",
        perturb_path,
    )

    grid_path = PROJECT_ROOT / "metadata/driver/celloracle_module6_8_grid_arrows.tsv.gz"
    grid = pd.read_csv(grid_path, sep="\t")
    egr1_grid = grid.loc[grid["tf"].astype(str).eq(TARGET_TF)]
    add_check(
        checks,
        "egr1_saved_grid_arrows",
        "celloracle",
        "pass" if len(egr1_grid) else "warning",
        int(len(egr1_grid)),
        "Saved native CellOracle grid arrows are retained for provenance; Figure 3C reaggregates cells on the exact common baseline grid.",
        grid_path,
    )

    perturb_report_path = PROJECT_ROOT / "metadata/driver/celloracle_module6_8_perturbation_report.json"
    perturb_report = read_json(perturb_report_path)
    parameters = perturb_report.get("parameters", {})
    per_tf = {
        row.get("tf"): row
        for row in perturb_report.get("result", {}).get("per_tf_reports", [])
    }
    same_run = TARGET_TF in per_tf and "SOX4" in per_tf and parameters.get("seed") == SEED
    add_check(
        checks,
        "egr1_sox4_simulation_parameter_equivalence",
        "celloracle",
        "pass" if same_run else "fail",
        parameters,
        "EGR1 and SOX4 were simulated in the same Module 6.8 run under one shared parameter contract.",
        perturb_report_path,
    )

    expression = state_expression_table(adata, TARGET_TF)
    expression_path = out_dir / "figure3_egr1_state_expression_detection.tsv"
    expression.to_csv(expression_path, sep="\t", index=False)
    for row in expression.to_dict(orient="records"):
        add_check(
            checks,
            f"egr1_expression_{row['celloracle_state']}",
            "expression",
            "pass" if row["n_cells"] > 0 else "warning",
            row,
            "Raw-count mean expression and detection rate by CellOracle state.",
            expression_path,
        )

    scenic_cols = [column for column in adata.obs.columns if "egr1" in column.lower() and "auc" in column.lower()]
    add_check(
        checks,
        "egr1_regulon_auc",
        "scenic",
        "pass" if scenic_cols else "warning",
        scenic_cols,
        "No EGR1 regulon AUC may be inferred from unrelated regulons; missing EGR1 AUC is retained as missing evidence.",
        h5ad_path,
    )
    if not scenic_cols:
        review_risks.append(
            {
                "flag": "egr1_scenic_auc_missing",
                "severity": "review_attention",
                "detail": "The prepared h5ad contains no EGR1-specific SCENIC/cisTarget regulon AUC column.",
            }
        )

    candidate = candidate_inventory(PROJECT_ROOT)
    candidate_path = out_dir / "figure3_egr1_candidate_preflight_metrics.tsv"
    candidate.to_csv(candidate_path, sep="\t", index=False)
    add_check(
        checks,
        "candidate_evidence_inventory",
        "selection",
        "pass" if candidate["tf"].eq(TARGET_TF).any() else "fail",
        {"n_candidates": int(len(candidate)), "columns": candidate.columns.tolist()},
        "CellOracle, integrated, phase, LODO, state-specificity, and perturbation-gene evidence gathered for AP-1/CEBPB/EGR1 candidates.",
        candidate_path,
    )

    networks = subset_network_inventory(PROJECT_ROOT)
    network_path = out_dir / "figure3_egr1_sctenifoldknk_network_inventory.tsv"
    networks.to_csv(network_path, sep="\t", index=False)
    add_check(
        checks,
        "egr1_in_sctenifoldknk_networks",
        "sctenifoldknk",
        "pass" if len(networks) and networks["egr1_in_network_genes"].all() else "fail",
        networks.to_dict(orient="records"),
        "Existing scTenifoldKnk networks and EGR1 result contracts inventoried by subset.",
        network_path,
    )
    low_parameter = networks.loc[
        pd.to_numeric(networks["nc_nNet"], errors="coerce").fillna(0).lt(10)
        | pd.to_numeric(networks["nc_nCells"], errors="coerce").fillna(0).lt(500)
    ]
    if len(low_parameter):
        review_risks.append(
            {
                "flag": "existing_sctenifoldknk_low_replication",
                "severity": "review_attention",
                "detail": "All existing Module 7.2 subsets use nc_nNet=1 and nc_nCells=100; Figure 3E needs a dedicated higher-replication stress-transition run.",
            }
        )

    checks_df = pd.DataFrame(checks)
    checks_path = out_dir / "figure3_egr1_preflight_report.tsv"
    checks_df.to_csv(checks_path, sep="\t", index=False)
    fail_count = int(checks_df["status"].eq("fail").sum())
    warning_count = int(checks_df["status"].eq("warning").sum())
    report = {
        "module": "Figure 3 EGR1 preflight",
        "target_tf": TARGET_TF,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "overall_status": "fail" if fail_count else ("warning" if warning_count else "pass"),
        "n_checks": int(len(checks_df)),
        "n_fail": fail_count,
        "n_warning": warning_count,
        "inputs": {
            "h5ad": str(h5ad_path.resolve()),
            "fitted_oracle": str(fitted_oracle.resolve()),
            "celloracle_perturbation": str(perturb_path.resolve()),
            "celloracle_grid_arrows": str(grid_path.resolve()),
            "celloracle_report": str(perturb_report_path.resolve()),
        },
        "celloracle_simulation_parameters": parameters,
        "egr1_celloracle_source": {
            "n_cells": int(len(egr1)),
            "delta_embedding_columns": delta_columns,
            "saved_grid_rows": int(len(egr1_grid)),
            "condition": per_tf.get(TARGET_TF, {}).get("condition"),
            "source_oracle": perturb_report.get("result", {}).get("input_oracle"),
        },
        "strict_main": {
            "n_cells": int(strict_mask.sum()),
            "pseudotime_column": PSEUDOTIME_COLUMN,
            "state_counts": {
                str(key): int(value)
                for key, value in adata.obs.loc[strict_mask, STATE_COLUMN].astype(str).value_counts().items()
            },
        },
        "state_expression_detection": expression.to_dict(orient="records"),
        "scenic_egr1_auc_columns": scenic_cols,
        "candidate_metric_table": str(candidate_path.resolve()),
        "network_inventory": networks.to_dict(orient="records"),
        "review_risk_flags": review_risks,
        "outputs": {
            "checks_tsv": str(checks_path.resolve()),
            "protected_fingerprints": str(protected_path.resolve()),
            "state_expression": str(expression_path.resolve()),
            "candidate_metrics": str(candidate_path.resolve()),
            "network_inventory": str(network_path.resolve()),
        },
        "runtime": {
            "python": sys.version,
            "executable": sys.executable,
            "platform": platform.platform(),
            "anndata": ad.__version__,
            "pandas": pd.__version__,
            "numpy": np.__version__,
        },
        "caveat": "Preflight verifies computational availability and provenance. It does not establish experimental causality for EGR1.",
    }
    report = json_safe(report)
    report_path = out_dir / "figure3_egr1_preflight_report.json"
    write_json(report, report_path)
    adata.file.close() if getattr(adata, "file", None) and adata.file.is_open else None
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--fitted-oracle", type=Path, default=DEFAULT_FITTED_ORACLE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = run(args.h5ad, args.fitted_oracle, args.out_dir)
    print(
        json.dumps(
            {
                "overall_status": report["overall_status"],
                "n_checks": report["n_checks"],
                "n_fail": report["n_fail"],
                "n_warning": report["n_warning"],
                "report": str((args.out_dir / "figure3_egr1_preflight_report.json").resolve()),
            },
            ensure_ascii=False,
        )
    )
    return 0 if report["n_fail"] == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())

