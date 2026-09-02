from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]

CNV_SUPPORTED_TIERS = {"module3_high_conf_malignant", "module3_cnv_supported_malignant"}
TRUE_VALUES = {"1", "true", "t", "yes", "y", "aneuploid"}
FALSE_VALUES = {"0", "false", "f", "no", "n", "diploid", "unknown", "not.defined", "not_run", "nan", ""}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 5.5: overlay CNV and malignant evidence along trajectory.")
    parser.add_argument(
        "--evidence-cells",
        type=Path,
        default=ROOT / "metadata/trajectory/trajectory_module5_2_stage_root_end_cells.tsv.gz",
    )
    parser.add_argument(
        "--module-scores",
        type=Path,
        default=ROOT / "metadata/trajectory/trajectory_module5_4_module_scores_by_cell.tsv.gz",
    )
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/trajectory")
    parser.add_argument("--figures-dir", type=Path, default=ROOT / "figures/trajectory")
    parser.add_argument("--n-bins", type=int, default=10)
    return parser.parse_args()


def pseudotime_inputs(metadata_dir: Path) -> list[dict[str, object]]:
    return [
        {
            "run_id": "main_strict",
            "path": metadata_dir / "trajectory_module5_3_main_strict_pseudotime_merged.tsv.gz",
            "methods": {
                "monocle3": "main_strict__monocle3_norm",
                "slingshot_scanvi": "main_strict__slingshot_scanvi_norm",
                "slingshot_hepatocyte_pca": "main_strict__slingshot_hepatocyte_pca_norm",
            },
        },
        {
            "run_id": "sensitivity_include_review",
            "path": metadata_dir / "trajectory_module5_3_sensitivity_include_review_pseudotime_merged.tsv.gz",
            "methods": {
                "monocle3": "sensitivity_include_review__monocle3_norm",
                "slingshot_scanvi": "sensitivity_include_review__slingshot_scanvi_norm",
                "slingshot_hepatocyte_pca": "sensitivity_include_review__slingshot_hepatocyte_pca_norm",
            },
        },
    ]


def normalize_label(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip().lower()


def boolish(value: object) -> float:
    if pd.isna(value):
        return np.nan
    if isinstance(value, (bool, np.bool_)):
        return float(value)
    if isinstance(value, (int, float, np.integer, np.floating)) and np.isfinite(value):
        return float(value != 0)
    label = normalize_label(value)
    if label in TRUE_VALUES:
        return 1.0
    if label in FALSE_VALUES:
        return 0.0
    return np.nan


def numeric_fraction(values: pd.Series) -> float:
    numeric = values.map(boolish).replace([np.inf, -np.inf], np.nan).dropna()
    if numeric.empty:
        return float("nan")
    return float(numeric.mean())


def assign_cnv_evidence_tier(row: dict[str, object] | pd.Series) -> str:
    call = normalize_label(row.get("malignant_hcc_call", ""))
    role = normalize_label(row.get("trajectory_root_end_role", ""))
    copykat_status = normalize_label(row.get("copykat_status", ""))
    copykat_pred = normalize_label(row.get("copykat_pred", ""))
    cnv_proxy_status = normalize_label(row.get("cnv_proxy_status", ""))

    if call == "malignant_hcc_high_conf":
        return "module3_high_conf_malignant"
    if call in {"malignant_hcc_cnv_support", "malignant_hcc_probable"}:
        return "module3_cnv_supported_malignant"
    if call in {"malignant_hcc_marker_proliferation_needs_cnv_review", "cnv_not_available"}:
        return "malignant_like_needs_review"
    if role == "end_malignant_review":
        return "malignant_like_needs_review"
    if copykat_status == "aneuploid" or copykat_pred == "aneuploid":
        return "copykat_aneuploid_without_module3_call"
    if cnv_proxy_status == "aneuploid_proxy":
        return "cnv_proxy_aneuploid_without_module3_call"
    return "no_cnv_evidence_or_reference"


def assign_pseudotime_bins(values: pd.Series, n_bins: int = 10) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    out = pd.Series(pd.NA, index=values.index, dtype="Int64")
    finite = numeric.replace([np.inf, -np.inf], np.nan).dropna()
    if finite.empty:
        return out
    n_effective = min(n_bins, int(finite.nunique()))
    if n_effective <= 1:
        out.loc[finite.index] = 0
        return out
    ranked = finite.rank(method="first")
    out.loc[finite.index] = pd.qcut(ranked, q=n_effective, labels=list(range(n_effective)), duplicates="drop").astype("Int64")
    return out


def mean_numeric(values: pd.Series) -> float:
    numeric = pd.to_numeric(values, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    if numeric.empty:
        return float("nan")
    return float(numeric.mean())


def series_or_nan(df: pd.DataFrame, column: str) -> pd.Series:
    if column in df.columns:
        return df[column]
    return pd.Series(np.nan, index=df.index)


def bin_cnv_overlay(df: pd.DataFrame, pseudotime_col: str = "pseudotime_norm", n_bins: int = 10) -> pd.DataFrame:
    if pseudotime_col not in df.columns:
        raise ValueError(f"Missing pseudotime column: {pseudotime_col}")

    work = df.copy()
    work[pseudotime_col] = pd.to_numeric(work[pseudotime_col], errors="coerce")
    if "cnv_evidence_tier" not in work.columns:
        work["cnv_evidence_tier"] = work.apply(assign_cnv_evidence_tier, axis=1)
    work["pseudotime_bin"] = assign_pseudotime_bins(work[pseudotime_col], n_bins=n_bins)
    work = work.dropna(subset=[pseudotime_col, "pseudotime_bin"])
    if work.empty:
        return pd.DataFrame()

    rows = []
    for bin_id, group in work.groupby("pseudotime_bin", observed=True, sort=True):
        tier = group["cnv_evidence_tier"].astype(str)
        copykat_aneuploid = (
            group["copykat_aneuploid"]
            if "copykat_aneuploid" in group.columns
            else group.get("copykat_status", pd.Series("", index=group.index)).astype(str).str.lower().eq("aneuploid")
        )
        cnv_proxy_aneuploid = (
            group["cnv_proxy_aneuploid"]
            if "cnv_proxy_aneuploid" in group.columns
            else group.get("cnv_proxy_status", pd.Series("", index=group.index)).astype(str).str.lower().eq("aneuploid_proxy")
        )
        row = {
            "pseudotime_bin": int(bin_id),
            "n_cells": int(group.shape[0]),
            "mean_pseudotime": mean_numeric(group[pseudotime_col]),
            "cnv_supported_fraction": numeric_fraction(tier.isin(CNV_SUPPORTED_TIERS)),
            "high_conf_fraction": numeric_fraction(tier.eq("module3_high_conf_malignant")),
            "cnv_support_fraction": numeric_fraction(tier.eq("module3_cnv_supported_malignant")),
            "review_fraction": numeric_fraction(tier.eq("malignant_like_needs_review")),
            "copykat_aneuploid_fraction": numeric_fraction(pd.Series(copykat_aneuploid, index=group.index)),
            "cnv_proxy_aneuploid_fraction": numeric_fraction(pd.Series(cnv_proxy_aneuploid, index=group.index)),
            "mean_cnv_proxy_burden": mean_numeric(series_or_nan(group, "cnv_proxy_burden")),
            "mean_cnv_proxy_high_bin_fraction": mean_numeric(series_or_nan(group, "cnv_proxy_high_bin_fraction")),
            "mean_cnv_proxy_max_abs_bin_log2": mean_numeric(series_or_nan(group, "cnv_proxy_max_abs_bin_log2")),
            "mean_hcc_malignant_score_z": mean_numeric(series_or_nan(group, "hcc_malignant_associated_score_z")),
            "mean_proliferation_score_z": mean_numeric(series_or_nan(group, "proliferation_score_z")),
            "mean_hcc_malignant_module": mean_numeric(series_or_nan(group, "HCC_Malignant_Associated")),
            "mean_proliferation_module": mean_numeric(series_or_nan(group, "Proliferation")),
            "mean_stressed_injured_module": mean_numeric(series_or_nan(group, "Stressed_Injured")),
            "mean_mature_hepatocyte_module": mean_numeric(series_or_nan(group, "Mature_Hepatocyte")),
        }
        rows.append(row)
    return pd.DataFrame(rows).sort_values("pseudotime_bin").reset_index(drop=True)


def finite_spearman(values: pd.Series, pseudotime: pd.Series) -> dict[str, object]:
    data = pd.DataFrame({"value": pd.to_numeric(values, errors="coerce"), "pseudotime": pd.to_numeric(pseudotime, errors="coerce")})
    data = data.replace([np.inf, -np.inf], np.nan).dropna()
    if data.shape[0] < 3 or data["value"].nunique() < 2 or data["pseudotime"].nunique() < 2:
        return {"n_cells": int(data.shape[0]), "spearman_rho": np.nan, "spearman_pvalue": np.nan}
    rho, pvalue = spearmanr(data["value"], data["pseudotime"])
    return {"n_cells": int(data.shape[0]), "spearman_rho": float(rho), "spearman_pvalue": float(pvalue)}


def minmax_scale(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    finite = numeric.replace([np.inf, -np.inf], np.nan).dropna()
    if finite.empty:
        return numeric * np.nan
    low = finite.min()
    high = finite.max()
    if high == low:
        return numeric * 0.0
    return (numeric - low) / (high - low)


def read_evidence_cells(path: Path) -> pd.DataFrame:
    desired = [
        "cell_id",
        "sample_id",
        "study_sample",
        "cnv_sample",
        "sample_source_class",
        "cnv_proxy_burden",
        "cnv_proxy_z",
        "cnv_proxy_high_bin_fraction",
        "cnv_proxy_max_abs_bin_log2",
        "cnv_proxy_status",
        "hcc_malignant_associated_score_z",
        "proliferation_score_z",
        "regenerative_progenitor_score_z",
        "malignant_hcc_evidence",
        "copykat_pred",
        "copykat_status",
        "malignant_hcc_call",
        "malignant_hcc_cnv_method",
        "malignant_hcc_evidence_copykat",
        "trajectory_role",
        "trajectory_include_main",
        "trajectory_include_cnv_strict",
        "sample_disease_stage",
        "cell_disease_stage",
        "trajectory_root_end_role",
        "trajectory_root_cell_selected",
        "trajectory_end_cell_selected",
        "global_umap_1",
        "global_umap_2",
    ]
    available = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
    usecols = [column for column in desired if column in available]
    missing = sorted(set(desired) - set(usecols))
    evidence = pd.read_csv(path, sep="\t", usecols=usecols)
    if "cell_id" not in evidence.columns:
        raise ValueError(f"Missing cell_id in {path}")
    evidence["cell_id"] = evidence["cell_id"].astype(str)
    evidence["cnv_evidence_tier"] = evidence.apply(assign_cnv_evidence_tier, axis=1)
    evidence["module3_cnv_supported"] = evidence["cnv_evidence_tier"].isin(CNV_SUPPORTED_TIERS)
    evidence["malignant_like_review"] = evidence["cnv_evidence_tier"].eq("malignant_like_needs_review")
    evidence["copykat_aneuploid"] = evidence.get("copykat_status", pd.Series("", index=evidence.index)).astype(str).str.lower().eq("aneuploid")
    evidence["cnv_proxy_aneuploid"] = evidence.get("cnv_proxy_status", pd.Series("", index=evidence.index)).astype(str).str.lower().eq("aneuploid_proxy")
    evidence.attrs["missing_columns"] = missing
    return evidence


def correlation_rows(df: pd.DataFrame, pseudotime_col: str, run_id: str, method: str) -> list[dict[str, object]]:
    features = {
        "module3_cnv_supported": "module3_cnv_supported",
        "malignant_like_review": "malignant_like_review",
        "copykat_aneuploid": "copykat_aneuploid",
        "cnv_proxy_aneuploid": "cnv_proxy_aneuploid",
        "cnv_proxy_burden": "cnv_proxy_burden",
        "cnv_proxy_z": "cnv_proxy_z",
        "cnv_proxy_high_bin_fraction": "cnv_proxy_high_bin_fraction",
        "cnv_proxy_max_abs_bin_log2": "cnv_proxy_max_abs_bin_log2",
        "hcc_malignant_associated_score_z": "hcc_malignant_associated_score_z",
        "proliferation_score_z": "proliferation_score_z",
        "regenerative_progenitor_score_z": "regenerative_progenitor_score_z",
        "HCC_Malignant_Associated": "HCC_Malignant_Associated",
        "Proliferation": "Proliferation",
        "Stressed_Injured": "Stressed_Injured",
        "Mature_Hepatocyte": "Mature_Hepatocyte",
    }
    rows = []
    for feature, column in features.items():
        if column not in df.columns:
            continue
        values = df[column].map(boolish) if df[column].dtype == bool else df[column]
        stats = finite_spearman(values, df[pseudotime_col])
        rows.append({"run_id": run_id, "method": method, "feature": feature, **stats})
    return rows


def evidence_tier_counts(df: pd.DataFrame, run_id: str, method: str) -> pd.DataFrame:
    counts = df["cnv_evidence_tier"].value_counts(dropna=False).rename_axis("cnv_evidence_tier").reset_index(name="n_cells")
    total = counts["n_cells"].sum()
    counts["fraction"] = counts["n_cells"] / total if total else np.nan
    counts.insert(0, "method", method)
    counts.insert(0, "run_id", run_id)
    return counts


def plot_cnv_overlay(bin_summary: pd.DataFrame, run_id: str, method: str, figures_dir: Path) -> str | None:
    if bin_summary.empty:
        return None
    figures_dir.mkdir(parents=True, exist_ok=True)
    sub = bin_summary.sort_values("pseudotime_bin")
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    fraction_specs = [
        ("cnv_supported_fraction", "CNV-supported malignant fraction"),
        ("copykat_aneuploid_fraction", "CopyKAT aneuploid fraction"),
        ("review_fraction", "malignant-like review fraction"),
        ("cnv_proxy_aneuploid_fraction", "CNV proxy aneuploid fraction"),
    ]
    for column, label in fraction_specs:
        if column in sub.columns:
            ax.plot(sub["mean_pseudotime"], sub[column], marker="o", linewidth=1.8, markersize=3, label=label)

    scaled_specs = [
        ("mean_hcc_malignant_module", "HCC malignant module, scaled"),
        ("mean_proliferation_module", "Proliferation module, scaled"),
        ("mean_cnv_proxy_burden", "CNV proxy burden, scaled"),
    ]
    for column, label in scaled_specs:
        if column in sub.columns and pd.to_numeric(sub[column], errors="coerce").notna().any():
            ax.plot(sub["mean_pseudotime"], minmax_scale(sub[column]), linestyle="--", linewidth=1.4, label=label)

    ax.set_xlabel("Normalized pseudotime")
    ax.set_ylabel("Fraction or min-max scaled mean")
    ax.set_ylim(-0.03, 1.03)
    ax.set_title(f"CNV and malignant evidence overlay: {run_id} / {method}")
    ax.legend(fontsize=7, ncol=2, frameon=False)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    path = figures_dir / f"trajectory_module5_5_cnv_overlay__{run_id}__{method}.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return str(path.resolve())


def main() -> int:
    args = parse_args()
    start = time.time()
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    args.figures_dir.mkdir(parents=True, exist_ok=True)

    evidence = read_evidence_cells(args.evidence_cells)
    module_scores = pd.read_csv(args.module_scores, sep="\t")
    module_scores["cell_id"] = module_scores["cell_id"].astype(str)

    overlay_rows = []
    bin_summary_rows = []
    tier_count_rows = []
    corr_rows = []
    figure_paths = []

    for run in pseudotime_inputs(args.metadata_dir):
        run_id = str(run["run_id"])
        pt = pd.read_csv(run["path"], sep="\t")
        pt["cell_id"] = pt["cell_id"].astype(str)
        scores = module_scores.loc[module_scores["run_id"].eq(run_id)].drop(columns=["run_id"], errors="ignore")

        for method, pt_col in run["methods"].items():
            print(f"OVERLAY {run_id} {method}", flush=True)
            method_pt = pt[["cell_id", pt_col]].rename(columns={pt_col: "pseudotime_norm"})
            method_pt = method_pt.dropna(subset=["pseudotime_norm"])
            merged = method_pt.merge(evidence, on="cell_id", how="left").merge(scores, on="cell_id", how="left")
            merged.insert(0, "method", method)
            merged.insert(0, "run_id", run_id)

            overlay_rows.append(merged)
            tier_count_rows.append(evidence_tier_counts(merged, run_id, method))
            corr_rows.extend(correlation_rows(merged, "pseudotime_norm", run_id, method))

            bins = bin_cnv_overlay(merged, pseudotime_col="pseudotime_norm", n_bins=args.n_bins)
            if not bins.empty:
                bins.insert(0, "method", method)
                bins.insert(0, "run_id", run_id)
                bin_summary_rows.append(bins)
                figure_path = plot_cnv_overlay(bins, run_id, method, args.figures_dir)
                if figure_path:
                    figure_paths.append(figure_path)

    overlay = pd.concat(overlay_rows, ignore_index=True)
    bin_summary = pd.concat(bin_summary_rows, ignore_index=True) if bin_summary_rows else pd.DataFrame()
    tier_counts = pd.concat(tier_count_rows, ignore_index=True) if tier_count_rows else pd.DataFrame()
    correlations = pd.DataFrame(corr_rows)
    if not correlations.empty:
        correlations["abs_spearman_rho"] = correlations["spearman_rho"].abs()
        correlations = correlations.sort_values(["run_id", "method", "abs_spearman_rho"], ascending=[True, True, False])

    overlay_path = args.metadata_dir / "trajectory_module5_5_cnv_malignant_overlay_by_cell.tsv.gz"
    bin_summary_path = args.metadata_dir / "trajectory_module5_5_pseudotime_bin_summary.tsv"
    tier_counts_path = args.metadata_dir / "trajectory_module5_5_evidence_tier_counts.tsv"
    correlations_path = args.metadata_dir / "trajectory_module5_5_evidence_correlations.tsv"
    report_path = args.metadata_dir / "trajectory_module5_5_report.json"

    overlay.to_csv(overlay_path, sep="\t", index=False, compression="gzip")
    bin_summary.to_csv(bin_summary_path, sep="\t", index=False)
    tier_counts.to_csv(tier_counts_path, sep="\t", index=False)
    correlations.to_csv(correlations_path, sep="\t", index=False)

    report = {
        "module": "5.5",
        "method": "overlay CNV, CopyKAT, CNV proxy, and malignant marker evidence along trajectory pseudotime",
        "n_bins": int(args.n_bins),
        "n_overlay_rows": int(overlay.shape[0]),
        "n_bin_rows": int(bin_summary.shape[0]),
        "n_correlation_rows": int(correlations.shape[0]),
        "missing_evidence_columns": evidence.attrs.get("missing_columns", []),
        "outputs": {
            "overlay_by_cell": str(overlay_path.resolve()),
            "pseudotime_bin_summary": str(bin_summary_path.resolve()),
            "evidence_tier_counts": str(tier_counts_path.resolve()),
            "evidence_correlations": str(correlations_path.resolve()),
            "figures": figure_paths,
        },
        "elapsed_seconds": round(time.time() - start, 3),
    }
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
