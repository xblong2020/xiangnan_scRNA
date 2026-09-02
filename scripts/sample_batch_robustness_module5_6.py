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

TRUE_VALUES = {"1", "true", "t", "yes", "y", "aneuploid"}
FALSE_VALUES = {"0", "false", "f", "no", "n", "diploid", "unknown", "not.defined", "not_run", "nan", ""}

GROUP_COLUMNS = ["sample_id", "study_sample", "cnv_sample", "dataset", "source_h5ad", "_scvi_batch"]
LOO_GROUP_COLUMNS = ["sample_id", "cnv_sample", "dataset", "_scvi_batch"]
FEATURES = [
    "module3_cnv_supported",
    "copykat_aneuploid",
    "cnv_proxy_aneuploid",
    "malignant_like_review",
    "HCC_Malignant_Associated",
    "Proliferation",
    "Stressed_Injured",
    "Mature_Hepatocyte",
    "cnv_proxy_burden",
    "cnv_proxy_high_bin_fraction",
    "cnv_proxy_max_abs_bin_log2",
    "hcc_malignant_associated_score_z",
    "proliferation_score_z",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 5.6: sample and batch robustness for trajectory evidence trends.")
    parser.add_argument(
        "--overlay",
        type=Path,
        default=ROOT / "metadata/trajectory/trajectory_module5_5_cnv_malignant_overlay_by_cell.tsv.gz",
    )
    parser.add_argument(
        "--stage-cells",
        type=Path,
        default=ROOT / "metadata/trajectory/trajectory_module5_2_stage_root_end_cells.tsv.gz",
    )
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/trajectory")
    parser.add_argument("--figures-dir", type=Path, default=ROOT / "figures/trajectory")
    parser.add_argument("--n-bins", type=int, default=10)
    parser.add_argument("--min-group-cells", type=int, default=30)
    parser.add_argument("--min-remaining-cells", type=int, default=100)
    return parser.parse_args()


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


def numeric_feature(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values.astype(float)
    numeric = pd.to_numeric(values, errors="coerce")
    if numeric.notna().sum() > 0:
        return numeric.astype(float)
    return values.map(boolish).astype(float)


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


def finite_spearman(values: pd.Series, pseudotime: pd.Series, min_cells: int = 3) -> dict[str, object]:
    data = pd.DataFrame({"value": numeric_feature(values), "pseudotime": pd.to_numeric(pseudotime, errors="coerce")})
    data = data.replace([np.inf, -np.inf], np.nan).dropna()
    if data.shape[0] < min_cells or data["value"].nunique() < 2 or data["pseudotime"].nunique() < 2:
        return {"n_cells": int(data.shape[0]), "spearman_rho": np.nan, "spearman_pvalue": np.nan}
    rho, pvalue = spearmanr(data["value"], data["pseudotime"])
    return {"n_cells": int(data.shape[0]), "spearman_rho": float(rho), "spearman_pvalue": float(pvalue)}


def trend_delta(df: pd.DataFrame, pseudotime_col: str, feature_col: str, n_bins: int = 10) -> dict[str, object]:
    if pseudotime_col not in df.columns or feature_col not in df.columns:
        return {
            "n_cells": 0,
            "n_bins_observed": 0,
            "early_mean": np.nan,
            "late_mean": np.nan,
            "late_minus_early_delta": np.nan,
        }
    data = pd.DataFrame(
        {
            "pseudotime": pd.to_numeric(df[pseudotime_col], errors="coerce"),
            "value": numeric_feature(df[feature_col]),
        }
    ).replace([np.inf, -np.inf], np.nan)
    data = data.dropna()
    if data.empty:
        return {
            "n_cells": 0,
            "n_bins_observed": 0,
            "early_mean": np.nan,
            "late_mean": np.nan,
            "late_minus_early_delta": np.nan,
        }
    data["bin"] = assign_pseudotime_bins(data["pseudotime"], n_bins=n_bins)
    data = data.dropna(subset=["bin"])
    if data.empty:
        return {
            "n_cells": 0,
            "n_bins_observed": 0,
            "early_mean": np.nan,
            "late_mean": np.nan,
            "late_minus_early_delta": np.nan,
        }
    means = data.groupby("bin", observed=True)["value"].mean().sort_index()
    early = float(means.iloc[0])
    late = float(means.iloc[-1])
    return {
        "n_cells": int(data.shape[0]),
        "n_bins_observed": int(means.shape[0]),
        "early_mean": early,
        "late_mean": late,
        "late_minus_early_delta": float(late - early),
    }


def summarize_group_trends(
    df: pd.DataFrame,
    group_col: str,
    pseudotime_col: str,
    feature_col: str,
    n_bins: int = 10,
    min_cells: int = 30,
) -> pd.DataFrame:
    if group_col not in df.columns:
        return pd.DataFrame()
    rows = []
    for group, group_df in df.dropna(subset=[group_col]).groupby(group_col, observed=True, sort=True):
        stats = trend_delta(group_df, pseudotime_col, feature_col, n_bins=n_bins)
        if stats["n_cells"] < min_cells:
            continue
        corr = finite_spearman(group_df[feature_col], group_df[pseudotime_col])
        rows.append(
            {
                "group": str(group),
                "n_cells": int(stats["n_cells"]),
                "pseudotime_min": float(pd.to_numeric(group_df[pseudotime_col], errors="coerce").min()),
                "pseudotime_max": float(pd.to_numeric(group_df[pseudotime_col], errors="coerce").max()),
                "pseudotime_mean": float(pd.to_numeric(group_df[pseudotime_col], errors="coerce").mean()),
                "n_bins_observed": int(stats["n_bins_observed"]),
                "early_mean": stats["early_mean"],
                "late_mean": stats["late_mean"],
                "late_minus_early_delta": stats["late_minus_early_delta"],
                "spearman_rho": corr["spearman_rho"],
                "spearman_pvalue": corr["spearman_pvalue"],
            }
        )
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(["late_minus_early_delta", "n_cells"], ascending=[False, False]).reset_index(drop=True)


def leave_one_group_out_delta(
    df: pd.DataFrame,
    group_col: str,
    pseudotime_col: str,
    feature_col: str,
    n_bins: int = 10,
    min_remaining_cells: int = 100,
) -> pd.DataFrame:
    if group_col not in df.columns:
        return pd.DataFrame()
    base = df.dropna(subset=[group_col]).copy()
    overall = trend_delta(base, pseudotime_col, feature_col, n_bins=n_bins)
    rows = []
    for group, group_df in base.groupby(group_col, observed=True, sort=True):
        keep = base[group_col].ne(group)
        remaining = base.loc[keep]
        if remaining.shape[0] < min_remaining_cells:
            continue
        loo = trend_delta(remaining, pseudotime_col, feature_col, n_bins=n_bins)
        rows.append(
            {
                "omitted_group": str(group),
                "omitted_n_cells": int(group_df.shape[0]),
                "remaining_n_cells": int(remaining.shape[0]),
                "overall_delta": overall["late_minus_early_delta"],
                "loo_delta": loo["late_minus_early_delta"],
                "delta_shift_from_overall": (
                    float(loo["late_minus_early_delta"] - overall["late_minus_early_delta"])
                    if pd.notna(loo["late_minus_early_delta"]) and pd.notna(overall["late_minus_early_delta"])
                    else np.nan
                ),
            }
        )
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values("delta_shift_from_overall").reset_index(drop=True)


def batch_centered_spearman(
    df: pd.DataFrame,
    batch_col: str,
    pseudotime_col: str,
    feature_col: str,
    min_cells: int = 30,
) -> dict[str, object]:
    if batch_col not in df.columns:
        return {"n_cells": 0, "n_batches": 0, "spearman_rho": np.nan, "spearman_pvalue": np.nan}
    data = pd.DataFrame(
        {
            "batch": df[batch_col],
            "pseudotime": pd.to_numeric(df[pseudotime_col], errors="coerce"),
            "value": numeric_feature(df[feature_col]),
        }
    ).replace([np.inf, -np.inf], np.nan)
    data = data.dropna()
    batch_sizes = data["batch"].value_counts()
    keep_batches = batch_sizes[batch_sizes >= 2].index
    data = data.loc[data["batch"].isin(keep_batches)].copy()
    if data.shape[0] < min_cells or data["batch"].nunique() < 1:
        return {"n_cells": int(data.shape[0]), "n_batches": int(data["batch"].nunique()), "spearman_rho": np.nan, "spearman_pvalue": np.nan}
    data["pseudotime_centered"] = data["pseudotime"] - data.groupby("batch", observed=True)["pseudotime"].transform("mean")
    data["value_centered"] = data["value"] - data.groupby("batch", observed=True)["value"].transform("mean")
    if data["value_centered"].nunique() < 2 or data["pseudotime_centered"].nunique() < 2:
        return {"n_cells": int(data.shape[0]), "n_batches": int(data["batch"].nunique()), "spearman_rho": np.nan, "spearman_pvalue": np.nan}
    rho, pvalue = spearmanr(data["value_centered"], data["pseudotime_centered"])
    return {
        "n_cells": int(data.shape[0]),
        "n_batches": int(data["batch"].nunique()),
        "spearman_rho": float(rho),
        "spearman_pvalue": float(pvalue),
    }


def largest_group_fraction(df: pd.DataFrame, group_col: str) -> dict[str, object]:
    if group_col not in df.columns:
        return {"n_cells": 0, "n_groups": 0, "largest_group": "", "largest_group_n_cells": 0, "largest_group_fraction": np.nan}
    counts = df[group_col].dropna().astype(str).value_counts()
    if counts.empty:
        return {"n_cells": 0, "n_groups": 0, "largest_group": "", "largest_group_n_cells": 0, "largest_group_fraction": np.nan}
    total = int(counts.sum())
    return {
        "n_cells": total,
        "n_groups": int(counts.shape[0]),
        "largest_group": str(counts.index[0]),
        "largest_group_n_cells": int(counts.iloc[0]),
        "largest_group_fraction": float(counts.iloc[0] / total),
    }


def robustness_label(overall_delta: float, positive_group_fraction: float, min_loo_delta: float) -> str:
    if pd.isna(overall_delta):
        return "insufficient"
    if overall_delta > 0 and positive_group_fraction >= 0.6 and (pd.isna(min_loo_delta) or min_loo_delta > 0):
        return "robust_positive"
    if overall_delta > 0 and positive_group_fraction >= 0.6:
        return "positive_with_group_sensitivity"
    if overall_delta > 0:
        return "overall_positive_group_mixed"
    if overall_delta < 0:
        return "overall_negative"
    return "flat"


def read_overlay_with_batch_metadata(overlay_path: Path, stage_cells_path: Path) -> pd.DataFrame:
    overlay = pd.read_csv(overlay_path, sep="\t")
    overlay["cell_id"] = overlay["cell_id"].astype(str)
    stage_columns = ["cell_id", "dataset", "source_h5ad", "_scvi_batch"]
    available = pd.read_csv(stage_cells_path, sep="\t", nrows=0).columns.tolist()
    usecols = [column for column in stage_columns if column in available]
    stage = pd.read_csv(stage_cells_path, sep="\t", usecols=usecols)
    stage["cell_id"] = stage["cell_id"].astype(str)
    add_cols = [column for column in usecols if column == "cell_id" or column not in overlay.columns]
    if len(add_cols) > 1:
        overlay = overlay.merge(stage[add_cols].drop_duplicates("cell_id"), on="cell_id", how="left")
    return overlay


def composition_rows(df: pd.DataFrame, group_col: str, run_id: str, method: str) -> list[dict[str, object]]:
    if group_col not in df.columns:
        return []
    counts = df[group_col].dropna().astype(str).value_counts().rename_axis("group").reset_index(name="n_cells")
    total = int(counts["n_cells"].sum())
    dominance = largest_group_fraction(df, group_col)
    rows = []
    for _, row in counts.iterrows():
        rows.append(
            {
                "run_id": run_id,
                "method": method,
                "group_type": group_col,
                "group": row["group"],
                "n_cells": int(row["n_cells"]),
                "fraction": float(row["n_cells"] / total) if total else np.nan,
                "n_groups": dominance["n_groups"],
                "largest_group": dominance["largest_group"],
                "largest_group_fraction": dominance["largest_group_fraction"],
            }
        )
    return rows


def plot_sample_robustness(group_trends: pd.DataFrame, run_id: str, method: str, figures_dir: Path) -> str | None:
    sub = group_trends.loc[
        group_trends["run_id"].eq(run_id)
        & group_trends["method"].eq(method)
        & group_trends["group_type"].eq("sample_id")
        & group_trends["feature"].isin(["module3_cnv_supported", "HCC_Malignant_Associated", "Proliferation"])
    ].copy()
    if sub.empty:
        return None
    figures_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, 5.2))
    for feature, feature_df in sub.groupby("feature", observed=True):
        ordered = feature_df.sort_values("pseudotime_mean")
        ax.scatter(ordered["pseudotime_mean"], ordered["late_minus_early_delta"], s=np.clip(ordered["n_cells"], 20, 220), alpha=0.7, label=feature)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Sample mean normalized pseudotime")
    ax.set_ylabel("Late minus early pseudotime delta")
    ax.set_title(f"Sample-level robustness: {run_id} / {method}")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    path = figures_dir / f"trajectory_module5_6_sample_robustness__{run_id}__{method}.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return str(path.resolve())


def main() -> int:
    args = parse_args()
    start = time.time()
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    args.figures_dir.mkdir(parents=True, exist_ok=True)

    overlay = read_overlay_with_batch_metadata(args.overlay, args.stage_cells)
    group_columns = [column for column in GROUP_COLUMNS if column in overlay.columns]
    loo_group_columns = [column for column in LOO_GROUP_COLUMNS if column in overlay.columns]
    features = [feature for feature in FEATURES if feature in overlay.columns]

    group_trend_rows = []
    loo_rows = []
    adjusted_rows = []
    composition = []
    overall_delta_lookup: dict[tuple[str, str, str], float] = {}
    dominance_lookup: dict[tuple[str, str, str], dict[str, object]] = {}

    for (run_id, method), sub in overlay.groupby(["run_id", "method"], observed=True, sort=True):
        print(f"ROBUSTNESS {run_id} {method}", flush=True)
        for group_col in group_columns:
            composition.extend(composition_rows(sub, group_col, str(run_id), str(method)))
            dominance_lookup[(str(run_id), str(method), group_col)] = largest_group_fraction(sub, group_col)

        for feature in features:
            overall = trend_delta(sub, "pseudotime_norm", feature, n_bins=args.n_bins)
            overall_delta_lookup[(str(run_id), str(method), feature)] = overall["late_minus_early_delta"]
            raw_corr = finite_spearman(sub[feature], sub["pseudotime_norm"])
            for group_col in group_columns:
                trends = summarize_group_trends(sub, group_col, "pseudotime_norm", feature, n_bins=args.n_bins, min_cells=args.min_group_cells)
                if not trends.empty:
                    trends.insert(0, "feature", feature)
                    trends.insert(0, "group_type", group_col)
                    trends.insert(0, "method", method)
                    trends.insert(0, "run_id", run_id)
                    group_trend_rows.append(trends)

                adjusted = batch_centered_spearman(sub, group_col, "pseudotime_norm", feature, min_cells=args.min_group_cells)
                dominance = largest_group_fraction(sub, group_col)
                adjusted_rows.append(
                    {
                        "run_id": run_id,
                        "method": method,
                        "batch_type": group_col,
                        "feature": feature,
                        "raw_n_cells": raw_corr["n_cells"],
                        "raw_spearman_rho": raw_corr["spearman_rho"],
                        "raw_spearman_pvalue": raw_corr["spearman_pvalue"],
                        "centered_n_cells": adjusted["n_cells"],
                        "centered_n_batches": adjusted["n_batches"],
                        "centered_spearman_rho": adjusted["spearman_rho"],
                        "centered_spearman_pvalue": adjusted["spearman_pvalue"],
                        "largest_group_fraction": dominance["largest_group_fraction"],
                    }
                )

            for group_col in loo_group_columns:
                loo = leave_one_group_out_delta(
                    sub,
                    group_col,
                    "pseudotime_norm",
                    feature,
                    n_bins=args.n_bins,
                    min_remaining_cells=args.min_remaining_cells,
                )
                if not loo.empty:
                    loo.insert(0, "feature", feature)
                    loo.insert(0, "group_type", group_col)
                    loo.insert(0, "method", method)
                    loo.insert(0, "run_id", run_id)
                    loo_rows.append(loo)

    group_trends = pd.concat(group_trend_rows, ignore_index=True) if group_trend_rows else pd.DataFrame()
    leave_one_out = pd.concat(loo_rows, ignore_index=True) if loo_rows else pd.DataFrame()
    adjusted = pd.DataFrame(adjusted_rows)
    composition_df = pd.DataFrame(composition)

    robustness_rows = []
    for (run_id, method, feature), overall_delta in overall_delta_lookup.items():
        for group_col in group_columns:
            if not group_trends.empty:
                feature_groups = group_trends.loc[
                    group_trends["run_id"].eq(run_id)
                    & group_trends["method"].eq(method)
                    & group_trends["group_type"].eq(group_col)
                    & group_trends["feature"].eq(feature)
                ]
            else:
                feature_groups = pd.DataFrame()
            positive_fraction = float(feature_groups["late_minus_early_delta"].gt(0).mean()) if not feature_groups.empty else np.nan
            median_group_delta = float(feature_groups["late_minus_early_delta"].median()) if not feature_groups.empty else np.nan

            if group_col in loo_group_columns and not leave_one_out.empty:
                loo_sub = leave_one_out.loc[
                    leave_one_out["run_id"].eq(run_id)
                    & leave_one_out["method"].eq(method)
                    & leave_one_out["group_type"].eq(group_col)
                    & leave_one_out["feature"].eq(feature)
                ]
            else:
                loo_sub = pd.DataFrame()
            min_loo_delta = float(loo_sub["loo_delta"].min()) if not loo_sub.empty else np.nan
            max_abs_loo_shift = float(loo_sub["delta_shift_from_overall"].abs().max()) if not loo_sub.empty else np.nan
            dominance = dominance_lookup.get((run_id, method, group_col), {})
            robustness_rows.append(
                {
                    "run_id": run_id,
                    "method": method,
                    "group_type": group_col,
                    "feature": feature,
                    "overall_delta": overall_delta,
                    "n_groups_tested": int(feature_groups.shape[0]),
                    "positive_group_fraction": positive_fraction,
                    "median_group_delta": median_group_delta,
                    "min_loo_delta": min_loo_delta,
                    "max_abs_loo_shift": max_abs_loo_shift,
                    "largest_group_fraction": dominance.get("largest_group_fraction", np.nan),
                    "robustness_label": robustness_label(overall_delta, positive_fraction, min_loo_delta),
                }
            )
    robustness = pd.DataFrame(robustness_rows)

    figure_paths = []
    if not group_trends.empty:
        for (run_id, method), _ in group_trends.groupby(["run_id", "method"], observed=True, sort=True):
            path = plot_sample_robustness(group_trends, str(run_id), str(method), args.figures_dir)
            if path:
                figure_paths.append(path)

    group_trends_path = args.metadata_dir / "trajectory_module5_6_group_trends.tsv"
    leave_one_out_path = args.metadata_dir / "trajectory_module5_6_leave_one_group_out.tsv"
    adjusted_path = args.metadata_dir / "trajectory_module5_6_batch_adjusted_correlations.tsv"
    composition_path = args.metadata_dir / "trajectory_module5_6_group_composition.tsv"
    robustness_path = args.metadata_dir / "trajectory_module5_6_robustness_summary.tsv"
    report_path = args.metadata_dir / "trajectory_module5_6_report.json"

    group_trends.to_csv(group_trends_path, sep="\t", index=False)
    leave_one_out.to_csv(leave_one_out_path, sep="\t", index=False)
    adjusted.to_csv(adjusted_path, sep="\t", index=False)
    composition_df.to_csv(composition_path, sep="\t", index=False)
    robustness.to_csv(robustness_path, sep="\t", index=False)

    report = {
        "module": "5.6",
        "method": "sample and batch robustness of trajectory-associated CNV and malignant evidence",
        "n_bins": int(args.n_bins),
        "min_group_cells": int(args.min_group_cells),
        "min_remaining_cells": int(args.min_remaining_cells),
        "group_columns": group_columns,
        "loo_group_columns": loo_group_columns,
        "features": features,
        "n_group_trend_rows": int(group_trends.shape[0]),
        "n_leave_one_out_rows": int(leave_one_out.shape[0]),
        "n_batch_adjusted_rows": int(adjusted.shape[0]),
        "outputs": {
            "group_trends": str(group_trends_path.resolve()),
            "leave_one_group_out": str(leave_one_out_path.resolve()),
            "batch_adjusted_correlations": str(adjusted_path.resolve()),
            "group_composition": str(composition_path.resolve()),
            "robustness_summary": str(robustness_path.resolve()),
            "figures": figure_paths,
        },
        "elapsed_seconds": round(time.time() - start, 3),
    }
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
