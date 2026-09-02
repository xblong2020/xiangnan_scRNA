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

try:
    from score_celloracle_perturbation_module6_9b import aggregate_tf_scores
except ModuleNotFoundError:
    from scripts.score_celloracle_perturbation_module6_9b import aggregate_tf_scores


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_H5AD = PROJECT_ROOT / "data/processed/driver/celloracle_module6_6/celloracle_module6_6_input.h5ad"
DEFAULT_METADATA_DIR = PROJECT_ROOT / "metadata/driver"
FATE_COL = "cellrank_fate_prob_cnv_supported_malignant"
PHASE_COL = "driver_main_strict__pseudotime_phase"


def attach_cell_metadata(cell_scores: pd.DataFrame, h5ad_path: Path) -> pd.DataFrame:
    adata = ad.read_h5ad(h5ad_path)
    cols = [
        "dataset",
        "sample_id",
        "celloracle_main_strict",
        "driver_primary_eligible",
        "trajectory_include_main",
        "trajectory_include_cnv_strict",
        PHASE_COL,
    ]
    missing = [col for col in cols if col not in adata.obs.columns]
    if missing:
        raise ValueError(f"Missing h5ad obs columns for Module 6.10: {missing}")
    metadata = adata.obs[cols].copy()
    metadata.insert(0, "cell_id", adata.obs_names.astype(str))
    metadata["driver_union"] = True
    metadata["phase_3level"] = metadata[PHASE_COL].astype(str).replace({"middle": "intermediate"})
    merged = cell_scores.merge(metadata, on="cell_id", how="left", validate="many_to_one", suffixes=("", "_meta"))
    if merged["dataset"].isna().any():
        raise ValueError("Some cell-level score rows could not be aligned to h5ad metadata")
    return merged


def aggregate_subset_scores(
    subset_scores: pd.DataFrame,
    subset_name: str,
    fate_high_threshold: float,
) -> pd.DataFrame:
    summary = aggregate_tf_scores(subset_scores, fate_high_threshold=fate_high_threshold)
    summary.insert(0, "subset", subset_name)
    summary.insert(1, "n_subset_cells", int(subset_scores["cell_id"].nunique()) if "cell_id" in subset_scores.columns else int(len(subset_scores)))
    return summary


def compute_core_subset_scores(cell_scores: pd.DataFrame, fate_high_threshold: float) -> pd.DataFrame:
    subsets = []
    subsets.append(aggregate_subset_scores(cell_scores, "driver_union_all_cells", fate_high_threshold))
    main = cell_scores.loc[cell_scores["celloracle_main_strict"].astype(bool)]
    subsets.append(aggregate_subset_scores(main, "main_strict_cells", fate_high_threshold))
    for phase in ["early", "intermediate", "late"]:
        phase_df = cell_scores.loc[cell_scores["phase_3level"].astype(str) == phase]
        if len(phase_df):
            subsets.append(aggregate_subset_scores(phase_df, f"phase_{phase}", fate_high_threshold))
    return pd.concat(subsets, axis=0, ignore_index=True)


def compute_leave_one_group_out(
    cell_scores: pd.DataFrame,
    group_col: str,
    fate_high_threshold: float,
    min_remaining_cells: int,
) -> pd.DataFrame:
    frames = []
    for group in sorted(cell_scores[group_col].astype(str).dropna().unique()):
        kept = cell_scores.loc[cell_scores[group_col].astype(str) != group]
        if kept["cell_id"].nunique() < min_remaining_cells:
            continue
        summary = aggregate_subset_scores(kept, f"leave_one_{group_col}_out", fate_high_threshold)
        summary.insert(2, "left_out_col", group_col)
        summary.insert(3, "left_out_group", group)
        summary.insert(4, "n_left_out_cells", int(cell_scores.loc[cell_scores[group_col].astype(str) == group, "cell_id"].nunique()))
        frames.append(summary)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, axis=0, ignore_index=True)


def summarize_rank_stability(subset_scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for tf, df in subset_scores.groupby("tf", sort=False):
        rows.append(
            {
                "tf": tf,
                "n_subsets": int(df["subset"].nunique()),
                "median_rank": float(df["quantitative_rank"].median()),
                "max_rank": int(df["quantitative_rank"].max()),
                "min_rank": int(df["quantitative_rank"].min()),
                "rank_iqr": float(df["quantitative_rank"].quantile(0.75) - df["quantitative_rank"].quantile(0.25)),
                "top3_fraction": float((df["quantitative_rank"] <= 3).mean()),
                "top5_fraction": float((df["quantitative_rank"] <= 5).mean()),
                "mean_quantitative_score": float(df["quantitative_perturbation_score"].mean()),
                "min_quantitative_score": float(df["quantitative_perturbation_score"].min()),
            }
        )
    return pd.DataFrame(rows).sort_values(["top5_fraction", "median_rank"], ascending=[False, True]).reset_index(drop=True)


def classify_control_tfs(tf_scores: pd.DataFrame, n_low: int = 3) -> pd.DataFrame:
    rows = []
    ranked = tf_scores.sort_values("quantitative_perturbation_score", ascending=True)
    for tf in ranked.head(n_low)["tf"]:
        rows.append({"tf": tf, "control_type": "low_score_tf", "control_rationale": "Bottom quantitative perturbation score among simulated TFs."})
    non_malignant = tf_scores.loc[pd.to_numeric(tf_scores["malignant_fate_direction_score"], errors="coerce") <= 0]
    for tf in non_malignant.sort_values("malignant_fate_direction_score").head(n_low)["tf"]:
        rows.append({"tf": tf, "control_type": "non_malignant_direction_tf", "control_rationale": "Perturbation does not shift high-CNV/malignant-like cells away from malignant direction."})
    if "state_specificity_ratio" in tf_scores.columns:
        score_median = pd.to_numeric(tf_scores["quantitative_perturbation_score"], errors="coerce").median()
        proxy_pool = tf_scores.loc[pd.to_numeric(tf_scores["quantitative_perturbation_score"], errors="coerce") <= score_median]
        if proxy_pool.empty:
            proxy_pool = tf_scores
        proxy = proxy_pool.sort_values(["state_specificity_ratio", "quantitative_perturbation_score"], ascending=[True, True]).head(1)
        for tf in proxy["tf"]:
            rows.append(
                {
                    "tf": tf,
                    "control_type": "housekeeping_like_proxy",
                    "control_rationale": "Proxy only: no canonical housekeeping TF was included in the 15 simulated perturbation TFs; selected as low state-specificity/low-impact broad-effect comparator.",
                }
            )
    controls = pd.DataFrame(rows).drop_duplicates(["tf", "control_type"]).reset_index(drop=True)
    return controls


def compute_dataset_sample_dominance(cell_scores: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    rows = []
    for tf, df in cell_scores.groupby("tf", sort=False):
        for group_col in group_cols:
            counts = df.drop_duplicates("cell_id")[group_col].astype(str).value_counts(dropna=False)
            if counts.empty:
                continue
            rows.append(
                {
                    "tf": tf,
                    "group_col": group_col,
                    "dominant_group": str(counts.index[0]),
                    "dominant_group_cells": int(counts.iloc[0]),
                    "total_cells": int(counts.sum()),
                    "dominant_fraction": float(counts.iloc[0] / counts.sum()),
                    "n_groups": int(len(counts)),
                }
            )
    return pd.DataFrame(rows)


def summarize_lodo_drop(full_scores: pd.DataFrame, lodo_scores: pd.DataFrame) -> pd.DataFrame:
    if lodo_scores.empty:
        return pd.DataFrame()
    base = full_scores.loc[full_scores["subset"] == "driver_union_all_cells", ["tf", "quantitative_perturbation_score", "quantitative_rank"]].rename(
        columns={
            "quantitative_perturbation_score": "full_quantitative_score",
            "quantitative_rank": "full_rank",
        }
    )
    merged = lodo_scores.merge(base, on="tf", how="left")
    merged["score_delta_vs_full"] = merged["quantitative_perturbation_score"] - merged["full_quantitative_score"]
    merged["rank_delta_vs_full"] = merged["quantitative_rank"] - merged["full_rank"]
    rows = []
    for tf, df in merged.groupby("tf", sort=False):
        worst = df.loc[df["score_delta_vs_full"].idxmin()]
        rows.append(
            {
                "tf": tf,
                "min_lodo_score": float(df["quantitative_perturbation_score"].min()),
                "max_lodo_rank": int(df["quantitative_rank"].max()),
                "top5_lodo_fraction": float((df["quantitative_rank"] <= 5).mean()),
                "worst_left_out_col": worst["left_out_col"],
                "worst_left_out_group": worst["left_out_group"],
                "worst_score_delta_vs_full": float(worst["score_delta_vs_full"]),
            }
        )
    return pd.DataFrame(rows).sort_values(["top5_lodo_fraction", "min_lodo_score"], ascending=[False, False]).reset_index(drop=True)


def flag_review_risks(
    dominance: pd.DataFrame,
    lodo_summary: pd.DataFrame,
    controls: pd.DataFrame,
    dominance_threshold: float,
) -> pd.DataFrame:
    rows = []
    for _, row in dominance.iterrows():
        if row["dominant_fraction"] >= dominance_threshold:
            rows.append(
                {
                    "risk_type": "single_group_dominance",
                    "tf": row["tf"],
                    "detail": f"{row['group_col']}={row['dominant_group']} contributes {row['dominant_fraction']:.2%} of cells.",
                    "severity": "review_attention",
                }
            )
    if not lodo_summary.empty:
        for _, row in lodo_summary.iterrows():
            if row["top5_lodo_fraction"] < 0.75:
                rows.append(
                    {
                        "risk_type": "lodo_rank_instability",
                        "tf": row["tf"],
                        "detail": f"Top5 fraction under leave-one-out is {row['top5_lodo_fraction']:.2f}; worst exclusion {row['worst_left_out_col']}={row['worst_left_out_group']}.",
                        "severity": "review_attention",
                    }
                )
    if "housekeeping_like_proxy" in set(controls["control_type"]):
        rows.append(
            {
                "risk_type": "control_limitation",
                "tf": "NA",
                "detail": "Canonical housekeeping TF was not included in the 15 simulated CellOracle TFs; a proxy control is reported.",
                "severity": "method_caveat",
            }
        )
    return pd.DataFrame(rows)


def build_main_union_comparison(core_subset_scores: pd.DataFrame) -> pd.DataFrame:
    union = core_subset_scores.loc[
        core_subset_scores["subset"] == "driver_union_all_cells",
        ["tf", "quantitative_rank", "quantitative_perturbation_score"],
    ].rename(
        columns={
            "quantitative_rank": "driver_union_rank",
            "quantitative_perturbation_score": "driver_union_score",
        }
    )
    main = core_subset_scores.loc[
        core_subset_scores["subset"] == "main_strict_cells",
        ["tf", "quantitative_rank", "quantitative_perturbation_score"],
    ].rename(
        columns={
            "quantitative_rank": "main_strict_rank",
            "quantitative_perturbation_score": "main_strict_score",
        }
    )
    comparison = union.merge(main, on="tf", how="outer")
    comparison["rank_delta_main_minus_union"] = comparison["main_strict_rank"] - comparison["driver_union_rank"]
    comparison["score_delta_main_minus_union"] = comparison["main_strict_score"] - comparison["driver_union_score"]
    return comparison.sort_values(["driver_union_rank", "main_strict_rank"]).reset_index(drop=True)


def build_phase_wide_summary(core_subset_scores: pd.DataFrame) -> pd.DataFrame:
    phases = ["phase_early", "phase_intermediate", "phase_late"]
    frames = []
    for phase in phases:
        sub = core_subset_scores.loc[
            core_subset_scores["subset"] == phase,
            ["tf", "quantitative_rank", "quantitative_perturbation_score"],
        ].rename(
            columns={
                "quantitative_rank": f"{phase}_rank",
                "quantitative_perturbation_score": f"{phase}_score",
            }
        )
        frames.append(sub)
    if not frames:
        return pd.DataFrame()
    phase_summary = frames[0]
    for frame in frames[1:]:
        phase_summary = phase_summary.merge(frame, on="tf", how="outer")
    return phase_summary.sort_values(["phase_late_rank", "phase_intermediate_rank", "phase_early_rank"]).reset_index(drop=True)


def run_module6_10(
    h5ad_path: Path,
    metadata_dir: Path,
    fate_high_threshold: float,
    min_remaining_cells: int,
    dominance_threshold: float,
) -> dict:
    metadata_dir.mkdir(parents=True, exist_ok=True)
    cell_scores = pd.read_csv(metadata_dir / "celloracle_module6_9b_cell_level_scores.tsv.gz", sep="\t")
    cell_scores = attach_cell_metadata(cell_scores, h5ad_path)
    full_tf_scores = pd.read_csv(metadata_dir / "celloracle_module6_9b_quantitative_tf_scores.tsv", sep="\t")

    core_subset_scores = compute_core_subset_scores(cell_scores, fate_high_threshold)
    lodo_dataset = compute_leave_one_group_out(cell_scores, "dataset", fate_high_threshold, min_remaining_cells)
    lodo_sample = compute_leave_one_group_out(cell_scores, "sample_id", fate_high_threshold, min_remaining_cells)
    lodo_scores = pd.concat([df for df in [lodo_dataset, lodo_sample] if len(df)], axis=0, ignore_index=True)
    stability = summarize_rank_stability(core_subset_scores)
    lodo_summary = summarize_lodo_drop(core_subset_scores, lodo_scores)
    controls = classify_control_tfs(full_tf_scores)
    dominance = compute_dataset_sample_dominance(cell_scores, ["dataset", "sample_id"])
    risk_flags = flag_review_risks(dominance, lodo_summary, controls, dominance_threshold)
    main_union_comparison = build_main_union_comparison(core_subset_scores)
    phase_wide_summary = build_phase_wide_summary(core_subset_scores)

    outputs = {
        "cell_scores_with_metadata": metadata_dir / "celloracle_module6_10_cell_scores_with_metadata.tsv.gz",
        "subset_scores": metadata_dir / "celloracle_module6_10_subset_scores.tsv",
        "main_union_comparison": metadata_dir / "celloracle_module6_10_main_strict_vs_driver_union.tsv",
        "phase_wide_summary": metadata_dir / "celloracle_module6_10_phase_wide_summary.tsv",
        "leave_one_group_out_scores": metadata_dir / "celloracle_module6_10_leave_one_group_out_scores.tsv.gz",
        "rank_stability": metadata_dir / "celloracle_module6_10_rank_stability.tsv",
        "lodo_summary": metadata_dir / "celloracle_module6_10_lodo_summary.tsv",
        "negative_controls": metadata_dir / "celloracle_module6_10_negative_controls.tsv",
        "dataset_sample_dominance": metadata_dir / "celloracle_module6_10_dataset_sample_dominance.tsv",
        "review_risk_flags": metadata_dir / "celloracle_module6_10_review_risk_flags.tsv",
    }
    cell_scores.to_csv(outputs["cell_scores_with_metadata"], sep="\t", index=False)
    core_subset_scores.to_csv(outputs["subset_scores"], sep="\t", index=False)
    main_union_comparison.to_csv(outputs["main_union_comparison"], sep="\t", index=False)
    phase_wide_summary.to_csv(outputs["phase_wide_summary"], sep="\t", index=False)
    lodo_scores.to_csv(outputs["leave_one_group_out_scores"], sep="\t", index=False)
    stability.to_csv(outputs["rank_stability"], sep="\t", index=False)
    lodo_summary.to_csv(outputs["lodo_summary"], sep="\t", index=False)
    controls.to_csv(outputs["negative_controls"], sep="\t", index=False)
    dominance.to_csv(outputs["dataset_sample_dominance"], sep="\t", index=False)
    risk_flags.to_csv(outputs["review_risk_flags"], sep="\t", index=False)

    top_main = core_subset_scores.loc[core_subset_scores["subset"] == "main_strict_cells"].head(5)["tf"].astype(str).tolist()
    top_union = core_subset_scores.loc[core_subset_scores["subset"] == "driver_union_all_cells"].head(5)["tf"].astype(str).tolist()
    return {
        "outputs": {k: str(v) for k, v in outputs.items()},
        "n_cell_score_rows": int(len(cell_scores)),
        "n_tfs": int(full_tf_scores["tf"].nunique()),
        "n_datasets": int(cell_scores["dataset"].astype(str).nunique()),
        "n_samples": int(cell_scores["sample_id"].astype(str).nunique()),
        "top5_driver_union": top_union,
        "top5_main_strict": top_main,
        "n_risk_flags": int(len(risk_flags)),
        "negative_controls": controls.to_dict(orient="records"),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 6.10 CellOracle perturbation robustness and negative controls")
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA_DIR)
    parser.add_argument("--fate-high-threshold", type=float, default=0.5)
    parser.add_argument("--min-remaining-cells", type=int, default=1000)
    parser.add_argument("--dominance-threshold", type=float, default=0.5)
    parser.add_argument("--report", type=Path, default=DEFAULT_METADATA_DIR / "celloracle_module6_10_report.json")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = datetime.now(timezone.utc)
    result = run_module6_10(
        h5ad_path=args.h5ad,
        metadata_dir=args.metadata_dir,
        fate_high_threshold=args.fate_high_threshold,
        min_remaining_cells=args.min_remaining_cells,
        dominance_threshold=args.dominance_threshold,
    )
    finished = datetime.now(timezone.utc)
    report = {
        "module": "6.10",
        "method": "CellOracle perturbation robustness and negative controls",
        "created_at_utc": finished.isoformat(),
        "elapsed_seconds": (finished - started).total_seconds(),
        "parameters": {
            "h5ad": str(args.h5ad),
            "metadata_dir": str(args.metadata_dir),
            "fate_high_threshold": args.fate_high_threshold,
            "min_remaining_cells": args.min_remaining_cells,
            "dominance_threshold": args.dominance_threshold,
        },
        "result": result,
        "python_runtime": {
            "executable": sys.executable,
            "version": sys.version,
            "platform": platform.platform(),
        },
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(
        {
            "n_tfs": result["n_tfs"],
            "n_datasets": result["n_datasets"],
            "n_samples": result["n_samples"],
            "top5_driver_union": result["top5_driver_union"],
            "top5_main_strict": result["top5_main_strict"],
            "n_risk_flags": result["n_risk_flags"],
            "report": str(args.report),
        },
        ensure_ascii=False,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
