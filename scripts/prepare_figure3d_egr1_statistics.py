#!/usr/bin/env python3
"""Calculate Figure 3D EGR1 perturbation-score statistics along pseudotime."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.neighbors import NearestNeighbors

try:
    from figure3_egr1_common import (
        PROJECT_ROOT,
        PSEUDOTIME_COLUMN,
        SEED,
        TARGET_TF,
        assign_fixed_stage,
        json_safe,
        write_json,
    )
except ModuleNotFoundError:
    from scripts.figure3_egr1_common import (
        PROJECT_ROOT,
        PSEUDOTIME_COLUMN,
        SEED,
        TARGET_TF,
        assign_fixed_stage,
        json_safe,
        write_json,
    )


DEFAULT_C_DIR = PROJECT_ROOT / "metadata/driver/figure3c_egr1"
DEFAULT_OUT_DIR = PROJECT_ROOT / "metadata/driver/figure3d_egr1"


def map_pseudotime_to_grid(
    cells: pd.DataFrame,
    grid: pd.DataFrame,
    x_col: str,
    y_col: str,
    k: int,
) -> pd.DataFrame:
    coordinates = cells[[x_col, y_col]].to_numpy(dtype=float)
    query = grid[["grid_x", "grid_y"]].to_numpy(dtype=float)
    k_eff = min(k, len(cells))
    nn = NearestNeighbors(n_neighbors=k_eff)
    nn.fit(coordinates)
    distances, indices = nn.kneighbors(query)
    pseudotime = pd.to_numeric(cells["pseudotime"], errors="coerce").to_numpy(dtype=float)
    out = grid.copy()
    means, sds = [], []
    for row_distances, row_indices in zip(distances, indices, strict=True):
        bandwidth = max(float(np.median(row_distances)), 1e-8)
        weights = np.exp(-(row_distances**2) / (2 * bandwidth**2))
        weights /= weights.sum()
        values = pseudotime[row_indices]
        mean = float(np.sum(values * weights))
        means.append(mean)
        sds.append(float(np.sqrt(np.sum(weights * (values - mean) ** 2))))
    out["pseudotime_grid"] = means
    out["pseudotime_sd"] = sds
    out["valid"] = (
        out["keep_score"].astype(str).str.lower().isin({"true", "1"})
        & np.isfinite(pd.to_numeric(out["inner_product_score_grid"], errors="coerce"))
        & np.isfinite(out["pseudotime_grid"])
    )
    out["fixed_stage"] = assign_fixed_stage(out["pseudotime_grid"])
    return out


def summarize_bins(frame: pd.DataFrame, space: str, n_bins: int) -> pd.DataFrame:
    data = frame.loc[frame["valid"]].copy()
    data["pseudotime_bin"] = pd.cut(
        data["pseudotime_grid"],
        bins=np.linspace(0, 1, n_bins + 1),
        labels=False,
        include_lowest=True,
    )
    rows = []
    for bin_index, group in data.groupby("pseudotime_bin", observed=True):
        scores = group["inner_product_score_grid"].to_numpy(dtype=float)
        rows.append(
            {
                "space": space,
                "pseudotime_bin": int(bin_index) + 1,
                "pseudotime_start": float(bin_index) / n_bins,
                "pseudotime_end": float(bin_index + 1) / n_bins,
                "pseudotime_center": float(group["pseudotime_grid"].mean()),
                "n_grid_points": int(len(group)),
                "score_mean": float(np.mean(scores)),
                "score_median": float(np.median(scores)),
                "score_q25": float(np.quantile(scores, 0.25)),
                "score_q75": float(np.quantile(scores, 0.75)),
                "positive_fraction": float(np.mean(scores > 0)),
                "absolute_score_mean": float(np.mean(np.abs(scores))),
                "absolute_score_median": float(np.median(np.abs(scores))),
            }
        )
    return pd.DataFrame(rows)


def bootstrap_spearman(frame: pd.DataFrame, rng: np.random.Generator, n_bootstrap: int) -> dict:
    data = frame.loc[frame["valid"], ["pseudotime_grid", "inner_product_score_grid"]].dropna()
    x = data["pseudotime_grid"].to_numpy(dtype=float)
    y = data["inner_product_score_grid"].to_numpy(dtype=float)
    rho = float(stats.spearmanr(x, y).statistic)
    boot = np.empty(n_bootstrap, dtype=float)
    for index in range(n_bootstrap):
        selected = rng.integers(0, len(data), size=len(data))
        boot[index] = stats.spearmanr(x[selected], y[selected]).statistic
    return {
        "rho": rho,
        "ci95_low": float(np.nanquantile(boot, 0.025)),
        "ci95_high": float(np.nanquantile(boot, 0.975)),
        "n_grid_points": int(len(data)),
        "n_bootstrap": int(n_bootstrap),
    }


def fixed_stage_rows(frame: pd.DataFrame, space: str) -> tuple[list[dict], list[dict]]:
    data = frame.loc[frame["valid"] & frame["fixed_stage"].notna()].copy()
    summary_rows = []
    for stage in ["early", "intermediate", "late"]:
        scores = data.loc[data["fixed_stage"].eq(stage), "inner_product_score_grid"].to_numpy(dtype=float)
        summary_rows.append(
            {
                "row_type": "stage_summary",
                "space": space,
                "comparison": stage,
                "metric": "score",
                "test": "",
                "n_1": int(len(scores)),
                "n_2": np.nan,
                "score_mean": float(np.mean(scores)),
                "score_median": float(np.median(scores)),
                "score_q25": float(np.quantile(scores, 0.25)),
                "score_q75": float(np.quantile(scores, 0.75)),
                "positive_fraction": float(np.mean(scores > 0)),
                "absolute_score_mean": float(np.mean(np.abs(scores))),
                "absolute_score_median": float(np.median(np.abs(scores))),
                "statistic": np.nan,
                "p_value": np.nan,
                "p_adjust": np.nan,
                "effect_median_difference": np.nan,
            }
        )

    test_rows = []
    stage_values = {
        stage: data.loc[data["fixed_stage"].eq(stage), "inner_product_score_grid"].to_numpy(dtype=float)
        for stage in ["early", "intermediate", "late"]
    }
    kruskal = stats.kruskal(stage_values["early"], stage_values["intermediate"], stage_values["late"])
    test_rows.append(
        {
            "row_type": "test",
            "space": space,
            "comparison": "early_vs_intermediate_vs_late",
            "metric": "signed_score",
            "test": "Kruskal-Wallis",
            "n_1": int(len(data)),
            "n_2": np.nan,
            "score_mean": np.nan,
            "score_median": np.nan,
            "score_q25": np.nan,
            "score_q75": np.nan,
            "positive_fraction": np.nan,
            "absolute_score_mean": np.nan,
            "absolute_score_median": np.nan,
            "statistic": float(kruskal.statistic),
            "p_value": float(kruskal.pvalue),
            "effect_median_difference": np.nan,
        }
    )
    for other in ["early", "late"]:
        for metric, transform in [("signed_score", lambda x: x), ("absolute_score", np.abs)]:
            left = transform(stage_values["intermediate"])
            right = transform(stage_values[other])
            test = stats.mannwhitneyu(left, right, alternative="two-sided")
            test_rows.append(
                {
                    "row_type": "test",
                    "space": space,
                    "comparison": f"intermediate_vs_{other}",
                    "metric": metric,
                    "test": "Mann-Whitney U",
                    "n_1": int(len(left)),
                    "n_2": int(len(right)),
                    "score_mean": np.nan,
                    "score_median": np.nan,
                    "score_q25": np.nan,
                    "score_q75": np.nan,
                    "positive_fraction": np.nan,
                    "absolute_score_mean": np.nan,
                    "absolute_score_median": np.nan,
                    "statistic": float(test.statistic),
                    "p_value": float(test.pvalue),
                    "effect_median_difference": float(np.median(left) - np.median(right)),
                }
            )
    p_values = np.array([row["p_value"] for row in test_rows], dtype=float)
    order = np.argsort(p_values)
    adjusted = np.empty_like(p_values)
    ranked = p_values[order] * len(p_values) / np.arange(1, len(p_values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1.0)
    adjusted[order] = ranked
    for row, value in zip(test_rows, adjusted, strict=True):
        row["p_adjust"] = float(value)
    return summary_rows, test_rows


def bootstrap_peak_stability(
    frame: pd.DataFrame,
    space: str,
    rng: np.random.Generator,
    n_bootstrap: int,
    n_bins: int,
) -> tuple[pd.DataFrame, dict]:
    data = frame.loc[frame["valid"] & frame["fixed_stage"].notna()].copy()
    data["pseudotime_bin"] = pd.cut(
        data["pseudotime_grid"],
        bins=np.linspace(0, 1, n_bins + 1),
        labels=False,
        include_lowest=True,
    )
    stage_peak = []
    bin_peak = []
    for _ in range(n_bootstrap):
        boot = data.iloc[rng.integers(0, len(data), size=len(data))]
        stage_abs = boot.groupby("fixed_stage", observed=True)["inner_product_score_grid"].apply(
            lambda values: float(np.median(np.abs(values)))
        )
        if len(stage_abs):
            stage_peak.append(str(stage_abs.idxmax()))
        bin_abs = boot.groupby("pseudotime_bin", observed=True)["inner_product_score_grid"].apply(
            lambda values: float(np.median(np.abs(values)))
        )
        if len(bin_abs):
            bin_peak.append(int(bin_abs.idxmax()) + 1)
    rows = []
    for scope, values in [("fixed_stage", stage_peak), ("pseudotime_bin", bin_peak)]:
        counts = pd.Series(values).value_counts(dropna=False)
        for label, count in counts.items():
            rows.append(
                {
                    "space": space,
                    "peak_scope": scope,
                    "peak_label": label,
                    "bootstrap_count": int(count),
                    "bootstrap_fraction": float(count / n_bootstrap),
                    "n_bootstrap": int(n_bootstrap),
                }
            )
    summary = {
        "intermediate_peak_fraction": float(np.mean(np.asarray(stage_peak) == "intermediate")),
        "modal_peak_stage": str(pd.Series(stage_peak).mode().iloc[0]),
        "modal_peak_bin": int(pd.Series(bin_peak).mode().iloc[0]),
    }
    return pd.DataFrame(rows), summary


def fit_two_change_points(frame: pd.DataFrame) -> dict:
    data = frame.loc[frame["valid"], ["pseudotime_grid", "inner_product_score_grid"]].sort_values("pseudotime_grid")
    x = data["pseudotime_grid"].to_numpy(dtype=float)
    y = data["inner_product_score_grid"].to_numpy(dtype=float)
    candidates = np.unique(np.quantile(x, np.linspace(0.20, 0.80, 25)))
    best = None
    for first in candidates:
        for second in candidates:
            if second <= first + 0.10:
                continue
            masks = [x < first, (x >= first) & (x < second), x >= second]
            if min(mask.sum() for mask in masks) < 10:
                continue
            sse = sum(float(np.sum((y[mask] - np.mean(y[mask])) ** 2)) for mask in masks)
            if best is None or sse < best["sse"]:
                best = {"change_point_1": float(first), "change_point_2": float(second), "sse": sse}
    return best or {"change_point_1": None, "change_point_2": None, "sse": None}


def observed_pattern(stage_table: pd.DataFrame) -> dict:
    summary = stage_table.loc[
        stage_table["row_type"].eq("stage_summary") & stage_table["space"].eq("CellOracle UMAP grid")
    ].copy()
    signed_peak = summary.loc[summary["score_median"].abs().idxmax(), "comparison"]
    absolute_peak = summary.loc[summary["absolute_score_median"].idxmax(), "comparison"]
    medians = dict(zip(summary["comparison"], summary["score_median"], strict=True))
    if absolute_peak == "intermediate":
        conclusion = (
            "Virtual EGR1 knockout showed its largest absolute grid-level perturbation score in the "
            "prespecified intermediate stage; the sign and pairwise statistics determine whether it opposed progression."
        )
    elif np.all(np.asarray(list(medians.values())) > 0):
        conclusion = "EGR1 knockout scores were positive across all three fixed stages, with the strongest absolute stage reported directly."
    elif np.all(np.asarray(list(medians.values())) < 0):
        conclusion = "EGR1 knockout scores were negative across all three fixed stages, with the strongest absolute stage reported directly."
    else:
        conclusion = "EGR1 knockout showed a stage-dependent mixed-sign perturbation-score pattern without imposing an intermediate-stage peak."
    return {
        "signed_absolute_median_peak_stage": str(signed_peak),
        "absolute_effect_peak_stage": str(absolute_peak),
        "stage_score_medians": {str(key): float(value) for key, value in medians.items()},
        "results_recommendation": conclusion,
    }


def run(c_dir: Path, out_dir: Path, k_neighbors: int, n_bins: int, n_bootstrap: int, seed: int) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    cells = pd.read_csv(c_dir / "figure3c_egr1_cells_with_tsne_projection.tsv.gz", sep="\t")
    if not cells["tf"].astype(str).eq(TARGET_TF).all():
        raise ValueError("Figure 3D input contains non-EGR1 rows")
    spaces = {
        "CellOracle UMAP grid": (
            pd.read_csv(c_dir / "figure3c_egr1_inner_product_grid_umap.tsv.gz", sep="\t"),
            "umap_1",
            "umap_2",
            "figure3d_egr1_pseudotime_inner_product_umap.tsv.gz",
        ),
        "Expression t-SNE grid": (
            pd.read_csv(c_dir / "figure3c_egr1_inner_product_grid_tsne.tsv.gz", sep="\t"),
            "tsne_1",
            "tsne_2",
            "figure3d_egr1_pseudotime_inner_product_tsne.tsv.gz",
        ),
    }
    rng = np.random.default_rng(seed)
    mapped, bin_frames, stage_rows, test_rows, peak_frames = {}, [], [], [], []
    correlations, peaks, change_points = {}, {}, {}
    for space, (grid, x_col, y_col, filename) in spaces.items():
        result = map_pseudotime_to_grid(cells, grid, x_col, y_col, k_neighbors)
        result["space"] = space
        result.to_csv(out_dir / filename, sep="\t", index=False, compression="gzip")
        mapped[space] = result
        bin_frames.append(summarize_bins(result, space, n_bins))
        summaries, tests = fixed_stage_rows(result, space)
        stage_rows.extend(summaries)
        test_rows.extend(tests)
        correlation_rng = np.random.default_rng(rng.integers(0, np.iinfo(np.int32).max))
        peak_rng = np.random.default_rng(rng.integers(0, np.iinfo(np.int32).max))
        correlations[space] = bootstrap_spearman(result, correlation_rng, n_bootstrap)
        peak_frame, peak_summary = bootstrap_peak_stability(
            result, space, peak_rng, n_bootstrap, n_bins
        )
        peak_frames.append(peak_frame)
        peaks[space] = peak_summary
        change_points[space] = fit_two_change_points(result)

    bins = pd.concat(bin_frames, ignore_index=True)
    bins_path = out_dir / "figure3d_egr1_pseudotime_bin_summary.tsv"
    bins.to_csv(bins_path, sep="\t", index=False)
    stage_table = pd.DataFrame(stage_rows + test_rows)
    stage_path = out_dir / "figure3d_egr1_stage_comparison.tsv"
    stage_table.to_csv(stage_path, sep="\t", index=False)
    peak_table = pd.concat(peak_frames, ignore_index=True)
    peak_path = out_dir / "figure3d_egr1_bootstrap_peak_summary.tsv"
    peak_table.to_csv(peak_path, sep="\t", index=False)
    pattern = observed_pattern(stage_table)
    report = {
        "module": "Figure 3D",
        "target_tf": TARGET_TF,
        "score_definition": "EGR1 perturbation vector dot baseline developmental vector",
        "pseudotime": PSEUDOTIME_COLUMN,
        "n_cells": int(len(cells)),
        "n_bins": int(n_bins),
        "k_neighbors": int(k_neighbors),
        "seed": int(seed),
        "fixed_stage_definition": {
            "early": "[0.00, 0.33)",
            "intermediate": "[0.33, 0.67)",
            "late": "[0.67, 1.00]",
        },
        "spearman_bootstrap": correlations,
        "bootstrap_peak_stability": peaks,
        "data_driven_change_point_sensitivity": change_points,
        "observed_umap_pattern": pattern,
        "outputs": {
            "umap_grid": str((out_dir / spaces["CellOracle UMAP grid"][3]).resolve()),
            "tsne_grid": str((out_dir / spaces["Expression t-SNE grid"][3]).resolve()),
            "bin_summary": str(bins_path.resolve()),
            "stage_comparison": str(stage_path.resolve()),
            "bootstrap_peak_summary": str(peak_path.resolve()),
        },
        "caveat": "The fixed early/intermediate/late analysis is primary. Change points, LOESS, and t-SNE are sensitivity analyses.",
    }
    report_path = out_dir / "figure3d_egr1_report.json"
    write_json(json_safe(report), report_path)
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figure3c-dir", type=Path, default=DEFAULT_C_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--k-neighbors", type=int, default=50)
    parser.add_argument("--n-bins", type=int, default=10)
    parser.add_argument("--n-bootstrap", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=SEED)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = run(
        args.figure3c_dir,
        args.out_dir,
        args.k_neighbors,
        args.n_bins,
        args.n_bootstrap,
        args.seed,
    )
    print(
        json.dumps(
            {
                "target_tf": report["target_tf"],
                "observed_umap_pattern": report["observed_umap_pattern"],
                "report": str((args.out_dir / "figure3d_egr1_report.json").resolve()),
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

