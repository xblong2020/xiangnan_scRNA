from __future__ import annotations

import argparse
import json
import math
import platform
import time
from collections.abc import Iterable, Sequence
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA_DIR = ROOT / "metadata/driver"
DEFAULT_FIGURE_DIR = ROOT / "figures/driver"
DEFAULT_DRIVER_H5AD = ROOT / "data/processed/driver/driver_hepatocyte_trajectory.module6_1.h5ad"
DEFAULT_RUN_METHOD_LONG = DEFAULT_METADATA_DIR / "driver_module6_1_run_method_long.tsv.gz"
DEFAULT_CELLRANK_FATE = DEFAULT_METADATA_DIR / "driver_module6_2_cellrank_fate_probabilities.tsv.gz"
DEFAULT_PYSCENIC_AUC = DEFAULT_METADATA_DIR / "driver_module6_3_pyscenic_regulon_auc.tsv.gz"
DEFAULT_CISTARGET_AUC = DEFAULT_METADATA_DIR / "driver_module6_3c_cistarget_regulon_auc.tsv.gz"
DEFAULT_TF_TARGETS = DEFAULT_METADATA_DIR / "module8_tf_target_signature_genes.tsv"
DEFAULT_SOX4_STATE_SPECIFIC = DEFAULT_METADATA_DIR / "sctenifoldknk_module7_3_malignant_like_state_specific_genes.tsv"
DEFAULT_PROCESSED_DIR = ROOT / "data/processed"

TF_GENES = ["HNF4A", "PPARA", "JUN", "FOS", "JUND", "ATF3", "CEBPB", "EGR1", "SOX4"]
AP1_TFS = ["JUN", "FOS", "JUND", "ATF3"]
CEBPB_EGR1_TFS = ["CEBPB", "EGR1"]
RETENTION_TFS = ["HNF4A", "PPARA"]
SOX4_TFS = ["SOX4"]

ORDER_FEATURES = [
    "A_hnf4a_ppara_loss",
    "B_transition_activation",
    "C_sox4_axis",
    "C_malignant_like_fate",
]
ORDER_COMPARISONS = [
    ("A_loss_before_B_transition", "A_hnf4a_ppara_loss", "B_transition_activation", "<"),
    ("B_transition_before_C_sox4", "B_transition_activation", "C_sox4_axis", "<"),
    ("C_sox4_before_or_equal_malignant_fate", "C_sox4_axis", "C_malignant_like_fate", "<="),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 9.1 temporal hierarchy evidence analysis.")
    parser.add_argument("--driver-h5ad", type=Path, default=DEFAULT_DRIVER_H5AD)
    parser.add_argument("--run-method-long", type=Path, default=DEFAULT_RUN_METHOD_LONG)
    parser.add_argument("--cellrank-fate", type=Path, default=DEFAULT_CELLRANK_FATE)
    parser.add_argument("--pyscenic-auc", type=Path, default=DEFAULT_PYSCENIC_AUC)
    parser.add_argument("--cistarget-auc", type=Path, default=DEFAULT_CISTARGET_AUC)
    parser.add_argument("--tf-targets", type=Path, default=DEFAULT_TF_TARGETS)
    parser.add_argument("--sox4-state-specific", type=Path, default=DEFAULT_SOX4_STATE_SPECIFIC)
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--processed-dir", type=Path, default=DEFAULT_PROCESSED_DIR)
    parser.add_argument("--velocity-h5ad", type=Path, default=None)
    parser.add_argument("--top-n-sox4-targets", type=int, default=50)
    parser.add_argument("--n-bins", type=int, default=20)
    parser.add_argument("--n-bootstrap", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20260615)
    return parser.parse_args()


def package_version(name: str) -> str:
    try:
        return version(name)
    except PackageNotFoundError:
        return "not_installed"


def read_tsv_or_empty(path: Path, **kwargs) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", **kwargs)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def zscore_series(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce").astype(float)
    finite = numeric.replace([np.inf, -np.inf], np.nan)
    mean = finite.mean(skipna=True)
    std = finite.std(skipna=True, ddof=0)
    if not np.isfinite(std) or std == 0:
        return pd.Series(np.nan, index=values.index, dtype=float)
    return (finite - mean) / std


def mean_zscore(frame: pd.DataFrame, columns: Sequence[str]) -> pd.Series:
    available = [column for column in columns if column in frame.columns]
    if not available:
        return pd.Series(np.nan, index=frame.index, dtype=float)
    z = pd.DataFrame({column: zscore_series(frame[column]) for column in available}, index=frame.index)
    return z.mean(axis=1, skipna=True)


def tf_feature_columns(frame: pd.DataFrame, tfs: Sequence[str]) -> list[str]:
    prefixes = ["expr", "regulon", "regulon_pyscenic", "regulon_cistarget"]
    columns = []
    for prefix in prefixes:
        for tf in tfs:
            column = f"{prefix}_{tf}"
            if column in frame.columns:
                columns.append(column)
    return list(dict.fromkeys(columns))


def feature_availability_row(axis: str, requested: Sequence[str], frame: pd.DataFrame) -> dict[str, object]:
    available = [column for column in requested if column in frame.columns]
    return {
        "axis": axis,
        "n_requested_features": len(requested),
        "n_available_features": len(available),
        "requested_features": ";".join(requested),
        "available_features": ";".join(available),
        "missing_features": ";".join([column for column in requested if column not in frame.columns]),
    }


def build_axis_scores(feature_values: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    retention_features = tf_feature_columns(feature_values, RETENTION_TFS)
    ap1_features = tf_feature_columns(feature_values, AP1_TFS)
    cebpb_egr1_features = tf_feature_columns(feature_values, CEBPB_EGR1_TFS)
    sox4_features = tf_feature_columns(feature_values, SOX4_TFS)
    if "sox4_target_signature" in feature_values.columns:
        sox4_features.append("sox4_target_signature")
    malignant_features = [
        column
        for column in [
            "module_HCC_Malignant_Associated",
            "module_Proliferation",
            "cellrank_fate_prob_cnv_supported_malignant",
        ]
        if column in feature_values.columns
    ]

    scores = pd.DataFrame(index=feature_values.index)
    scores["A_hnf4a_ppara_retention"] = mean_zscore(feature_values, retention_features)
    scores["A_hnf4a_ppara_loss"] = -scores["A_hnf4a_ppara_retention"]
    scores["B_ap1_activation"] = mean_zscore(feature_values, ap1_features)
    scores["B_cebpb_egr1_activation"] = mean_zscore(feature_values, cebpb_egr1_features)
    scores["B_transition_activation"] = mean_zscore(scores, ["B_ap1_activation", "B_cebpb_egr1_activation"])
    scores["C_sox4_axis"] = mean_zscore(feature_values, sox4_features)
    scores["C_malignant_like_fate"] = mean_zscore(feature_values, malignant_features)

    availability = pd.DataFrame(
        [
            feature_availability_row("A_hnf4a_ppara_retention", retention_features, feature_values),
            feature_availability_row("A_hnf4a_ppara_loss", retention_features, feature_values),
            feature_availability_row("B_ap1_activation", ap1_features, feature_values),
            feature_availability_row("B_cebpb_egr1_activation", cebpb_egr1_features, feature_values),
            feature_availability_row(
                "B_transition_activation",
                ["B_ap1_activation", "B_cebpb_egr1_activation"],
                scores,
            ),
            feature_availability_row("C_sox4_axis", sox4_features, feature_values),
            feature_availability_row("C_malignant_like_fate", malignant_features, feature_values),
        ]
    )
    return scores, availability


def assign_pseudotime_bins(values: pd.Series, n_bins: int = 20) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce").replace([np.inf, -np.inf], np.nan)
    out = pd.Series(pd.NA, index=values.index, dtype="Int64")
    finite = numeric.dropna()
    if finite.empty:
        return out
    n_effective = min(n_bins, int(finite.nunique()), int(finite.shape[0]))
    if n_effective <= 1:
        out.loc[finite.index] = 0
        return out
    ranked = finite.rank(method="first")
    out.loc[finite.index] = pd.qcut(ranked, q=n_effective, labels=False, duplicates="drop").astype("Int64")
    return out


def smooth_values(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    if len(y) < 4:
        return y.astype(float)
    try:
        from statsmodels.nonparametric.smoothers_lowess import lowess

        smoothed = lowess(y, x, frac=min(0.6, max(0.25, 4.0 / len(y))), return_sorted=False)
        return np.asarray(smoothed, dtype=float)
    except Exception:
        return pd.Series(y).rolling(window=3, min_periods=1, center=True).mean().to_numpy(dtype=float)


def bin_temporal_values(
    values: pd.Series,
    pseudotime: pd.Series,
    n_bins: int,
    feature: str,
) -> pd.DataFrame:
    bins = assign_pseudotime_bins(pseudotime, n_bins=n_bins)
    work = pd.DataFrame({"pseudotime": pseudotime, "score": values, "bin": bins}).replace([np.inf, -np.inf], np.nan).dropna()
    if work.empty:
        return pd.DataFrame(
            columns=["feature", "bin", "n_cells", "mean_pseudotime", "mean_score", "sem_score", "smoothed_score"]
        )
    summary = (
        work.groupby("bin", observed=True)
        .agg(
            n_cells=("score", "size"),
            mean_pseudotime=("pseudotime", "mean"),
            mean_score=("score", "mean"),
            sem_score=("score", "sem"),
        )
        .reset_index()
        .sort_values("mean_pseudotime")
    )
    summary["feature"] = feature
    summary["smoothed_score"] = smooth_values(
        summary["mean_pseudotime"].to_numpy(dtype=float),
        summary["mean_score"].to_numpy(dtype=float),
    )
    return summary[["feature", "bin", "n_cells", "mean_pseudotime", "mean_score", "sem_score", "smoothed_score"]]


def detect_onset_from_bins(binned: pd.DataFrame) -> tuple[float, str]:
    if binned.empty or binned.shape[0] < 3:
        return np.nan, "insufficient_bins"
    scores = binned["smoothed_score"].to_numpy(dtype=float)
    times = binned["mean_pseudotime"].to_numpy(dtype=float)
    finite = np.isfinite(scores) & np.isfinite(times)
    if finite.sum() < 3:
        return np.nan, "insufficient_finite_bins"
    scores = scores[finite]
    times = times[finite]
    baseline_n = max(1, int(math.ceil(0.2 * len(scores))))
    baseline = scores[:baseline_n]
    baseline_mean = float(np.nanmean(baseline))
    baseline_sd = float(np.nanstd(baseline))
    fallback_sd = float(np.nanstd(scores))
    threshold_sd = baseline_sd if np.isfinite(baseline_sd) and baseline_sd > 0 else fallback_sd
    threshold = baseline_mean + 0.5 * threshold_sd if np.isfinite(threshold_sd) else baseline_mean
    for idx in range(baseline_n, len(scores)):
        window = scores[idx : min(idx + 3, len(scores))]
        if len(window) >= 2 and np.all(window > threshold):
            return float(times[idx]), "tested"
    return np.nan, "no_onset_detected"


def summarize_temporal_trend(
    values: pd.Series,
    pseudotime: pd.Series,
    n_bins: int = 20,
    feature: str = "",
) -> dict[str, object]:
    values = pd.to_numeric(values, errors="coerce").replace([np.inf, -np.inf], np.nan)
    pseudotime = pd.to_numeric(pseudotime, errors="coerce").replace([np.inf, -np.inf], np.nan)
    finite = values.notna() & pseudotime.notna()
    feature_name = feature or str(getattr(values, "name", "feature"))
    if finite.sum() < 10 or values.loc[finite].nunique() < 2 or pseudotime.loc[finite].nunique() < 2:
        return {
            "feature": feature_name,
            "n_cells": int(finite.sum()),
            "spearman_rho": np.nan,
            "spearman_pvalue": np.nan,
            "late_vs_early_delta": np.nan,
            "onset_time": np.nan,
            "peak_time": np.nan,
            "max_slope_time": np.nan,
            "trend_direction": "insufficient",
            "trend_status": "insufficient",
        }

    rho, pvalue = spearmanr(values.loc[finite], pseudotime.loc[finite])
    binned = bin_temporal_values(values.loc[finite], pseudotime.loc[finite], n_bins=n_bins, feature=feature_name)
    onset_time, onset_status = detect_onset_from_bins(binned)
    if binned.empty:
        peak_time = np.nan
        max_slope_time = np.nan
        delta = np.nan
    else:
        peak_idx = int(np.nanargmax(binned["smoothed_score"].to_numpy(dtype=float)))
        peak_time = float(binned["mean_pseudotime"].iloc[peak_idx])
        if binned.shape[0] >= 2:
            dy = np.diff(binned["smoothed_score"].to_numpy(dtype=float))
            dx = np.diff(binned["mean_pseudotime"].to_numpy(dtype=float))
            slopes = np.divide(dy, dx, out=np.full_like(dy, np.nan, dtype=float), where=dx != 0)
            slope_idx = int(np.nanargmax(slopes)) if np.isfinite(slopes).any() else 0
            max_slope_time = float(binned["mean_pseudotime"].iloc[slope_idx + 1])
        else:
            max_slope_time = np.nan
        edge_n = max(1, int(math.ceil(0.2 * binned.shape[0])))
        delta = float(binned["mean_score"].iloc[-edge_n:].mean() - binned["mean_score"].iloc[:edge_n].mean())
    if np.isfinite(rho) and np.isfinite(delta) and rho > 0.2 and delta > 0:
        direction = "increasing"
    elif np.isfinite(rho) and np.isfinite(delta) and rho < -0.2 and delta < 0:
        direction = "decreasing"
    else:
        direction = "flat_or_mixed"
    return {
        "feature": feature_name,
        "n_cells": int(finite.sum()),
        "spearman_rho": float(rho),
        "spearman_pvalue": float(pvalue),
        "late_vs_early_delta": delta,
        "onset_time": onset_time,
        "peak_time": peak_time,
        "max_slope_time": max_slope_time,
        "trend_direction": direction,
        "trend_status": onset_status,
    }


def grouped_bootstrap_order_tests(
    data: pd.DataFrame,
    pseudotime_col: str,
    group_col: str,
    n_bootstrap: int = 500,
    random_state: int = 20260615,
    n_bins: int = 20,
    features: Sequence[str] = ORDER_FEATURES,
) -> pd.DataFrame:
    required = [pseudotime_col, *features]
    missing = [column for column in required if column not in data.columns]
    if missing:
        raise KeyError(f"Missing columns for bootstrap order test: {missing}")
    if group_col not in data.columns:
        work = data.copy()
        work[group_col] = np.arange(work.shape[0]).astype(str)
    else:
        work = data.copy()
    work[group_col] = work[group_col].astype(str).fillna("unknown")
    groups = [group for group, sub in work.groupby(group_col, observed=True) if sub.shape[0] > 0]
    rng = np.random.default_rng(random_state)
    results = {name: [] for name, _, _, _ in ORDER_COMPARISONS}
    deltas = {name: [] for name, _, _, _ in ORDER_COMPARISONS}
    successes = 0

    for _ in range(n_bootstrap):
        if not groups:
            break
        sampled_groups = rng.choice(groups, size=len(groups), replace=True)
        sampled = pd.concat([work.loc[work[group_col].eq(group)] for group in sampled_groups], ignore_index=True)
        onsets: dict[str, float] = {}
        for feature in features:
            summary = summarize_temporal_trend(sampled[feature], sampled[pseudotime_col], n_bins=n_bins, feature=feature)
            onsets[feature] = float(summary["onset_time"]) if pd.notna(summary["onset_time"]) else np.nan
        if not all(np.isfinite(onsets.get(feature, np.nan)) for feature in features):
            continue
        successes += 1
        for name, upstream, downstream, operator in ORDER_COMPARISONS:
            upstream_time = onsets[upstream]
            downstream_time = onsets[downstream]
            passed = upstream_time <= downstream_time if operator == "<=" else upstream_time < downstream_time
            results[name].append(bool(passed))
            deltas[name].append(float(downstream_time - upstream_time))

    rows = []
    for name, upstream, downstream, operator in ORDER_COMPARISONS:
        values = results[name]
        rows.append(
            {
                "comparison": name,
                "upstream_feature": upstream,
                "downstream_feature": downstream,
                "operator": operator,
                "n_bootstrap_requested": int(n_bootstrap),
                "n_bootstrap_successful": int(len(values)),
                "order_probability": float(np.mean(values)) if values else np.nan,
                "median_time_delta": float(np.median(deltas[name])) if deltas[name] else np.nan,
                "support_label": "supported" if values and float(np.mean(values)) >= 0.7 else "not_supported",
            }
        )
    return pd.DataFrame(rows)


def benjamini_hochberg(pvalues: Sequence[float]) -> list[float]:
    p = np.asarray([1.0 if pd.isna(value) else float(value) for value in pvalues], dtype=float)
    n = len(p)
    order = np.argsort(p)
    adjusted = np.empty(n, dtype=float)
    running_min = 1.0
    for rank, idx in enumerate(order[::-1], start=1):
        original_rank = n - rank + 1
        value = min(running_min, p[idx] * n / original_rank)
        running_min = value
        adjusted[idx] = value
    return np.clip(adjusted, 0, 1).tolist()


def cohen_d(group_a: pd.Series, group_b: pd.Series) -> float:
    a = pd.to_numeric(group_a, errors="coerce").dropna().astype(float)
    b = pd.to_numeric(group_b, errors="coerce").dropna().astype(float)
    if a.empty or b.empty:
        return np.nan
    var_a = a.var(ddof=1) if len(a) > 1 else 0.0
    var_b = b.var(ddof=1) if len(b) > 1 else 0.0
    pooled = math.sqrt(((len(a) - 1) * var_a + (len(b) - 1) * var_b) / max(len(a) + len(b) - 2, 1))
    if pooled == 0 or not np.isfinite(pooled):
        return np.nan
    return float((a.mean() - b.mean()) / pooled)


def compute_cellrank_sox4_association(
    cells: pd.DataFrame,
    sox4_col: str = "C_sox4_axis",
    fate_col: str = "cellrank_fate_prob_cnv_supported_malignant",
) -> pd.DataFrame:
    if sox4_col not in cells.columns or fate_col not in cells.columns:
        return pd.DataFrame(
            [
                {
                    "status": "missing_required_columns",
                    "n_cells": 0,
                    "spearman_rho": np.nan,
                    "spearman_pvalue": np.nan,
                    "spearman_p.adjust": np.nan,
                    "high_fate_effect_size": np.nan,
                    "high_fate_pvalue": np.nan,
                    "high_fate_p.adjust": np.nan,
                    "high_sox4_high_fate_odds_ratio": np.nan,
                    "high_sox4_high_fate_pvalue": np.nan,
                    "high_sox4_high_fate_p.adjust": np.nan,
                }
            ]
        )
    work = cells[[sox4_col, fate_col]].replace([np.inf, -np.inf], np.nan).dropna()
    if work.shape[0] < 10 or work[sox4_col].nunique() < 2 or work[fate_col].nunique() < 2:
        return pd.DataFrame(
            [
                {
                    "status": "insufficient",
                    "n_cells": int(work.shape[0]),
                    "spearman_rho": np.nan,
                    "spearman_pvalue": np.nan,
                    "spearman_p.adjust": np.nan,
                    "high_fate_effect_size": np.nan,
                    "high_fate_pvalue": np.nan,
                    "high_fate_p.adjust": np.nan,
                    "high_sox4_high_fate_odds_ratio": np.nan,
                    "high_sox4_high_fate_pvalue": np.nan,
                    "high_sox4_high_fate_p.adjust": np.nan,
                }
            ]
        )
    rho, spearman_p = spearmanr(work[sox4_col], work[fate_col])
    fate_cutoff = work[fate_col].quantile(0.75)
    sox4_cutoff = work[sox4_col].quantile(0.75)
    high_fate = work[fate_col] >= fate_cutoff
    high_sox4 = work[sox4_col] >= sox4_cutoff
    effect = cohen_d(work.loc[high_fate, sox4_col], work.loc[~high_fate, sox4_col])
    if high_fate.sum() > 0 and (~high_fate).sum() > 0:
        high_fate_p = float(mannwhitneyu(work.loc[high_fate, sox4_col], work.loc[~high_fate, sox4_col], alternative="greater").pvalue)
    else:
        high_fate_p = np.nan
    contingency = np.array(
        [
            [int((~high_sox4 & ~high_fate).sum()), int((~high_sox4 & high_fate).sum())],
            [int((high_sox4 & ~high_fate).sum()), int((high_sox4 & high_fate).sum())],
        ],
        dtype=int,
    )
    odds_ratio, fisher_p = fisher_exact(contingency, alternative="greater")
    adjusted = benjamini_hochberg([spearman_p, high_fate_p, fisher_p])
    return pd.DataFrame(
        [
            {
                "status": "tested",
                "n_cells": int(work.shape[0]),
                "fate_high_quantile": 0.75,
                "sox4_high_quantile": 0.75,
                "spearman_rho": float(rho),
                "spearman_pvalue": float(spearman_p),
                "spearman_p.adjust": adjusted[0],
                "high_fate_effect_size": effect,
                "high_fate_pvalue": high_fate_p,
                "high_fate_p.adjust": adjusted[1],
                "high_sox4_high_fate_odds_ratio": float(odds_ratio),
                "high_sox4_high_fate_pvalue": float(fisher_p),
                "high_sox4_high_fate_p.adjust": adjusted[2],
            }
        ]
    )


def audit_velocity_feasibility(paths: Path | Sequence[Path] | None) -> pd.DataFrame:
    if paths is None:
        return pd.DataFrame(
            [
                {
                    "path": "",
                    "exists": False,
                    "layers": "",
                    "has_spliced_unspliced": False,
                    "velocity_status": "not_testable_missing_velocity_h5ad",
                }
            ]
        )
    if isinstance(paths, (str, Path)):
        path_list = [Path(paths)]
    else:
        path_list = [Path(path) for path in paths]
    rows = []
    for path in path_list:
        if not path.exists():
            rows.append(
                {
                    "path": str(path),
                    "exists": False,
                    "layers": "",
                    "has_spliced_unspliced": False,
                    "velocity_status": "not_testable_missing_velocity_h5ad",
                }
            )
            continue
        try:
            with h5py.File(path, "r") as handle:
                layers = list(handle.get("layers", {}).keys()) if "layers" in handle else []
        except Exception as exc:
            rows.append(
                {
                    "path": str(path),
                    "exists": True,
                    "layers": "",
                    "has_spliced_unspliced": False,
                    "velocity_status": f"not_testable_read_error:{type(exc).__name__}",
                }
            )
            continue
        has_layers = {"spliced", "unspliced"}.issubset(set(layers))
        rows.append(
            {
                "path": str(path),
                "exists": True,
                "layers": ";".join(layers),
                "has_spliced_unspliced": bool(has_layers),
                "velocity_status": "testable_spliced_unspliced_present"
                if has_layers
                else "not_testable_missing_spliced_unspliced_layers",
            }
        )
    return pd.DataFrame(rows)


def discover_h5ad_files(processed_dir: Path) -> list[Path]:
    if not processed_dir.exists():
        return []
    return sorted(processed_dir.rglob("*.h5ad"))


def normalize_gene_symbols(values: Iterable[object]) -> list[str]:
    return ["" if pd.isna(value) else str(value).strip().upper().split(".")[0] for value in values]


def select_sox4_target_genes(tf_targets_path: Path, state_specific_path: Path, top_n: int) -> list[str]:
    genes: list[str] = []
    tf_targets = read_tsv_or_empty(tf_targets_path)
    if not tf_targets.empty and {"tf", "gene"}.issubset(tf_targets.columns):
        sub = tf_targets.loc[tf_targets["tf"].astype(str).str.upper().eq("SOX4")].copy()
        if "rank" in sub.columns:
            sub["rank"] = pd.to_numeric(sub["rank"], errors="coerce")
            sub = sub.sort_values("rank")
        genes.extend(normalize_gene_symbols(sub["gene"].head(top_n)))
    state_specific = read_tsv_or_empty(state_specific_path)
    if not state_specific.empty and {"tf", "gene"}.issubset(state_specific.columns):
        sub = state_specific.loc[state_specific["tf"].astype(str).str.upper().eq("SOX4")].copy()
        sort_cols = [column for column in ["malignant_like_specificity_ratio", "malignant_like_fdr"] if column in sub.columns]
        if sort_cols:
            for column in sort_cols:
                sub[column] = pd.to_numeric(sub[column], errors="coerce")
            ascending = [False if column == "malignant_like_specificity_ratio" else True for column in sort_cols]
            sub = sub.sort_values(sort_cols, ascending=ascending)
        genes.extend(normalize_gene_symbols(sub["gene"].head(top_n)))
    return [gene for gene in dict.fromkeys(genes) if gene]


def read_expression_features(h5ad_path: Path, genes: Sequence[str]) -> pd.DataFrame:
    import anndata as ad

    if not h5ad_path.exists():
        raise FileNotFoundError(h5ad_path)
    wanted = [gene.upper() for gene in dict.fromkeys(genes) if gene]
    backed = ad.read_h5ad(h5ad_path, backed="r")
    try:
        var_upper = pd.Index(normalize_gene_symbols(backed.var_names.astype(str)))
        gene_to_pos = {gene: idx for idx, gene in enumerate(var_upper) if gene not in {"", None}}
        available = [gene for gene in wanted if gene in gene_to_pos]
        if not available:
            return pd.DataFrame(index=backed.obs_names.astype(str))
        positions = [gene_to_pos[gene] for gene in available]
        sub = backed[:, positions].to_memory()
        x = sub.X
        if sparse.issparse(x):
            x = x.toarray()
        expr = pd.DataFrame(np.asarray(x, dtype=np.float32), index=sub.obs_names.astype(str), columns=available)
    finally:
        if getattr(backed, "isbacked", False):
            backed.file.close()
    return expr


def read_regulon_auc_features(path: Path, source: str, tfs: Sequence[str]) -> pd.DataFrame:
    auc = read_tsv_or_empty(path)
    if auc.empty or "cell_id" not in auc.columns:
        return pd.DataFrame()
    keep = ["cell_id"]
    rename = {}
    for tf in tfs:
        column = f"{tf}(+)"
        if column in auc.columns:
            keep.append(column)
            rename[column] = f"regulon_{source}_{tf}"
    if len(keep) == 1:
        return pd.DataFrame(index=auc["cell_id"].astype(str))
    out = auc[keep].rename(columns=rename)
    return out.set_index(out["cell_id"].astype(str)).drop(columns=["cell_id"])


def build_base_feature_values(
    driver_h5ad: Path,
    pyscenic_auc: Path,
    cistarget_auc: Path,
    tf_targets: Path,
    sox4_state_specific: Path,
    cellrank_fate: Path,
    top_n_sox4_targets: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    sox4_targets = select_sox4_target_genes(tf_targets, sox4_state_specific, top_n=top_n_sox4_targets)
    expression_genes = list(dict.fromkeys([*TF_GENES, *sox4_targets]))
    expr = read_expression_features(driver_h5ad, expression_genes)
    features = pd.DataFrame(index=expr.index)
    for tf in TF_GENES:
        if tf in expr.columns:
            features[f"expr_{tf}"] = expr[tf]
    sox4_available = [gene for gene in sox4_targets if gene in expr.columns]
    if sox4_available:
        features["sox4_target_signature"] = mean_zscore(expr, sox4_available)

    for path, source in [(pyscenic_auc, "pyscenic"), (cistarget_auc, "cistarget")]:
        regulons = read_regulon_auc_features(path, source, TF_GENES)
        if not regulons.empty:
            features = features.join(regulons, how="left")

    fate = read_tsv_or_empty(cellrank_fate, usecols=lambda col: col in {"cell_id", "cellrank_fate_prob_cnv_supported_malignant"})
    if not fate.empty and "cell_id" in fate.columns:
        fate = fate.drop_duplicates("cell_id").set_index(fate["cell_id"].astype(str)).drop(columns=["cell_id"])
        features = features.join(fate, how="left")

    feature_manifest = []
    for column in features.columns:
        feature_manifest.append(
            {
                "feature": column,
                "n_cells": int(features[column].notna().sum()),
                "feature_class": column.split("_", 1)[0],
            }
        )
    return features, pd.DataFrame(feature_manifest)


def build_temporal_cell_scores(
    run_method_long: pd.DataFrame,
    base_features: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    cell_score_rows = []
    availability_rows = []
    for (run_id, method), sub in run_method_long.groupby(["run_id", "method"], observed=True):
        sub = sub.copy()
        sub["cell_id"] = sub["cell_id"].astype(str)
        feature_values = base_features.reindex(sub["cell_id"]).copy()
        for module in ["HCC_Malignant_Associated", "Proliferation"]:
            if module in sub.columns:
                feature_values[f"module_{module}"] = pd.to_numeric(sub[module].to_numpy(), errors="coerce")
        scores, availability = build_axis_scores(feature_values)
        out = scores.reset_index(names="cell_id")
        out.insert(1, "run_id", run_id)
        out.insert(2, "method", method)
        metadata_cols = [
            "pseudotime_norm",
            "dataset",
            "sample_id",
            "study_sample",
            "cnv_sample",
            "sample_source_class",
            "cell_disease_stage",
            "trajectory_role",
            "trajectory_root_end_role",
            "consensus_pseudotime_phase",
            "driver_eligible",
        ]
        for column in metadata_cols:
            if column in sub.columns:
                out[column] = sub[column].to_numpy()
        for raw_feature in ["cellrank_fate_prob_cnv_supported_malignant"]:
            if raw_feature in feature_values.columns:
                out[raw_feature] = feature_values[raw_feature].to_numpy()
        cell_score_rows.append(out)
        availability.insert(0, "method", method)
        availability.insert(0, "run_id", run_id)
        availability_rows.append(availability)
    return pd.concat(cell_score_rows, ignore_index=True), pd.concat(availability_rows, ignore_index=True)


def summarize_all_temporal_trends(cell_scores: pd.DataFrame, n_bins: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    binned_rows = []
    onset_rows = []
    for (run_id, method), sub in cell_scores.groupby(["run_id", "method"], observed=True):
        pseudotime = pd.to_numeric(sub["pseudotime_norm"], errors="coerce")
        for feature in ORDER_FEATURES:
            binned = bin_temporal_values(pd.to_numeric(sub[feature], errors="coerce"), pseudotime, n_bins=n_bins, feature=feature)
            if not binned.empty:
                binned.insert(0, "method", method)
                binned.insert(0, "run_id", run_id)
                binned_rows.append(binned)
            summary = summarize_temporal_trend(sub[feature], pseudotime, n_bins=n_bins, feature=feature)
            onset_rows.append({"run_id": run_id, "method": method, **summary})
    trend_bins = pd.concat(binned_rows, ignore_index=True) if binned_rows else pd.DataFrame()
    onset_times = pd.DataFrame(onset_rows)
    return trend_bins, onset_times


def run_grouped_bootstrap_by_method(
    cell_scores: pd.DataFrame,
    n_bootstrap: int,
    random_state: int,
    n_bins: int,
) -> pd.DataFrame:
    rows = []
    for idx, ((run_id, method), sub) in enumerate(cell_scores.groupby(["run_id", "method"], observed=True)):
        tests = grouped_bootstrap_order_tests(
            sub,
            pseudotime_col="pseudotime_norm",
            group_col="cnv_sample",
            n_bootstrap=n_bootstrap,
            random_state=random_state + idx,
            n_bins=n_bins,
        )
        tests.insert(0, "method", method)
        tests.insert(0, "run_id", run_id)
        rows.append(tests)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def build_evidence_grade(
    bootstrap: pd.DataFrame,
    cellrank_assoc: pd.DataFrame,
    velocity: pd.DataFrame,
) -> pd.DataFrame:
    if bootstrap.empty:
        supported_methods = 0
        total_methods = 0
        mean_order_probability = np.nan
    else:
        method_support = (
            bootstrap.assign(is_supported=bootstrap["order_probability"].astype(float).ge(0.7))
            .groupby(["run_id", "method"], observed=True)["is_supported"]
            .all()
            .reset_index()
        )
        supported_methods = int(method_support["is_supported"].sum())
        total_methods = int(method_support.shape[0])
        mean_order_probability = float(pd.to_numeric(bootstrap["order_probability"], errors="coerce").mean())

    cellrank_supported = False
    if not cellrank_assoc.empty and cellrank_assoc.loc[0, "status"] == "tested":
        cellrank_supported = (
            float(cellrank_assoc.loc[0, "spearman_rho"]) > 0
            and float(cellrank_assoc.loc[0, "spearman_p.adjust"]) < 0.05
            and float(cellrank_assoc.loc[0, "high_fate_effect_size"]) > 0
        )
    velocity_statuses = set(velocity["velocity_status"].astype(str)) if not velocity.empty else {"not_testable_missing_velocity_h5ad"}
    velocity_testable = any(status == "testable_spliced_unspliced_present" for status in velocity_statuses)
    velocity_not_testable = not velocity_testable

    temporal_supported = total_methods > 0 and supported_methods / total_methods >= 0.5 and mean_order_probability >= 0.7
    if temporal_supported and cellrank_supported and velocity_not_testable:
        final = "temporal_supported_velocity_not_testable"
    elif temporal_supported and cellrank_supported:
        final = "temporal_supported"
    elif temporal_supported:
        final = "pseudotime_supported_cellrank_not_supported"
    else:
        final = "temporal_not_supported"
    return pd.DataFrame(
        [
            {
                "evidence_domain": "Module9.1_temporal_hierarchy",
                "supported_methods": supported_methods,
                "total_methods": total_methods,
                "mean_order_probability": mean_order_probability,
                "cellrank_sox4_fate_supported": bool(cellrank_supported),
                "velocity_testable": bool(velocity_testable),
                "velocity_status_summary": ";".join(sorted(velocity_statuses)),
                "final_support_label": final,
            }
        ]
    )


def configure_plot_style() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica"],
            "font.size": 8,
            "axes.labelsize": 8,
            "axes.titlesize": 9,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def save_figure(fig, path_base: Path) -> dict[str, str]:
    outputs = {}
    path_base.parent.mkdir(parents=True, exist_ok=True)
    for suffix in ["png", "pdf", "svg"]:
        path = path_base.with_suffix(f".{suffix}")
        fig.savefig(path, bbox_inches="tight")
        outputs[suffix] = str(path.resolve())
    return outputs


def write_figures(
    trend_bins: pd.DataFrame,
    onset_times: pd.DataFrame,
    bootstrap: pd.DataFrame,
    cell_scores: pd.DataFrame,
    cellrank_assoc: pd.DataFrame,
    evidence: pd.DataFrame,
    figure_dir: Path,
) -> dict[str, str]:
    configure_plot_style()
    import matplotlib.pyplot as plt

    outputs: dict[str, str] = {}
    colors = {
        "A_hnf4a_ppara_loss": "#0072B2",
        "B_transition_activation": "#D55E00",
        "C_sox4_axis": "#CC79A7",
        "C_malignant_like_fate": "#009E73",
    }
    if not trend_bins.empty:
        fig, ax = plt.subplots(figsize=(6.2, 4.0))
        primary = trend_bins.loc[trend_bins["run_id"].eq("main_strict")].copy()
        if primary.empty:
            primary = trend_bins.copy()
        summary = (
            primary.groupby(["feature", "bin"], observed=True)
            .agg(mean_pseudotime=("mean_pseudotime", "mean"), smoothed_score=("smoothed_score", "mean"))
            .reset_index()
        )
        for feature, sub in summary.groupby("feature", observed=True):
            sub = sub.sort_values("mean_pseudotime")
            ax.plot(
                sub["mean_pseudotime"],
                sub["smoothed_score"],
                marker="o",
                markersize=2.5,
                linewidth=1.8,
                label=feature,
                color=colors.get(feature),
            )
        ax.set_xlabel("Normalized pseudotime")
        ax.set_ylabel("Smoothed axis score")
        ax.set_title("Module 9.1 three-axis temporal trends")
        ax.legend(frameon=False, fontsize=6)
        ax.grid(alpha=0.2)
        fig.tight_layout()
        outputs.update({f"three_axis_pseudotime_trends_{k}": v for k, v in save_figure(fig, figure_dir / "module9_1_three_axis_pseudotime_trends").items()})
        plt.close(fig)

    if not onset_times.empty:
        fig, ax = plt.subplots(figsize=(6.2, 3.6))
        plot = onset_times.loc[onset_times["feature"].isin(ORDER_FEATURES)].copy()
        plot["run_method"] = plot["run_id"].astype(str) + " / " + plot["method"].astype(str)
        y_labels = plot["run_method"].drop_duplicates().tolist()
        y_pos = {label: idx for idx, label in enumerate(y_labels)}
        for feature, sub in plot.groupby("feature", observed=True):
            ax.scatter(
                sub["onset_time"],
                sub["run_method"].map(y_pos),
                label=feature,
                color=colors.get(feature),
                s=18,
                alpha=0.8,
            )
        ax.set_yticks(list(y_pos.values()))
        ax.set_yticklabels(y_labels)
        ax.set_xlabel("Onset pseudotime")
        ax.set_title("Onset order across trajectory methods")
        ax.legend(frameon=False, fontsize=6, ncol=2)
        ax.grid(axis="x", alpha=0.2)
        fig.tight_layout()
        outputs.update({f"onset_order_forest_{k}": v for k, v in save_figure(fig, figure_dir / "module9_1_onset_order_forest").items()})
        plt.close(fig)

    fate_col = "cellrank_fate_prob_cnv_supported_malignant"
    if fate_col in cell_scores.columns:
        assoc_cells = (
            cell_scores.loc[cell_scores["method"].eq("monocle3"), ["cell_id", "C_sox4_axis", fate_col]]
            .drop_duplicates("cell_id")
            .replace([np.inf, -np.inf], np.nan)
            .dropna()
        )
        if not assoc_cells.empty:
            sample = assoc_cells.sample(n=min(3000, assoc_cells.shape[0]), random_state=1)
            fig, ax = plt.subplots(figsize=(4.4, 3.6))
            ax.scatter(sample["C_sox4_axis"], sample[fate_col], s=7, alpha=0.35, linewidths=0)
            rho = cellrank_assoc.loc[0, "spearman_rho"] if not cellrank_assoc.empty else np.nan
            ax.set_xlabel("SOX4 axis score")
            ax.set_ylabel("CellRank CNV malignant fate probability")
            ax.set_title(f"SOX4 axis vs malignant fate (rho={float(rho):.2f})" if pd.notna(rho) else "SOX4 axis vs malignant fate")
            ax.grid(alpha=0.2)
            fig.tight_layout()
            outputs.update({f"sox4_cellrank_fate_association_{k}": v for k, v in save_figure(fig, figure_dir / "module9_1_sox4_cellrank_fate_association").items()})
            plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2))
    if not bootstrap.empty:
        boot_summary = bootstrap.groupby("comparison", observed=True)["order_probability"].mean().reset_index()
        axes[0].barh(boot_summary["comparison"], boot_summary["order_probability"], color="#0072B2")
        axes[0].axvline(0.7, color="black", linestyle="--", linewidth=1)
        axes[0].set_xlim(0, 1)
        axes[0].set_xlabel("Order probability")
    axes[0].set_title("Bootstrap order support")
    if not evidence.empty:
        labels = ["Temporal methods", "CellRank", "Velocity testable"]
        values = [
            float(evidence.loc[0, "supported_methods"]) / max(float(evidence.loc[0, "total_methods"]), 1.0),
            1.0 if bool(evidence.loc[0, "cellrank_sox4_fate_supported"]) else 0.0,
            1.0 if bool(evidence.loc[0, "velocity_testable"]) else 0.0,
        ]
        axes[1].bar(labels, values, color=["#009E73", "#CC79A7", "#999999"])
        axes[1].set_ylim(0, 1)
        axes[1].set_ylabel("Support fraction")
        axes[1].tick_params(axis="x", rotation=30)
    axes[1].set_title(str(evidence.loc[0, "final_support_label"]) if not evidence.empty else "Evidence summary")
    fig.tight_layout()
    outputs.update({f"temporal_evidence_summary_{k}": v for k, v in save_figure(fig, figure_dir / "module9_1_temporal_evidence_summary").items()})
    plt.close(fig)
    return outputs


def build_main_conclusions(evidence: pd.DataFrame, bootstrap: pd.DataFrame, cellrank: pd.DataFrame, velocity: pd.DataFrame) -> str:
    label = evidence.loc[0, "final_support_label"] if not evidence.empty else "not_available"
    supported_methods = int(evidence.loc[0, "supported_methods"]) if not evidence.empty else 0
    total_methods = int(evidence.loc[0, "total_methods"]) if not evidence.empty else 0
    mean_order = float(evidence.loc[0, "mean_order_probability"]) if not evidence.empty and pd.notna(evidence.loc[0, "mean_order_probability"]) else np.nan
    lines = [
        "# Module 9.1 Temporal Hierarchy Conclusions",
        "",
        "## Current Status",
        f"Final temporal support label: {label}.",
        "",
        "## Pseudotime Ordering",
        f"Supported trajectory methods: {supported_methods}/{total_methods}. Mean bootstrap order probability: {mean_order:.3f}.",
    ]
    if not bootstrap.empty:
        for comparison, sub in bootstrap.groupby("comparison", observed=True):
            lines.append(
                f"- {comparison}: mean order probability {pd.to_numeric(sub['order_probability'], errors='coerce').mean():.3f}."
            )
    lines.extend(["", "## CellRank SOX4 Fate Evidence"])
    if not cellrank.empty and cellrank.loc[0, "status"] == "tested":
        lines.append(
            "SOX4 axis was positively associated with CNV-supported malignant fate "
            f"(Spearman rho {float(cellrank.loc[0, 'spearman_rho']):.3f}, "
            f"FDR {float(cellrank.loc[0, 'spearman_p.adjust']):.3g})."
        )
    else:
        lines.append("CellRank SOX4 fate association was not testable.")
    lines.extend(["", "## RNA Velocity"])
    statuses = sorted(set(velocity["velocity_status"].astype(str))) if not velocity.empty else ["not_testable_missing_velocity_h5ad"]
    lines.append(f"Velocity feasibility status: {'; '.join(statuses)}.")
    lines.append("Current project H5AD files lack paired spliced/unspliced layers, so velocity timing remains an explicit evidence gap.")
    return "\n".join(lines) + "\n"


def run_module9_1(
    driver_h5ad: Path = DEFAULT_DRIVER_H5AD,
    run_method_long_path: Path = DEFAULT_RUN_METHOD_LONG,
    cellrank_fate_path: Path = DEFAULT_CELLRANK_FATE,
    pyscenic_auc_path: Path = DEFAULT_PYSCENIC_AUC,
    cistarget_auc_path: Path = DEFAULT_CISTARGET_AUC,
    tf_targets_path: Path = DEFAULT_TF_TARGETS,
    sox4_state_specific_path: Path = DEFAULT_SOX4_STATE_SPECIFIC,
    metadata_dir: Path = DEFAULT_METADATA_DIR,
    figure_dir: Path = DEFAULT_FIGURE_DIR,
    processed_dir: Path = DEFAULT_PROCESSED_DIR,
    velocity_h5ad: Path | None = None,
    top_n_sox4_targets: int = 50,
    n_bins: int = 20,
    n_bootstrap: int = 500,
    seed: int = 20260615,
) -> dict[str, object]:
    start = time.time()
    metadata_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    base_features, feature_manifest = build_base_feature_values(
        driver_h5ad=driver_h5ad,
        pyscenic_auc=pyscenic_auc_path,
        cistarget_auc=cistarget_auc_path,
        tf_targets=tf_targets_path,
        sox4_state_specific=sox4_state_specific_path,
        cellrank_fate=cellrank_fate_path,
        top_n_sox4_targets=top_n_sox4_targets,
    )
    run_method_long = read_tsv_or_empty(run_method_long_path)
    if run_method_long.empty:
        raise ValueError(f"No run-method pseudotime rows found: {run_method_long_path}")
    cell_scores, availability = build_temporal_cell_scores(run_method_long, base_features)
    trend_bins, onset_times = summarize_all_temporal_trends(cell_scores, n_bins=n_bins)
    bootstrap = run_grouped_bootstrap_by_method(
        cell_scores,
        n_bootstrap=n_bootstrap,
        random_state=seed,
        n_bins=n_bins,
    )
    cellrank_assoc_cells = (
        cell_scores.loc[cell_scores["method"].eq("monocle3"), ["cell_id", "C_sox4_axis", "cellrank_fate_prob_cnv_supported_malignant"]]
        .drop_duplicates("cell_id")
        .copy()
    )
    cellrank_assoc = compute_cellrank_sox4_association(cellrank_assoc_cells)

    if velocity_h5ad is not None:
        velocity_paths = [velocity_h5ad]
    else:
        velocity_paths = discover_h5ad_files(processed_dir)
    velocity = audit_velocity_feasibility(velocity_paths)
    evidence = build_evidence_grade(bootstrap, cellrank_assoc, velocity)
    figure_outputs = write_figures(trend_bins, onset_times, bootstrap, cell_scores, cellrank_assoc, evidence, figure_dir)

    outputs = {
        "temporal_cell_scores": str((metadata_dir / "module9_1_temporal_cell_scores.tsv.gz").resolve()),
        "axis_feature_availability": str((metadata_dir / "module9_1_axis_feature_availability.tsv").resolve()),
        "feature_manifest": str((metadata_dir / "module9_1_feature_manifest.tsv").resolve()),
        "pseudotime_trends": str((metadata_dir / "module9_1_pseudotime_trends.tsv").resolve()),
        "onset_times": str((metadata_dir / "module9_1_onset_times.tsv").resolve()),
        "bootstrap_order_tests": str((metadata_dir / "module9_1_bootstrap_order_tests.tsv").resolve()),
        "cellrank_sox4_fate_association": str((metadata_dir / "module9_1_cellrank_sox4_fate_association.tsv").resolve()),
        "velocity_feasibility": str((metadata_dir / "module9_1_velocity_feasibility.tsv").resolve()),
        "evidence_grade": str((metadata_dir / "module9_1_evidence_grade.tsv").resolve()),
        "report": str((metadata_dir / "module9_1_report.json").resolve()),
        "main_conclusions": str((metadata_dir / "module9_1_main_conclusions.md").resolve()),
        **figure_outputs,
    }
    cell_scores.to_csv(outputs["temporal_cell_scores"], sep="\t", index=False, compression="gzip")
    availability.to_csv(outputs["axis_feature_availability"], sep="\t", index=False)
    feature_manifest.to_csv(outputs["feature_manifest"], sep="\t", index=False)
    trend_bins.to_csv(outputs["pseudotime_trends"], sep="\t", index=False)
    onset_times.to_csv(outputs["onset_times"], sep="\t", index=False)
    bootstrap.to_csv(outputs["bootstrap_order_tests"], sep="\t", index=False)
    cellrank_assoc.to_csv(outputs["cellrank_sox4_fate_association"], sep="\t", index=False)
    velocity.to_csv(outputs["velocity_feasibility"], sep="\t", index=False)
    evidence.to_csv(outputs["evidence_grade"], sep="\t", index=False)
    Path(outputs["main_conclusions"]).write_text(
        build_main_conclusions(evidence, bootstrap, cellrank_assoc, velocity),
        encoding="utf-8",
    )

    report = {
        "module": "9.1",
        "method": "Temporal hierarchy evidence using pseudotime, regulon/expression axis scores, CellRank fate probabilities and velocity feasibility audit",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "driver_h5ad": str(driver_h5ad.resolve()),
            "run_method_long": str(run_method_long_path.resolve()),
            "cellrank_fate": str(cellrank_fate_path.resolve()),
            "pyscenic_auc": str(pyscenic_auc_path.resolve()),
            "cistarget_auc": str(cistarget_auc_path.resolve()),
            "tf_targets": str(tf_targets_path.resolve()),
            "sox4_state_specific": str(sox4_state_specific_path.resolve()),
            "velocity_h5ad": str(velocity_h5ad.resolve()) if velocity_h5ad else "",
        },
        "parameters": {
            "top_n_sox4_targets": int(top_n_sox4_targets),
            "n_bins": int(n_bins),
            "n_bootstrap": int(n_bootstrap),
            "seed": int(seed),
        },
        "outputs": outputs,
        "n_temporal_cell_score_rows": int(cell_scores.shape[0]),
        "n_onset_rows": int(onset_times.shape[0]),
        "n_bootstrap_rows": int(bootstrap.shape[0]),
        "n_velocity_files_audited": int(velocity.shape[0]),
        "final_support_label": str(evidence.loc[0, "final_support_label"]) if not evidence.empty else "not_available",
        "package_versions": {
            "python": platform.python_version(),
            "pandas": package_version("pandas"),
            "numpy": package_version("numpy"),
            "scipy": package_version("scipy"),
            "anndata": package_version("anndata"),
            "statsmodels": package_version("statsmodels"),
            "matplotlib": package_version("matplotlib"),
        },
        "elapsed_seconds": round(time.time() - start, 3),
    }
    Path(outputs["report"]).write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def main() -> None:
    args = parse_args()
    report = run_module9_1(
        driver_h5ad=args.driver_h5ad,
        run_method_long_path=args.run_method_long,
        cellrank_fate_path=args.cellrank_fate,
        pyscenic_auc_path=args.pyscenic_auc,
        cistarget_auc_path=args.cistarget_auc,
        tf_targets_path=args.tf_targets,
        sox4_state_specific_path=args.sox4_state_specific,
        metadata_dir=args.metadata_dir,
        figure_dir=args.figure_dir,
        processed_dir=args.processed_dir,
        velocity_h5ad=args.velocity_h5ad,
        top_n_sox4_targets=args.top_n_sox4_targets,
        n_bins=args.n_bins,
        n_bootstrap=args.n_bootstrap,
        seed=args.seed,
    )
    print(
        json.dumps(
            {
                "report": report["outputs"]["report"],
                "final_support_label": report["final_support_label"],
                "n_temporal_cell_score_rows": report["n_temporal_cell_score_rows"],
                "n_velocity_files_audited": report["n_velocity_files_audited"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
