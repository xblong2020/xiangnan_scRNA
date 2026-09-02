from __future__ import annotations

import argparse
import ast
import json
import math
import time
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from statsmodels.nonparametric.smoothers_lowess import lowess


ROOT = Path(__file__).resolve().parents[1]


def record_path(path: Path | str) -> str:
    value = Path(path)
    try:
        return str(value.resolve().relative_to(ROOT.resolve())).replace("\\", "/")
    except ValueError:
        return str(value)

PROGRAMS = {
    "identity_maintenance": {"HNF4A", "PPARA", "HLF", "AR", "PPARGC1A"},
    "AP1_stress_transition": {"JUN", "JUNB", "JUND", "FOS", "FOSB"},
    "CEBPB_EGR1_transition": {"CEBPB", "CEBPD", "EGR1"},
    "SOX4_malignant_stabilization": {"SOX4"},
}

AXES = {
    "HNF4A_PPARA_identity": ["HNF4A", "PPARA", "HLF"],
    "AP1_stress_transition": ["JUN", "JUNB", "JUND", "FOS", "FOSB"],
    "CEBPB_EGR1_transition": ["CEBPB", "CEBPD", "EGR1"],
    "SOX4_malignant_stabilization": ["SOX4"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Integrate formal canonical SCENIC 6.3b with CellRank and robustness.")
    parser.add_argument(
        "--input-h5ad",
        type=Path,
        default=ROOT / "data/processed/driver/driver_union_full_expression.module6_3b.formal.h5ad",
    )
    parser.add_argument(
        "--ctx-output",
        type=Path,
        default=ROOT / "metadata/driver/scenic_module6_3b/driver_module6_3b_canonical_regulons_seed777.tsv",
    )
    parser.add_argument(
        "--auc",
        type=Path,
        default=ROOT / "metadata/driver/scenic_module6_3b/driver_module6_3b_canonical_regulon_auc_9512.csv",
    )
    parser.add_argument(
        "--old-summary",
        type=Path,
        default=ROOT / "metadata/driver/driver_module6_3_pyscenic_regulons.tsv",
    )
    parser.add_argument(
        "--trajectory-h5ad",
        type=Path,
        default=ROOT / "data/processed/driver/driver_hepatocyte_trajectory.module6_1.h5ad",
    )
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/driver/scenic_module6_3b")
    parser.add_argument("--processed-dir", type=Path, default=ROOT / "data/processed/driver")
    parser.add_argument("--figures-dir", type=Path, default=ROOT / "figures/driver/scenic_module6_3b")
    parser.add_argument("--fate-key", default="cellrank_fate_prob_cnv_supported_malignant")
    parser.add_argument("--pseudotime-key", default="driver_main_strict__pseudotime_median")
    parser.add_argument("--phase-key", default="driver_main_strict__pseudotime_phase")
    parser.add_argument("--dataset-key", default="dataset")
    parser.add_argument("--sample-key", default="sample_id")
    parser.add_argument("--cnv-key", default="driver_primary_module3_cnv_supported")
    parser.add_argument("--min-group-cells", type=int, default=20)
    parser.add_argument("--top-n", type=int, default=20)
    return parser.parse_args()


def package_version(name: str) -> str:
    try:
        return version(name)
    except PackageNotFoundError:
        return "not_installed"


def bh_qvalues(pvalues: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvalues, errors="coerce")
    q = pd.Series(np.nan, index=p.index, dtype=float)
    valid = p.dropna()
    if valid.empty:
        return q
    order = valid.sort_values().index
    values = valid.loc[order].to_numpy(dtype=float)
    adjusted = values * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    q.loc[order] = np.clip(adjusted, 0.0, 1.0)
    return q


def parse_targets(value: object) -> set[str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return set()
    try:
        parsed = ast.literal_eval(str(value))
    except (ValueError, SyntaxError):
        return {part.strip() for part in str(value).replace(",", ";").split(";") if part.strip()}
    genes = set()
    for item in parsed:
        if isinstance(item, (list, tuple)) and item:
            genes.add(str(item[0]))
        elif isinstance(item, str):
            genes.add(item)
    return genes


def read_ctx(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", header=[0, 1], index_col=[0, 1])
    if isinstance(frame.columns, pd.MultiIndex):
        frame.columns = [str(bottom) if not str(bottom).startswith("Unnamed") else str(top) for top, bottom in frame.columns]
    return frame.reset_index()


def read_auc(path: Path, cell_ids: pd.Index) -> pd.DataFrame:
    auc = pd.read_csv(path)
    cell_col = next((column for column in ["Cell", "CellID", "cell_id", "cells"] if column in auc.columns), None)
    if cell_col:
        auc.index = auc.pop(cell_col).astype(str)
    else:
        auc.index = auc.index.astype(str)
    missing = pd.Index(cell_ids.astype(str)).difference(auc.index)
    if len(missing):
        raise ValueError(f"{len(missing)} cells are missing from AUCell output")
    auc = auc.loc[cell_ids.astype(str)].apply(pd.to_numeric, errors="coerce")
    auc = auc.loc[:, auc.nunique(dropna=True) > 1]
    if auc.empty:
        raise ValueError("No variable regulons in AUCell output")
    return auc


def safe_corr(x: pd.Series, y: pd.Series, method: str) -> tuple[float, float, int]:
    mask = x.notna() & y.notna() & np.isfinite(x) & np.isfinite(y)
    n = int(mask.sum())
    if n < 3 or x.loc[mask].nunique() < 2 or y.loc[mask].nunique() < 2:
        return np.nan, np.nan, n
    fn = spearmanr if method == "spearman" else pearsonr
    result = fn(x.loc[mask], y.loc[mask])
    return float(result.statistic if hasattr(result, "statistic") else result[0]), float(result.pvalue if hasattr(result, "pvalue") else result[1]), n


def cnv_positive_mask(cells: pd.DataFrame, cnv_key: str) -> pd.Series:
    if cnv_key in cells.columns:
        values = cells[cnv_key]
        if pd.api.types.is_bool_dtype(values):
            return values.fillna(False).astype(bool)
        return values.astype(str).str.lower().isin({"true", "1", "1.0", "yes", "cnv_supported_malignant"})
    terminal = cells.get("cellrank_terminal_state", pd.Series("", index=cells.index))
    return terminal.astype(str).eq("cnv_supported_malignant")


def association_table(auc: pd.DataFrame, cells: pd.DataFrame, fate_key: str, pseudotime_key: str, phase_key: str, cnv_key: str) -> pd.DataFrame:
    fate = pd.to_numeric(cells[fate_key], errors="coerce")
    pseudotime = pd.to_numeric(cells[pseudotime_key], errors="coerce")
    fate_mask = fate.notna()
    cnv = cnv_positive_mask(cells, cnv_key)
    rows = []
    for regulon in auc.columns:
        values = pd.to_numeric(auc[regulon], errors="coerce")
        sr, sp, sn = safe_corr(values.loc[fate_mask], fate.loc[fate_mask], "spearman")
        pr, pp, _ = safe_corr(values.loc[fate_mask], fate.loc[fate_mask], "pearson")
        ts, tsp, tn = safe_corr(values, pseudotime, "spearman")
        tp, tpp, _ = safe_corr(values, pseudotime, "pearson")
        phase_means = {}
        for phase in ["early", "middle", "late"]:
            phase_values = values.loc[fate_mask & cells[phase_key].astype(str).eq(phase)] if phase_key in cells.columns else pd.Series(dtype=float)
            phase_means[f"{phase}_mean_AUC"] = float(phase_values.mean()) if not phase_values.empty else np.nan
        pos = values.loc[fate_mask & cnv]
        neg = values.loc[fate_mask & ~cnv]
        tf = regulon[:-3] if regulon.endswith("(+)") else regulon
        rows.append(
            {
                "regulon": regulon,
                "TF": tf,
                "regulon_size": np.nan,
                "spearman_rho": sr,
                "spearman_p": sp,
                "n_cells": sn,
                "pearson_r": pr,
                "pearson_p": pp,
                "pseudotime_spearman_rho": ts,
                "pseudotime_spearman_p": tsp,
                "pseudotime_pearson_r": tp,
                "pseudotime_pearson_p": tpp,
                "pseudotime_n_cells": tn,
                "CNV_positive_mean_AUC": float(pos.mean()) if not pos.empty else np.nan,
                "CNV_negative_mean_AUC": float(neg.mean()) if not neg.empty else np.nan,
                "CNV_difference_AUC": float(pos.mean() - neg.mean()) if not pos.empty and not neg.empty else np.nan,
                **phase_means,
            }
        )
    out = pd.DataFrame(rows)
    out["spearman_FDR"] = bh_qvalues(out["spearman_p"])
    out["pearson_FDR"] = bh_qvalues(out["pearson_p"])
    out["pseudotime_FDR"] = bh_qvalues(out["pseudotime_spearman_p"])
    return out.sort_values(["spearman_rho", "pseudotime_spearman_rho"], ascending=[False, False]).reset_index(drop=True)


def sample_pseudobulk(auc: pd.DataFrame, cells: pd.DataFrame, fate_key: str, sample_key: str, dataset_key: str, phase_key: str, cnv_key: str) -> pd.DataFrame:
    fate = pd.to_numeric(cells[fate_key], errors="coerce")
    keep = fate.notna()
    metadata = cells.loc[keep, [column for column in [sample_key, dataset_key, phase_key] if column in cells.columns]].copy()
    metadata["fate_probability"] = fate.loc[keep]
    metadata["cnv_state"] = np.where(cnv_positive_mask(cells, cnv_key).loc[keep], "CNV_positive", "CNV_negative")
    joined = pd.concat([metadata, auc.loc[keep]], axis=1)
    group_cols = [column for column in [sample_key, dataset_key, phase_key, "cnv_state"] if column in joined.columns]
    mean = joined.groupby(group_cols, observed=True).mean(numeric_only=True).reset_index()
    n = joined.groupby(group_cols, observed=True).size().rename("n_cells").reset_index()
    out = mean.merge(n, on=group_cols, how="left")
    return out


def pseudobulk_correlations(pseudobulk: pd.DataFrame, regulons: list[str], sample_key: str, fate_key: str = "fate_probability") -> pd.DataFrame:
    rows = []
    if pseudobulk.empty:
        return pd.DataFrame()
    group = pseudobulk.groupby(sample_key, observed=True)
    fate_sample = group[fate_key].mean()
    for regulon in regulons:
        auc_sample = group[regulon].mean()
        joined = pd.concat([auc_sample.rename("auc"), fate_sample.rename("fate")], axis=1).dropna()
        rho, p, n = safe_corr(joined["auc"], joined["fate"], "spearman")
        rows.append({"regulon": regulon, "sample_pseudobulk_spearman_rho": rho, "sample_pseudobulk_spearman_p": p, "sample_pseudobulk_n": n})
    out = pd.DataFrame(rows)
    if not out.empty:
        out["sample_pseudobulk_FDR"] = bh_qvalues(out["sample_pseudobulk_spearman_p"])
    return out


def leave_one_out(auc: pd.DataFrame, cells: pd.DataFrame, fate_key: str, group_key: str, method: str) -> pd.DataFrame:
    fate = pd.to_numeric(cells[fate_key], errors="coerce")
    keep = fate.notna()
    rows = []
    groups = sorted(cells.loc[keep, group_key].dropna().astype(str).unique()) if group_key in cells.columns else []
    for omitted in groups:
        mask = keep & cells[group_key].astype(str).ne(omitted)
        for regulon in auc.columns:
            rho, p, n = safe_corr(auc[regulon].loc[mask], fate.loc[mask], method)
            rows.append({"regulon": regulon, "omitted_group": omitted, "n_cells": n, "rho": rho, "p": p})
    out = pd.DataFrame(rows)
    if not out.empty:
        out["FDR"] = bh_qvalues(out["p"])
        out["effect_direction"] = np.sign(out["rho"])
        out["rank"] = out.groupby("omitted_group", observed=True)["rho"].rank(method="min", ascending=False)
    return out


def robustness_table(association: pd.DataFrame, sample_loo: pd.DataFrame, dataset_loo: pd.DataFrame, sample_pb: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in association.itertuples(index=False):
        expected = 1 if pd.isna(row.spearman_rho) or row.spearman_rho >= 0 else -1
        sample = sample_loo.loc[sample_loo["regulon"].eq(row.regulon)] if not sample_loo.empty else pd.DataFrame()
        dataset = dataset_loo.loc[dataset_loo["regulon"].eq(row.regulon)] if not dataset_loo.empty else pd.DataFrame()
        sample_consistency = float((sample["rho"] * expected > 0).mean()) if not sample.empty else np.nan
        dataset_consistency = float((dataset["rho"] * expected > 0).mean()) if not dataset.empty else np.nan
        sample_stable = bool(not sample.empty and sample_consistency >= 0.8 and not (sample["rho"] * expected < 0).any())
        dataset_stable = bool(not dataset.empty and dataset_consistency >= 0.8 and not (dataset["rho"] * expected < 0).any())
        support_count = int(pd.Series([row.spearman_FDR < 0.05, row.pseudotime_FDR < 0.05, sample_stable, dataset_stable]).fillna(False).sum())
        tier = "Tier A" if row.spearman_FDR < 0.05 and row.pseudotime_FDR < 0.05 and sample_stable and dataset_stable else ("Tier B" if support_count >= 2 else "Tier C")
        rows.append(
            {
                "regulon": row.regulon,
                "sample_loo_direction_consistency": sample_consistency,
                "dataset_loo_direction_consistency": dataset_consistency,
                "sample_loo_stable": sample_stable,
                "dataset_loo_stable": dataset_stable,
                "seed_stability": "single_seed_777_not_run",
                "evidence_tier": tier,
            }
        )
    return pd.DataFrame(rows)


def biological_program(tf: str) -> str:
    for program, members in PROGRAMS.items():
        if tf in members:
            return program
    return "other_canonical_regulatory_program"


def evaluate_axes(associations: pd.DataFrame) -> pd.DataFrame:
    rows = []
    present = set(associations["TF"].astype(str))
    for axis, tfs in AXES.items():
        detected = sorted(set(tfs).intersection(present))
        if not detected:
            axis_status = "not_detected"
        elif all(tf in present for tf in tfs):
            axis_status = "detected_complete_axis"
        else:
            axis_status = "detected_partial_axis"
        for tf in tfs:
            sub = associations.loc[associations["TF"].eq(tf)]
            if sub.empty:
                rows.append(
                    {
                        "axis": axis,
                        "axis_status": axis_status,
                        "TF": tf,
                        "canonical_regulon_detected": False,
                        "interpretation": "not_detected_by_canonical_ctx",
                    }
                )
                continue
            row = sub.iloc[0]
            multi_layer = bool(
                pd.notna(row.get("spearman_FDR"))
                and row["spearman_FDR"] < 0.05
                and pd.notna(row.get("pseudotime_FDR"))
                and row["pseudotime_FDR"] < 0.05
                and bool(row.get("sample_loo_stable", False))
                and bool(row.get("dataset_loo_stable", False))
            )
            rows.append(
                {
                    "axis": axis,
                    "axis_status": axis_status,
                    "TF": tf,
                    "canonical_regulon_detected": True,
                    "regulon": row.get("regulon"),
                    "regulon_size": row.get("regulon_size"),
                    "motif_NES": row.get("motif_NES_max"),
                    "spearman_rho": row.get("spearman_rho"),
                    "spearman_FDR": row.get("spearman_FDR"),
                    "pseudotime_spearman_rho": row.get("pseudotime_spearman_rho"),
                    "pseudotime_FDR": row.get("pseudotime_FDR"),
                    "sample_loo_stable": row.get("sample_loo_stable"),
                    "dataset_loo_stable": row.get("dataset_loo_stable"),
                    "evidence_tier": row.get("evidence_tier"),
                    "interpretation": "supported_by_multiple_layers" if multi_layer else "partial_or_directionally_inconsistent_support",
                }
            )
    return pd.DataFrame(rows)


def save_figure(fig: plt.Figure, path_base: Path) -> list[str]:
    path_base.parent.mkdir(parents=True, exist_ok=True)
    outputs = []
    for suffix in [".png", ".pdf"]:
        path = path_base.parent / f"{path_base.name}{suffix}"
        fig.savefig(path, bbox_inches="tight")
        outputs.append(record_path(path))
    plt.close(fig)
    return outputs


def plot_volcano(assoc: pd.DataFrame, path_base: Path) -> list[str]:
    x = assoc["spearman_rho"].to_numpy(dtype=float)
    y = -np.log10(np.clip(assoc["spearman_FDR"].fillna(1).to_numpy(dtype=float), 1e-300, 1))
    fig, ax = plt.subplots(figsize=(5.4, 4.0))
    ax.scatter(x, y, s=18, c=np.where(assoc["spearman_FDR"].fillna(1) < 0.05, "#D55E00", "#999999"), alpha=0.85)
    top = assoc.sort_values("spearman_rho", ascending=False).head(8)
    for row in top.itertuples(index=False):
        ax.text(row.spearman_rho, -math.log10(max(float(row.spearman_FDR), 1e-300)), row.TF, fontsize=7)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_xlabel("Spearman rho with CellRank CNV fate probability")
    ax.set_ylabel("-log10(BH-FDR)")
    ax.set_title("Canonical SCENIC regulon-fate association")
    fig.tight_layout()
    return save_figure(fig, path_base)


def plot_phase_heatmap(assoc: pd.DataFrame, path_base: Path, top_n: int) -> list[str]:
    top = assoc.sort_values("spearman_rho", ascending=False).head(top_n)
    matrix = top.set_index("TF")[["early_mean_AUC", "middle_mean_AUC", "late_mean_AUC"]]
    fig, ax = plt.subplots(figsize=(4.8, max(3.5, 0.22 * matrix.shape[0])))
    im = ax.imshow(matrix.to_numpy(dtype=float), aspect="auto", cmap="magma")
    ax.set_yticks(range(matrix.shape[0]))
    ax.set_yticklabels(matrix.index)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["early", "middle", "late"])
    ax.set_title("Canonical regulon activity by trajectory phase")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="Mean AUC")
    fig.tight_layout()
    return save_figure(fig, path_base)


def plot_pseudotime(auc: pd.DataFrame, cells: pd.DataFrame, assoc: pd.DataFrame, path_base: Path) -> list[str]:
    wanted = [tf for tf in ["HNF4A", "PPARA", "JUN", "JUND", "FOS", "CEBPB", "EGR1", "SOX4"] if tf in set(assoc["TF"])]
    if not wanted:
        return []
    ncols = 4
    nrows = int(math.ceil(len(wanted) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(11.0, 2.7 * nrows), squeeze=False)
    pseudotime = pd.to_numeric(cells["driver_main_strict__pseudotime_median"], errors="coerce")
    for ax, tf in zip(axes.flat, wanted):
        regulon = f"{tf}(+)"
        mask = pseudotime.notna() & auc[regulon].notna()
        x = pseudotime.loc[mask].to_numpy(dtype=float)
        y = auc.loc[mask, regulon].to_numpy(dtype=float)
        ax.scatter(x, y, s=2, alpha=0.12, color="#0072B2")
        if len(x) > 10:
            smoothed = lowess(y, x, frac=0.25, return_sorted=True)
            ax.plot(smoothed[:, 0], smoothed[:, 1], color="#D55E00", lw=1.3)
        ax.set_title(tf)
        ax.set_xlabel("pseudotime")
        ax.set_ylabel("AUC")
    for ax in axes.flat[len(wanted) :]:
        ax.axis("off")
    fig.suptitle("Canonical regulon activity across pseudotime", y=1.01)
    fig.tight_layout()
    return save_figure(fig, path_base)


def plot_cnv_heatmap(assoc: pd.DataFrame, path_base: Path, top_n: int) -> list[str]:
    top = assoc.reindex(assoc["CNV_difference_AUC"].abs().sort_values(ascending=False).head(top_n).index)
    matrix = top.set_index("TF")[["CNV_negative_mean_AUC", "CNV_positive_mean_AUC"]]
    fig, ax = plt.subplots(figsize=(4.2, max(3.5, 0.22 * matrix.shape[0])))
    im = ax.imshow(matrix.to_numpy(dtype=float), aspect="auto", cmap="coolwarm")
    ax.set_yticks(range(matrix.shape[0]))
    ax.set_yticklabels(matrix.index)
    ax.set_xticks(range(2))
    ax.set_xticklabels(["CNV-negative", "CNV-positive"], rotation=25, ha="right")
    ax.set_title("Canonical regulon activity by CNV state")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="Mean AUC")
    fig.tight_layout()
    return save_figure(fig, path_base)


def plot_sample_robustness(loo: pd.DataFrame, assoc: pd.DataFrame, path_base: Path) -> list[str]:
    if loo.empty:
        return []
    top = assoc.sort_values("spearman_rho", ascending=False).head(20)["regulon"].tolist()
    matrix = loo.loc[loo["regulon"].isin(top)].pivot(index="regulon", columns="omitted_group", values="rho")
    matrix = matrix.reindex(top).dropna(how="all")
    fig, ax = plt.subplots(figsize=(max(6, 0.45 * matrix.shape[1]), max(3.5, 0.25 * matrix.shape[0])))
    im = ax.imshow(matrix.to_numpy(dtype=float), aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_yticks(range(matrix.shape[0]))
    ax.set_yticklabels([x.replace("(+)", "") for x in matrix.index])
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=70, ha="right")
    ax.set_title("Leave-one-sample-out Spearman rho")
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02, label="rho")
    fig.tight_layout()
    return save_figure(fig, path_base)


def plot_comparison(comparison: pd.DataFrame, path_base: Path) -> list[str]:
    if comparison.empty:
        return []
    sub = comparison.loc[comparison["canonical_present"] | comparison["old_present"]].copy()
    sub = sub.sort_values(["canonical_present", "canonical_rank"], ascending=[False, True]).head(40)
    values = np.column_stack([sub["old_rank"].fillna(0), sub["canonical_rank"].fillna(0)])
    fig, ax = plt.subplots(figsize=(5.5, max(4, 0.2 * sub.shape[0])))
    im = ax.imshow(values, aspect="auto", cmap="viridis")
    ax.set_yticks(range(sub.shape[0]))
    ax.set_yticklabels(sub["TF"])
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["old 6.3 rank", "canonical 6.3b rank"])
    ax.set_title("Exploratory versus canonical SCENIC TF rank")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="rank")
    fig.tight_layout()
    return save_figure(fig, path_base)


def plot_umap(auc: pd.DataFrame, assoc: pd.DataFrame, trajectory_path: Path, path_base: Path) -> list[str]:
    if not trajectory_path.exists():
        return []
    trajectory = ad.read_h5ad(trajectory_path)
    if "X_umap" not in trajectory.obsm:
        return []
    common = pd.Index(trajectory.obs_names.astype(str)).intersection(auc.index)
    if len(common) < 100:
        return []
    coords = np.asarray(trajectory[common].obsm["X_umap"])
    wanted = [tf for tf in ["HNF4A", "JUN", "FOS", "CEBPB", "EGR1", "SOX4"] if f"{tf}(+)" in auc.columns]
    if not wanted:
        return []
    ncols = 3
    nrows = int(math.ceil(len(wanted) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(10, 3.1 * nrows), squeeze=False)
    for ax, tf in zip(axes.flat, wanted):
        values = auc.loc[common, f"{tf}(+)"].to_numpy(dtype=float)
        sc = ax.scatter(coords[:, 0], coords[:, 1], c=values, s=2, cmap="magma", linewidths=0)
        ax.set_title(tf)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
    for ax in axes.flat[len(wanted) :]:
        ax.axis("off")
    fig.suptitle("Canonical regulon activity on trajectory UMAP", y=1.01)
    fig.tight_layout()
    trajectory.file.close() if getattr(trajectory, "isbacked", False) else None
    return save_figure(fig, path_base)


def main() -> None:
    start = time.time()
    args = parse_args()
    for path in [args.input_h5ad, args.ctx_output, args.auc, args.old_summary]:
        if not path.exists():
            raise FileNotFoundError(path)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    args.processed_dir.mkdir(parents=True, exist_ok=True)
    args.figures_dir.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(args.input_h5ad)
    cells = adata.obs.copy()
    cells.index = adata.obs_names.astype(str)
    auc = read_auc(args.auc, pd.Index(adata.obs_names.astype(str)))
    if args.fate_key not in cells.columns:
        raise KeyError(f"Missing fate key: {args.fate_key}")
    fate_mask = pd.to_numeric(cells[args.fate_key], errors="coerce").notna()
    ctx = read_ctx(args.ctx_output)
    summary = pd.read_csv(args.metadata_dir / "driver_module6_3b_canonical_regulons_seed777.csv")
    summary["TF"] = summary["TF"].astype(str)
    summary["regulon"] = summary["regulon"].astype(str)
    associations = association_table(auc, cells, args.fate_key, args.pseudotime_key, args.phase_key, args.cnv_key)
    associations = associations.merge(summary[["TF", "regulon", "regulon_size", "motif_NES_max"]], on=["TF", "regulon"], how="left", suffixes=("", "_summary"))
    if "regulon_size_summary" in associations.columns:
        summary_sizes = pd.to_numeric(associations["regulon_size_summary"], errors="coerce")
        association_sizes = pd.to_numeric(associations["regulon_size"], errors="coerce")
        associations["regulon_size"] = summary_sizes.fillna(association_sizes)
        associations = associations.drop(columns=["regulon_size_summary"])
    associations = associations.merge(summary[["regulon", "primary_analysis_eligible"]], on="regulon", how="left")
    pseudobulk = sample_pseudobulk(auc, cells, args.fate_key, args.sample_key, args.dataset_key, args.phase_key, args.cnv_key)
    sample_pb = pseudobulk_correlations(pseudobulk, list(auc.columns), args.sample_key)
    associations = associations.merge(sample_pb, on="regulon", how="left")
    sample_loo = leave_one_out(auc, cells, args.fate_key, args.sample_key, "spearman")
    dataset_loo = leave_one_out(auc, cells, args.fate_key, args.dataset_key, "spearman")
    robustness = robustness_table(associations, sample_loo, dataset_loo, sample_pb)
    associations = associations.merge(robustness, on="regulon", how="left")
    associations["biological_program_annotation"] = associations["TF"].map(biological_program)
    associations["rank"] = associations["spearman_rho"].rank(method="min", ascending=False).astype("Int64")
    associations = associations.sort_values(["spearman_rho", "pseudotime_spearman_rho"], ascending=[False, False]).reset_index(drop=True)
    axis_table = evaluate_axes(associations)
    grn_qc = {}
    grn_qc_path = args.metadata_dir / "driver_module6_3b_grn_qc_report.json"
    if grn_qc_path.exists():
        grn_qc = json.loads(grn_qc_path.read_text(encoding="utf-8"))
    resource_validation_path = ROOT / "metadata/driver/scenic_resources_v10/resource_validation.json"
    resource_validation = json.loads(resource_validation_path.read_text(encoding="utf-8")) if resource_validation_path.exists() else {}
    ranking_validation = resource_validation.get("resources", {}).get("ranking_10kb", {})
    tf_list_path = ROOT / "data/processed/driver/scenic_module6_3b_formal/driver_union_tfs_in_matrix.txt"
    n_input_tfs = len([line for line in tf_list_path.read_text(encoding="utf-8").splitlines() if line.strip()]) if tf_list_path.exists() else np.nan

    association_path = args.metadata_dir / "driver_module6_3b_regulon_cellrank_fate_association.tsv"
    pseudotime_path = args.metadata_dir / "driver_module6_3b_regulon_pseudotime_association.tsv"
    pseudobulk_path = args.metadata_dir / "driver_module6_3b_sample_pseudobulk_regulon_activity.tsv.gz"
    sample_loo_path = args.metadata_dir / "driver_module6_3b_regulon_leave_one_sample_out.tsv"
    dataset_loo_path = args.metadata_dir / "driver_module6_3b_regulon_leave_one_dataset_out.tsv"
    robustness_path = args.metadata_dir / "driver_module6_3b_regulon_robustness.tsv"
    axis_path = args.metadata_dir / "driver_module6_3b_three_axis_evaluation.tsv"
    auc_path = args.metadata_dir / "driver_module6_3b_canonical_regulon_auc_9512.tsv.gz"
    associations.to_csv(association_path, sep="\t", index=False)
    associations[[column for column in associations.columns if "pseudotime" in column or column in {"regulon", "TF", "regulon_size", "early_mean_AUC", "middle_mean_AUC", "late_mean_AUC"}]].to_csv(pseudotime_path, sep="\t", index=False)
    pseudobulk.to_csv(pseudobulk_path, sep="\t", index=False, compression="gzip")
    sample_loo.to_csv(sample_loo_path, sep="\t", index=False)
    dataset_loo.to_csv(dataset_loo_path, sep="\t", index=False)
    robustness.to_csv(robustness_path, sep="\t", index=False)
    axis_table.to_csv(axis_path, sep="\t", index=False)
    auc_out = auc.copy()
    auc_out.insert(0, "cell_id", auc_out.index.astype(str))
    auc_out.to_csv(auc_path, sep="\t", index=False, compression="gzip")

    comparison = pd.DataFrame()
    if args.old_summary.exists():
        old = pd.read_csv(args.old_summary, sep="\t")
        old["TF"] = old["tf"].astype(str)
        old["old_rank"] = old["cnv_fate_spearman_rho"].rank(method="min", ascending=False) if "cnv_fate_spearman_rho" in old.columns else np.nan
        old = old.rename(columns={"regulon": "old_regulon", "cnv_fate_spearman_rho": "old_spearman_rho"})
        new = associations[["TF", "regulon", "spearman_rho", "rank"]].rename(columns={"spearman_rho": "canonical_spearman_rho", "rank": "canonical_rank"})
        comparison = old[["TF", "old_regulon", "old_spearman_rho", "old_rank"]].merge(new, on="TF", how="outer", indicator=True)
        comparison["old_present"] = comparison["_merge"].isin(["both", "left_only"])
        comparison["canonical_present"] = comparison["_merge"].isin(["both", "right_only"])
        comparison["comparison_status"] = comparison["_merge"].map({"both": "shared_tf", "left_only": "old_only", "right_only": "canonical_only"})
        comparison = comparison.drop(columns=["_merge"])
    comparison_path = args.metadata_dir / "driver_module6_3_vs_6_3b_comparison.tsv"
    comparison.to_csv(comparison_path, sep="\t", index=False)

    master = associations.copy()
    top_motif = ctx.assign(NES=pd.to_numeric(ctx["NES"], errors="coerce")).sort_values("NES", ascending=False).drop_duplicates("TF")
    master = master.merge(top_motif[["TF", "MotifID"]].rename(columns={"MotifID": "motif"}), on="TF", how="left")
    master["motif_NES"] = master["motif_NES_max"]
    master["evidence_tier"] = master["evidence_tier"].fillna("Tier C")
    master.loc[master["regulon_size"].lt(10), "evidence_tier"] = "Tier C_small_regulon"
    master_path = args.metadata_dir / "driver_module6_3b_canonical_scenic_master_table.tsv"
    master.to_csv(master_path, sep="\t", index=False)

    adata.obsm["module6_3b_canonical_regulon_auc"] = auc.to_numpy(dtype=np.float32)
    adata.uns["module6_3b_canonical_regulon_auc_names"] = auc.columns.tolist()
    adata.uns["module6_3b_canonical_scenic"] = {
        "method": "GRNBoost2 -> cisTarget ctx -> AUCell",
        "n_cells": int(adata.n_obs),
        "n_regulons": int(auc.shape[1]),
        "n_input_tfs": int(n_input_tfs) if pd.notna(n_input_tfs) else None,
        "n_grn_edges": grn_qc.get("n_edges"),
        "n_cellrank_fate_cells": int(cells[args.fate_key].notna().sum()),
        "seed": 777,
    }
    output_h5ad = args.processed_dir / "driver_canonical_scenic.module6_3b.h5ad"
    adata.write_h5ad(output_h5ad, compression="gzip")

    figure_outputs = []
    figure_outputs += plot_volcano(associations, args.figures_dir / "driver_module6_3b_regulon_cellrank_fate_volcano")
    figure_outputs += plot_phase_heatmap(associations, args.figures_dir / "driver_module6_3b_regulon_phase_heatmap", args.top_n)
    figure_outputs += plot_pseudotime(auc, cells, associations, args.figures_dir / "driver_module6_3b_regulon_pseudotime_loess")
    figure_outputs += plot_cnv_heatmap(associations, args.figures_dir / "driver_module6_3b_regulon_cnv_state_heatmap", args.top_n)
    figure_outputs += plot_sample_robustness(sample_loo, associations, args.figures_dir / "driver_module6_3b_sample_robustness_heatmap")
    figure_outputs += plot_comparison(comparison, args.figures_dir / "driver_module6_3b_old_vs_canonical_comparison")
    figure_outputs += plot_umap(auc, associations, args.trajectory_h5ad, args.figures_dir / "driver_module6_3b_top_regulons_umap")

    top_positive = associations.sort_values("spearman_rho", ascending=False).head(15)
    top_negative = associations.sort_values("spearman_rho", ascending=True).head(15)
    status_path = ROOT / "metadata/driver/driver_module6_3b_canonical_scenic_status.json"
    status = {
        "module": "6.3b",
        "status": "MODULE6_3B_COMPLETE",
        "checkpoint_history": [
            {"stage": "INPUT_QC", "status": "INPUT_QC_COMPLETE"},
            {"stage": "RESOURCES", "status": "RESOURCES_COMPLETE"},
            {"stage": "GRNBoost2", "status": "GRN_COMPLETE"},
            {"stage": "GRN_QC", "status": "GRN_QC_COMPLETE"},
            {"stage": "CTX", "status": "CTX_COMPLETE"},
            {"stage": "CANONICAL_REGULON_QC", "status": "CTX_QC_COMPLETE_WITH_WARNINGS"},
            {"stage": "AUCELL", "status": "AUCELL_COMPLETE"},
            {"stage": "CELLRANK_ASSOCIATION", "status": "ASSOCIATION_COMPLETE"},
            {"stage": "ROBUSTNESS", "status": "ROBUSTNESS_COMPLETE"},
            {"stage": "REPORTING", "status": "MODULE6_3B_COMPLETE"},
        ],
        "completed": {
            "input_qc": record_path(ROOT / "metadata/driver/scenic_module6_3b_input_qc.tsv"),
            "grnboost2": record_path(ROOT / "metadata/driver/scenic_module6_3b/driver_module6_3b_grnboost2_seed777_adjacencies.tsv.gz"),
            "ctx": record_path(args.ctx_output),
            "aucell": record_path(args.auc),
            "cellrank_association": record_path(association_path),
            "robustness": record_path(robustness_path),
            "three_axis_evaluation": record_path(axis_path),
        },
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_regulons": int(auc.shape[1]),
        "n_input_tfs": int(n_input_tfs) if pd.notna(n_input_tfs) else None,
        "n_grn_edges": grn_qc.get("n_edges"),
        "n_primary_analysis_regulons": int(summary["primary_analysis_eligible"].sum()),
        "n_small_regulons_lt_10_targets": int((summary["regulon_size"] < 10).sum()),
        "n_cellrank_fate_cells": int(cells[args.fate_key].notna().sum()),
        "n_cellrank_fate_datasets": int(cells.loc[fate_mask, args.dataset_key].astype(str).nunique()),
        "n_cellrank_fate_samples": int(cells.loc[fate_mask, args.sample_key].astype(str).nunique()),
        "n_unknown_metadata_cells_outside_fate_subset": int((~fate_mask).sum()),
        "top_positive_regulons": top_positive[["TF", "regulon", "spearman_rho", "spearman_FDR"]].to_dict(orient="records"),
        "top_negative_regulons": top_negative[["TF", "regulon", "spearman_rho", "spearman_FDR"]].to_dict(orient="records"),
        "limitation": "Single GRNBoost2 seed was run; seed stability is marked not_run_single_seed_777 in robustness output.",
        "publication_readiness": {
            "extended_data_and_hypothesis_prioritization": "ready_after_manual_figure_review",
            "main_figure_computational_regulatory_activity": "ready_with_explicit_association_language",
            "main_figure_causal_mechanism": "not_ready_without_orthogonal_perturbation",
            "remaining_computational_gap": "multi_seed_GRNBoost2_stability",
        },
        "outputs": {"master_table": record_path(master_path), "three_axis_evaluation": record_path(axis_path), "figures": figure_outputs, "final_report": record_path(ROOT / "reports/module6_3b_canonical_scenic_final_report.md")},
    }
    status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")

    report_lines = [
        "# Module 6.3b canonical SCENIC final report",
        "",
        "## 1. Objective",
        "",
        "Formalize the regulatory analysis on the 9,512-cell full-expression driver union using GRNBoost2, matching hg38 mc_v10_clust cisTarget resources, canonical motif pruning, AUCell, CellRank association and sample-aware robustness.",
        "",
        "## 2. Input",
        "",
        f"- Cells: `{adata.n_obs}`; genes after all-zero cleanup: `{adata.n_vars}`; canonical regulons scored: `{len(auc.columns)}`.",
        f"- Primary regulons with at least 10 targets: `{int(summary['primary_analysis_eligible'].sum())}`; small regulons retained for supplementary review: `{int((summary['regulon_size'] < 10).sum())}`.",
        f"- CellRank intersection: `{cells[args.fate_key].notna().sum()}` cells.",
        f"- Known datasets in the CellRank intersection: `{cells.loc[fate_mask, args.dataset_key].astype(str).nunique()}`; known samples: `{cells.loc[fate_mask, args.sample_key].astype(str).nunique()}`; the remaining `{int((~fate_mask).sum())}` AUCell cells are outside the CellRank subset and have Unknown metadata by design.",
        "",
        "## 3. Resources",
        "",
        "- Genome: hg38, RefSeq 80.",
        "- Ranking family: mc_v10_clust, with 10 kb around TSS and proximal promoter resources.",
        "- Motif annotation: motifs-v10nr_clust-nr HGNC.",
        f"- 10 kb ranking coverage for formal genes: `{ranking_validation.get('input_gene_coverage', 'NA')}/{adata.n_vars}`; genes outside the RefSeq 80 ranking universe: `{ranking_validation.get('input_gene_missing', 'NA')}`.",
        "",
        "## 4. GRN and cisTarget",
        "",
        f"- GRN method: GRNBoost2, seed 777, full-expression input; edges=`{grn_qc.get('n_edges', 'NA')}`, represented TFs=`{grn_qc.get('n_tfs', n_input_tfs)}`, represented targets=`{grn_qc.get('n_targets', adata.n_vars)}`.",
        f"- Canonical regulons scored: `{auc.shape[1]}`.",
        "- Historical Module 6.3 adjacency and 6.3c results are excluded from formal 6.3b inference.",
        "",
        "## 5. AUCell and CellRank association",
        "",
        f"- AUCell coverage: `{auc.shape[0]} cells x {auc.shape[1]} regulons`.",
        f"- Fate association uses only the `{cells[args.fate_key].notna().sum()}` non-null CellRank fate cells.",
        f"- Top positive regulons: `{', '.join(top_positive['TF'].head(10))}`.",
        f"- Top negative regulons: `{', '.join(top_negative['TF'].head(10))}`.",
        "",
        "## 6. Pseudotime and robustness",
        "",
        "- Spearman and Pearson associations use BH-FDR.",
        "- Pseudotime associations include LOESS visualization and early/middle/late phase summaries.",
        "- Sample pseudobulk, leave-one-sample-out and leave-one-dataset-out tables are generated.",
        "- Seed stability remains a declared gap because only seed 777 was run; the limitation is explicitly recorded in the robustness table.",
        "",
        "## 7. Old versus canonical comparison",
        "",
        "The comparison table records shared TFs, old-only hits, canonical-only hits and rank changes. Old Module 6.3 is labeled exploratory co-expression without motif pruning; Module 6.3b is the formal canonical branch.",
        "",
        "## 8. Three-axis assessment",
        "",
        "The three-axis assessment was added after the unbiased canonical analysis. Missing canonical regulons and directionally inconsistent layers are reported explicitly in the axis table below.",
        "",
        "## 9. Limitations and publication readiness",
        "",
        "- The data are observational single-cell measurements; SCENIC and CellRank support regulatory activity and fate association rather than direct causality.",
        "- CellRank fate probabilities are computational inferences and remain sensitive to terminal-state definitions and sample composition.",
        "- GRNBoost2 stochasticity is incompletely assessed with one seed.",
        "- Motif databases encode prior sequence/regulatory knowledge and do not establish TF binding in these samples.",
        "- Canonical computational results are suitable for Extended Data and hypothesis prioritization. Main-figure causal claims require orthogonal perturbation or external binding evidence.",
        "- Publication readiness: canonical computational activity is suitable for Extended Data and association-based main-figure panels; causal mechanism is not ready without orthogonal perturbation, and multi-seed GRN stability remains an open computational check.",
        "",
        "## 10. Output files",
        "",
    ]
    axis_lines = [
        "",
        "### Axis-level results",
        "",
    ]
    for axis, sub in axis_table.groupby("axis", sort=False):
        detected = sub.loc[sub["canonical_regulon_detected"], "TF"].astype(str).tolist()
        absent = sub.loc[~sub["canonical_regulon_detected"], "TF"].astype(str).tolist()
        status_value = str(sub["axis_status"].iloc[0])
        detected_metrics = []
        for row in sub.loc[sub["canonical_regulon_detected"]].itertuples(index=False):
            fate_rho = "NA" if pd.isna(row.spearman_rho) else f"{row.spearman_rho:.3f}"
            time_rho = "NA" if pd.isna(row.pseudotime_spearman_rho) else f"{row.pseudotime_spearman_rho:.3f}"
            detected_metrics.append(
                f"{row.TF}(fate_rho={fate_rho},pseudotime_rho={time_rho},sample_LOSO={row.sample_loo_stable},dataset_LODO={row.dataset_loo_stable},tier={row.evidence_tier})"
            )
        axis_lines.append(
            f"- `{axis}`: status=`{status_value}`; detected=`{','.join(detected) if detected else 'none'}`; not detected=`{','.join(absent) if absent else 'none'}`; metrics=`{'; '.join(detected_metrics) if detected_metrics else 'none'}`."
        )
    insert_at = report_lines.index("## 9. Limitations and publication readiness")
    report_lines[insert_at:insert_at] = axis_lines
    for label, path in {
        "association": association_path,
        "pseudotime": pseudotime_path,
        "pseudobulk": pseudobulk_path,
        "sample_LOSO": sample_loo_path,
        "dataset_LODO": dataset_loo_path,
        "robustness": robustness_path,
        "three_axis_evaluation": axis_path,
        "AUC": auc_path,
        "master_table": master_path,
        "integrated_h5ad": output_h5ad,
        "grnboost2_report": grn_qc_path,
        "resource_validation": resource_validation_path,
        "environment": ROOT / "metadata/driver/scenic_module6_3b_environment.txt",
        "status": status_path,
    }.items():
        report_lines.append(f"- `{label}`: `{record_path(path)}`")
    (ROOT / "reports/module6_3b_canonical_scenic_final_report.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    report = {
        "module": "6.3b",
        "status": "MODULE6_3B_COMPLETE",
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_regulons": int(auc.shape[1]),
        "n_primary_analysis_regulons": int(summary["primary_analysis_eligible"].sum()),
        "n_small_regulons_lt_10_targets": int((summary["regulon_size"] < 10).sum()),
        "n_cellrank_fate_cells": int(cells[args.fate_key].notna().sum()),
        "n_cellrank_fate_datasets": int(cells.loc[fate_mask, args.dataset_key].astype(str).nunique()),
        "n_cellrank_fate_samples": int(cells.loc[fate_mask, args.sample_key].astype(str).nunique()),
        "n_unknown_metadata_cells_outside_fate_subset": int((~fate_mask).sum()),
        "outputs": {
            "master_table": record_path(master_path),
            "three_axis_evaluation": record_path(axis_path),
            "report": record_path(ROOT / "reports/module6_3b_canonical_scenic_final_report.md"),
            "grnboost2_report": record_path(grn_qc_path),
            "resource_validation": record_path(resource_validation_path),
            "environment": record_path(ROOT / "metadata/driver/scenic_module6_3b_environment.txt"),
            "figures": figure_outputs,
        },
        "package_versions": {name: package_version(name) for name in ["anndata", "pandas", "numpy", "scipy", "pyarrow", "pyscenic", "ctxcore", "matplotlib", "statsmodels"]},
        "elapsed_seconds": round(time.time() - start, 3),
    }
    (args.metadata_dir / "driver_module6_3b_final_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
