from __future__ import annotations

import argparse
import json
import math
import time
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Iterable

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import fisher_exact, pearsonr, t as t_dist


ROOT = Path(__file__).resolve().parents[1]

MODULE = "6.4"
BASE_GRN_URL = (
    "https://raw.githubusercontent.com/morris-lab/CellOracle/master/"
    "celloracle/data/promoter_base_GRN/"
    "hg38_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet"
)

PREFERRED_PRO_TFS = [
    "JUN",
    "FOS",
    "JUNB",
    "JUND",
    "ATF3",
    "IRF1",
    "MYC",
    "CEBPB",
    "KLF2",
    "EGR1",
    "SOX4",
    "MAFB",
]
PREFERRED_MAINTENANCE_TFS = ["HNF4A", "PPARA", "HLF"]
RESERVE_TFS = ["MAFF", "CEBPD", "KLF4", "STAT1", "IRF8", "SPI1", "PPARGC1A", "AR"]
MAINTENANCE_TFS = set(PREFERRED_MAINTENANCE_TFS + ["PPARGC1A", "AR"])
NONCANONICAL_EXCLUSIONS = {"TPI1"}
IMMUNE_RISK_TFS = {"SPI1", "IRF8"}
INTERFERON_CONTEXT_TFS = {"STAT1", "IRF1"}

BIOLOGY_SCORES = {
    "JUN": 5.0,
    "FOS": 5.0,
    "JUNB": 5.0,
    "JUND": 4.5,
    "ATF3": 4.5,
    "IRF1": 4.0,
    "MYC": 5.0,
    "CEBPB": 4.0,
    "KLF2": 3.5,
    "EGR1": 4.0,
    "SOX4": 4.0,
    "MAFB": 3.5,
    "HNF4A": 5.0,
    "PPARA": 4.5,
    "HLF": 4.0,
    "MAFF": 3.5,
    "CEBPD": 3.5,
    "KLF4": 3.5,
    "STAT1": 3.0,
    "IRF8": 2.5,
    "SPI1": 2.0,
    "PPARGC1A": 4.0,
    "AR": 3.5,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 6.4: select CellOracle input TFs for CNV malignant fate.")
    parser.add_argument(
        "--input-h5ad",
        type=Path,
        default=ROOT / "data/processed/driver/driver_cistarget_regulon_activity.module6_3c.h5ad",
    )
    parser.add_argument(
        "--regulon-summary",
        type=Path,
        default=ROOT / "metadata/driver/driver_module6_3c_cistarget_regulon_summary.tsv",
    )
    parser.add_argument(
        "--auc",
        type=Path,
        default=ROOT / "metadata/driver/driver_module6_3c_cistarget_regulon_auc.tsv.gz",
    )
    parser.add_argument(
        "--cellrank-drivers",
        type=Path,
        default=ROOT / "metadata/driver/driver_module6_2_cellrank_lineage_drivers.tsv.gz",
    )
    parser.add_argument(
        "--tf-catalog",
        type=Path,
        default=ROOT / "metadata/driver/scenic_resources/allTFs_hg38.txt",
    )
    parser.add_argument(
        "--celloracle-base-grn",
        type=Path,
        default=ROOT / "metadata/driver/scenic_resources/celloracle_hg38_promoter_base_grn.parquet",
    )
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/driver")
    parser.add_argument("--figures-dir", type=Path, default=ROOT / "figures/driver")
    parser.add_argument("--fate-key", default="cellrank_fate_prob_cnv_supported_malignant")
    parser.add_argument("--time-key", default="driver_main_strict__pseudotime_median")
    parser.add_argument("--phase-key", default="driver_main_strict__pseudotime_phase")
    parser.add_argument("--dataset-key", default="dataset")
    parser.add_argument("--sample-key", default="sample_id")
    parser.add_argument("--top-n-cellrank-drivers", type=int, default=200)
    parser.add_argument("--min-detection-rate", type=float, default=0.05)
    parser.add_argument("--min-detected-datasets", type=int, default=3)
    parser.add_argument("--min-base-grn-links", type=int, default=10)
    parser.add_argument("--min-group-cells", type=int, default=30)
    parser.add_argument("--main-panel-size", type=int, default=15)
    return parser.parse_args()


def configure_plot_style() -> None:
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


def package_version(name: str) -> str:
    try:
        return version(name)
    except PackageNotFoundError:
        return "not_installed"


def write_dataframe(path: Path, df: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = "gzip" if path.suffix == ".gz" else None
    df.to_csv(path, sep="\t", index=False, compression=compression)


def bh_qvalues(pvalues: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvalues, errors="coerce")
    q = pd.Series(np.nan, index=p.index, dtype=float)
    valid = p.dropna()
    if valid.empty:
        return q
    order = valid.sort_values().index
    ranked = valid.loc[order].to_numpy(dtype=float)
    m = float(len(ranked))
    adjusted = ranked * m / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    q.loc[order] = np.clip(adjusted, 0.0, 1.0)
    return q


def clip01(value: object) -> float:
    out = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(out) or not np.isfinite(out):
        return 0.0
    return float(np.clip(out, 0.0, 1.0))


def standardize(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    sd = numeric.std(ddof=0)
    if pd.isna(sd) or sd == 0:
        return pd.Series(np.nan, index=values.index, dtype=float)
    return (numeric - numeric.mean()) / sd


def parse_target_genes(value: object) -> set[str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return set()
    return {part.strip().upper() for part in str(value).replace(",", ";").split(";") if part.strip()}


def read_tf_catalog(path: Path) -> set[str]:
    if not path.exists():
        raise FileNotFoundError(path)
    values = []
    for line in path.read_text(encoding="utf-8").splitlines():
        item = line.strip()
        if item and not item.startswith("#"):
            values.append(item.upper())
    if not values:
        raise ValueError(f"TF catalog is empty: {path}")
    return set(values)


def read_base_grn(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix.lower() == ".parquet":
        return pd.read_parquet(path)
    if path.suffix.lower() in {".csv"}:
        return pd.read_csv(path)
    return pd.read_csv(path, sep="\t")


def read_auc_matrix(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    sep = "\t" if path.suffix == ".gz" or path.suffix == ".tsv" else ","
    auc = pd.read_csv(path, sep=sep)
    if auc.shape[1] == 1:
        auc = pd.read_csv(path)
    if auc.empty:
        raise ValueError(f"AUCell matrix is empty: {path}")
    return auc


def align_auc_to_cells(auc: pd.DataFrame, cell_ids: pd.Index) -> pd.DataFrame:
    work = auc.copy()
    cell_col = None
    for candidate in ["cell_id", "CellID", "Cell", "cells"]:
        if candidate in work.columns:
            cell_col = candidate
            break
    if cell_col is not None:
        work.index = work.pop(cell_col).astype(str)
    else:
        work.index = work.index.astype(str)
    missing = pd.Index(cell_ids.astype(str)).difference(pd.Index(work.index.astype(str)))
    if len(missing) > 0:
        raise ValueError(f"{len(missing)} h5ad cells are missing from AUCell matrix.")
    aligned = work.loc[cell_ids.astype(str)].apply(pd.to_numeric, errors="coerce")
    aligned = aligned.loc[:, aligned.nunique(dropna=True) > 1]
    if aligned.empty:
        raise ValueError("No variable AUCell regulon columns remain after alignment.")
    return aligned


def compute_base_grn_compatibility(base_grn: pd.DataFrame, candidate_tfs: Iterable[str], min_links: int = 10) -> pd.DataFrame:
    special = {"PEAK_ID", "GENE_SHORT_NAME"}
    col_map = {str(col).upper(): col for col in base_grn.columns if str(col).upper() not in special}
    gene_col = "gene_short_name" if "gene_short_name" in base_grn.columns else None
    rows = []
    for tf in candidate_tfs:
        tf_upper = str(tf).upper()
        source_col = col_map.get(tf_upper)
        if source_col is None:
            rows.append(
                {
                    "tf": tf_upper,
                    "tf_in_base_grn": False,
                    "base_grn_tf_name": "",
                    "base_grn_outgoing_links": 0,
                    "base_grn_target_genes": 0,
                    "base_grn_has_min_links": False,
                }
            )
            continue
        values = pd.to_numeric(base_grn[source_col], errors="coerce").fillna(0)
        mask = values > 0
        target_genes = int(base_grn.loc[mask, gene_col].astype(str).str.upper().nunique()) if gene_col else int(mask.sum())
        rows.append(
            {
                "tf": tf_upper,
                "tf_in_base_grn": True,
                "base_grn_tf_name": str(source_col),
                "base_grn_outgoing_links": int(mask.sum()),
                "base_grn_target_genes": target_genes,
                "base_grn_has_min_links": bool(mask.sum() >= min_links),
            }
        )
    return pd.DataFrame(rows)


def compute_expression_summary(
    adata: ad.AnnData,
    candidate_tfs: Iterable[str],
    fate_key: str,
    dataset_key: str,
    dataset_detection_min_rate: float = 0.01,
) -> pd.DataFrame:
    candidate_tfs = [str(tf).upper() for tf in candidate_tfs]
    var_map = {str(gene).upper(): gene for gene in adata.var_names.astype(str)}
    present = [tf for tf in candidate_tfs if tf in var_map]
    cells = adata.obs.copy()
    cells.index = adata.obs_names.astype(str)
    main_mask = cells[fate_key].notna().to_numpy() if fate_key in cells.columns else np.ones(adata.n_obs, dtype=bool)

    expression_vectors: dict[str, np.ndarray] = {}
    if present:
        subset = adata[:, [var_map[tf] for tf in present]]
        matrix = subset.layers["counts"] if "counts" in subset.layers else subset.X
        if sparse.issparse(matrix):
            matrix = matrix.tocsr()
            for i, tf in enumerate(present):
                expression_vectors[tf] = np.asarray(matrix[:, i].toarray()).ravel()
        else:
            array = np.asarray(matrix)
            for i, tf in enumerate(present):
                expression_vectors[tf] = np.asarray(array[:, i]).ravel()

    rows = []
    for tf in candidate_tfs:
        vec = expression_vectors.get(tf)
        if vec is None:
            rows.append(
                {
                    "tf": tf,
                    "tf_in_expression": False,
                    "detection_rate_main": 0.0,
                    "mean_expression_main": 0.0,
                    "detected_dataset_count": 0,
                }
            )
            continue
        main_values = vec[main_mask]
        detection_rate = float((main_values > 0).mean()) if main_values.size else 0.0
        mean_expr = float(np.mean(main_values)) if main_values.size else 0.0
        detected_datasets = 0
        if dataset_key in cells.columns:
            for _, idx in cells.loc[main_mask].groupby(dataset_key, observed=True).groups.items():
                positions = cells.index.get_indexer(pd.Index(idx))
                positions = positions[positions >= 0]
                if positions.size and float((vec[positions] > 0).mean()) >= dataset_detection_min_rate:
                    detected_datasets += 1
        rows.append(
            {
                "tf": tf,
                "tf_in_expression": True,
                "detection_rate_main": detection_rate,
                "mean_expression_main": mean_expr,
                "detected_dataset_count": int(detected_datasets),
            }
        )
    return pd.DataFrame(rows)


def compute_cellrank_target_overlap(
    candidates: pd.DataFrame,
    drivers: pd.DataFrame,
    universe: pd.Index,
    top_n_drivers: int = 200,
    lineage: str = "cnv_supported_malignant",
) -> pd.DataFrame:
    work = drivers.copy()
    if "lineage" in work.columns:
        work = work.loc[work["lineage"].astype(str).eq(lineage)].copy()
    work["gene_upper"] = work["gene"].astype(str).str.upper()
    if "rank_positive_corr" in work.columns:
        top = work.sort_values("rank_positive_corr").head(top_n_drivers)
    else:
        top = work.sort_values("corr", ascending=False).head(top_n_drivers)
    top_driver_genes = set(top["gene_upper"])
    universe_set = {str(gene).upper() for gene in universe.astype(str)}
    universe_size = len(universe_set)
    rows = []
    for _, row in candidates.iterrows():
        tf = str(row["tf"]).upper()
        targets = parse_target_genes(row.get("target_genes", ""))
        targets = targets.intersection(universe_set)
        overlap = targets.intersection(top_driver_genes)
        a = len(overlap)
        b = max(len(targets) - a, 0)
        c = max(len(top_driver_genes.intersection(universe_set)) - a, 0)
        d = max(universe_size - a - b - c, 0)
        odds, pval = fisher_exact([[a, b], [c, d]], alternative="greater") if universe_size else (np.nan, np.nan)
        self_rows = work.loc[work["gene_upper"].eq(tf)].copy()
        self_corr = np.nan
        self_qval = np.nan
        self_rank = np.nan
        if not self_rows.empty:
            self_rows = self_rows.sort_values("rank_positive_corr" if "rank_positive_corr" in self_rows.columns else "corr")
            first = self_rows.iloc[0]
            self_corr = pd.to_numeric(first.get("corr", np.nan), errors="coerce")
            self_qval = pd.to_numeric(first.get("qval", np.nan), errors="coerce")
            self_rank = pd.to_numeric(first.get("rank_positive_corr", np.nan), errors="coerce")
        rows.append(
            {
                "tf": tf,
                "cellrank_top_driver_overlap_n": int(a),
                "cellrank_top_driver_overlap_genes": ";".join(sorted(overlap)),
                "cellrank_target_overlap_oddsratio": float(odds) if pd.notna(odds) else np.nan,
                "cellrank_target_overlap_p": float(pval) if pd.notna(pval) else np.nan,
                "tf_self_cellrank_corr": float(self_corr) if pd.notna(self_corr) else np.nan,
                "tf_self_cellrank_qval": float(self_qval) if pd.notna(self_qval) else np.nan,
                "tf_self_cellrank_rank": float(self_rank) if pd.notna(self_rank) else np.nan,
            }
        )
    out = pd.DataFrame(rows)
    out["cellrank_target_overlap_q"] = bh_qvalues(out["cellrank_target_overlap_p"])
    return out


def assign_roles(candidates: pd.DataFrame) -> pd.DataFrame:
    out = candidates.copy()
    tf_upper = out["tf"].astype(str).str.upper()
    out["tf"] = tf_upper
    out["role"] = np.where(tf_upper.isin(MAINTENANCE_TFS), "anti_cnv_hepatocyte_maintenance", "pro_cnv_candidate")
    out["perturbation_mode"] = np.where(
        out["role"].eq("anti_cnv_hepatocyte_maintenance"),
        "OE/rescue_primary;KO_direction_check",
        "KO/LOF",
    )
    out["expected_cnv_direction"] = np.where(out["role"].eq("anti_cnv_hepatocyte_maintenance"), -1.0, 1.0)
    flags = []
    for tf in tf_upper:
        notes = []
        if tf in IMMUNE_RISK_TFS:
            notes.append("immune_contamination_risk")
        if tf in INTERFERON_CONTEXT_TFS:
            notes.append("interferon_or_inflammatory_context")
        if tf in NONCANONICAL_EXCLUSIONS:
            notes.append("noncanonical_tf_exclusion")
        flags.append(";".join(notes))
    out["biology_risk_flag"] = flags
    out["biology_score"] = [BIOLOGY_SCORES.get(tf, 2.0) for tf in tf_upper]
    return out


def apply_hard_filters(
    candidates: pd.DataFrame,
    tf_catalog: set[str],
    expression_summary: pd.DataFrame,
    base_grn_compatibility: pd.DataFrame,
    maintenance_tfs: set[str],
    noncanonical_exclusions: set[str],
    min_base_grn_links: int = 10,
    min_detection_rate: float = 0.05,
    min_detected_datasets: int = 3,
    min_nes: float = 3.0,
    min_targets: int = 10,
    maintenance_min_targets: int = 8,
) -> pd.DataFrame:
    out = candidates.copy()
    out["tf"] = out["tf"].astype(str).str.upper()
    out = out.merge(expression_summary, on="tf", how="left")
    out = out.merge(base_grn_compatibility, on="tf", how="left")
    for col in ["tf_in_expression", "tf_in_base_grn", "base_grn_has_min_links"]:
        if col in out.columns:
            out[col] = out[col].fillna(False).astype(bool)
    for col in ["detection_rate_main", "mean_expression_main", "detected_dataset_count", "base_grn_outgoing_links"]:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce").fillna(0)
    catalog_upper = {str(tf).upper() for tf in tf_catalog}
    maintenance_upper = {str(tf).upper() for tf in maintenance_tfs}
    noncanonical_upper = {str(tf).upper() for tf in noncanonical_exclusions}

    reasons_by_row = []
    passes = []
    min_targets_used = []
    for _, row in out.iterrows():
        tf = str(row["tf"]).upper()
        reasons = []
        tf_in_catalog_or_base = tf in catalog_upper or bool(row.get("tf_in_base_grn", False))
        target_threshold = maintenance_min_targets if tf in maintenance_upper else min_targets
        min_targets_used.append(target_threshold)
        if tf in noncanonical_upper:
            reasons.append("noncanonical_exclusion")
        if not tf_in_catalog_or_base:
            reasons.append("not_in_tf_catalog_or_base_grn")
        if not bool(row.get("tf_in_expression", False)):
            reasons.append("not_in_expression_matrix")
        if not bool(row.get("tf_in_base_grn", False)):
            reasons.append("not_in_celloracle_base_grn")
        if float(row.get("base_grn_outgoing_links", 0)) < min_base_grn_links:
            reasons.append(f"base_grn_links_lt_{min_base_grn_links}")
        if float(row.get("detection_rate_main", 0)) < min_detection_rate:
            reasons.append(f"detection_rate_lt_{min_detection_rate:g}")
        if int(row.get("detected_dataset_count", 0)) < min_detected_datasets:
            reasons.append(f"detected_datasets_lt_{min_detected_datasets}")
        if float(row.get("best_nes", 0)) < min_nes:
            reasons.append(f"best_nes_lt_{min_nes:g}")
        if int(row.get("n_targets", 0)) < target_threshold:
            reasons.append(f"n_targets_lt_{target_threshold}")
        passes.append(len(reasons) == 0)
        reasons_by_row.append(";".join(reasons))
    out["min_targets_required"] = min_targets_used
    out["hard_filter_pass"] = passes
    out["exclusion_reason"] = reasons_by_row
    return out


def compute_phase_trends(auc: pd.DataFrame, cells: pd.DataFrame, candidates: pd.DataFrame, phase_key: str) -> pd.DataFrame:
    rows = []
    for _, row in candidates.iterrows():
        regulon = row["regulon"]
        tf = str(row["tf"]).upper()
        means = {}
        if phase_key in cells.columns and regulon in auc.columns:
            for phase in ["early", "middle", "late"]:
                idx = cells.index[cells[phase_key].astype(str).eq(phase)]
                idx = pd.Index(idx).intersection(auc.index)
                values = pd.to_numeric(auc.loc[idx, regulon], errors="coerce").dropna() if len(idx) else pd.Series(dtype=float)
                means[f"phase_{phase}_mean_auc"] = float(values.mean()) if not values.empty else np.nan
        rows.append(
            {
                "tf": tf,
                **means,
                "phase_late_minus_early_auc": float(
                    means.get("phase_late_mean_auc", np.nan) - means.get("phase_early_mean_auc", np.nan)
                )
                if pd.notna(means.get("phase_late_mean_auc", np.nan)) and pd.notna(means.get("phase_early_mean_auc", np.nan))
                else np.nan,
            }
        )
    return pd.DataFrame(rows)


def safe_pearson(x: pd.Series, y: pd.Series) -> tuple[float, float, int]:
    mask = x.notna() & y.notna() & np.isfinite(x) & np.isfinite(y)
    if int(mask.sum()) < 3 or x.loc[mask].nunique() < 2 or y.loc[mask].nunique() < 2:
        return np.nan, np.nan, int(mask.sum())
    r, p = pearsonr(x.loc[mask], y.loc[mask])
    return float(r), float(p), int(mask.sum())


def compute_group_correlations(
    auc: pd.DataFrame,
    cells: pd.DataFrame,
    candidates: pd.DataFrame,
    fate_key: str,
    group_cols: list[str],
    min_group_cells: int,
) -> pd.DataFrame:
    rows = []
    fate = pd.to_numeric(cells[fate_key], errors="coerce")
    for group_col in group_cols:
        if group_col not in cells.columns:
            continue
        for group, idx in cells.groupby(group_col, observed=True, sort=True).groups.items():
            idx = pd.Index(idx).intersection(auc.index)
            for _, row in candidates.iterrows():
                regulon = row["regulon"]
                values = pd.to_numeric(auc.loc[idx, regulon], errors="coerce") if regulon in auc.columns else pd.Series(dtype=float)
                r, p, n = safe_pearson(values, fate.loc[idx])
                rows.append(
                    {
                        "tf": str(row["tf"]).upper(),
                        "regulon": regulon,
                        "group_type": group_col,
                        "group": str(group),
                        "n_cells": n,
                        "eligible_group": bool(n >= min_group_cells and pd.notna(r)),
                        "pearson_r": r,
                        "pearson_p": p,
                    }
                )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["pearson_q"] = bh_qvalues(out["pearson_p"])
    return out


def compute_leave_one_dataset_out(
    auc: pd.DataFrame,
    cells: pd.DataFrame,
    candidates: pd.DataFrame,
    fate_key: str,
    dataset_key: str,
    min_remaining_cells: int,
) -> pd.DataFrame:
    if dataset_key not in cells.columns:
        return pd.DataFrame()
    fate = pd.to_numeric(cells[fate_key], errors="coerce")
    rows = []
    for omitted, _ in cells.groupby(dataset_key, observed=True, sort=True).groups.items():
        idx = cells.index[~cells[dataset_key].astype(str).eq(str(omitted))]
        idx = pd.Index(idx).intersection(auc.index)
        for _, row in candidates.iterrows():
            regulon = row["regulon"]
            values = pd.to_numeric(auc.loc[idx, regulon], errors="coerce") if regulon in auc.columns else pd.Series(dtype=float)
            r, p, n = safe_pearson(values, fate.loc[idx])
            rows.append(
                {
                    "tf": str(row["tf"]).upper(),
                    "regulon": regulon,
                    "omitted_dataset": str(omitted),
                    "remaining_n_cells": n,
                    "eligible_loo": bool(n >= min_remaining_cells and pd.notna(r)),
                    "loo_pearson_r": r,
                    "loo_pearson_p": p,
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["loo_pearson_q"] = bh_qvalues(out["loo_pearson_p"])
    return out


def compute_dataset_adjusted_betas(
    auc: pd.DataFrame,
    cells: pd.DataFrame,
    candidates: pd.DataFrame,
    fate_key: str,
    dataset_key: str,
) -> pd.DataFrame:
    rows = []
    fate = pd.to_numeric(cells[fate_key], errors="coerce")
    for _, row in candidates.iterrows():
        regulon = row["regulon"]
        tf = str(row["tf"]).upper()
        values = pd.to_numeric(auc[regulon], errors="coerce") if regulon in auc.columns else pd.Series(np.nan, index=cells.index)
        mask = values.notna() & fate.notna() & np.isfinite(values) & np.isfinite(fate)
        if dataset_key in cells.columns:
            mask &= cells[dataset_key].notna()
        if int(mask.sum()) < 10 or values.loc[mask].nunique() < 2:
            rows.append({"tf": tf, "dataset_adjusted_beta": np.nan, "dataset_adjusted_p": np.nan, "dataset_adjusted_n": int(mask.sum())})
            continue
        y = standardize(fate.loc[mask]).to_numpy(dtype=float)
        x = standardize(values.loc[mask]).to_numpy(dtype=float)
        if np.isnan(y).any() or np.isnan(x).any():
            rows.append({"tf": tf, "dataset_adjusted_beta": np.nan, "dataset_adjusted_p": np.nan, "dataset_adjusted_n": int(mask.sum())})
            continue
        design_parts = [np.ones(mask.sum()), x]
        if dataset_key in cells.columns:
            dummies = pd.get_dummies(cells.loc[mask, dataset_key].astype(str), drop_first=True, dtype=float)
            for col in dummies.columns:
                design_parts.append(dummies[col].to_numpy(dtype=float))
        xmat = np.column_stack(design_parts)
        beta, _, rank, _ = np.linalg.lstsq(xmat, y, rcond=None)
        residual = y - xmat @ beta
        df = max(len(y) - rank, 1)
        sigma2 = float((residual @ residual) / df)
        try:
            cov = sigma2 * np.linalg.pinv(xmat.T @ xmat)
            se = float(math.sqrt(max(cov[1, 1], 0)))
            t_stat = float(beta[1] / se) if se > 0 else np.nan
            pval = float(2 * t_dist.sf(abs(t_stat), df)) if pd.notna(t_stat) else np.nan
        except Exception:
            pval = np.nan
        rows.append(
            {
                "tf": tf,
                "dataset_adjusted_beta": float(beta[1]),
                "dataset_adjusted_p": pval,
                "dataset_adjusted_n": int(mask.sum()),
            }
        )
    out = pd.DataFrame(rows)
    out["dataset_adjusted_q"] = bh_qvalues(out["dataset_adjusted_p"])
    return out


def summarize_robustness(
    candidates: pd.DataFrame,
    group_corr: pd.DataFrame,
    loo: pd.DataFrame,
    adjusted: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for _, row in candidates.iterrows():
        tf = str(row["tf"]).upper()
        expected = float(row["expected_cnv_direction"])
        sub = group_corr.loc[group_corr["tf"].eq(tf) & group_corr["eligible_group"]].copy() if not group_corr.empty else pd.DataFrame()
        eligible = int(sub.shape[0])
        direction_fraction = float((sub["pearson_r"] * expected > 0).mean()) if eligible else np.nan
        dataset_sub = sub.loc[sub["group_type"].eq("dataset")] if not sub.empty else pd.DataFrame()
        dataset_direction_fraction = float((dataset_sub["pearson_r"] * expected > 0).mean()) if not dataset_sub.empty else np.nan

        loo_sub = loo.loc[loo["tf"].eq(tf) & loo["eligible_loo"]].copy() if not loo.empty else pd.DataFrame()
        loo_flip = bool((loo_sub["loo_pearson_r"] * expected < 0).any()) if not loo_sub.empty else False
        loo_min_directional = float((loo_sub["loo_pearson_r"] * expected).min()) if not loo_sub.empty else np.nan
        adj = adjusted.loc[adjusted["tf"].eq(tf)].iloc[0].to_dict() if not adjusted.loc[adjusted["tf"].eq(tf)].empty else {}
        rows.append(
            {
                "tf": tf,
                "eligible_group_count": eligible,
                "direction_consistency_fraction": direction_fraction,
                "dataset_direction_consistency_fraction": dataset_direction_fraction,
                "loo_dataset_count": int(loo_sub.shape[0]),
                "loo_any_direction_flip": loo_flip,
                "loo_min_directional_r": loo_min_directional,
                "dataset_adjusted_beta": adj.get("dataset_adjusted_beta", np.nan),
                "dataset_adjusted_p": adj.get("dataset_adjusted_p", np.nan),
                "dataset_adjusted_q": adj.get("dataset_adjusted_q", np.nan),
                "dataset_adjusted_n": adj.get("dataset_adjusted_n", np.nan),
            }
        )
    return pd.DataFrame(rows)


def score_candidates(candidates: pd.DataFrame) -> pd.DataFrame:
    out = candidates.copy()
    motif_scores = []
    fate_scores = []
    overlap_scores = []
    robustness_scores = []
    compatibility_scores = []
    totals = []
    for _, row in out.iterrows():
        expected = float(row.get("expected_cnv_direction", 1.0))
        directional_fate = expected * float(row.get("cnv_fate_pearson_r", 0) if pd.notna(row.get("cnv_fate_pearson_r", np.nan)) else 0)
        directional_time = expected * float(row.get("pseudotime_pearson_r", 0) if pd.notna(row.get("pseudotime_pearson_r", np.nan)) else 0)
        directional_phase = expected * float(row.get("phase_late_minus_early_auc", 0) if pd.notna(row.get("phase_late_minus_early_auc", np.nan)) else 0)
        directional_adjusted = expected * float(row.get("dataset_adjusted_beta", 0) if pd.notna(row.get("dataset_adjusted_beta", np.nan)) else 0)
        tf_self_directional = expected * float(row.get("tf_self_cellrank_corr", 0) if pd.notna(row.get("tf_self_cellrank_corr", np.nan)) else 0)

        fate_scale = 0.25 if expected < 0 else 0.7
        time_scale = 0.6 if expected < 0 else 0.7
        motif = (
            8.0 * clip01(float(row.get("best_nes", 0)) / 6.0)
            + 6.0 * clip01(float(row.get("direct_annotation_fraction", 0)) / 0.6)
            + 6.0 * clip01(float(row.get("n_targets", 0)) / 40.0)
        )
        fate = (
            10.0 * clip01(directional_fate / fate_scale)
            + (5.0 if float(row.get("cnv_fate_pearson_q", 1.0)) <= 0.05 else 0.0)
            + 5.0 * clip01(directional_time / time_scale)
            + (5.0 if directional_phase > 0 else 0.0)
        )
        overlap = (
            8.0 * clip01(float(row.get("cellrank_top_driver_overlap_n", 0)) / 8.0)
            + (6.0 if float(row.get("cellrank_target_overlap_q", 1.0)) <= 0.05 else 0.0)
            + 6.0 * clip01(tf_self_directional / 0.2)
        )
        consistency = row.get("direction_consistency_fraction", np.nan)
        consistency = float(consistency) if pd.notna(consistency) else 0.0
        no_flip = not bool(row.get("loo_any_direction_flip", False))
        robustness = (
            8.0 * clip01(consistency)
            + (6.0 if no_flip else 0.0)
            + 6.0 * clip01(directional_adjusted / 0.2)
        )
        compatibility = (
            (5.0 if bool(row.get("tf_in_base_grn", False)) and float(row.get("base_grn_outgoing_links", 0)) >= 10 else 0.0)
            + 2.5 * clip01(float(row.get("detection_rate_main", 0)) / 0.2)
            + 2.5 * clip01(float(row.get("detected_dataset_count", 0)) / 3.0)
        )
        biology = float(row.get("biology_score", 2.0))
        total = motif + fate + overlap + robustness + compatibility + biology
        motif_scores.append(round(motif, 3))
        fate_scores.append(round(fate, 3))
        overlap_scores.append(round(overlap, 3))
        robustness_scores.append(round(robustness, 3))
        compatibility_scores.append(round(compatibility, 3))
        totals.append(round(total, 3))
    out["motif_score"] = motif_scores
    out["fate_score"] = fate_scores
    out["cellrank_overlap_score"] = overlap_scores
    out["robustness_score"] = robustness_scores
    out["compatibility_score"] = compatibility_scores
    out["total_score"] = totals
    out["tier"] = np.where(
        ~out["hard_filter_pass"].astype(bool),
        "Tier C",
        np.where(out["total_score"] >= 75, "Tier A", np.where(out["total_score"] >= 60, "Tier B", "Tier C")),
    )
    return out


def select_main_panel(
    scored: pd.DataFrame,
    preferred_pro: list[str],
    preferred_maintenance: list[str],
    reserve_order: list[str],
    min_panel_size: int = 12,
    max_panel_size: int = 15,
) -> pd.DataFrame:
    out = scored.copy()
    out["selected_for_main_panel"] = False
    out["selection_order"] = np.nan
    tf_index = {str(tf).upper(): idx for idx, tf in out["tf"].astype(str).str.upper().items()}
    selected: list[str] = []

    def add_tf(tf: str) -> None:
        tf_upper = str(tf).upper()
        if len(selected) >= max_panel_size or tf_upper in selected or tf_upper not in tf_index:
            return
        row = out.loc[tf_index[tf_upper]]
        if bool(row.get("hard_filter_pass", False)):
            selected.append(tf_upper)

    for tf in preferred_pro:
        add_tf(tf)
    for tf in preferred_maintenance:
        add_tf(tf)
    if len(selected) < min_panel_size:
        for tf in reserve_order:
            add_tf(tf)
            if len(selected) >= min_panel_size:
                break
    if len(selected) < min_panel_size:
        ranked = out.loc[out["hard_filter_pass"].astype(bool)].sort_values("total_score", ascending=False)
        for tf in ranked["tf"].astype(str):
            add_tf(tf)
            if len(selected) >= min_panel_size:
                break
    for order, tf in enumerate(selected, start=1):
        idx = tf_index[tf]
        out.loc[idx, "selected_for_main_panel"] = True
        out.loc[idx, "selection_order"] = order
    return out


def plot_evidence_heatmap(scored: pd.DataFrame, path_base: Path) -> list[str]:
    display = scored.loc[scored["selected_for_main_panel"] | scored["tier"].isin(["Tier A", "Tier B"])].copy()
    if display.empty:
        return []
    display = display.sort_values(["selected_for_main_panel", "total_score"], ascending=[False, False])
    score_cols = [
        "motif_score",
        "fate_score",
        "cellrank_overlap_score",
        "robustness_score",
        "compatibility_score",
        "biology_score",
    ]
    matrix = display[score_cols].to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(6.6, max(3.2, 0.26 * display.shape[0])))
    im = ax.imshow(matrix, aspect="auto", cmap="viridis")
    labels = display["tf"] + "  " + display["tier"] + np.where(display["selected_for_main_panel"], " *", "")
    ax.set_yticks(range(display.shape[0]))
    ax.set_yticklabels(labels)
    ax.set_xticks(range(len(score_cols)))
    ax.set_xticklabels([col.replace("_score", "").replace("_", " ") for col in score_cols], rotation=35, ha="right")
    ax.set_title("Module 6.4 CellOracle TF evidence scores")
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Score component")
    fig.tight_layout()
    outputs = []
    for suffix in [".png", ".pdf"]:
        out = path_base.parent / f"{path_base.name}{suffix}"
        fig.savefig(out, bbox_inches="tight")
        outputs.append(str(out))
    plt.close(fig)
    return outputs


def order_output_columns(scored: pd.DataFrame) -> pd.DataFrame:
    required = [
        "tf",
        "role",
        "perturbation_mode",
        "tier",
        "total_score",
        "motif_score",
        "fate_score",
        "cellrank_overlap_score",
        "robustness_score",
        "compatibility_score",
        "biology_score",
        "hard_filter_pass",
        "exclusion_reason",
        "selected_for_main_panel",
    ]
    extras = [col for col in scored.columns if col not in required]
    return scored.loc[:, required + extras]


def main() -> None:
    start = time.time()
    args = parse_args()
    configure_plot_style()
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    args.figures_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.input_h5ad)
    cells = adata.obs.copy()
    cells.index = adata.obs_names.astype(str)
    if args.fate_key not in cells.columns:
        raise KeyError(f"Missing CellRank fate key in h5ad obs: {args.fate_key}")

    candidates = pd.read_csv(args.regulon_summary, sep="\t")
    required_regulon_cols = {"tf", "regulon", "best_nes", "n_targets", "target_genes"}
    missing_cols = required_regulon_cols.difference(candidates.columns)
    if missing_cols:
        raise KeyError(f"Regulon summary missing columns: {sorted(missing_cols)}")
    candidates = assign_roles(candidates)

    tf_catalog = read_tf_catalog(args.tf_catalog)
    base_grn = read_base_grn(args.celloracle_base_grn)
    base_compat = compute_base_grn_compatibility(base_grn, candidates["tf"], min_links=args.min_base_grn_links)
    expression = compute_expression_summary(adata, candidates["tf"], args.fate_key, args.dataset_key)

    drivers = pd.read_csv(args.cellrank_drivers, sep="\t")
    overlap = compute_cellrank_target_overlap(candidates, drivers, adata.var_names, args.top_n_cellrank_drivers)

    auc_raw = read_auc_matrix(args.auc)
    auc = align_auc_to_cells(auc_raw, pd.Index(adata.obs_names.astype(str)))
    phase = compute_phase_trends(auc, cells, candidates, args.phase_key)
    group_corr = compute_group_correlations(
        auc,
        cells,
        candidates,
        fate_key=args.fate_key,
        group_cols=[args.dataset_key, args.sample_key],
        min_group_cells=args.min_group_cells,
    )
    loo = compute_leave_one_dataset_out(
        auc,
        cells,
        candidates,
        fate_key=args.fate_key,
        dataset_key=args.dataset_key,
        min_remaining_cells=500,
    )
    adjusted = compute_dataset_adjusted_betas(auc, cells, candidates, args.fate_key, args.dataset_key)
    robustness = summarize_robustness(candidates, group_corr, loo, adjusted)

    filtered = apply_hard_filters(
        candidates,
        tf_catalog=tf_catalog,
        expression_summary=expression,
        base_grn_compatibility=base_compat,
        maintenance_tfs=MAINTENANCE_TFS,
        noncanonical_exclusions=NONCANONICAL_EXCLUSIONS,
        min_base_grn_links=args.min_base_grn_links,
        min_detection_rate=args.min_detection_rate,
        min_detected_datasets=args.min_detected_datasets,
    )
    scored = (
        filtered.merge(overlap, on="tf", how="left")
        .merge(phase, on="tf", how="left")
        .merge(robustness, on="tf", how="left")
    )
    scored = score_candidates(scored)
    scored = select_main_panel(
        scored,
        preferred_pro=PREFERRED_PRO_TFS,
        preferred_maintenance=PREFERRED_MAINTENANCE_TFS,
        reserve_order=RESERVE_TFS,
        min_panel_size=args.main_panel_size,
        max_panel_size=15,
    )
    scored = scored.sort_values(
        ["selected_for_main_panel", "selection_order", "hard_filter_pass", "total_score"],
        ascending=[False, True, False, False],
    ).reset_index(drop=True)
    scored["rank_by_total_score"] = scored["total_score"].rank(method="first", ascending=False).astype(int)
    scored = order_output_columns(scored)

    selection_path = args.metadata_dir / "celloracle_tf_selection.module6_4.tsv"
    report_path = args.metadata_dir / "celloracle_tf_selection.module6_4.json"
    input_tf_path = args.metadata_dir / "celloracle_input_tfs.module6_4.txt"
    group_path = args.metadata_dir / "celloracle_tf_selection_group_robustness.module6_4.tsv"
    loo_path = args.metadata_dir / "celloracle_tf_selection_leave_one_dataset_out.module6_4.tsv"
    figure_outputs = plot_evidence_heatmap(scored, args.figures_dir / "celloracle_tf_evidence_heatmap.module6_4")

    write_dataframe(selection_path, scored)
    write_dataframe(group_path, group_corr)
    write_dataframe(loo_path, loo)
    selected = scored.loc[scored["selected_for_main_panel"], "tf"].astype(str).tolist()
    input_tf_path.write_text("\n".join(selected) + "\n", encoding="utf-8")

    hard_filter_counts = scored["hard_filter_pass"].value_counts(dropna=False).to_dict()
    report = {
        "module": MODULE,
        "method": "CellOracle input TF selection for CellRank CNV-supported malignant fate",
        "goal": "select 12-15 TFs for CellOracle KO/OE simulation",
        "input_h5ad": str(args.input_h5ad),
        "regulon_summary": str(args.regulon_summary),
        "auc": str(args.auc),
        "cellrank_drivers": str(args.cellrank_drivers),
        "tf_catalog": str(args.tf_catalog),
        "celloracle_base_grn": str(args.celloracle_base_grn),
        "celloracle_base_grn_source_url": BASE_GRN_URL,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_fate_non_null_cells": int(cells[args.fate_key].notna().sum()),
        "n_candidate_regulons": int(candidates.shape[0]),
        "n_hard_filter_pass": int(scored["hard_filter_pass"].sum()),
        "n_selected_main_panel": int(scored["selected_for_main_panel"].sum()),
        "selected_tfs": selected,
        "hard_filter_counts": {str(key): int(value) for key, value in hard_filter_counts.items()},
        "thresholds": {
            "min_detection_rate": args.min_detection_rate,
            "min_detected_datasets": args.min_detected_datasets,
            "min_base_grn_links": args.min_base_grn_links,
            "min_best_nes": 3.0,
            "min_targets": 10,
            "maintenance_min_targets": 8,
            "top_n_cellrank_drivers": args.top_n_cellrank_drivers,
            "min_group_cells": args.min_group_cells,
        },
        "scoring_components": {
            "motif_score": 20,
            "fate_score": 25,
            "cellrank_overlap_score": 20,
            "robustness_score": 20,
            "compatibility_score": 10,
            "biology_score": 5,
        },
        "notes": [
            "Pro-CNV TFs are intended for KO/LOF CellOracle simulation.",
            "Anti-CNV hepatocyte-maintenance TFs are intended for OE/rescue-primary simulation.",
            "CellOracle compatibility uses hg38 promoter base GRN because matched scATAC is not available in this workspace.",
            "TPI1 is treated as a noncanonical TF exclusion unless an external TF catalog/base GRN explicitly supports it.",
        ],
        "references": [
            {
                "name": "CellOracle",
                "url": "https://www.nature.com/articles/s41586-022-05688-9",
            },
            {
                "name": "SCENIC/cisTarget/AUCell",
                "url": "https://www.nature.com/articles/nmeth.4463",
            },
            {
                "name": "CellRank",
                "url": "https://www.nature.com/articles/s41592-021-01346-6",
            },
        ],
        "outputs": {
            "selection_table": str(selection_path),
            "input_tfs": str(input_tf_path),
            "group_robustness": str(group_path),
            "leave_one_dataset_out": str(loo_path),
            "figures": figure_outputs,
            "report": str(report_path),
        },
        "package_versions": {
            "anndata": package_version("anndata"),
            "pandas": package_version("pandas"),
            "numpy": package_version("numpy"),
            "scipy": package_version("scipy"),
            "matplotlib": package_version("matplotlib"),
            "pyarrow": package_version("pyarrow"),
            "celloracle": package_version("celloracle"),
        },
        "elapsed_seconds": round(time.time() - start, 3),
    }
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
