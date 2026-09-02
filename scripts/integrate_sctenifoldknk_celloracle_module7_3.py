from __future__ import annotations

import argparse
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA_DIR = PROJECT_ROOT / "metadata/driver"
DEFAULT_PERTURBATION = DEFAULT_METADATA_DIR / "sctenifoldknk_module7_2_all_perturbation_genes.tsv"
DEFAULT_DRIVER_UNION = DEFAULT_METADATA_DIR / "sctenifoldknk_module7_2_driver_union_all_perturbation_genes.tsv"
DEFAULT_MALIGNANT_LIKE = DEFAULT_METADATA_DIR / "sctenifoldknk_module7_2_malignant_like_perturbation_genes.tsv"
DEFAULT_MAIN_STRICT = DEFAULT_METADATA_DIR / "sctenifoldknk_module7_2_main_strict_perturbation_genes.tsv"
DEFAULT_TIERS = DEFAULT_METADATA_DIR / "celloracle_module6_11_candidate_tier_table.tsv"
DEFAULT_GRN = DEFAULT_METADATA_DIR / "celloracle_module6_7_grn_links_filtered.tsv.gz"
DEFAULT_PHASE = DEFAULT_METADATA_DIR / "celloracle_module6_10_phase_wide_summary.tsv"

AP1_TFS = {"JUN", "FOS", "ATF3", "JUND"}
MARKER_PANELS = {
    "HCC_Malignant_Associated": ["AFP", "GPC3", "SPP1", "MDK", "IGF2BP3", "MUC1", "CEACAM5"],
    "Stressed_Injured": ["HSPA1A", "HSPA1B", "HSP90AA1", "DNAJB1", "FOS", "JUN", "JUNB", "ATF3", "DDIT3", "SAA1", "SAA2", "MT1G", "MT2A", "IER3"],
    "Proliferation": ["MKI67", "TOP2A", "STMN1", "TYMS", "UBE2C", "PCNA", "MCM2", "MCM5", "HMGB2"],
}


def _fdr_col(df: pd.DataFrame) -> str:
    for col in ["p.adj", "p_adj", "padj", "FDR", "fdr", "qvalue", "q_value"]:
        if col in df.columns:
            return col
    raise ValueError("No FDR-adjusted p-value column found")


def _distance_col(df: pd.DataFrame) -> str:
    for col in ["distance", "Distance", "dist", "perturbation_score", "score"]:
        if col in df.columns:
            return col
    raise ValueError("No scTenifoldKnk distance/score column found")


def _minmax(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce").astype(float)
    if values.notna().sum() == 0:
        return pd.Series(0.0, index=series.index)
    low = values.min()
    high = values.max()
    if np.isclose(high, low):
        return pd.Series(1.0, index=series.index)
    return (values - low) / (high - low)


def _top_genes(df: pd.DataFrame, top_n: int, fdr_col: str, distance_col: str) -> set[str]:
    ordered = df.copy()
    ordered[fdr_col] = pd.to_numeric(ordered[fdr_col], errors="coerce").fillna(1.0)
    ordered[distance_col] = pd.to_numeric(ordered[distance_col], errors="coerce").fillna(0.0)
    ordered = ordered.sort_values([fdr_col, distance_col, "gene"], ascending=[True, False, True])
    return set(ordered.head(top_n)["gene"].astype(str))


def build_concordance_summary(
    perturbation_genes: pd.DataFrame,
    grn_links: pd.DataFrame,
    fdr_threshold: float = 0.05,
    top_n: int = 100,
) -> pd.DataFrame:
    if perturbation_genes.empty:
        return pd.DataFrame()
    fdr_col = _fdr_col(perturbation_genes)
    distance_col = _distance_col(perturbation_genes)
    rows = []
    grn = grn_links.copy()
    if "source" not in grn.columns or "target" not in grn.columns:
        raise ValueError("CellOracle GRN links must contain source and target columns")
    grn["source"] = grn["source"].astype(str)
    grn["target"] = grn["target"].astype(str)

    group_cols = ["tf"]
    if "subset" in perturbation_genes.columns:
        group_cols = ["subset", "tf"]
    for group_key, tf_df in perturbation_genes.groupby(group_cols, sort=False):
        if isinstance(group_key, tuple) and len(group_key) == 2:
            subset, tf = group_key
        elif isinstance(group_key, tuple) and len(group_key) == 1:
            subset, tf = None, group_key[0]
        else:
            subset, tf = None, group_key
        work = tf_df.copy()
        work[fdr_col] = pd.to_numeric(work[fdr_col], errors="coerce").fillna(1.0)
        work[distance_col] = pd.to_numeric(work[distance_col], errors="coerce").fillna(0.0)
        significant = work.loc[work[fdr_col] <= fdr_threshold]
        top = _top_genes(work, top_n=top_n, fdr_col=fdr_col, distance_col=distance_col)
        targets = set(grn.loc[grn["source"] == str(tf), "target"])
        overlap = top & targets
        union = top | targets
        row = {
            "tf": str(tf),
            "n_tested_genes": int(work["gene"].astype(str).nunique()),
            "n_significant_perturbed_genes": int(significant["gene"].astype(str).nunique()),
            "mean_distance_significant": float(significant[distance_col].mean()) if len(significant) else 0.0,
            "max_distance": float(work[distance_col].max()) if len(work) else 0.0,
            "n_celloracle_grn_targets": int(len(targets)),
            "n_grn_overlap_genes": int(len(overlap)),
            "top_gene_grn_target_jaccard": float(len(overlap) / len(union)) if union else 0.0,
        }
        if subset is not None:
            row["subset"] = str(subset)
        rows.append(row)
    summary = pd.DataFrame(rows)
    summary["scTenifoldKnk_score"] = (
        _minmax(summary["n_significant_perturbed_genes"])
        + _minmax(summary["mean_distance_significant"])
        + _minmax(summary["top_gene_grn_target_jaccard"])
    ) / 3
    summary = summary.sort_values(
        ["scTenifoldKnk_score", "n_significant_perturbed_genes", "top_gene_grn_target_jaccard"],
        ascending=[False, False, False],
    ).reset_index(drop=True)
    summary["scTenifoldKnk_rank"] = np.arange(1, len(summary) + 1)
    return summary


def build_marker_overlap_summary(
    perturbation_genes: pd.DataFrame,
    marker_panels: dict[str, list[str]] | None = None,
    fdr_threshold: float = 0.05,
) -> pd.DataFrame:
    if marker_panels is None:
        marker_panels = MARKER_PANELS
    if perturbation_genes.empty:
        return pd.DataFrame()
    fdr_col = _fdr_col(perturbation_genes)
    distance_col = _distance_col(perturbation_genes)
    group_cols = ["tf"]
    if "subset" in perturbation_genes.columns:
        group_cols = ["subset", "tf"]
    rows = []
    for group_key, df in perturbation_genes.groupby(group_cols, sort=False):
        if isinstance(group_key, tuple) and len(group_key) == 2:
            subset, tf = group_key
        elif isinstance(group_key, tuple) and len(group_key) == 1:
            subset, tf = "unknown", group_key[0]
        else:
            subset, tf = "unknown", group_key
        work = df.copy()
        work["gene"] = work["gene"].astype(str)
        work[fdr_col] = pd.to_numeric(work[fdr_col], errors="coerce").fillna(1.0)
        work[distance_col] = pd.to_numeric(work[distance_col], errors="coerce").fillna(0.0)
        significant = work.loc[work[fdr_col] <= fdr_threshold]
        for panel, genes in marker_panels.items():
            panel_set = set(genes)
            panel_rows = work.loc[work["gene"].isin(panel_set)]
            sig_panel = significant.loc[significant["gene"].isin(panel_set)]
            sig_genes = (
                sig_panel.sort_values([distance_col, "gene"], ascending=[False, True])["gene"]
                .drop_duplicates()
                .astype(str)
                .tolist()
            )
            rows.append(
                {
                    "subset": str(subset),
                    "tf": str(tf),
                    "marker_panel": panel,
                    "n_panel_genes_available": int(panel_rows["gene"].nunique()),
                    "n_significant_marker_genes": int(len(sig_genes)),
                    "significant_marker_genes": ";".join(sig_genes),
                    "mean_marker_distance": float(panel_rows[distance_col].mean()) if len(panel_rows) else 0.0,
                    "mean_significant_marker_distance": float(sig_panel[distance_col].mean()) if len(sig_panel) else 0.0,
                }
            )
    return pd.DataFrame(rows)


def build_state_specific_gene_table(
    perturbation_genes: pd.DataFrame,
    fdr_threshold: float = 0.05,
    ratio_threshold: float = 2.0,
) -> pd.DataFrame:
    if "subset" not in perturbation_genes.columns:
        raise ValueError("subset column is required to compute malignant-like specificity")
    fdr_col = _fdr_col(perturbation_genes)
    distance_col = _distance_col(perturbation_genes)
    work = perturbation_genes.copy()
    work[fdr_col] = pd.to_numeric(work[fdr_col], errors="coerce").fillna(1.0)
    work[distance_col] = pd.to_numeric(work[distance_col], errors="coerce").fillna(0.0)
    rows = []
    for (tf, gene), df in work.groupby(["tf", "gene"], sort=False):
        by_subset = df.set_index("subset")
        if "malignant_like" not in by_subset.index:
            continue
        malignant = by_subset.loc["malignant_like"]
        if isinstance(malignant, pd.DataFrame):
            malignant = malignant.iloc[0]
        malignant_fdr = float(malignant[fdr_col])
        malignant_distance = float(malignant[distance_col])
        if malignant_fdr > fdr_threshold:
            continue
        comparator_distances = []
        comparator_fdrs = []
        for subset in ["main_strict", "driver_union_all"]:
            if subset in by_subset.index:
                comp = by_subset.loc[subset]
                if isinstance(comp, pd.DataFrame):
                    comp = comp.iloc[0]
                comparator_distances.append(float(comp[distance_col]))
                comparator_fdrs.append(float(comp[fdr_col]))
        max_comparator_distance = max(comparator_distances) if comparator_distances else 0.0
        min_comparator_fdr = min(comparator_fdrs) if comparator_fdrs else 1.0
        ratio = malignant_distance / (max_comparator_distance + 1e-12)
        is_specific = min_comparator_fdr > fdr_threshold or ratio >= ratio_threshold
        if is_specific:
            rows.append(
                {
                    "tf": str(tf),
                    "gene": str(gene),
                    "malignant_like_distance": malignant_distance,
                    "malignant_like_fdr": malignant_fdr,
                    "max_comparator_distance": max_comparator_distance,
                    "min_comparator_fdr": min_comparator_fdr,
                    "malignant_like_specificity_ratio": float(ratio),
                    "specificity_rule": "malignant_sig_comparator_nonsig"
                    if min_comparator_fdr > fdr_threshold
                    else f"malignant_distance_ratio_ge_{ratio_threshold:g}",
                }
            )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(["malignant_like_specificity_ratio", "malignant_like_distance"], ascending=[False, False]).reset_index(drop=True)


def build_integrated_evidence_matrix(concordance: pd.DataFrame, candidate_tiers: pd.DataFrame) -> pd.DataFrame:
    merged = candidate_tiers.merge(concordance, on="tf", how="left")
    numeric_defaults = [
        "n_significant_perturbed_genes",
        "top_gene_grn_target_jaccard",
        "mean_distance_significant",
        "scTenifoldKnk_score",
    ]
    for col in numeric_defaults:
        if col not in merged.columns:
            merged[col] = 0.0
        merged[col] = pd.to_numeric(merged[col], errors="coerce").fillna(0.0)
    if "quantitative_perturbation_score" not in merged.columns:
        merged["quantitative_perturbation_score"] = 0.0
    merged["celloracle_score_scaled"] = _minmax(merged["quantitative_perturbation_score"])
    merged["sctenifold_gene_count_scaled"] = _minmax(merged["n_significant_perturbed_genes"])
    merged["sctenifold_distance_scaled"] = _minmax(merged["mean_distance_significant"])
    merged["sctenifold_grn_overlap_scaled"] = _minmax(merged["top_gene_grn_target_jaccard"])
    merged["integrated_module7_score"] = (
        0.4 * merged["celloracle_score_scaled"]
        + 0.25 * merged["sctenifold_gene_count_scaled"]
        + 0.2 * merged["sctenifold_distance_scaled"]
        + 0.15 * merged["sctenifold_grn_overlap_scaled"]
    )
    merged = merged.sort_values(["integrated_module7_score", "tf"], ascending=[False, True]).reset_index(drop=True)
    merged["module7_integrated_rank"] = np.arange(1, len(merged) + 1)
    return merged


def _panel_count(marker_summary: pd.DataFrame, tf: str, panel: str, subset: str = "driver_union_all") -> int:
    if marker_summary.empty:
        return 0
    rows = marker_summary.loc[
        marker_summary["tf"].astype(str).eq(str(tf))
        & marker_summary["marker_panel"].astype(str).eq(panel)
        & marker_summary["subset"].astype(str).eq(subset)
    ]
    if rows.empty:
        return 0
    return int(pd.to_numeric(rows["n_significant_marker_genes"], errors="coerce").fillna(0).max())


def _panel_genes(marker_summary: pd.DataFrame, tf: str, panel: str, subset: str = "driver_union_all") -> str:
    if marker_summary.empty:
        return ""
    rows = marker_summary.loc[
        marker_summary["tf"].astype(str).eq(str(tf))
        & marker_summary["marker_panel"].astype(str).eq(panel)
        & marker_summary["subset"].astype(str).eq(subset)
    ]
    if rows.empty:
        return ""
    return str(rows.iloc[0].get("significant_marker_genes", ""))


def build_biological_axis_summary(
    integrated_matrix: pd.DataFrame,
    marker_summary: pd.DataFrame,
    phase_summary: pd.DataFrame,
) -> pd.DataFrame:
    if integrated_matrix.empty:
        return pd.DataFrame()
    phase_cols = [col for col in ["tf", "phase_early_rank", "phase_early_score", "phase_late_rank", "phase_late_score"] if col in phase_summary.columns]
    phase = phase_summary[phase_cols].copy() if phase_cols else pd.DataFrame({"tf": []})
    merged = integrated_matrix.merge(phase, on="tf", how="left", suffixes=("", "_phase"))
    for col in ["phase_early_rank", "phase_early_score", "phase_late_rank", "phase_late_score"]:
        phase_col = f"{col}_phase"
        if phase_col in merged.columns:
            if col in merged.columns:
                merged[col] = merged[col].combine_first(merged[phase_col])
            else:
                merged[col] = merged[phase_col]
    rows = []
    for _, row in merged.iterrows():
        tf = str(row["tf"])
        hcc_count = _panel_count(marker_summary, tf, "HCC_Malignant_Associated")
        stress_count = _panel_count(marker_summary, tf, "Stressed_Injured")
        proliferation_count = _panel_count(marker_summary, tf, "Proliferation")
        if str(row.get("candidate_tier", "")) == "Tier 1":
            interpretation = "Tier 1 HCC rescue replication"
        elif tf in AP1_TFS:
            interpretation = "AP-1 early/stress/proliferation axis"
        elif str(row.get("candidate_tier", "")) == "Tier 2":
            interpretation = "malignant-like state-specific replication"
        else:
            interpretation = "supplementary/control calibration"
        rows.append(
            {
                "tf": tf,
                "candidate_tier": row.get("candidate_tier", ""),
                "module7_integrated_rank": row.get("module7_integrated_rank", np.nan),
                "integrated_module7_score": row.get("integrated_module7_score", np.nan),
                "phase_early_rank": row.get("phase_early_rank", np.nan),
                "phase_early_score": row.get("phase_early_score", np.nan),
                "phase_late_rank": row.get("phase_late_rank", np.nan),
                "hcc_marker_count": hcc_count,
                "hcc_marker_genes": _panel_genes(marker_summary, tf, "HCC_Malignant_Associated"),
                "stress_marker_count": stress_count,
                "stress_marker_genes": _panel_genes(marker_summary, tf, "Stressed_Injured"),
                "proliferation_marker_count": proliferation_count,
                "proliferation_marker_genes": _panel_genes(marker_summary, tf, "Proliferation"),
                "module7_axis_interpretation": interpretation,
            }
        )
    return pd.DataFrame(rows).sort_values(["module7_integrated_rank", "tf"], ascending=[True, True]).reset_index(drop=True)


def flag_control_outperformance(integrated_matrix: pd.DataFrame) -> pd.DataFrame:
    if integrated_matrix.empty:
        return pd.DataFrame(columns=["risk_type", "tf", "detail", "severity"])
    tiers = integrated_matrix["candidate_tier"].astype(str)
    tier1 = integrated_matrix.loc[tiers.eq("Tier 1")]
    controls = integrated_matrix.loc[tiers.eq("Negative/control")]
    if tier1.empty or controls.empty:
        return pd.DataFrame(columns=["risk_type", "tf", "detail", "severity"])
    cutoff = float(pd.to_numeric(tier1["integrated_module7_score"], errors="coerce").max())
    rows = []
    for _, row in controls.iterrows():
        score = float(row["integrated_module7_score"])
        if score > cutoff:
            rows.append(
                {
                    "risk_type": "control_outperforms_tier1",
                    "tf": row["tf"],
                    "detail": f"Negative/control TF score {score:.3f} exceeds max Tier 1 score {cutoff:.3f}.",
                    "severity": "review_attention",
                }
            )
    return pd.DataFrame(rows)


def run_integration(
    perturbation_path: Path,
    tier_path: Path,
    grn_path: Path,
    phase_path: Path,
    metadata_dir: Path,
    fdr_threshold: float,
    top_n: int,
    malignant_like_path: Path | None = None,
    main_strict_path: Path | None = None,
) -> dict:
    perturb = pd.read_csv(perturbation_path, sep="\t")
    frames = [perturb]
    if malignant_like_path is not None and malignant_like_path.exists():
        frames.append(pd.read_csv(malignant_like_path, sep="\t"))
    if main_strict_path is not None and main_strict_path.exists():
        frames.append(pd.read_csv(main_strict_path, sep="\t"))
    all_perturb = pd.concat(frames, ignore_index=True)
    tiers = pd.read_csv(tier_path, sep="\t")
    grn = pd.read_csv(grn_path, sep="\t")
    phase = pd.read_csv(phase_path, sep="\t") if phase_path.exists() else pd.DataFrame()
    concordance = build_concordance_summary(perturb, grn, fdr_threshold=fdr_threshold, top_n=top_n)
    matrix = build_integrated_evidence_matrix(concordance, tiers)
    marker_summary = build_marker_overlap_summary(all_perturb, fdr_threshold=fdr_threshold)
    state_specific = build_state_specific_gene_table(all_perturb, fdr_threshold=fdr_threshold)
    axis_summary = build_biological_axis_summary(matrix, marker_summary, phase)
    risks = flag_control_outperformance(matrix)

    outputs = {
        "concordance_summary": str(metadata_dir / "sctenifoldknk_module7_3_concordance_summary.tsv"),
        "integrated_evidence_matrix": str(metadata_dir / "sctenifoldknk_module7_3_integrated_evidence_matrix.tsv"),
        "marker_overlap_summary": str(metadata_dir / "sctenifoldknk_module7_3_marker_overlap_summary.tsv"),
        "malignant_like_state_specific_genes": str(metadata_dir / "sctenifoldknk_module7_3_malignant_like_state_specific_genes.tsv"),
        "biological_axis_summary": str(metadata_dir / "sctenifoldknk_module7_3_biological_axis_summary.tsv"),
        "review_risk_flags": str(metadata_dir / "sctenifoldknk_module7_3_review_risk_flags.tsv"),
        "report": str(metadata_dir / "sctenifoldknk_module7_3_report.json"),
    }
    metadata_dir.mkdir(parents=True, exist_ok=True)
    concordance.to_csv(outputs["concordance_summary"], sep="\t", index=False)
    matrix.to_csv(outputs["integrated_evidence_matrix"], sep="\t", index=False)
    marker_summary.to_csv(outputs["marker_overlap_summary"], sep="\t", index=False)
    state_specific.to_csv(outputs["malignant_like_state_specific_genes"], sep="\t", index=False)
    axis_summary.to_csv(outputs["biological_axis_summary"], sep="\t", index=False)
    risks.to_csv(outputs["review_risk_flags"], sep="\t", index=False)
    report = {
        "module": "7.3",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "driver_union_perturbation": str(perturbation_path),
            "malignant_like_perturbation": str(malignant_like_path) if malignant_like_path else None,
            "main_strict_perturbation": str(main_strict_path) if main_strict_path else None,
            "candidate_tiers": str(tier_path),
            "celloracle_grn": str(grn_path),
            "celloracle_phase_summary": str(phase_path),
        },
        "outputs": outputs,
        "n_tfs": int(concordance["tf"].nunique()) if not concordance.empty else 0,
        "n_marker_overlap_rows": int(len(marker_summary)),
        "n_malignant_like_state_specific_genes": int(len(state_specific)),
        "n_biological_axis_rows": int(len(axis_summary)),
        "fdr_threshold": fdr_threshold,
        "top_n": top_n,
        "python_runtime": {"version": platform.python_version(), "platform": platform.platform()},
    }
    Path(outputs["report"]).write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 7.3 integrate scTenifoldKnk perturbation with CellOracle evidence")
    parser.add_argument("--perturbation", type=Path, default=DEFAULT_DRIVER_UNION)
    parser.add_argument("--malignant-like", type=Path, default=DEFAULT_MALIGNANT_LIKE)
    parser.add_argument("--main-strict", type=Path, default=DEFAULT_MAIN_STRICT)
    parser.add_argument("--tiers", type=Path, default=DEFAULT_TIERS)
    parser.add_argument("--grn", type=Path, default=DEFAULT_GRN)
    parser.add_argument("--phase", type=Path, default=DEFAULT_PHASE)
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA_DIR)
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    parser.add_argument("--top-n", type=int, default=100)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    report = run_integration(
        args.perturbation,
        args.tiers,
        args.grn,
        args.phase,
        args.metadata_dir,
        args.fdr_threshold,
        args.top_n,
        malignant_like_path=args.malignant_like,
        main_strict_path=args.main_strict,
    )
    print(json.dumps({"report": report["outputs"]["report"], "n_tfs": report["n_tfs"]}, indent=2))


if __name__ == "__main__":
    main()
