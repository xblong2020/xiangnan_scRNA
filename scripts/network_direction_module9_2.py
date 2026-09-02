from __future__ import annotations

import argparse
import json
import math
import platform
import time
from collections.abc import Iterable, Mapping, Sequence
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA_DIR = ROOT / "metadata/driver"
DEFAULT_FIGURE_DIR = ROOT / "figures/driver"
DEFAULT_CELLORACLE_TF_DELTA = DEFAULT_METADATA_DIR / "celloracle_module6_8_tf_delta_summary.tsv"
DEFAULT_CELLORACLE_TOP_GENES = DEFAULT_METADATA_DIR / "celloracle_module6_8_top_gene_delta_by_state.tsv.gz"
DEFAULT_CELLORACLE_QUANT_SCORES = DEFAULT_METADATA_DIR / "celloracle_module6_9b_quantitative_tf_scores.tsv"
DEFAULT_TF_TARGETS = DEFAULT_METADATA_DIR / "module8_tf_target_signature_genes.tsv"

A_TFS = ["HNF4A", "PPARA"]
AP1_TFS = ["JUN", "FOS", "JUND", "ATF3"]
CEBPB_EGR1_TFS = ["CEBPB", "EGR1"]
B_TFS = AP1_TFS + CEBPB_EGR1_TFS
C_TFS = ["SOX4"]
CONTROL_TFS = ["HLF", "IRF1", "JUNB", "MAFB", "MAFF", "MYC"]

AXIS_TF_GROUPS: dict[str, list[str]] = {
    "A_upstream": A_TFS,
    "B_ap1": AP1_TFS,
    "B_cebpb_egr1": CEBPB_EGR1_TFS,
    "B_transition": B_TFS,
    "C_sox4": C_TFS,
    "control": CONTROL_TFS,
}

CORE_TARGET_AXES = ["A_upstream", "B_transition", "C_sox4", "control"]
EPS = 1e-12


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 9.2 network direction evidence analysis.")
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--celloracle-tf-delta", type=Path, default=DEFAULT_CELLORACLE_TF_DELTA)
    parser.add_argument("--celloracle-top-genes", type=Path, default=DEFAULT_CELLORACLE_TOP_GENES)
    parser.add_argument("--celloracle-quantitative-scores", type=Path, default=DEFAULT_CELLORACLE_QUANT_SCORES)
    parser.add_argument("--tf-targets", type=Path, default=DEFAULT_TF_TARGETS)
    parser.add_argument("--sctenifold-pattern", default="sctenifoldknk_module7_2_*_perturbation_genes.tsv")
    parser.add_argument("--restore-celloracle-dir", type=Path, default=None)
    parser.add_argument("--require-restore", action="store_true")
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    parser.add_argument("--skip-figures", action="store_true")
    return parser.parse_args()


def package_version(name: str) -> str:
    try:
        return version(name)
    except PackageNotFoundError:
        return "not_installed"


def read_tsv_or_empty(path: Path, **kwargs) -> pd.DataFrame:
    if path is None or not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", **kwargs)


def normalize_gene(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip().upper()


def unique_clean(values: Iterable[object]) -> list[str]:
    cleaned = sorted({normalize_gene(value) for value in values if normalize_gene(value)})
    return cleaned


def infer_source_axis(tf: object) -> str:
    tf_name = normalize_gene(tf)
    if tf_name in A_TFS:
        return "A_upstream"
    if tf_name in B_TFS:
        return "B_transition"
    if tf_name in C_TFS:
        return "C_sox4"
    if tf_name in CONTROL_TFS:
        return "control"
    return "other"


def build_signature_sets(signature_table: pd.DataFrame) -> dict[str, list[str]]:
    """Build raw Module 8 and canonical A/B/C signature gene sets."""
    if signature_table.empty:
        return {axis: [] for axis in CORE_TARGET_AXES}
    required = {"axis", "tf", "gene"}
    missing = required.difference(signature_table.columns)
    if missing:
        raise ValueError(f"signature table missing columns: {sorted(missing)}")

    table = signature_table.copy()
    table["axis"] = table["axis"].astype(str)
    table["tf_norm"] = table["tf"].map(normalize_gene)
    table["gene_norm"] = table["gene"].map(normalize_gene)
    if "signature_class" not in table.columns:
        table["signature_class"] = ""

    signature_sets: dict[str, list[str]] = {}
    for axis, group in table.groupby("axis", dropna=False):
        signature_sets[str(axis)] = unique_clean(group["gene_norm"])

    signature_sets["A_upstream"] = unique_clean(table.loc[table["tf_norm"].isin(A_TFS), "gene_norm"])
    signature_sets["B_transition"] = unique_clean(
        table.loc[
            table["tf_norm"].isin(B_TFS)
            | table["axis"].astype(str).str.contains("ap1|transition|cebpb|egr1", case=False, regex=True),
            "gene_norm",
        ]
    )
    signature_sets["C_sox4"] = unique_clean(
        table.loc[
            table["tf_norm"].isin(C_TFS) | table["axis"].astype(str).str.contains("sox4", case=False, regex=True),
            "gene_norm",
        ]
    )
    signature_sets["control"] = unique_clean(
        table.loc[
            table["tf_norm"].isin(CONTROL_TFS)
            | table["axis"].astype(str).str.contains("control", case=False, regex=True)
            | table["signature_class"].astype(str).str.contains("control", case=False, regex=True),
            "gene_norm",
        ]
    )
    return signature_sets


def _rank_perturbation_genes(table: pd.DataFrame) -> pd.DataFrame:
    if table.empty:
        return table.copy()
    required = {"tf", "gene"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"scTenifoldKnk table missing columns: {sorted(missing)}")

    ranked = table.copy()
    ranked["perturb_tf"] = ranked["tf"].map(normalize_gene)
    ranked["gene_norm"] = ranked["gene"].map(normalize_gene)
    if "subset" not in ranked.columns:
        ranked["subset"] = "unknown"
    ranked["subset"] = ranked["subset"].astype(str)
    if "distance" not in ranked.columns:
        ranked["distance"] = np.nan
    if "p.adj" not in ranked.columns:
        ranked["p.adj"] = np.nan
    ranked["distance"] = pd.to_numeric(ranked["distance"], errors="coerce")
    ranked["p.adj"] = pd.to_numeric(ranked["p.adj"], errors="coerce")

    grouped = (
        ranked.groupby(["subset", "perturb_tf", "gene_norm"], as_index=False)
        .agg(distance=("distance", "max"), p_adj=("p.adj", "min"))
        .sort_values(["subset", "perturb_tf", "p_adj", "distance"], ascending=[True, True, True, False])
    )
    grouped["rank"] = grouped.groupby(["subset", "perturb_tf"]).cumcount() + 1
    return grouped


def _safe_neg_log10(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce").clip(lower=1e-300)
    return -np.log10(numeric)


def compute_sctenifold_signature_impact(
    perturbation_genes: pd.DataFrame,
    signature_sets: Mapping[str, Sequence[str]],
    fdr_threshold: float = 0.05,
) -> pd.DataFrame:
    ranked = _rank_perturbation_genes(perturbation_genes)
    rows: list[dict[str, object]] = []
    if ranked.empty:
        return pd.DataFrame(
            columns=[
                "subset",
                "perturb_tf",
                "source_axis",
                "target_axis",
                "n_signature_genes",
                "n_overlap",
                "n_sig_fdr05",
                "sig_fraction",
                "overlap_fraction",
                "mean_rank",
                "median_rank",
                "mean_distance",
                "median_distance",
                "mean_neg_log10_fdr",
                "impact_score",
                "signature_gene_hits",
            ]
        )

    for (subset, perturb_tf), group in ranked.groupby(["subset", "perturb_tf"], sort=True):
        gene_index = group.set_index("gene_norm", drop=False)
        for target_axis, genes in signature_sets.items():
            signature_genes = sorted({normalize_gene(gene) for gene in genes if normalize_gene(gene)})
            if not signature_genes:
                continue
            overlap = gene_index.loc[gene_index.index.intersection(signature_genes)].copy()
            n_overlap = int(overlap["gene_norm"].nunique()) if not overlap.empty else 0
            sig = overlap.loc[overlap["p_adj"] <= fdr_threshold].copy() if not overlap.empty else overlap
            n_sig = int(sig["gene_norm"].nunique()) if not sig.empty else 0
            sig_fraction = n_sig / len(signature_genes)
            overlap_fraction = n_overlap / len(signature_genes)
            mean_rank = float(overlap["rank"].mean()) if not overlap.empty else math.nan
            median_rank = float(overlap["rank"].median()) if not overlap.empty else math.nan
            mean_distance = float(overlap["distance"].mean()) if not overlap.empty else 0.0
            median_distance = float(overlap["distance"].median()) if not overlap.empty else 0.0
            mean_neg_log10_fdr = float(_safe_neg_log10(sig["p_adj"]).mean()) if not sig.empty else 0.0
            impact_score = float(sig_fraction * (1.0 + mean_neg_log10_fdr / 10.0) * (1.0 + math.log1p(max(mean_distance, 0.0))))
            rows.append(
                {
                    "subset": subset,
                    "perturb_tf": perturb_tf,
                    "source_axis": infer_source_axis(perturb_tf),
                    "target_axis": target_axis,
                    "n_signature_genes": len(signature_genes),
                    "n_overlap": n_overlap,
                    "n_sig_fdr05": n_sig,
                    "sig_fraction": sig_fraction,
                    "overlap_fraction": overlap_fraction,
                    "mean_rank": mean_rank,
                    "median_rank": median_rank,
                    "mean_distance": mean_distance,
                    "median_distance": median_distance,
                    "mean_neg_log10_fdr": mean_neg_log10_fdr,
                    "impact_score": impact_score,
                    "signature_gene_hits": ";".join(sorted(sig["gene_norm"].unique())) if not sig.empty else "",
                }
            )
    return pd.DataFrame(rows)


def compute_celloracle_tf_axis_delta(
    tf_delta_summary: pd.DataFrame,
    axis_tf_groups: Mapping[str, Sequence[str]] = AXIS_TF_GROUPS,
) -> pd.DataFrame:
    if tf_delta_summary.empty:
        return pd.DataFrame(
            columns=[
                "perturb_tf",
                "source_axis",
                "target_axis",
                "celloracle_state",
                "n_target_tfs_expected",
                "n_target_tfs_observed",
                "mean_delta_x",
                "median_delta_x",
                "mean_abs_delta_x",
                "target_tfs_observed",
            ]
        )
    required = {"tf", "target_tf", "celloracle_state", "mean_delta_x"}
    missing = required.difference(tf_delta_summary.columns)
    if missing:
        raise ValueError(f"CellOracle TF delta table missing columns: {sorted(missing)}")

    table = tf_delta_summary.copy()
    table["perturb_tf"] = table["tf"].map(normalize_gene)
    table["target_tf_norm"] = table["target_tf"].map(normalize_gene)
    table["celloracle_state"] = table["celloracle_state"].astype(str)
    table["mean_delta_x"] = pd.to_numeric(table["mean_delta_x"], errors="coerce")
    if "median_delta_x" not in table.columns:
        table["median_delta_x"] = table["mean_delta_x"]
    table["median_delta_x"] = pd.to_numeric(table["median_delta_x"], errors="coerce")
    if "mean_abs_delta_x" not in table.columns:
        table["mean_abs_delta_x"] = table["mean_delta_x"].abs()
    table["mean_abs_delta_x"] = pd.to_numeric(table["mean_abs_delta_x"], errors="coerce")

    pooled = table.copy()
    pooled["celloracle_state"] = "all_states_pooled"
    table = pd.concat([table, pooled], ignore_index=True)

    rows: list[dict[str, object]] = []
    for (perturb_tf, state), group in table.groupby(["perturb_tf", "celloracle_state"], sort=True):
        for axis, target_tfs in axis_tf_groups.items():
            target_set = {normalize_gene(tf) for tf in target_tfs}
            observed = group.loc[group["target_tf_norm"].isin(target_set)].copy()
            rows.append(
                {
                    "perturb_tf": perturb_tf,
                    "source_axis": infer_source_axis(perturb_tf),
                    "target_axis": axis,
                    "celloracle_state": state,
                    "n_target_tfs_expected": len(target_set),
                    "n_target_tfs_observed": int(observed["target_tf_norm"].nunique()) if not observed.empty else 0,
                    "mean_delta_x": float(observed["mean_delta_x"].mean()) if not observed.empty else math.nan,
                    "median_delta_x": float(observed["median_delta_x"].median()) if not observed.empty else math.nan,
                    "mean_abs_delta_x": float(observed["mean_abs_delta_x"].mean()) if not observed.empty else math.nan,
                    "target_tfs_observed": ";".join(sorted(observed["target_tf_norm"].unique())) if not observed.empty else "",
                }
            )
    return pd.DataFrame(rows)


def compute_celloracle_top_gene_signature_overlap(
    top_gene_delta: pd.DataFrame,
    signature_sets: Mapping[str, Sequence[str]],
) -> pd.DataFrame:
    if top_gene_delta.empty:
        return pd.DataFrame(
            columns=[
                "perturb_tf",
                "source_axis",
                "target_axis",
                "celloracle_state",
                "n_signature_genes",
                "n_overlap",
                "overlap_fraction",
                "mean_delta_x",
                "mean_abs_delta_x",
                "signature_gene_hits",
            ]
        )
    required = {"tf", "celloracle_state", "gene", "mean_delta_x"}
    missing = required.difference(top_gene_delta.columns)
    if missing:
        raise ValueError(f"CellOracle top gene table missing columns: {sorted(missing)}")

    table = top_gene_delta.copy()
    table["perturb_tf"] = table["tf"].map(normalize_gene)
    table["gene_norm"] = table["gene"].map(normalize_gene)
    table["celloracle_state"] = table["celloracle_state"].astype(str)
    table["mean_delta_x"] = pd.to_numeric(table["mean_delta_x"], errors="coerce")
    if "abs_mean_delta_x" not in table.columns:
        table["abs_mean_delta_x"] = table["mean_delta_x"].abs()
    table["abs_mean_delta_x"] = pd.to_numeric(table["abs_mean_delta_x"], errors="coerce")

    rows: list[dict[str, object]] = []
    for (perturb_tf, state), group in table.groupby(["perturb_tf", "celloracle_state"], sort=True):
        gene_index = group.set_index("gene_norm", drop=False)
        for target_axis, genes in signature_sets.items():
            signature_genes = sorted({normalize_gene(gene) for gene in genes if normalize_gene(gene)})
            if not signature_genes:
                continue
            overlap = gene_index.loc[gene_index.index.intersection(signature_genes)].copy()
            rows.append(
                {
                    "perturb_tf": perturb_tf,
                    "source_axis": infer_source_axis(perturb_tf),
                    "target_axis": target_axis,
                    "celloracle_state": state,
                    "n_signature_genes": len(signature_genes),
                    "n_overlap": int(overlap["gene_norm"].nunique()) if not overlap.empty else 0,
                    "overlap_fraction": float(overlap["gene_norm"].nunique() / len(signature_genes)) if not overlap.empty else 0.0,
                    "mean_delta_x": float(overlap["mean_delta_x"].mean()) if not overlap.empty else math.nan,
                    "mean_abs_delta_x": float(overlap["abs_mean_delta_x"].mean()) if not overlap.empty else math.nan,
                    "signature_gene_hits": ";".join(sorted(overlap["gene_norm"].unique())) if not overlap.empty else "",
                }
            )
    return pd.DataFrame(rows)


def compute_celloracle_axis_impact_matrix(celloracle_axis_delta: pd.DataFrame) -> pd.DataFrame:
    if celloracle_axis_delta.empty:
        return pd.DataFrame(
            columns=[
                "subset",
                "perturb_tf",
                "source_axis",
                "target_axis",
                "impact_score",
                "signed_delta",
                "sig_fraction",
                "n_sig_fdr05",
                "n_signature_genes",
                "evidence_source",
            ]
        )
    required = {"perturb_tf", "source_axis", "target_axis", "celloracle_state", "mean_delta_x", "mean_abs_delta_x"}
    missing = required.difference(celloracle_axis_delta.columns)
    if missing:
        raise ValueError(f"CellOracle axis delta table missing columns: {sorted(missing)}")
    matrix = celloracle_axis_delta.copy()
    matrix["subset"] = "celloracle_" + matrix["celloracle_state"].astype(str)
    matrix["impact_score"] = pd.to_numeric(matrix["mean_abs_delta_x"], errors="coerce")
    matrix["signed_delta"] = pd.to_numeric(matrix["mean_delta_x"], errors="coerce")
    matrix["sig_fraction"] = np.nan
    matrix["n_sig_fdr05"] = np.nan
    matrix["n_signature_genes"] = matrix.get("n_target_tfs_expected", np.nan)
    matrix["evidence_source"] = "celloracle_abs_delta"
    return matrix[
        [
            "subset",
            "perturb_tf",
            "source_axis",
            "target_axis",
            "impact_score",
            "signed_delta",
            "sig_fraction",
            "n_sig_fdr05",
            "n_signature_genes",
            "evidence_source",
        ]
    ].copy()


def _mean_impact(
    impact: pd.DataFrame,
    source_axes: Sequence[str],
    target_axes: Sequence[str],
    subset: str | None = None,
) -> float:
    if impact.empty:
        return math.nan
    mask = impact["source_axis"].isin(source_axes) & impact["target_axis"].isin(target_axes)
    if subset is not None and "subset" in impact.columns:
        mask &= impact["subset"].eq(subset)
    values = pd.to_numeric(impact.loc[mask, "impact_score"], errors="coerce").dropna()
    if values.empty:
        return math.nan
    return float(values.mean())


def _directionality_index(forward: float, reverse: float) -> float:
    if math.isnan(forward) or math.isnan(reverse):
        return math.nan
    return float((forward - reverse) / (abs(forward) + abs(reverse) + EPS))


def compute_asymmetry_tests(signature_impact: pd.DataFrame, evidence_source: str = "sctenifold_signature") -> pd.DataFrame:
    if signature_impact.empty:
        return pd.DataFrame(
            [
                {
                    "comparison": "network_direction",
                    "subset": "all",
                    "forward_impact": math.nan,
                    "reverse_impact": math.nan,
                    "directionality_index": math.nan,
                    "forward_vs_control_ratio": math.nan,
                    "support_label": "insufficient_data",
                    "evidence_source": evidence_source,
                }
            ]
        )

    subsets = sorted(signature_impact["subset"].dropna().astype(str).unique()) if "subset" in signature_impact.columns else ["all"]
    rows: list[dict[str, object]] = []
    comparison_specs = [
        ("A_to_B_vs_C_to_A", ["A_upstream"], ["B_transition"], ["C_sox4"], ["A_upstream"], ["B_transition"]),
        ("A_to_C_vs_C_to_A", ["A_upstream"], ["C_sox4"], ["C_sox4"], ["A_upstream"], ["C_sox4"]),
        ("B_to_C_vs_C_to_B", ["B_transition"], ["C_sox4"], ["C_sox4"], ["B_transition"], ["C_sox4"]),
    ]

    for subset in subsets:
        subset_value = subset if "subset" in signature_impact.columns else None
        for name, f_source, f_target, r_source, r_target, control_target in comparison_specs:
            forward = _mean_impact(signature_impact, f_source, f_target, subset_value)
            reverse = _mean_impact(signature_impact, r_source, r_target, subset_value)
            control = _mean_impact(signature_impact, ["control"], control_target, subset_value)
            directionality = _directionality_index(forward, reverse)
            if math.isnan(forward) or math.isnan(reverse):
                label = "insufficient_data"
            elif forward > reverse * 1.25 and (math.isnan(control) or forward > control * 1.25):
                label = "forward_greater_than_reverse"
            elif forward > reverse * 1.25:
                label = "partial_forward_greater_than_reverse"
            else:
                label = "no_forward_asymmetry"
            rows.append(
                {
                    "comparison": name,
                    "subset": subset,
                    "forward_source_axis": ",".join(f_source),
                    "forward_target_axis": ",".join(f_target),
                    "reverse_source_axis": ",".join(r_source),
                    "reverse_target_axis": ",".join(r_target),
                    "forward_impact": forward,
                    "reverse_impact": reverse,
                    "control_impact": control,
                    "directionality_index": directionality,
                    "forward_reverse_ratio": forward / (reverse + EPS) if not math.isnan(forward) and not math.isnan(reverse) else math.nan,
                    "forward_vs_control_ratio": forward / (control + EPS) if not math.isnan(forward) and not math.isnan(control) else math.nan,
                    "support_label": label,
                    "evidence_source": evidence_source,
                }
            )

        sox4_self = _mean_impact(signature_impact, ["C_sox4"], ["C_sox4"], subset_value)
        sox4_reverse = _mean_impact(signature_impact, ["C_sox4"], ["A_upstream", "B_transition"], subset_value)
        control = _mean_impact(signature_impact, ["control"], ["C_sox4"], subset_value)
        directionality = _directionality_index(sox4_self, sox4_reverse)
        if math.isnan(sox4_self) or math.isnan(sox4_reverse):
            label = "insufficient_data"
        elif sox4_self > sox4_reverse * 2.0:
            label = "weak_reverse_upstream"
        elif sox4_self > sox4_reverse * 1.25:
            label = "partial_weak_reverse_upstream"
        else:
            label = "broad_reverse_impact"
        rows.append(
            {
                "comparison": "SOX4_self_vs_reverse_upstream",
                "subset": subset,
                "forward_source_axis": "C_sox4",
                "forward_target_axis": "C_sox4",
                "reverse_source_axis": "C_sox4",
                "reverse_target_axis": "A_upstream,B_transition",
                "forward_impact": sox4_self,
                "reverse_impact": sox4_reverse,
                "control_impact": control,
                "directionality_index": directionality,
                "forward_reverse_ratio": sox4_self / (sox4_reverse + EPS)
                if not math.isnan(sox4_self) and not math.isnan(sox4_reverse)
                else math.nan,
                "forward_vs_control_ratio": sox4_self / (control + EPS)
                if not math.isnan(sox4_self) and not math.isnan(control)
                else math.nan,
                "support_label": label,
                "evidence_source": evidence_source,
            }
        )

    consensus = pd.DataFrame(rows)
    if not consensus.empty:
        supported_labels = {
            "forward_greater_than_reverse",
            "partial_forward_greater_than_reverse",
            "weak_reverse_upstream",
            "partial_weak_reverse_upstream",
        }
        summary = (
            consensus.assign(is_supported=consensus["support_label"].isin(supported_labels))
            .groupby("comparison", as_index=False)
            .agg(
                subset=("subset", lambda values: "consensus"),
                forward_impact=("forward_impact", "mean"),
                reverse_impact=("reverse_impact", "mean"),
                control_impact=("control_impact", "mean"),
                directionality_index=("directionality_index", "mean"),
                forward_reverse_ratio=("forward_reverse_ratio", "mean"),
                forward_vs_control_ratio=("forward_vs_control_ratio", "mean"),
                n_supported_subsets=("is_supported", "sum"),
                n_tested_subsets=("is_supported", "size"),
            )
        )
        summary["support_label"] = np.where(
            summary["n_supported_subsets"] >= np.ceil(summary["n_tested_subsets"] / 2),
            "consensus_supported",
            "consensus_not_supported",
        )
        summary["evidence_source"] = evidence_source
        for col in ["forward_source_axis", "forward_target_axis", "reverse_source_axis", "reverse_target_axis"]:
            summary[col] = "consensus"
        consensus = pd.concat([consensus, summary[consensus.columns]], ignore_index=True, sort=False)
    return consensus


def combine_asymmetry_sources(*tables: pd.DataFrame) -> pd.DataFrame:
    frames = [table.copy() for table in tables if table is not None and not table.empty]
    if not frames:
        return pd.DataFrame()
    combined = pd.concat(frames, ignore_index=True, sort=False)
    consensus = combined.loc[combined["subset"].eq("consensus")].copy()
    if consensus.empty:
        return combined

    supported_labels = {"consensus_supported"}
    rows: list[dict[str, object]] = []
    for comparison, group in consensus.groupby("comparison", sort=True):
        n_sources = int(group["evidence_source"].nunique())
        n_supported = int(group["support_label"].isin(supported_labels).sum())
        forward_mean = float(pd.to_numeric(group["forward_impact"], errors="coerce").mean())
        reverse_mean = float(pd.to_numeric(group["reverse_impact"], errors="coerce").mean())
        control_mean = float(pd.to_numeric(group["control_impact"], errors="coerce").mean())
        if n_sources >= 2 and n_supported == n_sources:
            label = "joint_consensus_supported"
        elif n_supported > 0:
            label = "partial_joint_consensus_support"
        else:
            label = "joint_consensus_not_supported"
        rows.append(
            {
                "comparison": comparison,
                "subset": "joint_consensus",
                "forward_source_axis": "joint_consensus",
                "forward_target_axis": "joint_consensus",
                "reverse_source_axis": "joint_consensus",
                "reverse_target_axis": "joint_consensus",
                "forward_impact": forward_mean,
                "reverse_impact": reverse_mean,
                "control_impact": control_mean,
                "directionality_index": float(pd.to_numeric(group["directionality_index"], errors="coerce").mean()),
                "forward_reverse_ratio": forward_mean / (reverse_mean + EPS),
                "forward_vs_control_ratio": forward_mean / (control_mean + EPS),
                "support_label": label,
                "evidence_source": "combined",
                "n_supported_sources": n_supported,
                "n_tested_sources": n_sources,
            }
        )
    return pd.concat([combined, pd.DataFrame(rows)], ignore_index=True, sort=False)


def audit_restore_availability(restore_celloracle_dir: Path | None, require_restore: bool = False) -> pd.DataFrame:
    if restore_celloracle_dir is None or not restore_celloracle_dir.exists():
        status = "required_but_missing" if require_restore else "not_available_existing_outputs"
        return pd.DataFrame(
            [
                {
                    "restore_dir": str(restore_celloracle_dir) if restore_celloracle_dir is not None else "",
                    "restore_available": False,
                    "restore_status": status,
                    "n_tf_delta_files": 0,
                    "n_top_gene_files": 0,
                }
            ]
        )
    tf_delta_files = list(restore_celloracle_dir.rglob("*tf_delta_summary*.tsv"))
    top_gene_files = list(restore_celloracle_dir.rglob("*top_gene_delta*.tsv*"))
    available = bool(tf_delta_files or top_gene_files)
    status = "available" if available else ("required_but_missing" if require_restore else "not_available_existing_outputs")
    return pd.DataFrame(
        [
            {
                "restore_dir": str(restore_celloracle_dir),
                "restore_available": available,
                "restore_status": status,
                "n_tf_delta_files": len(tf_delta_files),
                "n_top_gene_files": len(top_gene_files),
            }
        ]
    )


def summarize_network_evidence(asymmetry_tests: pd.DataFrame, restore_audit: pd.DataFrame) -> pd.DataFrame:
    joint = asymmetry_tests.loc[asymmetry_tests["subset"].eq("joint_consensus")].copy() if not asymmetry_tests.empty else pd.DataFrame()
    source_consensus = asymmetry_tests.loc[asymmetry_tests["subset"].eq("consensus")].copy() if not asymmetry_tests.empty else pd.DataFrame()
    supported = joint["support_label"].eq("joint_consensus_supported") if not joint.empty else pd.Series(dtype=bool)
    partial = joint["support_label"].eq("partial_joint_consensus_support") if not joint.empty else pd.Series(dtype=bool)
    n_supported = int(supported.sum()) if not supported.empty else 0
    n_partial = int(partial.sum()) if not partial.empty else 0
    n_core = int(len(joint)) if not joint.empty else 0
    n_sctenifold_supported = int(
        source_consensus.loc[source_consensus["evidence_source"].eq("sctenifold_signature"), "support_label"].eq("consensus_supported").sum()
    ) if not source_consensus.empty else 0
    n_celloracle_supported = int(
        source_consensus.loc[source_consensus["evidence_source"].eq("celloracle_abs_delta"), "support_label"].eq("consensus_supported").sum()
    ) if not source_consensus.empty else 0
    restore_status = (
        str(restore_audit.loc[0, "restore_status"]) if not restore_audit.empty and "restore_status" in restore_audit.columns else "unknown"
    )
    restore_available = (
        bool(restore_audit.loc[0, "restore_available"]) if not restore_audit.empty and "restore_available" in restore_audit.columns else False
    )

    if n_core == 0:
        label = "network_direction_not_testable"
    elif n_supported >= 3:
        label = "ko_network_direction_supported"
    elif n_supported >= 2 or (n_supported + n_partial) >= 3:
        label = "partial_ko_network_direction_support"
    else:
        label = "ko_network_direction_not_supported"
    if not restore_available and label in {"ko_network_direction_supported", "partial_ko_network_direction_support"}:
        label = f"{label}_restore_not_available"

    return pd.DataFrame(
        [
            {
                "network_direction_label": label,
                "n_joint_consensus_comparisons": n_core,
                "n_supported_joint_consensus_comparisons": n_supported,
                "n_partial_joint_consensus_comparisons": n_partial,
                "n_sctenifold_supported_consensus_comparisons": n_sctenifold_supported,
                "n_celloracle_supported_consensus_comparisons": n_celloracle_supported,
                "restore_available": restore_available,
                "restore_status": restore_status,
                "evidence_scope": "KO network-direction evidence only" if not restore_available else "KO plus supplied restore evidence",
            }
        ]
    )


def load_sctenifold_tables(metadata_dir: Path, pattern: str) -> tuple[pd.DataFrame, list[str]]:
    paths = sorted(metadata_dir.glob(pattern))
    frames: list[pd.DataFrame] = []
    for path in paths:
        df = read_tsv_or_empty(path)
        if df.empty:
            continue
        if "subset" not in df.columns:
            subset = path.name.replace("sctenifoldknk_module7_2_", "").replace("_perturbation_genes.tsv", "")
            df["subset"] = subset
        frames.append(df)
    if not frames:
        return pd.DataFrame(), [str(path) for path in paths]
    return pd.concat(frames, ignore_index=True, sort=False), [str(path) for path in paths]


def build_directionality_matrix(signature_impact: pd.DataFrame) -> pd.DataFrame:
    if signature_impact.empty:
        return pd.DataFrame()
    canonical = signature_impact.loc[signature_impact["target_axis"].isin(CORE_TARGET_AXES)].copy()
    if canonical.empty:
        return pd.DataFrame()
    return (
        canonical.groupby(["subset", "perturb_tf", "source_axis", "target_axis"], as_index=False)
        .agg(
            impact_score=("impact_score", "mean"),
            sig_fraction=("sig_fraction", "mean"),
            n_sig_fdr05=("n_sig_fdr05", "sum"),
            n_signature_genes=("n_signature_genes", "mean"),
        )
        .sort_values(["subset", "source_axis", "perturb_tf", "target_axis"])
    )


def merge_celloracle_sctenifold_concordance(
    directionality_matrix: pd.DataFrame,
    celloracle_axis_delta: pd.DataFrame,
) -> pd.DataFrame:
    if directionality_matrix.empty or celloracle_axis_delta.empty:
        return pd.DataFrame()
    oracle = celloracle_axis_delta.loc[celloracle_axis_delta["celloracle_state"].eq("malignant_or_malignant_like")].copy()
    if oracle.empty:
        oracle = celloracle_axis_delta.loc[celloracle_axis_delta["celloracle_state"].eq("all_states_pooled")].copy()
    oracle["celloracle_abs_axis_delta"] = pd.to_numeric(oracle["mean_abs_delta_x"], errors="coerce")
    merged = directionality_matrix.merge(
        oracle[["perturb_tf", "target_axis", "celloracle_state", "mean_delta_x", "celloracle_abs_axis_delta"]],
        on=["perturb_tf", "target_axis"],
        how="inner",
    )
    return merged


def write_table(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def save_figures(
    directionality_matrix: pd.DataFrame,
    asymmetry_tests: pd.DataFrame,
    concordance: pd.DataFrame,
    figure_dir: Path,
) -> list[str]:
    import matplotlib.pyplot as plt

    figure_dir.mkdir(parents=True, exist_ok=True)
    saved: list[str] = []

    if not directionality_matrix.empty:
        subset = "driver_union_all" if "driver_union_all" in set(directionality_matrix["subset"]) else str(directionality_matrix["subset"].iloc[0])
        heat = directionality_matrix.loc[directionality_matrix["subset"].eq(subset)].copy()
        heat = heat.pivot_table(index="perturb_tf", columns="target_axis", values="impact_score", aggfunc="mean").fillna(0.0)
        if not heat.empty:
            fig, ax = plt.subplots(figsize=(7.5, max(3.5, 0.35 * len(heat))))
            image = ax.imshow(heat.to_numpy(), aspect="auto", cmap="viridis")
            ax.set_xticks(np.arange(len(heat.columns)))
            ax.set_xticklabels(heat.columns, rotation=35, ha="right")
            ax.set_yticks(np.arange(len(heat.index)))
            ax.set_yticklabels(heat.index)
            ax.set_title(f"Module 9.2 scTenifoldKnk signature impact ({subset})")
            fig.colorbar(image, ax=ax, label="impact_score")
            fig.tight_layout()
            for ext in ["png", "pdf", "svg"]:
                out = figure_dir / f"module9_2_network_direction_heatmap.{ext}"
                fig.savefig(out, dpi=300)
                saved.append(str(out))
            plt.close(fig)

    consensus = asymmetry_tests.loc[asymmetry_tests["subset"].eq("joint_consensus")].copy() if not asymmetry_tests.empty else pd.DataFrame()
    if not consensus.empty:
        fig, ax = plt.subplots(figsize=(8, 4.2))
        values = pd.to_numeric(consensus["directionality_index"], errors="coerce").fillna(0.0)
        ax.bar(np.arange(len(consensus)), values, color="#3f6fb5")
        ax.axhline(0.0, color="black", linewidth=0.8)
        ax.set_xticks(np.arange(len(consensus)))
        ax.set_xticklabels(consensus["comparison"], rotation=30, ha="right")
        ax.set_ylabel("directionality_index")
        ax.set_title("Module 9.2 forward vs reverse asymmetry")
        fig.tight_layout()
        for ext in ["png", "pdf", "svg"]:
            out = figure_dir / f"module9_2_asymmetry_panel.{ext}"
            fig.savefig(out, dpi=300)
            saved.append(str(out))
        plt.close(fig)

    if not concordance.empty:
        plot_df = concordance.dropna(subset=["impact_score", "celloracle_abs_axis_delta"]).copy()
        if not plot_df.empty:
            fig, ax = plt.subplots(figsize=(5.2, 4.4))
            ax.scatter(plot_df["impact_score"], plot_df["celloracle_abs_axis_delta"], color="#2a9d8f", alpha=0.8)
            for _, row in plot_df.iterrows():
                ax.annotate(str(row["perturb_tf"]), (row["impact_score"], row["celloracle_abs_axis_delta"]), fontsize=7, alpha=0.7)
            ax.set_xlabel("scTenifoldKnk signature impact")
            ax.set_ylabel("CellOracle abs axis delta")
            ax.set_title("CellOracle vs scTenifoldKnk concordance")
            fig.tight_layout()
            for ext in ["png", "pdf", "svg"]:
                out = figure_dir / f"module9_2_celloracle_sctenifold_concordance.{ext}"
                fig.savefig(out, dpi=300)
                saved.append(str(out))
            plt.close(fig)

    return saved


def write_conclusions(path: Path, evidence_grade: pd.DataFrame, asymmetry_tests: pd.DataFrame, restore_audit: pd.DataFrame) -> None:
    grade = evidence_grade.iloc[0].to_dict() if not evidence_grade.empty else {}
    consensus = asymmetry_tests.loc[asymmetry_tests["subset"].eq("joint_consensus")].copy() if not asymmetry_tests.empty else pd.DataFrame()
    if consensus.empty and not asymmetry_tests.empty:
        consensus = asymmetry_tests.loc[asymmetry_tests["subset"].eq("consensus")].copy()
    lines = [
        "# Module 9.2 Network Direction Evidence",
        "",
        f"- Network direction label: `{grade.get('network_direction_label', 'not_available')}`.",
        f"- Evidence scope: {grade.get('evidence_scope', 'not_available')}.",
        f"- Restore status: `{grade.get('restore_status', 'unknown')}`.",
        "",
        "## Joint consensus asymmetry tests",
        "",
    ]
    if consensus.empty:
        lines.append("- No joint consensus asymmetry tests were available.")
    else:
        for _, row in consensus.iterrows():
            lines.append(
                "- `{comparison}`: directionality_index={directionality:.3g}, "
                "forward_impact={forward:.3g}, reverse_impact={reverse:.3g}, label=`{label}`.".format(
                    comparison=row["comparison"],
                    directionality=float(row["directionality_index"]) if pd.notna(row["directionality_index"]) else math.nan,
                    forward=float(row["forward_impact"]) if pd.notna(row["forward_impact"]) else math.nan,
                    reverse=float(row["reverse_impact"]) if pd.notna(row["reverse_impact"]) else math.nan,
                    label=row["support_label"],
                )
            )
    lines.extend(
        [
            "",
            "## Interpretation guardrails",
            "",
            "- scTenifoldKnk is treated as perturbation impact magnitude, not signed activation or repression.",
            "- CellOracle TF-axis delta is treated as signed local model shift after KO.",
            "- Restore/fixed-mode evidence is not claimed unless restore outputs are supplied.",
        ]
    )
    if not restore_audit.empty and bool(restore_audit.loc[0, "restore_available"]):
        lines.append("- Supplied restore outputs were detected and audited.")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_module(args: argparse.Namespace) -> dict[str, object]:
    start = time.time()
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)

    signatures = read_tsv_or_empty(args.tf_targets)
    signature_sets = build_signature_sets(signatures)
    sctenifold, sctenifold_paths = load_sctenifold_tables(args.metadata_dir, args.sctenifold_pattern)
    celloracle_tf_delta = read_tsv_or_empty(args.celloracle_tf_delta)
    celloracle_top_genes = read_tsv_or_empty(args.celloracle_top_genes)
    quant_scores = read_tsv_or_empty(args.celloracle_quantitative_scores)

    signature_impact = compute_sctenifold_signature_impact(sctenifold, signature_sets, fdr_threshold=args.fdr_threshold)
    directionality_matrix = build_directionality_matrix(signature_impact)
    celloracle_axis_delta = compute_celloracle_tf_axis_delta(celloracle_tf_delta, AXIS_TF_GROUPS)
    celloracle_top_overlap = compute_celloracle_top_gene_signature_overlap(celloracle_top_genes, signature_sets)
    sctenifold_asymmetry = compute_asymmetry_tests(
        signature_impact.loc[signature_impact["target_axis"].isin(CORE_TARGET_AXES)].copy(),
        evidence_source="sctenifold_signature",
    )
    celloracle_impact = compute_celloracle_axis_impact_matrix(
        celloracle_axis_delta.loc[celloracle_axis_delta["target_axis"].isin(CORE_TARGET_AXES)].copy()
    )
    celloracle_asymmetry = compute_asymmetry_tests(celloracle_impact, evidence_source="celloracle_abs_delta")
    asymmetry_tests = combine_asymmetry_sources(sctenifold_asymmetry, celloracle_asymmetry)
    restore_audit = audit_restore_availability(args.restore_celloracle_dir, require_restore=args.require_restore)
    evidence_grade = summarize_network_evidence(asymmetry_tests, restore_audit)
    concordance = merge_celloracle_sctenifold_concordance(directionality_matrix, celloracle_axis_delta)

    outputs = {
        "signature_impact": args.metadata_dir / "module9_2_sctenifold_signature_impact.tsv",
        "celloracle_axis_delta": args.metadata_dir / "module9_2_celloracle_tf_axis_delta.tsv",
        "celloracle_top_overlap": args.metadata_dir / "module9_2_celloracle_top_gene_signature_overlap.tsv",
        "directionality_matrix": args.metadata_dir / "module9_2_network_directionality_matrix.tsv",
        "celloracle_impact_matrix": args.metadata_dir / "module9_2_celloracle_impact_matrix.tsv",
        "asymmetry_tests": args.metadata_dir / "module9_2_asymmetry_tests.tsv",
        "restore_availability": args.metadata_dir / "module9_2_restore_availability.tsv",
        "evidence_grade": args.metadata_dir / "module9_2_evidence_grade.tsv",
        "concordance": args.metadata_dir / "module9_2_celloracle_sctenifold_concordance.tsv",
        "conclusions": args.metadata_dir / "module9_2_main_conclusions.md",
        "report": args.metadata_dir / "module9_2_report.json",
    }
    write_table(signature_impact, outputs["signature_impact"])
    write_table(celloracle_axis_delta, outputs["celloracle_axis_delta"])
    write_table(celloracle_top_overlap, outputs["celloracle_top_overlap"])
    write_table(directionality_matrix, outputs["directionality_matrix"])
    write_table(celloracle_impact, outputs["celloracle_impact_matrix"])
    write_table(asymmetry_tests, outputs["asymmetry_tests"])
    write_table(restore_audit, outputs["restore_availability"])
    write_table(evidence_grade, outputs["evidence_grade"])
    write_table(concordance, outputs["concordance"])
    write_conclusions(outputs["conclusions"], evidence_grade, asymmetry_tests, restore_audit)

    figure_outputs: list[str] = []
    if not args.skip_figures:
        figure_outputs = save_figures(directionality_matrix, asymmetry_tests, concordance, args.figure_dir)

    report = {
        "module": "9.2_network_direction_evidence",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "runtime_seconds": round(time.time() - start, 3),
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "package_versions": {
            "numpy": package_version("numpy"),
            "pandas": package_version("pandas"),
            "matplotlib": package_version("matplotlib"),
        },
        "inputs": {
            "tf_targets": str(args.tf_targets),
            "celloracle_tf_delta": str(args.celloracle_tf_delta),
            "celloracle_top_genes": str(args.celloracle_top_genes),
            "celloracle_quantitative_scores": str(args.celloracle_quantitative_scores),
            "sctenifold_paths": sctenifold_paths,
            "restore_celloracle_dir": str(args.restore_celloracle_dir) if args.restore_celloracle_dir is not None else "",
        },
        "counts": {
            "n_signature_rows": int(len(signatures)),
            "n_signature_sets": int(len(signature_sets)),
            "n_sctenifold_rows": int(len(sctenifold)),
            "n_sctenifold_signature_impact_rows": int(len(signature_impact)),
            "n_celloracle_tf_delta_rows": int(len(celloracle_tf_delta)),
            "n_celloracle_axis_delta_rows": int(len(celloracle_axis_delta)),
            "n_celloracle_top_overlap_rows": int(len(celloracle_top_overlap)),
            "n_celloracle_impact_rows": int(len(celloracle_impact)),
            "n_quantitative_score_rows": int(len(quant_scores)),
            "n_asymmetry_rows": int(len(asymmetry_tests)),
            "n_concordance_rows": int(len(concordance)),
        },
        "evidence_grade": evidence_grade.iloc[0].to_dict() if not evidence_grade.empty else {},
        "outputs": {key: str(value) for key, value in outputs.items()},
        "figure_outputs": figure_outputs,
    }
    outputs["report"].write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    return report


def main() -> None:
    args = parse_args()
    report = run_module(args)
    print(json.dumps({"module": report["module"], "evidence_grade": report["evidence_grade"], "outputs": report["outputs"]}, indent=2))


if __name__ == "__main__":
    main()
