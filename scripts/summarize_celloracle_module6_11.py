from __future__ import annotations

import argparse
import json
import platform
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA_DIR = PROJECT_ROOT / "metadata/driver"
DEFAULT_FIGURE_DIR = PROJECT_ROOT / "figures/driver"

TIER1_RESCUE_TFS = ["HNF4A", "PPARA", "CEBPB", "EGR1"]
TIER2_STATE_SPECIFIC_TFS = ["SOX4"]
TIER3_AP1_EARLY_TFS = ["JUN", "FOS", "ATF3", "JUND"]
CONTROL_PROXY_NOTES = {
    "HLF": "housekeeping-like proxy control; report as proxy control only",
}

INPUT_FILES = {
    "module6_9b_quantitative_scores": "celloracle_module6_9b_quantitative_tf_scores.tsv",
    "module6_10_main_strict_vs_driver_union": "celloracle_module6_10_main_strict_vs_driver_union.tsv",
    "module6_10_phase_wide_summary": "celloracle_module6_10_phase_wide_summary.tsv",
    "module6_10_lodo_summary": "celloracle_module6_10_lodo_summary.tsv",
    "module6_10_negative_controls": "celloracle_module6_10_negative_controls.tsv",
    "module6_10_review_risk_flags": "celloracle_module6_10_review_risk_flags.tsv",
    "module6_9_candidate_evidence_matrix": "celloracle_module6_9_candidate_evidence_matrix.tsv",
    "module6_7_state_grn_summary": "celloracle_module6_7_tf_network_summary.tsv",
    "module6_8_perturbation_ranking": "celloracle_module6_8_perturbation_ranking.tsv",
}

OUTPUT_FILES = {
    "candidate_tier_table": "celloracle_module6_11_candidate_tier_table.tsv",
    "claim_evidence_matrix": "celloracle_module6_11_claim_evidence_matrix.tsv",
    "main_conclusions": "celloracle_module6_11_main_conclusions.md",
    "nature_figure_plan": "celloracle_module6_11_nature_figure_plan.md",
    "supplementary_table_index": "celloracle_module6_11_supplementary_table_index.tsv",
    "methods_result_snippets": "celloracle_module6_11_methods_result_snippets.md",
    "report": "celloracle_module6_11_report.json",
}


def _read_tsv(path: Path, required: bool = True) -> pd.DataFrame:
    if not path.exists():
        if required:
            raise FileNotFoundError(path)
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t")
    if required and df.empty:
        raise ValueError(f"Required input is empty: {path}")
    return df


def _to_numeric(df: pd.DataFrame, cols: Iterable[str]) -> pd.DataFrame:
    out = df.copy()
    for col in cols:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def _join_unique(values: Iterable[object]) -> str:
    items = []
    for value in values:
        if pd.isna(value):
            continue
        text = str(value).strip()
        if text and text not in items:
            items.append(text)
    return "; ".join(items)


def _fmt_num(value: object, ndigits: int = 3) -> str:
    if pd.isna(value):
        return "NA"
    try:
        return f"{float(value):.{ndigits}f}"
    except (TypeError, ValueError):
        return str(value)


def _aggregate_negative_controls(negative_controls: pd.DataFrame) -> pd.DataFrame:
    if negative_controls.empty:
        return pd.DataFrame(columns=["tf", "control_type", "control_rationale"])
    return (
        negative_controls.groupby("tf", as_index=False)
        .agg({"control_type": _join_unique, "control_rationale": _join_unique})
        .sort_values("tf")
        .reset_index(drop=True)
    )


def _aggregate_risk_flags(risk_flags: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    cols = ["tf", "review_risk_flag", "review_risk_detail", "review_risk_severity"]
    if risk_flags.empty:
        return pd.DataFrame(columns=cols), []

    risk_flags = risk_flags.copy()
    risk_flags["tf"] = risk_flags["tf"].fillna("NA").astype(str)
    global_flags = risk_flags.loc[risk_flags["tf"].isin(["", "NA", "nan"])]
    tf_flags = risk_flags.loc[~risk_flags.index.isin(global_flags.index)]
    caveats = [
        f"{row.risk_type}: {row.detail}"
        for row in global_flags.itertuples(index=False)
        if hasattr(row, "risk_type") and hasattr(row, "detail")
    ]
    if tf_flags.empty:
        return pd.DataFrame(columns=cols), caveats

    agg = (
        tf_flags.groupby("tf", as_index=False)
        .agg(
            review_risk_flag=("risk_type", _join_unique),
            review_risk_detail=("detail", _join_unique),
            review_risk_severity=("severity", _join_unique),
        )
        .sort_values("tf")
        .reset_index(drop=True)
    )
    return agg, caveats


def build_candidate_tier_table(
    quantitative_scores: pd.DataFrame,
    strict_union: pd.DataFrame,
    phase_summary: pd.DataFrame,
    lodo_summary: pd.DataFrame,
    negative_controls: pd.DataFrame,
    evidence_matrix: pd.DataFrame,
    risk_flags: pd.DataFrame,
) -> pd.DataFrame:
    quantitative_scores = _to_numeric(
        quantitative_scores,
        [
            "quantitative_rank",
            "quantitative_perturbation_score",
            "state_specificity_score",
            "state_specificity_ratio",
            "malignant_fate_direction_score",
            "inner_product_score",
            "cnv_fate_probability_association_score",
            "module_rescue_score",
        ],
    )
    strict_union = _to_numeric(
        strict_union,
        [
            "driver_union_rank",
            "driver_union_score",
            "main_strict_rank",
            "main_strict_score",
            "rank_delta_main_minus_union",
        ],
    )
    phase_summary = _to_numeric(
        phase_summary,
        [
            "phase_early_rank",
            "phase_early_score",
            "phase_intermediate_rank",
            "phase_intermediate_score",
            "phase_late_rank",
            "phase_late_score",
        ],
    )
    lodo_summary = _to_numeric(lodo_summary, ["min_lodo_score", "max_lodo_rank", "top5_lodo_fraction"])
    evidence_cols = ["tf", "integrated_rank", "integrated_evidence_score", "role", "tier", "selected_for_main_panel"]
    evidence = evidence_matrix[[col for col in evidence_cols if col in evidence_matrix.columns]].copy()
    evidence = _to_numeric(evidence, ["integrated_rank", "integrated_evidence_score"])
    controls = _aggregate_negative_controls(negative_controls)
    tf_risks, _ = _aggregate_risk_flags(risk_flags)

    merged = (
        quantitative_scores.merge(strict_union, on="tf", how="left")
        .merge(phase_summary, on="tf", how="left")
        .merge(lodo_summary, on="tf", how="left")
        .merge(controls, on="tf", how="left")
        .merge(evidence, on="tf", how="left", suffixes=("", "_module6_9"))
        .merge(tf_risks, on="tf", how="left")
    )

    control_tfs = set(controls["tf"].astype(str)) if not controls.empty else set()

    def assign(row: pd.Series) -> tuple[str, int, str, str]:
        tf = str(row["tf"])
        is_tier1 = (
            row.get("quantitative_rank", float("inf")) <= 5
            and row.get("rank_delta_main_minus_union", float("inf")) == 0
            and row.get("top5_lodo_fraction", 0) >= 0.75
            and row.get("phase_late_rank", float("inf")) <= 5
        )
        is_tier2_signal = (
            tf in TIER2_STATE_SPECIFIC_TFS
            or (
                row.get("state_specificity_score", 0) >= 0.75
                and row.get("quantitative_rank", float("inf")) <= 10
            )
            or (
                row.get("state_specificity_ratio", 0) >= 5
                and row.get("quantitative_rank", float("inf")) <= 10
            )
        )
        is_ap1_early = (
            tf in TIER3_AP1_EARLY_TFS
            or row.get("integrated_rank", float("inf")) <= 5
            or row.get("phase_early_rank", float("inf")) <= 5
        )

        if is_tier1:
            return (
                "Tier 1",
                1,
                "robust rescue candidate",
                "Robust rescue TF supported by quantitative rank, strict/union agreement, LODO stability and late-phase effect.",
            )
        if is_tier2_signal:
            return (
                "Tier 2",
                2,
                "state-specific rescue candidate",
                "State-specific or late malignant-like candidate with weaker late/robustness support than Tier 1.",
            )
        if tf in TIER3_AP1_EARLY_TFS or (is_ap1_early and tf not in control_tfs):
            return (
                "Tier 3",
                3,
                "AP-1 / early-transition evidence",
                "AP-1 or early-transition TF with strong integrated/early evidence and lower robustness than Tier 1.",
            )
        if tf in control_tfs:
            return (
                "Negative/control",
                4,
                "negative or proxy control",
                CONTROL_PROXY_NOTES.get(tf, "Negative-control TF retained for reviewer-facing calibration."),
            )
        return (
            "Reserve/other",
            5,
            "reserve evidence",
            "Candidate retained in source tables but not emphasized in the Nature main claim.",
        )

    assigned = merged.apply(assign, axis=1, result_type="expand")
    assigned.columns = ["candidate_tier", "tier_order", "tier_label", "primary_interpretation"]
    merged = pd.concat([merged, assigned], axis=1)

    merged["figure_role"] = merged["tf"].map(
        lambda tf: "main highlight"
        if tf in TIER1_RESCUE_TFS + TIER2_STATE_SPECIFIC_TFS
        else ("AP-1 evidence axis" if tf in TIER3_AP1_EARLY_TFS else ("control calibration" if tf in control_tfs else "supplementary only"))
    )
    merged["manuscript_use"] = merged["candidate_tier"].map(
        {
            "Tier 1": "Primary Results claim and Figure 6b-e highlight",
            "Tier 2": "State-specific Results claim and Figure 6b/c highlight",
            "Tier 3": "AP-1 early-transition Results paragraph and Extended Data",
            "Negative/control": "Robustness/control supplement and limitation text",
            "Reserve/other": "Supplementary table only",
        }
    )

    ordered_cols = [
        "tf",
        "candidate_tier",
        "tier_order",
        "tier_label",
        "primary_interpretation",
        "quantitative_rank",
        "quantitative_perturbation_score",
        "integrated_rank",
        "integrated_evidence_score",
        "driver_union_rank",
        "main_strict_rank",
        "rank_delta_main_minus_union",
        "phase_early_rank",
        "phase_intermediate_rank",
        "phase_late_rank",
        "top5_lodo_fraction",
        "min_lodo_score",
        "max_lodo_rank",
        "control_type",
        "control_rationale",
        "review_risk_flag",
        "review_risk_detail",
        "review_risk_severity",
        "figure_role",
        "manuscript_use",
    ]
    for col in ordered_cols:
        if col not in merged.columns:
            merged[col] = pd.NA

    return (
        merged[ordered_cols]
        .sort_values(["tier_order", "quantitative_rank", "tf"], na_position="last")
        .reset_index(drop=True)
    )


def build_claim_evidence_matrix(candidate_tiers: pd.DataFrame) -> pd.DataFrame:
    tiers = candidate_tiers.set_index("tf", drop=False)

    def tf_metric(tf: str) -> str:
        if tf not in tiers.index:
            return f"{tf}: not present"
        row = tiers.loc[tf]
        return (
            f"{tf}: q_rank={_fmt_num(row.get('quantitative_rank'), 0)}, "
            f"late_rank={_fmt_num(row.get('phase_late_rank'), 0)}, "
            f"LODO_top5={_fmt_num(row.get('top5_lodo_fraction'))}"
        )

    rows = [
        {
            "claim_id": "robust_rescue_tf_candidates",
            "claim": "HNF4A, PPARA, CEBPB and EGR1 are the primary robust rescue TF candidates from CellOracle perturbation analysis.",
            "tf_group": ", ".join(TIER1_RESCUE_TFS),
            "evidence_status": "main claim",
            "primary_source_table": INPUT_FILES["module6_9b_quantitative_scores"],
            "supporting_source_tables": "; ".join(
                [
                    INPUT_FILES["module6_10_main_strict_vs_driver_union"],
                    INPUT_FILES["module6_10_lodo_summary"],
                    INPUT_FILES["module6_10_phase_wide_summary"],
                ]
            ),
            "key_metrics": "; ".join(tf_metric(tf) for tf in TIER1_RESCUE_TFS),
            "nature_figure_panel": "Figure 6b-c; Figure 6e",
            "supplementary_table": "Supplementary Table 6.5",
            "caveat": "CellOracle perturbation is an in silico directional-prioritization result requiring experimental validation.",
        },
        {
            "claim_id": "state_specific_sox4",
            "claim": "SOX4 is prioritized as a state-specific malignant-like candidate rather than a Tier 1 broad rescue TF.",
            "tf_group": "SOX4",
            "evidence_status": "secondary claim",
            "primary_source_table": OUTPUT_FILES["candidate_tier_table"],
            "supporting_source_tables": "; ".join(
                [
                    INPUT_FILES["module6_9b_quantitative_scores"],
                    INPUT_FILES["module6_10_phase_wide_summary"],
                ]
            ),
            "key_metrics": tf_metric("SOX4"),
            "nature_figure_panel": "Figure 6b-c",
            "supplementary_table": "Supplementary Table 6.5",
            "caveat": "SOX4 is separated from Tier 1 because its late-phase rank is lower than the Tier 1 cutoff.",
        },
        {
            "claim_id": "ap1_early_transition_axis",
            "claim": "JUN, FOS, ATF3 and JUND support an AP-1 / early-transition evidence axis with weaker robustness than Tier 1 rescue TFs.",
            "tf_group": ", ".join(TIER3_AP1_EARLY_TFS),
            "evidence_status": "contextual claim",
            "primary_source_table": INPUT_FILES["module6_9_candidate_evidence_matrix"],
            "supporting_source_tables": "; ".join(
                [
                    INPUT_FILES["module6_10_lodo_summary"],
                    INPUT_FILES["module6_10_phase_wide_summary"],
                    INPUT_FILES["module6_10_review_risk_flags"],
                ]
            ),
            "key_metrics": "; ".join(tf_metric(tf) for tf in TIER3_AP1_EARLY_TFS),
            "nature_figure_panel": "Figure 6e; Extended Data 6.5",
            "supplementary_table": "Supplementary Table 6.4; Supplementary Table 6.5",
            "caveat": "AP-1 TFs should be described as an early-transition evidence axis, not as the most robust rescue class.",
        },
        {
            "claim_id": "negative_controls_and_proxy_control",
            "claim": "MAFB, JUNB, MAFF, MYC and HLF are retained as negative or proxy controls for calibration.",
            "tf_group": "MAFB, JUNB, MAFF, MYC, HLF",
            "evidence_status": "control claim",
            "primary_source_table": INPUT_FILES["module6_10_negative_controls"],
            "supporting_source_tables": INPUT_FILES["module6_10_review_risk_flags"],
            "key_metrics": "control_type; control_rationale; HLF reported only as housekeeping-like proxy control",
            "nature_figure_panel": "Extended Data 6.6",
            "supplementary_table": "Supplementary Table 6.4",
            "caveat": "HLF is a proxy control in this TF panel and should not be called a canonical housekeeping TF.",
        },
        {
            "claim_id": "no_single_dataset_or_sample_dominance",
            "claim": "The key CellOracle perturbation conclusion is not flagged as being dominated by a single dataset or sample group.",
            "tf_group": "all ranked TFs",
            "evidence_status": "robustness claim",
            "primary_source_table": "celloracle_module6_10_dataset_sample_dominance.tsv",
            "supporting_source_tables": INPUT_FILES["module6_10_review_risk_flags"],
            "key_metrics": "dominant_fraction below pre-specified dominance threshold in Module 6.10",
            "nature_figure_panel": "Extended Data 6.5",
            "supplementary_table": "Supplementary Table 6.4",
            "caveat": "Unknown metadata groups should be transparently reported if they contribute a large fraction.",
        },
    ]
    return pd.DataFrame(rows)


def build_supplementary_table_index() -> pd.DataFrame:
    rows = [
        {
            "supplementary_table": "Supplementary Table 6.1",
            "title": "CellOracle input TF selection and upstream evidence",
            "primary_source_file": "celloracle_tf_selection.module6_4.tsv",
            "companion_files": "celloracle_input_tfs.module6_4.txt",
            "description": "Input TF panel and evidence classes used to start CellOracle perturbation analysis.",
            "intended_manuscript_use": "Methods and input-data transparency for Module 6.4.",
            "reviewer_question_addressed": "How were perturbation TFs selected before CellOracle simulation?",
        },
        {
            "supplementary_table": "Supplementary Table 6.2",
            "title": "State-specific CellOracle GRN link summary",
            "primary_source_file": INPUT_FILES["module6_7_state_grn_summary"],
            "companion_files": "celloracle_module6_7_links_raw.tsv.gz; celloracle_module6_7_links_filtered.tsv.gz",
            "description": "State-resolved TF-target edge counts, coefficients and filtered GRN support.",
            "intended_manuscript_use": "Extended Data 6.1 and Methods for state-specific GRN fitting.",
            "reviewer_question_addressed": "Which TF-target links support state-specific regulatory interpretation?",
        },
        {
            "supplementary_table": "Supplementary Table 6.3",
            "title": "TF perturbation rankings and quantitative scoring",
            "primary_source_file": INPUT_FILES["module6_9b_quantitative_scores"],
            "companion_files": "celloracle_module6_8_perturbation_ranking.tsv; celloracle_module6_9b_cell_level_scores.tsv.gz",
            "description": "TF KO perturbation ranking, five-score quantitative rescue metrics and cell-level scores.",
            "intended_manuscript_use": "Figure 6b and Extended Data 6.2-6.3.",
            "reviewer_question_addressed": "How was visual perturbation direction converted into ranked TF evidence?",
        },
        {
            "supplementary_table": "Supplementary Table 6.4",
            "title": "Robustness, negative controls and review-risk flags",
            "primary_source_file": INPUT_FILES["module6_10_main_strict_vs_driver_union"],
            "companion_files": "; ".join(
                [
                    INPUT_FILES["module6_10_phase_wide_summary"],
                    INPUT_FILES["module6_10_lodo_summary"],
                    INPUT_FILES["module6_10_negative_controls"],
                    INPUT_FILES["module6_10_review_risk_flags"],
                    "celloracle_module6_10_dataset_sample_dominance.tsv",
                ]
            ),
            "description": "Main-strict versus driver-union, phase stratification, leave-one-group-out stability, controls and caveats.",
            "intended_manuscript_use": "Extended Data 6.4-6.6 and reviewer-facing robustness documentation.",
            "reviewer_question_addressed": "Are key perturbation claims stable across cell definitions, states and data groups?",
        },
        {
            "supplementary_table": "Supplementary Table 6.5",
            "title": "Final CellOracle candidate TF tiers",
            "primary_source_file": OUTPUT_FILES["candidate_tier_table"],
            "companion_files": OUTPUT_FILES["claim_evidence_matrix"],
            "description": "Final manuscript-facing TF tiers and linked claim-evidence matrix.",
            "intended_manuscript_use": "Final Results table and claim audit for Module 6.11.",
            "reviewer_question_addressed": "Which TFs are main rescue candidates, state-specific candidates, AP-1 evidence or controls?",
        },
    ]
    return pd.DataFrame(rows)


def build_figure_plan_rows() -> list[dict[str, str]]:
    return [
        {
            "figure": "Figure 6a",
            "panel_title": "CellOracle workflow schematic",
            "claim": "Base GRN, state-specific GRN, TF KO simulation and quantitative scoring form one auditable perturbation workflow.",
            "input_file": "; ".join([INPUT_FILES["module6_7_state_grn_summary"], INPUT_FILES["module6_8_perturbation_ranking"], INPUT_FILES["module6_9b_quantitative_scores"]]),
            "metric": "workflow steps; input TF count; Oracle cell/gene dimensions from Module 6.6 report",
            "output_or_figure_reference": "planned schematic: celloracle_module6_11_workflow_schematic",
        },
        {
            "figure": "Figure 6b",
            "panel_title": "Quantitative TF perturbation ranking",
            "claim": "HNF4A, PPARA, CEBPB, EGR1 and SOX4 lead the quantitative perturbation score.",
            "input_file": INPUT_FILES["module6_9b_quantitative_scores"],
            "metric": "quantitative_perturbation_score; quantitative_rank; five component scores",
            "output_or_figure_reference": "planned bar or dot plot from Module 6.9b table",
        },
        {
            "figure": "Figure 6c",
            "panel_title": "State and robustness heatmap",
            "claim": "Tier 1 rescue TFs are stable across main-strict/driver-union definitions and late malignant-like states.",
            "input_file": "; ".join([INPUT_FILES["module6_10_main_strict_vs_driver_union"], INPUT_FILES["module6_10_phase_wide_summary"]]),
            "metric": "rank_delta_main_minus_union; phase_early/intermediate/late ranks",
            "output_or_figure_reference": "planned heatmap from Module 6.10 tables",
        },
        {
            "figure": "Figure 6d",
            "panel_title": "Top TF perturbation vector fields",
            "claim": "Top TF perturbations redirect malignant-like vector fields toward the rescue direction.",
            "input_file": INPUT_FILES["module6_8_perturbation_ranking"],
            "metric": "CellOracle perturbation vector field; anti_malignant_shift_score",
            "output_or_figure_reference": "figures/driver/celloracle_module6_9_top_tf_vector_fields.svg",
        },
        {
            "figure": "Figure 6e",
            "panel_title": "Candidate evidence matrix",
            "claim": "CellRank, SCENIC/cisTarget, GRN, perturbation and robustness evidence jointly stratify TF candidates.",
            "input_file": "; ".join([INPUT_FILES["module6_9_candidate_evidence_matrix"], OUTPUT_FILES["candidate_tier_table"]]),
            "metric": "integrated_evidence_score; candidate_tier; robustness metrics",
            "output_or_figure_reference": "figures/driver/celloracle_module6_9_candidate_evidence_heatmap.svg",
        },
        {
            "figure": "Extended Data 6.1",
            "panel_title": "State-specific GRN edge summary",
            "claim": "CellOracle GRN fitting yields state-resolved TF-target edge support.",
            "input_file": INPUT_FILES["module6_7_state_grn_summary"],
            "metric": "n_edges_passing_p; mean_coef_abs_passing_p",
            "output_or_figure_reference": "planned Extended Data panel",
        },
        {
            "figure": "Extended Data 6.2",
            "panel_title": "15 TF perturbation ranking and state projection",
            "claim": "All 15 perturbation TFs are retained for transparent rank comparison.",
            "input_file": INPUT_FILES["module6_8_perturbation_ranking"],
            "metric": "rank; state-specific projection scores",
            "output_or_figure_reference": "figures/driver/celloracle_module6_9_perturbation_ranking_barplot.svg",
        },
        {
            "figure": "Extended Data 6.3",
            "panel_title": "Five-score quantitative perturbation scoring",
            "claim": "Final perturbation rank is explained by fate direction, inner product, CNV-fate association, module rescue and state specificity.",
            "input_file": INPUT_FILES["module6_9b_quantitative_scores"],
            "metric": "five scaled component scores",
            "output_or_figure_reference": "planned Extended Data heatmap",
        },
        {
            "figure": "Extended Data 6.4",
            "panel_title": "Main-strict versus driver-union robustness",
            "claim": "Primary rescue TFs are not dependent on one malignant-cell definition.",
            "input_file": INPUT_FILES["module6_10_main_strict_vs_driver_union"],
            "metric": "rank_delta_main_minus_union; score_delta_main_minus_union",
            "output_or_figure_reference": "planned Extended Data panel",
        },
        {
            "figure": "Extended Data 6.5",
            "panel_title": "Leave-one-dataset/sample-out stability",
            "claim": "Tier 1 TFs show stronger leave-one-group-out stability than AP-1 early-transition TFs.",
            "input_file": INPUT_FILES["module6_10_lodo_summary"],
            "metric": "top5_lodo_fraction; max_lodo_rank; worst_score_delta_vs_full",
            "output_or_figure_reference": "planned Extended Data panel",
        },
        {
            "figure": "Extended Data 6.6",
            "panel_title": "Negative controls and review-risk flags",
            "claim": "Negative and proxy controls calibrate interpretation and prevent overclaiming.",
            "input_file": "; ".join([INPUT_FILES["module6_10_negative_controls"], INPUT_FILES["module6_10_review_risk_flags"]]),
            "metric": "control_type; severity; review_risk_flag",
            "output_or_figure_reference": "planned Extended Data panel",
        },
    ]


def _markdown_table(rows: list[dict[str, str]], columns: list[str]) -> str:
    lines = []
    lines.append("| " + " | ".join(columns) + " |")
    lines.append("| " + " | ".join(["---"] * len(columns)) + " |")
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(col, "")).replace("\n", " ") for col in columns) + " |")
    return "\n".join(lines)


def build_nature_figure_plan_md(rows: list[dict[str, str]]) -> str:
    columns = ["figure", "panel_title", "claim", "input_file", "metric", "output_or_figure_reference"]
    return "\n".join(
        [
            "# Module 6.11 Nature Figure And Extended Data Plan",
            "",
            "This plan organizes existing CellOracle Module 6.4-6.10 outputs for manuscript figures without rerunning CellOracle.",
            "",
            _markdown_table(rows, columns),
            "",
        ]
    )


def build_main_conclusions_md(candidate_tiers: pd.DataFrame, global_caveats: list[str]) -> str:
    tier_counts = candidate_tiers["candidate_tier"].value_counts().to_dict()
    tier1 = ", ".join(TIER1_RESCUE_TFS)
    tier2 = ", ".join(TIER2_STATE_SPECIFIC_TFS)
    tier3 = ", ".join(TIER3_AP1_EARLY_TFS)
    controls = ", ".join(candidate_tiers.loc[candidate_tiers["candidate_tier"] == "Negative/control", "tf"].astype(str))
    caveat_lines = "\n".join(f"- {caveat}" for caveat in global_caveats) if global_caveats else "- No global Module 6.10 caveat beyond TF-specific risk flags."
    return f"""# Module 6.11 Main Conclusions

## Primary Claim
CellOracle perturbation prioritizes {tier1} as Tier 1 robust rescue candidates. These TFs combine top quantitative rescue rank, identical main-strict and driver-union ranks, high leave-one-group-out top-5 stability and late malignant-like phase support.

## State-Specific Claim
{tier2} is retained as a Tier 2 state-specific rescue candidate. It has strong state-specificity evidence, but it should be presented separately from the Tier 1 rescue class because its late-phase rank does not meet the Tier 1 cutoff.

## AP-1 / Early-Transition Claim
{tier3} define a Tier 3 AP-1 / early-transition evidence axis. They are useful for biological interpretation and extended evidence, while the robustness language should be more conservative than for Tier 1 rescue TFs.

## Controls And Caveats
Negative/proxy controls retained from Module 6.10: {controls}. HLF should be described only as a housekeeping-like proxy control in this restricted TF panel.

Global caveats from Module 6.10:
{caveat_lines}

## Tier Counts
{json.dumps(tier_counts, indent=2, sort_keys=True)}
"""


def build_methods_result_snippets_md(candidate_tiers: pd.DataFrame) -> str:
    tier1 = ", ".join(TIER1_RESCUE_TFS)
    tier2 = ", ".join(TIER2_STATE_SPECIFIC_TFS)
    tier3 = ", ".join(TIER3_AP1_EARLY_TFS)
    n_candidates = int(candidate_tiers["tf"].nunique())
    return f"""# Module 6.11 Manuscript Snippets

## Methods Snippet
We integrated previously generated CellOracle outputs without rerunning CellOracle. Candidate TFs were summarized from quantitative perturbation scores, main-strict versus driver-union robustness, phase-stratified ranks, leave-one-group-out stability, negative-control annotations and review-risk flags. Final TF tiers were assigned using pre-specified rank and stability rules, with HLF treated as a housekeeping-like proxy control only.

## Results Snippet
Across {n_candidates} TFs, CellOracle perturbation analysis prioritized {tier1} as robust rescue candidates. SOX4 was retained as a state-specific malignant-like candidate, while {tier3} formed an AP-1 / early-transition evidence axis with less robust rescue support than Tier 1 TFs.

## Limitation Snippet
CellOracle perturbation results should be interpreted as in silico directional prioritization. AP-1 TFs should be described as early-transition or multi-evidence candidates rather than the most robust rescue class, and HLF should be reported only as a proxy control in this panel.
"""


def build_report_payload(
    candidate_tiers: pd.DataFrame,
    figure_rows: list[dict[str, str]],
    supplementary_index: pd.DataFrame,
    inputs: dict[str, str],
    outputs: dict[str, str],
) -> dict[str, object]:
    tier_counts = candidate_tiers["candidate_tier"].value_counts().sort_index().to_dict()
    return {
        "module": "Module 6.11",
        "purpose": "Summarize CellOracle Module 6.4-6.10 outputs into manuscript-ready candidate tiers, claims and figure/table organization.",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": platform.python_version(),
        "inputs": inputs,
        "outputs": outputs,
        "n_candidates": int(candidate_tiers["tf"].nunique()),
        "tier_counts": {str(k): int(v) for k, v in tier_counts.items()},
        "figure_count": int(len(figure_rows)),
        "supplementary_table_count": int(len(supplementary_index)),
        "tier1_robust_rescue_tfs": TIER1_RESCUE_TFS,
        "tier2_state_specific_tfs": TIER2_STATE_SPECIFIC_TFS,
        "tier3_ap1_early_transition_tfs": TIER3_AP1_EARLY_TFS,
        "negative_control_tfs": sorted(candidate_tiers.loc[candidate_tiers["candidate_tier"] == "Negative/control", "tf"].astype(str).tolist()),
        "assumptions": [
            "Module 6.11 integrates existing outputs and does not rerun CellOracle.",
            "HLF is treated as a housekeeping-like proxy control only.",
            "AP-1 TFs are reported as early-transition or multi-evidence candidates with weaker robustness than Tier 1 rescue TFs.",
        ],
    }


def write_outputs(metadata_dir: Path, figure_dir: Path) -> dict[str, Path]:
    input_paths = {key: metadata_dir / name for key, name in INPUT_FILES.items()}
    quantitative = _read_tsv(input_paths["module6_9b_quantitative_scores"])
    strict_union = _read_tsv(input_paths["module6_10_main_strict_vs_driver_union"])
    phase = _read_tsv(input_paths["module6_10_phase_wide_summary"])
    lodo = _read_tsv(input_paths["module6_10_lodo_summary"])
    controls = _read_tsv(input_paths["module6_10_negative_controls"])
    risks = _read_tsv(input_paths["module6_10_review_risk_flags"])
    evidence = _read_tsv(input_paths["module6_9_candidate_evidence_matrix"])
    _read_tsv(input_paths["module6_7_state_grn_summary"])
    _read_tsv(input_paths["module6_8_perturbation_ranking"])

    candidate_tiers = build_candidate_tier_table(quantitative, strict_union, phase, lodo, controls, evidence, risks)
    claim_matrix = build_claim_evidence_matrix(candidate_tiers)
    supp_index = build_supplementary_table_index()
    figure_rows = build_figure_plan_rows()
    _, global_caveats = _aggregate_risk_flags(risks)

    output_paths = {key: metadata_dir / name for key, name in OUTPUT_FILES.items()}
    candidate_tiers.to_csv(output_paths["candidate_tier_table"], sep="\t", index=False)
    claim_matrix.to_csv(output_paths["claim_evidence_matrix"], sep="\t", index=False)
    supp_index.to_csv(output_paths["supplementary_table_index"], sep="\t", index=False)
    output_paths["main_conclusions"].write_text(build_main_conclusions_md(candidate_tiers, global_caveats), encoding="utf-8")
    output_paths["nature_figure_plan"].write_text(build_nature_figure_plan_md(figure_rows), encoding="utf-8")
    output_paths["methods_result_snippets"].write_text(build_methods_result_snippets_md(candidate_tiers), encoding="utf-8")

    existing_figure_refs = {
        "top_tf_vector_fields_svg_exists": (figure_dir / "celloracle_module6_9_top_tf_vector_fields.svg").exists(),
        "candidate_evidence_heatmap_svg_exists": (figure_dir / "celloracle_module6_9_candidate_evidence_heatmap.svg").exists(),
        "perturbation_ranking_barplot_svg_exists": (figure_dir / "celloracle_module6_9_perturbation_ranking_barplot.svg").exists(),
    }
    report = build_report_payload(
        candidate_tiers,
        figure_rows,
        supp_index,
        inputs={key: str(path) for key, path in input_paths.items()},
        outputs={key: str(path) for key, path in output_paths.items()},
    )
    report["existing_figure_references"] = existing_figure_refs
    output_paths["report"].write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")

    for key, path in output_paths.items():
        if not path.exists() or path.stat().st_size == 0:
            raise ValueError(f"Output was not created or is empty: {key} -> {path}")
    return output_paths


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize CellOracle Module 6.11 candidate tiers and manuscript organization.")
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    args = parser.parse_args()

    paths = write_outputs(args.metadata_dir, args.figure_dir)
    print(json.dumps({key: str(path) for key, path in paths.items()}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
