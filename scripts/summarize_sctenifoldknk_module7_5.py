from __future__ import annotations

import argparse
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA_DIR = PROJECT_ROOT / "metadata/driver"
DEFAULT_FIGURE_DIR = PROJECT_ROOT / "figures/driver"

MAIN_TIER1_TFS = {"HNF4A", "PPARA", "CEBPB", "EGR1"}
MAIN_STATE_TFS = {"SOX4"}
MAIN_AP1_TFS = {"ATF3", "FOS", "JUN", "JUND"}
SUPPLEMENT_CONTROL_TFS = {"HLF", "MYC", "MAFF", "JUNB", "MAFB"}
DISPLAY_GROUP_ORDER = {
    "main_tier1": 1,
    "main_state_specific": 2,
    "main_ap1_axis": 3,
    "supplement_control": 4,
    "supplement_reserve": 5,
}


def read_tsv_or_empty(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def assign_module7_display_group(tf: str) -> str:
    tf_upper = str(tf).upper()
    if tf_upper in MAIN_TIER1_TFS:
        return "main_tier1"
    if tf_upper in MAIN_STATE_TFS:
        return "main_state_specific"
    if tf_upper in MAIN_AP1_TFS:
        return "main_ap1_axis"
    if tf_upper in SUPPLEMENT_CONTROL_TFS:
        return "supplement_control"
    return "supplement_reserve"


def manuscript_use_for_display_group(display_group: str) -> str:
    return {
        "main_tier1": "main_text_tier1_replication",
        "main_state_specific": "main_text_state_specific_replication",
        "main_ap1_axis": "main_text_ap1_axis",
        "supplement_control": "supplementary_control_calibration",
        "supplement_reserve": "supplementary_reserve_context",
    }.get(display_group, "supplementary_reserve_context")


def build_tf_level_replication_matrix(
    integrated: pd.DataFrame,
    biological_axis: pd.DataFrame,
    state_specific: pd.DataFrame,
) -> pd.DataFrame:
    if integrated.empty:
        return pd.DataFrame()

    matrix = integrated.copy()
    matrix["tf"] = matrix["tf"].astype(str)

    axis_cols = [
        "tf",
        "hcc_marker_count",
        "hcc_marker_genes",
        "stress_marker_count",
        "stress_marker_genes",
        "proliferation_marker_count",
        "proliferation_marker_genes",
        "module7_axis_interpretation",
    ]
    if not biological_axis.empty:
        available_axis_cols = [col for col in axis_cols if col in biological_axis.columns]
        if "tf" in available_axis_cols:
            matrix = matrix.merge(biological_axis[available_axis_cols], on="tf", how="left")

    if not state_specific.empty and {"tf", "gene"}.issubset(state_specific.columns):
        state_counts = (
            state_specific.assign(tf=state_specific["tf"].astype(str))
            .groupby("tf")["gene"]
            .nunique()
            .reset_index(name="malignant_like_state_specific_gene_count")
        )
    else:
        state_counts = pd.DataFrame(columns=["tf", "malignant_like_state_specific_gene_count"])
    matrix = matrix.merge(state_counts, on="tf", how="left")

    matrix["malignant_like_state_specific_gene_count"] = (
        pd.to_numeric(matrix["malignant_like_state_specific_gene_count"], errors="coerce").fillna(0).astype(int)
    )
    for col in ["hcc_marker_count", "stress_marker_count", "proliferation_marker_count"]:
        if col in matrix.columns:
            matrix[col] = pd.to_numeric(matrix[col], errors="coerce").fillna(0).astype(int)

    matrix["display_group"] = matrix["tf"].map(assign_module7_display_group)
    matrix["display_order"] = matrix["display_group"].map(DISPLAY_GROUP_ORDER).fillna(99).astype(int)
    matrix["main_text_priority"] = matrix["display_group"].str.startswith("main_")
    matrix["manuscript_use"] = matrix["display_group"].map(manuscript_use_for_display_group)
    matrix["manuscript_panel"] = matrix["display_group"].map(
        {
            "main_tier1": "Figure 7A TF-level replication matrix",
            "main_state_specific": "Figure 7A state-specific replication",
            "main_ap1_axis": "Figure 7B AP-1 concordance and pathway panel",
            "supplement_control": "Supplementary Figure 7 control calibration",
            "supplement_reserve": "Supplementary Table 7 reserve context",
        }
    )

    sort_cols = ["display_order"]
    for col in ["module7_integrated_rank", "quantitative_rank", "tf"]:
        if col in matrix.columns:
            sort_cols.append(col)
    matrix = matrix.sort_values(sort_cols, kind="mergesort").reset_index(drop=True)

    preferred_cols = [
        "tf",
        "display_group",
        "main_text_priority",
        "manuscript_use",
        "manuscript_panel",
        "candidate_tier",
        "quantitative_rank",
        "quantitative_perturbation_score",
        "scTenifoldKnk_rank",
        "scTenifoldKnk_score",
        "module7_integrated_rank",
        "integrated_module7_score",
        "n_significant_perturbed_genes",
        "mean_distance_significant",
        "n_celloracle_grn_targets",
        "n_grn_overlap_genes",
        "top_gene_grn_target_jaccard",
        "hcc_marker_count",
        "hcc_marker_genes",
        "stress_marker_count",
        "stress_marker_genes",
        "proliferation_marker_count",
        "proliferation_marker_genes",
        "malignant_like_state_specific_gene_count",
        "module7_axis_interpretation",
        "primary_interpretation",
        "review_risk_flag",
        "review_risk_detail",
    ]
    ordered_cols = [col for col in preferred_cols if col in matrix.columns]
    remaining_cols = [col for col in matrix.columns if col not in ordered_cols and col != "display_order"]
    return matrix[ordered_cols + remaining_cols]


def build_pathway_level_enrichment_matrix(enrichment_summary: pd.DataFrame, tf_matrix: pd.DataFrame) -> pd.DataFrame:
    if enrichment_summary.empty:
        return pd.DataFrame()

    matrix = enrichment_summary.copy()
    matrix["tf"] = matrix["tf"].astype(str)
    tf_meta_cols = [col for col in ["tf", "display_group", "candidate_tier", "manuscript_panel"] if col in tf_matrix.columns]
    if tf_meta_cols:
        matrix = matrix.merge(tf_matrix[tf_meta_cols].drop_duplicates("tf"), on="tf", how="left")
    if "display_group" not in matrix.columns:
        matrix["display_group"] = matrix["tf"].map(assign_module7_display_group)
    else:
        matrix["display_group"] = matrix["display_group"].fillna(matrix["tf"].map(assign_module7_display_group))

    matrix["display_order"] = matrix["display_group"].map(DISPLAY_GROUP_ORDER).fillna(99).astype(int)
    matrix["p.adjust_sort"] = pd.to_numeric(matrix.get("p.adjust", 1.0), errors="coerce").fillna(1.0)
    matrix["abs_NES_sort"] = pd.to_numeric(matrix.get("NES", 0.0), errors="coerce").abs().fillna(0.0)
    matrix = matrix.sort_values(
        ["display_order", "tf", "analysis", "database", "p.adjust_sort", "abs_NES_sort"],
        ascending=[True, True, True, True, True, False],
        kind="mergesort",
    )
    matrix["pathway_rank_within_tf_database"] = matrix.groupby(["tf", "analysis", "database"]).cumcount() + 1
    matrix["manuscript_use"] = matrix["display_group"].apply(
        lambda group: "main_pathway_panel" if str(group).startswith("main_") else "supplementary_pathway_calibration"
    )
    matrix["pathway_focus_axis"] = matrix.apply(classify_pathway_focus_axis, axis=1)

    drop_cols = ["display_order", "p.adjust_sort", "abs_NES_sort"]
    return matrix.drop(columns=[col for col in drop_cols if col in matrix.columns])


def classify_pathway_focus_axis(row: pd.Series) -> str:
    group = str(row.get("display_group", ""))
    term = str(row.get("term_name", "")).lower()
    if group == "main_ap1_axis":
        if any(token in term for token in ["heat", "stress", "response", "prolifer", "cell cycle", "apoptosis"]):
            return "AP-1 stress/proliferation"
        return "AP-1 broader pathway context"
    if group == "main_tier1":
        if any(token in term for token in ["metabolism", "fatty", "bile", "oxid", "ppar", "hepato", "ketone"]):
            return "Tier 1 liver/metabolic context"
        return "Tier 1 replication context"
    if group == "main_state_specific":
        return "SOX4 malignant-like state context"
    return "control/calibration context"


def plot_celloracle_sctenifold_concordance(tf_matrix: pd.DataFrame, figure_dir: Path) -> dict:
    if tf_matrix.empty:
        return {}
    required = {"tf", "quantitative_perturbation_score", "scTenifoldKnk_score", "display_group"}
    if not required.issubset(tf_matrix.columns):
        return {}

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure_dir.mkdir(parents=True, exist_ok=True)
    plot_df = tf_matrix.copy()
    plot_df["quantitative_perturbation_score"] = pd.to_numeric(
        plot_df["quantitative_perturbation_score"], errors="coerce"
    )
    plot_df["scTenifoldKnk_score"] = pd.to_numeric(plot_df["scTenifoldKnk_score"], errors="coerce")
    plot_df["n_significant_perturbed_genes"] = pd.to_numeric(
        plot_df.get("n_significant_perturbed_genes", 20), errors="coerce"
    ).fillna(20)
    plot_df = plot_df.dropna(subset=["quantitative_perturbation_score", "scTenifoldKnk_score"])
    if plot_df.empty:
        return {}

    palette = {
        "main_tier1": "#1f77b4",
        "main_state_specific": "#2ca02c",
        "main_ap1_axis": "#d62728",
        "supplement_control": "#7f7f7f",
        "supplement_reserve": "#9467bd",
    }
    labels = {
        "main_tier1": "Tier 1 focus",
        "main_state_specific": "SOX4 state-specific",
        "main_ap1_axis": "AP-1 axis",
        "supplement_control": "Control/calibration",
        "supplement_reserve": "Reserve",
    }
    min_count = float(plot_df["n_significant_perturbed_genes"].min())
    max_count = float(plot_df["n_significant_perturbed_genes"].max())
    denom = max(max_count - min_count, 1.0)
    plot_df["point_size"] = 70 + (plot_df["n_significant_perturbed_genes"] - min_count) / denom * 230

    fig, ax = plt.subplots(figsize=(7.2, 5.2), dpi=160)
    for group, group_df in plot_df.groupby("display_group", sort=False):
        ax.scatter(
            group_df["quantitative_perturbation_score"],
            group_df["scTenifoldKnk_score"],
            s=group_df["point_size"],
            c=palette.get(group, "#333333"),
            label=labels.get(group, group),
            alpha=0.82,
            edgecolor="white",
            linewidth=0.8,
        )

    for _, row in plot_df.iterrows():
        ax.annotate(
            str(row["tf"]),
            (row["quantitative_perturbation_score"], row["scTenifoldKnk_score"]),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=8,
        )

    corr = plot_df["quantitative_perturbation_score"].corr(plot_df["scTenifoldKnk_score"], method="spearman")
    ax.set_xlabel("CellOracle quantitative perturbation score")
    ax.set_ylabel("scTenifoldKnk replication score")
    ax.set_title("CellOracle vs scTenifoldKnk TF concordance")
    ax.text(
        0.02,
        0.92,
        f"Spearman rho = {corr:.2f}" if pd.notna(corr) else "Spearman rho = NA",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 2},
    )
    ax.grid(True, color="#dddddd", linewidth=0.6, alpha=0.7)
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    fig.tight_layout()

    output_stem = figure_dir / "sctenifoldknk_module7_5_celloracle_vs_sctenifold_concordance"
    outputs = {
        "concordance_png": str(output_stem.with_suffix(".png")),
        "concordance_pdf": str(output_stem.with_suffix(".pdf")),
        "concordance_svg": str(output_stem.with_suffix(".svg")),
    }
    for path in outputs.values():
        fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return outputs


def build_module7_report_payload(
    concordance: pd.DataFrame,
    enrichment_summary: pd.DataFrame,
    risk_flags: pd.DataFrame,
    inputs: dict,
    outputs: dict,
    tf_matrix: pd.DataFrame | None = None,
    pathway_matrix: pd.DataFrame | None = None,
    figure_outputs: dict | None = None,
) -> dict:
    tf_matrix = tf_matrix if tf_matrix is not None else pd.DataFrame()
    pathway_matrix = pathway_matrix if pathway_matrix is not None else pd.DataFrame()
    figure_outputs = figure_outputs or {}
    return {
        "module": "7.5",
        "method": "Integrated scTenifoldKnk replication and enrichment summary for CellOracle candidates",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": inputs,
        "outputs": outputs,
        "n_tfs_with_sctenifoldknk_results": int(concordance["tf"].nunique()) if "tf" in concordance.columns else 0,
        "n_enrichment_rows": int(len(enrichment_summary)),
        "n_tf_level_matrix_rows": int(len(tf_matrix)),
        "n_pathway_level_matrix_rows": int(len(pathway_matrix)),
        "n_figure_outputs": int(len(figure_outputs)),
        "n_review_risk_flags": int(len(risk_flags)),
        "python_runtime": {"version": platform.python_version(), "platform": platform.platform()},
    }


def build_main_conclusions(tf_matrix: pd.DataFrame, enrichment_summary: pd.DataFrame, risk_flags: pd.DataFrame) -> str:
    lines = ["# Module 7.5 scTenifoldKnk Integration Conclusions", ""]
    if not tf_matrix.empty and {"tf", "candidate_tier", "module7_integrated_rank"}.issubset(tf_matrix.columns):
        top = tf_matrix.sort_values("module7_integrated_rank").head(5)
        lines.append("## Top Integrated TFs")
        lines.append(", ".join(top["tf"].astype(str).tolist()))
        lines.append("")
        tier1 = tf_matrix.loc[tf_matrix["display_group"].astype(str).eq("main_tier1"), "tf"].astype(str).tolist()
        if tier1:
            lines.append("## Main Tier 1 Replication Set")
            lines.append(", ".join(tier1))
            lines.append("")
        sox4 = tf_matrix.loc[tf_matrix["tf"].astype(str).eq("SOX4")]
        if not sox4.empty and "malignant_like_state_specific_gene_count" in sox4.columns:
            lines.append("## SOX4 State-Specific Evidence")
            lines.append(
                f"SOX4 malignant-like state-specific perturbation genes: "
                f"{int(sox4.iloc[0]['malignant_like_state_specific_gene_count'])}"
            )
            lines.append("")
        ap1 = tf_matrix.loc[tf_matrix["display_group"].astype(str).eq("main_ap1_axis"), "tf"].astype(str).tolist()
        if ap1:
            lines.append("## AP-1 Axis")
            lines.append(", ".join(ap1))
            lines.append("")
    if not enrichment_summary.empty:
        lines.append("## Enrichment Scope")
        databases = sorted(enrichment_summary["database"].astype(str).unique()) if "database" in enrichment_summary.columns else []
        lines.append(f"Pathway databases summarized: {', '.join(databases)}")
        lines.append("")
    if len(risk_flags):
        lines.append("## Review Risk Flags")
        lines.append(f"{len(risk_flags)} Module 7 review-risk flags were generated.")
    else:
        lines.append("## Review Risk Flags")
        lines.append("No Module 7 review-risk flags were generated.")
    lines.append("")
    return "\n".join(lines)


def run_module7_summary(metadata_dir: Path, figure_dir: Path) -> dict:
    concordance_path = metadata_dir / "sctenifoldknk_module7_3_concordance_summary.tsv"
    integrated_path = metadata_dir / "sctenifoldknk_module7_3_integrated_evidence_matrix.tsv"
    axis_path = metadata_dir / "sctenifoldknk_module7_3_biological_axis_summary.tsv"
    state_specific_path = metadata_dir / "sctenifoldknk_module7_3_malignant_like_state_specific_genes.tsv"
    enrichment_path = metadata_dir / "sctenifoldknk_module7_4_top_enrichment_summary.tsv"
    risks_path = metadata_dir / "sctenifoldknk_module7_3_review_risk_flags.tsv"

    concordance = read_tsv_or_empty(concordance_path)
    integrated = read_tsv_or_empty(integrated_path)
    biological_axis = read_tsv_or_empty(axis_path)
    state_specific = read_tsv_or_empty(state_specific_path)
    enrichment = read_tsv_or_empty(enrichment_path)
    risks = read_tsv_or_empty(risks_path)

    tf_matrix = build_tf_level_replication_matrix(integrated, biological_axis, state_specific)
    pathway_matrix = build_pathway_level_enrichment_matrix(enrichment, tf_matrix)
    figure_outputs = plot_celloracle_sctenifold_concordance(tf_matrix, figure_dir)

    outputs = {
        "report": str(metadata_dir / "sctenifoldknk_module7_5_report.json"),
        "tf_level_replication_matrix": str(metadata_dir / "sctenifoldknk_module7_5_tf_level_replication_matrix.tsv"),
        "pathway_level_enrichment_matrix": str(metadata_dir / "sctenifoldknk_module7_5_pathway_level_enrichment_matrix.tsv"),
        "main_conclusions": str(metadata_dir / "sctenifoldknk_module7_5_main_conclusions.md"),
        "supplementary_table_index": str(metadata_dir / "sctenifoldknk_module7_5_supplementary_table_index.tsv"),
        **figure_outputs,
    }
    inputs = {
        "concordance": str(concordance_path),
        "integrated_evidence": str(integrated_path),
        "biological_axis": str(axis_path),
        "state_specific_genes": str(state_specific_path),
        "enrichment_summary": str(enrichment_path),
        "risk_flags": str(risks_path),
    }
    metadata_dir.mkdir(parents=True, exist_ok=True)
    tf_matrix.to_csv(outputs["tf_level_replication_matrix"], sep="\t", index=False)
    pathway_matrix.to_csv(outputs["pathway_level_enrichment_matrix"], sep="\t", index=False)
    report = build_module7_report_payload(
        concordance,
        enrichment,
        risks,
        inputs=inputs,
        outputs=outputs,
        tf_matrix=tf_matrix,
        pathway_matrix=pathway_matrix,
        figure_outputs=figure_outputs,
    )
    Path(outputs["report"]).write_text(json.dumps(report, indent=2), encoding="utf-8")
    Path(outputs["main_conclusions"]).write_text(build_main_conclusions(tf_matrix, enrichment, risks), encoding="utf-8")
    supp = pd.DataFrame(
        [
            {
                "table_id": "Supplementary Table 7.1",
                "primary_source_file": "sctenifoldknk_module7_5_tf_level_replication_matrix.tsv",
                "description": "TF-level replication matrix integrating CellOracle tiering, scTenifoldKnk perturbation, GRN overlap, marker evidence and state specificity.",
                "intended_manuscript_use": "Main Figure 7A and supplementary control calibration.",
            },
            {
                "table_id": "Supplementary Table 7.2",
                "primary_source_file": "sctenifoldknk_module7_5_pathway_level_enrichment_matrix.tsv",
                "description": "Pathway-level enrichment matrix labeling KEGG, Reactome, GO BP ORA and preranked GSEA terms for main and supplementary TF groups.",
                "intended_manuscript_use": "Main Figure 7B pathway panel and supplementary pathway calibration.",
            },
            {
                "table_id": "Figure 7.1",
                "primary_source_file": Path(outputs.get("concordance_png", "")).name,
                "description": "CellOracle quantitative perturbation score versus scTenifoldKnk replication score concordance plot.",
                "intended_manuscript_use": "CellOracle vs scTenifoldKnk concordance figure.",
            },
        ]
    )
    supp.to_csv(outputs["supplementary_table_index"], sep="\t", index=False)
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 7.5 summarize scTenifoldKnk replication and enrichment")
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    report = run_module7_summary(args.metadata_dir, args.figure_dir)
    print(
        json.dumps(
            {
                "report": report["outputs"]["report"],
                "n_tfs": report["n_tfs_with_sctenifoldknk_results"],
                "tf_level_matrix": report["outputs"]["tf_level_replication_matrix"],
                "pathway_level_matrix": report["outputs"]["pathway_level_enrichment_matrix"],
                "n_figures": report["n_figure_outputs"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
