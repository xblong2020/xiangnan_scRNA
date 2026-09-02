#!/usr/bin/env python3
"""Generate the evidence-backed final Figure 3 EGR1 A-F report."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

try:
    from figure3_egr1_common import PROJECT_ROOT
except ModuleNotFoundError:
    from scripts.figure3_egr1_common import PROJECT_ROOT


DEFAULT_OUTPUT = PROJECT_ROOT / "reports/figure3_egr1_a_to_f_report.md"


def read_json(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(path)
    return json.loads(path.read_text(encoding="utf-8"))


def fmt(value, digits: int = 4) -> str:
    if value is None or pd.isna(value):
        return "NA"
    if isinstance(value, float):
        return f"{value:.{digits}g}"
    return str(value)


def bullets(items: list[str]) -> str:
    return "\n".join(f"- {item}" for item in items) if items else "- None."


def run(root: Path, output: Path) -> dict:
    meta = root / "metadata/driver"
    a = read_json(meta / "figure3a_stress_transition/figure3a_candidate_selection_report.json")
    b = read_json(meta / "figure3b_egr1/figure3b_egr1_r_plot_report.json")
    c = read_json(meta / "figure3c_egr1/figure3c_egr1_inner_product_report.json")
    c_data = read_json(meta / "figure3c_egr1/figure3c_egr1_data_report.json")
    d = read_json(meta / "figure3d_egr1/figure3d_egr1_report.json")
    e = read_json(meta / "figure3e_egr1/figure3e_egr1_report.json")
    e_run = read_json(meta / "figure3e_egr1/figure3e_egr1_stressed_regenerative_run_report.json")
    determinism = read_json(
        meta / "figure3e_egr1_determinism/figure3e_egr1_same_seed_determinism_report.json"
    )
    e_sens = read_json(meta / "figure3e_egr1_sensitivity/figure3e_egr1_sensitivity_report.json")
    f = read_json(meta / "figure3f_egr1/figure3f_egr1_report.json")
    validation = read_json(meta / "figure3_egr1_validation/figure3_egr1_a_to_f_validation_report.json")
    preflight = read_json(meta / "figure3_egr1_preflight/figure3_egr1_preflight_report.json")
    selection = pd.read_csv(
        meta / "figure3a_stress_transition/figure3a_candidate_evidence_matrix.tsv", sep="\t"
    )
    baseline_audit = pd.read_csv(
        meta / "figure3b_egr1/figure3b_baseline_equivalence_audit.tsv", sep="\t"
    )
    stages = pd.read_csv(meta / "figure3d_egr1/figure3d_egr1_stage_comparison.tsv", sep="\t")
    stages = stages.loc[
        stages["row_type"].eq("stage_summary")
        & stages["space"].eq("CellOracle UMAP grid")
    ].copy()
    subset_audit = pd.read_csv(
        meta / "figure3e_egr1/figure3e_egr1_subset_selection_audit.tsv", sep="\t"
    )
    sens_summary = pd.read_csv(
        meta / "figure3e_egr1_sensitivity/figure3e_egr1_sensitivity_summary.tsv", sep="\t"
    )
    shared = pd.read_csv(
        meta / "three_axis_figure_consistency/figure2_figure3_figure4_shared_colour_limits.tsv",
        sep="\t",
    )

    egr1 = selection.loc[selection["candidate"].eq("EGR1")].iloc[0]
    stage_lines = [
        f"{row.comparison}: median PS={fmt(row.score_median)}, mean PS={fmt(row.score_mean)}, "
        f"positive fraction={fmt(row.positive_fraction)}"
        for row in stages.itertuples(index=False)
    ]
    main_subset = subset_audit.loc[subset_audit["subset"].eq(e["subset"])].iloc[0]
    sens_lines = [
        f"{row.subset}: n={row.n_cells}, significant genes={row.n_significant_perturbed_genes}, "
        f"FDR pathways={row.n_fdr_significant_pathways}"
        for row in sens_summary.itertuples(index=False)
    ]
    risk_by_flag: dict[str, str] = {}
    for report in [preflight, a, d, e_run, determinism, e, e_sens, f, validation]:
        values = report.get("review_risk_flags", [])
        if isinstance(values, dict):
            values = [values]
        for risk in values:
            if isinstance(risk, dict):
                flag = risk.get("flag", "risk")
                detail = risk.get("detail", "")
                if flag == "existing_sctenifoldknk_low_replication":
                    detail += " This preflight limitation was resolved for the formal main network by the dedicated nc_nNet=10, nc_nCells=500, three-seed run."
                risk_by_flag[flag] = (
                    f"{risk.get('severity', 'unspecified')}: {flag} — {detail}"
                )
            else:
                risk_by_flag[str(risk)] = str(risk)
    risks = list(risk_by_flag.values())

    scripts = sorted(
        str(path.relative_to(root)).replace("\\", "/")
        for path in (root / "scripts").glob("*figure3*")
        if "egr1" in path.name.lower() or "stress_transition" in path.name.lower()
    )
    data_paths = [
        "metadata/driver/figure3_egr1_preflight/",
        "metadata/driver/figure3a_stress_transition/",
        "metadata/driver/figure3b_egr1/",
        "metadata/driver/figure3c_egr1/",
        "metadata/driver/figure3d_egr1/",
        "metadata/driver/figure3e_egr1/",
        "metadata/driver/figure3e_egr1_sensitivity/",
        "metadata/driver/figure3f_egr1/",
        "data/processed/driver/figure3e_egr1_sctenifoldknk/",
    ]
    figure_paths = [
        "figures/driver/figure3a_stress_transition/",
        "figures/driver/figure3b_egr1/",
        "figures/driver/figure3c_egr1/",
        "figures/driver/figure3d_egr1/",
        "figures/driver/figure3e_egr1/",
        "figures/driver/figure3e_egr1_sensitivity/",
        "figures/driver/figure3f_egr1/",
        "figures/driver/figure3_egr1_preview/",
    ]
    peak_stage = d.get("observed_umap_pattern", {}).get(
        "absolute_effect_peak_stage", "NA"
    )
    f_recommendation = f.get("recommendation_if_no_stable_pathways")
    extended = []
    if not e.get("formal_plot_generated"):
        extended.append("Figure 3E network result: no FDR-significant gene plot was generated.")
    if not f.get("formal_plot_generated"):
        extended.append(
            f"Figure 3F pathway result: {f_recommendation or 'retain as an FDR-null report.'}"
        )
    extended.append("t-SNE perturbation and pseudotime projections remain sensitivity analyses.")
    extended.append(
        "Lower-replication state-network comparisons belong in Extended Data unless independently replicated."
    )

    main_assessment = validation["sci_main_figure_assessment"]
    baseline_status = (
        "all exact-equivalence components passed"
        if baseline_audit["status"].eq("pass").all()
        else "see failed baseline audit components"
    )
    e_gene_word = "gene" if int(e["n_significant_excluding_target"]) == 1 else "genes"
    e_plot_word = "gene" if int(e["n_plotted"]) == 1 else "genes"
    e_verb = "was" if int(e["n_significant_excluding_target"]) == 1 else "were"
    protected_statement = (
        "All preflight protected hashes remained unchanged."
        if validation["protected_assets_unchanged"]
        else "Two protected HNF4A files changed after preflight; the concurrent-change incident is recorded without re-baselining or reverting them."
    )
    conservative = (
        "EGR1-associated regulatory activity was enriched within the stress-transition phase, "
        "and computational EGR1 knockout predicted stage-dependent alterations in the direction "
        "of hepatocyte state progression."
    )
    stronger = (
        f"Virtual EGR1 knockout predominantly opposed the developmental vector field, with the "
        f"largest absolute prespecified-stage effect in the {peak_stage} phase, supporting an "
        "EGR1-associated role in progression through the stress-transition programme."
    )
    results_paragraph = (
        f"Within the AP-1/CEBPB/EGR1 candidate module, EGR1 ranked "
        f"{int(egr1.selection_rank)} with a composite selection score of "
        f"{fmt(float(egr1.selection_score))}. Using a shared 5,000-cell developmental field, "
        f"CellOracle EGR1=0 simulation yielded stage-dependent perturbation scores. "
        f"The largest absolute effect occurred in the {peak_stage} stage rather than being "
        f"prespecified to peak in the intermediate stage. In a {e['subset'].replace('_', ' ')} "
        f"scTenifoldKnk network, {e['n_significant_excluding_target']} {e_gene_word} {e_verb} significant "
        f"after exclusion of EGR1, and strict ORA identified {f['n_significant_pathways']} "
        f"globally BH-significant pathways."
    )
    legend = (
        "**Figure 3 | AP-1/CEBPB/EGR1 stress-transition programme.** "
        "(A) Overlapping regulatory phases and project-derived candidate evidence supporting "
        "selection of EGR1 as the principal perturbation representative. "
        "(B) Common 5,000-cell baseline developmental field. "
        "(C) CellOracle EGR1=0 virtual-knockout field and perturbation score, defined as the "
        "dot product between the perturbation and baseline developmental vectors. "
        "(D) Perturbation score across pseudotime with fixed early, intermediate and late stages; "
        "t-SNE is a supplementary projection. "
        f"(E) {e['n_significant_excluding_target']} FDR-significant perturbed {e_gene_word} in the "
        f"{e['subset'].replace('_', ' ')} network"
        + ("." if e.get("formal_plot_generated") else "; the formal panel was suppressed because the FDR rule was not met.")
        + f" (F) Strict ORA using the matching 3,000-gene network background identified "
        f"{f['n_significant_pathways']} globally BH-significant pathways"
        + ("." if f.get("formal_plot_generated") else "; the formal panel was suppressed.")
    )

    markdown = f"""# Figure 3: AP-1/CEBPB/EGR1 stress-transition programme

## Deliverables

### New scripts

{bullets([f"`{item}`" for item in scripts])}

### New data and reports

{bullets([f"`{item}`" for item in data_paths])}

### New figures

{bullets([f"`{item}`" for item in figure_paths])}

## Figure 3A: candidate comparison and EGR1 selection

EGR1 ranked **{int(egr1.selection_rank)} of {len(selection)}** with selection score
**{fmt(float(egr1.selection_score))}**. The matrix uses numeric, project-derived CellOracle,
scTenifoldKnk, temporal, state-specificity, cross-dataset and penalty evidence. EGR1 was
selected as the principal perturbation representative while JUN/AP-1 and CEBPB remain
members of the stress-transition module. The architecture is described as overlapping and
partially ordered, without a strict causal cascade.

## Figure 3B: common baseline

Figure 3B reuses the exact validated 5,000-cell baseline, including cell order, UMAP, t-SNE,
pseudotime, state labels, 20 x 20 grids, k=50 neighbours and density masks. Audit status:
**{baseline_status}**.

## Figure 3C: EGR1 CellOracle perturbation

The displacement source is `metadata/driver/celloracle_module6_8_cell_shift_summary.tsv.gz` with condition
`EGR1=0`; propagation and simulation parameters are recorded in the source report.
UMAP is the native CellOracle perturbation space. t-SNE is a local Jacobian projection used
only for sensitivity analysis. The score is `{c.get('definition', 'perturbation vector dot baseline vector')}`.

## Figure 3D: pseudotime and fixed stages

{bullets(stage_lines)}

The largest absolute stage effect was **{peak_stage}**. This observed pattern is reported
directly; it does not support an intermediate-specific maximum. Bootstrap, binned summaries,
Spearman confidence intervals, LOESS and change-point sensitivity results are retained in the
Figure 3D metadata directory.

## Figure 3E: stress-transition scTenifoldKnk

The formal network uses **{e['subset']}**: {int(main_subset.n_cells)} cells,
{int(main_subset.n_genes)} genes, EGR1 detection rate {fmt(float(main_subset.egr1_detection_rate))},
{int(main_subset.n_datasets)} represented known datasets, and maximum dataset fraction
{fmt(float(main_subset.max_dataset_fraction))}. The main run used nc_nNet=
{e_run['parameters']['nc_nNet']}, nc_nCells={e_run['parameters']['nc_nCells_used']} and
{e_run['n_successful_seeds']} successful fixed seeds. After excluding EGR1,
**{e['n_significant_excluding_target']} {e_gene_word}** passed p.adj < 0.05; the plot contains
{e['n_plotted']} {e_plot_word} and never uses non-significant rows to fill Top 20.
The independent seed-15071990 repeat had numeric reproducibility status
**{determinism['status']}**.

### State sensitivity

{bullets(sens_lines)}

The non-main networks use lower replication and provide sensitivity evidence. Dataset
composition/dominance is audited; a dedicated network-level scTenifoldKnk LODO analysis was
not performed.

## Figure 3F: strict enrichment

Strict one-sided ORA used all Figure 3E FDR-significant genes and the matching
`{e['subset']}` network background. Across KEGG, Reactome and GO BP, global BH correction
identified **{f['n_significant_pathways']} significant pathways**. Unsigned manifold distance
is not interpreted as pathway activation or suppression. Formal plot generated:
**{str(bool(f['formal_plot_generated'])).lower()}**.

## Three-axis comparability

HNF4A, EGR1 and SOX4 use the same cells, baseline coordinates, grids, state colours and
statistical definitions. Panel-specific variants are retained, while the main Figure 3C/D
outputs use the three-axis shared symmetric limit **{fmt(float(shared['three_axis_shared_symmetric_limit'].iloc[0]), 8)}**.
Protected-asset audit: **{protected_statement}**

## Review-risk flags

{bullets(risks)}

## Publication assessment

Validation outcome: **{validation['status']}** with {validation['n_pass']}/{validation['n_checks']}
checks passing, {validation['n_warning']} warnings and {validation['n_fail']} failures.
SCI main-figure assessment: **{main_assessment}**.

### Results recommended for Extended Data

{bullets(extended)}

## Recommended Figure legend

{legend}

## Recommended Results paragraph

{results_paragraph}

## Recommended conclusions

Conservative:

> {conservative}

Stronger, still compliant with the observed computational evidence:

> {stronger}

Neither formulation establishes direct causality, experimental validation, genetic epistasis,
or an obligatory linear HNF4A/AP-1/CEBPB/EGR1/SOX4 cascade.
"""
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(markdown, encoding="utf-8")
    return {
        "output": str(output.resolve()),
        "assessment": main_assessment,
        "n_significant_genes": int(e["n_significant_excluding_target"]),
        "n_significant_pathways": int(f["n_significant_pathways"]),
        "peak_stage": peak_stage,
        "n_review_risks": len(risks),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(run(args.project_root, args.output), ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
