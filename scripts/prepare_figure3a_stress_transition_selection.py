#!/usr/bin/env python3
"""Build the Figure 3A AP-1/CEBPB/EGR1 evidence matrix from project results."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

try:
    from figure3_egr1_common import (
        PENALTY_WEIGHTS,
        PROJECT_ROOT,
        SELECTION_WEIGHTS,
        TARGET_TF,
        compute_selection_score,
        json_safe,
        minmax_scale,
        write_json,
    )
except ModuleNotFoundError:
    from scripts.figure3_egr1_common import (
        PENALTY_WEIGHTS,
        PROJECT_ROOT,
        SELECTION_WEIGHTS,
        TARGET_TF,
        compute_selection_score,
        json_safe,
        minmax_scale,
        write_json,
    )


DEFAULT_H5AD = PROJECT_ROOT / "data/processed/driver/celloracle_module6_6/celloracle_module6_6_input.h5ad"
DEFAULT_OUT_DIR = PROJECT_ROOT / "metadata/driver/figure3a_stress_transition"
AP1_MEMBERS = ["JUN", "JUNB", "JUND", "FOS", "ATF3"]
AP1_DIAGRAM_MEMBERS = ["JUN", "JUNB", "JUND", "FOS", "FOSB", "ATF3"]
INDIVIDUAL_CANDIDATES = AP1_MEMBERS + ["CEBPB", "EGR1"]


def expression_detection_by_state(adata: ad.AnnData, genes: list[str]) -> pd.DataFrame:
    matrix = adata.layers["counts"] if "counts" in adata.layers else adata.X
    rows = []
    states = adata.obs["celloracle_state"].astype(str).to_numpy()
    for gene in genes:
        if gene not in adata.var_names:
            continue
        ix = int(adata.var_names.get_loc(gene))
        column = matrix[:, ix]
        values = column.toarray().ravel() if sparse.issparse(column) else np.asarray(column).ravel()
        for state in sorted(pd.unique(states)):
            sub = values[states == state]
            rows.append(
                {
                    "tf": gene,
                    "celloracle_state": state,
                    "n_cells": int(len(sub)),
                    "mean_expression": float(np.mean(sub)),
                    "detection_rate": float(np.mean(sub > 0)),
                }
            )
    return pd.DataFrame(rows)


def sctenifold_subset_metrics(root: Path) -> pd.DataFrame:
    frames = []
    for subset in ["driver_union_all", "main_strict", "malignant_like"]:
        path = root / f"metadata/driver/sctenifoldknk_module7_2_{subset}_perturbation_genes.tsv"
        data = pd.read_csv(path, sep="\t")
        data["p.adj"] = pd.to_numeric(data["p.adj"], errors="coerce")
        data["distance"] = pd.to_numeric(data["distance"], errors="coerce")
        data = data.loc[
            data["tf"].astype(str).isin(INDIVIDUAL_CANDIDATES)
            & ~data["gene"].astype(str).eq(data["tf"].astype(str))
        ]
        summary = (
            data.groupby("tf", as_index=False)
            .agg(
                n_tested=("gene", "size"),
                n_significant=("p.adj", lambda values: int(np.sum(pd.to_numeric(values, errors="coerce") < 0.05))),
                mean_distance=("distance", "mean"),
                max_distance=("distance", "max"),
            )
            .assign(subset=subset)
        )
        summary["subset_rank"] = summary["n_significant"].rank(method="min", ascending=False)
        frames.append(summary)
    return pd.concat(frames, ignore_index=True)


def build_member_metrics(root: Path, h5ad_path: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    quantitative = pd.read_csv(root / "metadata/driver/celloracle_module6_9b_quantitative_tf_scores.tsv", sep="\t")
    stability = pd.read_csv(root / "metadata/driver/celloracle_module6_10_rank_stability.tsv", sep="\t")
    phase = pd.read_csv(root / "metadata/driver/celloracle_module6_10_phase_wide_summary.tsv", sep="\t")
    selection = pd.read_csv(root / "metadata/driver/celloracle_tf_selection.module6_4.tsv", sep="\t")
    integrated = pd.read_csv(root / "metadata/driver/sctenifoldknk_module7_3_integrated_evidence_matrix.tsv", sep="\t")
    state_shift = pd.read_csv(root / "metadata/driver/celloracle_module6_8_state_shift_summary.tsv", sep="\t")
    sct = sctenifold_subset_metrics(root)

    adata = ad.read_h5ad(h5ad_path)
    expression = expression_detection_by_state(adata, INDIVIDUAL_CANDIDATES)
    adata.file.close() if getattr(adata, "file", None) and adata.file.is_open else None

    transition_states = {"stressed_injured", "regenerative_progenitor"}
    shift_rows = []
    for tf, frame in state_shift.loc[state_shift["tf"].isin(INDIVIDUAL_CANDIDATES)].groupby("tf"):
        values = pd.to_numeric(frame["embedding_shift_norm_mean"], errors="coerce")
        weights = pd.to_numeric(frame["n_cells"], errors="coerce")
        transition = frame["celloracle_state"].astype(str).isin(transition_states)
        transition_mean = float(np.average(values[transition], weights=weights[transition]))
        other_mean = float(np.average(values[~transition], weights=weights[~transition]))
        shift_rows.append(
            {
                "tf": tf,
                "transition_embedding_shift_mean": transition_mean,
                "non_transition_embedding_shift_mean": other_mean,
                "transition_shift_ratio": transition_mean / max(other_mean, 1e-12),
            }
        )
    shift_specificity = pd.DataFrame(shift_rows)

    expression_rows = []
    for tf, frame in expression.groupby("tf"):
        transition = frame["celloracle_state"].isin(transition_states)
        transition_detection = float(
            np.average(frame.loc[transition, "detection_rate"], weights=frame.loc[transition, "n_cells"])
        )
        non_transition_detection = float(
            np.average(frame.loc[~transition, "detection_rate"], weights=frame.loc[~transition, "n_cells"])
        )
        breadth = float(np.mean(frame["detection_rate"] >= 0.01))
        transition_excess = max(0.0, transition_detection - non_transition_detection)
        expression_rows.append(
            {
                "tf": tf,
                "transition_detection_rate": transition_detection,
                "non_transition_detection_rate": non_transition_detection,
                "detection_breadth_fraction": breadth,
                "generic_stress_raw": breadth * (1.0 - min(transition_excess, 1.0)),
            }
        )
    expression_summary = pd.DataFrame(expression_rows)

    sct_wide = sct.pivot(index="tf", columns="subset", values=["n_significant", "subset_rank"])
    sct_wide.columns = [f"{metric}_{subset}" for metric, subset in sct_wide.columns]
    sct_wide = sct_wide.reset_index()
    count_cols = [column for column in sct_wide.columns if column.startswith("n_significant_")]
    rank_cols = [column for column in sct_wide.columns if column.startswith("subset_rank_")]
    sct_wide["sct_significant_mean"] = sct_wide[count_cols].mean(axis=1)
    sct_wide["sct_significant_cv"] = sct_wide[count_cols].std(axis=1) / sct_wide["sct_significant_mean"].clip(lower=1)
    sct_wide["sct_rank_mean"] = sct_wide[rank_cols].mean(axis=1)

    keep_quantitative = [
        "tf",
        "quantitative_rank",
        "quantitative_perturbation_score",
        "proliferation_module_rescue_score",
        "state_specificity_ratio",
    ]
    keep_stability = [
        "tf",
        "max_rank",
        "top5_fraction",
        "mean_quantitative_score",
        "min_quantitative_score",
    ]
    keep_phase = [
        "tf",
        "phase_early_score",
        "phase_intermediate_score",
        "phase_late_score",
        "phase_early_rank",
        "phase_intermediate_rank",
        "phase_late_rank",
    ]
    keep_selection = [
        "tf",
        "biology_score",
        "dataset_direction_consistency_fraction",
        "loo_min_directional_r",
        "detected_dataset_count",
    ]
    keep_integrated = [
        "tf",
        "integrated_rank",
        "module7_integrated_rank",
        "integrated_module7_score",
        "scTenifoldKnk_rank",
        "n_significant_perturbed_genes",
        "top_gene_grn_target_jaccard",
    ]
    member = (
        pd.DataFrame({"tf": INDIVIDUAL_CANDIDATES})
        .merge(quantitative[keep_quantitative], on="tf", how="left")
        .merge(stability[keep_stability], on="tf", how="left")
        .merge(phase[keep_phase], on="tf", how="left")
        .merge(selection[keep_selection], on="tf", how="left")
        .merge(integrated[keep_integrated], on="tf", how="left")
        .merge(shift_specificity, on="tf", how="left")
        .merge(expression_summary, on="tf", how="left")
        .merge(sct_wide, on="tf", how="left")
    )

    phase_centers = np.array([1 / 6, 0.5, 5 / 6])
    phase_values = member[["phase_early_score", "phase_intermediate_score", "phase_late_score"]].to_numpy(float)
    phase_sums = phase_values.sum(axis=1)
    member["temporal_center"] = np.divide(
        phase_values @ phase_centers,
        phase_sums,
        out=np.full(len(member), 0.5),
        where=phase_sums > 0,
    )
    member["temporal_proximity_raw"] = np.clip(1.0 - 2.0 * np.abs(member["temporal_center"] - 0.5), 0, 1)
    member["intermediate_dominance_raw"] = np.divide(
        member["phase_intermediate_score"],
        member[["phase_early_score", "phase_intermediate_score", "phase_late_score"]].max(axis=1).clip(lower=1e-12),
    )

    member["celloracle_robustness"] = (
        minmax_scale(member["quantitative_perturbation_score"])
        + minmax_scale(member["mean_quantitative_score"])
        + minmax_scale(member["top5_fraction"])
    ) / 3
    member["sctenifoldknk_robustness"] = (
        minmax_scale(member["integrated_module7_score"])
        + minmax_scale(member["sct_significant_mean"])
        + (1.0 - minmax_scale(member["sct_significant_cv"]))
    ) / 3
    member["transition_state_specificity"] = (
        minmax_scale(member["transition_shift_ratio"])
        + minmax_scale(member["transition_detection_rate"] - member["non_transition_detection_rate"])
    ) / 2
    member["temporal_positioning"] = (
        member["temporal_proximity_raw"] + member["intermediate_dominance_raw"]
    ) / 2
    member["cross_dataset_stability"] = (
        minmax_scale(member["top5_fraction"])
        + (1.0 - minmax_scale(member["max_rank"]))
        + minmax_scale(member["dataset_direction_consistency_fraction"])
        + minmax_scale(member["loo_min_directional_r"])
    ) / 4
    member["pathway_interpretability"] = (
        minmax_scale(member["n_significant_perturbed_genes"])
        + minmax_scale(member["top_gene_grn_target_jaccard"])
    ) / 2
    member["cross_method_concordance"] = 1.0 - np.abs(
        member["celloracle_robustness"] - member["sctenifoldknk_robustness"]
    )
    member["generic_stress_penalty"] = minmax_scale(member["generic_stress_raw"])
    member["proliferation_dependency_penalty"] = minmax_scale(
        member["proliferation_module_rescue_score"].abs()
    )
    member["literature_overlap_penalty"] = minmax_scale(member["biology_score"])
    member["literature_overlap"] = member["literature_overlap_penalty"]

    provenance = {
        "celloracle_robustness": "Module 6.9b quantitative score, Module 6.10 mean score, and top-5 subset stability.",
        "sctenifoldknk_robustness": "Module 7 integrated score, mean FDR-significant gene count across three existing subsets, and inverse count CV.",
        "transition_state_specificity": "CellOracle shift ratio and expression-detection excess in stressed_injured + regenerative_progenitor versus other states.",
        "temporal_positioning": "Proximity of phase-score center of mass to pseudotime 0.5 and intermediate-phase dominance.",
        "cross_dataset_stability": "LODO/LOSO top-5 fraction, worst rank, dataset direction consistency, and minimum directional correlation.",
        "pathway_interpretability": "FDR-significant perturbation-gene count and CellOracle GRN overlap Jaccard.",
        "generic_stress_penalty": "Data-derived expression breadth with reduced penalty for transition-state detection excess.",
        "proliferation_dependency_penalty": "Absolute CellOracle proliferation-module perturbation score, used as a dependency proxy.",
        "literature_overlap_penalty": "Existing Module 6.4 biology_score used as a project-recorded literature/biology-prior proxy; no star rating was assigned.",
    }
    return member, expression, provenance


def collapse_candidates(member: pd.DataFrame) -> pd.DataFrame:
    metric_columns = [
        "celloracle_robustness",
        "sctenifoldknk_robustness",
        "transition_state_specificity",
        "temporal_positioning",
        "cross_dataset_stability",
        "pathway_interpretability",
        "cross_method_concordance",
        "generic_stress_penalty",
        "proliferation_dependency_penalty",
        "literature_overlap_penalty",
        "literature_overlap",
    ]
    group_map = {tf: "JUN/AP-1" for tf in AP1_MEMBERS}
    group_map.update({"CEBPB": "CEBPB", "EGR1": "EGR1"})
    grouped = member.assign(candidate=member["tf"].map(group_map)).groupby("candidate", as_index=False)
    candidate = grouped[metric_columns].mean()
    candidate["member_tfs"] = candidate["candidate"].map(
        {
            "JUN/AP-1": ";".join(AP1_MEMBERS),
            "CEBPB": "CEBPB",
            "EGR1": "EGR1",
        }
    )
    candidate["missing_diagram_members"] = candidate["candidate"].map(
        {"JUN/AP-1": "FOSB (not retained in the 3,000-gene CellOracle/scTenifoldKnk network)", "CEBPB": "", "EGR1": ""}
    )
    # Re-standardize every scored field across the three displayed candidates.
    for column in metric_columns:
        candidate[column] = minmax_scale(candidate[column])
    candidate["generic_stress_risk"] = candidate["generic_stress_penalty"]
    candidate["proliferation_dependency"] = candidate["proliferation_dependency_penalty"]
    candidate["celloracle_evidence"] = candidate["celloracle_robustness"]
    candidate["sctenifoldknk_evidence"] = candidate["sctenifoldknk_robustness"]
    candidate["leave_one_dataset_out_stability"] = candidate["cross_dataset_stability"]
    scored = compute_selection_score(candidate)
    scored["selection_score_scaled"] = minmax_scale(scored["selection_score"])
    scored["final_role"] = scored["candidate"].map(
        {
            "JUN/AP-1": "Stress-transition programme members",
            "CEBPB": "Stress-transition programme member",
            "EGR1": "Principal perturbation representative",
        }
    )
    return scored


def run(h5ad_path: Path, out_dir: Path) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    member, expression, provenance = build_member_metrics(PROJECT_ROOT, h5ad_path)
    candidate = collapse_candidates(member)
    candidate_path = out_dir / "figure3a_candidate_evidence_matrix.tsv"
    member_path = out_dir / "figure3a_candidate_member_metrics.tsv"
    expression_path = out_dir / "figure3a_candidate_expression_by_state.tsv"
    candidate.to_csv(candidate_path, sep="\t", index=False)
    member.to_csv(member_path, sep="\t", index=False)
    expression.to_csv(expression_path, sep="\t", index=False)

    egr1 = candidate.loc[candidate["candidate"].eq(TARGET_TF)].iloc[0]
    leader = candidate.iloc[0]
    egr1_is_first = int(egr1["selection_rank"]) == 1
    score_gap = float(leader["selection_score"] - egr1["selection_score"])
    review_risks = []
    if not egr1_is_first:
        review_risks.append(
            {
                "flag": "egr1_not_highest_selection_score",
                "severity": "review_attention",
                "detail": (
                    f"EGR1 ranked {int(egr1['selection_rank'])}; {leader['candidate']} ranked first "
                    f"with a selection-score gap of {score_gap:.4f}."
                ),
            }
        )
    review_risks.extend(
        [
            {
                "flag": "egr1_scenic_auc_missing",
                "severity": "review_attention",
                "detail": "No EGR1-specific SCENIC regulon AUC was available in the prepared h5ad.",
            },
            {
                "flag": "literature_overlap_proxy",
                "severity": "methodological_caveat",
                "detail": "Literature overlap uses the existing Module 6.4 biology_score as a recorded project prior rather than a new systematic bibliography.",
            },
        ]
    )
    report = {
        "module": "Figure 3A",
        "target_tf": TARGET_TF,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "programme": "AP-1/CEBPB/EGR1 stress-transition programme",
        "ap1_diagram_members": AP1_DIAGRAM_MEMBERS,
        "ap1_scored_members": AP1_MEMBERS,
        "selection_weights": SELECTION_WEIGHTS,
        "penalty_weights": PENALTY_WEIGHTS,
        "metric_provenance": provenance,
        "candidate_ranking": candidate[
            ["candidate", "selection_rank", "selection_score", "selection_score_scaled", "final_role"]
        ].to_dict(orient="records"),
        "egr1_selection": {
            "actual_rank": int(egr1["selection_rank"]),
            "selection_score": float(egr1["selection_score"]),
            "highest_candidate": str(leader["candidate"]),
            "score_gap_to_highest": score_gap,
            "empirically_highest": bool(egr1_is_first),
            "non_statistical_reason_if_not_first": (
                None
                if egr1_is_first
                else "EGR1 was prespecified as the individual perturbation representative because it has complete true CellOracle displacement, scTenifoldKnk availability, an individual-TF interpretation, and lower AP-1-family redundancy; the empirical matrix rank is preserved."
            ),
        },
        "claim": "EGR1 was selected as the principal perturbation representative of the AP-1/CEBPB/EGR1 stress-transition programme.",
        "architecture_language": "overlapping regulatory phases; partially ordered regulatory architecture",
        "review_risk_flags": review_risks,
        "outputs": {
            "candidate_evidence_matrix": str(candidate_path.resolve()),
            "candidate_member_metrics": str(member_path.resolve()),
            "candidate_expression_by_state": str(expression_path.resolve()),
        },
        "caveat": "The candidate score supports transparent perturbation-target selection and does not establish EGR1 as a causal HCC driver or a unique module member.",
    }
    report_path = out_dir / "figure3a_candidate_selection_report.json"
    write_json(json_safe(report), report_path)
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = run(args.h5ad, args.out_dir)
    print(
        json.dumps(
            {
                "target_tf": report["target_tf"],
                "actual_rank": report["egr1_selection"]["actual_rank"],
                "highest_candidate": report["egr1_selection"]["highest_candidate"],
                "report": str((args.out_dir / "figure3a_candidate_selection_report.json").resolve()),
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

