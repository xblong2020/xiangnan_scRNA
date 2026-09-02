from __future__ import annotations

import pandas as pd
from pathlib import Path

from scripts.audit_sctenifoldknk_reproducibility_v2 import (
    build_three_axis_validation_matrix,
    classify_rerun_decision,
    compute_gene_reproducibility,
)
from scripts.run_sctenifoldknk_reproducibility_audit_v2 import build_rerun_jobs, select_targets_to_rerun
from scripts.integrate_sctenifoldknk_reproducibility_audit_v2 import (
    build_axis_summary,
    compute_pathway_reproducibility,
    decide_panel_position,
)


def _seed_table(seed: int, genes: list[str], z_values: list[float]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "tf": ["TF1"] * len(genes),
            "gene": genes,
            "distance": [0.5, 0.4, 0.3][: len(genes)],
            "p.adj": [0.01, 0.02, 0.20][: len(genes)],
            "p.value": [0.001, 0.002, 0.20][: len(genes)],
            "Z": z_values,
            "FC": [2.0] * len(genes),
            "seed": [seed] * len(genes),
        }
    )


def test_compute_gene_reproducibility_reports_recurrence_and_direction() -> None:
    tables = {
        1: _seed_table(1, ["A", "B", "C"], [1.0, 1.0, -1.0]),
        2: _seed_table(2, ["A", "B", "D"], [2.0, 2.0, -1.0]),
        3: _seed_table(3, ["A", "B", "E"], [3.0, -3.0, -1.0]),
    }

    summary = compute_gene_reproducibility(tables, target_tf="TF1", fdr_threshold=0.05)

    assert int(summary.loc[summary["gene"].eq("A"), "n_significant_seeds"].iloc[0]) == 3
    assert int(summary.loc[summary["gene"].eq("B"), "n_significant_seeds"].iloc[0]) == 3
    assert int(summary.loc[summary["gene"].eq("C"), "n_significant_seeds"].iloc[0]) == 0
    assert summary.loc[summary["gene"].eq("A"), "direction_concordant"].iloc[0]
    assert not summary.loc[summary["gene"].eq("B"), "direction_concordant"].iloc[0]


def test_classify_rerun_decision_requires_all_three_contract_dimensions() -> None:
    needs_rerun = classify_rerun_decision(
        target_tf="HNF4A",
        nc_nnet=1,
        nc_ncells=100,
        seeds=[11],
        successful_seeds=1,
        input_cells=8098,
    )
    complete = classify_rerun_decision(
        target_tf="EGR1",
        nc_nnet=10,
        nc_ncells=500,
        seeds=[15071990, 15071991, 15071992],
        successful_seeds=3,
        input_cells=646,
    )

    assert needs_rerun["needs_rerun"] is True
    assert "nc_nNet" in " ".join(needs_rerun["reasons"])
    assert complete["needs_rerun"] is False


def test_validation_matrix_preserves_parameter_asymmetry_and_evidence_grade() -> None:
    summary = pd.DataFrame(
        {
            "tf": ["HNF4A", "EGR1", "SOX4"],
            "axis": ["Identity", "Stress", "Malignant state"],
            "nc_nNet": [10, 10, 10],
            "nc_nCells": [500, 500, 500],
            "seeds": ["1,2,3"] * 3,
            "successful_seeds": [3, 3, 3],
            "input_cells": [8098, 646, 552],
            "input_genes": [3000] * 3,
            "significant_genes": [10, 2, 8],
            "fdr_genes": [10, 2, 8],
            "three_of_three_recurrent_genes": [5, 1, 4],
            "two_of_three_recurrent_genes": [3, 1, 2],
            "gene_jaccard": [0.5, 0.1, 0.4],
            "direction_concordance": [0.9, 0.8, 0.85],
            "enriched_pathways": [2, 0, 3],
            "fdr_pathways": [1, 0, 2],
            "programme_level_recurrence": ["identity", "limited", "stress/malignant"],
            "reproducibility_grade": ["A", "C", "B"],
            "evidence_strength": ["strong", "limited", "moderate"],
            "main_figure_suitability": ["KEEP_MAIN", "MOVE_EXTENDED", "REVISE"],
            "extended_data_suitability": ["supplementary", "preferred", "required"],
            "limitation": ["", "", ""],
        }
    )

    matrix = build_three_axis_validation_matrix(summary)

    assert list(matrix.columns) == ["Dimension", "HNF4A", "EGR1", "SOX4"]
    assert matrix.loc[matrix["Dimension"].eq("nc_nNet"), "HNF4A"].iloc[0] == "10"
    assert matrix.loc[matrix["Dimension"].eq("Input cells"), "SOX4"].iloc[0] == "552"
    assert matrix.loc[matrix["Dimension"].eq("Main Figure suitability"), "EGR1"].iloc[0] == "MOVE_EXTENDED"


def test_build_rerun_jobs_only_targets_underpowered_historical_runs() -> None:
    audit = {
        "historical_runs": {
            "HNF4A": {"decision": {"needs_rerun": True}, "subset": "normal_reference"},
            "EGR1": {"decision": {"needs_rerun": False}, "subset": "stressed_regenerative"},
            "SOX4": {"decision": {"needs_rerun": True}, "subset": "malignant_like"},
        }
    }

    targets = select_targets_to_rerun(audit, requested_targets=["HNF4A", "EGR1", "SOX4"])
    jobs = build_rerun_jobs(
        project_root=Path("C:/project"),
        audit=audit,
        targets=targets,
        seeds=[15071990, 15071991, 15071992],
    )

    assert sorted({job["target_tf"] for job in jobs}) == ["HNF4A", "SOX4"]
    assert len(jobs) == 6
    assert all(job["nc_nNet"] == 10 and job["nc_nCells"] == 500 for job in jobs)
    assert all("sctenifoldknk_reproducibility_audit_v2" in str(job["output_dir"]) for job in jobs)


def test_compute_pathway_reproducibility_counts_seed_recurrence() -> None:
    seed1 = pd.DataFrame(
        {
            "tf": ["TF1", "TF1"],
            "seed": [1, 1],
            "database": ["KEGG", "KEGG"],
            "term": ["Pathway A", "Pathway B"],
            "p.adjust": [0.01, 0.20],
            "pvalue": [0.001, 0.2],
            "overlap_count": [3, 1],
            "programme_annotation": ["identity", ""],
        }
    )
    seed2 = seed1.copy()
    seed2["seed"] = 2
    seed2.loc[seed2["term"].eq("Pathway B"), "p.adjust"] = 0.01
    seed3 = seed1.copy()
    seed3["seed"] = 3

    summary = compute_pathway_reproducibility([seed1, seed2, seed3], target_tf="TF1")

    pathway_a = summary.loc[summary["term"].eq("Pathway A")].iloc[0]
    pathway_b = summary.loc[summary["term"].eq("Pathway B")].iloc[0]
    assert int(pathway_a["n_significant_seeds"]) == 3
    assert int(pathway_b["n_significant_seeds"]) == 1
    assert bool(pathway_a["recurrent_3_of_3"])
    assert not bool(pathway_b["recurrent_2_of_3"])


def test_build_axis_summary_assigns_conservative_grade_and_panel_position() -> None:
    seed_tables = {
        1: _seed_table(1, ["A", "B", "C"], [1.0, 1.0, -1.0]).assign(tf="HNF4A"),
        2: _seed_table(2, ["A", "B", "D"], [2.0, 2.0, -1.0]).assign(tf="HNF4A"),
        3: _seed_table(3, ["A", "B", "E"], [3.0, -3.0, -1.0]).assign(tf="HNF4A"),
    }
    pathway = pd.DataFrame(
        {
                "tf": ["HNF4A"],
            "seed": [1],
            "database": ["KEGG"],
            "term": ["Pathway A"],
            "p.adjust": [0.01],
            "pvalue": [0.001],
            "overlap_count": [3],
            "programme_annotation": ["identity"],
        }
    )
    summary = build_axis_summary(
        tf="HNF4A",
        axis="Identity",
        subset="normal_reference",
        seed_tables=seed_tables,
        pathway_results=[pathway.assign(seed=seed, **{"p.adjust": 0.01 if seed == 1 else 0.2}) for seed in [1, 2, 3]],
        input_cells=8098,
        input_genes=3000,
    )

    assert summary["reproducibility_grade"] == "B"
    assert decide_panel_position("HNF4A", summary)["figure2f"] == "MOVE_EXTENDED"
