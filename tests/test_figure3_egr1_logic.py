from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from scripts.audit_figure3e_egr1_determinism import run as audit_determinism
from scripts.figure3_egr1_common import (
    assign_fixed_stage,
    assert_figure3_path_isolated,
    choose_stress_transition_subset,
    compute_selection_score,
    select_significant_genes,
)


def test_fixed_stage_contract_uses_prespecified_boundaries() -> None:
    stages = assign_fixed_stage(pd.Series([0.0, 0.3299, 0.33, 0.6699, 0.67, 1.0]))
    assert stages.tolist() == ["early", "early", "intermediate", "intermediate", "late", "late"]


def test_stress_regenerative_is_selected_when_stressed_alone_is_too_small() -> None:
    audits = pd.DataFrame(
        [
            {"subset": "stressed_injured", "n_cells": 247, "n_datasets": 6, "egr1_detection_rate": 0.2, "max_dataset_fraction": 0.4},
            {"subset": "stressed_regenerative", "n_cells": 646, "n_datasets": 6, "egr1_detection_rate": 0.2, "max_dataset_fraction": 0.4},
            {"subset": "intermediate_pseudotime", "n_cells": 1600, "n_datasets": 6, "egr1_detection_rate": 0.2, "max_dataset_fraction": 0.4},
        ]
    )
    selected, reason = choose_stress_transition_subset(audits)
    assert selected == "stressed_regenerative"
    assert "prespecified priority" in reason


def test_significant_gene_selection_never_fills_with_non_significant_rows() -> None:
    table = pd.DataFrame(
        {
            "tf": ["EGR1"] * 4,
            "gene": ["EGR1", "A", "B", "C"],
            "distance": [9.0, 3.0, 5.0, 7.0],
            "p.adj": [0.001, 0.01, 0.20, 0.04],
        }
    )
    shown, total = select_significant_genes(table, top_n=20)
    assert total == 2
    assert shown["gene"].tolist() == ["C", "A"]
    assert shown["p.adj"].lt(0.05).all()
    assert "EGR1" not in shown["gene"].tolist()


def test_selection_score_preserves_empirical_ranking() -> None:
    frame = pd.DataFrame(
        {
            "candidate": ["EGR1", "CEBPB"],
            "celloracle_robustness": [0.8, 0.9],
            "sctenifoldknk_robustness": [0.8, 0.9],
            "transition_state_specificity": [0.8, 0.9],
            "temporal_positioning": [0.8, 0.9],
            "cross_dataset_stability": [0.8, 0.9],
            "pathway_interpretability": [0.8, 0.9],
            "generic_stress_penalty": [0.0, 0.0],
            "proliferation_dependency_penalty": [0.0, 0.0],
            "literature_overlap_penalty": [0.0, 0.0],
        }
    )
    scored = compute_selection_score(frame)
    assert scored.iloc[0]["candidate"] == "CEBPB"
    assert scored.loc[scored["candidate"].eq("EGR1"), "selection_rank"].item() == 2


def test_output_paths_cannot_cross_protected_tf_namespaces() -> None:
    assert_figure3_path_isolated(Path("metadata/driver/figure3c_egr1/report.json"))
    with pytest.raises(ValueError):
        assert_figure3_path_isolated(Path("metadata/driver/figure2c_sox4/report.json"))
    with pytest.raises(ValueError):
        assert_figure3_path_isolated(Path("metadata/driver/figure2c_hnf4a/report.json"))


def test_same_seed_audit_accepts_identical_gene_level_results(tmp_path: Path) -> None:
    frame = pd.DataFrame(
        {
            "tf": ["EGR1", "EGR1", "EGR1"],
            "gene": ["A", "B", "C"],
            "distance": [1.0, 2.0, 3.0],
            "p.adj": [0.01, 0.02, 0.03],
            "p.value": [0.001, 0.002, 0.003],
            "Z": [-2.0, 0.0, 2.0],
            "FC": [0.5, 1.0, 1.5],
        }
    )
    output = tmp_path / "figure3_egr1_determinism"
    output.mkdir()
    original = output / "original.tsv"
    repeat = output / "repeat.tsv"
    frame.to_csv(original, sep="\t", index=False)
    frame.to_csv(repeat, sep="\t", index=False)
    report = audit_determinism(original, repeat, output)
    assert report["status"] == "pass"
    assert report["numeric_values_reproducible"] is True
