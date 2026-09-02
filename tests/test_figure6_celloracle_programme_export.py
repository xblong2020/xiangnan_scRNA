import numpy as np
import pandas as pd

from scripts.figure6_00_export_celloracle_programme_deltas import (
    build_frozen_programmes,
    score_programmes,
    select_sox4_targets,
)


def test_select_sox4_targets_is_ranked_and_deduplicated():
    tf_targets = pd.DataFrame(
        {"tf": ["SOX4", "SOX4", "HNF4A"], "gene": ["G2", "G1", "X"], "rank": [2, 1, 1]}
    )
    state = pd.DataFrame(
        {
            "tf": ["SOX4", "SOX4"],
            "gene": ["G2", "G3"],
            "malignant_like_specificity_ratio": [4.0, 3.0],
            "malignant_like_fdr": [0.01, 0.02],
        }
    )
    assert select_sox4_targets(tf_targets, state, top_n=50) == ["G1", "G2", "G3"]


def test_score_programmes_uses_available_frozen_genes_only():
    matrix = np.array([[1.0, 3.0, 9.0], [2.0, 4.0, 8.0]])
    scores, manifest = score_programmes(matrix, ["A", "B", "C"], {"identity": ["A", "B", "MISSING"]})
    np.testing.assert_allclose(scores["identity"], [2.0, 3.0])
    assert manifest[0]["n_requested"] == 3
    assert manifest[0]["n_available"] == 2
    assert manifest[0]["genes_missing"] == "MISSING"


def test_build_frozen_programmes_preserves_module5_gene_sets():
    availability = pd.DataFrame(
        {
            "run_id": ["main_strict"] * 4,
            "module": ["Mature_Hepatocyte", "Stressed_Injured", "Proliferation", "HCC_Malignant_Associated"],
            "genes_available": ["ALB;HNF4A", "FOS;JUN", "PCNA", "SPP1"],
        }
    )
    targets = pd.DataFrame({"tf": ["SOX4"], "gene": ["TG"], "rank": [1]})
    state = pd.DataFrame(columns=["tf", "gene"])
    programmes = build_frozen_programmes(availability, targets, state)
    assert programmes["identity_program_change"] == ["ALB", "HNF4A"]
    assert programmes["stress_transition_change"] == ["FOS", "JUN"]
    assert programmes["sox4_programme_change"] == ["SOX4", "TG"]

