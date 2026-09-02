import unittest

import numpy as np
import pandas as pd

from scripts.simulate_celloracle_perturbation_module6_8 import (
    build_perturbation_condition,
    compute_malignant_axis,
    compute_sparse_embedding_shift,
    select_available_tfs,
    summarize_cell_shifts,
    summarize_gene_delta_by_state,
    summarize_perturbation_ranking,
)


class Module68CellOraclePerturbationLogicTest(unittest.TestCase):
    def test_select_available_tfs_reports_missing_and_keeps_order(self):
        selected, missing = select_available_tfs(
            input_tfs=["JUN", "FOS", "MYC"],
            available_genes={"JUN", "MYC"},
            active_regulatory_genes={"JUN"},
        )

        self.assertEqual(selected, ["JUN"])
        self.assertEqual(missing, ["FOS", "MYC"])

    def test_build_perturbation_condition_defaults_to_knockout(self):
        condition = build_perturbation_condition("JUN", mode="knockout", expression_value=None)

        self.assertEqual(condition, {"JUN": 0.0})

    def test_compute_malignant_axis_uses_normal_to_malignant_centroids(self):
        embedding = np.array([[0.0, 0.0], [0.0, 2.0], [3.0, 0.0], [3.0, 2.0]])
        states = pd.Series(["normal_reference", "normal_reference", "malignant_or_malignant_like", "malignant_or_malignant_like"])

        axis = compute_malignant_axis(
            embedding,
            states,
            start_state="normal_reference",
            end_state="malignant_or_malignant_like",
        )

        np.testing.assert_allclose(axis, np.array([1.0, 0.0]))

    def test_compute_sparse_embedding_shift_matches_neighbor_weighting(self):
        embedding = np.array([[0.0, 0.0], [1.0, 0.0]])
        corrcoef = np.array([[0.0, 1.0], [0.0, 0.0]])
        indices = np.array([[1], [0]])

        shift = compute_sparse_embedding_shift(
            embedding=embedding,
            corrcoef=corrcoef,
            neighbor_indices=indices,
            sigma_corr=0.05,
        )

        np.testing.assert_allclose(shift, np.zeros((2, 2)), atol=1e-8)

    def test_summarize_cell_shifts_projects_to_malignant_axis(self):
        obs = pd.DataFrame({"celloracle_state": ["a", "a", "b"]}, index=["c1", "c2", "c3"])
        delta_embedding = np.array([[1.0, 0.0], [0.0, 2.0], [-1.0, 0.0]])
        delta_x = np.array([[1.0, -1.0], [0.0, 2.0], [3.0, 4.0]])

        summary = summarize_cell_shifts(
            tf="JUN",
            obs=obs,
            delta_embedding=delta_embedding,
            delta_x=delta_x,
            malignant_axis=np.array([1.0, 0.0]),
        )

        self.assertEqual(list(summary["cell_id"]), ["c1", "c2", "c3"])
        self.assertAlmostEqual(float(summary.loc[0, "embedding_shift_norm"]), 1.0)
        self.assertAlmostEqual(float(summary.loc[2, "malignant_axis_projection"]), -1.0)
        self.assertAlmostEqual(float(summary.loc[0, "mean_abs_delta_x"]), 1.0)

    def test_summarize_gene_delta_by_state_returns_top_genes(self):
        delta_x = np.array([[2.0, -1.0, 0.0], [4.0, -3.0, 1.0], [0.0, 5.0, -6.0]])
        states = pd.Series(["a", "a", "b"])

        summary = summarize_gene_delta_by_state(
            tf="JUN",
            genes=["G1", "G2", "G3"],
            states=states,
            delta_x=delta_x,
            top_n=1,
        )

        summary = summary.set_index(["celloracle_state", "gene"])
        self.assertAlmostEqual(float(summary.loc[("a", "G1"), "mean_delta_x"]), 3.0)
        self.assertAlmostEqual(float(summary.loc[("b", "G3"), "mean_delta_x"]), -6.0)

    def test_summarize_perturbation_ranking_prefers_anti_malignant_shift(self):
        state_summary = pd.DataFrame(
            {
                "tf": ["JUN", "FOS"],
                "celloracle_state": ["malignant_or_malignant_like", "malignant_or_malignant_like"],
                "n_cells": [10, 10],
                "mean_abs_delta_x_mean": [0.1, 0.2],
                "embedding_shift_norm_mean": [0.3, 0.4],
                "malignant_axis_projection_mean": [-0.5, -0.2],
            }
        )

        ranking = summarize_perturbation_ranking(state_summary)

        self.assertEqual(list(ranking["tf"]), ["JUN", "FOS"])
        self.assertAlmostEqual(float(ranking.loc[0, "anti_malignant_shift_score"]), 0.5)


if __name__ == "__main__":
    unittest.main()
