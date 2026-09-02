import unittest

import numpy as np
import pandas as pd

from scripts.visualize_celloracle_perturbation_module6_9 import (
    build_candidate_evidence_matrix,
    minmax_scale,
    pivot_state_projection,
    select_top_tfs,
)


class Module69CellOracleVisualizationLogicTest(unittest.TestCase):
    def test_minmax_scale_handles_constant_and_missing_values(self):
        values = pd.Series([2.0, 2.0, np.nan])

        scaled = minmax_scale(values)

        self.assertEqual(list(scaled), [0.5, 0.5, 0.0])

    def test_build_candidate_evidence_matrix_merges_inputs_and_scores(self):
        ranking = pd.DataFrame(
            {
                "tf": ["JUN", "FOS"],
                "rank": [1, 2],
                "anti_malignant_shift_score": [0.8, 0.2],
                "weighted_mean_abs_delta_x": [0.1, 0.3],
                "malignant_axis_projection_mean": [-0.8, -0.2],
            }
        )
        selection = pd.DataFrame(
            {
                "tf": ["JUN", "FOS"],
                "total_score": [90.0, 80.0],
                "motif_score": [20.0, 10.0],
                "fate_score": [25.0, 15.0],
                "cellrank_overlap_score": [10.0, 5.0],
                "robustness_score": [12.0, 8.0],
                "compatibility_score": [8.0, 4.0],
                "biology_score": [5.0, 4.0],
                "cnv_fate_pearson_r": [0.7, 0.3],
                "tf_self_cellrank_corr": [0.4, 0.1],
                "base_grn_target_genes": [100, 50],
                "selected_for_main_panel": [True, True],
            }
        )
        grn = pd.DataFrame(
            {
                "tf": ["JUN", "FOS"],
                "celloracle_state": ["malignant_or_malignant_like", "malignant_or_malignant_like"],
                "n_edges_passing_p": [500, 200],
                "mean_coef_abs_passing_p": [0.1, 0.2],
            }
        )

        evidence = build_candidate_evidence_matrix(ranking, selection, grn)

        self.assertEqual(list(evidence["tf"]), ["JUN", "FOS"])
        self.assertGreater(float(evidence.loc[0, "integrated_evidence_score"]), float(evidence.loc[1, "integrated_evidence_score"]))
        self.assertIn("anti_malignant_shift_score_scaled", evidence.columns)

    def test_pivot_state_projection_orders_states_and_tfs(self):
        state = pd.DataFrame(
            {
                "tf": ["FOS", "JUN", "JUN"],
                "celloracle_state": ["normal_reference", "malignant_or_malignant_like", "normal_reference"],
                "malignant_axis_projection_mean": [0.1, -0.5, 0.0],
            }
        )

        pivot = pivot_state_projection(state, tf_order=["JUN", "FOS"])

        self.assertEqual(list(pivot.index), ["JUN", "FOS"])
        self.assertEqual(list(pivot.columns), ["normal_reference", "malignant_or_malignant_like"])
        self.assertAlmostEqual(float(pivot.loc["JUN", "malignant_or_malignant_like"]), -0.5)

    def test_select_top_tfs_uses_integrated_score_then_rank(self):
        evidence = pd.DataFrame(
            {
                "tf": ["A", "B", "C"],
                "integrated_evidence_score": [0.2, 0.9, 0.9],
                "rank": [3, 2, 1],
            }
        )

        self.assertEqual(select_top_tfs(evidence, n=2), ["C", "B"])


if __name__ == "__main__":
    unittest.main()
