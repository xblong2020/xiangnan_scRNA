import unittest

import numpy as np
import pandas as pd

from scripts.score_celloracle_perturbation_module6_9b import (
    aggregate_tf_scores,
    compute_local_gradients,
    compute_state_specificity,
    minmax_scale,
    score_cell_level_vectors,
)


class Module69bCellOracleQuantitativeScoringTest(unittest.TestCase):
    def test_compute_local_gradients_recovers_linear_field(self):
        embedding = np.array(
            [
                [0.0, 0.0],
                [1.0, 0.0],
                [0.0, 1.0],
                [1.0, 1.0],
                [2.0, 1.0],
            ]
        )
        values = 2 * embedding[:, 0] - embedding[:, 1]

        gradients = compute_local_gradients(embedding, values, k=4)

        np.testing.assert_allclose(np.nanmean(gradients, axis=0), np.array([2.0, -1.0]), atol=1e-8)

    def test_score_cell_level_vectors_uses_gradient_inner_products(self):
        cell_shift = pd.DataFrame(
            {
                "tf": ["JUN", "JUN"],
                "cell_id": ["c1", "c2"],
                "celloracle_state": ["malignant_or_malignant_like", "normal_reference"],
                "delta_embedding_1": [-1.0, 1.0],
                "delta_embedding_2": [0.0, 0.0],
                "embedding_shift_norm": [1.0, 1.0],
                "malignant_axis_projection": [-1.0, 1.0],
            }
        )
        metadata = pd.DataFrame(
            {
                "cell_id": ["c1", "c2"],
                "cellrank_fate_prob_cnv_supported_malignant": [0.9, 0.1],
                "driver_main_strict__pseudotime_mean": [0.8, 0.2],
                "driver_main_strict__module_HCC_Malignant_Associated": [1.0, 0.0],
                "driver_main_strict__module_Proliferation": [0.5, 0.0],
                "driver_main_strict__module_Mature_Hepatocyte": [0.0, 1.0],
            }
        )
        gradients = {
            "cnv_fate": np.array([[1.0, 0.0], [1.0, 0.0]]),
            "pseudotime": np.array([[1.0, 0.0], [1.0, 0.0]]),
            "hcc_module": np.array([[1.0, 0.0], [1.0, 0.0]]),
            "proliferation_module": np.array([[1.0, 0.0], [1.0, 0.0]]),
            "mature_module": np.array([[-1.0, 0.0], [-1.0, 0.0]]),
        }

        scored = score_cell_level_vectors(cell_shift, metadata, gradients)

        first = scored.iloc[0]
        self.assertAlmostEqual(float(first["cnv_fate_probability_association_cell_score"]), 1.0)
        self.assertAlmostEqual(float(first["inner_product_cell_score"]), 1.0)
        self.assertAlmostEqual(float(first["module_rescue_cell_score"]), 1.0)

    def test_compute_state_specificity_rewards_malignant_specific_effect(self):
        cell_scores = pd.DataFrame(
            {
                "tf": ["JUN", "JUN", "FOS", "FOS"],
                "celloracle_state": [
                    "malignant_or_malignant_like",
                    "normal_reference",
                    "malignant_or_malignant_like",
                    "normal_reference",
                ],
                "malignant_fate_direction_cell_score": [2.0, 0.1, 1.0, 1.0],
                "embedding_shift_norm": [2.0, 0.1, 1.0, 1.0],
            }
        )

        specificity = compute_state_specificity(cell_scores).set_index("tf")

        self.assertGreater(float(specificity.loc["JUN", "state_specificity_ratio"]), 1.0)
        self.assertAlmostEqual(float(specificity.loc["FOS", "state_specificity_ratio"]), 1.0)
        self.assertGreater(float(specificity.loc["JUN", "state_specificity_score"]), float(specificity.loc["FOS", "state_specificity_score"]))

    def test_aggregate_tf_scores_builds_composite_rank(self):
        cell_scores = pd.DataFrame(
            {
                "tf": ["JUN", "JUN", "FOS", "FOS"],
                "celloracle_state": [
                    "malignant_or_malignant_like",
                    "normal_reference",
                    "malignant_or_malignant_like",
                    "normal_reference",
                ],
                "cellrank_fate_prob_cnv_supported_malignant": [0.9, 0.1, 0.9, 0.1],
                "malignant_fate_direction_cell_score": [2.0, 0.0, 0.5, 0.0],
                "inner_product_cell_score": [2.0, 0.0, 0.5, 0.0],
                "cnv_fate_probability_association_cell_score": [2.0, 0.0, 0.5, 0.0],
                "hcc_malignant_module_rescue_cell_score": [2.0, 0.0, 0.5, 0.0],
                "proliferation_module_rescue_cell_score": [2.0, 0.0, 0.5, 0.0],
                "mature_hepatocyte_module_rescue_cell_score": [2.0, 0.0, 0.5, 0.0],
                "module_rescue_cell_score": [2.0, 0.0, 0.5, 0.0],
                "embedding_shift_norm": [1.0, 0.1, 1.0, 0.1],
            }
        )

        summary = aggregate_tf_scores(cell_scores, fate_high_threshold=0.5)

        self.assertEqual(list(summary["tf"]), ["JUN", "FOS"])
        self.assertGreater(float(summary.loc[0, "quantitative_perturbation_score"]), float(summary.loc[1, "quantitative_perturbation_score"]))

    def test_minmax_scale_constant_values(self):
        scaled = minmax_scale(pd.Series([3.0, 3.0]))

        self.assertEqual(list(scaled), [0.5, 0.5])


if __name__ == "__main__":
    unittest.main()
