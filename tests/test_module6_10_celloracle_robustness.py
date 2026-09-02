import unittest

import pandas as pd

from scripts.robustness_celloracle_perturbation_module6_10 import (
    aggregate_subset_scores,
    classify_control_tfs,
    compute_dataset_sample_dominance,
    compute_leave_one_group_out,
    summarize_rank_stability,
)


class Module610CellOracleRobustnessTest(unittest.TestCase):
    def _scores(self):
        return pd.DataFrame(
            {
                "tf": ["A", "A", "A", "B", "B", "B", "C", "C", "C"],
                "cell_id": [f"c{i}" for i in range(9)],
                "celloracle_state": ["malignant_or_malignant_like"] * 9,
                "dataset": ["D1", "D1", "D2", "D1", "D2", "D2", "D1", "D2", "D2"],
                "sample_id": ["S1", "S1", "S2", "S1", "S2", "S3", "S1", "S2", "S3"],
                "celloracle_main_strict": [True, True, False, True, False, False, True, False, False],
                "driver_main_strict__pseudotime_phase": ["early", "late", "late", "early", "middle", "late", "early", "middle", "late"],
                "cellrank_fate_prob_cnv_supported_malignant": [0.9] * 9,
                "malignant_fate_direction_cell_score": [3, 3, 2, 2, 2, 2, -1, -1, -1],
                "inner_product_cell_score": [3, 3, 2, 2, 2, 2, -1, -1, -1],
                "cnv_fate_probability_association_cell_score": [3, 3, 2, 2, 2, 2, -1, -1, -1],
                "hcc_malignant_module_rescue_cell_score": [3, 3, 2, 2, 2, 2, -1, -1, -1],
                "proliferation_module_rescue_cell_score": [3, 3, 2, 2, 2, 2, -1, -1, -1],
                "mature_hepatocyte_module_rescue_cell_score": [3, 3, 2, 2, 2, 2, -1, -1, -1],
                "module_rescue_cell_score": [3, 3, 2, 2, 2, 2, -1, -1, -1],
                "embedding_shift_norm": [1] * 9,
            }
        )

    def test_aggregate_subset_scores_ranks_tf_with_larger_rescue(self):
        summary = aggregate_subset_scores(self._scores(), subset_name="all", fate_high_threshold=0.5)

        self.assertEqual(list(summary["tf"]), ["A", "B", "C"])
        self.assertEqual(summary.loc[0, "subset"], "all")

    def test_leave_one_group_out_keeps_expected_groups(self):
        lodo = compute_leave_one_group_out(
            self._scores(),
            group_col="dataset",
            fate_high_threshold=0.5,
            min_remaining_cells=1,
        )

        self.assertEqual(set(lodo["left_out_group"]), {"D1", "D2"})
        self.assertIn("A", set(lodo["tf"]))

    def test_rank_stability_flags_stable_top_tf(self):
        subset_scores = pd.DataFrame(
            {
                "subset": ["s1", "s1", "s2", "s2"],
                "tf": ["A", "B", "A", "B"],
                "quantitative_rank": [1, 2, 1, 2],
                "quantitative_perturbation_score": [0.9, 0.1, 0.8, 0.2],
            }
        )

        stability = summarize_rank_stability(subset_scores).set_index("tf")

        self.assertEqual(int(stability.loc["A", "top5_fraction"]), 1)
        self.assertEqual(float(stability.loc["A", "median_rank"]), 1.0)

    def test_control_classification_selects_low_and_non_malignant_controls(self):
        tf_scores = pd.DataFrame(
            {
                "tf": ["A", "B", "C"],
                "quantitative_rank": [1, 2, 3],
                "quantitative_perturbation_score": [0.9, 0.5, 0.1],
                "malignant_fate_direction_score": [0.5, 0.1, -0.2],
                "state_specificity_ratio": [1.5, 1.2, 1.0],
            }
        )

        controls = classify_control_tfs(tf_scores, n_low=1)

        self.assertIn("low_score_tf", set(controls["control_type"]))
        self.assertIn("non_malignant_direction_tf", set(controls["control_type"]))

    def test_dataset_sample_dominance_reports_largest_fraction(self):
        dominance = compute_dataset_sample_dominance(self._scores(), group_cols=["dataset", "sample_id"])

        row = dominance.loc[(dominance["tf"] == "A") & (dominance["group_col"] == "dataset")].iloc[0]
        self.assertEqual(row["dominant_group"], "D1")
        self.assertAlmostEqual(float(row["dominant_fraction"]), 2 / 3)


if __name__ == "__main__":
    unittest.main()
