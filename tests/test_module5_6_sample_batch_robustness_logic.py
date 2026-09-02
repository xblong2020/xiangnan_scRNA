import unittest

import pandas as pd

from scripts.sample_batch_robustness_module5_6 import (
    batch_centered_spearman,
    largest_group_fraction,
    leave_one_group_out_delta,
    summarize_group_trends,
)


class Module56SampleBatchRobustnessLogicTest(unittest.TestCase):
    def test_summarize_group_trends_reports_late_minus_early_delta(self):
        df = pd.DataFrame(
            {
                "sample": ["A", "A", "A", "A", "B", "B", "B", "B"],
                "pseudotime": [0.1, 0.2, 0.8, 0.9, 0.1, 0.2, 0.8, 0.9],
                "cnv_supported": [False, False, True, True, False, True, True, True],
            }
        )

        summary = summarize_group_trends(df, "sample", "pseudotime", "cnv_supported", n_bins=2, min_cells=3)

        self.assertEqual(set(summary["group"]), {"A", "B"})
        self.assertAlmostEqual(float(summary.loc[summary["group"].eq("A"), "late_minus_early_delta"].iloc[0]), 1.0)
        self.assertAlmostEqual(float(summary.loc[summary["group"].eq("B"), "late_minus_early_delta"].iloc[0]), 0.5)

    def test_leave_one_group_out_delta_tracks_dependence_on_each_sample(self):
        df = pd.DataFrame(
            {
                "sample": ["A", "A", "A", "A", "B", "B", "B", "B"],
                "pseudotime": [0.1, 0.2, 0.8, 0.9, 0.1, 0.2, 0.8, 0.9],
                "cnv_supported": [False, False, True, True, False, False, False, False],
            }
        )

        loo = leave_one_group_out_delta(df, "sample", "pseudotime", "cnv_supported", n_bins=2, min_remaining_cells=3)

        self.assertAlmostEqual(float(loo.loc[loo["omitted_group"].eq("B"), "loo_delta"].iloc[0]), 1.0)
        self.assertAlmostEqual(float(loo.loc[loo["omitted_group"].eq("A"), "loo_delta"].iloc[0]), 0.0)
        self.assertLess(float(loo.loc[loo["omitted_group"].eq("A"), "delta_shift_from_overall"].iloc[0]), 0.0)

    def test_batch_centered_spearman_uses_within_batch_signal(self):
        df = pd.DataFrame(
            {
                "batch": ["b1", "b1", "b1", "b1", "b2", "b2", "b2", "b2"],
                "pseudotime": [0.1, 0.2, 0.8, 0.9, 0.1, 0.2, 0.8, 0.9],
                "score": [0.0, 0.1, 1.0, 1.1, 5.0, 5.1, 6.0, 6.1],
            }
        )

        result = batch_centered_spearman(df, "batch", "pseudotime", "score", min_cells=3)

        self.assertEqual(result["n_batches"], 2)
        self.assertGreater(result["spearman_rho"], 0.9)

    def test_largest_group_fraction_reports_dominance(self):
        df = pd.DataFrame({"sample": ["A", "A", "A", "B"]})

        result = largest_group_fraction(df, "sample")

        self.assertEqual(result["n_groups"], 2)
        self.assertAlmostEqual(result["largest_group_fraction"], 0.75)


if __name__ == "__main__":
    unittest.main()
