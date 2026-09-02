import unittest

import numpy as np
import pandas as pd

from scripts.analyze_trajectory_trends_module5_4 import (
    assign_pseudotime_bins,
    compute_module_scores,
    summarize_trend,
)


class Module54TrendLogicTest(unittest.TestCase):
    def test_assign_pseudotime_bins_orders_quantiles(self):
        values = pd.Series(np.linspace(0, 1, 10), name="pt")
        bins = assign_pseudotime_bins(values, n_bins=5)

        self.assertEqual(sorted(bins.dropna().unique().tolist()), [0, 1, 2, 3, 4])
        self.assertEqual(int(bins.iloc[0]), 0)
        self.assertEqual(int(bins.iloc[-1]), 4)

    def test_compute_module_scores_uses_only_available_genes(self):
        expr = pd.DataFrame(
            {
                "ALB": [2.0, 4.0],
                "APOA1": [4.0, 8.0],
                "MKI67": [0.0, 10.0],
            },
            index=["cell_a", "cell_b"],
        )
        panels = {"Mature": ["ALB", "APOA1", "MISSING"], "Cycle": ["MKI67"]}

        scores, availability = compute_module_scores(expr, panels)

        self.assertEqual(scores.loc["cell_a", "Mature"], 3.0)
        self.assertEqual(scores.loc["cell_b", "Mature"], 6.0)
        self.assertEqual(availability["Mature"]["n_available"], 2)
        self.assertEqual(availability["Mature"]["genes_available"], "ALB;APOA1")

    def test_summarize_trend_detects_increasing_expression(self):
        pseudotime = pd.Series([0, 1, 2, 3, 4], dtype=float)
        values = pd.Series([1, 2, 3, 4, 5], dtype=float)

        summary = summarize_trend(values, pseudotime)

        self.assertAlmostEqual(summary["spearman_rho"], 1.0)
        self.assertGreater(summary["delta_last_first_bin"], 0)
        self.assertEqual(summary["trend_direction"], "increasing")


if __name__ == "__main__":
    unittest.main()
