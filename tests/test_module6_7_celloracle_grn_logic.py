import unittest

import pandas as pd

from scripts.fit_celloracle_grn_module6_7 import (
    merge_links_dict,
    summarize_grn_links,
    summarize_tf_network,
)


class Module67CellOracleGrnLogicTest(unittest.TestCase):
    def test_merge_links_dict_adds_state_and_standard_columns(self):
        links = {
            "state_a": pd.DataFrame(
                {
                    "source": ["JUN", "FOS"],
                    "target": ["A", "B"],
                    "coef_mean": [0.2, -0.3],
                    "coef_abs": [0.2, 0.3],
                    "p": [0.001, 0.2],
                    "-logp": [3.0, 0.7],
                }
            )
        }

        merged = merge_links_dict(links)

        self.assertEqual(list(merged["celloracle_state"].unique()), ["state_a"])
        self.assertIn("coef_abs", merged.columns)
        self.assertEqual(len(merged), 2)

    def test_summarize_grn_links_counts_edges_per_state(self):
        merged = pd.DataFrame(
            {
                "celloracle_state": ["a", "a", "b"],
                "source": ["JUN", "FOS", "JUN"],
                "target": ["X", "Y", "Z"],
                "coef_mean": [0.1, -0.2, 0.3],
                "coef_abs": [0.1, 0.2, 0.3],
                "p": [0.001, 0.02, 0.5],
            }
        )

        summary = summarize_grn_links(merged, p_threshold=0.05).set_index("celloracle_state")

        self.assertEqual(int(summary.loc["a", "n_edges_total"]), 2)
        self.assertEqual(int(summary.loc["a", "n_edges_passing_p"]), 2)
        self.assertEqual(int(summary.loc["b", "n_edges_passing_p"]), 0)

    def test_summarize_tf_network_tracks_selected_tf_edges_and_sign(self):
        merged = pd.DataFrame(
            {
                "celloracle_state": ["a", "a", "a", "b"],
                "source": ["JUN", "JUN", "FOS", "JUN"],
                "target": ["X", "Y", "Z", "Q"],
                "coef_mean": [0.4, -0.1, 0.2, -0.5],
                "coef_abs": [0.4, 0.1, 0.2, 0.5],
                "p": [0.001, 0.2, 0.01, 0.03],
            }
        )

        summary = summarize_tf_network(merged, input_tfs=["JUN", "FOS"], p_threshold=0.05).set_index(
            ["celloracle_state", "tf"]
        )

        self.assertEqual(int(summary.loc[("a", "JUN"), "n_edges_total"]), 2)
        self.assertEqual(int(summary.loc[("a", "JUN"), "n_edges_passing_p"]), 1)
        self.assertEqual(int(summary.loc[("a", "JUN"), "n_positive_edges_passing_p"]), 1)
        self.assertEqual(int(summary.loc[("b", "JUN"), "n_negative_edges_passing_p"]), 1)
        self.assertEqual(int(summary.loc[("a", "FOS"), "n_edges_passing_p"]), 1)


if __name__ == "__main__":
    unittest.main()
