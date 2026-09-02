import unittest
from pathlib import Path

import pandas as pd


class Figure1HRssLogicTests(unittest.TestCase):
    def test_target_specific_regulon_ranks_first_by_rss(self) -> None:
        script_path = Path("scripts/plot_figure1h_rss_ranking.py")
        self.assertTrue(script_path.exists())

        from scripts.plot_figure1h_rss_ranking import calculate_rss, select_top_regulons

        auc = pd.DataFrame(
            {
                "Target(+)": [0.90, 0.80, 0.01, 0.02],
                "Other(+)": [0.10, 0.10, 0.70, 0.80],
            },
            index=["cell_1", "cell_2", "cell_3", "cell_4"],
        )
        labels = pd.Series(["malignant", "malignant", "normal", "normal"], index=auc.index)

        rss = calculate_rss(auc, labels, min_cells_per_state=2)
        top = select_top_regulons(rss, target_state="malignant", top_n=1)

        self.assertEqual(top.iloc[0]["regulon"], "Target(+)")
        self.assertGreater(top.iloc[0]["rss"], 0.5)


if __name__ == "__main__":
    unittest.main()
