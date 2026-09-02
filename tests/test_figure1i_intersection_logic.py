import unittest
from pathlib import Path

import pandas as pd


class Figure1IIntersectionLogicTests(unittest.TestCase):
    def test_figure_uses_ggsci_lancet_colors(self) -> None:
        from scripts.plot_figure1i_cellrank_cistarget_intersection import LANCET_COLORS

        self.assertEqual(LANCET_COLORS["cellrank"], "#00468B")
        self.assertEqual(LANCET_COLORS["cistarget"], "#42B540")
        self.assertEqual(LANCET_COLORS["shared"], "#ED0000")

    def test_intersection_keeps_significant_positive_tf_drivers(self) -> None:
        script_path = Path("scripts/plot_figure1i_cellrank_cistarget_intersection.py")
        self.assertTrue(script_path.exists())

        from scripts.plot_figure1i_cellrank_cistarget_intersection import compute_intersection

        cellrank = pd.DataFrame(
            {
                "gene": ["JUN", "FOS", "TPI1", "HNF4A"],
                "corr": [0.31, 0.28, 0.45, -0.20],
                "qval": [0.001, 0.002, 0.001, 0.001],
                "rank_positive_corr": [20, 40, 10, 60],
            }
        )
        regulons = pd.DataFrame(
            {
                "tf": ["JUN", "FOS", "TPI1", "HNF4A"],
                "regulon": ["JUN(+)", "FOS(+)", "TPI1(+)", "HNF4A(+)"],
                "best_nes": [4.0, 4.5, 5.0, 4.0],
                "n_targets": [20, 18, 12, 16],
                "n_motifs": [12, 8, 1, 6],
            }
        )

        result = compute_intersection(
            cellrank,
            regulons,
            known_tfs={"JUN", "FOS", "HNF4A"},
            top_n=50,
            qvalue_cutoff=0.05,
            min_corr=0.0,
            min_nes=3.0,
            min_motifs=3,
        )

        self.assertEqual(result["cellrank_tf_set"], {"JUN", "FOS"})
        self.assertEqual(result["cistarget_tf_set"], {"JUN", "FOS", "HNF4A"})
        self.assertEqual(result["overlap_set"], {"JUN", "FOS"})
        self.assertEqual(result["candidates"]["tf"].tolist(), ["JUN", "FOS"])


if __name__ == "__main__":
    unittest.main()
