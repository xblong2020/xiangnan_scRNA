import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd


class Figure1GTsneLogicTests(unittest.TestCase):
    def test_tsne_defaults_are_bounded_for_the_full_cytotrace2_cell_set(self) -> None:
        from scripts.compute_figure1g_hepatocyte_tsne import parse_args

        with patch("sys.argv", ["compute_figure1g_hepatocyte_tsne.py"]):
            args = parse_args()

        self.assertEqual(args.max_iter, 500)
        self.assertEqual(args.angle, 0.7)

    def test_select_latent_rows_preserves_score_cell_order(self) -> None:
        script_path = Path("scripts/compute_figure1g_hepatocyte_tsne.py")
        self.assertTrue(script_path.exists())

        from scripts.compute_figure1g_hepatocyte_tsne import select_latent_rows

        latent = pd.DataFrame(
            {"SCVI_1": [3.0, 1.0, 2.0], "SCVI_2": [30.0, 10.0, 20.0]},
            index=pd.Index(["cell_3", "cell_1", "cell_2"], name="cell_id"),
        )
        selected = select_latent_rows(latent, ["cell_2", "cell_1"])

        self.assertEqual(selected.index.tolist(), ["cell_2", "cell_1"])
        self.assertEqual(selected["SCVI_1"].tolist(), [2.0, 1.0])


if __name__ == "__main__":
    unittest.main()
