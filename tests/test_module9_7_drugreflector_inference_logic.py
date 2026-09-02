import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from scripts.drugreflector_inference_module9_7 import (
    build_consensus_table,
    build_gene_coverage_rows,
    flatten_top_compounds,
    read_vscore_series,
)


class Module97DrugReflectorInferenceLogicTest(unittest.TestCase):
    def test_read_vscore_series_normalizes_and_averages_duplicate_genes(self):
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "vscore.tsv"
            pd.DataFrame(
                [
                    {"gene": " hnf4a ", "v_score": 1.0},
                    {"gene": "HNF4A", "v_score": 0.5},
                    {"gene": "SOX4", "v_score": -1.0},
                ]
            ).to_csv(path, sep="\t", index=False)

            series = read_vscore_series(path, "test_signature")

        self.assertEqual(series.name, "test_signature")
        self.assertEqual(set(series.index), {"HNF4A", "SOX4"})
        self.assertAlmostEqual(float(series.loc["HNF4A"]), 0.75)

    def test_build_gene_coverage_rows_reports_fold_and_union_overlap(self):
        series = pd.Series([1.0, -1.0, 0.5], index=["HNF4A", "SOX4", "ALB"], name="sig")
        rows = build_gene_coverage_rows(
            "primary",
            series,
            [{"HNF4A", "JUN"}, {"SOX4", "FOS"}],
        )

        union = [row for row in rows if row["fold"] == "union"][0]
        self.assertEqual(union["n_query_genes"], 3)
        self.assertEqual(union["n_model_landmark_genes"], 4)
        self.assertEqual(union["n_overlap_genes"], 2)
        self.assertAlmostEqual(union["query_gene_coverage_fraction"], 2 / 3)

    def test_flatten_top_compounds_adds_one_based_rank(self):
        top = {
            "sig": pd.DataFrame(
                {
                    "compound": ["drug_a", "drug_b"],
                    "rank": [1, 0],
                    "logit": [0.1, 0.8],
                    "prob": [0.2, 0.7],
                }
            )
        }

        frame = flatten_top_compounds(top, "primary")

        self.assertEqual(frame["compound"].tolist(), ["drug_b", "drug_a"])
        self.assertEqual(frame["rank_1based"].tolist(), [1, 2])
        self.assertEqual(frame["source_label"].unique().tolist(), ["primary"])

    def test_build_consensus_table_prioritizes_compounds_in_both_lists(self):
        primary = pd.DataFrame(
            {
                "compound": ["drug_a", "drug_b"],
                "rank_0based": [0, 1],
                "rank_1based": [1, 2],
                "logit": [2.0, 1.0],
                "prob": [0.7, 0.2],
            }
        )
        sensitivity = pd.DataFrame(
            {
                "compound": ["drug_b", "drug_c"],
                "rank_0based": [0, 1],
                "rank_1based": [1, 2],
                "logit": [3.0, 1.5],
                "prob": [0.8, 0.1],
            }
        )

        consensus = build_consensus_table(primary, sensitivity)

        self.assertEqual(consensus.iloc[0]["compound"], "drug_b")
        self.assertTrue(bool(consensus.iloc[0]["in_both_top_lists"]))
        self.assertEqual(int(consensus.iloc[0]["primary_rank_1based"]), 2)
        self.assertEqual(int(consensus.iloc[0]["sensitivity_rank_1based"]), 1)


if __name__ == "__main__":
    unittest.main()
