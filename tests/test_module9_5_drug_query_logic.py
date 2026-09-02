import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts.drug_query_module9_5 import (
    build_clue_query_payload,
    build_l1000fwd_payload,
    normalize_l1000fwd_results,
    parse_gmt,
)


class Module95DrugQueryLogicTest(unittest.TestCase):
    def test_parse_gmt_reads_signature_name_description_and_genes(self):
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "sig.gmt"
            path.write_text("sig_name\tdescription\tA\tB\tC\n", encoding="utf-8")

            parsed = parse_gmt(path)

        self.assertEqual(parsed["name"], "sig_name")
        self.assertEqual(parsed["description"], "description")
        self.assertEqual(parsed["genes"], ["A", "B", "C"])

    def test_l1000fwd_payload_uses_up_and_down_genes(self):
        payload = build_l1000fwd_payload(["HNF4A", "PPARA"], ["SOX4", "JUN"])

        self.assertEqual(payload["up_genes"], ["HNF4A", "PPARA"])
        self.assertEqual(payload["down_genes"], ["SOX4", "JUN"])

    def test_normalize_l1000fwd_results_labels_similar_as_candidate_and_opposite_as_counter(self):
        topn = {
            "similar": [
                {"sig_id": "SIG_A", "score": 0.91, "pvals": 0.001, "qvals": 0.01},
            ],
            "opposite": [
                {"sig_id": "SIG_B", "score": -0.88, "pvals": 0.002, "qvals": 0.02},
            ],
        }

        normalized = normalize_l1000fwd_results(topn)
        by_sig = normalized.set_index("sig_id")

        self.assertEqual(by_sig.loc["SIG_A", "candidate_direction"], "similar_to_reversal_signature")
        self.assertEqual(by_sig.loc["SIG_B", "candidate_direction"], "opposite_to_reversal_signature")
        self.assertGreater(float(by_sig.loc["SIG_A", "final_rank_score"]), float(by_sig.loc["SIG_B", "final_rank_score"]))

    def test_normalize_current_l1000fwd_format_parses_sig_id_metadata(self):
        topn = {
            "similar": [
                {
                    "sig_id": "CPC013_VCAP_24H:BRD-K76938712-001-01-2:10",
                    "scores": 0.083,
                    "zscores": -1.68,
                    "pvals": 3.7e-8,
                    "qvals": 1.9e-4,
                }
            ],
            "opposite": [],
        }

        row = normalize_l1000fwd_results(topn).iloc[0]

        self.assertEqual(row["compound_id"], "BRD-K76938712")
        self.assertEqual(row["cell_line"], "VCAP")
        self.assertEqual(row["time"], "24H")
        self.assertEqual(row["dose"], "10")
        self.assertAlmostEqual(float(row["raw_score"]), 0.083)

    def test_clue_payload_contains_stringified_gmt_and_query_metadata(self):
        payload = build_clue_query_payload(["HNF4A", "PPARA"], ["SOX4", "JUN"], query_name="module9_5")

        self.assertEqual(payload["tool_id"], "sig_query")
        self.assertIn("HNF4A", payload["uptag"])
        self.assertIn("SOX4", payload["dntag"])
        self.assertEqual(payload["name"], "module9_5")


if __name__ == "__main__":
    unittest.main()
