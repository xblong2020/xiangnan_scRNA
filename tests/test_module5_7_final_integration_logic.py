import unittest
from tempfile import TemporaryDirectory
from pathlib import Path

import pandas as pd

from scripts.final_trajectory_evidence_module5_7 import (
    build_conclusion_table,
    classify_final_evidence,
    first_last_delta,
    make_markdown_table,
)


class Module57FinalIntegrationLogicTest(unittest.TestCase):
    def test_first_last_delta_uses_ordered_bins(self):
        df = pd.DataFrame(
            {
                "pseudotime_bin": [2, 0, 1],
                "cnv_supported_fraction": [0.8, 0.1, 0.3],
            }
        )

        result = first_last_delta(df, "cnv_supported_fraction")

        self.assertAlmostEqual(result["early_value"], 0.1)
        self.assertAlmostEqual(result["late_value"], 0.8)
        self.assertAlmostEqual(result["delta"], 0.7)

    def test_classify_final_evidence_marks_supported_with_sample_caveat(self):
        row = {
            "cnv_supported_delta": 0.6,
            "hcc_malignant_module_delta": 1.0,
            "proliferation_module_delta": 0.4,
            "sample_robustness_label": "overall_positive_group_mixed",
            "dataset_robustness_label": "robust_positive",
        }

        self.assertEqual(classify_final_evidence(row), "supported_with_sample_composition_caveat")

    def test_classify_final_evidence_marks_consensus_supported(self):
        row = {
            "cnv_supported_delta": 0.5,
            "hcc_malignant_module_delta": 0.8,
            "proliferation_module_delta": 0.2,
            "sample_robustness_label": "robust_positive",
            "dataset_robustness_label": "robust_positive",
        }

        self.assertEqual(classify_final_evidence(row), "consensus_supported")

    def test_build_conclusion_table_merges_bin_and_robustness_evidence(self):
        bins = pd.DataFrame(
            {
                "run_id": ["main", "main"],
                "method": ["m1", "m1"],
                "pseudotime_bin": [0, 1],
                "cnv_supported_fraction": [0.0, 0.6],
                "review_fraction": [0.0, 0.1],
                "mean_hcc_malignant_module": [0.2, 1.2],
                "mean_proliferation_module": [0.1, 0.5],
            }
        )
        robustness = pd.DataFrame(
            {
                "run_id": ["main", "main"],
                "method": ["m1", "m1"],
                "group_type": ["sample_id", "dataset"],
                "feature": ["module3_cnv_supported", "HCC_Malignant_Associated"],
                "robustness_label": ["overall_positive_group_mixed", "robust_positive"],
                "positive_group_fraction": [0.25, 1.0],
                "min_loo_delta": [0.2, 0.8],
            }
        )
        correlations = pd.DataFrame(
            {
                "run_id": ["main"],
                "method": ["m1"],
                "feature": ["module3_cnv_supported"],
                "spearman_rho": [0.42],
            }
        )

        table = build_conclusion_table(bins, robustness, correlations)

        self.assertEqual(table.shape[0], 1)
        self.assertAlmostEqual(table.loc[0, "cnv_supported_delta"], 0.6)
        self.assertEqual(table.loc[0, "final_evidence_label"], "supported_with_sample_composition_caveat")
        self.assertAlmostEqual(table.loc[0, "module3_cnv_supported_spearman_rho"], 0.42)

    def test_make_markdown_table_does_not_require_optional_tabulate(self):
        conclusion = pd.DataFrame(
            {
                "run_id": ["main"],
                "method": ["m1"],
                "final_evidence_label": ["consensus_supported"],
                "cnv_supported_delta": [0.6],
                "hcc_malignant_module_delta": [1.0],
                "proliferation_module_delta": [0.4],
                "sample_robustness_label": ["robust_positive"],
                "dataset_robustness_label": ["robust_positive"],
                "sample_positive_group_fraction": [0.8],
                "sample_centered_cnv_spearman_rho": [0.2],
            }
        )
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "table.md"

            make_markdown_table(conclusion, path)

            text = path.read_text(encoding="utf-8")
        self.assertIn("final_evidence_label", text)
        self.assertIn("consensus_supported", text)


if __name__ == "__main__":
    unittest.main()
