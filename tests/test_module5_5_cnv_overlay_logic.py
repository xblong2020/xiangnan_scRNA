import unittest

import pandas as pd

from scripts.overlay_cnv_malignant_evidence_module5_5 import (
    assign_cnv_evidence_tier,
    bin_cnv_overlay,
    numeric_fraction,
)


class Module55CnvOverlayLogicTest(unittest.TestCase):
    def test_assign_cnv_evidence_tier_prioritizes_module3_calls(self):
        self.assertEqual(
            assign_cnv_evidence_tier(
                {
                    "malignant_hcc_call": "malignant_hcc_high_conf",
                    "copykat_status": "diploid",
                    "trajectory_root_end_role": "intermediate_trajectory",
                }
            ),
            "module3_high_conf_malignant",
        )
        self.assertEqual(
            assign_cnv_evidence_tier(
                {
                    "malignant_hcc_call": "malignant_hcc_cnv_support",
                    "copykat_status": "aneuploid",
                    "trajectory_root_end_role": "end_malignant_cnv",
                }
            ),
            "module3_cnv_supported_malignant",
        )

    def test_assign_cnv_evidence_tier_marks_review_and_reference(self):
        self.assertEqual(
            assign_cnv_evidence_tier(
                {
                    "malignant_hcc_call": "malignant_hcc_marker_proliferation_needs_cnv_review",
                    "copykat_status": "not.defined",
                    "trajectory_root_end_role": "end_malignant_review",
                }
            ),
            "malignant_like_needs_review",
        )
        self.assertEqual(
            assign_cnv_evidence_tier(
                {
                    "malignant_hcc_call": "not_module3_candidate",
                    "copykat_status": "Unknown",
                    "trajectory_root_end_role": "root_reference",
                }
            ),
            "no_cnv_evidence_or_reference",
        )

    def test_numeric_fraction_handles_string_booleans(self):
        values = pd.Series(["true", "False", True, "yes", "0"])
        self.assertEqual(numeric_fraction(values), 0.6)

    def test_bin_cnv_overlay_summarizes_fraction_and_scores(self):
        df = pd.DataFrame(
            {
                "pseudotime": [0.1, 0.2, 0.8, 0.9],
                "cnv_evidence_tier": [
                    "no_cnv_evidence_or_reference",
                    "module3_cnv_supported_malignant",
                    "module3_cnv_supported_malignant",
                    "module3_high_conf_malignant",
                ],
                "copykat_aneuploid": [False, True, True, True],
                "cnv_proxy_burden": [0.0, 0.5, 0.7, 0.8],
                "HCC_Malignant_Associated": [0.0, 1.0, 2.0, 3.0],
            }
        )

        summary = bin_cnv_overlay(df, pseudotime_col="pseudotime", n_bins=2)

        self.assertEqual(summary.shape[0], 2)
        self.assertLess(summary.loc[0, "cnv_supported_fraction"], summary.loc[1, "cnv_supported_fraction"])
        self.assertLess(summary.loc[0, "mean_hcc_malignant_module"], summary.loc[1, "mean_hcc_malignant_module"])


if __name__ == "__main__":
    unittest.main()
