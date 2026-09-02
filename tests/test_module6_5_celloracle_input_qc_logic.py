import unittest

import pandas as pd

from scripts.qc_celloracle_inputs_module6_5 import (
    check_required_columns,
    classify_celloracle_environment,
    summarize_base_grn_for_tfs,
    summarize_tf_input,
)


class Module65CellOracleInputQcLogicTest(unittest.TestCase):
    def test_check_required_columns_reports_missing_columns(self):
        result = check_required_columns(
            available=["dataset", "sample_id", "cellrank_fate_prob_cnv_supported_malignant"],
            required=[
                "dataset",
                "sample_id",
                "driver_main_strict__pseudotime_mean",
            ],
        )

        self.assertFalse(result["all_present"])
        self.assertEqual(result["missing"], ["driver_main_strict__pseudotime_mean"])

    def test_summarize_base_grn_for_tfs_counts_links_and_targets(self):
        base_grn = pd.DataFrame(
            {
                "peak_id": ["p1", "p2", "p3"],
                "gene_short_name": ["A", "B", "B"],
                "JUN": [1.0, 0.0, 1.0],
                "HNF4A": [0.0, 1.0, 0.0],
            }
        )

        summary = summarize_base_grn_for_tfs(base_grn, ["JUN", "TPI1"]).set_index("tf")

        self.assertTrue(bool(summary.loc["JUN", "tf_in_base_grn"]))
        self.assertEqual(int(summary.loc["JUN", "base_grn_outgoing_links"]), 2)
        self.assertEqual(int(summary.loc["JUN", "base_grn_target_genes"]), 2)
        self.assertFalse(bool(summary.loc["TPI1", "tf_in_base_grn"]))

    def test_summarize_tf_input_merges_selection_expression_and_base_grn(self):
        selection = pd.DataFrame(
            {
                "tf": ["JUN", "HNF4A"],
                "role": ["pro_cnv_candidate", "anti_cnv_hepatocyte_maintenance"],
                "perturbation_mode": ["KO", "OE"],
                "hard_filter_pass": [True, True],
                "selected_for_main_panel": [True, True],
            }
        )
        base_grn_summary = pd.DataFrame(
            {
                "tf": ["JUN", "HNF4A"],
                "tf_in_base_grn": [True, True],
                "base_grn_outgoing_links": [20, 10],
                "base_grn_target_genes": [18, 9],
            }
        )

        result = summarize_tf_input(
            tfs=["JUN", "HNF4A", "FOS"],
            genes=pd.Index(["JUN", "FOS", "A"]),
            selection=selection,
            base_grn_summary=base_grn_summary,
        ).set_index("tf")

        self.assertTrue(bool(result.loc["JUN", "tf_in_expression"]))
        self.assertEqual(result.loc["HNF4A", "perturbation_mode"], "OE")
        self.assertFalse(bool(result.loc["FOS", "tf_in_selection_table"]))
        self.assertFalse(bool(result.loc["HNF4A", "tf_in_expression"]))

    def test_classify_celloracle_environment_requires_import_success(self):
        ok = classify_celloracle_environment(package_present=True, import_ok=True, import_error="")
        failed = classify_celloracle_environment(
            package_present=True,
            import_ok=False,
            import_error="ModuleNotFoundError: No module named 'pysam'",
        )

        self.assertEqual(ok, "usable")
        self.assertEqual(failed, "package_present_but_import_failed")


if __name__ == "__main__":
    unittest.main()
