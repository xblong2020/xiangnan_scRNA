import unittest

import numpy as np
import pandas as pd

from scripts.select_celloracle_tfs_module6_4 import (
    apply_hard_filters,
    bh_qvalues,
    compute_base_grn_compatibility,
    compute_cellrank_target_overlap,
    select_main_panel,
)


class Module64CellOracleTfSelectionLogicTest(unittest.TestCase):
    def test_bh_qvalues_are_monotonic_in_original_order(self):
        q = bh_qvalues(pd.Series([0.03, 0.001, 0.2], index=["a", "b", "c"]))

        self.assertAlmostEqual(q.loc["b"], 0.003)
        self.assertAlmostEqual(q.loc["a"], 0.045)
        self.assertAlmostEqual(q.loc["c"], 0.2)

    def test_base_grn_compatibility_counts_candidate_outgoing_links(self):
        base_grn = pd.DataFrame(
            {
                "peak_id": ["p1", "p2", "p3"],
                "gene_short_name": ["A", "B", "B"],
                "JUN": [1.0, 0.0, 1.0],
                "HNF4A": [0.0, 1.0, 0.0],
            }
        )

        compat = compute_base_grn_compatibility(base_grn, ["JUN", "TPI1"], min_links=2).set_index("tf")

        self.assertTrue(bool(compat.loc["JUN", "tf_in_base_grn"]))
        self.assertEqual(int(compat.loc["JUN", "base_grn_outgoing_links"]), 2)
        self.assertEqual(int(compat.loc["JUN", "base_grn_target_genes"]), 2)
        self.assertFalse(bool(compat.loc["TPI1", "tf_in_base_grn"]))

    def test_target_overlap_uses_fisher_enrichment_against_gene_universe(self):
        candidates = pd.DataFrame(
            {
                "tf": ["JUN", "HNF4A"],
                "target_genes": ["A;B;C", "D;E"],
            }
        )
        drivers = pd.DataFrame(
            {
                "gene": ["A", "B", "F", "JUN"],
                "lineage": ["cnv_supported_malignant"] * 4,
                "corr": [0.4, 0.3, 0.2, 0.1],
                "qval": [0.001, 0.002, 0.01, 0.05],
                "rank_positive_corr": [1, 2, 3, 4],
            }
        )
        universe = pd.Index(list("ABCDEFGHIJ"))

        overlap = compute_cellrank_target_overlap(candidates, drivers, universe, top_n_drivers=3).set_index("tf")

        self.assertEqual(int(overlap.loc["JUN", "cellrank_top_driver_overlap_n"]), 2)
        self.assertGreater(float(overlap.loc["JUN", "cellrank_target_overlap_oddsratio"]), 1.0)
        self.assertAlmostEqual(float(overlap.loc["JUN", "tf_self_cellrank_corr"]), 0.1)

    def test_hard_filters_exclude_noncanonical_tf_and_allow_maintenance_min_targets(self):
        candidates = pd.DataFrame(
            {
                "tf": ["TPI1", "HNF4A", "JUN"],
                "best_nes": [4.0, 3.1, 5.0],
                "n_targets": [12, 8, 20],
            }
        )
        expression = pd.DataFrame(
            {
                "tf": ["TPI1", "HNF4A", "JUN"],
                "tf_in_expression": [True, True, True],
                "detection_rate_main": [0.8, 0.2, 0.5],
                "detected_dataset_count": [5, 3, 4],
            }
        )
        compat = pd.DataFrame(
            {
                "tf": ["TPI1", "HNF4A", "JUN"],
                "tf_in_base_grn": [False, True, True],
                "base_grn_outgoing_links": [0, 20, 50],
            }
        )

        filtered = apply_hard_filters(
            candidates,
            tf_catalog={"TPI1", "HNF4A", "JUN"},
            expression_summary=expression,
            base_grn_compatibility=compat,
            maintenance_tfs={"HNF4A"},
            noncanonical_exclusions={"TPI1"},
            min_base_grn_links=10,
        ).set_index("tf")

        self.assertFalse(bool(filtered.loc["TPI1", "hard_filter_pass"]))
        self.assertIn("noncanonical_exclusion", filtered.loc["TPI1", "exclusion_reason"])
        self.assertTrue(bool(filtered.loc["HNF4A", "hard_filter_pass"]))
        self.assertTrue(bool(filtered.loc["JUN", "hard_filter_pass"]))

    def test_select_main_panel_prefers_expected_roles_and_keeps_size_limit(self):
        scored = pd.DataFrame(
            {
                "tf": ["JUN", "FOS", "HNF4A", "PPARA", "MAFF"],
                "role": [
                    "pro_cnv_candidate",
                    "pro_cnv_candidate",
                    "anti_cnv_hepatocyte_maintenance",
                    "anti_cnv_hepatocyte_maintenance",
                    "pro_cnv_candidate",
                ],
                "total_score": [90.0, 80.0, 85.0, 70.0, 65.0],
                "hard_filter_pass": [True, True, True, True, True],
            }
        )

        selected = select_main_panel(
            scored,
            preferred_pro=["JUN", "FOS"],
            preferred_maintenance=["HNF4A"],
            reserve_order=["MAFF", "PPARA"],
            min_panel_size=3,
            max_panel_size=4,
        ).set_index("tf")

        self.assertTrue(bool(selected.loc["JUN", "selected_for_main_panel"]))
        self.assertTrue(bool(selected.loc["HNF4A", "selected_for_main_panel"]))
        self.assertFalse(bool(selected.loc["PPARA", "selected_for_main_panel"]))
        self.assertEqual(int(selected["selected_for_main_panel"].sum()), 3)


if __name__ == "__main__":
    unittest.main()
