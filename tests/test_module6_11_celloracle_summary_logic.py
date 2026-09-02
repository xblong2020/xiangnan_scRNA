import unittest

import pandas as pd

from scripts.summarize_celloracle_module6_11 import (
    build_candidate_tier_table,
    build_claim_evidence_matrix,
    build_figure_plan_rows,
    build_report_payload,
    build_supplementary_table_index,
)


class Module611CellOracleSummaryLogicTest(unittest.TestCase):
    def _inputs(self):
        quantitative = pd.DataFrame(
            {
                "tf": ["HNF4A", "SOX4", "JUN", "MAFB", "HLF"],
                "quantitative_rank": [1, 5, 6, 15, 10],
                "quantitative_perturbation_score": [0.9, 0.75, 0.7, 0.05, 0.49],
                "state_specificity_score": [0.4, 0.81, 0.55, -0.1, 0.2],
                "state_specificity_ratio": [1.4, 6.4, 1.8, 3.0, 1.5],
                "malignant_fate_direction_score": [0.8, 0.7, 0.5, -0.2, 0.1],
                "inner_product_score": [0.8, 0.7, 0.5, -0.2, 0.1],
                "cnv_fate_probability_association_score": [0.8, 0.7, 0.5, -0.2, 0.1],
                "module_rescue_score": [0.8, 0.7, 0.5, -0.2, 0.1],
            }
        )
        strict_union = pd.DataFrame(
            {
                "tf": ["HNF4A", "SOX4", "JUN", "MAFB", "HLF"],
                "driver_union_rank": [1, 5, 6, 15, 10],
                "driver_union_score": [0.9, 0.75, 0.7, 0.05, 0.49],
                "main_strict_rank": [1, 5, 7, 15, 10],
                "main_strict_score": [0.9, 0.74, 0.65, 0.05, 0.49],
                "rank_delta_main_minus_union": [0, 0, 1, 0, 0],
            }
        )
        phase = pd.DataFrame(
            {
                "tf": ["HNF4A", "SOX4", "JUN", "MAFB", "HLF"],
                "phase_early_rank": [15, 5, 1, 10, 12],
                "phase_intermediate_rank": [3, 4, 7, 15, 9],
                "phase_late_rank": [1, 6, 5, 15, 9],
            }
        )
        lodo = pd.DataFrame(
            {
                "tf": ["HNF4A", "SOX4", "JUN", "MAFB", "HLF"],
                "min_lodo_score": [0.8, 0.7, 0.4, 0.01, 0.2],
                "max_lodo_rank": [2, 6, 9, 15, 11],
                "top5_lodo_fraction": [1.0, 0.92, 0.5, 0.0, 0.0],
            }
        )
        controls = pd.DataFrame(
            {
                "tf": ["MAFB", "HLF"],
                "control_type": ["low_score_tf", "housekeeping_like_proxy_tf"],
                "control_rationale": ["lowest score", "panel-level proxy control"],
            }
        )
        evidence = pd.DataFrame(
            {
                "tf": ["HNF4A", "SOX4", "JUN", "MAFB", "HLF"],
                "integrated_rank": [9, 11, 1, 14, 15],
                "integrated_evidence_score": [0.55, 0.52, 0.77, 0.38, 0.26],
                "role": ["rescue", "state", "ap1", "control", "control"],
                "tier": ["candidate"] * 5,
                "selected_for_main_panel": [1, 1, 1, 0, 0],
            }
        )
        risks = pd.DataFrame(
            {
                "risk_type": ["lodo_instability", "control_proxy"],
                "tf": ["JUN", "NA"],
                "detail": ["weaker than Tier 1", "HLF is proxy only"],
                "severity": ["medium", "medium"],
            }
        )
        return quantitative, strict_union, phase, lodo, controls, evidence, risks

    def test_candidate_tier_assignment_is_stable(self):
        table = build_candidate_tier_table(*self._inputs()).set_index("tf")

        self.assertEqual(table.loc["HNF4A", "candidate_tier"], "Tier 1")
        self.assertEqual(table.loc["SOX4", "candidate_tier"], "Tier 2")
        self.assertEqual(table.loc["JUN", "candidate_tier"], "Tier 3")
        self.assertEqual(table.loc["MAFB", "candidate_tier"], "Negative/control")
        self.assertEqual(table.loc["HLF", "candidate_tier"], "Negative/control")
        self.assertIn("proxy control", table.loc["HLF", "primary_interpretation"])

    def test_claim_evidence_matrix_has_source_tables_and_metrics(self):
        table = build_candidate_tier_table(*self._inputs())
        claims = build_claim_evidence_matrix(table)

        self.assertGreaterEqual(len(claims), 5)
        self.assertFalse(claims["primary_source_table"].isna().any())
        self.assertFalse(claims["key_metrics"].isna().any())
        self.assertTrue((claims["key_metrics"].astype(str).str.len() > 0).all())

    def test_figure_plan_has_required_panel_fields(self):
        rows = build_figure_plan_rows()

        figure6 = [row for row in rows if row["figure"].startswith("Figure 6")]
        self.assertEqual(len(figure6), 5)
        for row in figure6:
            self.assertTrue(row["claim"])
            self.assertTrue(row["input_file"])
            self.assertTrue(row["metric"])
            self.assertTrue(row["output_or_figure_reference"])

    def test_supplementary_index_has_required_manuscript_use(self):
        index = build_supplementary_table_index()

        self.assertEqual(len(index), 5)
        self.assertFalse(index["primary_source_file"].isna().any())
        self.assertTrue((index["description"].astype(str).str.len() > 0).all())
        self.assertTrue((index["intended_manuscript_use"].astype(str).str.len() > 0).all())

    def test_report_integrity_contains_counts(self):
        table = build_candidate_tier_table(*self._inputs())
        figures = build_figure_plan_rows()
        supp = build_supplementary_table_index()
        report = build_report_payload(table, figures, supp, inputs={"a": "b"}, outputs={"c": "d"})

        self.assertEqual(report["n_candidates"], 5)
        self.assertEqual(report["figure_count"], len(figures))
        self.assertEqual(report["supplementary_table_count"], len(supp))
        self.assertIn("Tier 1", report["tier_counts"])
        self.assertIn("HLF", report["negative_control_tfs"])


if __name__ == "__main__":
    unittest.main()
