import unittest

import numpy as np
import pandas as pd

from scripts.drug_reversal_signature_module9_4 import (
    build_primary_signature,
    compute_c_fate_correlations,
    is_housekeeping_or_qc_gene,
    records_from_malignant_tf_markers,
    resolve_signature_records,
)


class Module94DrugReversalSignatureLogicTest(unittest.TestCase):
    def test_conflict_resolution_prioritizes_hepatocyte_rescue_over_lower_weight_down(self):
        records = pd.DataFrame(
            [
                {
                    "gene": "ALB",
                    "desired_direction": "down",
                    "component": "ap1_stress_proliferation",
                    "source_file": "ap1.tsv",
                    "source_metric": "rank",
                    "source_rank": 1,
                    "evidence_weight": 0.8,
                },
                {
                    "gene": "ALB",
                    "desired_direction": "up",
                    "component": "mature_hepatocyte",
                    "source_file": "fixed",
                    "source_metric": "fixed_marker",
                    "source_rank": 1,
                    "evidence_weight": 0.9,
                },
            ]
        )

        resolved = resolve_signature_records(records).set_index("gene")

        self.assertEqual(resolved.loc["ALB", "desired_direction"], "up")
        self.assertFalse(bool(resolved.loc["ALB", "conflict_flag"]))
        self.assertTrue(bool(resolved.loc["ALB", "include_sensitivity"]))

    def test_unresolved_equal_priority_conflict_is_sensitivity_only(self):
        records = pd.DataFrame(
            [
                {
                    "gene": "AMBIG",
                    "desired_direction": "down",
                    "component": "ap1_stress_proliferation",
                    "source_file": "down.tsv",
                    "source_metric": "rank",
                    "source_rank": 1,
                    "evidence_weight": 0.8,
                },
                {
                    "gene": "AMBIG",
                    "desired_direction": "up",
                    "component": "tier1_rescue",
                    "source_file": "up.tsv",
                    "source_metric": "rank",
                    "source_rank": 1,
                    "evidence_weight": 0.7,
                },
            ]
        )

        resolved = resolve_signature_records(records).set_index("gene")

        self.assertTrue(bool(resolved.loc["AMBIG", "conflict_flag"]))
        self.assertFalse(bool(resolved.loc["AMBIG", "include_primary"]))
        self.assertTrue(bool(resolved.loc["AMBIG", "include_sensitivity"]))

    def test_ap1_and_sox4_tf_identities_remain_down_when_pathways_overlap_rescue(self):
        records = pd.DataFrame(
            [
                {
                    "gene": "JUN",
                    "desired_direction": "down",
                    "component": "ap1_stress_proliferation",
                    "source_file": "ap1.tsv",
                    "source_metric": "rank",
                    "source_rank": 1,
                    "evidence_weight": 0.8,
                },
                {
                    "gene": "JUN",
                    "desired_direction": "up",
                    "component": "hnf4a_ppara_rescue",
                    "source_file": "rescue.tsv",
                    "source_metric": "rank",
                    "source_rank": 1,
                    "evidence_weight": 1.0,
                },
                {
                    "gene": "SOX4",
                    "desired_direction": "down",
                    "component": "sox4_state_specific",
                    "source_file": "sox4.tsv",
                    "source_metric": "rank",
                    "source_rank": 1,
                    "evidence_weight": 1.0,
                },
                {
                    "gene": "SOX4",
                    "desired_direction": "up",
                    "component": "hnf4a_ppara_rescue",
                    "source_file": "rescue.tsv",
                    "source_metric": "rank",
                    "source_rank": 1,
                    "evidence_weight": 1.0,
                },
            ]
        )

        resolved = resolve_signature_records(records).set_index("gene")

        self.assertEqual(resolved.loc["JUN", "desired_direction"], "down")
        self.assertEqual(resolved.loc["SOX4", "desired_direction"], "down")
        self.assertFalse(bool(resolved.loc["JUN", "conflict_flag"]))
        self.assertFalse(bool(resolved.loc["SOX4", "conflict_flag"]))

    def test_primary_signature_limits_each_direction_and_sorts_by_weight_then_rank(self):
        records = pd.DataFrame(
            [
                {
                    "gene": f"UP{i}",
                    "desired_direction": "up",
                    "component": "hnf4a_ppara_rescue",
                    "source_file": "up.tsv",
                    "source_metric": "rank",
                    "source_rank": i,
                    "evidence_weight": 1.0 - (i * 0.01),
                }
                for i in range(5)
            ]
            + [
                {
                    "gene": f"DOWN{i}",
                    "desired_direction": "down",
                    "component": "c_malignant_like_fate",
                    "source_file": "down.tsv",
                    "source_metric": "spearman_rho",
                    "source_rank": i,
                    "evidence_weight": 1.0 - (i * 0.01),
                }
                for i in range(5)
            ]
        )
        resolved = resolve_signature_records(records)

        primary = build_primary_signature(resolved, max_genes_per_direction=3, min_genes_per_direction=2)

        self.assertEqual(primary["up"], ["UP0", "UP1", "UP2"])
        self.assertEqual(primary["down"], ["DOWN0", "DOWN1", "DOWN2"])

    def test_primary_signature_keeps_fixed_tf_anchor_before_generic_records(self):
        records = pd.DataFrame(
            [
                {
                    "gene": "ATF3",
                    "desired_direction": "down",
                    "component": "ap1_stress_proliferation",
                    "source_file": "fixed_malignant_tf_markers",
                    "source_metric": "fixed_tf_marker",
                    "source_rank": 4,
                    "evidence_weight": 0.8,
                },
            ]
            + [
                {
                    "gene": f"GENERIC{i}",
                    "desired_direction": "down",
                    "component": "c_malignant_like_fate",
                    "source_file": "expr",
                    "source_metric": "spearman_rho",
                    "source_rank": i,
                    "evidence_weight": 1.0,
                }
                for i in range(5)
            ]
        )
        resolved = resolve_signature_records(records)

        primary = build_primary_signature(resolved, max_genes_per_direction=3, min_genes_per_direction=1)

        self.assertIn("ATF3", primary["down"])

    def test_housekeeping_qc_genes_are_excluded_from_primary(self):
        self.assertTrue(is_housekeeping_or_qc_gene("MT-CO1"))
        self.assertTrue(is_housekeeping_or_qc_gene("RPL13A"))
        self.assertTrue(is_housekeeping_or_qc_gene("RPS3"))
        self.assertTrue(is_housekeeping_or_qc_gene("MALAT1"))
        self.assertFalse(is_housekeeping_or_qc_gene("HNF4A"))

        records = pd.DataFrame(
            [
                {
                    "gene": "MT-CO1",
                    "desired_direction": "down",
                    "component": "c_malignant_like_fate",
                    "source_file": "expr",
                    "source_metric": "spearman_rho",
                    "source_rank": 1,
                    "evidence_weight": 1.0,
                }
            ]
        )

        resolved = resolve_signature_records(records).set_index("gene")

        self.assertTrue(bool(resolved.loc["MT-CO1", "housekeeping_or_qc_flag"]))
        self.assertFalse(bool(resolved.loc["MT-CO1", "include_primary"]))
        self.assertTrue(bool(resolved.loc["MT-CO1", "include_sensitivity"]))

    def test_fixed_malignant_tf_markers_anchor_sox4_and_ap1_down(self):
        records = pd.DataFrame(records_from_malignant_tf_markers())

        resolved = resolve_signature_records(records).set_index("gene")

        for gene in ["SOX4", "JUN", "FOS", "JUND", "ATF3"]:
            self.assertEqual(resolved.loc[gene, "desired_direction"], "down")
            self.assertTrue(bool(resolved.loc[gene, "include_primary"]))
        self.assertNotIn("CEBPB", resolved.index)
        self.assertNotIn("EGR1", resolved.index)

    def test_c_fate_correlation_keeps_positive_significant_top_genes(self):
        expression = pd.DataFrame(
            {
                "POS": [0, 1, 2, 3, 4, 5],
                "NEG": [5, 4, 3, 2, 1, 0],
                "FLAT": [1, 1, 1, 1, 1, 1],
                "WEAK": [0, 1, 0, 1, 0, 1],
            },
            index=[f"c{i}" for i in range(6)],
        )
        fate = pd.Series([0, 1, 2, 3, 4, 5], index=expression.index)

        result = compute_c_fate_correlations(expression, fate, rho_threshold=0.15, q_threshold=0.05, top_n=10)

        self.assertEqual(result["gene"].tolist(), ["POS"])
        self.assertGreater(float(result.loc[0, "spearman_rho"]), 0.9)
        self.assertLess(float(result.loc[0, "q_value"]), 0.05)


if __name__ == "__main__":
    unittest.main()
