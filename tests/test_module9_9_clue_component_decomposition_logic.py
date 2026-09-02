import unittest

import pandas as pd

from scripts.clue_component_decomposition_module9_9 import (
    aggregate_clue_components,
    add_component_percentiles,
    build_evidence_adjusted_priority,
    enrich_compound_metadata,
)


class Module99ClueComponentDecompositionLogicTest(unittest.TestCase):
    def test_aggregate_uses_up_for_rescue_and_negative_down_for_malignant(self):
        signatures = pd.DataFrame(
            [
                {
                    "compound": "BRD-K1",
                    "cell_line": "HEPG2",
                    "cs_up": 0.4,
                    "cs_down": -0.2,
                    "cs_combined": 0.3,
                },
                {
                    "compound": "BRD-K1",
                    "cell_line": "HA1E",
                    "cs_up": 0.2,
                    "cs_down": -0.4,
                    "cs_combined": 0.3,
                },
            ]
        )

        row = aggregate_clue_components(signatures).iloc[0]

        self.assertAlmostEqual(float(row["clue_rescue_component"]), 0.3)
        self.assertAlmostEqual(float(row["clue_malignant_suppression_component"]), 0.3)
        self.assertAlmostEqual(float(row["clue_combined_component"]), 0.3)
        self.assertAlmostEqual(float(row["clue_hepg2_rescue_component"]), 0.4)
        self.assertAlmostEqual(
            float(row["clue_hepg2_malignant_suppression_component"]), 0.2
        )

    def test_component_percentiles_reward_balanced_branch_support(self):
        aggregate = pd.DataFrame(
            [
                {
                    "compound": "balanced",
                    "clue_rescue_component": 0.5,
                    "clue_malignant_suppression_component": 0.5,
                    "clue_combined_component": 0.5,
                },
                {
                    "compound": "one_sided",
                    "clue_rescue_component": 0.8,
                    "clue_malignant_suppression_component": -0.2,
                    "clue_combined_component": 0.2,
                },
                {
                    "compound": "low",
                    "clue_rescue_component": -0.5,
                    "clue_malignant_suppression_component": -0.5,
                    "clue_combined_component": -0.5,
                },
            ]
        )

        result = add_component_percentiles(aggregate).set_index("compound")

        self.assertGreater(
            float(result.loc["balanced", "clue_branch_balance_percentile"]),
            float(result.loc["one_sided", "clue_branch_balance_percentile"]),
        )

    def test_evidence_adjusted_priority_favors_cross_platform_balance(self):
        drugreflector = pd.DataFrame(
            [
                {
                    "compound": "cross_supported",
                    "decomposition_score": 0.95,
                    "both_biological_branches_top_200": True,
                },
                {
                    "compound": "model_only",
                    "decomposition_score": 0.99,
                    "both_biological_branches_top_200": True,
                },
            ]
        )
        clue = pd.DataFrame(
            [
                {
                    "compound": "cross_supported",
                    "clue_branch_balance_percentile": 0.95,
                    "clue_combined_percentile": 0.90,
                    "clue_liver_context_percentile": 0.85,
                }
            ]
        )

        priority = build_evidence_adjusted_priority(drugreflector, clue)

        self.assertEqual(priority.iloc[0]["compound"], "cross_supported")
        self.assertEqual(
            priority.set_index("compound").loc["cross_supported", "evidence_tier"],
            "A_cross_platform_balanced",
        )

    def test_evidence_adjusted_priority_fills_existing_missing_l1000_fields(self):
        drugreflector = pd.DataFrame(
            [
                {
                    "compound": "drug",
                    "decomposition_score": 0.9,
                    "both_biological_branches_top_200": True,
                    "l1000_support_score": pd.NA,
                }
            ]
        )
        clue = pd.DataFrame(
            [
                {
                    "compound": "drug",
                    "clue_branch_balance_percentile": 0.8,
                    "clue_combined_percentile": 0.8,
                    "clue_liver_context_percentile": 0.8,
                }
            ]
        )
        l1000 = pd.DataFrame(
            [
                {
                    "compound": "drug",
                    "l1000_similar_best_rank": 2,
                    "l1000_opposite_best_rank": pd.NA,
                    "l1000_support_score": 0.5,
                }
            ]
        )

        priority = build_evidence_adjusted_priority(drugreflector, clue, l1000)

        self.assertAlmostEqual(float(priority.iloc[0]["l1000_support_score"]), 0.5)
        self.assertNotIn("l1000_support_score_external", priority.columns)

    def test_enrich_compound_metadata_fills_missing_names_without_overwriting_existing(self):
        priority = pd.DataFrame(
            [
                {"compound": "BRD-K1", "pert_iname": "existing"},
                {"compound": "BRD-K2", "pert_iname": pd.NA},
            ]
        )
        metadata = pd.DataFrame(
            [
                {"compound": "BRD-K1", "pert_iname": "replacement", "pert_type": "trt_cp"},
                {"compound": "BRD-K2", "pert_iname": "filled", "pert_type": "trt_cp"},
            ]
        )

        enriched = enrich_compound_metadata(priority, metadata).set_index("compound")

        self.assertEqual(enriched.loc["BRD-K1", "pert_iname"], "existing")
        self.assertEqual(enriched.loc["BRD-K2", "pert_iname"], "filled")
        self.assertEqual(enriched.loc["BRD-K2", "pert_type"], "trt_cp")


if __name__ == "__main__":
    unittest.main()
