import unittest

import pandas as pd

from scripts.landmark_signature_decomposition_module9_9 import (
    build_decomposed_profiles,
    build_final_priority_table,
    select_landmark_signature,
)


class Module99LandmarkDecompositionLogicTest(unittest.TestCase):
    def setUp(self):
        self.signature = pd.DataFrame(
            [
                {
                    "gene": "SOX4",
                    "desired_direction": "down",
                    "component": "sox4_state_specific",
                    "final_weight": 1.0,
                    "include_primary": True,
                    "include_sensitivity": True,
                    "conflict_flag": False,
                    "housekeeping_or_qc_flag": False,
                },
                {
                    "gene": "JUN",
                    "desired_direction": "down",
                    "component": "ap1_stress_proliferation",
                    "final_weight": 0.8,
                    "include_primary": False,
                    "include_sensitivity": True,
                    "conflict_flag": False,
                    "housekeeping_or_qc_flag": False,
                },
                {
                    "gene": "HNF4A",
                    "desired_direction": "up",
                    "component": "hnf4a_ppara_rescue",
                    "final_weight": 1.2,
                    "include_primary": True,
                    "include_sensitivity": True,
                    "conflict_flag": False,
                    "housekeeping_or_qc_flag": False,
                },
                {
                    "gene": "ALB",
                    "desired_direction": "up",
                    "component": "mature_hepatocyte",
                    "final_weight": 0.9,
                    "include_primary": True,
                    "include_sensitivity": True,
                    "conflict_flag": False,
                    "housekeeping_or_qc_flag": False,
                },
            ]
        )

    def test_select_landmark_signature_keeps_only_auditable_landmark_genes(self):
        selected = select_landmark_signature(
            self.signature,
            landmark_genes={"SOX4", "JUN", "HNF4A"},
            sensitivity=True,
        )

        self.assertEqual(set(selected["gene"]), {"SOX4", "JUN", "HNF4A"})
        self.assertNotIn("ALB", selected["gene"].tolist())

    def test_combined_balanced_has_equal_up_and_down_absolute_weight(self):
        selected = select_landmark_signature(
            self.signature,
            landmark_genes={"SOX4", "JUN", "HNF4A"},
            sensitivity=True,
        )
        profiles = build_decomposed_profiles(selected)
        combined = profiles.loc[profiles["profile"].eq("combined_balanced")]

        up_mass = combined.loc[combined["v_score"].gt(0), "v_score"].abs().sum()
        down_mass = combined.loc[combined["v_score"].lt(0), "v_score"].abs().sum()

        self.assertAlmostEqual(float(up_mass), 1.0)
        self.assertAlmostEqual(float(down_mass), 1.0)
        self.assertTrue(
            profiles.loc[profiles["profile"].eq("malignant_only"), "v_score"].lt(0).all()
        )
        self.assertTrue(
            profiles.loc[profiles["profile"].eq("rescue_only"), "v_score"].gt(0).all()
        )

    def test_final_priority_rewards_support_in_both_biological_branches(self):
        predictions = pd.DataFrame(
            [
                {"compound": "balanced", "profile": "malignant_only", "rank_1based": 10},
                {"compound": "balanced", "profile": "rescue_only", "rank_1based": 12},
                {"compound": "balanced", "profile": "combined_balanced", "rank_1based": 8},
                {"compound": "one_sided", "profile": "malignant_only", "rank_1based": 1},
                {"compound": "one_sided", "profile": "rescue_only", "rank_1based": 500},
                {"compound": "one_sided", "profile": "combined_balanced", "rank_1based": 20},
            ]
        )

        priority = build_final_priority_table(predictions, n_compounds=1000)

        self.assertEqual(priority.iloc[0]["compound"], "balanced")
        self.assertGreater(
            float(priority.set_index("compound").loc["balanced", "branch_balance_score"]),
            float(priority.set_index("compound").loc["one_sided", "branch_balance_score"]),
        )
        self.assertEqual(
            int(priority.set_index("compound").loc["balanced", "n_profiles_top_200"]),
            3,
        )


if __name__ == "__main__":
    unittest.main()
