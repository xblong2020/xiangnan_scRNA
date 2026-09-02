import unittest

import numpy as np
import pandas as pd

from scripts.statistical_mediation_module9_3 import (
    assess_outcome_availability,
    benjamini_hochberg,
    bootstrap_indirect_effects,
    build_bulk_axis_scores,
    build_scrna_pseudobulk_scores,
    compute_path_coefficients,
    fit_cox_path_model,
)


class Module93MediationLogicTest(unittest.TestCase):
    def test_bulk_axis_scores_invert_tier1_rescue_to_hnf4a_ppara_loss(self):
        scores = pd.DataFrame(
            {
                "dataset_id": ["D1"] * 9,
                "sample": ["s1", "s2", "s3"] * 3,
                "axis": ["tier1_rescue"] * 3
                + ["ap1_stress_proliferation"] * 3
                + ["sox4_state_specific"] * 3,
                "signature_score": [3.0, 2.0, 1.0, 1.0, 2.0, 3.0, 0.5, 1.0, 2.0],
                "sample_type": ["tumor"] * 9,
            }
        )

        axis_scores = build_bulk_axis_scores(scores)

        indexed = axis_scores.set_index("sample")
        self.assertLess(indexed.loc["s1", "A_hnf4a_ppara_loss"], indexed.loc["s3", "A_hnf4a_ppara_loss"])
        self.assertGreater(indexed.loc["s3", "B_transition_activation"], indexed.loc["s1", "B_transition_activation"])
        self.assertIn("C_sox4_axis", axis_scores.columns)

    def test_synthetic_path_model_recovers_positive_sequential_mediation(self):
        rng = np.random.default_rng(10)
        n = 120
        a = rng.normal(size=n)
        b = 0.8 * a + rng.normal(scale=0.15, size=n)
        c = 0.7 * b + 0.2 * a + rng.normal(scale=0.15, size=n)
        y = 0.6 * c + 0.2 * b + rng.normal(scale=0.15, size=n)
        data = pd.DataFrame(
            {
                "dataset_id": "synthetic",
                "A_hnf4a_ppara_loss": a,
                "B_transition_activation": b,
                "C_sox4_axis": c,
                "outcome_value": y,
            }
        )

        coefficients = compute_path_coefficients(data, outcome_col="outcome_value", outcome_type="continuous")
        effects = bootstrap_indirect_effects(data, outcome_col="outcome_value", n_bootstrap=120, random_state=3)
        sequential = effects.set_index("indirect_path").loc["A_to_B_to_C_to_outcome"]

        self.assertGreater(float(coefficients.set_index("path").loc["A_to_B", "coef"]), 0.5)
        self.assertGreater(float(coefficients.set_index("path").loc["B_to_C", "coef"]), 0.4)
        self.assertGreater(float(sequential["effect"]), 0.15)
        self.assertGreater(float(sequential["ci_low"]), 0.0)

    def test_cox_wrapper_reports_not_testable_when_events_are_insufficient(self):
        low_event = pd.DataFrame(
            {
                "time": [1, 2, 3, 4, 5, 6],
                "event": [0, 0, 0, 0, 0, 1],
                "A_hnf4a_ppara_loss": np.linspace(-1, 1, 6),
                "B_transition_activation": np.linspace(-1, 1, 6),
                "C_sox4_axis": np.linspace(-1, 1, 6),
            }
        )

        result = fit_cox_path_model(low_event, time_col="time", event_col="event")

        self.assertEqual(result.loc[0, "status"], "not_testable_insufficient_events")

    def test_cox_wrapper_returns_hazard_ratio_for_testable_data(self):
        rng = np.random.default_rng(4)
        n = 50
        c = rng.normal(size=n)
        time = rng.exponential(scale=np.exp(-0.25 * c), size=n) + 0.1
        event = rng.binomial(1, 0.65, size=n)
        data = pd.DataFrame(
            {
                "time": time,
                "event": event,
                "A_hnf4a_ppara_loss": rng.normal(size=n),
                "B_transition_activation": rng.normal(size=n),
                "C_sox4_axis": c,
            }
        )

        result = fit_cox_path_model(data, time_col="time", event_col="event")

        self.assertEqual(result.loc[0, "status"], "tested")
        self.assertTrue(np.isfinite(float(result.loc[0, "hazard_ratio"])))

    def test_scrna_pseudobulk_aggregates_axis_scores_and_malignant_fraction(self):
        cells = pd.DataFrame(
            {
                "run_id": ["r1"] * 4,
                "method": ["monocle3"] * 4,
                "dataset": ["D"] * 4,
                "cnv_sample": ["s1", "s1", "s1", "s2"],
                "A_hnf4a_ppara_loss": [1.0, 2.0, 3.0, 10.0],
                "B_transition_activation": [2.0, 3.0, 4.0, 11.0],
                "C_sox4_axis": [3.0, 4.0, 5.0, 12.0],
                "C_malignant_like_fate": [0.0, 1.0, 1.0, 0.5],
                "cell_disease_stage": [
                    "stage_0_reference_hepatocyte",
                    "stage_4_cnv_supported_malignant",
                    "stage_4_malignant_like_review",
                    "stage_1_stressed_injured",
                ],
            }
        )
        stage_counts = pd.DataFrame(
            {
                "dataset": ["D", "D"],
                "cnv_sample": ["s1", "s2"],
                "sample_disease_stage": ["primary_hcc_tumor", "reference_adjacent_liver"],
            }
        )

        pseudo = build_scrna_pseudobulk_scores(cells, stage_counts)
        s1 = pseudo.set_index("cnv_sample").loc["s1"]

        self.assertAlmostEqual(float(s1["A_hnf4a_ppara_loss"]), 2.0)
        self.assertAlmostEqual(float(s1["malignant_fraction"]), 2 / 3)
        self.assertGreater(float(s1["sample_stage_ordinal"]), float(pseudo.set_index("cnv_sample").loc["s2", "sample_stage_ordinal"]))

    def test_outcome_availability_and_fdr(self):
        data = pd.DataFrame({"stage": [1, 2, 3, np.nan], "grade": [np.nan] * 4, "OS_event": [0, 1, 0, 1]})
        availability = assess_outcome_availability(data, ["stage", "grade", "OS_event"]).set_index("outcome")
        adjusted = benjamini_hochberg([0.01, 0.02, 0.5])

        self.assertEqual(availability.loc["stage", "status"], "available")
        self.assertEqual(availability.loc["grade", "status"], "not_testable_no_observed_values")
        self.assertLessEqual(adjusted[0], adjusted[1])


if __name__ == "__main__":
    unittest.main()
