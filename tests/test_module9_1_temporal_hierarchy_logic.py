import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import h5py
import numpy as np
import pandas as pd

from scripts.temporal_hierarchy_module9_1 import (
    audit_velocity_feasibility,
    build_axis_scores,
    compute_cellrank_sox4_association,
    grouped_bootstrap_order_tests,
    summarize_temporal_trend,
)


class Module91TemporalHierarchyLogicTest(unittest.TestCase):
    def test_axis_scores_invert_hnf4a_ppara_retention_to_loss(self):
        feature_values = pd.DataFrame(
            {
                "expr_HNF4A": [10.0, 5.0, 1.0],
                "expr_PPARA": [9.0, 5.0, 1.0],
                "regulon_HNF4A": [8.0, 5.0, 1.0],
                "regulon_PPARA": [7.0, 5.0, 1.0],
                "expr_JUN": [1.0, 5.0, 9.0],
                "expr_FOS": [1.0, 5.0, 9.0],
                "expr_CEBPB": [1.0, 5.0, 9.0],
                "expr_EGR1": [1.0, 5.0, 9.0],
                "expr_SOX4": [1.0, 3.0, 10.0],
                "sox4_target_signature": [1.0, 3.0, 10.0],
                "module_HCC_Malignant_Associated": [0.0, 2.0, 9.0],
                "module_Proliferation": [0.0, 2.0, 9.0],
                "cellrank_fate_prob_cnv_supported_malignant": [0.0, 0.5, 1.0],
            },
            index=["c1", "c2", "c3"],
        )

        scores, availability = build_axis_scores(feature_values)

        self.assertGreater(scores.loc["c1", "A_hnf4a_ppara_retention"], scores.loc["c3", "A_hnf4a_ppara_retention"])
        self.assertLess(scores.loc["c1", "A_hnf4a_ppara_loss"], scores.loc["c3", "A_hnf4a_ppara_loss"])
        self.assertTrue(availability.set_index("axis").loc["A_hnf4a_ppara_loss", "n_available_features"] >= 4)

    def test_onset_summary_recovers_ordered_synthetic_axis_timing(self):
        pseudotime = pd.Series(np.linspace(0.0, 1.0, 120), name="pseudotime")
        early = pd.Series((pseudotime > 0.18).astype(float), name="A_hnf4a_ppara_loss")
        middle = pd.Series((pseudotime > 0.42).astype(float), name="B_transition_activation")
        late = pd.Series((pseudotime > 0.68).astype(float), name="C_sox4_axis")

        a = summarize_temporal_trend(early, pseudotime, n_bins=20, feature="A_hnf4a_ppara_loss")
        b = summarize_temporal_trend(middle, pseudotime, n_bins=20, feature="B_transition_activation")
        c = summarize_temporal_trend(late, pseudotime, n_bins=20, feature="C_sox4_axis")

        self.assertLess(a["onset_time"], b["onset_time"])
        self.assertLess(b["onset_time"], c["onset_time"])
        self.assertEqual(a["trend_status"], "tested")

    def test_grouped_bootstrap_reports_high_order_probability_for_ordered_samples(self):
        rows = []
        rng = np.random.default_rng(2)
        for sample_idx in range(8):
            for x in np.linspace(0.0, 1.0, 40):
                rows.append(
                    {
                        "cell_id": f"s{sample_idx}_{x:.3f}",
                        "cnv_sample": f"s{sample_idx}",
                        "pseudotime": x,
                        "A_hnf4a_ppara_loss": float(x > 0.15) + rng.normal(0, 0.02),
                        "B_transition_activation": float(x > 0.40) + rng.normal(0, 0.02),
                        "C_sox4_axis": float(x > 0.62) + rng.normal(0, 0.02),
                        "C_malignant_like_fate": float(x > 0.78) + rng.normal(0, 0.02),
                    }
                )
        data = pd.DataFrame(rows)

        tests = grouped_bootstrap_order_tests(
            data,
            pseudotime_col="pseudotime",
            group_col="cnv_sample",
            n_bootstrap=80,
            random_state=7,
            n_bins=16,
        ).set_index("comparison")

        self.assertGreater(float(tests.loc["A_loss_before_B_transition", "order_probability"]), 0.8)
        self.assertGreater(float(tests.loc["B_transition_before_C_sox4", "order_probability"]), 0.8)

    def test_cellrank_association_detects_high_sox4_in_high_fate_cells(self):
        fate = np.r_[np.linspace(0.0, 0.3, 30), np.linspace(0.7, 1.0, 30)]
        sox4 = np.r_[np.linspace(-1.0, -0.2, 30), np.linspace(0.5, 2.0, 30)]
        cells = pd.DataFrame(
            {
                "cell_id": [f"c{i}" for i in range(60)],
                "C_sox4_axis": sox4,
                "cellrank_fate_prob_cnv_supported_malignant": fate,
            }
        )

        result = compute_cellrank_sox4_association(cells)

        self.assertEqual(result.loc[0, "status"], "tested")
        self.assertGreater(float(result.loc[0, "spearman_rho"]), 0.8)
        self.assertGreater(float(result.loc[0, "high_fate_effect_size"]), 1.0)

    def test_velocity_audit_marks_h5ad_without_spliced_unspliced_as_not_testable(self):
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "counts_only.h5ad"
            with h5py.File(path, "w") as handle:
                handle.create_group("layers").create_dataset("counts", data=np.ones((2, 2)))

            audit = audit_velocity_feasibility(path)

        self.assertEqual(audit.loc[0, "velocity_status"], "not_testable_missing_spliced_unspliced_layers")
        self.assertFalse(bool(audit.loc[0, "has_spliced_unspliced"]))


if __name__ == "__main__":
    unittest.main()
