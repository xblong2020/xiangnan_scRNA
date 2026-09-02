import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from scripts.network_direction_module9_2 import (
    AXIS_TF_GROUPS,
    audit_restore_availability,
    build_signature_sets,
    compute_asymmetry_tests,
    compute_celloracle_axis_impact_matrix,
    compute_celloracle_tf_axis_delta,
    compute_sctenifold_signature_impact,
)


class Module92NetworkDirectionLogicTest(unittest.TestCase):
    def test_signature_impact_deduplicates_signature_genes_and_counts_fdr_hits(self):
        signatures = pd.DataFrame(
            {
                "axis": ["tier1_rescue", "tier1_rescue", "tier1_rescue", "sox4_state_specific"],
                "tf": ["HNF4A", "PPARA", "HNF4A", "SOX4"],
                "gene": ["G1", "G1", "G2", "S1"],
                "signature_class": ["rescue", "rescue", "rescue", "state"],
            }
        )
        perturbation_genes = pd.DataFrame(
            {
                "tf": ["HNF4A", "HNF4A", "HNF4A", "SOX4"],
                "gene": ["G1", "G2", "NOISE", "S1"],
                "distance": [0.8, 0.4, 0.1, 1.0],
                "p.adj": [0.01, 0.20, 0.01, 0.001],
                "subset": ["driver_union_all"] * 4,
            }
        )

        signature_sets = build_signature_sets(signatures)
        impact = compute_sctenifold_signature_impact(perturbation_genes, signature_sets, fdr_threshold=0.05)
        hnf4a_rescue = impact.set_index(["perturb_tf", "target_axis", "subset"]).loc[
            ("HNF4A", "tier1_rescue", "driver_union_all")
        ]

        self.assertEqual(int(hnf4a_rescue["n_signature_genes"]), 2)
        self.assertEqual(int(hnf4a_rescue["n_overlap"]), 2)
        self.assertEqual(int(hnf4a_rescue["n_sig_fdr05"]), 1)
        self.assertAlmostEqual(float(hnf4a_rescue["sig_fraction"]), 0.5)
        self.assertGreater(float(hnf4a_rescue["impact_score"]), 0.0)

    def test_celloracle_axis_delta_averages_signed_target_tf_deltas(self):
        celloracle = pd.DataFrame(
            {
                "tf": ["HNF4A", "HNF4A", "HNF4A", "SOX4"],
                "target_tf": ["JUN", "FOS", "PPARA", "HNF4A"],
                "celloracle_state": ["malignant_or_malignant_like"] * 4,
                "mean_delta_x": [1.0, 3.0, -0.5, -0.25],
                "mean_abs_delta_x": [1.0, 3.0, 0.5, 0.25],
            }
        )

        axis_delta = compute_celloracle_tf_axis_delta(celloracle, AXIS_TF_GROUPS)
        hnf4a_to_b = axis_delta.set_index(["perturb_tf", "target_axis", "celloracle_state"]).loc[
            ("HNF4A", "B_transition", "malignant_or_malignant_like")
        ]

        self.assertEqual(int(hnf4a_to_b["n_target_tfs_observed"]), 2)
        self.assertAlmostEqual(float(hnf4a_to_b["mean_delta_x"]), 2.0)
        self.assertAlmostEqual(float(hnf4a_to_b["mean_abs_delta_x"]), 2.0)

    def test_celloracle_axis_impact_matrix_uses_abs_delta_as_impact_score(self):
        axis_delta = pd.DataFrame(
            {
                "perturb_tf": ["HNF4A", "SOX4"],
                "source_axis": ["A_upstream", "C_sox4"],
                "target_axis": ["B_transition", "A_upstream"],
                "celloracle_state": ["malignant_or_malignant_like", "malignant_or_malignant_like"],
                "mean_delta_x": [-2.0, 0.25],
                "mean_abs_delta_x": [2.0, 0.25],
            }
        )

        impact = compute_celloracle_axis_impact_matrix(axis_delta)
        hnf4a_to_b = impact.set_index(["perturb_tf", "target_axis", "subset"]).loc[
            ("HNF4A", "B_transition", "celloracle_malignant_or_malignant_like")
        ]

        self.assertEqual(hnf4a_to_b["evidence_source"], "celloracle_abs_delta")
        self.assertAlmostEqual(float(hnf4a_to_b["impact_score"]), 2.0)
        self.assertAlmostEqual(float(hnf4a_to_b["signed_delta"]), -2.0)

    def test_asymmetry_tests_detect_forward_greater_than_reverse(self):
        impact = pd.DataFrame(
            {
                "perturb_tf": ["HNF4A", "PPARA", "SOX4", "SOX4", "HLF"],
                "source_axis": ["A_upstream", "A_upstream", "C_sox4", "C_sox4", "control"],
                "target_axis": ["B_transition", "C_sox4", "A_upstream", "B_transition", "B_transition"],
                "subset": ["driver_union_all"] * 5,
                "impact_score": [2.0, 2.2, 0.2, 0.3, 0.1],
                "sig_fraction": [0.8, 0.9, 0.1, 0.1, 0.05],
            }
        )

        tests = compute_asymmetry_tests(impact)
        a_to_b = tests.loc[
            tests["comparison"].eq("A_to_B_vs_C_to_A") & tests["subset"].eq("driver_union_all")
        ].iloc[0]

        self.assertEqual(a_to_b["support_label"], "forward_greater_than_reverse")
        self.assertGreater(float(a_to_b["directionality_index"]), 0.7)
        self.assertGreater(float(a_to_b["forward_vs_control_ratio"]), 1.0)

    def test_sox4_reverse_test_marks_weak_reverse_when_upstream_axes_are_low(self):
        impact = pd.DataFrame(
            {
                "perturb_tf": ["SOX4", "SOX4", "SOX4"],
                "source_axis": ["C_sox4", "C_sox4", "C_sox4"],
                "target_axis": ["C_sox4", "A_upstream", "B_transition"],
                "subset": ["driver_union_all"] * 3,
                "impact_score": [2.0, 0.05, 0.1],
                "sig_fraction": [0.9, 0.02, 0.05],
            }
        )

        tests = compute_asymmetry_tests(impact)
        sox4_reverse = tests.loc[
            tests["comparison"].eq("SOX4_self_vs_reverse_upstream") & tests["subset"].eq("driver_union_all")
        ].iloc[0]

        self.assertEqual(sox4_reverse["support_label"], "weak_reverse_upstream")
        self.assertGreater(float(sox4_reverse["directionality_index"]), 0.8)

    def test_restore_audit_marks_missing_restore_outputs_as_not_available(self):
        with TemporaryDirectory() as tmpdir:
            audit = audit_restore_availability(Path(tmpdir) / "missing_restore_dir")

        self.assertEqual(audit.loc[0, "restore_status"], "not_available_existing_outputs")
        self.assertFalse(bool(audit.loc[0, "restore_available"]))


if __name__ == "__main__":
    unittest.main()
