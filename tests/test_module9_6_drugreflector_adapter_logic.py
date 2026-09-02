import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from scripts.drugreflector_adapter_module9_6 import (
    audit_drugreflector_dependencies,
    build_vscore_frame,
)


class Module96DrugReflectorAdapterLogicTest(unittest.TestCase):
    def test_build_vscore_frame_assigns_signed_scores_from_direction(self):
        signature = pd.DataFrame(
            [
                {
                    "gene": "HNF4A",
                    "desired_direction": "up",
                    "component": "hnf4a_ppara_rescue",
                    "final_weight": 1.0,
                    "include_primary": True,
                    "include_sensitivity": True,
                    "conflict_flag": False,
                    "housekeeping_or_qc_flag": False,
                    "source_file": "a.tsv",
                    "source_metric": "metric_a",
                },
                {
                    "gene": "SOX4",
                    "desired_direction": "down",
                    "component": "sox4_state_specific",
                    "final_weight": 0.8,
                    "include_primary": True,
                    "include_sensitivity": True,
                    "conflict_flag": False,
                    "housekeeping_or_qc_flag": False,
                    "source_file": "b.tsv",
                    "source_metric": "metric_b",
                },
            ]
        )

        vscore = build_vscore_frame(signature, primary_only=False).set_index("gene")

        self.assertEqual(float(vscore.loc["HNF4A", "v_score"]), 1.0)
        self.assertEqual(float(vscore.loc["SOX4", "v_score"]), -0.8)

    def test_primary_only_excludes_conflict_and_housekeeping_genes(self):
        signature = pd.DataFrame(
            [
                {
                    "gene": "HNF4A",
                    "desired_direction": "up",
                    "component": "hnf4a_ppara_rescue",
                    "final_weight": 1.0,
                    "include_primary": True,
                    "include_sensitivity": True,
                    "conflict_flag": False,
                    "housekeeping_or_qc_flag": False,
                    "source_file": "a.tsv",
                    "source_metric": "metric_a",
                },
                {
                    "gene": "AMBIG",
                    "desired_direction": "down",
                    "component": "c_malignant_like_fate",
                    "final_weight": 0.9,
                    "include_primary": False,
                    "include_sensitivity": True,
                    "conflict_flag": True,
                    "housekeeping_or_qc_flag": False,
                    "source_file": "b.tsv",
                    "source_metric": "metric_b",
                },
                {
                    "gene": "MT-CO1",
                    "desired_direction": "down",
                    "component": "c_malignant_like_fate",
                    "final_weight": 0.7,
                    "include_primary": False,
                    "include_sensitivity": True,
                    "conflict_flag": False,
                    "housekeeping_or_qc_flag": True,
                    "source_file": "c.tsv",
                    "source_metric": "metric_c",
                },
            ]
        )

        primary = build_vscore_frame(signature, primary_only=True)

        self.assertEqual(primary["gene"].tolist(), ["HNF4A"])

    def test_sensitivity_keeps_non_primary_auditable_genes(self):
        signature = pd.DataFrame(
            [
                {
                    "gene": "HNF4A",
                    "desired_direction": "up",
                    "component": "hnf4a_ppara_rescue",
                    "final_weight": 1.0,
                    "include_primary": True,
                    "include_sensitivity": True,
                    "conflict_flag": False,
                    "housekeeping_or_qc_flag": False,
                    "source_file": "a.tsv",
                    "source_metric": "metric_a",
                },
                {
                    "gene": "AMBIG",
                    "desired_direction": "down",
                    "component": "c_malignant_like_fate",
                    "final_weight": 0.9,
                    "include_primary": False,
                    "include_sensitivity": True,
                    "conflict_flag": True,
                    "housekeeping_or_qc_flag": False,
                    "source_file": "b.tsv",
                    "source_metric": "metric_b",
                },
            ]
        )

        sensitivity = build_vscore_frame(signature, primary_only=False)

        self.assertEqual(set(sensitivity["gene"].tolist()), {"AMBIG", "HNF4A"})
        self.assertTrue(bool(sensitivity.set_index("gene").loc["AMBIG", "conflict_flag"]))

    def test_dependency_audit_marks_model_missing_when_packages_and_checkpoints_absent(self):
        with TemporaryDirectory() as tmpdir:
            audit = audit_drugreflector_dependencies(
                candidate_checkpoint_dirs=[Path(tmpdir) / "missing_ckpts"],
                package_presence={"torch": False, "drugreflector": False, "zenodo_get": False},
                import_results={
                    "torch": {"import_ok": False, "error": "missing"},
                    "drugreflector": {"import_ok": False, "error": "missing"},
                    "zenodo_get": {"import_ok": False, "error": "missing"},
                },
            )

        self.assertEqual(audit["status"], "adapter_ready_model_missing")
        self.assertFalse(audit["packages"]["torch"]["present"])
        self.assertFalse(audit["checkpoint_summary"]["all_required_present"])

    def test_dependency_audit_marks_runtime_blocked_when_installed_package_import_fails(self):
        with TemporaryDirectory() as tmpdir:
            audit = audit_drugreflector_dependencies(
                candidate_checkpoint_dirs=[Path(tmpdir) / "missing_ckpts"],
                package_presence={"torch": True, "drugreflector": True, "zenodo_get": True},
                import_results={
                    "torch": {"import_ok": False, "error": "dll load failed"},
                    "drugreflector": {"import_ok": False, "error": "missing module"},
                    "zenodo_get": {"import_ok": True, "error": ""},
                },
            )

        self.assertEqual(audit["status"], "adapter_ready_runtime_blocked")
        self.assertFalse(audit["packages"]["torch"]["import_ok"])
        self.assertIn("dll", audit["packages"]["torch"]["import_error"])

    def test_dependency_audit_rejects_partial_checkpoint_file(self):
        with TemporaryDirectory() as tmpdir:
            checkpoint_dir = Path(tmpdir)
            (checkpoint_dir / "model_fold_0.pt").write_bytes(b"partial checkpoint")
            audit = audit_drugreflector_dependencies(
                candidate_checkpoint_dirs=[checkpoint_dir],
                package_presence={"torch": True, "drugreflector": True, "zenodo_get": True},
                import_results={
                    "torch": {"import_ok": True, "error": ""},
                    "drugreflector": {"import_ok": True, "error": ""},
                    "zenodo_get": {"import_ok": True, "error": ""},
                },
            )

        expected = {
            row["file_name"]: row
            for row in audit["checkpoint_summary"]["expected_checkpoints"]
        }
        self.assertTrue(expected["model_fold_0.pt"]["present"])
        self.assertFalse(expected["model_fold_0.pt"]["size_ok"])
        self.assertFalse(expected["model_fold_0.pt"]["valid"])
        self.assertEqual(audit["checkpoint_summary"]["n_valid_expected_checkpoints"], 0)
        self.assertFalse(audit["checkpoint_summary"]["all_required_present"])


if __name__ == "__main__":
    unittest.main()
