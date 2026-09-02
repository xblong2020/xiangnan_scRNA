import unittest

from scripts.define_trajectory_stage_root_end_module5_2 import (
    assign_cell_disease_stage,
    assign_root_end_role,
    assign_sample_disease_stage,
)


class Module52StageRootEndLogicTest(unittest.TestCase):
    def test_sample_stage_uses_source_class_progression(self):
        self.assertEqual(assign_sample_disease_stage({"sample_source_class": "non_hcc_liver"}), "reference_non_hcc_liver")
        self.assertEqual(assign_sample_disease_stage({"sample_source_class": "normal_adjacent"}), "reference_adjacent_liver")
        self.assertEqual(assign_sample_disease_stage({"sample_source_class": "cirrhotic_or_chronic_liver"}), "chronic_liver")
        self.assertEqual(assign_sample_disease_stage({"sample_source_class": "tumor"}), "primary_hcc_tumor")
        self.assertEqual(assign_sample_disease_stage({"sample_source_class": "pvtt_tumor"}), "pvtt_tumor")
        self.assertEqual(
            assign_sample_disease_stage({"sample_source_class": "metastatic_tumor_lymphnode"}),
            "metastatic_lymphnode_tumor",
        )

    def test_cell_stage_prioritizes_cnv_supported_malignancy(self):
        row = {
            "trajectory_role": "malignant_cnv_supported",
            "sample_source_class": "tumor",
            "malignant_hcc_call": "malignant_hcc_cnv_support",
        }
        self.assertEqual(assign_cell_disease_stage(row), "stage_4_cnv_supported_malignant")

    def test_cell_stage_maps_reference_and_intermediate_states(self):
        self.assertEqual(
            assign_cell_disease_stage({"trajectory_role": "normal_reference"}),
            "stage_0_reference_hepatocyte",
        )
        self.assertEqual(
            assign_cell_disease_stage({"trajectory_role": "regenerative_progenitor"}),
            "stage_2_regenerative_progenitor",
        )
        self.assertEqual(
            assign_cell_disease_stage({"trajectory_role": "proliferating_candidate"}),
            "stage_3_proliferating_candidate",
        )

    def test_root_end_role_marks_strict_root_and_malignant_endpoint(self):
        root = {
            "trajectory_role": "normal_reference",
            "sample_source_class": "normal_adjacent",
            "trajectory_include_main": True,
        }
        end = {
            "trajectory_role": "malignant_cnv_supported",
            "sample_source_class": "pvtt_tumor",
            "trajectory_include_cnv_strict": True,
        }
        review = {
            "trajectory_role": "malignant_like_scanvi_review",
            "sample_source_class": "tumor",
            "trajectory_include_main": True,
        }
        self.assertEqual(assign_root_end_role(root), "root_reference")
        self.assertEqual(assign_root_end_role(end), "end_malignant_cnv")
        self.assertEqual(assign_root_end_role(review), "end_malignant_review")


if __name__ == "__main__":
    unittest.main()
