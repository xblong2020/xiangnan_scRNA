import unittest

import anndata as ad
import numpy as np
import pandas as pd

from scripts.build_trajectory_object_module5_1 import assign_trajectory_role, sample_source_class
from scripts.build_trajectory_object_module5_1 import attach_obs


class Module51TrajectoryLogicTest(unittest.TestCase):
    def test_cnv_supported_calls_override_hepatocyte_state(self):
        row = {
            "malignant_hcc_call": "malignant_hcc_cnv_support",
            "hepatocyte_state_label": "regenerative_progenitor_like_hepatocyte",
            "scanvi_unified_final_strict_label": "Cholangiocyte_Progenitor",
        }
        self.assertEqual(assign_trajectory_role(row), "malignant_cnv_supported")

    def test_scanvi_malignant_review_is_preserved(self):
        row = {
            "malignant_hcc_call": "non_malignant_or_unresolved",
            "hepatocyte_state_label": "malignant_hepatocyte_candidate_needs_cnv",
            "scanvi_unified_final_strict_label": "malignant_like_hepatocyte_needs_review",
        }
        self.assertEqual(assign_trajectory_role(row), "malignant_like_scanvi_review")

    def test_non_cnv_hepatocyte_states_map_to_reference_roles(self):
        self.assertEqual(
            assign_trajectory_role(
                {
                    "malignant_hcc_call": "not_module3_candidate",
                    "hepatocyte_state_label": "normal_hepatocyte_like",
                    "scanvi_unified_final_strict_label": "normal_hepatocyte_like",
                }
            ),
            "normal_reference",
        )
        self.assertEqual(
            assign_trajectory_role(
                {
                    "malignant_hcc_call": "not_module3_candidate",
                    "hepatocyte_state_label": "proliferating_hepatocyte_candidate",
                    "scanvi_unified_final_strict_label": "Cholangiocyte_Progenitor",
                }
            ),
            "proliferating_candidate",
        )

    def test_sample_source_class_uses_dataset_specific_rules(self):
        self.assertEqual(sample_source_class("GSE149614", "HCC07T"), "tumor")
        self.assertEqual(sample_source_class("GSE149614", "HCC07N"), "normal_adjacent")
        self.assertEqual(sample_source_class("GSE202379", "SITTA2"), "non_hcc_liver")
        self.assertEqual(sample_source_class("GSE185477", "C41_NST"), "normal_adjacent")

    def test_attach_obs_accepts_missing_values_in_categorical_columns(self):
        adata = ad.AnnData(np.ones((2, 1)))
        adata.obs_names = ["cell_a", "cell_b"]
        cells = pd.DataFrame(
            {
                "cell_id": ["cell_a", "cell_b"],
                "trajectory_role": pd.Categorical(["normal_reference", None]),
            }
        )

        attach_obs(adata, cells)

        self.assertEqual(list(adata.obs["trajectory_role"].astype(str)), ["normal_reference", "Unknown"])


if __name__ == "__main__":
    unittest.main()
