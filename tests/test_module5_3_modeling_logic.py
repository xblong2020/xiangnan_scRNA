import unittest

import pandas as pd

from scripts.prepare_trajectory_modeling_module5_3 import (
    build_model_mask,
    select_start_end_clusters,
    stratified_model_sample,
)


class Module53ModelingLogicTest(unittest.TestCase):
    def test_model_mask_main_strict_excludes_review_endpoint(self):
        cells = pd.DataFrame(
            {
                "trajectory_include_cnv_strict": [True, True, False, False],
                "trajectory_include_main": [True, True, True, False],
                "trajectory_root_end_role": [
                    "root_reference",
                    "end_malignant_cnv",
                    "end_malignant_review",
                    "excluded_from_main_trajectory",
                ],
            }
        )

        self.assertEqual(build_model_mask(cells, "main_strict").tolist(), [True, True, False, False])
        self.assertEqual(build_model_mask(cells, "include_review").tolist(), [True, True, True, False])

    def test_select_start_end_clusters_uses_root_malignant_and_progenitor_enrichment(self):
        cells = pd.DataFrame(
            {
                "leiden_trajectory": ["0"] * 6 + ["1"] * 6 + ["2"] * 6,
                "trajectory_root_end_role": ["root_reference"] * 5
                + ["intermediate_trajectory"]
                + ["end_malignant_cnv"] * 4
                + ["intermediate_trajectory"] * 2
                + ["intermediate_trajectory"] * 6,
                "cell_disease_stage": ["stage_0_reference_hepatocyte"] * 6
                + ["stage_4_cnv_supported_malignant"] * 4
                + ["stage_1_stressed_injured"] * 2
                + ["stage_2_regenerative_progenitor"] * 5
                + ["stage_3_proliferating_candidate"],
            }
        )

        selected = select_start_end_clusters(cells, cluster_key="leiden_trajectory", min_cells=3)

        self.assertEqual(selected["start_cluster"], "0")
        self.assertEqual(selected["malignant_end_cluster"], "1")
        self.assertEqual(selected["progenitor_end_cluster"], "2")
        self.assertEqual(selected["end_clusters"], ["1", "2"])

    def test_stratified_sample_keeps_selected_root_and_end_cells(self):
        cells = pd.DataFrame(
            {
                "cell_id": [f"cell_{idx}" for idx in range(20)],
                "trajectory_root_end_role": ["root_reference"] * 10 + ["end_malignant_cnv"] * 10,
                "trajectory_root_cell_selected": [True] + [False] * 19,
                "trajectory_end_cell_selected": [False] * 19 + [True],
            }
        )

        sampled = stratified_model_sample(cells, max_cells=6, seed=7, strata_cols=["trajectory_root_end_role"])

        self.assertLessEqual(sampled.shape[0], 6)
        self.assertIn("cell_0", set(sampled["cell_id"]))
        self.assertIn("cell_19", set(sampled["cell_id"]))
        self.assertTrue({"root_reference", "end_malignant_cnv"}.issubset(set(sampled["trajectory_root_end_role"])))


if __name__ == "__main__":
    unittest.main()
