from __future__ import annotations

import unittest

import pandas as pd

from scripts.prepare_figure1c_cytotrace2_inputs import original_cell_id, stratified_sample


class Figure1CInputPrepTests(unittest.TestCase):
    def test_original_cell_id_strips_study_sample_prefix(self) -> None:
        study_sample = "GSE149614__GSE149614_raw_counts"
        cell_id = "GSE149614__GSE149614_raw_counts__HCC01T_AAACCTGAGGGCATGT"
        self.assertEqual(original_cell_id(cell_id, study_sample), "HCC01T_AAACCTGAGGGCATGT")

    def test_stratified_sample_preserves_groups_and_target_size(self) -> None:
        df = pd.DataFrame(
            {
                "cell_id": [f"cell_{i}" for i in range(12)],
                "major_celltype": ["A"] * 6 + ["B"] * 6,
                "study_sample": ["S1"] * 3 + ["S2"] * 3 + ["S1"] * 3 + ["S2"] * 3,
            }
        )
        sampled = stratified_sample(
            df,
            group_cols=["major_celltype", "study_sample"],
            target_n=8,
            min_per_group=1,
            seed=1,
        )
        self.assertEqual(sampled.shape[0], 8)
        counts = sampled.groupby(["major_celltype", "study_sample"]).size().to_dict()
        self.assertGreaterEqual(counts[("A", "S1")], 1)
        self.assertGreaterEqual(counts[("A", "S2")], 1)
        self.assertGreaterEqual(counts[("B", "S1")], 1)
        self.assertGreaterEqual(counts[("B", "S2")], 1)


if __name__ == "__main__":
    unittest.main()
