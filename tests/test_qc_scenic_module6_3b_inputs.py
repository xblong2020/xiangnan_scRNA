import unittest

import numpy as np
import pandas as pd
from scipy import sparse

from scripts.qc_scenic_module6_3b_inputs import (
    classify_gene_space,
    make_metric,
    summarize_metadata_columns,
)
from scripts.prepare_driver_union_full_expression_module6_3b import filter_genes


class ScenicModule63bInputQcTest(unittest.TestCase):
    def test_make_metric_marks_expected_value_pass(self):
        row = make_metric("n_cells", 9512, 9512, "loom cell count")

        self.assertEqual(row["status"], "PASS")
        self.assertEqual(row["value"], 9512)

    def test_make_metric_marks_failed_numeric_gate(self):
        row = make_metric("n_genes", 2000, 12000, "full-expression gene count")

        self.assertEqual(row["status"], "FAIL")

    def test_classify_gene_space_distinguishes_full_expression_from_hvg(self):
        self.assertEqual(classify_gene_space(12000, 2000), "full_expression")
        self.assertEqual(classify_gene_space(2000, 2000), "hvg_like")

    def test_summarize_metadata_columns_reports_missing_columns(self):
        result = summarize_metadata_columns(
            ["dataset", "sample_id", "driver_main_strict__pseudotime_phase"],
            ["dataset", "sample_id", "driver_main_strict__eligible"],
        )

        self.assertEqual(result["present"], 2)
        self.assertEqual(result["missing"], ["driver_main_strict__pseudotime_phase"])

    def test_filter_genes_does_not_keep_zero_expression_tf(self):
        matrix = sparse.csr_matrix(np.array([[1, 0, 0], [1, 0, 0]], dtype=float))

        kept = filter_genes(
            matrix,
            pd.Index(["JUN", "FOS", "GENE1"]),
            tf_list={"JUN", "FOS"},
            min_cells=1,
            min_mean=0.1,
            max_genes=10,
        )

        self.assertIn("JUN", kept)
        self.assertNotIn("FOS", kept)


if __name__ == "__main__":
    unittest.main()
