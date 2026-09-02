import unittest

import numpy as np
import pandas as pd
from scipy import sparse

from scripts.build_celloracle_object_module6_6 import (
    build_celloracle_state,
    compute_gene_metrics,
    select_celloracle_genes,
    validate_embedding_alignment,
)


class Module66CellOracleObjectLogicTest(unittest.TestCase):
    def test_build_celloracle_state_collapses_malignant_and_preserves_intermediate_states(self):
        obs = pd.DataFrame(
            {
                "trajectory_role": [
                    "normal_reference",
                    "malignant_cnv_supported",
                    "malignant_like_scanvi_review",
                    "regenerative_progenitor",
                    "stressed_injured",
                    "proliferating_candidate",
                    "unexpected_state",
                ],
            }
        )

        state = build_celloracle_state(obs)

        self.assertEqual(state.iloc[0], "normal_reference")
        self.assertEqual(state.iloc[1], "malignant_or_malignant_like")
        self.assertEqual(state.iloc[2], "malignant_or_malignant_like")
        self.assertEqual(state.iloc[3], "regenerative_progenitor")
        self.assertEqual(state.iloc[4], "stressed_injured")
        self.assertEqual(state.iloc[5], "proliferating_candidate")
        self.assertEqual(state.iloc[6], "other_trajectory")

    def test_compute_gene_metrics_reports_mean_detection_and_dispersion(self):
        matrix = sparse.csr_matrix(
            np.array(
                [
                    [1, 0, 4],
                    [1, 0, 0],
                    [1, 3, 0],
                ],
                dtype=float,
            )
        )

        metrics = compute_gene_metrics(matrix, pd.Index(["A", "B", "C"])).set_index("gene")

        self.assertAlmostEqual(metrics.loc["A", "detection_rate"], 1.0)
        self.assertAlmostEqual(metrics.loc["B", "mean_counts"], 1.0)
        self.assertGreater(metrics.loc["C", "dispersion_score"], metrics.loc["A", "dispersion_score"])

    def test_select_celloracle_genes_respects_cap_and_forces_input_tfs(self):
        metrics = pd.DataFrame(
            {
                "gene": ["A", "B", "C", "JUN", "HNF4A"],
                "detection_rate": [0.9, 0.9, 0.9, 0.01, 0.0],
                "mean_counts": [2, 2, 2, 0.1, 0],
                "dispersion_score": [10, 8, 6, 0, 0],
            }
        )

        selected = select_celloracle_genes(
            metrics,
            input_tfs=["JUN", "HNF4A", "MISSING"],
            max_genes=4,
            min_detection_rate=0.05,
        )

        self.assertEqual(len(selected), 4)
        self.assertIn("JUN", selected)
        self.assertIn("HNF4A", selected)
        self.assertNotIn("MISSING", selected)
        self.assertIn("A", selected)

    def test_validate_embedding_alignment_reorders_cells_and_rejects_missing(self):
        embedding = pd.DataFrame(
            {"UMAP1": [1.0, 2.0], "UMAP2": [3.0, 4.0]},
            index=["cell_b", "cell_a"],
        )

        aligned = validate_embedding_alignment(["cell_a", "cell_b"], embedding)

        self.assertEqual(aligned.shape, (2, 2))
        self.assertAlmostEqual(aligned[0, 0], 2.0)
        with self.assertRaises(ValueError):
            validate_embedding_alignment(["cell_a", "cell_c"], embedding)


if __name__ == "__main__":
    unittest.main()
