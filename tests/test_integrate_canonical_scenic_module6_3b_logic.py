import unittest

import numpy as np
import pandas as pd

from scripts.integrate_canonical_scenic_module6_3b import (
    bh_qvalues,
    cnv_positive_mask,
    evaluate_axes,
    parse_targets,
    safe_corr,
)


class CanonicalScenicModule63bLogicTest(unittest.TestCase):
    def test_bh_qvalues_preserve_index_and_adjust(self):
        q = bh_qvalues(pd.Series([0.03, 0.001, 0.2], index=["a", "b", "c"]))

        self.assertAlmostEqual(q.loc["b"], 0.003)
        self.assertAlmostEqual(q.loc["a"], 0.045)
        self.assertAlmostEqual(q.loc["c"], 0.2)

    def test_parse_targets_extracts_gene_symbols_from_ctx_literal(self):
        targets = parse_targets("[('JUN', 0.9), ('FOS', 0.8)]")

        self.assertEqual(targets, {"JUN", "FOS"})

    def test_safe_corr_reports_sample_count_and_finite_statistics(self):
        rho, pvalue, n = safe_corr(
            pd.Series([1.0, 2.0, np.nan, 4.0]),
            pd.Series([1.0, 3.0, 4.0, 5.0]),
            "spearman",
        )

        self.assertEqual(n, 3)
        self.assertTrue(np.isfinite(rho))
        self.assertTrue(np.isfinite(pvalue))

    def test_cnv_positive_mask_accepts_boolean_and_label_columns(self):
        cells = pd.DataFrame(
            {
                "cnv_bool": [True, False, None],
                "cnv_label": ["cnv_supported_malignant", "reference", "1"],
            }
        )

        self.assertEqual(cnv_positive_mask(cells, "cnv_bool").tolist(), [True, False, False])
        self.assertEqual(cnv_positive_mask(cells, "cnv_label").tolist(), [True, False, True])

    def test_evaluate_axes_marks_missing_canonical_regulons(self):
        associations = pd.DataFrame(
            {
                "TF": ["HNF4A", "JUN"],
                "regulon": ["HNF4A(+)", "JUN(+)"],
                "regulon_size": [20, 30],
                "motif_NES_max": [4.0, 5.0],
                "spearman_rho": [-0.2, -0.1],
                "spearman_FDR": [0.01, 0.01],
                "pseudotime_spearman_rho": [-0.2, 0.4],
                "pseudotime_FDR": [0.01, 0.01],
                "sample_loo_stable": [True, True],
                "dataset_loo_stable": [True, True],
                "evidence_tier": ["Tier A", "Tier A"],
            }
        )

        axes = evaluate_axes(associations)

        self.assertFalse(bool(axes.loc[(axes["axis"] == "HNF4A_PPARA_identity") & axes["TF"].eq("PPARA"), "canonical_regulon_detected"].iloc[0]))
        self.assertEqual(axes.loc[axes["TF"].eq("PPARA"), "interpretation"].iloc[0], "not_detected_by_canonical_ctx")


if __name__ == "__main__":
    unittest.main()
