import unittest

import numpy as np
import pandas as pd

from scripts.integrate_cistarget_regulon_module6_3c import (
    align_auc_to_cells,
    correlate_auc_with_fate,
    summarize_pruned_regulons,
)


class Module63cCisTargetIntegrationLogicTest(unittest.TestCase):
    def test_align_auc_to_cells_uses_cell_id_column_and_order(self):
        auc = pd.DataFrame(
            {
                "CellID": ["cell_b", "cell_a"],
                "JUN(+)": [0.8, 0.2],
                "SRP9(+)": [0.4, 0.1],
            }
        )

        aligned = align_auc_to_cells(auc, pd.Index(["cell_a", "cell_b"]))

        self.assertEqual(aligned.index.tolist(), ["cell_a", "cell_b"])
        self.assertAlmostEqual(aligned.loc["cell_a", "JUN(+)"], 0.2)
        self.assertAlmostEqual(aligned.loc["cell_b", "SRP9(+)"], 0.4)

    def test_correlate_auc_with_fate_uses_only_non_null_fate_cells(self):
        auc = pd.DataFrame({"JUN(+)": [1.0, 2.0, 100.0, 3.0]}, index=["c1", "c2", "c3", "c4"])
        cells = pd.DataFrame(
            {
                "cellrank_fate_prob_cnv_supported_malignant": [0.1, 0.2, np.nan, 0.3],
                "driver_main_strict__pseudotime_median": [0.1, 0.2, np.nan, 0.4],
            },
            index=auc.index,
        )

        result = correlate_auc_with_fate(
            auc,
            cells,
            fate_key="cellrank_fate_prob_cnv_supported_malignant",
            time_key="driver_main_strict__pseudotime_median",
        )

        self.assertEqual(int(result.loc[0, "n_fate_cells"]), 3)
        self.assertAlmostEqual(result.loc[0, "cnv_fate_pearson_r"], 1.0)

    def test_summarize_pruned_regulons_extracts_target_counts_and_best_nes(self):
        motifs = pd.DataFrame(
            {
                ("Enrichment", "TF"): ["JUN", "JUN", "SRP9"],
                ("Enrichment", "NES"): [4.0, 5.0, 3.5],
                ("Enrichment", "TargetGenes"): [
                    "[('A', 0.5), ('B', 0.2)]",
                    "[('B', 0.4), ('C', 0.1)]",
                    "[('D', 0.9)]",
                ],
                ("Motif", "MotifID"): ["m1", "m2", "m3"],
                ("Motif", "Annotation"): ["gene is directly annotated", "motif similar", "gene is directly annotated"],
            }
        )

        summary = summarize_pruned_regulons(motifs)

        jun = summary.set_index("regulon").loc["JUN(+)"]
        self.assertEqual(int(jun["n_targets"]), 3)
        self.assertAlmostEqual(jun["best_nes"], 5.0)
        self.assertEqual(jun["top_motif_id"], "m2")

    def test_summarize_pruned_regulons_accepts_pyscenic_tf_motif_index(self):
        motifs = pd.DataFrame(
            {
                "NES": [4.0],
                "TargetGenes": ["[('A', 0.5), ('B', 0.2)]"],
                "Annotation": ["gene is directly annotated"],
            },
            index=pd.MultiIndex.from_tuples([("JUN", "m1")], names=["TF", "MotifID"]),
        )

        summary = summarize_pruned_regulons(motifs)

        self.assertEqual(summary.loc[0, "regulon"], "JUN(+)")
        self.assertEqual(summary.loc[0, "top_motif_id"], "m1")


if __name__ == "__main__":
    unittest.main()
