import unittest

import pandas as pd

from scripts.drugreflector_metadata_crossvalidation_module9_8 import (
    aggregate_l1000fwd_candidates,
    build_clue_crosswalk,
    build_crossvalidation_table,
    collapse_pert_metadata,
    compute_entity_overlap,
)


class Module98DrugReflectorMetadataCrossvalidationLogicTest(unittest.TestCase):
    def test_collapse_pert_metadata_prefers_higher_priority_nonmissing_values(self):
        phase1 = pd.DataFrame(
            [
                {
                    "pert_id": "BRD-K1",
                    "pert_iname": "old_name",
                    "pert_type": "trt_cp",
                    "is_touchstone": "1",
                    "inchi_key_prefix": pd.NA,
                    "inchi_key": pd.NA,
                    "canonical_smiles": "CCC",
                    "pubchem_cid": "1",
                    "source_dataset": "phase1",
                    "source_priority": 1,
                }
            ]
        )
        phase2 = phase1.copy()
        phase2["pert_iname"] = "new_name"
        phase2["canonical_smiles"] = pd.NA
        phase2["source_dataset"] = "phase2"
        phase2["source_priority"] = 2

        collapsed = collapse_pert_metadata([phase1, phase2]).iloc[0]

        self.assertEqual(collapsed["pert_iname"], "new_name")
        self.assertEqual(collapsed["canonical_smiles"], "CCC")
        self.assertTrue(bool(collapsed["metadata_conflict_flag"]))
        self.assertIn("pert_iname", collapsed["metadata_conflict_columns"])

    def test_aggregate_l1000fwd_candidates_separates_similar_and_opposite(self):
        candidates = pd.DataFrame(
            [
                {
                    "compound_id": "BRD-K1",
                    "result_group": "similar",
                    "rank_within_group": 2,
                    "raw_score": 0.8,
                    "cell_line": "A549",
                    "sig_id": "S1",
                },
                {
                    "compound_id": "BRD-K1",
                    "result_group": "opposite",
                    "rank_within_group": 4,
                    "raw_score": -0.7,
                    "cell_line": "PC3",
                    "sig_id": "S2",
                },
            ]
        )

        row = aggregate_l1000fwd_candidates(candidates).iloc[0]

        self.assertEqual(int(row["l1000_similar_signature_count"]), 1)
        self.assertEqual(int(row["l1000_opposite_signature_count"]), 1)
        self.assertAlmostEqual(float(row["l1000_support_score"]), 0.25)

    def test_build_crossvalidation_table_labels_supported_and_discordant_hits(self):
        consensus = pd.DataFrame(
            [
                {
                    "compound": "BRD-K1",
                    "primary_rank_1based": 1,
                    "sensitivity_rank_1based": 2,
                    "in_both_top_lists": True,
                    "mean_rank_1based": 1.5,
                },
                {
                    "compound": "BRD-K2",
                    "primary_rank_1based": 3,
                    "sensitivity_rank_1based": 4,
                    "in_both_top_lists": True,
                    "mean_rank_1based": 3.5,
                },
            ]
        )
        metadata = pd.DataFrame(
            [
                {"compound": "BRD-K1", "pert_iname": "drug1", "metadata_conflict_flag": False},
                {"compound": "BRD-K2", "pert_iname": "drug2", "metadata_conflict_flag": False},
            ]
        )
        l1000 = pd.DataFrame(
            [
                {
                    "compound": "BRD-K1",
                    "l1000_similar_signature_count": 1,
                    "l1000_opposite_signature_count": 0,
                    "l1000_similar_best_rank": 2,
                    "l1000_opposite_best_rank": pd.NA,
                    "l1000_support_score": 0.5,
                },
                {
                    "compound": "BRD-K2",
                    "l1000_similar_signature_count": 0,
                    "l1000_opposite_signature_count": 1,
                    "l1000_similar_best_rank": pd.NA,
                    "l1000_opposite_best_rank": 1,
                    "l1000_support_score": -1.0,
                },
            ]
        )

        result = build_crossvalidation_table(consensus, metadata, l1000, "skipped_missing_api_key")
        by_id = result.set_index("compound")

        self.assertEqual(by_id.loc["BRD-K1", "crossvalidation_status"], "supported_by_l1000fwd_similar")
        self.assertEqual(by_id.loc["BRD-K2", "crossvalidation_status"], "discordant_l1000fwd_opposite")
        self.assertEqual(by_id.loc["BRD-K1", "clue_query_status"], "skipped_missing_api_key")

    def test_compute_entity_overlap_counts_id_name_and_structure_layers(self):
        drugreflector = pd.DataFrame(
            [
                {"compound": "BRD-K1", "pert_iname": "Drug A", "inchi_key": "AAAA"},
                {"compound": "BRD-K2", "pert_iname": "Drug-B", "inchi_key": "BBBB"},
            ]
        )
        l1000 = pd.DataFrame(
            [
                {"compound": "BRD-K1", "pert_iname": "other", "inchi_key": "CCCC"},
                {"compound": "BRD-K9", "pert_iname": "drug b", "inchi_key": "AAAA"},
            ]
        )

        overlap = compute_entity_overlap(drugreflector, l1000)

        self.assertEqual(overlap["n_exact_brd_id_overlap"], 1)
        self.assertEqual(overlap["n_normalized_name_overlap"], 1)
        self.assertEqual(overlap["n_inchi_key_overlap"], 1)

    def test_build_clue_crosswalk_uses_name_alias_and_cell_scores(self):
        consensus_metadata = pd.DataFrame(
            [
                {
                    "compound": "BRD-OLD",
                    "pert_iname": "Drug A",
                    "inchi_key": pd.NA,
                }
            ]
        )
        all_metadata = pd.DataFrame(
            [
                {
                    "compound": "BRD-NEW",
                    "pert_iname": "drug-a",
                    "inchi_key": pd.NA,
                }
            ]
        )
        clue_summary = pd.DataFrame(
            [{"clue_compound": "BRD-NEW", "clue_tau": 92.0}]
        )
        clue_cell = pd.DataFrame(
            [
                {
                    "clue_compound": "BRD-NEW",
                    "cell_line": "HEPG2",
                    "clue_tau": 88.0,
                }
            ]
        )

        crosswalk, _ = build_clue_crosswalk(
            consensus_metadata,
            clue_summary,
            clue_cell,
            all_metadata,
        )
        row = crosswalk.iloc[0]

        self.assertEqual(row["clue_match_type"], "normalized_name")
        self.assertEqual(row["clue_matched_ids"], "BRD-NEW")
        self.assertEqual(float(row["clue_tau"]), 92.0)
        self.assertEqual(float(row["clue_hepg2_tau"]), 88.0)
        self.assertTrue(bool(row["clue_strong_support"]))


if __name__ == "__main__":
    unittest.main()
