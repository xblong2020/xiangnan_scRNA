import unittest

import pandas as pd

from scripts.export_sctenifoldknk_input_module7_1 import (
    build_cell_subset_mask,
    build_export_report,
    validate_tf_coverage,
)
from scripts.integrate_sctenifoldknk_celloracle_module7_3 import (
    build_biological_axis_summary,
    build_concordance_summary,
    build_integrated_evidence_matrix,
    build_marker_overlap_summary,
    build_state_specific_gene_table,
    flag_control_outperformance,
)
from scripts.summarize_sctenifoldknk_enrichment_module7_4 import (
    build_enrichment_summary,
    compute_preranked_metric,
    compute_simple_gsea,
    filter_gene_sets_to_background,
    parse_gmt_lines,
    run_ora_for_tf,
    summarize_mapping_stats,
)
from scripts.summarize_sctenifoldknk_module7_5 import (
    assign_module7_display_group,
    build_module7_report_payload,
    build_pathway_level_enrichment_matrix,
    build_tf_level_replication_matrix,
)


class Module7ScTenifoldKnkLogicTest(unittest.TestCase):
    def test_validate_tf_coverage_requires_all_celloracle_tfs(self):
        genes = ["HNF4A", "PPARA", "CEBPB", "EGR1", "SOX4", "JUN"]
        tfs = ["HNF4A", "PPARA", "CEBPB", "EGR1", "SOX4", "JUN"]

        summary = validate_tf_coverage(genes, tfs)

        self.assertEqual(summary["n_input_tfs"], 6)
        self.assertEqual(summary["n_retained_tfs"], 6)
        self.assertEqual(summary["missing_tfs"], [])

    def test_validate_tf_coverage_raises_for_missing_tf(self):
        with self.assertRaisesRegex(ValueError, "Missing required scTenifoldKnk knockout TFs"):
            validate_tf_coverage(["HNF4A", "PPARA"], ["HNF4A", "SOX4"])

    def test_build_cell_subset_mask_supports_main_subsets(self):
        obs = pd.DataFrame(
            {
                "celloracle_state": ["normal_reference", "malignant_or_malignant_like", "malignant_or_malignant_like"],
                "celloracle_main_strict": [True, False, True],
            },
            index=["c1", "c2", "c3"],
        )

        self.assertEqual(build_cell_subset_mask(obs, "driver_union_all").tolist(), [True, True, True])
        self.assertEqual(build_cell_subset_mask(obs, "malignant_like").tolist(), [False, True, True])
        self.assertEqual(build_cell_subset_mask(obs, "main_strict").tolist(), [True, False, True])

    def test_export_report_records_genes_by_cells_orientation(self):
        report = build_export_report(
            subset_name="driver_union_all",
            n_genes=3000,
            n_cells=9512,
            retained_tfs=["HNF4A", "PPARA"],
            output_files={"counts_mtx": "counts.mtx"},
        )

        self.assertEqual(report["matrix_orientation"], "genes_x_cells")
        self.assertEqual(report["shape"]["n_genes"], 3000)
        self.assertEqual(report["shape"]["n_cells"], 9512)

    def test_preranked_metric_uses_fdr_and_signed_distance(self):
        df = pd.DataFrame(
            {
                "gene": ["A", "B", "C"],
                "distance": [2.0, 1.0, -3.0],
                "p.adj": [0.01, 1.0, 0.001],
            }
        )

        ranked = compute_preranked_metric(df).set_index("gene")

        self.assertGreater(ranked.loc["A", "preranked_metric"], 0)
        self.assertEqual(ranked.loc["B", "preranked_metric"], 0)
        self.assertLess(ranked.loc["C", "preranked_metric"], 0)

    def test_mapping_stats_reports_background_and_input_mapping_rates(self):
        stats = summarize_mapping_stats(
            tf="HNF4A",
            database="KEGG",
            n_background=3000,
            n_input=120,
            n_mapped_background=2800,
            n_mapped_input=100,
        )

        self.assertAlmostEqual(stats["background_mapping_rate"], 2800 / 3000)
        self.assertAlmostEqual(stats["input_mapping_rate"], 100 / 120)

    def test_enrichment_summary_keeps_top_terms_per_tf_database(self):
        enrichment = pd.DataFrame(
            {
                "tf": ["HNF4A", "HNF4A", "HNF4A"],
                "database": ["KEGG", "KEGG", "Reactome"],
                "term_id": ["hsa00010", "hsa04110", "R-HSA-1"],
                "term_name": ["Glycolysis", "Cell cycle", "Metabolism"],
                "p.adjust": [0.02, 0.001, 0.03],
                "gene_count": [5, 8, 3],
            }
        )

        summary = build_enrichment_summary(enrichment, top_n=1)

        self.assertEqual(len(summary), 2)
        self.assertIn("Cell cycle", set(summary["term_name"]))

    def test_parse_and_filter_gene_sets_uses_background_universe(self):
        gene_sets = parse_gmt_lines(["Pathway A\tNA\tA\tB\tZ", "Pathway B\tNA\tX"])
        filtered = filter_gene_sets_to_background(gene_sets, background={"A", "B", "C"}, min_size=2, max_size=10)

        self.assertEqual(list(filtered), ["Pathway A"])
        self.assertEqual(filtered["Pathway A"], {"A", "B"})

    def test_ora_uses_custom_background_and_reports_overlap(self):
        gene_sets = {"Pathway A": {"A", "B"}, "Pathway B": {"C", "D"}}
        rows = run_ora_for_tf(
            tf="HNF4A",
            subset="driver_union_all",
            significant_genes={"A", "B"},
            background={"A", "B", "C", "D"},
            gene_sets=gene_sets,
            database="KEGG",
        )

        self.assertEqual(rows.iloc[0]["term_name"], "Pathway A")
        self.assertEqual(int(rows.iloc[0]["gene_count"]), 2)
        self.assertEqual(rows.iloc[0]["overlap_genes"], "A;B")

    def test_simple_gsea_prioritizes_genes_at_top_of_ranked_list(self):
        ranked = pd.Series({"A": 3.0, "B": 2.0, "C": -1.0, "D": -2.0})
        rows = compute_simple_gsea(
            tf="JUN",
            subset="driver_union_all",
            ranked_metric=ranked,
            gene_sets={"Top pathway": {"A", "B"}, "Bottom pathway": {"C", "D"}},
            database="Reactome_GSEA",
            permutations=0,
        ).set_index("term_name")

        self.assertGreater(rows.loc["Top pathway", "NES"], 0)
        self.assertLess(rows.loc["Bottom pathway", "NES"], 0)

    def test_concordance_summary_scores_celloracle_target_overlap(self):
        perturb = pd.DataFrame(
            {
                "tf": ["HNF4A", "HNF4A", "HNF4A", "SOX4"],
                "gene": ["A", "B", "C", "D"],
                "p.adj": [0.001, 0.02, 0.2, 0.01],
                "distance": [3.0, 2.0, 0.5, 4.0],
            }
        )
        grn = pd.DataFrame(
            {
                "source": ["HNF4A", "HNF4A", "SOX4"],
                "target": ["A", "X", "D"],
                "celloracle_state": ["malignant_or_malignant_like"] * 3,
                "p": [0.001, 0.001, 0.001],
            }
        )

        summary = build_concordance_summary(perturb, grn, fdr_threshold=0.05, top_n=2).set_index("tf")

        self.assertEqual(int(summary.loc["HNF4A", "n_significant_perturbed_genes"]), 2)
        self.assertAlmostEqual(float(summary.loc["HNF4A", "top_gene_grn_target_jaccard"]), 1 / 3)
        self.assertEqual(int(summary.loc["SOX4", "n_grn_overlap_genes"]), 1)

    def test_integrated_evidence_matrix_preserves_celloracle_tiers(self):
        concordance = pd.DataFrame(
            {
                "tf": ["HNF4A", "HLF"],
                "n_significant_perturbed_genes": [120, 5],
                "top_gene_grn_target_jaccard": [0.25, 0.0],
                "mean_distance_significant": [2.5, 0.5],
                "scTenifoldKnk_rank": [1, 2],
            }
        )
        tiers = pd.DataFrame(
            {
                "tf": ["HNF4A", "HLF"],
                "candidate_tier": ["Tier 1", "Negative/control"],
                "quantitative_rank": [1, 10],
                "quantitative_perturbation_score": [0.9, 0.49],
            }
        )

        matrix = build_integrated_evidence_matrix(concordance, tiers).set_index("tf")

        self.assertEqual(matrix.loc["HNF4A", "candidate_tier"], "Tier 1")
        self.assertGreater(matrix.loc["HNF4A", "integrated_module7_score"], matrix.loc["HLF", "integrated_module7_score"])

    def test_control_outperformance_flags_stronger_controls(self):
        matrix = pd.DataFrame(
            {
                "tf": ["HNF4A", "PPARA", "HLF"],
                "candidate_tier": ["Tier 1", "Tier 1", "Negative/control"],
                "integrated_module7_score": [0.8, 0.7, 0.9],
            }
        )

        flags = flag_control_outperformance(matrix)

        self.assertEqual(len(flags), 1)
        self.assertEqual(flags.loc[0, "risk_type"], "control_outperforms_tier1")

    def test_marker_overlap_counts_hcc_stress_and_proliferation_genes(self):
        perturb = pd.DataFrame(
            {
                "subset": ["driver_union_all", "driver_union_all", "driver_union_all", "driver_union_all"],
                "tf": ["HNF4A", "HNF4A", "JUN", "JUN"],
                "gene": ["AFP", "GPC3", "FOS", "MKI67"],
                "p.adj": [0.001, 0.2, 0.001, 0.001],
                "distance": [3.0, 1.0, 2.0, 2.5],
            }
        )

        summary = build_marker_overlap_summary(perturb, fdr_threshold=0.05).set_index(["tf", "marker_panel"])

        self.assertEqual(int(summary.loc[("HNF4A", "HCC_Malignant_Associated"), "n_significant_marker_genes"]), 1)
        self.assertEqual(int(summary.loc[("JUN", "Stressed_Injured"), "n_significant_marker_genes"]), 1)
        self.assertEqual(int(summary.loc[("JUN", "Proliferation"), "n_significant_marker_genes"]), 1)

    def test_state_specific_gene_table_finds_malignant_like_specific_genes(self):
        perturb = pd.DataFrame(
            {
                "subset": ["malignant_like", "main_strict", "driver_union_all", "malignant_like", "main_strict"],
                "tf": ["SOX4", "SOX4", "SOX4", "SOX4", "SOX4"],
                "gene": ["SPP1", "SPP1", "SPP1", "ALB", "ALB"],
                "p.adj": [0.001, 0.5, 0.2, 0.001, 0.001],
                "distance": [4.0, 0.5, 0.8, 2.0, 2.0],
            }
        )

        state = build_state_specific_gene_table(perturb, fdr_threshold=0.05, ratio_threshold=2.0)

        self.assertEqual(list(state["gene"]), ["SPP1"])
        self.assertGreater(float(state.loc[0, "malignant_like_specificity_ratio"]), 2.0)

    def test_biological_axis_summary_combines_tier_phase_and_marker_evidence(self):
        integrated = pd.DataFrame(
            {
                "tf": ["HNF4A", "JUN"],
                "candidate_tier": ["Tier 1", "Tier 3"],
                "integrated_module7_score": [0.9, 0.5],
                "module7_integrated_rank": [1, 2],
            }
        )
        marker = pd.DataFrame(
            {
                "subset": ["driver_union_all", "driver_union_all", "driver_union_all"],
                "tf": ["HNF4A", "JUN", "JUN"],
                "marker_panel": ["HCC_Malignant_Associated", "Stressed_Injured", "Proliferation"],
                "n_significant_marker_genes": [2, 3, 1],
                "significant_marker_genes": ["AFP;GPC3", "FOS;JUN;ATF3", "MKI67"],
            }
        )
        phase = pd.DataFrame(
            {
                "tf": ["HNF4A", "JUN"],
                "phase_early_rank": [15, 1],
                "phase_early_score": [0.2, 0.8],
            }
        )

        axis = build_biological_axis_summary(integrated, marker, phase).set_index("tf")

        self.assertEqual(axis.loc["HNF4A", "module7_axis_interpretation"], "Tier 1 HCC rescue replication")
        self.assertEqual(axis.loc["JUN", "module7_axis_interpretation"], "AP-1 early/stress/proliferation axis")
        self.assertEqual(int(axis.loc["JUN", "stress_marker_count"]), 3)

    def test_module7_display_group_prioritizes_main_text_axes(self):
        self.assertEqual(assign_module7_display_group("HNF4A"), "main_tier1")
        self.assertEqual(assign_module7_display_group("SOX4"), "main_state_specific")
        self.assertEqual(assign_module7_display_group("JUN"), "main_ap1_axis")
        self.assertEqual(assign_module7_display_group("HLF"), "supplement_control")
        self.assertEqual(assign_module7_display_group("IRF1"), "supplement_reserve")

    def test_tf_level_replication_matrix_adds_manuscript_roles(self):
        integrated = pd.DataFrame(
            {
                "tf": ["HNF4A", "SOX4", "HLF"],
                "candidate_tier": ["Tier 1", "Tier 2", "Negative/control"],
                "quantitative_rank": [1, 5, 10],
                "quantitative_perturbation_score": [0.9, 0.7, 0.4],
                "module7_integrated_rank": [6, 11, 13],
                "integrated_module7_score": [0.55, 0.39, 0.34],
                "n_significant_perturbed_genes": [14, 33, 14],
                "n_grn_overlap_genes": [42, 0, 7],
                "top_gene_grn_target_jaccard": [0.04, 0.0, 0.02],
            }
        )
        axis = pd.DataFrame(
            {
                "tf": ["HNF4A", "SOX4", "HLF"],
                "hcc_marker_count": [2, 1, 0],
                "stress_marker_count": [0, 0, 0],
                "proliferation_marker_count": [0, 0, 0],
                "module7_axis_interpretation": ["Tier 1 HCC rescue replication", "state axis", "control"],
            }
        )
        state_specific = pd.DataFrame({"tf": ["SOX4", "SOX4"], "gene": ["A", "B"]})

        matrix = build_tf_level_replication_matrix(integrated, axis, state_specific).set_index("tf")

        self.assertEqual(matrix.loc["HNF4A", "display_group"], "main_tier1")
        self.assertEqual(matrix.loc["SOX4", "malignant_like_state_specific_gene_count"], 2)
        self.assertEqual(matrix.loc["HLF", "manuscript_use"], "supplementary_control_calibration")

    def test_pathway_level_enrichment_matrix_labels_main_and_supplement(self):
        enrichment = pd.DataFrame(
            {
                "analysis": ["ORA", "GSEA"],
                "subset": ["driver_union_all", "driver_union_all"],
                "tf": ["JUN", "HLF"],
                "database": ["GO_BP_ORA", "Reactome_GSEA"],
                "term_name": ["Cellular Response To Heat", "Metabolism"],
                "p.adjust": [0.001, 1.0],
                "NES": [float("nan"), 0.5],
                "gene_count": [6, 10],
            }
        )
        tf_matrix = pd.DataFrame(
            {
                "tf": ["JUN", "HLF"],
                "display_group": ["main_ap1_axis", "supplement_control"],
                "candidate_tier": ["Tier 3", "Negative/control"],
            }
        )

        matrix = build_pathway_level_enrichment_matrix(enrichment, tf_matrix).set_index("tf")

        self.assertEqual(matrix.loc["JUN", "manuscript_use"], "main_pathway_panel")
        self.assertEqual(matrix.loc["HLF", "manuscript_use"], "supplementary_pathway_calibration")
        self.assertIn("pathway_rank_within_tf_database", matrix.columns)

    def test_module7_report_payload_counts_outputs(self):
        report = build_module7_report_payload(
            concordance=pd.DataFrame({"tf": ["HNF4A", "PPARA"]}),
            enrichment_summary=pd.DataFrame({"tf": ["HNF4A"], "database": ["KEGG"]}),
            risk_flags=pd.DataFrame({"risk_type": ["x"]}),
            inputs={"a": "b"},
            outputs={"c": "d"},
        )

        self.assertEqual(report["n_tfs_with_sctenifoldknk_results"], 2)
        self.assertEqual(report["n_enrichment_rows"], 1)
        self.assertEqual(report["n_review_risk_flags"], 1)


if __name__ == "__main__":
    unittest.main()
