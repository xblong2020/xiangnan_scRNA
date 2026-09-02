import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from scripts.external_validation_module8 import (
    MODULE8_REQUIRED_MANIFEST_COLUMNS,
    audit_dataset_leakage,
    build_axis_level_evidence_grade,
    build_external_dataset_manifest,
    build_signature_registry,
    build_tf_target_signature_genes,
    compute_bulk_signature_scores,
    compute_adjusted_survival_association,
    compute_bulk_clinical_variable_association,
    compute_bulk_tumor_normal_association,
    compute_exploratory_survival_association,
    collapse_scores_to_axis,
    compute_group_recurrence,
    compute_signature_scores,
    discover_local_scrna_sources,
    flag_external_control_outperformance,
    infer_comparison_group,
    infer_icgc_sample_type,
    load_tcga_clinical_table,
    prepare_tcga_clinical_covariates,
    sample_to_patient_id,
    infer_tcga_sample_type,
    normalize_gene_symbols,
    strip_ensembl_version,
    summarize_scrna_dataset_scores,
    write_module8_nature_figures,
)


class Module8ExternalValidationLogicTest(unittest.TestCase):
    def test_manifest_has_required_columns_and_excludes_discovery_datasets(self):
        manifest = build_external_dataset_manifest(discovery_dataset_ids={"GSE151530", "GSE149614"})

        self.assertTrue(set(MODULE8_REQUIRED_MANIFEST_COLUMNS).issubset(manifest.columns))
        audited = audit_dataset_leakage(manifest, discovery_dataset_ids={"GSE151530", "GSE149614"})
        leaked = audited.loc[audited["dataset_id"].isin(["GSE151530", "GSE149614"])]

        self.assertTrue(leaked["included"].eq(False).all())
        self.assertTrue(leaked["exclusion_reason"].str.contains("discovery_or_lodo_dataset").all())

    def test_signature_registry_freezes_module7_axis_membership(self):
        tf_matrix = pd.DataFrame(
            {
                "tf": ["HNF4A", "PPARA", "JUN", "SOX4", "HLF"],
                "display_group": [
                    "main_tier1",
                    "main_tier1",
                    "main_ap1_axis",
                    "main_state_specific",
                    "supplement_control",
                ],
                "candidate_tier": ["Tier 1", "Tier 1", "Tier 3", "Tier 2", "Negative/control"],
            }
        )

        registry = build_signature_registry(tf_matrix).set_index("tf")

        self.assertEqual(registry.loc["HNF4A", "axis"], "tier1_rescue")
        self.assertEqual(registry.loc["JUN", "axis"], "ap1_stress_proliferation")
        self.assertEqual(registry.loc["SOX4", "axis"], "sox4_state_specific")
        self.assertEqual(registry.loc["HLF", "axis"], "control_calibration")

    def test_tf_target_signature_genes_selects_top_disturbed_genes(self):
        registry = pd.DataFrame(
            {
                "tf": ["HNF4A"],
                "axis": ["tier1_rescue"],
                "signature_class": ["main"],
            }
        )
        perturb = pd.DataFrame(
            {
                "tf": ["HNF4A", "HNF4A", "HNF4A"],
                "gene": ["A", "B", "C"],
                "distance": [5.0, 2.0, 1.0],
                "p.adj": [0.01, 0.02, 0.5],
            }
        )

        genes = build_tf_target_signature_genes(registry, perturb, top_n=2)

        self.assertEqual(genes["gene"].tolist(), ["A", "B"])
        self.assertEqual(genes["signature_source"].unique().tolist(), ["sctenifold_top_disturbed"])

    def test_signature_scores_use_gene_z_scores_and_report_detection_rate(self):
        expression = pd.DataFrame(
            {
                "cell1": [10.0, 8.0, 1.0],
                "cell2": [9.0, 7.0, 1.0],
                "cell3": [1.0, 2.0, 10.0],
            },
            index=["A", "B", "C"],
        )
        signatures = pd.DataFrame(
            {
                "axis": ["axis1", "axis1"],
                "tf": ["TF1", "TF1"],
                "gene": ["A", "B"],
            }
        )

        scores, detection = compute_signature_scores(expression, signatures, dataset_id="EXT1")

        self.assertGreater(scores.loc[scores["cell_id"].eq("cell1"), "signature_score"].iloc[0], 0)
        self.assertLess(scores.loc[scores["cell_id"].eq("cell3"), "signature_score"].iloc[0], 0)
        self.assertEqual(float(detection.loc[0, "gene_detection_rate"]), 1.0)

    def test_group_recurrence_reports_effect_direction_and_adjusted_pvalue(self):
        scores = pd.DataFrame(
            {
                "dataset_id": ["EXT1"] * 6,
                "axis": ["axis1"] * 6,
                "cell_id": [f"c{i}" for i in range(6)],
                "signature_score": [2.0, 1.8, 1.5, -1.0, -1.2, -1.5],
            }
        )
        metadata = pd.DataFrame(
            {
                "cell_id": [f"c{i}" for i in range(6)],
                "comparison_group": ["malignant", "malignant", "malignant", "reference", "reference", "reference"],
            }
        )

        recurrence = compute_group_recurrence(
            scores,
            metadata,
            dataset_id="EXT1",
            modality="scRNA-seq",
            positive_group="malignant",
            reference_group="reference",
        )

        self.assertEqual(recurrence.loc[0, "direction"], "positive")
        self.assertGreater(float(recurrence.loc[0, "effect_size"]), 0)
        self.assertIn("p.adjust", recurrence.columns)

    def test_collapse_scores_to_axis_averages_tf_scores_per_cell(self):
        scores = pd.DataFrame(
            {
                "dataset_id": ["EXT1", "EXT1", "EXT1"],
                "axis": ["axis1", "axis1", "axis2"],
                "tf": ["TF1", "TF2", "TF3"],
                "cell_id": ["c1", "c1", "c1"],
                "signature_score": [1.0, 3.0, -1.0],
            }
        )

        collapsed = collapse_scores_to_axis(scores).set_index("axis")

        self.assertEqual(float(collapsed.loc["axis1", "signature_score"]), 2.0)
        self.assertEqual(float(collapsed.loc["axis2", "signature_score"]), -1.0)

    def test_summarize_scrna_dataset_scores_keeps_required_statistics(self):
        scores = pd.DataFrame(
            {
                "dataset_id": ["EXT1", "EXT1"],
                "axis": ["axis1", "axis1"],
                "tf": ["TF1", "TF1"],
                "cell_id": ["c1", "c2"],
                "signature_score": [1.0, -1.0],
            }
        )
        metadata = pd.DataFrame({"cell_id": ["c1", "c2"], "comparison_group": ["malignant", "reference"]})
        detection = pd.DataFrame(
            {
                "dataset_id": ["EXT1"],
                "axis": ["axis1"],
                "tf": ["TF1"],
                "gene_detection_rate": [0.8],
            }
        )
        recurrence = pd.DataFrame(
            {
                "axis": ["axis1"],
                "effect_size": [2.0],
                "pvalue": [0.1],
                "p.adjust": [0.1],
                "direction": ["positive"],
                "missingness_rate": [0.0],
            }
        )

        summary = summarize_scrna_dataset_scores(scores, metadata, detection, recurrence)

        self.assertIn("mean_signature_score", summary.columns)
        self.assertEqual(set(summary["cell_state"]), {"malignant", "reference"})
        self.assertEqual(float(summary.loc[summary["cell_state"].eq("malignant"), "gene_detection_rate"].iloc[0]), 0.8)

    def test_axis_evidence_grade_requires_modality_recurrence(self):
        recurrence = pd.DataFrame(
            {
                "axis": ["axis1", "axis1", "axis2"],
                "modality": ["scRNA-seq", "bulk_RNA-seq", "scRNA-seq"],
                "dataset_id": ["EXT1", "TCGA-LIHC", "EXT1"],
                "direction": ["positive", "positive", "negative"],
                "effect_size": [1.0, 0.8, -0.5],
                "p.adjust": [0.01, 0.02, 0.2],
                "included": [True, True, True],
                "dataset_leakage_status": ["clean", "clean", "clean"],
            }
        )

        grades = build_axis_level_evidence_grade(recurrence).set_index("axis")

        self.assertEqual(grades.loc["axis1", "evidence_grade"], "B")
        self.assertEqual(grades.loc["axis2", "evidence_grade"], "D")
        self.assertEqual(int(grades.loc["axis1", "n_replicated_cohorts"]), 2)

    def test_control_outperformance_flags_control_axis_above_main_axes(self):
        grades = pd.DataFrame(
            {
                "axis": ["tier1_rescue", "control_calibration"],
                "axis_group": ["main", "control"],
                "recurrence_score": [0.4, 0.9],
            }
        )

        flags = flag_external_control_outperformance(grades)

        self.assertEqual(len(flags), 1)
        self.assertEqual(flags.loc[0, "risk_type"], "review_risk_external_control_outperformance")

    def test_normalize_gene_symbols_strips_suffix_and_uppercases(self):
        genes = ["hnf4a", "JUN.1", np.nan, " sox4 "]

        self.assertEqual(normalize_gene_symbols(genes), ["HNF4A", "JUN", "", "SOX4"])

    def test_strip_ensembl_version_removes_suffix(self):
        self.assertEqual(strip_ensembl_version("ENSG00000167578.15"), "ENSG00000167578")

    def test_infer_tcga_sample_type_uses_barcode_portion(self):
        samples = ["TCGA-XX-0001-01A", "TCGA-XX-0002-11A", "OTHER"]

        self.assertEqual(infer_tcga_sample_type(samples), ["tumor", "normal", "unknown"])

    def test_infer_icgc_sample_type_uses_suffix(self):
        samples = ["SP112215-DO50816-T", "SP135251-DO45167-N", "OTHER"]

        self.assertEqual(infer_icgc_sample_type(samples), ["tumor", "normal", "unknown"])

    def test_sample_to_patient_id_supports_tcga_and_icgc(self):
        self.assertEqual(sample_to_patient_id("TCGA-DD-AAE4-01A"), "TCGA-DD-AAE4")
        self.assertEqual(sample_to_patient_id("SP112215-DO50816-T"), "DO50816")

    def test_bulk_signature_scores_use_sample_z_scores(self):
        expression = pd.DataFrame(
            {
                "S1": [10.0, 8.0, 1.0],
                "S2": [9.0, 7.0, 1.0],
                "S3": [1.0, 2.0, 10.0],
            },
            index=["A", "B", "C"],
        )
        signatures = pd.DataFrame({"axis": ["axis1", "axis1"], "tf": ["TF1", "TF1"], "gene": ["A", "B"]})

        scores, detection = compute_bulk_signature_scores(expression, signatures, dataset_id="TCGA-LIHC")

        self.assertGreater(scores.loc[scores["sample"].eq("S1"), "signature_score"].iloc[0], 0)
        self.assertLess(scores.loc[scores["sample"].eq("S3"), "signature_score"].iloc[0], 0)
        self.assertEqual(float(detection.loc[0, "gene_detection_rate"]), 1.0)

    def test_bulk_tumor_normal_association_reports_required_statistics(self):
        scores = pd.DataFrame(
            {
                "dataset_id": ["TCGA-LIHC"] * 6,
                "axis": ["axis1"] * 6,
                "sample": [f"S{i}" for i in range(6)],
                "signature_score": [2.0, 1.8, 1.5, -1.0, -1.2, -1.5],
                "sample_type": ["tumor", "tumor", "tumor", "normal", "normal", "normal"],
            }
        )

        assoc = compute_bulk_tumor_normal_association(scores)

        self.assertEqual(assoc.loc[0, "direction"], "positive")
        self.assertGreater(float(assoc.loc[0, "effect_size"]), 0)
        self.assertIn("missingness_rate", assoc.columns)

    def test_exploratory_survival_association_marks_missing_covariates(self):
        scores = pd.DataFrame(
            {
                "dataset_id": ["TCGA-LIHC"] * 4,
                "axis": ["axis1"] * 4,
                "sample": ["S1", "S2", "S3", "S4"],
                "signature_score": [2.0, 1.0, -1.0, -2.0],
            }
        )
        survival = pd.DataFrame({"sample": ["S1", "S2", "S3", "S4"], "OS": [1, 0, 1, 0], "OS.time": [10, 20, 30, 40]})

        summary = compute_exploratory_survival_association(scores, survival)

        self.assertEqual(summary.loc[0, "status"], "missing_covariates_exploratory")
        self.assertIn("hazard_ratio", summary.columns)

    def test_adjusted_survival_association_uses_available_covariates(self):
        rng = np.random.default_rng(0)
        n = 80
        scores = pd.DataFrame(
            {
                "dataset_id": ["ICGC-LIRI-JP"] * n,
                "axis": ["axis1"] * n,
                "sample": [f"SP{i}-DO{i}-T" for i in range(n)],
                "signature_score": rng.normal(size=n),
            }
        )
        survival = pd.DataFrame(
            {
                "id": [f"DO{i}" for i in range(n)],
                "fustat": rng.binomial(1, 0.55, n),
                "futime": rng.exponential(scale=100, size=n) + 1,
            }
        )
        clinical = pd.DataFrame(
            {
                "Id": [f"DO{i}" for i in range(n)],
                "Age": rng.integers(0, 2, n),
                "Gender": rng.integers(0, 2, n),
                "Stage": rng.integers(0, 3, n),
            }
        )

        summary = compute_adjusted_survival_association(
            scores,
            survival,
            clinical,
            dataset_id="ICGC-LIRI-JP",
            survival_id_col="id",
            event_col="fustat",
            time_col="futime",
            clinical_id_col="Id",
            covariates=["Age", "Gender", "Stage"],
        )

        self.assertEqual(summary.loc[0, "status"], "adjusted_cox")
        self.assertIn("Age", summary.loc[0, "covariates_used"])
        self.assertIn("hazard_ratio", summary.columns)

    def test_adjusted_survival_association_uses_tumor_samples_when_available(self):
        scores = pd.DataFrame(
            {
                "dataset_id": ["TCGA-LIHC"] * 8,
                "axis": ["axis1"] * 8,
                "sample": [
                    "TCGA-A-0001-01A",
                    "TCGA-A-0002-01A",
                    "TCGA-A-0003-01A",
                    "TCGA-A-0004-01A",
                    "TCGA-A-0005-11A",
                    "TCGA-A-0006-11A",
                    "TCGA-A-0007-11A",
                    "TCGA-A-0008-11A",
                ],
                "sample_type": ["tumor", "tumor", "tumor", "tumor", "normal", "normal", "normal", "normal"],
                "signature_score": [0.1, 0.2, 1.0, 1.2, 10.0, 9.0, 8.0, 7.0],
            }
        )
        survival = pd.DataFrame(
            {
                "Id": [f"TCGA-A-000{i}" for i in range(1, 9)],
                "fustat": [0, 1, 0, 1, 0, 1, 0, 1],
                "futime": [10, 20, 30, 40, 10, 20, 30, 40],
            }
        )
        clinical = pd.DataFrame(
            {
                "Id": [f"TCGA-A-000{i}" for i in range(1, 9)],
                "age": [50, 55, 60, 65, 50, 55, 60, 65],
                "gender": [0, 1, 0, 1, 0, 1, 0, 1],
                "stage": [1, 1, 2, 2, 1, 1, 2, 2],
            }
        )

        summary = compute_adjusted_survival_association(
            scores,
            survival,
            clinical,
            dataset_id="TCGA-LIHC",
            survival_id_col="Id",
            event_col="fustat",
            time_col="futime",
            clinical_id_col="Id",
            covariates=["age", "gender", "stage"],
        )

        self.assertEqual(int(summary.loc[0, "n_samples"]), 4)

    def test_bulk_clinical_variable_association_compares_binary_one_vs_zero_groups(self):
        scores = pd.DataFrame(
            {
                "dataset_id": ["TCGA-LIHC"] * 4,
                "axis": ["axis1"] * 4,
                "sample": ["TCGA-A-0001-01A", "TCGA-A-0002-01A", "TCGA-A-0003-01A", "TCGA-A-0004-01A"],
                "signature_score": [0.1, 0.2, 1.0, 1.2],
            }
        )
        clinical = pd.DataFrame({"Id": ["TCGA-A-0001", "TCGA-A-0002", "TCGA-A-0003", "TCGA-A-0004"], "stage": [0, 0, 1, 1]})

        assoc = compute_bulk_clinical_variable_association(scores, clinical, dataset_id="TCGA-LIHC", variables=["stage"])

        self.assertEqual(assoc.loc[0, "clinical_variable"], "stage_1_vs_0")
        self.assertGreater(float(assoc.loc[0, "effect_size"]), 0)

    def test_bulk_clinical_variable_association_uses_tumor_samples_when_available(self):
        scores = pd.DataFrame(
            {
                "dataset_id": ["TCGA-LIHC"] * 4,
                "axis": ["axis1"] * 4,
                "sample": ["TCGA-A-0001-01A", "TCGA-A-0002-01A", "TCGA-A-0003-11A", "TCGA-A-0004-11A"],
                "sample_type": ["tumor", "tumor", "normal", "normal"],
                "signature_score": [0.1, 1.0, 10.0, 0.0],
            }
        )
        clinical = pd.DataFrame({"Id": ["TCGA-A-0001", "TCGA-A-0002", "TCGA-A-0003", "TCGA-A-0004"], "stage": [0, 1, 1, 0]})

        assoc = compute_bulk_clinical_variable_association(scores, clinical, dataset_id="TCGA-LIHC", variables=["stage"])

        self.assertEqual(int(assoc.loc[0, "n_samples"]), 2)
        self.assertGreater(float(assoc.loc[0, "effect_size"]), 0)

    def test_load_tcga_clinical_table_reads_tab_text_with_xls_suffix(self):
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "clinical.xls"
            path.write_text(
                "Id\tfutime\tfustat\tage\tgender\tgrade\tstage\tT\tM\tN\n"
                "TCGA-DD-A1EA\t2415\t0\t68\tMALE\tG2\tStage II\tT2\tM0\tN0\n"
                "TCGA-CC-A5UC\t347\t1\t63\tFEMALE\tG3\tStage IIIA\tT3\tM0\tN0\n",
                encoding="utf-8",
            )

            clinical = load_tcga_clinical_table(path)

        self.assertEqual(clinical.loc[0, "Id"], "TCGA-DD-A1EA")
        self.assertIn("stage", clinical.columns)
        self.assertEqual(len(clinical), 2)

    def test_prepare_tcga_clinical_covariates_encodes_full_clinical_fields(self):
        clinical = pd.DataFrame(
            {
                "Id": ["TCGA-DD-A1EA", "TCGA-CC-A5UC"],
                "age": [68, 63],
                "gender": ["MALE", "FEMALE"],
                "grade": ["G2", "G3"],
                "stage": ["Stage II", "Stage IIIA"],
                "T": ["T2", "T3"],
                "M": ["M0", "M0"],
                "N": ["N0", "N1"],
            }
        )

        prepared = prepare_tcga_clinical_covariates(clinical)

        self.assertEqual(prepared["gender"].tolist(), [1, 0])
        self.assertEqual(prepared["stage"].tolist(), [2, 3])
        self.assertEqual(prepared["grade"].tolist(), [2, 3])
        self.assertEqual(prepared["T"].tolist(), [2, 3])
        self.assertEqual(prepared["N"].tolist(), [0, 1])

    def test_write_module8_nature_figures_exports_figure8_panels(self):
        manifest = pd.DataFrame(
            {
                "dataset_id": ["GSE156625", "TCGA-LIHC", "ICGC-LIRI-JP", "GSE151530"],
                "modality": ["scRNA-seq", "bulk RNA-seq", "bulk RNA-seq", "scRNA-seq"],
                "included": [True, True, True, False],
                "exclusion_reason": ["", "", "", "discovery_or_lodo_dataset"],
            }
        )
        axis_grade = pd.DataFrame(
            {
                "axis": ["tier1_rescue", "ap1_stress_proliferation", "sox4_state_specific", "control_calibration"],
                "direction_consistency": [1.0, 1.0, 1.0, 0.8],
                "effect_size_median": [0.4, 0.3, 0.5, 0.2],
                "fdr_meta": [0.001, 0.002, 0.003, 0.02],
                "n_replicated_cohorts": [6, 6, 6, 5],
                "evidence_grade": ["B", "B", "B", "B"],
                "recurrence_score": [2.4, 2.1, 2.9, 2.0],
            }
        )
        scrna_recurrence = pd.DataFrame(
            {
                "dataset_id": ["GSE156625", "GSE156625", "CNP0000650", "CNP0000650"],
                "axis": ["tier1_rescue", "sox4_state_specific", "tier1_rescue", "sox4_state_specific"],
                "effect_size": [0.1, 0.2, 0.3, 0.4],
                "direction": ["positive", "positive", "positive", "positive"],
                "p.adjust": [0.01, 0.02, 0.03, 0.04],
            }
        )
        bulk_assoc = pd.DataFrame(
            {
                "dataset_id": ["TCGA-LIHC", "TCGA-LIHC", "ICGC-LIRI-JP", "ICGC-LIRI-JP"],
                "axis": ["tier1_rescue", "sox4_state_specific", "tier1_rescue", "sox4_state_specific"],
                "clinical_variable": ["tumor_vs_normal", "grade_high_vs_low", "tumor_vs_normal", "Stage_1_vs_0"],
                "effect_size": [0.4, 0.2, 0.5, 0.3],
                "hazard_ratio": [np.nan, np.nan, np.nan, np.nan],
                "p.adjust": [0.001, 0.02, 0.002, 0.03],
                "direction": ["positive", "positive", "positive", "positive"],
                "status": ["tested", "tested", "tested", "tested"],
            }
        )
        survival = pd.DataFrame(
            {
                "dataset_id": ["TCGA-LIHC", "ICGC-LIRI-JP"],
                "axis": ["tier1_rescue", "tier1_rescue"],
                "hazard_ratio": [4.0, 4.4],
                "p.adjust": [0.01, 0.02],
                "status": ["adjusted_cox", "adjusted_cox"],
            }
        )
        with TemporaryDirectory() as tmpdir:
            outputs = write_module8_nature_figures(
                manifest=manifest,
                axis_grade=axis_grade,
                figure_dir=Path(tmpdir),
                scrna_recurrence=scrna_recurrence,
                bulk_assoc=bulk_assoc,
                survival_summary=survival,
            )

            expected_prefixes = [
                "figure8_dataset_flow",
                "scrna_heatmap",
                "bulk_forestplot",
                "clinical_heatmap",
                "axis_recurrence_summary",
                "figure8_multipanel",
            ]
            for prefix in expected_prefixes:
                for suffix in ["png", "pdf", "svg"]:
                    key = f"{prefix}_{suffix}"
                    self.assertIn(key, outputs)
                    self.assertGreater(Path(outputs[key]).stat().st_size, 1000)

    def test_infer_comparison_group_uses_known_hcc_metadata_columns(self):
        metadata = pd.DataFrame(
            {
                "NormalvsTumor": ["T", "N", ""],
                "Type": ["T cell", "Malignant cell", "CAF"],
                "cell_type": ["C1_Tcell", "C10_Tumor", "C4_NK"],
                "tissue_source": ["Tumor", "Adjacent liver", "Tumor"],
            }
        )

        groups = infer_comparison_group(metadata)

        self.assertEqual(groups.tolist(), ["malignant", "reference", "reference"])

    def test_infer_comparison_group_does_not_misclassify_malignant_cell_as_t_cell(self):
        metadata = pd.DataFrame({"Type": ["Malignant cell", "T cell", "HPC-like"]})

        groups = infer_comparison_group(metadata)

        self.assertEqual(groups.tolist(), ["malignant", "reference", "reference"])

    def test_discover_local_scrna_sources_reports_supported_formats(self):
        sources = discover_local_scrna_sources(r"G:\wanyi_HCC_scRNA\HCCscRNA")

        self.assertIn("dataset_id", sources.columns)
        self.assertIn("input_format", sources.columns)


if __name__ == "__main__":
    unittest.main()
