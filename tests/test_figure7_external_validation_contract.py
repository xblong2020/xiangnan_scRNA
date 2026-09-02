from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"


REQUIRED = [
    "figure7_plot_theme.R",
    "figure7_00_preflight_audit.R",
    "figure7_01_prepare_bulk_expression.R",
    "figure7_02_calculate_bulk_axis_scores.R",
    "plot_figure7a_cohort_flow.R",
    "plot_figure7b_bulk_signature_mapping.R",
    "figure7_03_analyze_tumour_normal.R",
    "plot_figure7c_tumour_normal_forest.R",
    "figure7_04_analyze_clinical_associations.R",
    "plot_figure7d_clinical_heatmap.R",
    "figure7_05_fit_multivariable_cox.R",
    "plot_figure7e_multivariable_cox_forest.R",
    "figure7_06_evaluate_incremental_prediction.R",
    "plot_figure7f_incremental_prediction.R",
    "figure7_07_prepare_survival_groups.R",
    "plot_figure7g_survival_curves.R",
    "figure7_08_run_sensitivity_analyses.R",
    "plot_figure7h_sensitivity_summary.R",
    "validate_figure7_external_validation.R",
    "run_figure7_external_validation.ps1",
]


def test_required_figure7_entrypoints_exist():
    missing = [name for name in REQUIRED if not (SCRIPTS / name).is_file()]
    assert not missing, missing


def test_lancet_colour_contract_is_centralized():
    text = (SCRIPTS / "figure7_plot_theme.R").read_text(encoding="utf-8")
    assert 'ggsci::pal_lancet("lanonc")(9)' in text
    assert re.search(r"identity_loss\s*=\s*lancet_palette\[1\]", text)
    assert re.search(r"stress_transition\s*=\s*lancet_palette\[3\]", text)
    assert re.search(r"sox4_stabilization\s*=\s*lancet_palette\[2\]", text)
    for name in REQUIRED:
        if name.startswith("plot_figure7") and name.endswith(".R"):
            assert "figure7_plot_theme.R" in (SCRIPTS / name).read_text(encoding="utf-8")


def test_no_outcome_optimized_cutpoint_or_combat():
    text = (SCRIPTS / "figure7_core.R").read_text(encoding="utf-8").lower()
    operational = text.split("stage_validate <-", maxsplit=1)[0]
    assert "surv_cutpoint" not in operational
    assert "maximally selected" not in operational
    assert "combat(" not in operational
    assert "weights_locked = true" in text


def test_random_signature_and_validation_counts_are_prespecified():
    text = (SCRIPTS / "figure7_core.R").read_text(encoding="utf-8")
    assert "n_random = 500L" in text
    assert "n_boot = 500L" in text
    assert "repeats = 10L, folds = 5L" in text
    assert "FIGURE7_SEED <- 20260805L" in text


def test_figure7_paths_are_namespaced():
    text = (SCRIPTS / "figure7_core.R").read_text(encoding="utf-8")
    assert '"figure7_external_validation"' in text
    assert "figure7_protected_figure1_6_hashes_before.tsv" in text
    assert "figure7_validation_report.tsv" in text
