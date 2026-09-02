from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
V2_SCRIPT = ROOT / "scripts" / "figure7_v2_core.R"


def test_figure7_v2_contract_isolated_and_prespecified():
    text = V2_SCRIPT.read_text(encoding="utf-8")
    assert "figure7_external_validation_v2" in text
    assert "FIGURE7_V2_N_RANDOM <- 1000L" in text
    assert "JUNB" in text
    assert "unsigned_associated_target_programme_score" in text
    assert "ICGC_Age_Gender_Stage_fustat_futime_encoding_unverified" in text


def test_figure7_v2_requires_patient_level_and_locked_validation_rules():
    text = V2_SCRIPT.read_text(encoding="utf-8")
    assert "patient_level_scores" in text
    assert "repeats = 10L, folds = 5L" in text
    assert "actual 3-year predicted vs observed calibration" in text
    assert "coefficients_refit_in_ICGC" in text
