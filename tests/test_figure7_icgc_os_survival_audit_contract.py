from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "figure7_icgc_os_survival_audit.R"


def test_survival_audit_is_isolated_and_uses_frozen_scores():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "Figure7_ICGC_OS_Audit" in text
    assert "primary_frozen_programme" in (ROOT / "scripts" / "figure7_icgc_os_audit.py").read_text(encoding="utf-8")
    assert "Surv(os_time_days, os_event)" in text
    assert "scale(frozen_axis_score)" in text
    assert "age omitted" in text
    assert "SUPPLEMENTARY_ONLY" in text


def test_survival_audit_records_ph_and_nonlinearity_checks():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "cox.zph" in text
    assert "natural-spline" in text
    assert "global_PH_P" in text
    assert "ICGC_OS_FINAL_DECISION.md" in text
