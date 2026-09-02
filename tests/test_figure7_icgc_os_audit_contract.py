from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "figure7_icgc_os_audit.py"


def test_icgc_os_audit_is_read_only_and_namespaced():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "Figure7_ICGC_OS_Audit" in text
    assert "ICGC_OS_DERIVATION_SPEC.md" in text
    assert "ICGC_OS_UNBLOCK_GATE.md" in text
    assert "ICGC_expression_clinical_mapping.tsv" in text
    assert "ICGC_survival_QC_summary.tsv" in text
    assert "coxph" not in text.lower()
    assert "survfit" not in text.lower()


def test_icgc_os_audit_records_the_three_allowed_decisions():
    text = SCRIPT.read_text(encoding="utf-8")
    for status in (
        "UNBLOCKED_FOR_EXTERNAL_SURVIVAL_VALIDATION",
        "ESTIMABLE_BUT_NOT_VALIDATED",
        "REMAIN_BLOCKED",
    ):
        assert status in text
