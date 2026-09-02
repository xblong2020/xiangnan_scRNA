from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "metadata" / "driver" / "figure8_transcriptomic_reversal_v2_mainfigure"
DATA = ROOT / "data" / "processed" / "driver" / "figure8_transcriptomic_reversal_v2_mainfigure"
FIG = ROOT / "figures" / "driver" / "figure8_transcriptomic_reversal_v2_mainfigure"
REPORT = ROOT / "reports" / "figure8_transcriptomic_reversal_v2_mainfigure"


def read_tsv(path: Path) -> pd.DataFrame:
    assert path.exists(), f"Missing v2 artifact: {path}"
    return pd.read_csv(path, sep="\t")


def test_final_landmark_and_signature_contract() -> None:
    coverage = read_tsv(META / "figure8_v2_978_landmark_coverage.tsv")
    assert len(coverage) == 978
    assert coverage["gene"].nunique() == 978
    scores = read_tsv(META / "figure8_v2_gene_level_rescue_vscore.tsv")
    assert len(scores) == 978
    assert scores["final_rescue_vscore"].between(-1, 1).all()
    assert not scores["evidence_sources"].str.contains("compound|DrugReflector_rank", case=False, regex=True).any()


def test_final_random_and_mapping_contract() -> None:
    random = read_tsv(META / "figure8_v2_matched_random_benchmark.tsv.gz")
    assert random["signature_id"].nunique() >= 2000
    summary = read_tsv(META / "figure8_v2_random_specificity_summary.tsv")
    assert summary["empirical_p_two_sided"].between(0, 1).all()
    mapping = read_tsv(META / "figure8_v2_cross_framework_concordance.tsv")
    assert not mapping.get("fuzzy_match_used", pd.Series(False, index=mapping.index)).astype(bool).any()


def test_final_prism_and_context_contract() -> None:
    cells = read_tsv(META / "figure8_v2_cell_line_metadata_audit.tsv")
    by_name = cells.set_index(cells["cell_line"].str.upper())
    assert by_name.loc["HCC515", "verified_context"] == "lung_adenocarcinoma_non_liver"
    assert by_name.loc["HA1E", "verified_context"] == "kidney_derived_non_liver"
    assert "hepatoblastoma" in by_name.loc["HEPG2", "verified_context"]
    prism = read_tsv(META / "figure8_v2_prism_viability.tsv")
    assert "normal_cell_safety_established" in prism.columns
    assert not prism["normal_cell_safety_established"].fillna(False).astype(bool).any()


def test_final_figure_report_and_validation_contract() -> None:
    for suffix in ("pdf", "svg", "png", "tiff"):
        path = FIG / f"figure8_v2_mainfigure_a_to_g.{suffix}"
        assert path.exists() and path.stat().st_size > 0
    assert (REPORT / "figure8_v2_mainfigure_reanalysis_report.md").exists()
    validation = read_tsv(META / "figure8_v2_validation_report.tsv")
    assert (validation["status"] == "fail").sum() == 0
    assert len(validation) >= 30

