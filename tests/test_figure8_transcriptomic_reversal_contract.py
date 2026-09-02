from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "metadata" / "driver" / "figure8_transcriptomic_reversal"


def read_tsv(name: str) -> pd.DataFrame:
    path = META / name
    assert path.exists(), f"Missing Figure 8 output: {path}"
    return pd.read_csv(path, sep="\t")


def test_signature_variants_and_primary_contract() -> None:
    manifest = read_tsv("figure8_signature_variant_manifest.tsv")
    assert manifest["signature_id"].nunique() == 15
    primary = manifest.loc[manifest["signature_id"] == "primary_three_axis"].iloc[0]
    assert int(primary["up_gene_count"]) == 150
    assert int(primary["down_gene_count"]) == 150
    assert int(primary["landmark_up_count"] + primary["landmark_down_count"]) == 47


def test_random_and_cross_method_null_contracts() -> None:
    random_manifest = read_tsv("figure8_random_signature_manifest.tsv")
    assert random_manifest["signature_id"].nunique() == 1000

    overlap = read_tsv("figure8e_method_overlap.tsv")
    three_way = overlap.loc[
        overlap["DrugReflector"].astype(bool)
        & overlap["L1000FWD"].astype(bool)
        & overlap["CLUE"].astype(bool)
    ]
    assert len(three_way) == 1
    assert int(three_way.iloc[0]["count"]) == 0


def test_external_connectivity_is_not_direct_expression() -> None:
    external = read_tsv("figure8g_external_signature_scores.tsv.gz")
    assert not external["direct_axis_score_available"].astype(bool).any()
    assert external["identity_rescue_score"].isna().all()
    assert external["stress_suppression_score"].isna().all()
    assert external["sox4_suppression_score"].isna().all()


def test_integrated_coverage_and_validation_contracts() -> None:
    ranking = read_tsv("figure8h_candidate_ranking_full.tsv")
    assert ranking["evidence_coverage"].between(0, 1, inclusive="both").all()
    assert ranking["toxicity_unknown"].astype(bool).all()

    validation = read_tsv("figure8_validation_report.tsv")
    assert (validation["status"] == "fail").sum() == 0
    assert (validation["status"] == "review_risk").sum() == 1

