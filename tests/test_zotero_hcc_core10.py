from pathlib import Path

from zotero_hcc_classifier.core10 import (
    CandidateRecord,
    Core10Rule,
    has_hcc_anchor,
    is_excluded_from_core10,
    is_review_record,
    normalize_pmid,
    deduplicate_candidates,
    build_pubmed_query,
    score_candidate,
    select_core10,
)
from zotero_hcc_classifier.core10_preview import run_core10_preview
from zotero_hcc_classifier.core10_apply import _created_item_key, decide_import
from zotero_hcc_classifier.pdf_acquisition import legal_pdf_url, validate_pdf
from zotero_hcc_classifier.zotero_client import build_attachment_payload


def record(**overrides):
    base = {
        "title": "HCC paper",
        "abstract": "hepatocellular carcinoma study",
        "year": 2024,
        "doi": "",
        "pmid": "",
        "pmcid": "",
        "article_type": "Original Research",
        "citation_count": 0,
        "authors": ["A Author"],
        "journal": "Journal of Hepatology",
        "source": "PubMed",
    }
    base.update(overrides)
    return CandidateRecord(**base)


def rule_for(name):
    return Core10Rule(
        level2="02_Identity-Loss",
        level3=name,
        keywords=("HNF4A", "hepatocyte identity"),
        evidence_terms=("knockdown", "rescue", "functional assay", "in vivo"),
        seminal_terms=("identity",),
    )


def test_pubmed_query_contains_hcc_anchor_and_leaf_terms():
    query = build_pubmed_query("01_HNF4A", ["HNF4A", "hepatocyte identity"])
    assert "hepatocellular carcinoma" in query
    assert "HNF4A" in query


def test_deduplication_prefers_pmid_then_doi_then_normalized_title():
    records = [
        record(pmid="1", doi="10/a", title="Same"),
        record(pmid="1", doi="10/b", title="Same"),
        record(title="Same"),
    ]
    result = deduplicate_candidates(records)
    assert len(result) == 1
    assert result[0].pmid == "1"


def test_candidate_filter_rejects_non_hcc_broad_search_hits():
    assert has_hcc_anchor(record(title="Survival in esophageal squamous cell carcinoma", abstract="Cox regression")) is False
    assert has_hcc_anchor(record(title="HNF4A identity loss in HCC", abstract="hepatocellular carcinoma")) is True


def test_review_and_reply_article_types_are_normalized_before_core_selection():
    assert is_review_record(record(title="HCC identity overview (Review)", article_type="Original Research")) is True
    assert is_review_record(record(title="AASLD guidelines for the treatment of hepatocellular carcinoma", article_type="Original Research")) is True
    assert is_excluded_from_core10(record(title="Reply to: HCC clinical study", article_type="Original Research")) is True
    assert is_excluded_from_core10(record(title="Abstract 123: HCC study", article_type="Original Research")) is True


def test_pmid_urls_are_normalized_to_numeric_identifiers():
    assert normalize_pmid("/pubmed.ncbi.nlm.nih.gov/29907753") == "29907753"
    assert normalize_pmid("29907753") == "29907753"


def test_core_score_has_five_bounded_components_and_year_normalized_influence():
    scored = score_candidate(
        record(
            title="HNF4A identity loss in HCC",
            abstract="knockdown and rescue with functional assay in vivo",
            year=2024,
            citation_count=20,
        ),
        rule=rule_for("01_HNF4A"),
        current_year=2026,
    )
    assert 0 <= scored.total <= 100
    assert 0 <= scored.topic_score <= 40
    assert 0 <= scored.evidence_score <= 25
    assert 0 <= scored.influence_score <= 15
    assert 0 <= scored.quality_score <= 10
    assert 0 <= scored.recency_score <= 10


def test_abstract_topic_hits_receive_substantial_topic_weight():
    rule = rule_for("01_HNF4A")
    title_only = score_candidate(record(title="HCC paper", abstract="hepatocellular carcinoma", year=2024), rule=rule, current_year=2026)
    abstract_core = score_candidate(record(title="HCC paper", abstract="hepatocellular carcinoma HNF4A", year=2024), rule=rule, current_year=2026)
    assert abstract_core.topic_score - title_only.topic_score >= 6


def test_select_core10_keeps_at_most_three_reviews_and_requires_75():
    selected = select_core10(
        [
            record(title=f"Original {i} HNF4A identity loss in HCC", abstract="knockdown rescue functional assay in vivo", year=2024, citation_count=100 - i)
            for i in range(8)
        ]
        + [record(title=f"Review {i} HNF4A identity loss in HCC", abstract="review of HNF4A identity", article_type="Review", citation_count=200 - i) for i in range(4)]
        + [record(title="Weak HCC paper", abstract="background only", year=2010, citation_count=0)],
        limit=10,
    )
    assert all(item.total >= 75 for item in selected)
    assert sum(item.record.article_type == "Review" for item in selected) <= 3


def test_core10_preview_writes_report_from_acquired_candidates(tmp_path):
    (tmp_path / "data").mkdir(parents=True)
    (tmp_path / "data" / "items_before.json").write_text(json_text({"items": []}), encoding="utf-8")

    def fake_acquire(rules, **kwargs):
        return {rules[0].path: [record(title="HNF4A identity loss in HCC", abstract="knockdown rescue functional assay in vivo", year=2024, citation_count=100)]}, {"queries": [], "errors": []}

    result = run_core10_preview(output_root=tmp_path, acquire_fn=fake_acquire)
    assert result["qa"]["category_count"] == 32
    assert (tmp_path / "reports" / "CORE10_PREVIEW.md").exists()


def json_text(value):
    import json

    return json.dumps(value)


def test_existing_doi_is_reused_and_existing_collections_and_tags_are_preserved():
    existing = {
        "key": "EXIST1",
        "version": 4,
        "data": {
            "title": "Existing HCC paper",
            "DOI": "10/a",
            "collections": ["OLD"],
            "tags": [{"tag": "Project::Existing"}],
        },
    }
    action = decide_import(
        existing_items=[existing],
        candidate=record(doi="10/a", title="Existing HCC paper"),
        target_collection="TARGET",
        target_tag="Core10::01_HNF4A",
    )
    assert action.kind == "REUSE"
    assert action.item_key == "EXIST1"
    assert action.target_collections == ["OLD", "TARGET"]
    assert action.target_tags == ["Project::Existing", "Core10::01_HNF4A"]


def test_pdf_validator_rejects_html_disguised_as_pdf(tmp_path):
    path = tmp_path / "bad.pdf"
    path.write_bytes(b"<html>paywall</html>")
    assert validate_pdf(path) is False


def test_pmcid_provides_a_legal_pmc_pdf_candidate():
    assert legal_pdf_url(record(pmcid="PMC123456")) == "https://pmc.ncbi.nlm.nih.gov/articles/PMC123456/pdf/"


def test_created_item_key_accepts_zotero_success_string_response():
    assert _created_item_key({"successful": {"0": "ABCD1234"}}) == "ABCD1234"


def test_attachment_payload_sets_parent_item_and_pdf_fields():
    payload = build_attachment_payload("PARENT1", "paper.pdf")
    assert payload["parentItem"] == "PARENT1"
    assert payload["linkMode"] == "imported_file"
    assert payload["contentType"] == "application/pdf"
    assert payload["filename"] == "paper.pdf"
