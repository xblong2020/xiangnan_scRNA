import json
from pathlib import Path

import pytest

import zotero_hcc_classifier.classifier as classifier_module
from zotero_hcc_classifier.classifier import classify_item, classify_items, is_hcc_candidate, run_qa
from zotero_hcc_classifier.config import load_collection_spec, load_config
from zotero_hcc_classifier.dry_run import run_dry_run
from zotero_hcc_classifier.apply_classification import apply_rows
from zotero_hcc_classifier.rollback import rollback_changes
from zotero_hcc_classifier.zotero_client import (
    ZoteroClient,
    ZoteroAPIError,
    build_collection_create_payload,
    build_membership_patch,
    build_write_headers,
    merge_collections,
)


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "zotero_hcc_classifier" / "config.yaml"
RULES_PATH = ROOT / "zotero_hcc_classifier" / "classification_rules.yaml"


def item(title, abstract="", tags=None, extra="", publication=""):
    return {
        "key": "ITEM1",
        "version": 7,
        "data": {
            "itemType": "journalArticle",
            "title": title,
            "abstractNote": abstract,
            "tags": [{"tag": tag} for tag in (tags or [])],
            "extra": extra,
            "publicationTitle": publication,
            "date": "2024",
            "DOI": "fixture-doi-sentinel",
            "collections": ["ORIGINAL"],
        },
    }


def test_rules_define_exactly_one_root_eight_modules_and_thirty_two_leaf_collections():
    spec = load_collection_spec(RULES_PATH)
    assert spec.root_name == "HCC_Hepatocyte-State-Plasticity"
    assert len(spec.modules) == 8
    assert sum(len(module.leaves) for module in spec.modules) == 32


def test_config_forbids_metadata_tag_and_delete_mutations():
    config = load_config(CONFIG_PATH)
    assert config["classification"]["modify_tags"] is False
    assert config["classification"]["modify_metadata"] is False
    assert config["classification"]["delete_items"] is False
    assert config["classification"]["delete_collections"] is False


def test_write_headers_include_server_id_and_local_api_key():
    headers = build_write_headers("server-123", "secret-key")
    assert headers["Zotero-Server-ID"] == "server-123"
    assert headers["Zotero-API-Key"] == "secret-key"


def test_localhost_transport_uses_ipv4_loopback_for_zotero():
    client = ZoteroClient("http://localhost:23119/api")
    assert client.base_url == "http://127.0.0.1:23119/api"


def test_merge_collections_preserves_existing_membership_and_deduplicates():
    assert merge_collections(["A", "B"], ["B", "C"]) == ["A", "B", "C"]


def test_patch_payload_contains_only_collections():
    assert build_membership_patch(["A", "C"]) == {"collections": ["A", "C"]}


def test_collection_create_payload_is_v3_object_array():
    assert build_collection_create_payload("Root", False) == [{"name": "Root", "parentCollection": False}]


def test_hcc_anchor_requires_hcc_context_for_generic_liver_disease():
    assert is_hcc_candidate(item("NAFLD in chronic liver disease", "hepatitis and cirrhosis")) is False
    assert is_hcc_candidate(item("HCC single-cell atlas", "hepatocellular carcinoma")) is True
    assert is_hcc_candidate(item("HCC metabolism", "SOX4 was cited as background only")) is True
    assert is_hcc_candidate(item("SHCC1 metabolism", "generic liver disease")) is False


def test_reference_only_sox4_mention_does_not_trigger_sox4_leaf():
    result = classify_item(item("HCC metabolism", "SOX4 was cited as background only"))
    assert "04_Malignant-Plasticity/01_SOX4" not in result.recommended_paths


def test_one_paper_can_receive_multiple_method_and_validation_labels():
    result = classify_item(
        item(
            "SOX4 scRNA-seq inferCNV SCENIC TCGA xenograft in HCC",
            "Primary results show malignant plasticity and an orthotopic xenograft model.",
        )
    )
    assert set(result.recommended_paths) >= {
        "01_Disease-State/04_Single-cell-atlas",
        "04_Malignant-Plasticity/01_SOX4",
        "05_Malignant-Evolution/01_CNV-aneuploidy",
        "06_Regulatory-Network/01_SCENIC-pySCENIC",
        "07_External-Validation/01_TCGA-LIHC",
        "08_Experimental-Validation/04_Xenograft-in-vivo",
    }


def test_batch_classification_loads_rules_once(monkeypatch):
    calls = 0
    original = classifier_module.load_collection_spec

    def counted(path):
        nonlocal calls
        calls += 1
        return original(path)

    monkeypatch.setattr(classifier_module, "load_collection_spec", counted)
    classify_items([item("HCC scRNA-seq", "hepatocellular carcinoma") for _ in range(3)])
    assert calls == 1


def test_qa_counts_reviewable_paths_and_ignores_low_confidence_suggestions():
    low_rows = [
        {"ItemKey": "ITEM1", "RecommendedPath": f"M/P{i}", "Confidence": 0.68, "EvidenceTerms": ""}
        for i in range(9)
    ]
    review_rows = [
        {"ItemKey": "ITEM2", "RecommendedPath": f"M/P{i}", "Confidence": 0.70, "EvidenceTerms": ""}
        for i in range(9)
    ]
    low_result = run_qa(low_rows)
    review_result = run_qa(review_rows)
    assert low_result["stop_formal_apply"] is False
    assert review_result["stop_formal_apply"] is True


def test_dry_run_writes_required_snapshot_preview_and_ambiguous_artifacts(tmp_path):
    result = run_dry_run(FakeClient(), output_root=tmp_path)
    assert result["candidate_count"] == 1
    assert (tmp_path / "data" / "collections_before.json").exists()
    assert (tmp_path / "data" / "items_before.json").exists()
    assert (tmp_path / "data" / "classification_preview.csv").exists()
    assert (tmp_path / "reports" / "CLASSIFICATION_PREVIEW.md").exists()


def test_apply_rejects_medium_confidence_and_preserves_existing_collections(tmp_path):
    client = FakeClient()
    result = apply_rows(
        client,
        rows=[
            {
                "ItemKey": "HIGH1",
                "RecommendedPath": "01_Disease-State/04_Single-cell-atlas",
                "Confidence": 0.9,
                "Action": "APPLY",
            },
            {
                "ItemKey": "MED1",
                "RecommendedPath": "04_Malignant-Plasticity/01_SOX4",
                "Confidence": 0.75,
                "Action": "REVIEW",
            },
        ],
        output_root=tmp_path,
        collection_keys={"01_Disease-State/04_Single-cell-atlas": "NEW1"},
    )
    assert result.modified_item_keys == ["HIGH1"]
    assert result.skipped_item_keys == ["MED1"]
    assert client.patched["HIGH1"] == ["ORIGINAL", "NEW1"]


def test_apply_reloads_item_and_retries_additive_patch_once_on_412(tmp_path):
    client = ConflictClient()
    result = apply_rows(
        client,
        rows=[{"ItemKey": "HIGH1", "RecommendedPath": "01_Disease-State/04_Single-cell-atlas", "Confidence": 0.9, "Action": "APPLY"}],
        output_root=tmp_path,
        collection_keys={"01_Disease-State/04_Single-cell-atlas": "NEW1"},
    )
    assert result.modified_item_keys == ["HIGH1"]
    assert result.errors == []
    assert client.patch_attempts == 2
    assert client.patched == ["ORIGINAL", "NEW1"]


def test_rollback_targets_only_recorded_modified_items(tmp_path):
    (tmp_path / "data").mkdir(parents=True)
    (tmp_path / "data" / "changes_applied.json").write_text(json.dumps({
        "modified_items": [{"item_key": "A", "original_collections": ["ORIGINAL"], "added_collections": ["NEW"]}],
    }), encoding="utf-8")
    client = RollbackClient()
    result = rollback_changes(client, output_root=tmp_path)
    assert result["restored_item_keys"] == ["A"]
    assert result["errors"] == []
    assert client.patched == {"A": ["ORIGINAL"]}
    assert client.get_calls == ["A"]


class FakeClient:
    server_id = "server-123"
    api_version = "3"
    schema_version = "44"
    library_version = 101

    def __init__(self):
        self.patched = {}
        self.items = [
            item("HCC single-cell atlas", "hepatocellular carcinoma scRNA-seq"),
        ]
        self.items[0]["key"] = "ITEM1"
        self.collections = [
            {"key": "ORIGINAL", "version": 3, "data": {"name": "Existing", "parentCollection": False}}
        ]

    def get_environment(self):
        return {
            "reachable": True,
            "api_version": self.api_version,
            "schema_version": self.schema_version,
            "server_id": self.server_id,
        }

    def get_collections(self):
        return self.collections

    def get_items_top(self):
        return self.items

    def patch_item_collections(self, item_key, collections, version, api_key=None):
        self.patched[item_key] = list(collections)
        return {"key": item_key, "version": version + 1, "data": {"collections": list(collections)}}


class ConflictClient:
    def __init__(self):
        self.item = {"key": "HIGH1", "version": 7, "data": {"collections": ["ORIGINAL"]}}
        self.patch_attempts = 0
        self.patched = None

    def get_item(self, item_key):
        return {"key": self.item["key"], "version": self.item["version"], "data": {"collections": list(self.item["data"]["collections"])}}

    def patch_item_collections(self, item_key, collections, version, api_key=None):
        self.patch_attempts += 1
        if self.patch_attempts == 1:
            self.item["version"] += 1
            raise ZoteroAPIError(412, "version changed")
        self.patched = list(collections)
        self.item["data"]["collections"] = list(collections)
        self.item["version"] += 1
        return {"key": item_key, "version": self.item["version"]}


class RollbackClient:
    def __init__(self):
        self.patched = {}
        self.get_calls = []

    def get_item(self, item_key):
        self.get_calls.append(item_key)
        return {"key": item_key, "version": 3, "data": {"collections": ["ORIGINAL", "NEW"]}}

    def patch_item_collections(self, item_key, collections, version, api_key=None):
        self.patched[item_key] = list(collections)
        return {"key": item_key, "version": version + 1}
