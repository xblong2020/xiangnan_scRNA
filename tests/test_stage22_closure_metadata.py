from __future__ import annotations

import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BUILDER = (
    ROOT
    / "22-SCI生信研究Introduction撰写器"
    / "code_repository_archival_closure_v1"
    / "01_git_audit"
    / "build_stage22_closure.py"
)

spec = importlib.util.spec_from_file_location("stage22_closure_builder", BUILDER)
assert spec and spec.loader
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)


def test_publication_state_preserves_github_release_and_zenodo_pending():
    state = builder.resolve_publication_state(
        release_ready=True,
        github_release_status="PUBLISHED",
        zenodo_status="MANUAL_ACTION_REQUIRED",
        permanent_identifier=None,
    )

    assert state["repository_archival_status"] == "OPEN_PENDING_EXTERNAL_ARCHIVAL"
    assert state["stage22_repository_blocker"] == "PENDING_ZENODO_ARCHIVAL"
    assert state["final_gate"] == "MANUAL_ZENODO_ACTION_REQUIRED"
    assert state["permanent_identifier"] is None


def test_public_audit_root_is_redacted():
    assert builder.public_audit_root(
        r"C:UsersAdministratorOneDrive文档湘南学院单细胞"
    ) == "REDACTED_LOCAL_PATH"


def test_closure_json_keys_are_unique_case_insensitively():
    closure = ROOT / "22-SCI生信研究Introduction撰写器" / "code_repository_archival_closure_v1"

    def hook(pairs):
        seen = set()
        result = {}
        for key, value in pairs:
            normalized = key.lower()
            assert normalized not in seen, key
            seen.add(normalized)
            result[key] = value
        return result

    for path in closure.rglob("*.json"):
        json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=hook)
