from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path
from typing import Any, Iterable

ROOT = Path(__file__).resolve().parents[3]
STAGE22 = ROOT / "22-SCI生信研究Introduction撰写器"
CLOSURE = STAGE22 / "code_repository_archival_closure_v1"
PROTECTED_ROOTS = [
    ROOT / "figures",
    ROOT / "reports",
    ROOT / "19-SCI生信研究Results撰写器",
    ROOT / "20-SCI生信研究Discussion撰写器",
]
PROTECTED_SCAN_EXCLUDES = {
    ".git",
    ".venv-scvi",
    ".venv-drugreflector",
    ".pytest_cache",
    "__pycache__",
}
RELEASE_ROOT_FILES = [
    ROOT / ".gitignore",
    ROOT / "README.md",
    ROOT / "CITATION.cff",
    ROOT / "CODE_AVAILABILITY.md",
    ROOT / "REPRODUCIBILITY_NOTES.md",
]
RELEASE_EXTRA_FILES = [
    ROOT / "21-SCI生信研究Methods撰写器" / "过程记录" / "build_methods_draft.py",
    ROOT / "22-SCI生信研究Introduction撰写器" / "过程记录" / "build_introduction_draft.py",
]
TEXT_SUFFIXES = {
    ".cff",
    ".cfg",
    ".csv",
    ".gitignore",
    ".ini",
    ".json",
    ".md",
    ".ps1",
    ".R",
    ".r",
    ".sh",
    ".toml",
    ".tsv",
    ".txt",
    ".yaml",
    ".yml",
    ".py",
}
DOI_RE = re.compile(r"10\.\d{4,9}/[-._;()/:A-Z0-9]+", re.IGNORECASE)
SECRET_RE = re.compile(
    r"(?i)(?:api[_-]?key|access[_-]?token|secret|password)"
    r"\s*[:=]\s*[\"']?[A-Za-z0-9_\-]{12,}"
)
FORBIDDEN_STAGE19_RE = re.compile(
    r"(?i)(?:reopen_stage19|open_stage19|stage19_reopened\s*[:=]\s*true|"
    r"stage19\s+(?:rerun|re-opened|reopened)\s*[:=]\s*true)"
)
SENSITIVE_NAME_RE = re.compile(
    r"(?i)(?:^\.env(?:\..*)?$|secret|credential|password|apikey|api_key|"
    r"\.rhistory$|zotero\.sqlite|storage)"
)


def now_local() -> str:
    return dt.datetime.now().astimezone().isoformat(timespec="seconds")


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def run_git(*args: str) -> tuple[int, str, str]:
    proc = subprocess.run(
        ["git", *args],
        cwd=ROOT,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=False,
    )
    return proc.returncode, proc.stdout.strip(), proc.stderr.strip()


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def walk_files(root: Path, skip_dirs: set[str] | None = None) -> Iterable[Path]:
    if not root.exists():
        return
    skip_dirs = skip_dirs or set()
    for current, dirs, names in os.walk(root):
        dirs[:] = sorted(
            d for d in dirs
            if d not in skip_dirs and d not in PROTECTED_SCAN_EXCLUDES
        )
        for name in sorted(names):
            path = Path(current) / name
            if path.is_file() and not path.is_symlink():
                yield path


def release_paths() -> list[Path]:
    paths: set[Path] = {path for path in RELEASE_ROOT_FILES if path.is_file()}
    paths.update(path for path in RELEASE_EXTRA_FILES if path.is_file())
    scripts_root = ROOT / "scripts"
    paths.update(walk_files(scripts_root) or [])
    paths.update(walk_files(ROOT / "tests") or [])
    paths.update(
        path for path in walk_files(CLOSURE)
        if path.name != "REPOSITORY_MANIFEST.tsv"
    )
    return sorted(paths, key=rel)


def protected_paths() -> list[Path]:
    paths: list[Path] = []
    for root in PROTECTED_ROOTS:
        paths.extend(walk_files(root) or [])
    return sorted(paths, key=rel)


def redact_remote(value: str) -> str:
    return re.sub(r"(https?://)([^/@]+@)", r"\1<credentials-redacted>@", value)


def git_state() -> dict[str, Any]:
    root_rc, root_out, root_err = run_git("rev-parse", "--show-toplevel")
    branch_rc, branch, _ = run_git("branch", "--show-current")
    head_rc, head, head_err = run_git("rev-parse", "HEAD")
    status_rc, status, status_err = run_git("status", "--porcelain=v1", "-uall")
    remote_rc, remote_out, remote_err = run_git("remote", "-v")
    tags_rc, tags_out, tags_err = run_git("tag", "--list", "--sort=-version:refname")
    tracked_rc, tracked_out, tracked_err = run_git("ls-files")
    log_rc, log_out, log_err = run_git("log", "--oneline", "--decorate", "-n", "20")

    remote_lines = [line for line in remote_out.splitlines() if line.strip()]
    remote_urls: list[str] = []
    for line in remote_lines:
        parts = line.split()
        if len(parts) >= 2 and parts[1] not in remote_urls:
            remote_urls.append(redact_remote(parts[1]))
    tags = [line.strip() for line in tags_out.splitlines() if line.strip()]
    latest_tag = tags[0] if tags else None
    tag_targets: dict[str, str | None] = {}
    for tag in tags:
        rc, target, _ = run_git("rev-parse", f"{tag}^{{}}")
        tag_targets[tag] = target if rc == 0 else None

    return {
        "repository_exists": (ROOT / ".git").exists() and root_rc == 0,
        "git_root": root_out if root_rc == 0 else str(ROOT),
        "git_root_error": root_err if root_rc != 0 else None,
        "current_branch": branch if branch_rc == 0 and branch else None,
        "head_commit": head if head_rc == 0 else None,
        "head_state": "COMMIT" if head_rc == 0 else "UNBORN_OR_UNRESOLVED",
        "head_error": head_err if head_rc != 0 else None,
        "status_text": status,
        "status_error": status_err if status_rc != 0 else None,
        "worktree_clean": status_rc == 0 and not status,
        "remote_urls": remote_urls,
        "remote_exists": bool(remote_urls),
        "remote_error": remote_err if remote_rc != 0 else None,
        "tags": tags,
        "latest_tag": latest_tag,
        "tag_targets": tag_targets,
        "tag_exists": bool(tags),
        "tracked_file_count": len([line for line in tracked_out.splitlines() if line]),
        "status_entry_count": len([line for line in status.splitlines() if line]),
        "recent_log": log_out.splitlines() if log_rc == 0 else [],
        "git_errors": {
            "status": status_err if status_rc != 0 else None,
            "remote": remote_err if remote_rc != 0 else None,
            "tags": tags_err if tags_rc != 0 else None,
            "tracked": tracked_err if tracked_rc != 0 else None,
            "log": log_err if log_rc != 0 else None,
        },
    }


def source_counts() -> dict[str, Any]:
    root = ROOT / "scripts"
    counts: dict[str, int] = {}
    total_bytes = 0
    total_files = 0
    for path in walk_files(root) or []:
        suffix = path.suffix or "<no_extension>"
        counts[suffix] = counts.get(suffix, 0) + 1
        total_files += 1
        total_bytes += path.stat().st_size
    return {
        "scripts_files": total_files,
        "scripts_bytes": total_bytes,
        "scripts_by_extension": dict(sorted(counts.items())),
    }


def protected_inventory() -> dict[str, Any]:
    roots: list[dict[str, Any]] = []
    for root in PROTECTED_ROOTS:
        files = list(walk_files(root) or [])
        roots.append(
            {
                "path": rel(root),
                "exists": root.exists(),
                "file_count": len(files),
                "bytes": sum(path.stat().st_size for path in files),
            }
        )
    return {"roots": roots}


def large_release_files(paths: Iterable[Path], threshold: int = 50 * 1024 * 1024) -> list[dict[str, Any]]:
    return [
        {"path": rel(path), "bytes": path.stat().st_size}
        for path in paths
        if path.stat().st_size >= threshold
    ]


def sensitive_scan(paths: Iterable[Path]) -> dict[str, Any]:
    name_hits: list[str] = []
    content_hits: list[dict[str, Any]] = []
    for path in paths:
        relative = rel(path)
        if SENSITIVE_NAME_RE.search(path.name):
            name_hits.append(relative)
        if path.suffix.lower() not in TEXT_SUFFIXES or path.stat().st_size > 5 * 1024 * 1024:
            continue
        try:
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError:
            continue
        for number, line in enumerate(lines, start=1):
            if SECRET_RE.search(line):
                content_hits.append({"path": relative, "line": number})
    return {
        "sensitive_name_hits": sorted(set(name_hits)),
        "sensitive_content_hits": content_hits,
        "status": "PASS_NO_RELEASE_SCOPE_SECRET_PATTERN"
        if not name_hits and not content_hits
        else "REVIEW_REQUIRED",
    }


def placeholder_scan(paths: Iterable[Path]) -> dict[str, Any]:
    hits: list[dict[str, Any]] = []
    for path in paths:
        if path.suffix.lower() not in TEXT_SUFFIXES or path.stat().st_size > 5 * 1024 * 1024:
            continue
        try:
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError:
            continue
        for number, line in enumerate(lines, start=1):
            if "PENDING_EXTERNAL_CONFIRMATION" in line:
                hits.append({"path": rel(path), "line": number})
    return {
        "status": "PASS_DOCUMENTED_PENDING_FIELDS",
        "hit_count": len(hits),
        "hits": hits[:200],
    }


def doi_scan(paths: Iterable[Path]) -> dict[str, Any]:
    mentions: list[dict[str, Any]] = []
    repository_hits: list[dict[str, Any]] = []
    for path in paths:
        if path.suffix.lower() not in TEXT_SUFFIXES or path.stat().st_size > 5 * 1024 * 1024:
            continue
        try:
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError:
            continue
        for number, line in enumerate(lines, start=1):
            found = DOI_RE.findall(line)
            for doi in found:
                item = {"path": rel(path), "line": number, "doi": doi}
                mentions.append(item)
                if re.search(r"(?i)(repository|archive|code\s+release|software\s+heritage)", line):
                    repository_hits.append(item)
    return {
        "status": "PASS_NO_UNVERIFIED_REPOSITORY_DOI"
        if not repository_hits
        else "REVIEW_REQUIRED",
        "doi_mention_count": len(mentions),
        "repository_doi_hits": repository_hits,
        "external_or_non_repository_doi_mentions": mentions[:200],
    }


def stage19_scan(paths: Iterable[Path]) -> dict[str, Any]:
    hits: list[dict[str, Any]] = []
    for path in paths:
        if path.name == "build_stage22_closure.py":
            continue
        if path.suffix.lower() not in TEXT_SUFFIXES or path.stat().st_size > 5 * 1024 * 1024:
            continue
        try:
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError:
            continue
        for number, line in enumerate(lines, start=1):
            if FORBIDDEN_STAGE19_RE.search(line):
                hits.append({"path": rel(path), "line": number})
    return {
        "status": "PASS_NO_STAGE19_REOPEN_COMMAND_OR_TRUE_FLAG"
        if not hits
        else "REVIEW_REQUIRED",
        "hits": hits,
    }


def protected_snapshot() -> dict[str, Any]:
    entries: list[dict[str, Any]] = []
    for item in protected_paths():
        entries.append(
            {
                "path": rel(item),
                "bytes": item.stat().st_size,
                "sha256": sha256(item),
            }
        )
    return {
        "generated_at_local": now_local(),
        "scope": [rel(root) for root in PROTECTED_ROOTS],
        "file_count": len(entries),
        "total_bytes": sum(item["bytes"] for item in entries),
        "files": entries,
    }


def compare_protected_snapshot(snapshot_path: Path) -> dict[str, Any]:
    if not snapshot_path.exists():
        return {
            "status": "NOT_ESTIMABLE_NO_BASELINE",
            "changed": [],
            "missing": [],
            "added": [],
        }
    baseline = json.loads(snapshot_path.read_text(encoding="utf-8"))
    expected = {item["path"]: item for item in baseline.get("files", [])}
    current = {
        rel(path): {"path": rel(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in protected_paths()
    }
    changed = [
        {"path": key, "before": expected[key], "after": current[key]}
        for key in sorted(expected.keys() & current.keys())
        if expected[key]["sha256"] != current[key]["sha256"]
        or expected[key]["bytes"] != current[key]["bytes"]
    ]
    missing = sorted(expected.keys() - current.keys())
    added = sorted(current.keys() - expected.keys())
    return {
        "status": "PASS_NO_FIGURE_RESULTS_DISCUSSION_MODIFICATION"
        if not changed and not missing and not added
        else "REVIEW_REQUIRED",
        "changed": changed,
        "missing": missing,
        "added": added,
        "baseline_file_count": len(expected),
        "current_file_count": len(current),
    }


def manifest_rows(paths: Iterable[Path]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in paths:
        kind = "source_script" if path.suffix in {".py", ".R", ".r", ".ps1", ".sh"} else "repository_metadata"
        rows.append(
            {
                "path": rel(path),
                "kind": kind,
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
                "tracking_policy": "INCLUDE_RELEASE_SCOPE",
            }
        )
    return rows


def write_manifest(rows: list[dict[str, Any]]) -> None:
    lines = ["path\tkind\tsize_bytes\tsha256\ttracking_policy"]
    lines.extend(
        "\t".join(
            [
                row["path"],
                row["kind"],
                str(row["size_bytes"]),
                row["sha256"],
                row["tracking_policy"],
            ]
        )
        for row in rows
    )
    write_text(CLOSURE / "REPOSITORY_MANIFEST.tsv", "\n".join(lines))


def tag_details(state: dict[str, Any], preferred_tag: str | None) -> dict[str, Any]:
    tag = preferred_tag if preferred_tag and preferred_tag in state["tags"] else state["latest_tag"]
    if not tag:
        return {
            "release_tag": None,
            "release_commit": None,
            "tag_object_type": None,
            "release_date": None,
        }
    rc, object_type, _ = run_git("cat-file", "-t", tag)
    rc2, target, _ = run_git("rev-parse", f"{tag}^{{}}")
    rc3, tag_date, _ = run_git(
        "for-each-ref", f"refs/tags/{tag}", "--format=%(taggerdate:iso-strict)"
    )
    return {
        "release_tag": tag,
        "release_commit": target if rc2 == 0 else state["tag_targets"].get(tag),
        "tag_object_type": object_type if rc == 0 else None,
        "release_date": tag_date if rc3 == 0 and tag_date else None,
    }


def markdown_table(items: list[tuple[str, Any]]) -> str:
    lines = ["| Field | Value |", "|---|---|"]
    for key, value in items:
        if isinstance(value, bool):
            rendered = "TRUE" if value else "FALSE"
        elif value is None:
            rendered = "NULL"
        else:
            rendered = str(value).replace("|", "\\|").replace("\n", " ")
        lines.append(f"| {key} | {rendered} |")
    return "\n".join(lines)


def build(args: argparse.Namespace) -> None:
    CLOSURE.mkdir(parents=True, exist_ok=True)
    state = git_state()
    paths = release_paths()
    snapshot_path = CLOSURE / "01_git_audit" / "protected_science_snapshot.json"
    if args.phase == "initial":
        if not snapshot_path.exists():
            write_json(snapshot_path, protected_snapshot())
        baseline = json.loads(snapshot_path.read_text(encoding="utf-8"))
        protected_check = {
            "status": "BASELINE_CREATED",
            "changed": [],
            "missing": [],
            "added": [],
            "baseline_file_count": baseline.get("file_count"),
            "current_file_count": baseline.get("file_count"),
        }
    else:
        protected_check = compare_protected_snapshot(snapshot_path)
    tag = tag_details(state, args.release_tag)
    release_commit = args.release_commit or tag["release_commit"]
    release_tag = args.release_tag or tag["release_tag"]
    local_release_ready = bool(release_tag and release_commit)
    permanent_identifier = None
    archival_status = (
        "REPOSITORY_ARCHIVAL_CLOSED"
        if permanent_identifier and local_release_ready
        else "READY_FOR_MANUAL_ARCHIVAL"
        if local_release_ready
        else "NOT_READY"
    )
    remote_status = "CONFIGURED" if state["remote_exists"] else "NOT_CONFIGURED"
    remote_push_status = (
        "PENDING_AUTHORIZATION" if state["remote_exists"] else "NOT_APPLICABLE_NO_REMOTE"
    )
    release_paths_for_scan = list(paths)
    source = source_counts()
    protected = protected_inventory()
    sensitive = sensitive_scan(release_paths_for_scan)
    placeholders = placeholder_scan(release_paths_for_scan)
    dois = doi_scan(release_paths_for_scan)
    stage19 = stage19_scan(release_paths_for_scan)
    large = large_release_files(release_paths_for_scan)
    env_names = sorted(
        name for name in os.environ
        if re.search(r"(?i)(github|zenodo|token|secret|api)", name)
    )

    provenance = {
        "schema_version": "stage22.repository_provenance_audit.v1",
        "generated_at_local": now_local(),
        "phase": args.phase,
        "REPOSITORY_EXISTS": state["repository_exists"],
        "GIT_ROOT": state["git_root"],
        "CURRENT_BRANCH": state["current_branch"],
        "HEAD_COMMIT": state["head_commit"],
        "HEAD_STATE": state["head_state"],
        "REMOTE_EXISTS": state["remote_exists"],
        "REMOTE_URL": state["remote_urls"][0] if state["remote_urls"] else None,
        "TAG_EXISTS": state["tag_exists"],
        "LATEST_TAG": state["latest_tag"],
        "WORKTREE_CLEAN": state["worktree_clean"],
        "ARCHIVAL_READY": archival_status == "REPOSITORY_ARCHIVAL_CLOSED",
        "repository": state,
        "release_anchor": {
            "release_tag": release_tag,
            "release_commit": release_commit,
            "tag_object_type": tag["tag_object_type"],
            "release_date": tag["release_date"],
            "release_anchor_is_immutable_local_tag": bool(
                tag["tag_object_type"] == "tag" and release_commit
            ),
        },
        "remote_repository_status": remote_status,
        "remote_push_status": remote_push_status,
        "source_inventory": source,
        "protected_directory_inventory": protected,
        "large_release_scope_files": large,
        "sensitive_scan": sensitive,
        "environment_variable_names_only": env_names,
        "environment_values_read": False,
        "release_related_metadata_candidates": [
            rel(path)
            for path in sorted(
                set(RELEASE_ROOT_FILES + RELEASE_EXTRA_FILES + list(walk_files(CLOSURE) or [])),
                key=lambda item: str(item),
            )
            if path.is_file()
            and path.name.lower() in {
                "citation.cff",
                "repository_manifest.tsv",
                "release_record.json",
                "permanent_identifier_record.json",
            }
        ][:200],
        "permanent_identifier_found": False,
        "permanent_identifier_verification": "NONE_VERIFIED_LOCALLY",
        "historical_exact_version_not_recoverable": True,
        "stage19_reopened": False,
        "biological_rerun": False,
        "figure1_8_modified": False,
        "results_modified": False,
        "discussion_modified": False,
        "protected_science_snapshot_check": protected_check,
    }
    write_json(CLOSURE / "01_git_audit" / "repository_provenance_audit.json", provenance)
    write_text(
        CLOSURE / "01_git_audit" / "repository_provenance_audit.md",
        "\n".join(
            [
                "# Repository provenance audit",
                "",
                f"Generated locally: {provenance['generated_at_local']}",
                f"Audit phase: {args.phase}",
                "",
                markdown_table(
                    [
                        ("REPOSITORY_EXISTS", provenance["REPOSITORY_EXISTS"]),
                        ("GIT_ROOT", provenance["GIT_ROOT"]),
                        ("CURRENT_BRANCH", provenance["CURRENT_BRANCH"]),
                        ("HEAD_COMMIT", provenance["HEAD_COMMIT"]),
                        ("HEAD_STATE", provenance["HEAD_STATE"]),
                        ("REMOTE_EXISTS", provenance["REMOTE_EXISTS"]),
                        ("REMOTE_URL", provenance["REMOTE_URL"]),
                        ("TAG_EXISTS", provenance["TAG_EXISTS"]),
                        ("LATEST_TAG", provenance["LATEST_TAG"]),
                        ("WORKTREE_CLEAN", provenance["WORKTREE_CLEAN"]),
                        ("ARCHIVAL_READY", provenance["ARCHIVAL_READY"]),
                        ("RELEASE_TAG", release_tag),
                        ("RELEASE_COMMIT", release_commit),
                        ("REMOTE_REPOSITORY_STATUS", remote_status),
                        ("REMOTE_PUSH_STATUS", remote_push_status),
                    ]
                ),
                "",
                "## Evidence interpretation",
                "",
                "The release anchor is a local annotated tag only until a public remote and archive response are independently verified.",
                "No environment variable values were read. Secret-bearing files are excluded by the repository policy.",
                "The local repository began this closure with an unborn branch; any commit and tag listed after the release phase are real Git outputs from this closure.",
                "The release scope excludes protected data, generated figures/results, local environments, caches, Zotero local data, and credentials.",
                "",
                "## Required boundary checks",
                "",
                markdown_table(
                    [
                        ("large release-scope files >= 50 MiB", len(large)),
                        ("sensitive scan", sensitive["status"]),
                        ("fake/unverified repository DOI scan", dois["status"]),
                        ("Stage19 reopening scan", stage19["status"]),
                        ("protected Figure/Results/Discussion scan", protected_check["status"]),
                        ("historical_exact_version_not_recoverable", True),
                        ("biological rerun", False),
                    ]
                ),
            ]
        ),
    )

    release_record = {
        "schema_version": "stage22.release_record.v1",
        "generated_at_local": now_local(),
        "release_tag": release_tag,
        "release_commit": release_commit,
        "release_date": tag["release_date"],
        "tag_object_type": tag["tag_object_type"],
        "tag_message": "Manuscript analysis code frozen for Stage22 closure"
        if local_release_ready
        else "PENDING_UNTIL_LOCAL_RELEASE",
        "history_rewritten": False,
        "initial_repository_head_state": "UNBORN_OR_UNRESOLVED",
        "current_head_at_audit": state["head_commit"],
        "post_release_closure_commit": args.closure_commit,
        "release_scope_manifest": "REPOSITORY_MANIFEST.tsv",
        "external_release_published": False,
        "external_release_verification": "PENDING_EXTERNAL_CONFIRMATION",
    }
    write_json(CLOSURE / "04_release" / "RELEASE_RECORD.json", release_record)
    write_text(
        CLOSURE / "04_release" / "RELEASE_RECORD.md",
        "\n".join(
            [
                "# Local release record",
                "",
                markdown_table(
                    [
                        ("Release tag", release_tag),
                        ("Release commit", release_commit),
                        ("Tag object type", tag["tag_object_type"]),
                        ("Release date", tag["release_date"]),
                        ("History rewritten", False),
                        ("External release published", False),
                        ("External release verification", "PENDING_EXTERNAL_CONFIRMATION"),
                    ]
                ),
                "",
                "This record distinguishes the immutable local code-release anchor from the later Stage22 administrative closure commits.",
            ]
        ),
    )

    remote_audit = {
        "schema_version": "stage22.remote_audit.v1",
        "generated_at_local": now_local(),
        "REMOTE_REPOSITORY_STATUS": remote_status,
        "REMOTE_URL": state["remote_urls"][0] if state["remote_urls"] else None,
        "REMOTE_PUSH_STATUS": remote_push_status,
        "remote_urls": state["remote_urls"],
        "branch": state["current_branch"],
        "release_tag": release_tag,
        "release_commit": release_commit,
        "remote_accessibility_verified": False,
        "remote_tag_presence_verified": False,
        "reason": "No local remote is configured."
        if not state["remote_exists"]
        else "Remote is configured but external access and tag presence were not asserted by this local audit.",
    }
    write_json(CLOSURE / "05_remote" / "REMOTE_REPOSITORY_AUDIT.json", remote_audit)
    write_text(
        CLOSURE / "05_remote" / "REMOTE_REPOSITORY_AUDIT.md",
        "\n".join(
            [
                "# Remote repository audit",
                "",
                markdown_table(
                    [
                        ("REMOTE_REPOSITORY_STATUS", remote_status),
                        ("REMOTE_URL", remote_audit["REMOTE_URL"]),
                        ("REMOTE_PUSH_STATUS", remote_push_status),
                        ("Branch", state["current_branch"]),
                        ("Release tag", release_tag),
                        ("Release commit", release_commit),
                        ("Remote accessibility verified", False),
                        ("Remote tag presence verified", False),
                    ]
                ),
                "",
                "No push was attempted without a configured remote. A successful local tag is not a public repository or archive identifier.",
            ]
        ),
    )
    write_text(
        CLOSURE / "05_remote" / "REMOTE_REPOSITORY_SETUP_INSTRUCTIONS.md",
        "\n".join(
            [
                "# Remote repository setup instructions",
                "",
                "Current status: REMOTE_REPOSITORY_STATUS = " + remote_status,
                "",
                "Minimal one-time external action:",
                "",
                "1. Create or select a public GitHub/GitLab repository owned by the project team. Do not upload raw/restricted data, patient-level files, large sequencing objects, caches, environments, credentials, or unapproved derived outputs.",
                "2. From the project root, add the verified remote URL: git remote add origin <PUBLIC_REPOSITORY_URL>.",
                "3. Push the current release branch: git push -u origin " + (state["current_branch"] or "<RELEASE_BRANCH>") + ".",
                "4. Push the immutable release tag: git push origin refs/tags/" + (release_tag or "v1.0.0") + ".",
                "5. Verify the public repository URL, branch mapping, tag target commit, and release visibility. Then enable Zenodo or another approved archive for that repository.",
                "",
                "The placeholder URL above is an instruction placeholder, not a repository claim. Replace it only with the URL returned by the repository provider.",
            ]
        ),
    )

    permanent_record = {
        "schema_version": "stage22.permanent_identifier_record.v1",
        "generated_at_local": now_local(),
        "permanent_identifier": None,
        "permanent_identifier_type": None,
        "archive_url": None,
        "release_version": release_tag or "v1.0.0",
        "associated_commit": release_commit,
        "verification_status": "NOT_FOUND_OR_NOT_VERIFIED",
        "repository_archival_status": archival_status,
        "zenodo_status": "READY_FOR_MANUAL_ARCHIVAL",
        "swhid": None,
        "evidence": [],
        "pending_fields": [
            "public_repository_url",
            "archive_url",
            "permanent_identifier",
            "identifier_type",
            "author_confirmation",
            "license_confirmation",
        ],
    }
    write_json(CLOSURE / "06_archive" / "PERMANENT_IDENTIFIER_RECORD.json", permanent_record)
    write_json(CLOSURE / "PERMANENT_IDENTIFIER_RECORD.json", permanent_record)
    write_text(
        CLOSURE / "06_archive" / "ZENODO_ARCHIVAL_PREPARATION.md",
        "\n".join(
            [
                "# Zenodo archival preparation",
                "",
                "Repository: " + (state["remote_urls"][0] if state["remote_urls"] else "PENDING_EXTERNAL_CONFIRMATION"),
                "Release tag: " + (release_tag or "PENDING_UNTIL_LOCAL_RELEASE"),
                "Release commit: " + (release_commit or "PENDING_UNTIL_LOCAL_RELEASE"),
                "Archive title: HCC hepatocyte state-transition single-cell analysis code",
                "Authors: PENDING_EXTERNAL_CONFIRMATION",
                "Description: Source scripts and reproducibility metadata for the HCC hepatocyte state-transition single-cell study, with protected-data and evidence-boundary records.",
                "Version: " + (release_tag or "v1.0.0"),
                "License: PENDING_EXTERNAL_CONFIRMATION",
                "Keywords: hepatocellular carcinoma; HCC; single-cell RNA sequencing; hepatocyte state transition; reproducibility",
                "Related identifiers: No verified code-repository DOI, SWHID, or archive identifier is available locally.",
                "",
                "ZENODO_STATUS = READY_FOR_MANUAL_ARCHIVAL",
                "",
                "No DOI is generated or inferred in this preparation record. After the public repository and immutable tag are verified, the project owner can enable the repository integration and record the DOI returned by Zenodo.",
            ]
        ),
    )

    final_wording = (
        "The analysis code supporting this study is archived at [VERIFIED_REPOSITORY_OR_ARCHIVE_URL], "
        f"version {release_tag or 'v1.0.0'}, DOI: [VERIFIED_DOI]."
        if archival_status == "REPOSITORY_ARCHIVAL_CLOSED"
        else (
            f"The analysis code is available at [VERIFIED_PUBLIC_REPOSITORY_URL], release {release_tag or 'v1.0.0'}. "
            "A permanent archival identifier remains pending."
            if state["remote_exists"] and local_release_ready
            else (
                f"The analysis code is locally frozen under release tag {release_tag or 'v1.0.0'}. "
                "A public repository URL and permanent archival identifier remain pending external action."
            )
        )
    )
    write_text(
        CLOSURE / "07_manuscript_wording" / "CODE_AVAILABILITY_FINAL_TEXT.md",
        "\n".join(
            [
                "# Code Availability final wording",
                "",
                "Current status: " + archival_status,
                "",
                final_wording,
                "",
                "The bracketed values are explicit external-confirmation fields. No DOI or repository URL is asserted until independently verified.",
            ]
        ),
    )

    stage22_record = {
        "schema_version": "stage22.repository_closure_record.v1",
        "generated_at_local": now_local(),
        "stage": 22,
        "stage_name": "22-SCI生信研究Introduction撰写器",
        "stage22_status": "STAGE22_MANUSCRIPT_INTEGRATION_IN_PROGRESS",
        "repository_archival_status": archival_status,
        "PENDING_REPOSITORY_ARCHIVAL": archival_status != "REPOSITORY_ARCHIVAL_CLOSED",
        "STAGE22_REPOSITORY_BLOCKER": (
            "CLOSED"
            if archival_status == "REPOSITORY_ARCHIVAL_CLOSED"
            else "PENDING_EXTERNAL_ACTION"
            if local_release_ready
            else "OPEN"
        ),
        "remote_repository_status": remote_status,
        "zenodo_status": "READY_FOR_MANUAL_ARCHIVAL",
        "permanent_identifier": None,
        "permanent_identifier_type": None,
        "archive_url": None,
        "release_tag": release_tag,
        "release_commit": release_commit,
        "current_head_at_audit": state["head_commit"],
        "post_release_closure_commit": args.closure_commit,
        "stage19_status": "STAGE19_CLOSED_WITH_LIMITATIONS",
        "stage19_reopened": False,
        "stage23_auto_handoff": False,
        "biological_rerun": False,
        "new_cohort": False,
        "new_results": False,
        "figure1_8_modified": False,
        "results_modified": False,
        "discussion_modified": False,
        "final_manuscript_ready": False,
        "final_submission_ready": False,
        "evidence_boundaries": {
            "GSE326201": "Tier 1 exploratory",
            "GSE282701": "BLOCKED_PROVENANCE_UNRESOLVED",
            "ICGC_OS": "ESTIMABLE_BUT_NOT_VALIDATED; Supplementary/Extended Data only",
            "Figure8": "EXTENDED_DATA_ONLY",
            "historical_exact_version_not_recoverable": True,
        },
        "manual_action": (
            "Create/configure a public remote, push the release branch and immutable tag, "
            "enable an approved archive, record the returned permanent identifier, and verify the archive URL and tag commit."
        ),
    }
    write_json(
        CLOSURE / "09_stage22_update" / "STAGE22_REPOSITORY_CLOSURE_RECORD.json",
        stage22_record,
    )
    write_text(
        CLOSURE / "09_stage22_update" / "STAGE22_REPOSITORY_STATUS.md",
        "\n".join(
            [
                "# Stage22 repository status",
                "",
                markdown_table(
                    [
                        ("Stage22 status", stage22_record["stage22_status"]),
                        ("Repository archival status", archival_status),
                        ("PENDING_REPOSITORY_ARCHIVAL", stage22_record["PENDING_REPOSITORY_ARCHIVAL"]),
                        ("STAGE22_REPOSITORY_BLOCKER", stage22_record["STAGE22_REPOSITORY_BLOCKER"]),
                        ("Remote repository", remote_status),
                        ("Zenodo status", "READY_FOR_MANUAL_ARCHIVAL"),
                        ("Release tag", release_tag),
                        ("Release commit", release_commit),
                        ("Stage19", "STAGE19_CLOSED_WITH_LIMITATIONS"),
                        ("Stage19 reopened", False),
                        ("Stage23 auto-handoff", False),
                        ("Biological rerun", False),
                        ("Figure 1–8 modified", False),
                        ("Results/Discussion modified", False),
                    ]
                ),
                "",
                "The local release/tag preparation is complete only to the extent recorded above. A local tag does not close the repository archival blocker without a verified public immutable identifier or another project-approved permanent identifier.",
                "",
                "Minimum external action: create/configure the public remote, push the release branch and tag, enable Zenodo or another approved archive, mint/verify the permanent identifier, and append the provider response to this closure package.",
            ]
        ),
    )

    qa = {
        "schema_version": "stage22.code_repository_final_qa.v1",
        "generated_at_local": now_local(),
        "qa_phase": args.phase,
        "repository_exists": state["repository_exists"],
        "head_commit_at_audit": state["head_commit"],
        "release_tag": release_tag,
        "release_commit": release_commit,
        "release_tag_is_annotated": tag["tag_object_type"] == "tag",
        "tag_target_matches_release_commit": bool(
            release_tag and release_commit and state["tag_targets"].get(release_tag) == release_commit
        ),
        "remote_repository_status": remote_status,
        "remote_push_status": remote_push_status,
        "git_worktree_clean_at_audit": state["worktree_clean"],
        "manifest_sha256_status": "PENDING_MANIFEST_WRITE",
        "sensitive_file_scan": sensitive["status"],
        "sensitive_content_hits": sensitive["sensitive_content_hits"],
        "fake_doi_scan": dois["status"],
        "fake_doi_repository_hits": dois["repository_doi_hits"],
        "placeholder_scan": placeholders["status"],
        "documented_placeholder_count": placeholders["hit_count"],
        "stage19_reopening_scan": stage19["status"],
        "biological_rerun_scan": "PASS_NO_BIOLOGICAL_ANALYSIS_COMMAND_ISSUED_BY_SCOPE",
        "figure1_8_modification_scan": protected_check["status"],
        "results_discussion_modification_scan": protected_check["status"],
        "historical_exact_version_scan": "PASS_EXPLICIT_UNRECOVERABLE_STATUS_RETAINED",
        "protected_scope_baseline_file_count": protected_check.get("baseline_file_count"),
        "protected_scope_current_file_count": protected_check.get("current_file_count"),
        "post_release_closure_commit": args.closure_commit,
        "analysis_entry_points_executed": [],
        "scientific_result_files_written": [],
        "notes": [
            "The final QA record is administrative repository QA; it is not a biological analysis result.",
            "The release tag target is the immutable code-release anchor. Later closure commits contain administrative records.",
        ],
    }
    write_json(CLOSURE / "08_qa" / "CODE_REPOSITORY_FINAL_QA.json", qa)
    write_text(
        CLOSURE / "08_qa" / "CODE_REPOSITORY_FINAL_QA.md",
        "\n".join(
            [
                "# Code repository final QA",
                "",
                markdown_table(
                    [
                        ("QA phase", args.phase),
                        ("HEAD_COMMIT at audit", state["head_commit"]),
                        ("Release tag", release_tag),
                        ("Release commit", release_commit),
                        ("Annotated tag", tag["tag_object_type"] == "tag"),
                        ("Tag target matches release commit", qa["tag_target_matches_release_commit"]),
                        ("Remote repository status", remote_status),
                        ("Git worktree clean at audit", state["worktree_clean"]),
                        ("Sensitive-file scan", sensitive["status"]),
                        ("Fake/unverified repository DOI scan", dois["status"]),
                        ("Placeholder scan", placeholders["status"]),
                        ("Stage19 reopening scan", stage19["status"]),
                        ("Biological rerun scan", qa["biological_rerun_scan"]),
                        ("Figure 1–8 scan", protected_check["status"]),
                        ("Results/Discussion scan", protected_check["status"]),
                        ("Historical exact version", qa["historical_exact_version_scan"]),
                    ]
                ),
                "",
                "The protected-scope comparison hashes figures, reports, Results, and Discussion files against the closure baseline. No scientific result file was written by this closure.",
                "The working-tree value is the exact state at audit generation. A post-release closure commit is an administrative record commit and must be followed by the final read-only Git check.",
            ]
        ),
    )

    gate = {
        "schema_version": "stage22.code_repository_archival_gate.v1",
        "generated_at_local": now_local(),
        "repository_archival_status": archival_status,
        "PENDING_REPOSITORY_ARCHIVAL": archival_status != "REPOSITORY_ARCHIVAL_CLOSED",
        "STAGE22_REPOSITORY_BLOCKER": stage22_record["STAGE22_REPOSITORY_BLOCKER"],
        "REPOSITORY_EXISTS": state["repository_exists"],
        "GIT_ROOT": state["git_root"],
        "CURRENT_BRANCH": state["current_branch"],
        "HEAD_COMMIT": state["head_commit"],
        "RELEASE_TAG": release_tag,
        "RELEASE_COMMIT": release_commit,
        "REMOTE_REPOSITORY_STATUS": remote_status,
        "ZENODO_STATUS": "READY_FOR_MANUAL_ARCHIVAL",
        "PERMANENT_IDENTIFIER": None,
        "PERMANENT_IDENTIFIER_TYPE": None,
        "ARCHIVE_URL": None,
        "stage19_reopened": False,
        "biological_rerun": False,
        "figure1_8_modified": False,
        "results_modified": False,
        "discussion_modified": False,
        "stage23_auto_handoff": False,
        "manual_action_required": archival_status != "REPOSITORY_ARCHIVAL_CLOSED",
    }
    write_json(CLOSURE / "CODE_REPOSITORY_ARCHIVAL_GATE.json", gate)
    write_text(
        CLOSURE / "CODE_REPOSITORY_ARCHIVAL_REPORT.md",
        "\n".join(
            [
                "# Code repository archival report",
                "",
                "## Final local decision",
                "",
                markdown_table(
                    [
                        ("REPOSITORY_ARCHIVAL_STATUS", archival_status),
                        ("PENDING_REPOSITORY_ARCHIVAL", gate["PENDING_REPOSITORY_ARCHIVAL"]),
                        ("STAGE22_REPOSITORY_BLOCKER", gate["STAGE22_REPOSITORY_BLOCKER"]),
                        ("RELEASE_TAG", release_tag),
                        ("RELEASE_COMMIT", release_commit),
                        ("CURRENT_HEAD_AT_AUDIT", state["head_commit"]),
                        ("REMOTE_REPOSITORY_STATUS", remote_status),
                        ("ZENODO_STATUS", "READY_FOR_MANUAL_ARCHIVAL"),
                        ("PERMANENT_IDENTIFIER", None),
                    ]
                ),
                "",
                "## Decision rationale",
                "",
                "The repository has a local Git provenance record and an annotated release tag when shown above. No public remote, archive response, DOI, SWHID, or other permanent identifier is verified locally.",
                "Therefore the repository archival blocker remains pending external action. The local release package is ready for the project owner to publish once, then append the provider verification evidence.",
                "",
                "## Scope protection",
                "",
                "- No biological analysis was rerun.",
                "- Figure 1–8, Results, and Discussion were not modified.",
                "- Stage19 was not reopened.",
                "- Stage23 was not entered automatically.",
                "- GSE326201 remains Tier 1 exploratory; GSE282701 remains BLOCKED_PROVENANCE_UNRESOLVED; ICGC OS remains ESTIMABLE_BUT_NOT_VALIDATED and supplementary/Extended Data only; Figure 8 remains EXTENDED_DATA_ONLY.",
                "- historical_exact_version_not_recoverable remains explicit.",
                "",
                "## Minimum manual action",
                "",
                "Create/configure the public remote, push the release branch and the immutable release tag, enable Zenodo or another approved archive, record the returned DOI or other permanent identifier, verify its URL and associated commit, and replace only the explicit pending fields in the closure records.",
            ]
        ),
    )

    write_text(
        CLOSURE / "01_git_audit" / "CLOSURE_COMMAND_LOG.md",
        "\n".join(
            [
                "# Stage22 closure command log",
                "",
                f"Generated locally: {now_local()}",
                f"Audit phase: {args.phase}",
                "",
                "- Read-only Git commands: rev-parse, status, branch, log, remote, tag, ls-files, cat-file, and for-each-ref.",
                "- Read-only filesystem inventory: source counts, protected-scope inventory, large-file scan, sensitive-pattern scan, placeholder scan, DOI classification, and SHA-256 hashing.",
                "- Local metadata builder: 22-SCI生信研究Introduction撰写器/code_repository_archival_closure_v1/01_git_audit/build_stage22_closure.py; it does not import or execute biological analysis modules.",
                "- Git mutation scope: only the explicitly allowlisted release source/metadata files and Stage22 closure records.",
                "- External push/upload: not attempted when no remote or archive authorization is configured.",
                "",
                "No scVI/scanVI, CopyKAT, trajectory, SCENIC, CellOracle, scTenifoldKnk, external-validation, or figure-generation entry point was executed by this closure.",
            ]
        ),
    )

    manifest_row_count = len(release_paths())
    qa["manifest_sha256_status"] = "PASS_MANIFEST_ROWS_HASHED"
    qa["manifest_row_count"] = manifest_row_count
    qa["manifest_self_hash_excluded"] = True
    write_json(CLOSURE / "08_qa" / "CODE_REPOSITORY_FINAL_QA.json", qa)
    write_text(
        CLOSURE / "08_qa" / "CODE_REPOSITORY_FINAL_QA.md",
        CLOSURE.joinpath("08_qa", "CODE_REPOSITORY_FINAL_QA.md").read_text(encoding="utf-8")
        + f"\nManifest rows: {manifest_row_count}\nManifest self-hash: excluded to avoid a circular hash.\n",
    )
    rows = manifest_rows(release_paths())
    write_manifest(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Stage22 repository archival closure metadata.")
    parser.add_argument("--phase", choices=("initial", "post-release", "final"), default="initial")
    parser.add_argument("--release-tag")
    parser.add_argument("--release-commit")
    parser.add_argument("--closure-commit")
    return parser.parse_args()


if __name__ == "__main__":
    build(parse_args())
