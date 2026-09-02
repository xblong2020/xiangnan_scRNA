"""Export visible Codex project conversations as ChatGPT-ready Markdown files.

The exporter reads only user-facing ``event_msg`` records. System prompts,
developer instructions, tool calls, tool outputs, and hidden reasoning remain
outside the generated discussion documents.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


PROJECT_ROOT = Path(__file__).resolve().parents[1]
VISIBLE_EVENT_TYPES = {"user_message": "user", "agent_message": "assistant"}
SECRET_PATTERNS = (
    re.compile(r"\bsk-[A-Za-z0-9_-]{20,}\b"),
    re.compile(r"\bgh[pousr]_[A-Za-z0-9]{20,}\b"),
    re.compile(r"\bxox[baprs]-[A-Za-z0-9-]{20,}\b"),
    re.compile(r"\bAIza[0-9A-Za-z_-]{20,}\b"),
    re.compile(r"(?i)\bBearer\s+[A-Za-z0-9._~-]{12,}\b"),
    re.compile(
        r"(?im)(\b(?:api[_ -]?key|access[_ -]?token|refresh[_ -]?token|"
        r"authorization|password|passwd|secret)\b\s*[:=]\s*)([^\s,;]+)"
    ),
)
WINDOWS_ABSOLUTE_PATH_PATTERN = re.compile(
    r"(?i)(?<![A-Z0-9])[A-Z]:(?:\\\\+|\\|/(?!/))[^\r\n`\"'<>|()\[\]{}；。！？，、;]*"
)


def normalise_path(value: str | Path) -> str:
    """Return a platform-normalized absolute path for workspace matching."""
    return os.path.normcase(os.path.normpath(str(Path(value).resolve())))


def path_is_project_scoped(cwd: str | None, project_root: Path) -> bool:
    if not cwd:
        return False
    project = normalise_path(project_root)
    candidate = normalise_path(cwd)
    return candidate == project or candidate.startswith(project + os.sep)


def read_json_lines(path: Path) -> tuple[list[dict[str, Any]], int]:
    """Read a JSONL file while tolerating a partial trailing live-session line."""
    records: list[dict[str, Any]] = []
    parse_errors = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                parse_errors += 1
                continue
            if isinstance(record, dict):
                records.append(record)
    return records, parse_errors


def load_thread_titles(session_index_path: Path) -> dict[str, str]:
    """Choose the latest title recorded for every Codex thread ID."""
    if not session_index_path.is_file():
        return {}

    titles: dict[str, tuple[str, str]] = {}
    records, _ = read_json_lines(session_index_path)
    for record in records:
        session_id = record.get("id")
        title = record.get("thread_name")
        updated_at = str(record.get("updated_at", ""))
        if not isinstance(session_id, str) or not isinstance(title, str) or not title.strip():
            continue
        if session_id not in titles or updated_at >= titles[session_id][0]:
            titles[session_id] = (updated_at, title.strip())
    return {session_id: item[1] for session_id, item in titles.items()}


def discover_project_sessions(codex_home: Path, project_root: Path) -> list[dict[str, Any]]:
    """Locate unique Codex JSONL sessions whose metadata uses the project cwd."""
    candidates: list[Path] = []
    for relative_root in ("sessions", "archived_sessions"):
        root = codex_home / relative_root
        if root.is_dir():
            candidates.extend(root.rglob("*.jsonl"))

    found: dict[str, dict[str, Any]] = {}
    for source_path in sorted(candidates, key=lambda item: (item.stat().st_mtime, str(item)), reverse=True):
        records, parse_errors = read_json_lines(source_path)
        metadata = [
            record.get("payload", {})
            for record in records
            if record.get("type") == "session_meta" and isinstance(record.get("payload"), dict)
        ]
        if not any(path_is_project_scoped(item.get("cwd"), project_root) for item in metadata):
            continue

        session_id = next(
            (
                item.get("id")
                for item in metadata
                if isinstance(item.get("id"), str) and item.get("id")
            ),
            source_path.stem,
        )
        if session_id in found:
            continue
        found[session_id] = {
            "session_id": session_id,
            "source_path": source_path,
            "records": records,
            "parse_errors": parse_errors,
            "metadata": metadata[0] if metadata else {},
        }

    return sorted(
        found.values(),
        key=lambda item: (
            str(item["metadata"].get("timestamp", "")),
            str(item["source_path"]),
        ),
    )


def extract_visible_messages(records: Iterable[dict[str, Any]]) -> list[dict[str, str]]:
    """Extract frontend-visible user/assistant message events in source order."""
    messages: list[dict[str, str]] = []
    for record in records:
        if record.get("type") != "event_msg":
            continue
        payload = record.get("payload")
        if not isinstance(payload, dict):
            continue
        role = VISIBLE_EVENT_TYPES.get(payload.get("type"))
        message = payload.get("message")
        if role is None or not isinstance(message, str) or not message.strip():
            continue
        messages.append(
            {
                "role": role,
                "timestamp": str(record.get("timestamp", "")),
                "phase": str(payload.get("phase", "")),
                "message": message.strip(),
            }
        )
    return messages


def redact_for_external_discussion(text: str, project_root: Path) -> tuple[str, int]:
    """Remove credentials and machine-specific absolute paths from visible text."""
    redactions = 0
    cleaned = text.replace("\r\n", "\n").replace("\r", "\n")
    project_variants = {str(project_root), project_root.as_posix()}
    for variant in sorted(project_variants, key=len, reverse=True):
        if not variant:
            continue
        pattern = re.compile(re.escape(variant), flags=re.IGNORECASE)
        cleaned, count = pattern.subn("$PROJECT_ROOT", cleaned)
        redactions += count

    cleaned, count = WINDOWS_ABSOLUTE_PATH_PATTERN.subn("<LOCAL_PATH>", cleaned)
    redactions += count
    for pattern in SECRET_PATTERNS:
        if pattern.groups:
            cleaned, count = pattern.subn(lambda match: f"{match.group(1)}[REDACTED]", cleaned)
        else:
            cleaned, count = pattern.subn("[REDACTED]", cleaned)
        redactions += count
    return cleaned, redactions


def fallback_title(messages: list[dict[str, str]], session_id: str) -> str:
    for message in messages:
        if message["role"] != "user":
            continue
        first_line = next((line.strip() for line in message["message"].splitlines() if line.strip()), "")
        first_line = re.sub(r"^#+\s*", "", first_line)
        if first_line:
            return first_line[:80]
    return f"Project session {session_id[:8]}"


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def markdown_timestamp(value: str) -> str:
    return value.replace("T", " ").replace("Z", " UTC") if value else "时间未记录"


def render_session_markdown(
    *,
    title: str,
    session_id: str,
    messages: list[dict[str, str]],
    source_created_at: str,
    redaction_count: int,
) -> str:
    user_count = sum(message["role"] == "user" for message in messages)
    assistant_count = sum(message["role"] == "assistant" for message in messages)
    lines = [
        f"# {title}",
        "",
        "> 本文为可供 ChatGPT 继续讨论的项目会话转录。仅保留用户与助手可见消息；系统提示、开发者指令、工具调用、工具输出和内部推理均已排除。",
        "> 已自动掩盖凭据与本机绝对路径；`$PROJECT_ROOT` 表示本项目根目录。",
        "",
        "## 会话信息",
        "",
        f"- 会话 ID：`{session_id}`",
        f"- 开始时间：{markdown_timestamp(source_created_at)}",
        f"- 可见消息：用户 {user_count} 条；助手 {assistant_count} 条",
        f"- 自动掩盖：{redaction_count} 处",
        "",
        "## 建议讨论提示",
        "",
        "请基于以下会话记录继续讨论本项目。先区分已完成证据、待验证假设和未决问题，再给出可执行的下一步；涉及生物学结论时，请明确证据来源与限制。",
        "",
        "## 对话记录",
        "",
    ]
    for message in messages:
        if message["role"] == "user":
            role = "用户"
        elif message.get("phase") == "final":
            role = "助手（答复）"
        else:
            role = "助手（进度）"
        lines.extend(
            [
                f"### {role}｜{markdown_timestamp(message['timestamp'])}",
                "",
                message["message"],
                "",
                "---",
                "",
            ]
        )
    return "\n".join(lines).rstrip() + "\n"


def render_index_markdown(sessions: list[dict[str, Any]]) -> str:
    lines = [
        "# 项目会话：ChatGPT 讨论索引",
        "",
        "> 每份 Markdown 是一个本项目 Codex 会话的可见对话转录，并已排除系统、工具和内部推理内容。",
        "",
        "| 序号 | 开始日期 | 会话标题 | 用户/助手消息 | 文档 |",
        "| ---: | --- | --- | ---: | --- |",
    ]
    for index, session in enumerate(sessions, start=1):
        title = str(session["title"]).replace("|", "\\|")
        date = str(session["source_created_at"])[:10] or "未记录"
        message_counts = f"{session['user_message_count']}/{session['assistant_message_count']}"
        lines.append(
            f"| {index} | {date} | {title} | {message_counts} | [{session['output_file']}]({session['output_file']}) |"
        )
    lines.extend(
        [
            "",
            "使用每个会话文档中的“建议讨论提示”作为 ChatGPT 的首条任务说明，并将该文档作为上下文输入。",
            "",
        ]
    )
    return "\n".join(lines)


def ensure_output_directory(output_dir: Path) -> None:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(
            f"Output directory already contains files: {output_dir}. Choose a new output directory."
        )
    output_dir.mkdir(parents=True, exist_ok=True)


def export_project_sessions(
    *,
    codex_home: Path,
    project_root: Path,
    output_dir: Path,
) -> dict[str, Any]:
    """Export every Codex session scoped to ``project_root`` into individual MD files."""
    codex_home = Path(codex_home)
    project_root = Path(project_root).resolve()
    output_dir = Path(output_dir)
    ensure_output_directory(output_dir)

    titles = load_thread_titles(codex_home / "session_index.jsonl")
    discovered = discover_project_sessions(codex_home, project_root)
    exported: list[dict[str, Any]] = []
    for position, session in enumerate(discovered, start=1):
        raw_messages = extract_visible_messages(session["records"])
        redaction_count = 0
        messages: list[dict[str, str]] = []
        for message in raw_messages:
            cleaned_message, count = redact_for_external_discussion(message["message"], project_root)
            redaction_count += count
            messages.append({**message, "message": cleaned_message})

        session_id = str(session["session_id"])
        title = titles.get(session_id) or fallback_title(messages, session_id)
        title, title_redaction_count = redact_for_external_discussion(title, project_root)
        redaction_count += title_redaction_count
        source_created_at = str(session["metadata"].get("timestamp", ""))
        date_part = re.sub(r"[^0-9]", "", source_created_at)[:8] or "undated"
        output_file = f"session_{position:02d}_{date_part}_{session_id[:8]}.md"
        output_path = output_dir / output_file
        output_path.write_text(
            render_session_markdown(
                title=title,
                session_id=session_id,
                messages=messages,
                source_created_at=source_created_at,
                redaction_count=redaction_count,
            ),
            encoding="utf-8",
        )

        exported.append(
            {
                "session_id": session_id,
                "title": title,
                "source_created_at": source_created_at,
                "source_file_name": session["source_path"].name,
                "source_sha256": file_sha256(session["source_path"]),
                "output_file": output_file,
                "user_message_count": sum(message["role"] == "user" for message in messages),
                "assistant_message_count": sum(message["role"] == "assistant" for message in messages),
                "redaction_count": redaction_count,
                "input_parse_errors": session["parse_errors"],
            }
        )

    (output_dir / "index.md").write_text(render_index_markdown(exported), encoding="utf-8")
    report = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "project_root": "$PROJECT_ROOT",
        "session_count": len(exported),
        "scope": "visible user and assistant event messages only",
        "redaction_policy": "credentials and local absolute paths are masked",
        "ordering": "session metadata timestamp ascending",
        "environment": {
            "python_version": sys.version.split()[0],
            "platform": platform.platform(),
            "script": Path(__file__).name,
        },
        "sessions": exported,
    }
    (output_dir / "export_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export visible Codex project sessions into ChatGPT-ready Markdown files."
    )
    parser.add_argument(
        "--codex-home",
        type=Path,
        default=Path.home() / ".codex",
        help="Codex local data directory (default: ~/.codex).",
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=PROJECT_ROOT,
        help="Project root used to select matching session metadata.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="New empty export directory. Defaults to a timestamped docs/chatgpt_discussions folder.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve()
    output_dir = args.output_dir
    if output_dir is None:
        stamp = datetime.now().strftime("%Y%m%dT%H%M%S")
        output_dir = project_root / "docs" / "chatgpt_discussions" / f"export_{stamp}"
    report = export_project_sessions(
        codex_home=args.codex_home,
        project_root=project_root,
        output_dir=output_dir,
    )
    print(
        json.dumps(
            {
                "output_dir": str(Path(output_dir)),
                "session_count": report["session_count"],
                "report": str(Path(output_dir) / "export_report.json"),
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
