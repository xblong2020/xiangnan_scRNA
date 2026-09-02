from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def record_path(path: Path | str) -> str:
    value = Path(path)
    try:
        return str(value.resolve().relative_to(ROOT.resolve())).replace("\\", "/")
    except ValueError:
        return str(value)
