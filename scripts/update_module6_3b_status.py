from __future__ import annotations

import argparse
import json
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_STATUS = ROOT / "metadata/driver/driver_module6_3b_canonical_scenic_status.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Update the Module 6.3b checkpoint status JSON.")
    parser.add_argument("--status", required=True)
    parser.add_argument("--stage", default=None)
    parser.add_argument("--error", default=None)
    parser.add_argument("--recommended-resume-point", default=None)
    parser.add_argument("--completed", action="append", default=[], help="KEY=PATH; may be supplied multiple times")
    parser.add_argument("--status-json", type=Path, default=DEFAULT_STATUS)
    return parser.parse_args()


def relative_path(path: str) -> str:
    value = Path(path)
    try:
        return str(value.resolve().relative_to(ROOT.resolve())).replace("\\", "/")
    except ValueError:
        return str(value)


def main() -> None:
    args = parse_args()
    args.status_json.parent.mkdir(parents=True, exist_ok=True)
    if args.status_json.exists():
        current = json.loads(args.status_json.read_text(encoding="utf-8"))
    else:
        current = {"module": "6.3b", "completed": {}, "history": []}
    current["module"] = "6.3b"
    current["status"] = args.status
    current.setdefault("completed", {})
    current.setdefault("history", [])
    for item in args.completed:
        if "=" not in item:
            raise ValueError(f"completed item must be KEY=PATH: {item}")
        key, path = item.split("=", 1)
        current["completed"][key] = relative_path(path)
    event = {"status": args.status, "timestamp_epoch": time.time()}
    if args.stage:
        current["running_stage"] = args.stage
        event["stage"] = args.stage
    else:
        current.pop("running_stage", None)
    if args.error:
        current["error"] = args.error
        event["error"] = args.error
    else:
        current.pop("error", None)
    if args.recommended_resume_point:
        current["recommended_resume_point"] = args.recommended_resume_point
        event["recommended_resume_point"] = args.recommended_resume_point
    current["history"].append(event)
    current["history"] = current["history"][-30:]
    args.status_json.write_text(json.dumps(current, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
