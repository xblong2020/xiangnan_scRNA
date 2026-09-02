from __future__ import annotations

import argparse
import subprocess
import time
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RSCRIPT = Path(r"C:\Program Files\R\R-4.6.0\bin\Rscript.exe")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run CopyKAT one manifest run at a time and stop after prediction is written.")
    parser.add_argument("--manifest", type=Path, default=ROOT / "metadata/copykat_module3/copykat_input_manifest.tsv")
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/copykat_module3")
    parser.add_argument("--runs", default="", help="Optional comma-separated run_id allowlist.")
    parser.add_argument("--poll-seconds", type=int, default=20)
    parser.add_argument("--stable-checks", type=int, default=2)
    parser.add_argument("--max-seconds-per-run", type=int, default=7200)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def resolve_path(value: object) -> Path:
    path = Path(str(value))
    if path.exists():
        return path
    text = str(value).replace("\\", "/")
    for anchor in ("data/processed/", "metadata/"):
        idx = text.find(anchor)
        if idx >= 0:
            candidate = ROOT / Path(text[idx:])
            if candidate.exists():
                return candidate
    return path


def prediction_path(run_id: str, run_dir: Path) -> Path:
    return run_dir / f"{run_id}_copykat_prediction.txt"


def write_status(path: Path, rows: list[dict[str, object]]) -> None:
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def kill_process_tree(proc: subprocess.Popen) -> None:
    if proc.poll() is not None:
        return
    subprocess.run(["taskkill", "/PID", str(proc.pid), "/T", "/F"], capture_output=True, text=True)


def file_is_stable(path: Path, poll_seconds: int, stable_checks: int) -> bool:
    if not path.exists() or path.stat().st_size <= 0:
        return False
    last_size = path.stat().st_size
    for _ in range(stable_checks):
        time.sleep(poll_seconds)
        if not path.exists() or path.stat().st_size <= 0:
            return False
        size = path.stat().st_size
        if size != last_size:
            last_size = size
            continue
    return True


def run_one(row: pd.Series, args: argparse.Namespace) -> dict[str, object]:
    run_id = str(row["run_id"])
    run_dir = resolve_path(row["run_dir"])
    pred_path = prediction_path(run_id, run_dir)
    result_path = run_dir / f"{run_id}_copykat_result.rds"
    if pred_path.exists() and pred_path.stat().st_size > 0 and not args.overwrite:
        return {
            "run_id": run_id,
            "status": "exists" if result_path.exists() else "exists_prediction_only",
            "run_dir": str(run_dir.resolve()),
            "prediction_path": str(pred_path.resolve()),
            "result_rds": str(result_path.resolve()) if result_path.exists() else "",
            "n_cells": int(row.get("n_cells", 0)),
            "n_candidate": int(row.get("n_candidate", 0)),
            "n_reference": int(row.get("n_reference", 0)),
            "elapsed_seconds": 0,
            "message": "",
        }

    stdout_path = run_dir / f"{run_id}_watchdog.stdout.log"
    stderr_path = run_dir / f"{run_id}_watchdog.stderr.log"
    cmd = [
        str(RSCRIPT),
        str(ROOT / "scripts/run_copykat_module3.R"),
        "--runs",
        run_id,
        "--overwrite",
    ]
    start = time.time()
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
        proc = subprocess.Popen(cmd, cwd=ROOT, stdout=stdout, stderr=stderr)
        stable = False
        while proc.poll() is None:
            elapsed = time.time() - start
            if elapsed > args.max_seconds_per_run:
                kill_process_tree(proc)
                return {
                    "run_id": run_id,
                    "status": "timeout",
                    "run_dir": str(run_dir.resolve()),
                    "prediction_path": str(pred_path.resolve()) if pred_path.exists() else "",
                    "result_rds": str(result_path.resolve()) if result_path.exists() else "",
                    "n_cells": int(row.get("n_cells", 0)),
                    "n_candidate": int(row.get("n_candidate", 0)),
                    "n_reference": int(row.get("n_reference", 0)),
                    "elapsed_seconds": round(elapsed, 1),
                    "message": f"exceeded {args.max_seconds_per_run}s",
                }
            if pred_path.exists() and pred_path.stat().st_size > 0:
                stable = file_is_stable(pred_path, args.poll_seconds, args.stable_checks)
                if stable:
                    kill_process_tree(proc)
                    break
            time.sleep(args.poll_seconds)

        return_code = proc.poll()

    elapsed = round(time.time() - start, 1)
    if pred_path.exists() and pred_path.stat().st_size > 0:
        status = "ok_prediction_only" if return_code not in (0, None) or not result_path.exists() else "ok"
        message = f"return_code={return_code}; stopped_after_prediction={stable}"
    else:
        status = "error"
        message = f"return_code={return_code}; no prediction"
    return {
        "run_id": run_id,
        "status": status,
        "run_dir": str(run_dir.resolve()),
        "prediction_path": str(pred_path.resolve()) if pred_path.exists() else "",
        "result_rds": str(result_path.resolve()) if result_path.exists() else "",
        "n_cells": int(row.get("n_cells", 0)),
        "n_candidate": int(row.get("n_candidate", 0)),
        "n_reference": int(row.get("n_reference", 0)),
        "elapsed_seconds": elapsed,
        "message": message,
    }


def main() -> int:
    args = parse_args()
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(args.manifest, sep="\t", encoding="utf-8")
    allow = {x.strip() for x in args.runs.split(",") if x.strip()}
    if allow:
        manifest = manifest.loc[manifest["run_id"].astype(str).isin(allow)].copy()
    if manifest.empty:
        raise SystemExit("No runs selected")

    status_path = args.metadata_dir / "copykat_run_status.tsv"
    rows: list[dict[str, object]] = []
    for _, row in manifest.iterrows():
        run_id = str(row["run_id"])
        print(f"WATCHDOG {run_id} cells={row.get('n_cells')} candidates={row.get('n_candidate')}", flush=True)
        status = run_one(row, args)
        print(f"DONE {run_id} status={status['status']} elapsed={status['elapsed_seconds']}", flush=True)
        rows.append(status)
        write_status(status_path, rows)
    has_error = any(row["status"] in {"error", "timeout"} for row in rows)
    return 1 if has_error else 0


if __name__ == "__main__":
    raise SystemExit(main())
