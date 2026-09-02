from __future__ import annotations

import argparse
import gzip
import json
import os
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run formal full-expression GRNBoost2 for Module 6.3b.")
    parser.add_argument(
        "--loom",
        type=Path,
        default=ROOT / "data/processed/driver/scenic_module6_3b_formal/driver_union_full_expression_counts.loom",
    )
    parser.add_argument(
        "--tf-list",
        type=Path,
        default=ROOT / "data/processed/driver/scenic_module6_3b_formal/driver_union_tfs_in_matrix.txt",
    )
    parser.add_argument("--work-dir", type=Path, default=Path("C:/SCENIC63b_work/grn"))
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/driver/scenic_module6_3b")
    parser.add_argument("--num-workers", type=int, default=32)
    parser.add_argument("--seed", type=int, default=777)
    parser.add_argument("--timeout-seconds", type=int, default=0, help="0 means no subprocess timeout.")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def compress_tsv(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with source.open("rb") as src, gzip.open(destination, "wb", compresslevel=6) as dst:
        while True:
            chunk = src.read(1024 * 1024)
            if not chunk:
                break
            dst.write(chunk)


def main() -> None:
    start = time.time()
    args = parse_args()
    if not args.loom.exists():
        raise FileNotFoundError(args.loom)
    if not args.tf_list.exists():
        raise FileNotFoundError(args.tf_list)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    raw_output = args.work_dir / "driver_module6_3b_grnboost2_seed777_adjacencies.tsv"
    final_output = args.metadata_dir / "driver_module6_3b_grnboost2_seed777_adjacencies.tsv.gz"
    log_path = args.work_dir / "driver_module6_3b_grnboost2_seed777.log"
    command = [
        sys.executable,
        str(ROOT / "scripts/run_pyscenic_cli_with_numpy_compat.py"),
        "grn",
        str(args.loom),
        str(args.tf_list),
        "--method",
        "grnboost2",
        "--output",
        str(raw_output),
        "--num_workers",
        str(args.num_workers),
        "--client_or_address",
        "local",
        "--seed",
        str(args.seed),
        "--cell_id_attribute",
        "CellID",
        "--gene_attribute",
        "Gene",
        "--sparse",
    ]
    env = os.environ.copy()
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["NUMEXPR_NUM_THREADS"] = "1"
    report_path = args.metadata_dir / "driver_module6_3b_grnboost2_report.json"
    report = {
        "module": "6.3b",
        "method": "GRNBoost2",
        "seed": args.seed,
        "input_loom": str(args.loom),
        "input_tf_list": str(args.tf_list),
        "work_dir": str(args.work_dir),
        "command": command,
        "num_workers": args.num_workers,
        "status": "RUNNING",
        "started_at_epoch": time.time(),
        "outputs": {
            "raw_adjacencies": str(raw_output),
            "adjacencies": str(final_output),
            "log": str(log_path),
        },
    }
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    if final_output.exists() and not args.force:
        report["status"] = "ALREADY_PRESENT"
        report["elapsed_seconds"] = round(time.time() - start, 3)
        report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
        return
    if raw_output.exists() and not args.force:
        raise FileExistsError(f"Raw GRN output exists; use --force after inspecting it: {raw_output}")

    try:
        with log_path.open("w", encoding="utf-8") as log:
            log.write("COMMAND: " + " ".join(command) + "\n")
            log.flush()
            subprocess.run(
                command,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
                timeout=args.timeout_seconds or None,
                env=env,
            )
        compress_tsv(raw_output, final_output)
        report["status"] = "GRN_COMPLETE"
        report["raw_size_bytes"] = raw_output.stat().st_size
        report["compressed_size_bytes"] = final_output.stat().st_size
    except subprocess.TimeoutExpired as exc:
        report["status"] = "GRN_TIMEOUT"
        report["error"] = str(exc)
        raise
    except Exception as exc:
        report["status"] = "GRN_FAILED"
        report["error"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        report["elapsed_seconds"] = round(time.time() - start, 3)
        report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
