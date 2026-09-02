from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def record_path(path: Path | str) -> str:
    value = Path(path)
    try:
        return str(value.resolve().relative_to(ROOT.resolve())).replace("\\", "/")
    except ValueError:
        return str(value)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run AUCell on all formal Module 6.3b driver-union cells.")
    parser.add_argument(
        "--loom",
        type=Path,
        default=ROOT / "data/processed/driver/scenic_module6_3b_formal/driver_union_full_expression_counts.loom",
    )
    parser.add_argument(
        "--regulons",
        type=Path,
        default=ROOT / "metadata/driver/scenic_module6_3b/driver_module6_3b_canonical_regulons_seed777.tsv",
    )
    parser.add_argument("--work-dir", type=Path, default=Path("C:/SCENIC63b_work/aucell"))
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/driver/scenic_module6_3b")
    parser.add_argument("--num-workers", type=int, default=8)
    parser.add_argument("--seed", type=int, default=777)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> None:
    start = time.time()
    args = parse_args()
    for path in [args.loom, args.regulons]:
        if not path.exists():
            raise FileNotFoundError(path)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    output = args.metadata_dir / "driver_module6_3b_canonical_regulon_auc_9512.csv"
    log_path = args.work_dir / "driver_module6_3b_aucell_seed777.log"
    report_path = args.metadata_dir / "driver_module6_3b_aucell_report.json"
    command = [
        sys.executable,
        str(ROOT / "scripts/run_pyscenic_cli_with_numpy_compat.py"),
        "aucell",
        str(args.loom),
        str(args.regulons),
        "-o",
        str(output),
        "--num_workers",
        str(args.num_workers),
        "--seed",
        str(args.seed),
        "--cell_id_attribute",
        "CellID",
        "--gene_attribute",
        "Gene",
        "--rank_threshold",
        "5000",
        "--auc_threshold",
        "0.05",
        "--nes_threshold",
        "3.0",
    ]
    report = {
        "module": "6.3b",
        "method": "pySCENIC AUCell on all formal driver-union cells",
        "seed": args.seed,
        "inputs": {"loom": record_path(args.loom), "regulons": record_path(args.regulons)},
        "command": command,
        "status": "RUNNING",
        "outputs": {"auc_csv": record_path(output), "log": record_path(log_path)},
    }
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    if output.exists() and not args.force:
        report["status"] = "ALREADY_PRESENT"
        report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
        return
    env = os.environ.copy()
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["NUMEXPR_NUM_THREADS"] = "1"
    try:
        with log_path.open("w", encoding="utf-8") as log:
            log.write("COMMAND: " + " ".join(command) + "\n")
            log.flush()
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, env=env)
        report["status"] = "AUCELL_COMPLETE"
        report["auc_size_bytes"] = output.stat().st_size
    except Exception as exc:
        report["status"] = "AUCELL_FAILED"
        report["error"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        report["elapsed_seconds"] = round(time.time() - start, 3)
        report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
