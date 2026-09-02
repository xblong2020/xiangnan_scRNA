from __future__ import annotations

import argparse
import json
import os
import shutil
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
    parser = argparse.ArgumentParser(description="Run formal cisTarget ctx motif pruning for Module 6.3b.")
    parser.add_argument(
        "--adjacency",
        type=Path,
        default=ROOT / "metadata/driver/scenic_module6_3b/driver_module6_3b_grnboost2_seed777_adjacencies.tsv.gz",
    )
    parser.add_argument(
        "--expression-loom",
        type=Path,
        default=ROOT / "data/processed/driver/scenic_module6_3b_formal/driver_union_full_expression_counts.loom",
    )
    parser.add_argument(
        "--ranking-10kb",
        type=Path,
        default=Path("C:/SCENIC63b_work/db/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"),
    )
    parser.add_argument(
        "--ranking-proximal",
        type=Path,
        default=Path("C:/SCENIC63b_work/db/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"),
    )
    parser.add_argument(
        "--annotations",
        type=Path,
        default=Path("C:/SCENIC63b_work/db/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"),
    )
    parser.add_argument("--work-dir", type=Path, default=Path("C:/SCENIC63b_work/ctx"))
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/driver/scenic_module6_3b")
    parser.add_argument("--num-workers", type=int, default=8)
    parser.add_argument("--chunk-size", type=int, default=20)
    parser.add_argument("--seed", type=int, default=777)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> None:
    start = time.time()
    args = parse_args()
    for path in [args.adjacency, args.expression_loom, args.ranking_10kb, args.ranking_proximal, args.annotations]:
        if not path.exists():
            raise FileNotFoundError(path)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    raw_output = args.work_dir / "driver_module6_3b_cistarget_ctx_seed777.tsv"
    final_output = args.metadata_dir / "driver_module6_3b_canonical_regulons_seed777.tsv"
    log_path = args.work_dir / "driver_module6_3b_cistarget_ctx_seed777.log"
    report_path = args.metadata_dir / "driver_module6_3b_ctx_report.json"
    command = [
        sys.executable,
        str(ROOT / "scripts/run_pyscenic_cli_with_numpy_compat.py"),
        "ctx",
        str(args.adjacency),
        str(args.ranking_10kb),
        str(args.ranking_proximal),
        "--annotations_fname",
        str(args.annotations),
        "--expression_mtx_fname",
        str(args.expression_loom),
        "--output",
        str(raw_output),
        "--mode",
        "dask_multiprocessing",
        "--num_workers",
        str(args.num_workers),
        "--chunk_size",
        str(args.chunk_size),
        "--cell_id_attribute",
        "CellID",
        "--gene_attribute",
        "Gene",
        "--thresholds",
        "0.75",
        "0.90",
        "--top_n_targets",
        "50",
        "--top_n_regulators",
        "5",
        "10",
        "50",
        "--min_genes",
        "20",
        "--rank_threshold",
        "5000",
        "--auc_threshold",
        "0.05",
        "--nes_threshold",
        "3.0",
    ]
    env = os.environ.copy()
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["NUMEXPR_NUM_THREADS"] = "1"
    report = {
        "module": "6.3b",
        "method": "pySCENIC ctx cisTarget motif pruning",
        "seed": args.seed,
        "inputs": {
            "adjacency": record_path(args.adjacency),
            "expression_loom": record_path(args.expression_loom),
            "ranking_10kb": str(args.ranking_10kb),
            "ranking_proximal": str(args.ranking_proximal),
            "annotations": str(args.annotations),
        },
        "command": command,
        "status": "RUNNING",
        "outputs": {"raw_ctx": str(raw_output), "canonical_regulons": record_path(final_output), "log": str(log_path)},
    }
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    if final_output.exists() and not args.force:
        report["status"] = "ALREADY_PRESENT"
        report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
        return
    if raw_output.exists() and not args.force:
        raise FileExistsError(f"Raw ctx output exists; use --force after inspection: {raw_output}")
    try:
        with log_path.open("w", encoding="utf-8") as log:
            log.write("COMMAND: " + " ".join(command) + "\n")
            log.flush()
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, env=env)
        shutil.copy2(raw_output, final_output)
        report["status"] = "CTX_COMPLETE"
        report["raw_size_bytes"] = raw_output.stat().st_size
        report["canonical_regulons_size_bytes"] = final_output.stat().st_size
    except Exception as exc:
        report["status"] = "CTX_FAILED"
        report["error"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        report["elapsed_seconds"] = round(time.time() - start, 3)
        report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
