from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
import time
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run GRNBoost2 directly and write output before Dask shutdown.")
    parser.add_argument("--loom", type=Path, required=True)
    parser.add_argument("--tf-list", type=Path, required=True)
    parser.add_argument("--raw-output", type=Path, required=True)
    parser.add_argument("--final-output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--num-workers", type=int, default=8)
    parser.add_argument("--seed", type=int, default=777)
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
    for path in [args.loom, args.tf_list]:
        if not path.exists():
            raise FileNotFoundError(path)
    args.raw_output.parent.mkdir(parents=True, exist_ok=True)
    args.final_output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    report = {
        "module": "6.3b",
        "method": "direct Arboreto GRNBoost2 with write-before-shutdown",
        "seed": args.seed,
        "input_loom": str(args.loom),
        "input_tf_list": str(args.tf_list),
        "num_workers": args.num_workers,
        "status": "RUNNING",
        "outputs": {"raw": str(args.raw_output), "final": str(args.final_output)},
    }
    args.report.write_text(json.dumps(report, indent=2), encoding="utf-8")
    for name, value in {"object": object, "bool": bool, "int": int, "float": float}.items():
        if name not in np.__dict__:
            setattr(np, name, value)
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"

    cluster = None
    client = None
    try:
        from scripts.run_pyscenic_cli_with_numpy_compat import patch_arboreto_create_graph

        patch_arboreto_create_graph()
        from pyscenic.cli.utils import load_exp_matrix
        from arboreto.algo import grnboost2
        from dask.distributed import Client, LocalCluster

        tfs = [line.strip() for line in args.tf_list.read_text(encoding="utf-8").splitlines() if line.strip()]
        expression, genes, cells = load_exp_matrix(str(args.loom), False, True, "CellID", "Gene")
        report["n_cells"] = len(cells)
        report["n_genes"] = len(genes)
        report["n_tfs"] = len(tfs)
        cluster = LocalCluster(n_workers=args.num_workers, threads_per_worker=1, diagnostics_port=None)
        client = Client(cluster)
        network = grnboost2(
            expression_data=expression,
            gene_names=genes,
            tf_names=tfs,
            client_or_address=client,
            seed=args.seed,
            verbose=True,
        )
        network.to_csv(args.raw_output, sep="\t", index=False)
        compress_tsv(args.raw_output, args.final_output)
        report["n_edges"] = int(network.shape[0])
        report["raw_size_bytes"] = args.raw_output.stat().st_size
        report["final_size_bytes"] = args.final_output.stat().st_size
        report["status"] = "GRN_COMPLETE"
    except Exception as exc:
        report["status"] = "GRN_FAILED"
        report["error"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        report["elapsed_seconds"] = round(time.time() - start, 3)
        args.report.write_text(json.dumps(report, indent=2), encoding="utf-8")
        if client is not None:
            try:
                client.close(timeout=2)
            except Exception as exc:
                report["client_close_error"] = f"{type(exc).__name__}: {exc}"
        if cluster is not None:
            try:
                cluster.close(timeout=2)
            except Exception as exc:
                report["cluster_close_error"] = f"{type(exc).__name__}: {exc}"
        args.report.write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
