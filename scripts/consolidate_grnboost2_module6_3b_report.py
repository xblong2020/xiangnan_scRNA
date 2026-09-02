from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REPORT_DIR = ROOT / "metadata/driver/scenic_module6_3b"


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def main() -> None:
    direct_path = REPORT_DIR / "driver_module6_3b_grnboost2_direct_report.json"
    direct = load(direct_path)
    qc = load(REPORT_DIR / "driver_module6_3b_grn_qc_report.json")
    old = load(REPORT_DIR / "driver_module6_3b_grnboost2_report.json")
    if direct.get("status") != "GRN_COMPLETE":
        raise SystemExit("The direct GRNBoost2 report is not complete.")
    report = {
        "module": "6.3b",
        "method": "GRNBoost2",
        "status": "GRN_COMPLETE",
        "authoritative_run": "direct_write_before_shutdown",
        "seed": direct.get("seed"),
        "num_workers": direct.get("num_workers"),
        "input_loom": direct.get("input_loom"),
        "input_tf_list": direct.get("input_tf_list"),
        "n_cells": direct.get("n_cells"),
        "n_genes": direct.get("n_genes"),
        "n_tfs": direct.get("n_tfs"),
        "n_edges": direct.get("n_edges"),
        "elapsed_seconds": direct.get("elapsed_seconds"),
        "outputs": {
            "adjacencies": direct.get("outputs", {}).get("final"),
            "raw_adjacencies": direct.get("outputs", {}).get("raw"),
            "direct_report": str(direct_path),
            "qc": qc.get("outputs", {}).get("qc"),
            "top_edges": qc.get("outputs", {}).get("top_edges"),
            "qc_report": qc.get("outputs", {}).get("report"),
        },
        "attempt_history": [
            {
                "attempt": "pyscenic_cli_grn",
                "status": old.get("status", "GRN_FAILED"),
                "error": old.get("error", "CLI shutdown failure before output write"),
                "output_written": False,
            },
            {
                "attempt": "direct_write_before_shutdown",
                "status": direct.get("status"),
                "cluster_close_error": direct.get("cluster_close_error"),
                "output_written": True,
            },
        ],
    }
    (REPORT_DIR / "driver_module6_3b_grnboost2_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
