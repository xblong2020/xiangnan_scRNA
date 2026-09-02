from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def record_path(path: Path | str) -> str:
    value = Path(path)
    try:
        return str(value.resolve().relative_to(ROOT.resolve())).replace("\\", "/")
    except ValueError:
        return str(value)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="QC formal Module 6.3b GRNBoost2 adjacency output.")
    parser.add_argument(
        "--adjacency",
        type=Path,
        default=ROOT / "metadata/driver/scenic_module6_3b/driver_module6_3b_grnboost2_seed777_adjacencies.tsv.gz",
    )
    parser.add_argument("--tf-list", type=Path, default=ROOT / "data/processed/driver/scenic_module6_3b_formal/driver_union_tfs_in_matrix.txt")
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/driver/scenic_module6_3b")
    parser.add_argument("--reports-dir", type=Path, default=ROOT / "reports")
    return parser.parse_args()


def main() -> None:
    start = time.time()
    args = parse_args()
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    args.reports_dir.mkdir(parents=True, exist_ok=True)
    if not args.adjacency.exists():
        raise FileNotFoundError(args.adjacency)
    frame = pd.read_csv(args.adjacency, sep="\t")
    required = {"TF", "target", "importance"}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise KeyError(f"Missing adjacency columns: {missing}")
    tf_list = {line.strip() for line in args.tf_list.read_text(encoding="utf-8").splitlines() if line.strip()}
    importance = pd.to_numeric(frame["importance"], errors="coerce")
    duplicate_edges = int(frame.duplicated(subset=["TF", "target"]).sum())
    represented_tfs = set(frame["TF"].astype(str))
    targets = set(frame["target"].astype(str))
    invalid_tfs = represented_tfs - tf_list
    finite = np.isfinite(importance.to_numpy()).all()
    positive = bool((importance > 0).all())
    qc = pd.DataFrame(
        [
            {"metric": "n_edges", "value": int(len(frame)), "status": "PASS" if len(frame) > 0 else "FAIL"},
            {"metric": "n_tfs_represented", "value": int(len(represented_tfs)), "status": "PASS" if len(represented_tfs) > 0 else "FAIL"},
            {"metric": "n_targets_represented", "value": int(len(targets)), "status": "PASS" if len(targets) > 0 else "FAIL"},
            {"metric": "duplicate_tf_target_edges", "value": duplicate_edges, "status": "PASS" if duplicate_edges == 0 else "FAIL"},
            {"metric": "nonfinite_importance", "value": int((~np.isfinite(importance.to_numpy())).sum()), "status": "PASS" if finite else "FAIL"},
            {"metric": "nonpositive_importance", "value": int((importance <= 0).sum()), "status": "PASS" if positive else "FAIL"},
            {"metric": "tfs_not_in_input_list", "value": len(invalid_tfs), "status": "PASS" if not invalid_tfs else "FAIL"},
            {"metric": "importance_min", "value": float(importance.min()), "status": "PASS"},
            {"metric": "importance_median", "value": float(importance.median()), "status": "PASS"},
            {"metric": "importance_max", "value": float(importance.max()), "status": "PASS"},
        ]
    )
    qc_path = args.metadata_dir / "driver_module6_3b_grn_qc.tsv"
    qc.to_csv(qc_path, sep="\t", index=False)
    top_path = args.metadata_dir / "driver_module6_3b_grn_top_edges.tsv"
    frame.sort_values("importance", ascending=False).head(1000).to_csv(top_path, sep="\t", index=False)
    passed = bool(qc["status"].eq("PASS").all())
    lines = [
        "# Module 6.3b GRNBoost2 QC",
        "",
        f"- Status: **{'PASS' if passed else 'FAIL'}**",
        f"- Adjacency: `{record_path(args.adjacency)}`",
        f"- Edges: `{len(frame)}`; TFs: `{len(represented_tfs)}`; targets: `{len(targets)}`",
        "",
        "| Metric | Value | Status |",
        "|---|---:|---|",
    ]
    for row in qc.itertuples(index=False):
        lines.append(f"| {row.metric} | {row.value} | {row.status} |")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This adjacency is accepted as formal GRNBoost2 only when all QC metrics pass. It must remain distinct from the historical Module 6.3 co-expression adjacency.",
            "",
            f"- QC table: `{record_path(qc_path)}`",
            f"- Top edges: `{record_path(top_path)}`",
        ]
    )
    report_path = args.reports_dir / "module6_3b_grn_qc_report.md"
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    report = {
        "module": "6.3b",
        "method": "GRNBoost2 adjacency QC",
        "adjacency": record_path(args.adjacency),
        "status": "GRN_QC_COMPLETE" if passed else "GRN_QC_FAILED",
        "n_edges": int(len(frame)),
        "n_tfs": int(len(represented_tfs)),
        "n_targets": int(len(targets)),
        "invalid_tfs": sorted(invalid_tfs),
        "outputs": {"qc": record_path(qc_path), "top_edges": record_path(top_path), "report": record_path(report_path)},
        "elapsed_seconds": round(time.time() - start, 3),
    }
    (args.metadata_dir / "driver_module6_3b_grn_qc_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    if not passed:
        raise SystemExit("GRNBoost2 QC failed; inspect the report before ctx.")


if __name__ == "__main__":
    main()
