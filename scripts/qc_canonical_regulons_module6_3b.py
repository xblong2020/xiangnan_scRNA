from __future__ import annotations

import argparse
import ast
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
    parser = argparse.ArgumentParser(description="Parse and QC formal Module 6.3b cisTarget canonical regulons.")
    parser.add_argument(
        "--ctx-output",
        type=Path,
        default=ROOT / "metadata/driver/scenic_module6_3b/driver_module6_3b_canonical_regulons_seed777.tsv",
    )
    parser.add_argument(
        "--tf-list",
        type=Path,
        default=ROOT / "data/processed/driver/scenic_module6_3b_formal/driver_union_tfs_in_matrix.txt",
    )
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/driver/scenic_module6_3b")
    parser.add_argument("--reports-dir", type=Path, default=ROOT / "reports")
    parser.add_argument("--min-nes", type=float, default=3.0)
    return parser.parse_args()


def flatten_columns(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    if isinstance(out.columns, pd.MultiIndex):
        out.columns = [str(bottom) if not str(bottom).startswith("Unnamed") else str(top) for top, bottom in out.columns]
    return out


def read_ctx(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", header=[0, 1], index_col=[0, 1])
    frame = flatten_columns(frame).reset_index()
    frame.columns = [str(column) for column in frame.columns]
    return frame


def parse_targets(value: object) -> list[tuple[str, float]]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    try:
        parsed = ast.literal_eval(str(value))
    except (ValueError, SyntaxError):
        return [(part.strip(), np.nan) for part in str(value).replace(",", ";").split(";") if part.strip()]
    result = []
    for item in parsed:
        if isinstance(item, (list, tuple)) and item:
            result.append((str(item[0]), float(item[1]) if len(item) > 1 else np.nan))
        elif isinstance(item, str):
            result.append((item, np.nan))
    return result


def main() -> None:
    start = time.time()
    args = parse_args()
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    args.reports_dir.mkdir(parents=True, exist_ok=True)
    if not args.ctx_output.exists():
        raise FileNotFoundError(args.ctx_output)
    frame = read_ctx(args.ctx_output)
    required = {"TF", "MotifID", "NES", "Annotation", "Context", "TargetGenes"}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise KeyError(f"ctx output missing columns: {missing}")
    frame["TF"] = frame["TF"].astype(str)
    frame["NES"] = pd.to_numeric(frame["NES"], errors="coerce")
    tf_set = {line.strip() for line in args.tf_list.read_text(encoding="utf-8").splitlines() if line.strip()}
    frame = frame.loc[frame["NES"].ge(args.min_nes)].copy()
    target_rows = []
    for row in frame.itertuples(index=False):
        for target, weight in parse_targets(row.TargetGenes):
            target_rows.append(
                {
                    "TF": row.TF,
                    "target": target,
                    "importance": weight,
                    "MotifID": row.MotifID,
                    "NES": row.NES,
                    "Annotation": row.Annotation,
                    "Context": row.Context,
                }
            )
    targets = pd.DataFrame(target_rows)
    if targets.empty:
        raise ValueError("No target genes parsed from ctx output after NES filtering.")
    targets = targets.sort_values(["TF", "target", "NES"], ascending=[True, True, False])
    targets = targets.drop_duplicates(subset=["TF", "target"], keep="first")
    summary_rows = []
    for tf, sub in frame.groupby("TF", sort=True):
        target_sub = targets.loc[targets["TF"].eq(tf)]
        annotations = sub["Annotation"].astype(str)
        summary_rows.append(
            {
                "TF": tf,
                "regulon": f"{tf}(+)",
                "motif_count": int(sub.shape[0]),
                "regulon_size": int(target_sub["target"].nunique()),
                "motif_NES_max": float(sub["NES"].max()),
                "motif_NES_median": float(sub["NES"].median()),
                "direct_annotation_fraction": float(annotations.str.contains("direct", case=False, na=False).mean()),
                "target_genes": ";".join(sorted(target_sub["target"].astype(str).unique())),
                "tf_in_input_list": tf in tf_set,
            }
        )
    summary = pd.DataFrame(summary_rows).sort_values(["motif_NES_max", "regulon_size"], ascending=[False, False])
    summary["primary_analysis_eligible"] = summary["regulon_size"].ge(10)
    motif_path = args.metadata_dir / "driver_module6_3b_motif_enrichment_seed777.csv"
    summary_path = args.metadata_dir / "driver_module6_3b_canonical_regulons_seed777.csv"
    target_path = args.metadata_dir / "driver_module6_3b_regulon_targets.tsv.gz"
    frame.to_csv(motif_path, index=False)
    summary.to_csv(summary_path, index=False)
    targets.to_csv(target_path, sep="\t", index=False, compression="gzip")

    duplicate_regulons = int(summary["target_genes"].duplicated().sum())
    invalid_tfs = sorted(set(summary.loc[~summary["tf_in_input_list"], "TF"]))
    qc = pd.DataFrame(
        [
            {"metric": "ctx_motif_rows_after_nes", "value": int(frame.shape[0]), "status": "PASS" if frame.shape[0] > 0 else "FAIL"},
            {"metric": "canonical_regulon_count", "value": int(summary.shape[0]), "status": "PASS" if summary.shape[0] > 0 else "FAIL"},
            {"metric": "target_edge_count", "value": int(targets.shape[0]), "status": "PASS" if targets.shape[0] > 0 else "FAIL"},
            {"metric": "invalid_tfs", "value": len(invalid_tfs), "status": "PASS" if not invalid_tfs else "FAIL"},
            {"metric": "duplicate_regulon_target_sets", "value": duplicate_regulons, "status": "PASS"},
            {"metric": "regulons_with_less_than_10_targets", "value": int((summary["regulon_size"] < 10).sum()), "status": "WARN"},
            {"metric": "primary_analysis_eligible_regulons", "value": int(summary["primary_analysis_eligible"].sum()), "status": "PASS" if int(summary["primary_analysis_eligible"].sum()) > 0 else "FAIL"},
            {"metric": "minimum_regulon_size", "value": int(summary["regulon_size"].min()), "status": "PASS"},
            {"metric": "median_regulon_size", "value": float(summary["regulon_size"].median()), "status": "PASS"},
            {"metric": "maximum_regulon_size", "value": int(summary["regulon_size"].max()), "status": "PASS"},
        ]
    )
    qc_path = args.metadata_dir / "driver_module6_3b_canonical_regulon_qc.tsv"
    qc.to_csv(qc_path, sep="\t", index=False)
    passed = bool(~qc["status"].eq("FAIL").any())
    warning = bool(qc["status"].eq("WARN").any())
    lines = [
        "# Module 6.3b canonical regulon QC",
        "",
        f"- Status: **{'PASS_WITH_WARNINGS' if passed and warning else ('PASS' if passed else 'FAIL')}**",
            f"- Source ctx output: `{record_path(args.ctx_output)}`",
        f"- Motif rows after NES >= {args.min_nes}: `{frame.shape[0]}`",
        f"- Canonical regulons: `{summary.shape[0]}`",
        "",
        "| Metric | Value | Status |",
        "|---|---:|---|",
    ]
    for row in qc.itertuples(index=False):
        lines.append(f"| {row.metric} | {row.value} | {row.status} |")
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- Canonical regulon summary: `{record_path(summary_path)}`",
            f"- Motif enrichment: `{record_path(motif_path)}`",
            f"- Regulon targets: `{record_path(target_path)}`",
            f"- QC table: `{record_path(qc_path)}`",
            "",
            "The summary is generated from the formal GRNBoost2-derived ctx output. Historical Module 6.3/6.3c regulons are not merged into this table.",
        ]
    )
    report_path = args.reports_dir / "module6_3b_canonical_regulon_qc.md"
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    report = {
        "module": "6.3b",
        "method": "canonical regulon QC after cisTarget ctx",
        "status": "CTX_QC_COMPLETE_WITH_WARNINGS" if passed and warning else ("CTX_QC_COMPLETE" if passed else "CTX_QC_FAILED"),
        "n_motif_rows": int(frame.shape[0]),
        "n_canonical_regulons": int(summary.shape[0]),
        "n_target_edges": int(targets.shape[0]),
        "invalid_tfs": invalid_tfs,
        "outputs": {
            "canonical_regulons": record_path(summary_path),
            "motif_enrichment": record_path(motif_path),
            "regulon_targets": record_path(target_path),
            "qc": record_path(qc_path),
            "report": record_path(report_path),
        },
        "elapsed_seconds": round(time.time() - start, 3),
    }
    (args.metadata_dir / "driver_module6_3b_canonical_regulon_qc_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    if not passed:
        raise SystemExit("Canonical regulon QC failed; inspect the report before AUCell.")


if __name__ == "__main__":
    main()
