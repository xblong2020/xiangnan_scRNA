from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import time
from pathlib import Path
from urllib.request import Request, urlopen

import anndata as ad
import pandas as pd
import pyarrow.feather as feather
from ctxcore.rnkdb import FeatherRankingDatabase


ROOT = Path(__file__).resolve().parents[1]

RESOURCES = {
    "ranking_10kb": {
        "filename": "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        "url": "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        "sha1_url": "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather.sha1sum.txt",
        "kind": "ranking",
    },
    "ranking_proximal": {
        "filename": "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        "url": "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        "sha1_url": "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather.sha1sum.txt",
        "kind": "ranking",
    },
    "motif_annotation": {
        "filename": "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        "url": "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        "kind": "annotation",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download and validate matching hg38 mc_v10_clust SCENIC resources.")
    parser.add_argument("--resource-dir", type=Path, default=ROOT / "metadata/driver/scenic_resources_v10")
    parser.add_argument(
        "--gene-input",
        type=Path,
        default=ROOT / "data/processed/driver/driver_union_full_expression.module6_3b.formal.h5ad",
    )
    parser.add_argument("--reports-dir", type=Path, default=ROOT / "reports")
    parser.add_argument("--download", action="store_true", help="Download missing or incomplete resources with resume support.")
    parser.add_argument("--skip-proximal", action="store_true")
    return parser.parse_args()


def sha1_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha1()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def head_size(url: str) -> int | None:
    try:
        with urlopen(Request(url, method="HEAD"), timeout=60) as response:
            value = response.headers.get("content-length")
            return int(value) if value else None
    except Exception:
        return None


def expected_sha1(url: str) -> str | None:
    try:
        with urlopen(url, timeout=60) as response:
            text = response.read().decode("utf-8", "replace").strip()
        return text.split()[0].lower() if text else None
    except Exception:
        return None


def download_resume(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    command = [
        "curl.exe",
        "-L",
        "-C",
        "-",
        "--fail",
        "--retry",
        "5",
        "--retry-delay",
        "5",
        "--connect-timeout",
        "30",
        "-o",
        str(destination),
        url,
    ]
    subprocess.run(command, check=True)


def validate_ranking(path: Path, input_genes: set[str]) -> dict[str, object]:
    result: dict[str, object] = {
        "path": str(path),
        "exists": path.exists(),
        "size_bytes": path.stat().st_size if path.exists() else 0,
    }
    if not path.exists():
        result["status"] = "MISSING"
        return result
    with path.open("rb") as handle:
        result["head_magic"] = handle.read(6).decode("ascii", "replace")
        handle.seek(-6, 2)
        result["tail_magic"] = handle.read(6).decode("ascii", "replace")
    result["arrow_magic_ok"] = result["head_magic"] == "ARROW1" and result["tail_magic"] == "ARROW1"
    try:
        table = feather.read_table(path, memory_map=True)
        result["pyarrow_readable"] = True
        result["schema_rows"] = int(table.num_rows)
        result["schema_columns"] = int(table.num_columns)
        result["column_examples"] = table.column_names[:5]
    except Exception as exc:
        result["pyarrow_readable"] = False
        result["pyarrow_error"] = f"{type(exc).__name__}: {exc}"
    try:
        db = FeatherRankingDatabase(str(path), path.stem)
        db_genes = set(map(str, db.genes))
        result["ctxcore_readable"] = True
        result["ctxcore_total_genes"] = int(db.total_genes)
        result["input_gene_coverage"] = int(len(input_genes & db_genes))
        result["input_gene_missing"] = int(len(input_genes - db_genes))
    except Exception as exc:
        result["ctxcore_readable"] = False
        result["ctxcore_error"] = f"{type(exc).__name__}: {exc}"
    readable = result.get("arrow_magic_ok") and result.get("pyarrow_readable") and result.get("ctxcore_readable")
    result["gene_coverage_warning"] = bool(result.get("input_gene_missing", 0) > 0)
    result["status"] = "PASS_WITH_GENE_GAPS" if readable and result["gene_coverage_warning"] else ("PASS" if readable else "FAIL")
    return result


def validate_annotation(path: Path) -> dict[str, object]:
    result: dict[str, object] = {
        "path": str(path),
        "exists": path.exists(),
        "size_bytes": path.stat().st_size if path.exists() else 0,
    }
    if not path.exists():
        result["status"] = "MISSING"
        return result
    try:
        frame = pd.read_csv(path, sep="\t", nrows=5)
        result["readable"] = True
        result["columns"] = list(frame.columns)
        result["required_columns_present"] = {"#motif_id", "gene_name"}.issubset(set(frame.columns))
        result["status"] = "PASS" if result["required_columns_present"] else "FAIL"
    except Exception as exc:
        result["readable"] = False
        result["error"] = f"{type(exc).__name__}: {exc}"
        result["status"] = "FAIL"
    return result


def main() -> None:
    start = time.time()
    args = parse_args()
    args.resource_dir.mkdir(parents=True, exist_ok=True)
    args.reports_dir.mkdir(parents=True, exist_ok=True)
    input_genes = set(map(str, ad.read_h5ad(args.gene_input, backed="r").var_names))
    selected = {key: value for key, value in RESOURCES.items() if not (args.skip_proximal and key == "ranking_proximal")}
    records: dict[str, object] = {
        "module": "6.3b",
        "resource_set": "hg38 refseq_r80 mc_v10_clust + motifs-v10nr_clust-nr HGNC",
        "download_enabled": args.download,
        "input_gene_count": len(input_genes),
        "resources": {},
    }
    for key, spec in selected.items():
        destination = args.resource_dir / spec["filename"]
        expected_size = head_size(spec["url"])
        expected_hash = expected_sha1(spec["sha1_url"]) if spec.get("sha1_url") else None
        if args.download and (not destination.exists() or (expected_size and destination.stat().st_size != expected_size)):
            download_resume(spec["url"], destination)
        item = (
            validate_ranking(destination, input_genes)
            if spec["kind"] == "ranking"
            else validate_annotation(destination)
        )
        item["url"] = spec["url"]
        item["expected_size_bytes"] = expected_size
        item["expected_sha1"] = expected_hash
        if destination.exists():
            item["sha1"] = sha1_file(destination)
            item["sha1_match"] = expected_hash is None or item["sha1"] == expected_hash
        else:
            item["sha1_match"] = False
        if item.get("status") == "PASS" and not item.get("sha1_match", True):
            item["status"] = "FAIL"
        records["resources"][key] = item
    statuses = [item.get("status") for item in records["resources"].values()]
    records["status"] = "RESOURCES_COMPLETE" if all(status in {"PASS", "PASS_WITH_GENE_GAPS"} for status in statuses) else "RESOURCES_INCOMPLETE"
    records["outputs"] = {"resource_dir": str(args.resource_dir)}
    records["elapsed_seconds"] = round(time.time() - start, 3)
    validation_path = args.resource_dir / "resource_validation.json"
    validation_path.write_text(json.dumps(records, indent=2), encoding="utf-8")

    lines = [
        "# Module 6.3b cisTarget resource validation",
        "",
        f"- Overall status: **{records['status']}**",
        "- Resource set: `hg38 refseq_r80 mc_v10_clust` gene-based rankings with `motifs-v10nr_clust-nr.hgnc` annotation.",
        "- Gene coverage gaps are retained as an explicit limitation when formal genes are outside the RefSeq 80 ranking universe.",
        f"- Input genes checked: `{len(input_genes)}`",
        "",
        "| Resource | Size | Arrow/reader | Gene coverage | SHA1 | Status |",
        "|---|---:|---|---:|---|---|",
    ]
    for key, item in records["resources"].items():
        reader = item.get("pyarrow_readable", item.get("readable", False))
        coverage = item.get("input_gene_coverage", "n/a")
        sha_status = item.get("sha1_match", False)
        lines.append(f"| {key} | {item.get('size_bytes', 0)} | {reader} | {coverage} | {sha_status} | {item.get('status')} |")
    lines.extend(
        [
            "",
            "## Files",
            "",
            f"- Validation JSON: `{validation_path}`",
            f"- Resource directory: `{args.resource_dir}`",
            "",
            "## Matching rule",
            "",
            "The two ranking databases and the motif annotation use the v10 clust resource family. The old mc9nr ranking database and v9 annotation remain historical references and are not used for formal 6.3b ctx.",
        ]
    )
    (args.reports_dir / "module6_3b_cistarget_resource_validation.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    if records["status"] != "RESOURCES_COMPLETE":
        raise SystemExit("cisTarget resource validation failed; inspect resource_validation.json before ctx.")


if __name__ == "__main__":
    main()
