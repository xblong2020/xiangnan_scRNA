from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path
from typing import Any, Sequence

import pandas as pd
import requests


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA_DIR = ROOT / "metadata/driver"
DEFAULT_UP_GMT = DEFAULT_METADATA_DIR / "module9_4_drug_reversal_up.gmt"
DEFAULT_DOWN_GMT = DEFAULT_METADATA_DIR / "module9_4_drug_reversal_down.gmt"
DEFAULT_GENE_INFO = DEFAULT_METADATA_DIR / "GSE92742_Broad_LINCS_gene_info.txt.gz"
DEFAULT_JOBS_URL = "https://api.clue.io/api/jobs"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 9.8 CLUE sig_gutc_tool query.")
    parser.add_argument("--up-gmt", type=Path, default=DEFAULT_UP_GMT)
    parser.add_argument("--down-gmt", type=Path, default=DEFAULT_DOWN_GMT)
    parser.add_argument("--gene-info", type=Path, default=DEFAULT_GENE_INFO)
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA_DIR)
    parser.add_argument("--jobs-url", default=DEFAULT_JOBS_URL)
    parser.add_argument("--api-key-env", default="CLUE_API_KEY")
    parser.add_argument("--query-name", default="module9_8_drug_reversal")
    parser.add_argument("--existing-job-id", default="")
    parser.add_argument("--poll-seconds", type=float, default=10.0)
    parser.add_argument("--max-polls", type=int, default=90)
    parser.add_argument("--timeout-seconds", type=int, default=180)
    return parser.parse_args()


def read_gmt_genes(path: Path) -> list[str]:
    parts = path.read_text(encoding="utf-8").strip().split("\t")
    if len(parts) < 3:
        raise ValueError(f"invalid GMT: {path}")
    return list(dict.fromkeys(gene.strip().upper() for gene in parts[2:] if gene.strip()))


def map_symbols_to_entrez(
    genes: Sequence[str],
    gene_info: pd.DataFrame,
) -> tuple[list[str], list[str], int]:
    required = {"pr_gene_symbol", "pr_gene_id", "pr_is_bing"}
    missing = required.difference(gene_info.columns)
    if missing:
        raise ValueError(f"gene info missing columns: {sorted(missing)}")
    lookup = (
        gene_info.dropna(subset=["pr_gene_symbol", "pr_gene_id"])
        .drop_duplicates("pr_gene_symbol")
        .set_index("pr_gene_symbol")
    )
    mapped = [
        str(int(lookup.loc[gene, "pr_gene_id"]))
        for gene in genes
        if gene in lookup.index
    ]
    missing_genes = [gene for gene in genes if gene not in lookup.index]
    n_bing = int(
        pd.to_numeric(
            lookup.loc[[gene for gene in genes if gene in lookup.index], "pr_is_bing"],
            errors="coerce",
        )
        .fillna(0)
        .sum()
    )
    return mapped, missing_genes, n_bing


def paired_gmt_text(query_name: str, gene_ids: Sequence[str]) -> str:
    return "\t".join([query_name, "module9_8_clue_entrez", *gene_ids]) + "\n"


def safe_job_status(job: dict[str, Any], results_archive: Path) -> dict[str, Any]:
    return {
        "job_id": job.get("job_id"),
        "status": job.get("status"),
        "download_status": job.get("download_status"),
        "created": job.get("created"),
        "last_modified": job.get("last_modified"),
        "tool_id": job.get("tool_id"),
        "tool_version": job.get("tool_version"),
        "error_message": job.get("errorMessage", ""),
        "api_key_present": True,
        "secret_recorded": False,
        "results_archive": str(results_archive.resolve()),
    }


def submit_job(
    jobs_url: str,
    api_key: str,
    query_name: str,
    up_gmt: str,
    down_gmt: str,
    timeout_seconds: int,
) -> str:
    payload = {
        "tool_id": "sig_gutc_tool",
        "data_type": "L1000",
        "name": query_name,
        "dataset": "Touchstone",
        "ignoreWarnings": True,
        "uptag-cmapfile": up_gmt.strip(),
        "dntag-cmapfile": down_gmt.strip(),
    }
    response = requests.post(
        jobs_url,
        headers={"user_key": api_key, "Accept": "application/json"},
        json=payload,
        timeout=timeout_seconds,
    )
    result = response.json()
    if response.status_code != 200:
        errors = result.get("both") or result.get("system") or {}
        raise RuntimeError(f"CLUE submission failed with HTTP {response.status_code}: {errors}")
    job_id = result.get("result", {}).get("job_id")
    if not job_id:
        raise ValueError("CLUE submission response missing job_id")
    return str(job_id)


def list_jobs(jobs_url: str, api_key: str, timeout_seconds: int) -> list[dict[str, Any]]:
    response = requests.get(
        jobs_url,
        headers={"user_key": api_key, "Accept": "application/json"},
        timeout=timeout_seconds,
    )
    response.raise_for_status()
    payload = response.json()
    if not isinstance(payload, list):
        raise ValueError("CLUE jobs endpoint returned a non-list response")
    return payload


def wait_for_job(
    jobs_url: str,
    api_key: str,
    job_id: str,
    poll_seconds: float,
    max_polls: int,
    timeout_seconds: int,
) -> dict[str, Any]:
    for _ in range(max_polls):
        jobs = list_jobs(jobs_url, api_key, timeout_seconds)
        job = next((item for item in jobs if str(item.get("job_id")) == job_id), None)
        if job is not None and str(job.get("status", "")).lower() in {
            "completed",
            "failed",
            "error",
        }:
            return job
        time.sleep(poll_seconds)
    raise TimeoutError(f"CLUE job {job_id} did not finish after {max_polls} polls")


def download_results(job: dict[str, Any], destination: Path, timeout_seconds: int) -> None:
    url = str(job.get("download_url", ""))
    if not url:
        raise ValueError("completed CLUE job has no download_url")
    if url.startswith("//"):
        url = "https:" + url
    response = requests.get(url, timeout=max(timeout_seconds, 600))
    response.raise_for_status()
    destination.write_bytes(response.content)


def main() -> None:
    args = parse_args()
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    api_key = os.environ.get(args.api_key_env, "")
    if not api_key:
        raise EnvironmentError(f"{args.api_key_env} is not configured")

    gene_info = pd.read_csv(args.gene_info, sep="\t", compression="infer")
    up_genes = read_gmt_genes(args.up_gmt)
    down_genes = read_gmt_genes(args.down_gmt)
    up_ids, up_missing, up_bing = map_symbols_to_entrez(up_genes, gene_info)
    down_ids, down_missing, down_bing = map_symbols_to_entrez(down_genes, gene_info)
    up_text = paired_gmt_text(args.query_name, up_ids)
    down_text = paired_gmt_text(args.query_name, down_ids)

    up_output = args.metadata_dir / "module9_8_clue_up_entrez.gmt"
    down_output = args.metadata_dir / "module9_8_clue_down_entrez.gmt"
    archive_output = args.metadata_dir / "module9_8_clue_results.tar.gz"
    status_output = args.metadata_dir / "module9_8_clue_job_status.json"
    qc_output = args.metadata_dir / "module9_8_clue_input_qc.json"
    up_output.write_text(up_text, encoding="utf-8")
    down_output.write_text(down_text, encoding="utf-8")

    job_id = args.existing_job_id or submit_job(
        args.jobs_url,
        api_key,
        args.query_name,
        up_text,
        down_text,
        args.timeout_seconds,
    )
    job = wait_for_job(
        args.jobs_url,
        api_key,
        job_id,
        args.poll_seconds,
        args.max_polls,
        args.timeout_seconds,
    )
    if str(job.get("status", "")).lower() != "completed":
        raise RuntimeError(f"CLUE job ended with status {job.get('status')}")
    download_results(job, archive_output, args.timeout_seconds)

    status_output.write_text(
        json.dumps(safe_job_status(job, archive_output), indent=2, sort_keys=True),
        encoding="utf-8",
    )
    qc = {
        "query_name": args.query_name,
        "n_up_input_symbols": len(up_genes),
        "n_down_input_symbols": len(down_genes),
        "n_up_entrez_ids": len(up_ids),
        "n_down_entrez_ids": len(down_ids),
        "n_up_bing_genes": up_bing,
        "n_down_bing_genes": down_bing,
        "up_unmapped_symbols": up_missing,
        "down_unmapped_symbols": down_missing,
        "api_key_present": True,
        "secret_recorded": False,
    }
    qc_output.write_text(json.dumps(qc, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"job_id": job_id, "status": job.get("status"), **qc}, sort_keys=True))


if __name__ == "__main__":
    main()
