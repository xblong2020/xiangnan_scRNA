from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path
from typing import Any, Sequence

import pandas as pd
import requests

try:
    from scripts.clue_query_module9_8 import (
        list_jobs,
        map_symbols_to_entrez,
        safe_job_status,
    )
except ModuleNotFoundError:
    from clue_query_module9_8 import list_jobs, map_symbols_to_entrez, safe_job_status


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA_DIR = ROOT / "metadata/driver"
DEFAULT_PROFILES = DEFAULT_METADATA_DIR / "module9_9_landmark_decomposition_signature_long.tsv"
DEFAULT_GENE_INFO = DEFAULT_METADATA_DIR / "GSE92742_Broad_LINCS_gene_info.txt.gz"
DEFAULT_JOBS_URL = "https://api.clue.io/api/jobs"
QUERY_PROFILES = ["malignant_only", "rescue_only", "combined_balanced"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 9.9 CLUE landmark profile queries.")
    parser.add_argument("--profiles", type=Path, default=DEFAULT_PROFILES)
    parser.add_argument("--gene-info", type=Path, default=DEFAULT_GENE_INFO)
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA_DIR)
    parser.add_argument("--jobs-url", default=DEFAULT_JOBS_URL)
    parser.add_argument("--api-key-env", default="CLUE_API_KEY")
    parser.add_argument("--poll-seconds", type=float, default=10.0)
    parser.add_argument("--max-polls", type=int, default=120)
    parser.add_argument("--timeout-seconds", type=int, default=180)
    return parser.parse_args()


def gmt_text(query_name: str, gene_ids: Sequence[str]) -> str:
    return "\t".join([query_name, "module9_9_clue_entrez", *gene_ids]) + "\n"


def build_clue_profile_payload(
    profile: str,
    up_entrez: Sequence[str],
    down_entrez: Sequence[str],
) -> dict[str, Any]:
    if not up_entrez and not down_entrez:
        raise ValueError(f"{profile} has no mapped genes")
    query_name = f"module9_9_{profile}"
    payload: dict[str, Any] = {
        "tool_id": "sig_gutc_tool",
        "data_type": "L1000",
        "name": query_name,
        "dataset": "Touchstone",
        "ignoreWarnings": True,
        "es_tail": "both" if up_entrez and down_entrez else ("up" if up_entrez else "down"),
    }
    if up_entrez:
        payload["uptag-cmapfile"] = gmt_text(query_name, up_entrez).strip()
    if down_entrez:
        payload["dntag-cmapfile"] = gmt_text(query_name, down_entrez).strip()
    return payload


def submit_profile(
    jobs_url: str,
    api_key: str,
    payload: dict[str, Any],
    timeout_seconds: int,
) -> str:
    response = requests.post(
        jobs_url,
        headers={"user_key": api_key, "Accept": "application/json"},
        json=payload,
        timeout=timeout_seconds,
    )
    result = response.json()
    if response.status_code != 200:
        errors = {
            key: value
            for key, value in result.items()
            if key in {"system", "up", "down", "both"}
        }
        raise RuntimeError(
            f"CLUE submission failed for {payload['name']} "
            f"with HTTP {response.status_code}: {errors}"
        )
    job_id = result.get("result", {}).get("job_id")
    if not job_id:
        raise ValueError(f"CLUE response missing job_id for {payload['name']}")
    return str(job_id)


def wait_for_profiles(
    jobs_url: str,
    api_key: str,
    job_ids: dict[str, str],
    poll_seconds: float,
    max_polls: int,
    timeout_seconds: int,
) -> dict[str, dict[str, Any]]:
    pending = dict(job_ids)
    completed: dict[str, dict[str, Any]] = {}
    for _ in range(max_polls):
        jobs = list_jobs(jobs_url, api_key, timeout_seconds)
        jobs_by_id = {str(job.get("job_id")): job for job in jobs}
        for profile, job_id in list(pending.items()):
            job = jobs_by_id.get(job_id)
            if job is None:
                continue
            status = str(job.get("status", "")).lower()
            if status in {"completed", "failed", "error"}:
                completed[profile] = job
                pending.pop(profile)
        if not pending:
            return completed
        time.sleep(poll_seconds)
    raise TimeoutError(f"CLUE profiles did not finish: {sorted(pending)}")


def download_archive(job: dict[str, Any], destination: Path, timeout_seconds: int) -> None:
    url = str(job.get("download_url", ""))
    if not url:
        raise ValueError(f"job {job.get('job_id')} missing download_url")
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

    profiles = pd.read_csv(args.profiles, sep="\t")
    gene_info = pd.read_csv(args.gene_info, sep="\t", compression="infer")
    job_ids: dict[str, str] = {}
    qc_rows: list[dict[str, object]] = []
    for profile in QUERY_PROFILES:
        frame = profiles.loc[profiles["profile"].eq(profile)].copy()
        up_symbols = frame.loc[frame["v_score"].gt(0), "gene"].astype(str).tolist()
        down_symbols = frame.loc[frame["v_score"].lt(0), "gene"].astype(str).tolist()
        up_ids, up_missing, up_bing = map_symbols_to_entrez(up_symbols, gene_info)
        down_ids, down_missing, down_bing = map_symbols_to_entrez(down_symbols, gene_info)
        payload = build_clue_profile_payload(profile, up_ids, down_ids)
        query_name = payload["name"]
        if up_ids:
            (args.metadata_dir / f"{query_name}_up_entrez.gmt").write_text(
                gmt_text(query_name, up_ids), encoding="utf-8"
            )
        if down_ids:
            (args.metadata_dir / f"{query_name}_down_entrez.gmt").write_text(
                gmt_text(query_name, down_ids), encoding="utf-8"
            )
        query_status = "not_submitted_clue_requires_paired_sets"
        job_id = ""
        if up_ids and down_ids:
            job_id = submit_profile(
                args.jobs_url, api_key, payload, args.timeout_seconds
            )
            job_ids[profile] = job_id
            query_status = "submitted"
        qc_rows.append(
            {
                "profile": profile,
                "n_up_symbols": len(up_symbols),
                "n_down_symbols": len(down_symbols),
                "n_up_entrez": len(up_ids),
                "n_down_entrez": len(down_ids),
                "n_up_bing": up_bing,
                "n_down_bing": down_bing,
                "up_unmapped": ",".join(up_missing),
                "down_unmapped": ",".join(down_missing),
                "job_id": job_id,
                "query_status": query_status,
            }
        )

    if not job_ids:
        raise ValueError("no paired CLUE profile was available for submission")
    completed = wait_for_profiles(
        args.jobs_url,
        api_key,
        job_ids,
        args.poll_seconds,
        args.max_polls,
        args.timeout_seconds,
    )
    safe_statuses: dict[str, Any] = {}
    for profile, job in completed.items():
        if str(job.get("status", "")).lower() != "completed":
            raise RuntimeError(f"CLUE {profile} ended with status {job.get('status')}")
        archive = args.metadata_dir / f"module9_9_clue_{profile}_results.tar.gz"
        download_archive(job, archive, args.timeout_seconds)
        safe_statuses[profile] = safe_job_status(job, archive)

    pd.DataFrame(qc_rows).to_csv(
        args.metadata_dir / "module9_9_clue_profile_qc.tsv", sep="\t", index=False
    )
    (args.metadata_dir / "module9_9_clue_jobs.json").write_text(
        json.dumps(
            {
                "status": "completed",
                "profiles": safe_statuses,
                "api_key_present": True,
                "secret_recorded": False,
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "status": "completed",
                "job_ids": job_ids,
                "secret_recorded": False,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
