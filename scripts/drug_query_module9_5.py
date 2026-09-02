from __future__ import annotations

import argparse
import json
import os
import platform
import re
import time
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Sequence
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA_DIR = ROOT / "metadata/driver"
DEFAULT_UP_GMT = DEFAULT_METADATA_DIR / "module9_4_drug_reversal_up.gmt"
DEFAULT_DOWN_GMT = DEFAULT_METADATA_DIR / "module9_4_drug_reversal_down.gmt"
DEFAULT_L1000FWD_BASE_URL = "https://maayanlab.cloud/l1000fwd"
DEFAULT_CLUE_JOBS_URL = "https://api.clue.io/api/jobs"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 9.5 external drug signature query adapter.")
    parser.add_argument("--up-gmt", type=Path, default=DEFAULT_UP_GMT)
    parser.add_argument("--down-gmt", type=Path, default=DEFAULT_DOWN_GMT)
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA_DIR)
    parser.add_argument("--query-name", default="module9_5_drug_reversal_query")
    parser.add_argument("--l1000fwd-base-url", default=DEFAULT_L1000FWD_BASE_URL)
    parser.add_argument("--clue-jobs-url", default=DEFAULT_CLUE_JOBS_URL)
    parser.add_argument("--skip-l1000fwd", action="store_true")
    parser.add_argument("--submit-clue", action="store_true")
    parser.add_argument("--clue-api-key-env", default="CLUE_API_KEY")
    parser.add_argument("--timeout-seconds", type=int, default=60)
    parser.add_argument("--poll-seconds", type=float, default=5.0)
    parser.add_argument("--max-clue-polls", type=int, default=24)
    return parser.parse_args()


def package_version(name: str) -> str:
    try:
        return version(name)
    except PackageNotFoundError:
        return "not_installed"


def parse_gmt(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8").strip()
    if not text:
        raise ValueError(f"empty GMT file: {path}")
    parts = text.split("\t")
    if len(parts) < 3:
        raise ValueError(f"GMT must contain name, description, and at least one gene: {path}")
    genes = [gene.strip().upper() for gene in parts[2:] if gene.strip()]
    return {"name": parts[0], "description": parts[1], "genes": list(dict.fromkeys(genes))}


def build_l1000fwd_payload(up_genes: Sequence[str], down_genes: Sequence[str]) -> dict[str, list[str]]:
    return {
        "up_genes": [str(gene).strip().upper() for gene in up_genes if str(gene).strip()],
        "down_genes": [str(gene).strip().upper() for gene in down_genes if str(gene).strip()],
    }


def stringified_gmt(name: str, genes: Sequence[str]) -> str:
    return "\t".join([name, "module9_5_drug_reversal_query", *[str(gene).strip().upper() for gene in genes if str(gene).strip()]])


def build_clue_query_payload(up_genes: Sequence[str], down_genes: Sequence[str], query_name: str) -> dict[str, Any]:
    return {
        "tool_id": "sig_query",
        "name": query_name,
        "uptag": stringified_gmt(f"{query_name}_up", up_genes),
        "dntag": stringified_gmt(f"{query_name}_down", down_genes),
    }


def json_request(
    url: str,
    method: str = "GET",
    payload: dict[str, Any] | None = None,
    headers: dict[str, str] | None = None,
    timeout_seconds: int = 60,
) -> Any:
    body = None if payload is None else json.dumps(payload).encode("utf-8")
    request_headers = {"Content-Type": "application/json"}
    if headers:
        request_headers.update(headers)
    request = Request(url, data=body, headers=request_headers, method=method)
    with urlopen(request, timeout=timeout_seconds) as response:
        data = response.read().decode("utf-8")
        if not data:
            return {}
        return json.loads(data)


def l1000fwd_sig_search(base_url: str, payload: dict[str, list[str]], timeout_seconds: int) -> dict[str, Any]:
    base = base_url.rstrip("/")
    result = json_request(f"{base}/sig_search", method="POST", payload=payload, timeout_seconds=timeout_seconds)
    if isinstance(result, dict):
        return result
    raise ValueError("L1000FWD sig_search returned a non-object response")


def extract_l1000fwd_result_id(response: dict[str, Any]) -> str:
    for key in ["result_id", "id"]:
        value = response.get(key)
        if value:
            return str(value)
    if len(response) == 1:
        value = next(iter(response.values()))
        if value:
            return str(value)
    raise ValueError(f"could not find L1000FWD result id in response keys: {sorted(response)}")


def l1000fwd_topn(base_url: str, result_id: str, timeout_seconds: int) -> dict[str, Any]:
    base = base_url.rstrip("/")
    result = json_request(f"{base}/result/topn/{result_id}", method="GET", timeout_seconds=timeout_seconds)
    if isinstance(result, dict):
        return result
    raise ValueError("L1000FWD topn returned a non-object response")


def normalize_l1000fwd_results(topn: dict[str, Any]) -> pd.DataFrame:
    def parse_sig_id(sig_id: object) -> dict[str, object]:
        text = "" if sig_id is None else str(sig_id)
        match = re.search(
            r"^(?P<batch>[^_]+)_(?P<cell>[^_]+)_(?P<time>[^:]+):"
            r"(?P<pert_id>BRD-[A-Z]\d+)(?:-[^:]+)?:?(?P<dose>[^:]*)$",
            text,
        )
        if not match:
            return {"cell_line": None, "time": None, "pert_id": None, "dose": None}
        return {
            "cell_line": match.group("cell"),
            "time": match.group("time"),
            "pert_id": match.group("pert_id"),
            "dose": match.group("dose") or None,
        }

    rows = []
    for direction, label, sign in [
        ("similar", "similar_to_reversal_signature", 1.0),
        ("opposite", "opposite_to_reversal_signature", -1.0),
    ]:
        entries = topn.get(direction, [])
        if isinstance(entries, dict):
            entries = list(entries.values())
        for rank, entry in enumerate(entries, start=1):
            if not isinstance(entry, dict):
                continue
            sig_id = entry.get("sig_id") or entry.get("signature_id") or entry.get("id")
            parsed_sig = parse_sig_id(sig_id)
            score = pd.to_numeric(pd.Series([entry.get("score")]), errors="coerce").iloc[0]
            if pd.isna(score):
                score = pd.to_numeric(pd.Series([entry.get("similarity")]), errors="coerce").iloc[0]
            if pd.isna(score):
                score = pd.to_numeric(pd.Series([entry.get("scores")]), errors="coerce").iloc[0]
            magnitude = abs(float(score)) if pd.notna(score) else 0.0
            rows.append(
                {
                    "source_database": "L1000FWD",
                    "candidate_direction": label,
                    "result_group": direction,
                    "rank_within_group": rank,
                    "sig_id": sig_id,
                    "compound_id": entry.get("pert_id") or parsed_sig["pert_id"],
                    "compound_name": entry.get("pert_desc") or entry.get("perturbagen") or entry.get("compound"),
                    "cell_line": entry.get("cell_id") or parsed_sig["cell_line"],
                    "dose": entry.get("pert_dose") or parsed_sig["dose"],
                    "time": entry.get("pert_time") or parsed_sig["time"],
                    "raw_score": score,
                    "p_value": entry.get("pvals") if entry.get("pvals") is not None else entry.get("p_value"),
                    "q_value": entry.get("qvals") if entry.get("qvals") is not None else entry.get("q_value"),
                    "final_rank_score": sign * magnitude,
                    "evidence_notes": json.dumps(entry, sort_keys=True),
                }
            )
    if not rows:
        return pd.DataFrame(
            columns=[
                "source_database",
                "candidate_direction",
                "result_group",
                "rank_within_group",
                "sig_id",
                "compound_id",
                "compound_name",
                "cell_line",
                "dose",
                "time",
                "raw_score",
                "p_value",
                "q_value",
                "final_rank_score",
                "evidence_notes",
            ]
        )
    return pd.DataFrame(rows).sort_values(["final_rank_score", "result_group", "rank_within_group"], ascending=[False, True, True])


def run_l1000fwd_query(
    base_url: str,
    payload: dict[str, list[str]],
    timeout_seconds: int,
) -> tuple[dict[str, Any], dict[str, Any], pd.DataFrame]:
    search_response = l1000fwd_sig_search(base_url, payload, timeout_seconds)
    result_id = extract_l1000fwd_result_id(search_response)
    topn = l1000fwd_topn(base_url, result_id, timeout_seconds)
    normalized = normalize_l1000fwd_results(topn)
    return {"result_id": result_id, "search_response": search_response}, topn, normalized


def run_clue_query(
    jobs_url: str,
    payload: dict[str, Any],
    api_key: str,
    timeout_seconds: int,
    poll_seconds: float,
    max_polls: int,
) -> dict[str, Any]:
    submit = json_request(
        jobs_url,
        method="POST",
        payload=payload,
        headers={"user_key": api_key},
        timeout_seconds=timeout_seconds,
    )
    job_id = None
    if isinstance(submit, dict):
        job_id = submit.get("id") or submit.get("job_id")
    status_payload = {"submit_response": submit, "job_id": job_id}
    if not job_id:
        status_payload["status"] = "submitted_no_job_id"
        return status_payload
    status_url = f"{jobs_url.rstrip('/')}/{job_id}"
    last_status = None
    for _ in range(max_polls):
        time.sleep(poll_seconds)
        last_status = json_request(status_url, method="GET", headers={"user_key": api_key}, timeout_seconds=timeout_seconds)
        if isinstance(last_status, dict) and str(last_status.get("status", "")).lower() in {"done", "success", "completed", "failed", "error"}:
            break
    status_payload["status_response"] = last_status
    return status_payload


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def main() -> None:
    args = parse_args()
    metadata_dir = args.metadata_dir
    metadata_dir.mkdir(parents=True, exist_ok=True)

    up = parse_gmt(args.up_gmt)
    down = parse_gmt(args.down_gmt)
    l1000_payload = build_l1000fwd_payload(up["genes"], down["genes"])
    clue_payload = build_clue_query_payload(up["genes"], down["genes"], query_name=args.query_name)

    outputs = {
        "l1000_payload": metadata_dir / "module9_5_l1000fwd_query_payload.json",
        "l1000_search": metadata_dir / "module9_5_l1000fwd_search_response.json",
        "l1000_topn": metadata_dir / "module9_5_l1000fwd_topn_response.json",
        "l1000_candidates": metadata_dir / "module9_5_l1000fwd_candidate_ranking.tsv",
        "clue_payload": metadata_dir / "module9_5_clue_query_payload.json",
        "clue_status": metadata_dir / "module9_5_clue_query_status.json",
        "report": metadata_dir / "module9_5_report.json",
    }
    write_json(outputs["l1000_payload"], l1000_payload)
    write_json(outputs["clue_payload"], clue_payload)

    l1000_status: dict[str, Any] = {"status": "skipped"}
    topn: dict[str, Any] = {}
    candidates = normalize_l1000fwd_results({})
    if not args.skip_l1000fwd:
        try:
            l1000_status, topn, candidates = run_l1000fwd_query(args.l1000fwd_base_url, l1000_payload, args.timeout_seconds)
            l1000_status["status"] = "completed"
        except (HTTPError, URLError, TimeoutError, ValueError, OSError) as exc:
            l1000_status = {"status": "failed", "error_type": type(exc).__name__, "error": str(exc)}
    write_json(outputs["l1000_search"], l1000_status)
    write_json(outputs["l1000_topn"], topn)
    candidates.to_csv(outputs["l1000_candidates"], sep="\t", index=False)

    clue_key = os.environ.get(args.clue_api_key_env, "")
    clue_status: dict[str, Any] = {
        "status": "skipped_missing_api_key" if not clue_key else "skipped_submit_flag_not_set",
        "api_key_env": args.clue_api_key_env,
        "api_key_present": bool(clue_key),
        "secret_recorded": False,
    }
    if args.submit_clue and clue_key:
        try:
            clue_status = run_clue_query(
                args.clue_jobs_url,
                clue_payload,
                clue_key,
                args.timeout_seconds,
                args.poll_seconds,
                args.max_clue_polls,
            )
            clue_status.update({"status": "submitted", "api_key_env": args.clue_api_key_env, "api_key_present": True, "secret_recorded": False})
        except (HTTPError, URLError, TimeoutError, ValueError, OSError) as exc:
            clue_status = {
                "status": "failed",
                "api_key_env": args.clue_api_key_env,
                "api_key_present": True,
                "secret_recorded": False,
                "error_type": type(exc).__name__,
                "error": str(exc),
            }
    write_json(outputs["clue_status"], clue_status)

    report = {
        "module": "module9_5_drug_query",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "up_gmt": str(args.up_gmt.resolve()),
            "down_gmt": str(args.down_gmt.resolve()),
            "n_up_genes": len(up["genes"]),
            "n_down_genes": len(down["genes"]),
        },
        "outputs": {key: str(path.resolve()) for key, path in outputs.items()},
        "l1000fwd": {
            "base_url": args.l1000fwd_base_url,
            "status": l1000_status.get("status"),
            "result_id": l1000_status.get("result_id"),
            "n_candidates": int(len(candidates)),
            "n_similar": int(candidates["result_group"].eq("similar").sum()) if "result_group" in candidates.columns else 0,
            "n_opposite": int(candidates["result_group"].eq("opposite").sum()) if "result_group" in candidates.columns else 0,
        },
        "clue": clue_status,
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "pandas": package_version("pandas"),
        },
    }
    write_json(outputs["report"], report)
    print(json.dumps({"l1000fwd": report["l1000fwd"], "clue_status": clue_status.get("status")}, sort_keys=True))


if __name__ == "__main__":
    main()
