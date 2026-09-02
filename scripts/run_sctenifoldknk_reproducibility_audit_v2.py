#!/usr/bin/env python3
"""Run only the prespecified versioned HNF4A/SOX4 scTenifoldKnk reruns."""

from __future__ import annotations

import argparse
import json
import os
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

try:
    from audit_sctenifoldknk_reproducibility_v2 import FORMAL_SEEDS, PROJECT_ROOT, sha256_file
except ModuleNotFoundError:
    from scripts.audit_sctenifoldknk_reproducibility_v2 import FORMAL_SEEDS, PROJECT_ROOT, sha256_file


DEFAULT_DATA_DIR = PROJECT_ROOT / "data/processed/driver/sctenifoldknk_reproducibility_audit_v2"
DEFAULT_METADATA_DIR = PROJECT_ROOT / "metadata/driver/sctenifoldknk_reproducibility_audit_v2"
DEFAULT_REPORT_DIR = PROJECT_ROOT / "reports/sctenifoldknk_reproducibility_audit_v2"
DEFAULT_RSCRIPT = Path(r"C:\Program Files\R\R-4.5.0\bin\x64\Rscript.exe")
RERUN_TFS = ["HNF4A", "SOX4"]


def select_targets_to_rerun(audit: dict[str, Any], requested_targets: list[str] | None = None) -> list[str]:
    requested = requested_targets or ["HNF4A", "EGR1", "SOX4"]
    historical = audit.get("historical_runs", {})
    return [tf for tf in requested if tf in RERUN_TFS and bool(historical.get(tf, {}).get("decision", {}).get("needs_rerun"))]


def build_rerun_jobs(
    project_root: Path,
    audit: dict[str, Any],
    targets: list[str],
    seeds: list[int] | None = None,
    data_dir: Path | None = None,
    metadata_dir: Path | None = None,
) -> list[dict[str, Any]]:
    seeds = seeds or FORMAL_SEEDS
    data_dir = data_dir or project_root / "data/processed/driver/sctenifoldknk_reproducibility_audit_v2"
    metadata_dir = metadata_dir or project_root / "metadata/driver/sctenifoldknk_reproducibility_audit_v2"
    definitions = {
        "HNF4A": {
            "axis": "Identity",
            "subset": "normal_reference",
            "input_dir": project_root / "data/processed/driver/figure2e_hnf4a_sctenifoldknk/normal_reference",
            "matrix_file": "figure2e_hnf4a_normal_reference_counts_genes_x_cells.mtx",
            "genes_file": "figure2e_hnf4a_normal_reference_genes.tsv",
        },
        "SOX4": {
            "axis": "Malignant state",
            "subset": "malignant_like",
            "input_dir": project_root / "data/processed/driver/sctenifoldknk_module7_1/malignant_like",
            "matrix_file": "sctenifoldknk_counts_genes_x_cells.mtx",
            "genes_file": "sctenifoldknk_genes.tsv",
        },
    }
    jobs = []
    for tf in targets:
        if tf not in definitions:
            continue
        definition = definitions[tf]
        for seed in seeds:
            jobs.append(
                {
                    "run_id": f"{tf}_{definition['subset']}_seed{seed}",
                    "target_tf": tf,
                    "axis": definition["axis"],
                    "subset": definition["subset"],
                    "seed": int(seed),
                    "input_dir": definition["input_dir"],
                    "matrix_file": definition["matrix_file"],
                    "genes_file": definition["genes_file"],
                    "output_dir": data_dir / tf / definition["subset"] / f"seed_{seed}",
                    "metadata_dir": metadata_dir / tf / definition["subset"] / f"seed_{seed}",
                    "nc_nNet": 10,
                    "nc_nCells": 500,
                    "nc_nComp": 3,
                    "ma_nDim": 2,
                    "ncores": 8,
                    "qc": "false",
                    "qc_min_cells": 3,
                }
            )
    return jobs


def _historical_audit_path(project_root: Path) -> Path:
    return project_root / "metadata/driver/sctenifoldknk_reproducibility_audit_v2/historical_audit.json"


def _rscript_command(rscript: Path, script: Path, job: dict[str, Any], project_root: Path = PROJECT_ROOT) -> list[str]:
    return _rscript_command_with_mount(rscript, script, job, project_root, None)


def _rscript_command_with_mount(
    rscript: Path,
    script: Path,
    job: dict[str, Any],
    project_root: Path,
    r_mount: str | None,
) -> list[str]:
    def r_path(path: Path) -> str:
        resolved = Path(path).resolve()
        if r_mount is None:
            return str(resolved)
        relative = os.path.relpath(resolved, Path(project_root).resolve())
        return os.path.join(r_mount, relative)

    return [
        str(rscript),
        r_path(script),
        "--input-dir",
        r_path(job["input_dir"]),
        "--matrix-file",
        str(job["matrix_file"]),
        "--genes-file",
        str(job["genes_file"]),
        "--output-dir",
        r_path(job["output_dir"]),
        "--metadata-dir",
        r_path(job["metadata_dir"]),
        "--target-tf",
        str(job["target_tf"]),
        "--axis",
        str(job["axis"]),
        "--subset",
        str(job["subset"]),
        "--seed",
        str(job["seed"]),
        "--nc-nnet",
        str(job["nc_nNet"]),
        "--nc-ncells",
        str(job["nc_nCells"]),
        "--nc-ncomp",
        str(job["nc_nComp"]),
        "--ma-ndim",
        str(job["ma_nDim"]),
        "--ncores",
        str(job["ncores"]),
        "--qc",
        str(job["qc"]),
        "--qc-min-cells",
        str(job["qc_min_cells"]),
    ]


def mount_project_drive(project_root: Path) -> tuple[str | None, bool]:
    if os.name != "nt":
        return None, False
    for letter in "RSTUVWXYZ":
        drive_root = f"{letter}:\\"
        if Path(drive_root).exists():
            continue
        completed = subprocess.run(["subst", f"{letter}:", str(Path(project_root).resolve())], capture_output=True, text=True, check=False)
        if completed.returncode == 0:
            return drive_root, True
    raise RuntimeError("No unused Windows drive letter is available for an ASCII R workspace mount")


def unmount_project_drive(r_mount: str | None, created: bool) -> None:
    if not created or not r_mount or os.name != "nt":
        return
    subprocess.run(["subst", r_mount[:2], "/D"], capture_output=True, text=True, check=False)


def _find_output(metadata_dir: Path, target_tf: str, subset: str, seed: int, suffix: str) -> Path:
    return metadata_dir / target_tf / subset / f"seed_{seed}" / f"{target_tf}_{subset}_seed{seed}_{suffix}"


def _log_diagnostics(path: Path) -> tuple[str, str]:
    if not path.exists():
        return "", ""
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    warning_lines = [line for line in lines if "warning" in line.lower() or "deprecated" in line.lower()]
    error_lines = [line for line in lines if "error" in line.lower() or "execution halted" in line.lower()]
    return " | ".join(dict.fromkeys(warning_lines)), " | ".join(dict.fromkeys(error_lines))


def run_jobs(
    project_root: Path,
    audit: dict[str, Any],
    rscript: Path = DEFAULT_RSCRIPT,
    targets: list[str] | None = None,
    seeds: list[int] | None = None,
    data_dir: Path = DEFAULT_DATA_DIR,
    metadata_dir: Path = DEFAULT_METADATA_DIR,
    report_dir: Path = DEFAULT_REPORT_DIR,
    dry_run: bool = False,
) -> dict[str, Any]:
    selected = select_targets_to_rerun(audit, targets)
    jobs = build_rerun_jobs(project_root, audit, selected, seeds=seeds, data_dir=data_dir, metadata_dir=metadata_dir)
    script = project_root / "scripts/run_sctenifoldknk_reproducibility_audit_v2.R"
    report_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)
    log_dir = report_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict[str, Any]] = []
    r_mount, mount_created = mount_project_drive(project_root)
    try:
        for job in jobs:
            job["output_dir"].mkdir(parents=True, exist_ok=True)
            job["metadata_dir"].mkdir(parents=True, exist_ok=True)
            command = _rscript_command_with_mount(rscript, script, job, project_root, r_mount)
            stdout_path = log_dir / f"{job['run_id']}_stdout.log"
            stderr_path = log_dir / f"{job['run_id']}_stderr.log"
            start = datetime.now(timezone.utc)
            status = "dry_run" if dry_run else "failed"
            exit_status = None
            if not dry_run:
                with stdout_path.open("w", encoding="utf-8", errors="replace") as stdout, stderr_path.open("w", encoding="utf-8", errors="replace") as stderr:
                    completed = subprocess.run(command, cwd=r_mount or project_root, stdout=stdout, stderr=stderr, check=False)
                exit_status = int(completed.returncode)
                status = "success" if exit_status == 0 else "failed"
            end = datetime.now(timezone.utc)
            report_path = _find_output(job["metadata_dir"].parent.parent.parent, job["target_tf"], job["subset"], job["seed"], "run_report.json")
            report = json.loads(report_path.read_text(encoding="utf-8")) if report_path.exists() else {}
            warning_text, error_text = _log_diagnostics(stderr_path)
            prefix = f"{job['target_tf']}_{job['subset']}_seed{job['seed']}"
            expected_outputs = {
                "result_rds": job["output_dir"] / f"{prefix}_result.rds",
                "perturbation_genes": job["metadata_dir"] / f"{prefix}_perturbation_genes.tsv",
                "network_adjacency_summary": job["metadata_dir"] / f"{prefix}_network_adjacency_summary.tsv",
                "result_contract": job["metadata_dir"] / f"{prefix}_result_contract.tsv",
                "run_report": report_path,
            }
            manifest_rows.append(
                {
                "run_id": job["run_id"],
                "target_tf": job["target_tf"],
                "axis": job["axis"],
                "subset": job["subset"],
                "seed": job["seed"],
                "status": status,
                "exit_status": exit_status,
                "started_at": start.isoformat(),
                "finished_at": end.isoformat(),
                "elapsed_seconds": (end - start).total_seconds(),
                "command": subprocess.list2cmdline(command),
                "script_path": str(script.resolve()),
                "script_sha256": sha256_file(script),
                "input_matrix_path": str((job["input_dir"] / job["matrix_file"]).resolve()),
                "input_matrix_sha256": sha256_file(job["input_dir"] / job["matrix_file"]),
                "input_genes_path": str((job["input_dir"] / job["genes_file"]).resolve()),
                "input_genes_sha256": sha256_file(job["input_dir"] / job["genes_file"]),
                "output_dir": str(job["output_dir"].resolve()),
                "metadata_dir": str(job["metadata_dir"].resolve()),
                "result_rds": str(expected_outputs["result_rds"].resolve()),
                "perturbation_genes": str(expected_outputs["perturbation_genes"].resolve()),
                "network_adjacency_summary": str(expected_outputs["network_adjacency_summary"].resolve()),
                "result_contract": str(expected_outputs["result_contract"].resolve()),
                "run_report": str(expected_outputs["run_report"].resolve()),
                "result_rds_sha256": sha256_file(expected_outputs["result_rds"]),
                "perturbation_genes_sha256": sha256_file(expected_outputs["perturbation_genes"]),
                "network_adjacency_summary_sha256": sha256_file(expected_outputs["network_adjacency_summary"]),
                "result_contract_sha256": sha256_file(expected_outputs["result_contract"]),
                "run_report_sha256": sha256_file(expected_outputs["run_report"]),
                "nc_nNet": job["nc_nNet"],
                "nc_nCells": job["nc_nCells"],
                "nc_nComp": job["nc_nComp"],
                "ma_nDim": job["ma_nDim"],
                "ncores": job["ncores"],
                "qc": job["qc"],
                "qc_min_cells": job["qc_min_cells"],
                "r_version": report.get("runtime", {}).get("r", ""),
                "scTenifoldKnk_version": report.get("runtime", {}).get("scTenifoldKnk", ""),
                "dependency_versions": json.dumps(report.get("runtime", {}), ensure_ascii=False),
                "hostname": platform.node(),
                "stdout_log": str(stdout_path.resolve()),
                "stderr_log": str(stderr_path.resolve()),
                "warnings": warning_text,
                "errors": error_text if status in {"success", "dry_run"} else (error_text or "see stderr_log"),
                }
            )
    finally:
        unmount_project_drive(r_mount, mount_created)
    manifest_path = metadata_dir / "run_manifest.tsv"
    new_manifest = pd.DataFrame(manifest_rows)
    if manifest_path.exists() and manifest_path.stat().st_size:
        existing_manifest = pd.read_csv(manifest_path, sep="\t")
        for column in ["result_rds", "perturbation_genes", "network_adjacency_summary", "result_contract", "run_report"]:
            if column not in existing_manifest.columns:
                continue
            existing_manifest[column] = existing_manifest[column].map(
                lambda value: str((project_root / str(value)[3:]).resolve())
                if isinstance(value, str) and value[:3].lower() in {"r:/", "r:\\"} else value
                )
        if "stderr_log" in existing_manifest.columns:
            diagnostics = existing_manifest["stderr_log"].map(lambda value: _log_diagnostics(Path(str(value))))
            existing_manifest["warnings"] = [item[0] for item in diagnostics]
            existing_manifest["errors"] = [item[1] for item in diagnostics]
        for path_column, hash_column in [
            ("result_rds", "result_rds_sha256"),
            ("perturbation_genes", "perturbation_genes_sha256"),
            ("network_adjacency_summary", "network_adjacency_summary_sha256"),
            ("result_contract", "result_contract_sha256"),
            ("run_report", "run_report_sha256"),
        ]:
            if path_column in existing_manifest.columns:
                existing_manifest[hash_column] = existing_manifest[path_column].map(
                    lambda value: sha256_file(Path(str(value))) if isinstance(value, str) else None
                )
        if "run_id" in existing_manifest.columns and not new_manifest.empty:
            existing_manifest = existing_manifest.loc[~existing_manifest["run_id"].isin(new_manifest["run_id"])]
            if existing_manifest.empty:
                new_manifest = new_manifest.reset_index(drop=True)
            else:
                new_manifest = pd.concat([existing_manifest, new_manifest], ignore_index=True)
        elif new_manifest.empty:
            new_manifest = existing_manifest
    new_manifest.to_csv(manifest_path, sep="\t", index=False)
    return {"targets": selected, "n_jobs": len(jobs), "manifest": str(manifest_path.resolve()), "rows": manifest_rows}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--audit", type=Path, default=None)
    parser.add_argument("--rscript", type=Path, default=DEFAULT_RSCRIPT)
    parser.add_argument("--targets", default="HNF4A,SOX4")
    parser.add_argument("--seeds", default=",".join(map(str, FORMAL_SEEDS)))
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA_DIR)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    audit_path = args.audit or _historical_audit_path(args.project_root)
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    result = run_jobs(
        project_root=args.project_root,
        audit=audit,
        rscript=args.rscript,
        targets=[x.strip().upper() for x in args.targets.split(",") if x.strip()],
        seeds=[int(x) for x in args.seeds.split(",") if x.strip()],
        data_dir=args.data_dir,
        metadata_dir=args.metadata_dir,
        report_dir=args.report_dir,
        dry_run=args.dry_run,
    )
    print(json.dumps({key: value for key, value in result.items() if key != "rows"}, ensure_ascii=False, indent=2))
    return 0 if all(row["status"] in {"success", "dry_run"} for row in result["rows"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())
