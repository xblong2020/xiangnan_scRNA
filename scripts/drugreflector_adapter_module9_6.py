from __future__ import annotations

import argparse
import hashlib
import json
import platform
from datetime import datetime, timezone
import importlib
from importlib.metadata import PackageNotFoundError, version
from importlib.util import find_spec
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA_DIR = ROOT / "metadata/driver"
DEFAULT_SIGNATURE = DEFAULT_METADATA_DIR / "module9_4_drug_reversal_signature.tsv"
DEFAULT_CHECKPOINT_DIRS = [
    DEFAULT_METADATA_DIR / "drugreflector_checkpoints",
    ROOT / "data" / "processed" / "drugreflector_checkpoints",
    ROOT / "models" / "drugreflector",
]
OUTPUT_STEM = "module9_6_drugreflector"
EXPECTED_CHECKPOINT_COUNT = 3
PACKAGE_DISTRIBUTIONS = {
    "torch": ("torch",),
    "drugreflector": ("drugreflector",),
    "zenodo_get": ("zenodo-get", "zenodo_get"),
}
EXPECTED_CHECKPOINTS = {
    "model_fold_0.pt": {
        "size_bytes": 91355041,
        "md5": "0a27e253713c37f4874318b5ba0c27a9",
    },
    "model_fold_1.pt": {
        "size_bytes": 91355041,
        "md5": "0e785196fd046d946f84e4480c81ff53",
    },
    "model_fold_2.pt": {
        "size_bytes": 91355041,
        "md5": "d8e36f6a8f9fa7a22feda7acdd0bee86",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 9.6 DrugReflector adapter audit.")
    parser.add_argument("--signature", type=Path, default=DEFAULT_SIGNATURE)
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA_DIR)
    parser.add_argument("--primary-only", action="store_true")
    parser.add_argument("--seed", type=int, default=20260616)
    return parser.parse_args()


def package_version(name: str) -> str:
    for distribution in PACKAGE_DISTRIBUTIONS.get(name, (name,)):
        try:
            return version(distribution)
        except PackageNotFoundError:
            continue
    return "not_installed"


def read_tsv_or_empty(path: Path, **kwargs) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", **kwargs)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def normalize_gene_symbol(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip().upper()


def display_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def package_presence_map() -> dict[str, bool]:
    return {name: bool(find_spec(name)) for name in ["torch", "drugreflector", "zenodo_get"]}


def package_import_results() -> dict[str, dict[str, object]]:
    results: dict[str, dict[str, object]] = {}
    for name in ["torch", "drugreflector", "zenodo_get"]:
        if not find_spec(name):
            results[name] = {"import_ok": False, "error": "not_installed"}
            continue
        try:
            importlib.import_module(name)
            results[name] = {"import_ok": True, "error": ""}
        except Exception as exc:  # pragma: no cover - depends on local runtime DLLs/packages.
            results[name] = {"import_ok": False, "error": f"{type(exc).__name__}: {exc}"}
    return results


def build_vscore_frame(signature: pd.DataFrame, primary_only: bool) -> pd.DataFrame:
    required = {
        "gene",
        "desired_direction",
        "component",
        "final_weight",
        "include_primary",
        "include_sensitivity",
        "conflict_flag",
        "housekeeping_or_qc_flag",
        "source_file",
        "source_metric",
    }
    missing = required.difference(signature.columns)
    if missing:
        raise ValueError(f"signature missing required columns: {sorted(missing)}")

    frame = signature.copy()
    frame["gene"] = frame["gene"].map(normalize_gene_symbol)
    frame = frame.loc[frame["gene"].ne("")].copy()
    frame["final_weight"] = pd.to_numeric(frame["final_weight"], errors="coerce").fillna(0.0)
    frame["include_primary"] = frame["include_primary"].astype(bool)
    frame["include_sensitivity"] = frame["include_sensitivity"].astype(bool)
    frame["conflict_flag"] = frame["conflict_flag"].astype(bool)
    frame["housekeeping_or_qc_flag"] = frame["housekeeping_or_qc_flag"].astype(bool)
    frame["desired_direction"] = frame["desired_direction"].astype(str)

    keep_mask = frame["include_primary"] if primary_only else frame["include_sensitivity"]
    frame = frame.loc[keep_mask].copy()
    if primary_only:
        frame = frame.loc[~frame["conflict_flag"] & ~frame["housekeeping_or_qc_flag"]].copy()

    sign = frame["desired_direction"].map({"up": 1.0, "down": -1.0}).fillna(0.0)
    frame["v_score"] = sign * frame["final_weight"]
    frame = frame[
        [
            "gene",
            "v_score",
            "desired_direction",
            "component",
            "final_weight",
            "include_primary",
            "include_sensitivity",
            "conflict_flag",
            "housekeeping_or_qc_flag",
            "source_file",
            "source_metric",
        ]
    ].copy()
    return frame.sort_values(["v_score", "gene"], ascending=[False, True]).reset_index(drop=True)


def find_checkpoint_files(candidate_dirs: Sequence[Path]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    patterns = ("*.pt", "*.pth", "*.ckpt")
    for base_dir in candidate_dirs:
        if not base_dir.exists():
            continue
        seen: set[Path] = set()
        for pattern in patterns:
            for path in base_dir.rglob(pattern):
                resolved = path.resolve()
                if resolved in seen or not path.is_file():
                    continue
                seen.add(resolved)
                rows.append(
                    {
                        "path": display_path(resolved),
                        "size_bytes": int(path.stat().st_size),
                        "parent_dir": display_path(base_dir.resolve()),
                    }
                )
    return sorted(rows, key=lambda row: (str(row["parent_dir"]), str(row["path"])))


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def audit_expected_checkpoints(candidate_dirs: Sequence[Path]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for filename, expected in EXPECTED_CHECKPOINTS.items():
        found_path: Path | None = None
        for base_dir in candidate_dirs:
            candidate = base_dir / filename
            if candidate.is_file():
                found_path = candidate
                break

        if found_path is None:
            rows.append(
                {
                    "file_name": filename,
                    "present": False,
                    "path": "",
                    "size_bytes": 0,
                    "expected_size_bytes": int(expected["size_bytes"]),
                    "size_ok": False,
                    "md5": "",
                    "expected_md5": str(expected["md5"]),
                    "md5_ok": False,
                    "valid": False,
                }
            )
            continue

        size_bytes = int(found_path.stat().st_size)
        actual_md5 = md5sum(found_path)
        size_ok = size_bytes == int(expected["size_bytes"])
        md5_ok = actual_md5 == str(expected["md5"])
        rows.append(
            {
                "file_name": filename,
                "present": True,
                "path": display_path(found_path.resolve()),
                "size_bytes": size_bytes,
                "expected_size_bytes": int(expected["size_bytes"]),
                "size_ok": size_ok,
                "md5": actual_md5,
                "expected_md5": str(expected["md5"]),
                "md5_ok": md5_ok,
                "valid": size_ok and md5_ok,
            }
        )
    return rows


def audit_drugreflector_dependencies(
    candidate_checkpoint_dirs: Sequence[Path] | None = None,
    package_presence: dict[str, bool] | None = None,
    import_results: dict[str, dict[str, object]] | None = None,
) -> dict[str, object]:
    dirs = list(candidate_checkpoint_dirs or DEFAULT_CHECKPOINT_DIRS)
    presence = dict(package_presence or package_presence_map())
    imports = dict(import_results or package_import_results())
    checkpoint_files = find_checkpoint_files(dirs)
    expected_checkpoints = audit_expected_checkpoints(dirs)
    packages = {}
    for name in ["torch", "drugreflector", "zenodo_get"]:
        present = bool(presence.get(name, False))
        import_result = imports.get(name, {"import_ok": False, "error": "not_checked"})
        packages[name] = {
            "present": present,
            "version": package_version(name) if present else "not_installed",
            "import_ok": bool(import_result.get("import_ok", False)),
            "import_error": str(import_result.get("error", "")),
        }
    checkpoint_summary = {
        "candidate_dirs": [display_path(path) for path in dirs],
        "n_checkpoint_files": len(checkpoint_files),
        "n_valid_expected_checkpoints": int(sum(bool(row["valid"]) for row in expected_checkpoints)),
        "expected_checkpoint_count": EXPECTED_CHECKPOINT_COUNT,
        "all_required_present": all(bool(row["valid"]) for row in expected_checkpoints),
        "checkpoint_files": checkpoint_files,
        "expected_checkpoints": expected_checkpoints,
    }
    package_present = all(packages[name]["present"] for name in packages)
    runtime_ready = all(packages[name]["present"] and packages[name]["import_ok"] for name in packages)
    status = "adapter_ready_model_missing"
    if package_present and not runtime_ready:
        status = "adapter_ready_runtime_blocked"
    elif runtime_ready and checkpoint_summary["all_required_present"]:
        status = "adapter_ready_model_available"
    return {
        "status": status,
        "packages": packages,
        "checkpoint_summary": checkpoint_summary,
    }


def summarize_vscore_frame(frame: pd.DataFrame, label: str) -> list[dict[str, object]]:
    if frame.empty:
        return [
            {"scope": label, "metric": "n_genes", "value": 0, "status": "empty"},
        ]
    down = int(frame["desired_direction"].eq("down").sum())
    up = int(frame["desired_direction"].eq("up").sum())
    return [
        {"scope": label, "metric": "n_genes", "value": int(len(frame)), "status": "reported"},
        {"scope": label, "metric": "n_up_genes", "value": up, "status": "reported"},
        {"scope": label, "metric": "n_down_genes", "value": down, "status": "reported"},
        {"scope": label, "metric": "v_score_min", "value": float(frame["v_score"].min()), "status": "reported"},
        {"scope": label, "metric": "v_score_max", "value": float(frame["v_score"].max()), "status": "reported"},
        {"scope": label, "metric": "v_score_mean", "value": float(frame["v_score"].mean()), "status": "reported"},
        {"scope": label, "metric": "n_duplicate_gene_rows", "value": int(frame["gene"].duplicated().sum()), "status": "reported"},
        {"scope": label, "metric": "n_conflict_rows", "value": int(frame["conflict_flag"].sum()), "status": "reported"},
        {"scope": label, "metric": "n_housekeeping_qc_rows", "value": int(frame["housekeeping_or_qc_flag"].sum()), "status": "reported"},
    ]


def build_qc_table(signature: pd.DataFrame, primary: pd.DataFrame, sensitivity: pd.DataFrame, dependency_audit: dict[str, object]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    rows.extend(summarize_vscore_frame(primary, "primary"))
    rows.extend(summarize_vscore_frame(sensitivity, "sensitivity"))
    rows.extend(
        [
            {"scope": "signature", "metric": "n_signature_rows", "value": int(len(signature)), "status": "reported"},
            {"scope": "signature", "metric": "n_signature_genes", "value": int(signature["gene"].nunique()), "status": "reported"},
            {"scope": "signature", "metric": "n_include_primary", "value": int(signature["include_primary"].astype(bool).sum()), "status": "reported"},
            {"scope": "signature", "metric": "n_include_sensitivity", "value": int(signature["include_sensitivity"].astype(bool).sum()), "status": "reported"},
            {"scope": "dependency", "metric": "adapter_status", "value": dependency_audit["status"], "status": "reported"},
            {"scope": "dependency", "metric": "n_checkpoint_files", "value": dependency_audit["checkpoint_summary"]["n_checkpoint_files"], "status": "reported"},
            {"scope": "dependency", "metric": "n_valid_expected_checkpoints", "value": dependency_audit["checkpoint_summary"]["n_valid_expected_checkpoints"], "status": "reported"},
        ]
    )
    return pd.DataFrame(rows)


def build_runbook(
    report_status: str,
    signature_path: Path,
    primary_path: Path,
    sensitivity_path: Path,
    dependency_audit_path: Path,
) -> str:
    return "\n".join(
        [
            "# Module 9.6 DrugReflector Runbook",
            "",
            f"Current adapter status: `{report_status}`.",
            "",
            "## Current inputs",
            f"- Signature source: `{display_path(signature_path)}`",
            f"- Primary v-score: `{display_path(primary_path)}`",
            f"- Sensitivity v-score: `{display_path(sensitivity_path)}`",
            f"- Dependency audit: `{display_path(dependency_audit_path)}`",
            "",
            "## What is missing now",
            "- `torch`, `drugreflector`, and/or `zenodo_get` must be installed and importable in the selected Python runtime.",
            "- DrugReflector model checkpoints must match the expected file names, sizes, and MD5 checksums.",
            "",
            "## Recommended next execution steps",
            "1. Create or activate a dedicated DrugReflector-capable Python environment.",
            "2. Install runtime dependencies. Example command:",
            "   `python -m pip install torch zenodo-get git+https://github.com/Cellarity/drugreflector.git`",
            "3. Download the three DrugReflector checkpoints from Zenodo DOI `10.5281/zenodo.16912444` into one of:",
            f"   `{display_path(DEFAULT_CHECKPOINT_DIRS[0])}`",
            f"   `{display_path(DEFAULT_CHECKPOINT_DIRS[1])}`",
            f"   `{display_path(DEFAULT_CHECKPOINT_DIRS[2])}`",
            "4. Re-run this adapter audit to confirm `adapter_ready_model_available`.",
            "5. Use the primary v-score table as the first inference input; keep the sensitivity table for robustness reruns.",
            "",
            "## Suggested real inference handoff",
            "- Input orientation: positive scores indicate desired rescue activation; negative scores indicate desired malignant-program suppression.",
            "- Start with the primary v-score file, then compare top candidates against the sensitivity v-score file to assess robustness.",
            "- Keep Module 9.5 L1000FWD/CLUE payloads alongside DrugReflector outputs for cross-method triangulation.",
        ]
    ) + "\n"


def main() -> None:
    args = parse_args()
    np.random.seed(args.seed)
    metadata_dir = args.metadata_dir
    metadata_dir.mkdir(parents=True, exist_ok=True)

    signature = read_tsv_or_empty(args.signature)
    if signature.empty:
        raise FileNotFoundError(f"signature not found or empty: {args.signature}")

    primary = build_vscore_frame(signature, primary_only=True)
    sensitivity = build_vscore_frame(signature, primary_only=False)
    dependency_audit = audit_drugreflector_dependencies()

    outputs = {
        "vscore_primary": metadata_dir / f"{OUTPUT_STEM}_vscore_primary.tsv",
        "vscore_sensitivity": metadata_dir / f"{OUTPUT_STEM}_vscore_sensitivity.tsv",
        "input_qc": metadata_dir / f"{OUTPUT_STEM}_input_qc.tsv",
        "dependency_audit": metadata_dir / f"{OUTPUT_STEM}_dependency_audit.json",
        "runbook": metadata_dir / f"{OUTPUT_STEM}_runbook.md",
        "report": metadata_dir / f"{OUTPUT_STEM}_report.json",
    }

    primary.to_csv(outputs["vscore_primary"], sep="\t", index=False)
    sensitivity.to_csv(outputs["vscore_sensitivity"], sep="\t", index=False)
    qc = build_qc_table(signature, primary, sensitivity, dependency_audit)
    qc.to_csv(outputs["input_qc"], sep="\t", index=False)
    outputs["dependency_audit"].write_text(json.dumps(dependency_audit, indent=2, sort_keys=True), encoding="utf-8")
    outputs["runbook"].write_text(
        build_runbook(
            dependency_audit["status"],
            args.signature,
            outputs["vscore_primary"],
            outputs["vscore_sensitivity"],
            outputs["dependency_audit"],
        ),
        encoding="utf-8",
    )

    report = {
        "module": "module9_6_drugreflector_adapter",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": dependency_audit["status"],
        "inputs": {
            "signature": str(args.signature.resolve()),
            "primary_only_flag": bool(args.primary_only),
            "seed": args.seed,
        },
        "outputs": {key: str(path.resolve()) for key, path in outputs.items()},
        "summary": {
            "n_primary_vscore_genes": int(len(primary)),
            "n_sensitivity_vscore_genes": int(len(sensitivity)),
            "n_primary_up": int(primary["desired_direction"].eq("up").sum()),
            "n_primary_down": int(primary["desired_direction"].eq("down").sum()),
            "n_sensitivity_conflict_rows": int(sensitivity["conflict_flag"].sum()),
            "n_sensitivity_housekeeping_qc_rows": int(sensitivity["housekeeping_or_qc_flag"].sum()),
            "n_checkpoint_files": int(dependency_audit["checkpoint_summary"]["n_checkpoint_files"]),
            "n_valid_expected_checkpoints": int(dependency_audit["checkpoint_summary"]["n_valid_expected_checkpoints"]),
        },
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "pandas": package_version("pandas"),
            "numpy": package_version("numpy"),
            "torch": package_version("torch") if dependency_audit["packages"]["torch"]["present"] else "not_installed",
            "drugreflector": package_version("drugreflector") if dependency_audit["packages"]["drugreflector"]["present"] else "not_installed",
            "zenodo_get": package_version("zenodo_get") if dependency_audit["packages"]["zenodo_get"]["present"] else "not_installed",
        },
    }
    outputs["report"].write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(report["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
