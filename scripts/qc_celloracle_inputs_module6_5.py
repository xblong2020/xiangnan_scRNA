from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import anndata as ad
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_H5AD = PROJECT_ROOT / "data/processed/driver/driver_cistarget_regulon_activity.module6_3c.h5ad"
DEFAULT_TF_LIST = PROJECT_ROOT / "metadata/driver/celloracle_input_tfs.module6_4.txt"
DEFAULT_SELECTION = PROJECT_ROOT / "metadata/driver/celloracle_tf_selection.module6_4.tsv"
DEFAULT_BASE_GRN = PROJECT_ROOT / "metadata/driver/scenic_resources/celloracle_hg38_promoter_base_grn.parquet"
DEFAULT_QC_TSV = PROJECT_ROOT / "metadata/driver/celloracle_module6_5_input_qc.tsv"
DEFAULT_QC_JSON = PROJECT_ROOT / "metadata/driver/celloracle_module6_5_input_qc.json"
DEFAULT_CELLORACLE_PYTHON = Path("C:/co/Scripts/python.exe")

REQUIRED_OBS_COLUMNS = [
    "dataset",
    "sample_id",
    "cellrank_fate_prob_cnv_supported_malignant",
    "driver_main_strict__pseudotime_mean",
    "driver_main_strict__pseudotime_phase",
    "driver_primary_cnv_evidence_tier",
]


def read_tf_list(path: Path) -> list[str]:
    values = []
    seen = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        tf = line.strip()
        if not tf or tf in seen:
            continue
        values.append(tf)
        seen.add(tf)
    return values


def check_required_columns(available: Iterable[str], required: Iterable[str]) -> dict:
    available_set = set(available)
    missing = [col for col in required if col not in available_set]
    return {
        "required": list(required),
        "missing": missing,
        "all_present": len(missing) == 0,
    }


def summarize_base_grn_for_tfs(base_grn: pd.DataFrame, tfs: list[str]) -> pd.DataFrame:
    gene_col = "gene_short_name" if "gene_short_name" in base_grn.columns else None
    rows = []
    for tf in tfs:
        if tf in base_grn.columns:
            mask = base_grn[tf].fillna(0).astype(float) != 0
            target_genes = int(base_grn.loc[mask, gene_col].nunique()) if gene_col else int(mask.sum())
            rows.append(
                {
                    "tf": tf,
                    "tf_in_base_grn": True,
                    "base_grn_outgoing_links": int(mask.sum()),
                    "base_grn_target_genes": target_genes,
                }
            )
        else:
            rows.append(
                {
                    "tf": tf,
                    "tf_in_base_grn": False,
                    "base_grn_outgoing_links": 0,
                    "base_grn_target_genes": 0,
                }
            )
    return pd.DataFrame(rows)


def summarize_tf_input(
    tfs: list[str],
    genes: pd.Index,
    selection: pd.DataFrame,
    base_grn_summary: pd.DataFrame,
) -> pd.DataFrame:
    genes_set = set(genes.astype(str))
    out = pd.DataFrame({"tf": tfs})
    out["tf_in_expression"] = out["tf"].isin(genes_set)

    selection_cols = [
        "tf",
        "role",
        "perturbation_mode",
        "tier",
        "total_score",
        "hard_filter_pass",
        "selected_for_main_panel",
    ]
    available_selection_cols = [col for col in selection_cols if col in selection.columns]
    selection_small = selection[available_selection_cols].drop_duplicates("tf")
    selection_small["tf_in_selection_table"] = True

    out = out.merge(selection_small, how="left", on="tf")
    out = out.merge(base_grn_summary, how="left", on="tf")

    out["tf_in_selection_table"] = out["tf_in_selection_table"].where(out["tf_in_selection_table"].notna(), False).astype(bool)
    out["tf_in_base_grn"] = out["tf_in_base_grn"].where(out["tf_in_base_grn"].notna(), False).astype(bool)
    out["base_grn_outgoing_links"] = out["base_grn_outgoing_links"].fillna(0).astype(int)
    out["base_grn_target_genes"] = out["base_grn_target_genes"].fillna(0).astype(int)
    if "hard_filter_pass" in out:
        out["hard_filter_pass"] = out["hard_filter_pass"].where(out["hard_filter_pass"].notna(), False).astype(bool)
    if "selected_for_main_panel" in out:
        out["selected_for_main_panel"] = out["selected_for_main_panel"].where(out["selected_for_main_panel"].notna(), False).astype(bool)
    return out


def classify_celloracle_environment(package_present: bool, import_ok: bool, import_error: str) -> str:
    if package_present and import_ok:
        return "usable"
    if package_present and not import_ok:
        return "package_present_but_import_failed"
    if not package_present and import_error:
        return "not_installed_import_failed"
    return "not_installed"


def probe_celloracle_environment(python_exe: Path) -> dict:
    if not python_exe.exists():
        return {
            "python_executable": str(python_exe),
            "python_exists": False,
            "package_present": False,
            "package_version": None,
            "import_ok": False,
            "import_error": "python_executable_not_found",
            "status": "python_executable_not_found",
        }

    code = r"""
import importlib
import importlib.metadata as md
import json
out = {"package_present": False, "package_version": None, "import_ok": False, "import_error": ""}
try:
    out["package_version"] = md.version("celloracle")
    out["package_present"] = True
except Exception as e:
    out["import_error"] = f"{type(e).__name__}: {e}"
try:
    mod = importlib.import_module("celloracle")
    out["import_ok"] = True
    out["import_version"] = getattr(mod, "__version__", out["package_version"])
except Exception as e:
    out["import_error"] = f"{type(e).__name__}: {e}"
print(json.dumps(out, sort_keys=True))
"""
    proc = subprocess.run(
        [str(python_exe), "-c", code],
        text=True,
        capture_output=True,
        check=False,
        timeout=60,
    )
    try:
        result = json.loads(proc.stdout.strip().splitlines()[-1])
    except Exception:
        result = {
            "package_present": False,
            "package_version": None,
            "import_ok": False,
            "import_error": (proc.stderr or proc.stdout).strip(),
        }
    result["python_executable"] = str(python_exe)
    result["python_exists"] = True
    result["returncode"] = proc.returncode
    result["status"] = classify_celloracle_environment(
        bool(result.get("package_present")),
        bool(result.get("import_ok")),
        str(result.get("import_error", "")),
    )
    return result


def build_qc_report(
    h5ad_path: Path,
    tf_list_path: Path,
    selection_path: Path,
    base_grn_path: Path,
    celloracle_python: Path,
) -> tuple[pd.DataFrame, dict]:
    tfs = read_tf_list(tf_list_path)
    adata = ad.read_h5ad(h5ad_path, backed="r")
    obs_columns = list(adata.obs.columns)
    var_names = adata.var_names.to_series().astype(str)
    shape = tuple(int(x) for x in adata.shape)
    layers = list(adata.layers.keys())
    obsm_keys = list(adata.obsm.keys())
    uns_keys = list(adata.uns.keys())
    n_fate_non_null = int(adata.obs["cellrank_fate_prob_cnv_supported_malignant"].notna().sum()) if "cellrank_fate_prob_cnv_supported_malignant" in adata.obs else 0
    n_dataset = int(adata.obs["dataset"].nunique()) if "dataset" in adata.obs else 0
    n_sample = int(adata.obs["sample_id"].nunique()) if "sample_id" in adata.obs else 0
    adata.file.close()

    selection = pd.read_csv(selection_path, sep="\t")
    base_grn = pd.read_parquet(base_grn_path)
    base_grn_summary = summarize_base_grn_for_tfs(base_grn, tfs)
    tf_qc = summarize_tf_input(tfs, pd.Index(var_names), selection, base_grn_summary)

    obs_check = check_required_columns(obs_columns, REQUIRED_OBS_COLUMNS)
    env = probe_celloracle_environment(celloracle_python)
    input_ready = (
        len(tfs) >= 1
        and shape[0] > 0
        and shape[1] > 0
        and "counts" in layers
        and obs_check["all_present"]
        and bool(tf_qc["tf_in_expression"].all())
        and bool(tf_qc["tf_in_base_grn"].all())
        and bool(tf_qc["tf_in_selection_table"].all())
    )

    report = {
        "module": "6.5",
        "method": "CellOracle environment and input QC",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "project_root": str(PROJECT_ROOT),
        "input_h5ad": str(h5ad_path),
        "tf_list": str(tf_list_path),
        "selection_table": str(selection_path),
        "celloracle_base_grn": str(base_grn_path),
        "h5ad_shape": {"n_cells": shape[0], "n_genes": shape[1]},
        "layers": layers,
        "obsm_keys": obsm_keys,
        "uns_keys": uns_keys,
        "required_obs_check": obs_check,
        "n_fate_non_null_cells": n_fate_non_null,
        "n_datasets": n_dataset,
        "n_samples": n_sample,
        "n_input_tfs": len(tfs),
        "all_input_tfs_in_expression": bool(tf_qc["tf_in_expression"].all()),
        "all_input_tfs_in_base_grn": bool(tf_qc["tf_in_base_grn"].all()),
        "all_input_tfs_in_selection_table": bool(tf_qc["tf_in_selection_table"].all()),
        "all_input_tfs_hard_filter_pass": bool(tf_qc.get("hard_filter_pass", pd.Series([False])).all()),
        "input_data_ready_for_celloracle": bool(input_ready),
        "native_celloracle_environment_ready": bool(env.get("import_ok")),
        "celloracle_environment": env,
        "windows_install_assessment": (
            "not_stable_with_pip_on_current_windows_environment"
            if not env.get("import_ok")
            else "usable"
        ),
        "python_runtime": {
            "executable": sys.executable,
            "version": sys.version,
            "platform": platform.platform(),
        },
        "outputs": {
            "tf_qc": str(DEFAULT_QC_TSV),
            "report": str(DEFAULT_QC_JSON),
        },
    }
    return tf_qc, report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 6.5 CellOracle input QC")
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--tf-list", type=Path, default=DEFAULT_TF_LIST)
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--base-grn", type=Path, default=DEFAULT_BASE_GRN)
    parser.add_argument("--celloracle-python", type=Path, default=DEFAULT_CELLORACLE_PYTHON)
    parser.add_argument("--out-tsv", type=Path, default=DEFAULT_QC_TSV)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_QC_JSON)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    tf_qc, report = build_qc_report(
        h5ad_path=args.h5ad,
        tf_list_path=args.tf_list,
        selection_path=args.selection,
        base_grn_path=args.base_grn,
        celloracle_python=args.celloracle_python,
    )
    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    tf_qc.to_csv(args.out_tsv, sep="\t", index=False)
    args.out_json.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(
        {
            "input_data_ready_for_celloracle": report["input_data_ready_for_celloracle"],
            "native_celloracle_environment_ready": report["native_celloracle_environment_ready"],
            "n_input_tfs": report["n_input_tfs"],
            "tf_qc": str(args.out_tsv),
            "report": str(args.out_json),
        },
        ensure_ascii=False,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
