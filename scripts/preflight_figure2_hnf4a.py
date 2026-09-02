#!/usr/bin/env python3
"""Audit all inputs required for HNF4A Figure 2B-F without modifying sources."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


PROJECT_ROOT = Path(__file__).resolve().parents[1]
TARGET_TF = "HNF4A"


def as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.astype(str).str.lower().isin({"true", "1", "yes", "y"})


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer")


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def expression_vector(adata: ad.AnnData, gene: str) -> np.ndarray:
    idx = int(adata.var_names.get_loc(gene))
    matrix = adata.X[:, idx]
    if sparse.issparse(matrix):
        return np.asarray(matrix.toarray()).ravel()
    return np.asarray(matrix).ravel()


def regulon_vector(adata: ad.AnnData, tf: str) -> tuple[np.ndarray | None, str | None]:
    key = "module6_3c_cistarget_regulon_auc"
    names_key = "module6_3c_cistarget_regulon_auc_names"
    if key not in adata.obsm or names_key not in adata.uns:
        return None, None
    names = [str(x) for x in adata.uns[names_key]]
    candidates = [tf, f"{tf}(+)"]
    for candidate in candidates:
        if candidate in names:
            return np.asarray(adata.obsm[key])[:, names.index(candidate)], candidate
    return None, None


def run(root: Path, target_tf: str, out_json: Path, out_tsv: Path) -> dict:
    h5ad_path = root / "data/processed/driver/celloracle_module6_6/celloracle_module6_6_input.h5ad"
    fitted_oracle = root / "data/processed/driver/celloracle_module6_7/celloracle_module6_7_fitted.celloracle.oracle"
    network_summary_path = root / "metadata/driver/celloracle_module6_7_tf_network_summary.tsv"
    perturb_path = root / "metadata/driver/celloracle_module6_8_cell_shift_summary.tsv.gz"
    grid_path = root / "metadata/driver/celloracle_module6_8_grid_arrows.tsv.gz"
    perturb_report_path = root / "metadata/driver/celloracle_module6_8_perturbation_report.json"
    baseline_cells_path = root / "metadata/driver/figure2b_sox4/figure2b_sox4_plot_cells.tsv.gz"
    baseline_grid_path = root / "metadata/driver/figure2b_sox4/figure2b_sox4_baseline_grid_umap.tsv.gz"
    sct_paths = {
        "main_strict": root / "metadata/driver/sctenifoldknk_module7_2_main_strict_perturbation_genes.tsv",
        "malignant_like": root / "metadata/driver/sctenifoldknk_module7_2_malignant_like_perturbation_genes.tsv",
        "normal_reference": root / "metadata/driver/figure2e_hnf4a_sctenifoldknk/figure2e_hnf4a_normal_reference_perturbation_genes.tsv",
    }

    required_paths = [h5ad_path, fitted_oracle, network_summary_path, perturb_path, grid_path, perturb_report_path,
                      baseline_cells_path, baseline_grid_path]
    missing_paths = [str(p) for p in required_paths if not p.exists()]
    if missing_paths:
        raise FileNotFoundError("Missing required preflight inputs: " + "; ".join(missing_paths))

    hash_path = root / "metadata/driver/figure2_hnf4a_sox4_reference_hashes.tsv"
    if not hash_path.exists():
        sox4_files = sorted(
            p for parent in [root / "scripts", root / "metadata/driver", root / "figures/driver"]
            for p in parent.rglob("*sox4*") if p.is_file()
        )
        pd.DataFrame([
            {"path": str(p.relative_to(root)).replace("\\", "/"), "size_bytes": p.stat().st_size,
             "sha256": file_sha256(p)}
            for p in sox4_files
        ]).to_csv(hash_path, sep="\t", index=False)

    adata = ad.read_h5ad(h5ad_path)
    required_obs = [
        "driver_main_strict__eligible",
        "driver_main_strict__pseudotime_rank",
        "celloracle_state",
    ]
    obs_presence = {col: col in adata.obs.columns for col in required_obs}
    obsm_presence = {"X_celloracle_umap": "X_celloracle_umap" in adata.obsm}
    if target_tf not in adata.var_names:
        raise ValueError(f"{target_tf} is absent from the CellOracle expression matrix")

    expr = expression_vector(adata, target_tf)
    regulon, regulon_name = regulon_vector(adata, target_tf)
    states = adata.obs["celloracle_state"].astype(str)
    state_rows = []
    for state in sorted(states.unique()):
        mask = states.eq(state).to_numpy()
        values = expr[mask]
        reg_values = regulon[mask] if regulon is not None else None
        state_rows.append({
            "record_type": "state_expression",
            "state_or_subset": state,
            "n_cells": int(mask.sum()),
            "mean_expression": float(np.mean(values)),
            "detection_rate": float(np.mean(values > 0)),
            "mean_regulon_activity": float(np.mean(reg_values)) if reg_values is not None else np.nan,
            "available": True,
            "detail": regulon_name or "regulon unavailable",
        })
    state_stats = pd.DataFrame(state_rows)

    network = read_table(network_summary_path)
    in_network = bool((network["tf"].astype(str) == target_tf).any())
    perturb = read_table(perturb_path)
    target_perturb = perturb.loc[perturb["tf"].astype(str).eq(target_tf)].copy()
    grid = read_table(grid_path)
    target_grid = grid.loc[grid["tf"].astype(str).eq(target_tf)].copy()
    delta_columns = ["delta_embedding_1", "delta_embedding_2"]
    delta_available = bool(
        len(target_perturb) == adata.n_obs
        and all(c in target_perturb.columns for c in delta_columns)
        and np.isfinite(target_perturb[delta_columns].to_numpy(dtype=float)).all()
    )

    strict_mask = as_bool(adata.obs["driver_main_strict__eligible"]) & pd.to_numeric(
        adata.obs["driver_main_strict__pseudotime_rank"], errors="coerce"
    ).notna()
    baseline = read_table(baseline_cells_path)
    strict_ids = set(adata.obs_names[strict_mask].astype(str))
    baseline_ids = set(baseline["cell_id"].astype(str))
    strict_complete = bool(len(strict_ids) == 5000 and strict_ids == baseline_ids)

    with perturb_report_path.open(encoding="utf-8-sig") as handle:
        perturb_report = json.load(handle)
    parameters = perturb_report.get("parameters", {})
    perturbed_tfs = perturb_report.get("result", {}).get("perturbed_tfs", [])

    sct_summary = []
    for subset, path in sct_paths.items():
        available = path.exists()
        n_target_rows = 0
        n_significant = 0
        if available:
            dat = read_table(path)
            tf_dat = dat.loc[dat["tf"].astype(str).eq(target_tf)] if "tf" in dat.columns else dat
            n_target_rows = int(len(tf_dat))
            if "p.adj" in tf_dat.columns:
                n_significant = int((pd.to_numeric(tf_dat["p.adj"], errors="coerce") < 0.05).sum())
        if subset == "normal_reference":
            n_cells = int(states.eq("normal_reference").sum())
        elif subset == "malignant_like":
            n_cells = int(states.eq("malignant_or_malignant_like").sum())
        else:
            n_cells = int(strict_mask.sum())
        sct_summary.append({
            "record_type": "sctenifold_subset",
            "state_or_subset": subset,
            "n_cells": n_cells,
            "mean_expression": None,
            "detection_rate": None,
            "mean_regulon_activity": None,
            "available": available,
            "detail": f"target_rows={n_target_rows};fdr_lt_0.05={n_significant};path={path}",
        })

    normal_row = state_stats.loc[state_stats["state_or_subset"].eq("normal_reference")].iloc[0]
    selected_subset = "normal_reference"
    selection_reason = (
        "normal_reference selected as identity-high non-malignant state: "
        f"n={int(normal_row.n_cells)}, mean_expression={normal_row.mean_expression:.6g}, "
        f"detection_rate={normal_row.detection_rate:.6g}, "
        f"mean_regulon_activity={normal_row.mean_regulon_activity:.6g}"
    )

    audit_rows = state_rows + sct_summary + [
        {
            "record_type": "preflight_check",
            "state_or_subset": "HNF4A_in_fitted_network",
            "n_cells": np.nan,
            "mean_expression": np.nan,
            "detection_rate": np.nan,
            "mean_regulon_activity": np.nan,
            "available": in_network,
            "detail": str(network_summary_path),
        },
        {
            "record_type": "preflight_check",
            "state_or_subset": "HNF4A_celloracle_delta_embedding",
            "n_cells": len(target_perturb),
            "mean_expression": np.nan,
            "detection_rate": np.nan,
            "mean_regulon_activity": np.nan,
            "available": delta_available,
            "detail": str(perturb_path),
        },
        {
            "record_type": "preflight_check",
            "state_or_subset": "HNF4A_grid_arrows",
            "n_cells": len(target_grid),
            "mean_expression": np.nan,
            "detection_rate": np.nan,
            "mean_regulon_activity": np.nan,
            "available": len(target_grid) > 0,
            "detail": str(grid_path),
        },
        {
            "record_type": "preflight_check",
            "state_or_subset": "strict_5000_baseline_complete",
            "n_cells": len(baseline_ids),
            "mean_expression": np.nan,
            "detection_rate": np.nan,
            "mean_regulon_activity": np.nan,
            "available": strict_complete,
            "detail": str(baseline_cells_path),
        },
    ]
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(audit_rows).to_csv(out_tsv, sep="\t", index=False, na_rep="NA")

    report = {
        "module": "Figure 2 HNF4A preflight",
        "target_tf": target_tf,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "h5ad": str(h5ad_path.resolve()),
            "fitted_celloracle_object": str(fitted_oracle.resolve()),
            "celloracle_network_summary": str(network_summary_path.resolve()),
            "celloracle_cell_shift": str(perturb_path.resolve()),
            "celloracle_grid_arrows": str(grid_path.resolve()),
            "sox4_baseline_cells": str(baseline_cells_path.resolve()),
            "sox4_baseline_grid_umap": str(baseline_grid_path.resolve()),
        },
        "checks": {
            "target_in_fitted_network": in_network,
            "target_in_celloracle_simulation": target_tf in perturbed_tfs,
            "delta_embedding_available": delta_available,
            "virtual_knockout_summary_available": len(target_perturb) > 0,
            "grid_arrows_available": len(target_grid) > 0,
            "strict_5000_complete": strict_complete,
            "required_obs_columns": obs_presence,
            "required_obsm": obsm_presence,
        },
        "celloracle_parameters": {
            "n_propagation": parameters.get("n_propagation"),
            "seed": parameters.get("seed"),
            "grid_steps": parameters.get("grid_steps"),
            "grid_neighbors": parameters.get("grid_neighbors"),
            "source_report": str(perturb_report_path.resolve()),
        },
        "state_statistics": state_stats.drop(columns="record_type").to_dict(orient="records"),
        "sctenifold_subsets": sct_summary,
        "selected_identity_high_subset": selected_subset,
        "selection_reason": selection_reason,
        "normal_reference_network_required": not sct_paths["normal_reference"].exists(),
        "scientific_language": [
            "predicted perturbation",
            "virtual knockout",
            "computationally inferred state shift",
            "network perturbation evidence",
        ],
        "outputs": {"json": str(out_json.resolve()), "tsv": str(out_tsv.resolve())},
        "sox4_reference_hashes": str(hash_path.resolve()),
    }
    out_json.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--target-tf", default=TARGET_TF)
    parser.add_argument(
        "--out-json",
        type=Path,
        default=PROJECT_ROOT / "metadata/driver/figure2_hnf4a_preflight_report.json",
    )
    parser.add_argument(
        "--out-tsv",
        type=Path,
        default=PROJECT_ROOT / "metadata/driver/figure2_hnf4a_preflight_report.tsv",
    )
    args = parser.parse_args()
    print(json.dumps(run(args.project_root, args.target_tf, args.out_json, args.out_tsv), ensure_ascii=False))


if __name__ == "__main__":
    main()
