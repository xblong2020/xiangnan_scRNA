from __future__ import annotations

import argparse
import json
import time
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

PSEUDOTIME_RUNS = {
    "main_strict": {
        "path_name": "trajectory_module5_3_main_strict_pseudotime_merged.tsv.gz",
        "methods": {
            "monocle3": "main_strict__monocle3_norm",
            "slingshot_scanvi": "main_strict__slingshot_scanvi_norm",
            "slingshot_hepatocyte_pca": "main_strict__slingshot_hepatocyte_pca_norm",
        },
        "include_col": "trajectory_include_cnv_strict",
    },
    "sensitivity_include_review": {
        "path_name": "trajectory_module5_3_sensitivity_include_review_pseudotime_merged.tsv.gz",
        "methods": {
            "monocle3": "sensitivity_include_review__monocle3_norm",
            "slingshot_scanvi": "sensitivity_include_review__slingshot_scanvi_norm",
            "slingshot_hepatocyte_pca": "sensitivity_include_review__slingshot_hepatocyte_pca_norm",
        },
        "include_col": "trajectory_include_main",
    },
}

MODULE_COLUMNS = [
    "Mature_Hepatocyte",
    "Stressed_Injured",
    "Regenerative_Progenitor",
    "HCC_Malignant_Associated",
    "Proliferation",
    "Cholangiocyte",
    "Immune",
    "Endothelial",
    "Stromal_HSC_Pericyte",
]

CORE_CELL_COLUMNS = [
    "cell_id",
    "dataset",
    "sample_id",
    "study_sample",
    "cnv_sample",
    "source_h5ad",
    "_scvi_batch",
    "leiden_scvi",
    "leiden_hep",
    "leiden_trajectory",
    "sample_source_class",
    "sample_disease_stage",
    "cell_disease_stage",
    "trajectory_role",
    "trajectory_root_end_role",
    "trajectory_include_main",
    "trajectory_include_cnv_strict",
    "scanvi_unified_final_label",
    "scanvi_unified_final_strict_label",
    "hepatocyte_state_label",
    "malignant_hcc_call",
    "malignant_hcc_cnv_method",
    "copykat_pred",
    "copykat_status",
    "cnv_proxy_status",
    "cnv_proxy_burden",
    "cnv_proxy_z",
    "cnv_proxy_high_bin_fraction",
    "cnv_proxy_max_abs_bin_log2",
    "hcc_malignant_associated_score_z",
    "proliferation_score_z",
    "regenerative_progenitor_score_z",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 6.1: prepare driver analysis AnnData object and manifests.")
    parser.add_argument(
        "--input-h5ad",
        type=Path,
        default=ROOT / "data/processed/trajectory/trajectory_hepatocyte_cnv_scanvi.stage_root_end.module5_2.h5ad",
    )
    parser.add_argument("--metadata-trajectory-dir", type=Path, default=ROOT / "metadata/trajectory")
    parser.add_argument("--processed-dir", type=Path, default=ROOT / "data/processed/driver")
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/driver")
    parser.add_argument("--output-name", default="driver_hepatocyte_trajectory.module6_1.h5ad")
    parser.add_argument("--n-bins", type=int, default=10)
    parser.add_argument("--primary-run-id", default="main_strict")
    parser.add_argument("--compression", default="gzip")
    return parser.parse_args()


def write_dataframe(path: Path, df: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = "gzip" if path.suffix == ".gz" else None
    df.to_csv(path, sep="\t", index=False, compression=compression)


def bool_series(values: pd.Series) -> pd.Series:
    labels = values.astype("object").map(lambda value: "false" if pd.isna(value) else str(value).strip().lower())
    return labels.isin({"true", "1", "yes", "y"})


def assign_bins(values: pd.Series, n_bins: int) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    out = pd.Series(pd.NA, index=values.index, dtype="Int64")
    finite = numeric.replace([np.inf, -np.inf], np.nan).dropna()
    if finite.empty:
        return out
    n_effective = min(n_bins, int(finite.nunique()))
    if n_effective <= 1:
        out.loc[finite.index] = 0
        return out
    ranked = finite.rank(method="first")
    out.loc[finite.index] = pd.qcut(ranked, q=n_effective, labels=list(range(n_effective)), duplicates="drop").astype("Int64")
    return out


def driver_phase(bin_values: pd.Series, n_bins: int) -> pd.Series:
    out = pd.Series(pd.NA, index=bin_values.index, dtype="object")
    numeric = pd.to_numeric(bin_values, errors="coerce")
    early_cut = max(0, int(np.floor(n_bins * 0.3)) - 1)
    late_cut = int(np.ceil(n_bins * 0.7))
    out.loc[numeric.le(early_cut)] = "early"
    out.loc[numeric.gt(early_cut) & numeric.lt(late_cut)] = "middle"
    out.loc[numeric.ge(late_cut)] = "late"
    return out


def minmax_rank(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    finite = numeric.replace([np.inf, -np.inf], np.nan).dropna()
    out = pd.Series(np.nan, index=values.index, dtype=float)
    if finite.empty:
        return out
    if finite.nunique() == 1:
        out.loc[finite.index] = 0.0
        return out
    ranks = finite.rank(method="average")
    out.loc[finite.index] = (ranks - 1.0) / float(len(finite) - 1)
    return out


def add_pseudotime_columns(cells: pd.DataFrame, trajectory_dir: Path, n_bins: int) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    cells = cells.copy()
    records: list[dict[str, object]] = []
    for run_id, config in PSEUDOTIME_RUNS.items():
        path = trajectory_dir / str(config["path_name"])
        pt = pd.read_csv(path, sep="\t")
        pt["cell_id"] = pt["cell_id"].astype(str)
        pt = pt.set_index("cell_id")
        aligned = pt.reindex(cells["cell_id"].astype(str))

        method_cols = []
        for method, source_col in config["methods"].items():
            out_col = f"driver_{run_id}__{method}_pseudotime"
            cells[out_col] = pd.to_numeric(aligned[source_col], errors="coerce").to_numpy()
            method_cols.append(out_col)
            records.append(
                {
                    "run_id": run_id,
                    "method": method,
                    "source_column": source_col,
                    "n_cells_with_pseudotime": int(cells[out_col].notna().sum()),
                }
            )

        values = cells[method_cols].apply(pd.to_numeric, errors="coerce")
        cells[f"driver_{run_id}__pseudotime_mean"] = values.mean(axis=1, skipna=True)
        cells[f"driver_{run_id}__pseudotime_median"] = values.median(axis=1, skipna=True)
        cells[f"driver_{run_id}__pseudotime_std"] = values.std(axis=1, skipna=True)
        cells[f"driver_{run_id}__pseudotime_n_methods"] = values.notna().sum(axis=1).astype(int)
        cells[f"driver_{run_id}__pseudotime_rank"] = minmax_rank(cells[f"driver_{run_id}__pseudotime_median"])
        bin_col = f"driver_{run_id}__pseudotime_bin{n_bins}"
        cells[bin_col] = assign_bins(cells[f"driver_{run_id}__pseudotime_median"], n_bins)
        cells[f"driver_{run_id}__pseudotime_phase"] = driver_phase(cells[bin_col], n_bins)
        include_col = str(config["include_col"])
        include = bool_series(cells[include_col]) if include_col in cells.columns else pd.Series(False, index=cells.index)
        cells[f"driver_{run_id}__eligible"] = include & cells[f"driver_{run_id}__pseudotime_median"].notna()

    primary = "driver_main_strict__eligible"
    cells["driver_primary_eligible"] = cells[primary] if primary in cells.columns else False
    return cells, records


def add_module_scores(cells: pd.DataFrame, module_scores_path: Path) -> pd.DataFrame:
    if not module_scores_path.exists():
        return cells
    scores = pd.read_csv(module_scores_path, sep="\t")
    scores["cell_id"] = scores["cell_id"].astype(str)
    cells = cells.copy()
    for run_id, sub in scores.groupby("run_id", observed=True, sort=True):
        available = [column for column in MODULE_COLUMNS if column in sub.columns]
        if not available:
            continue
        aligned = sub.set_index("cell_id")[available].reindex(cells["cell_id"].astype(str))
        for column in available:
            cells[f"driver_{run_id}__module_{column}"] = pd.to_numeric(aligned[column], errors="coerce").to_numpy()
    return cells


def build_run_method_table(cells: pd.DataFrame, n_bins: int) -> pd.DataFrame:
    rows = []
    common = [
        "cell_id",
        "dataset",
        "sample_id",
        "study_sample",
        "cnv_sample",
        "sample_source_class",
        "sample_disease_stage",
        "cell_disease_stage",
        "trajectory_role",
        "trajectory_root_end_role",
        "malignant_hcc_call",
        "copykat_status",
        "cnv_proxy_status",
    ]
    for run_id, config in PSEUDOTIME_RUNS.items():
        run_base_cols = [
            f"driver_{run_id}__pseudotime_median",
            f"driver_{run_id}__pseudotime_rank",
            f"driver_{run_id}__pseudotime_bin{n_bins}",
            f"driver_{run_id}__pseudotime_phase",
            f"driver_{run_id}__eligible",
        ]
        module_cols = [f"driver_{run_id}__module_{column}" for column in MODULE_COLUMNS if f"driver_{run_id}__module_{column}" in cells.columns]
        for method in config["methods"]:
            method_col = f"driver_{run_id}__{method}_pseudotime"
            keep = [column for column in common + run_base_cols + module_cols if column in cells.columns]
            sub = cells.loc[cells[method_col].notna(), keep].copy()
            sub.insert(1, "run_id", run_id)
            sub.insert(2, "method", method)
            sub.insert(3, "pseudotime_norm", cells.loc[cells[method_col].notna(), method_col].to_numpy())
            sub = sub.rename(
                columns={
                    f"driver_{run_id}__pseudotime_median": "consensus_pseudotime_median",
                    f"driver_{run_id}__pseudotime_rank": "consensus_pseudotime_rank",
                    f"driver_{run_id}__pseudotime_bin{n_bins}": f"consensus_pseudotime_bin{n_bins}",
                    f"driver_{run_id}__pseudotime_phase": "consensus_pseudotime_phase",
                    f"driver_{run_id}__eligible": "driver_eligible",
                }
            )
            for module in MODULE_COLUMNS:
                prefixed = f"driver_{run_id}__module_{module}"
                if prefixed in sub.columns:
                    sub = sub.rename(columns={prefixed: module})
            rows.append(sub)
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, axis=0, ignore_index=True)


def add_primary_overlay_evidence(cells: pd.DataFrame, overlay_path: Path, primary_run_id: str) -> pd.DataFrame:
    if not overlay_path.exists():
        return cells
    desired = [
        "run_id",
        "method",
        "cell_id",
        "cnv_evidence_tier",
        "module3_cnv_supported",
        "malignant_like_review",
        "copykat_aneuploid",
        "cnv_proxy_aneuploid",
    ]
    overlay = pd.read_csv(overlay_path, sep="\t", usecols=lambda col: col in desired)
    overlay = overlay.loc[overlay["run_id"].astype(str).eq(primary_run_id)].copy()
    if overlay.empty:
        return cells
    summary = (
        overlay.groupby("cell_id", observed=True)
        .agg(
            driver_primary_cnv_evidence_tier=("cnv_evidence_tier", first_nonempty),
            driver_primary_module3_cnv_supported=("module3_cnv_supported", any_bool),
            driver_primary_malignant_like_review=("malignant_like_review", any_bool),
            driver_primary_copykat_aneuploid=("copykat_aneuploid", any_bool),
            driver_primary_cnv_proxy_aneuploid=("cnv_proxy_aneuploid", any_bool),
            driver_primary_n_methods_with_overlay=("method", "nunique"),
        )
        .reset_index()
    )
    aligned = summary.set_index("cell_id").reindex(cells["cell_id"].astype(str))
    cells = cells.copy()
    for column in summary.columns:
        if column == "cell_id":
            continue
        cells[column] = aligned[column].to_numpy()
    bool_cols = [
        "driver_primary_module3_cnv_supported",
        "driver_primary_malignant_like_review",
        "driver_primary_copykat_aneuploid",
        "driver_primary_cnv_proxy_aneuploid",
    ]
    for column in bool_cols:
        cells[column] = bool_series(cells[column])
    cells["driver_primary_cnv_evidence_tier"] = cells["driver_primary_cnv_evidence_tier"].fillna("not_in_primary_driver_set")
    cells["driver_primary_n_methods_with_overlay"] = cells["driver_primary_n_methods_with_overlay"].fillna(0).astype(int)
    return cells


def first_nonempty(values: pd.Series) -> str:
    nonempty = values.dropna().astype(str)
    nonempty = nonempty.loc[~nonempty.str.lower().isin({"", "nan", "none", "<na>"})]
    if nonempty.empty:
        return ""
    return str(nonempty.iloc[0])


def any_bool(values: pd.Series) -> bool:
    return bool(bool_series(values).any())


def build_candidate_sets(cells: pd.DataFrame, n_bins: int) -> pd.DataFrame:
    rows = []
    definitions = []
    for run_id in PSEUDOTIME_RUNS:
        eligible = bool_series(cells[f"driver_{run_id}__eligible"])
        bin_col = f"driver_{run_id}__pseudotime_bin{n_bins}"
        phase_col = f"driver_{run_id}__pseudotime_phase"
        definitions.extend(
            [
                (f"{run_id}__all_eligible", run_id, eligible),
                (f"{run_id}__early_phase", run_id, eligible & cells[phase_col].astype(str).eq("early")),
                (f"{run_id}__middle_phase", run_id, eligible & cells[phase_col].astype(str).eq("middle")),
                (f"{run_id}__late_phase", run_id, eligible & cells[phase_col].astype(str).eq("late")),
                (f"{run_id}__early_vs_late_contrast", run_id, eligible & cells[phase_col].astype(str).isin({"early", "late"})),
                (f"{run_id}__late_bin", run_id, eligible & pd.to_numeric(cells[bin_col], errors="coerce").eq(n_bins - 1)),
            ]
        )
    if "driver_primary_module3_cnv_supported" in cells.columns:
        primary = bool_series(cells["driver_primary_eligible"])
        cnv = bool_series(cells["driver_primary_module3_cnv_supported"])
        late = cells[f"driver_main_strict__pseudotime_phase"].astype(str).eq("late")
        definitions.extend(
            [
                ("main_strict__cnv_supported", "main_strict", primary & cnv),
                ("main_strict__late_cnv_supported", "main_strict", primary & late & cnv),
                ("main_strict__late_without_module3_cnv", "main_strict", primary & late & ~cnv),
            ]
        )

    for name, run_id, mask in definitions:
        sub = cells.loc[mask].copy()
        run_phase_col = f"driver_{run_id}__pseudotime_phase"
        rows.append(
            {
                "candidate_set": name,
                "run_id": run_id,
                "n_cells": int(sub.shape[0]),
                "n_datasets": int(sub["dataset"].astype(str).nunique()) if "dataset" in sub.columns else 0,
                "n_samples": int(sub["sample_id"].astype(str).nunique()) if "sample_id" in sub.columns else 0,
                "n_cnv_samples": int(sub["cnv_sample"].astype(str).nunique()) if "cnv_sample" in sub.columns else 0,
                "module3_cnv_supported_fraction": float(bool_series(sub.get("driver_primary_module3_cnv_supported", pd.Series(False))).mean())
                if sub.shape[0] > 0
                else np.nan,
                "late_phase_fraction": float(sub.get(run_phase_col, pd.Series("", index=sub.index)).astype(str).eq("late").mean())
                if sub.shape[0] > 0
                else np.nan,
            }
        )
    return pd.DataFrame(rows)


def build_feature_manifest(adata: ad.AnnData, cells: pd.DataFrame, n_bins: int) -> pd.DataFrame:
    rows = []
    rows.append({"feature_group": "expression_gene", "feature": "var_names", "n_features": int(adata.n_vars), "source": "module5_2_h5ad"})
    for layer in adata.layers.keys():
        rows.append({"feature_group": "expression_layer", "feature": str(layer), "n_features": int(adata.n_vars), "source": "module5_2_h5ad"})
    for key in adata.obsm.keys():
        rows.append({"feature_group": "embedding", "feature": str(key), "n_features": int(adata.obsm[key].shape[1]), "source": "module5_2_h5ad"})
    for run_id in PSEUDOTIME_RUNS:
        rows.append(
            {
                "feature_group": "pseudotime",
                "feature": f"driver_{run_id}__pseudotime_median",
                "n_features": 1,
                "source": "module5_3_consensus",
            }
        )
        rows.append(
            {
                "feature_group": "pseudotime_bin",
                "feature": f"driver_{run_id}__pseudotime_bin{n_bins}",
                "n_features": 1,
                "source": "module5_3_consensus",
            }
        )
        for module in MODULE_COLUMNS:
            column = f"driver_{run_id}__module_{module}"
            if column in cells.columns:
                rows.append({"feature_group": "module_score", "feature": column, "n_features": 1, "source": "module5_4"})
    for column in [
        "dataset",
        "sample_id",
        "study_sample",
        "cnv_sample",
        "_scvi_batch",
        "sample_source_class",
        "sample_disease_stage",
        "cell_disease_stage",
        "trajectory_role",
        "trajectory_root_end_role",
    ]:
        if column in cells.columns:
            rows.append({"feature_group": "covariate", "feature": column, "n_features": 1, "source": "module5_2_obs"})
    return pd.DataFrame(rows)


def attach_obs(adata: ad.AnnData, cells: pd.DataFrame) -> None:
    aligned = cells.set_index("cell_id").loc[adata.obs_names.astype(str)]
    for column in aligned.columns:
        if column == "cell_id":
            continue
        series = aligned[column]
        if pd.api.types.is_bool_dtype(series):
            adata.obs[column] = series.fillna(False).astype(bool).to_numpy()
        elif pd.api.types.is_integer_dtype(series):
            adata.obs[column] = series.to_numpy()
        elif pd.api.types.is_float_dtype(series):
            adata.obs[column] = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
        else:
            text = series.astype("object").where(series.notna(), "Unknown").astype(str)
            unique = text.loc[~text.str.lower().isin({"unknown", "nan", "none", "<na>"})].nunique()
            if unique <= max(1000, int(len(series) * 0.5)):
                adata.obs[column] = pd.Categorical(text.to_numpy())
            else:
                adata.obs[column] = text.where(~text.eq("Unknown"), "").to_numpy()


def main() -> None:
    start = time.time()
    args = parse_args()
    args.processed_dir.mkdir(parents=True, exist_ok=True)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.input_h5ad)
    cells = adata.obs.copy()
    cells.insert(0, "cell_id", adata.obs_names.astype(str))
    cells, pseudotime_records = add_pseudotime_columns(cells, args.metadata_trajectory_dir, args.n_bins)
    cells = add_module_scores(cells, args.metadata_trajectory_dir / "trajectory_module5_4_module_scores_by_cell.tsv.gz")
    cells = add_primary_overlay_evidence(
        cells,
        args.metadata_trajectory_dir / "trajectory_module5_5_cnv_malignant_overlay_by_cell.tsv.gz",
        args.primary_run_id,
    )

    driver_mask = cells[[f"driver_{run_id}__eligible" for run_id in PSEUDOTIME_RUNS]].any(axis=1)
    driver_cells = cells.loc[driver_mask].copy()
    driver_adata = adata[driver_cells.index.to_numpy(), :].copy()
    attach_obs(driver_adata, driver_cells)

    driver_adata.uns["module6_1_driver_object"] = {
        "module": "6.1",
        "method": "driver analysis object prepared from module 5 trajectory consensus and evidence overlays",
        "input_h5ad": str(args.input_h5ad),
        "primary_run_id": args.primary_run_id,
        "pseudotime_runs": list(PSEUDOTIME_RUNS.keys()),
        "n_bins": int(args.n_bins),
        "driver_object_scope": "union of cells with any module 5.3 normalized pseudotime",
        "sample_composition_caveat": "module 5.7 supported the trajectory evidence with sample composition caveat",
    }

    output_h5ad = args.processed_dir / args.output_name
    driver_adata.write_h5ad(output_h5ad, compression=args.compression)

    core_cols = [column for column in CORE_CELL_COLUMNS if column in driver_cells.columns]
    driver_cols = [column for column in driver_cells.columns if column.startswith("driver_")]
    cell_table = driver_cells[core_cols + driver_cols].copy()
    cells_path = args.metadata_dir / "driver_module6_1_cells.tsv.gz"
    write_dataframe(cells_path, cell_table)

    long_table = build_run_method_table(driver_cells, args.n_bins)
    long_path = args.metadata_dir / "driver_module6_1_run_method_long.tsv.gz"
    write_dataframe(long_path, long_table)

    candidate_sets = build_candidate_sets(driver_cells, args.n_bins)
    candidate_sets_path = args.metadata_dir / "driver_module6_1_candidate_sets.tsv"
    write_dataframe(candidate_sets_path, candidate_sets)

    feature_manifest = build_feature_manifest(driver_adata, driver_cells, args.n_bins)
    feature_manifest_path = args.metadata_dir / "driver_module6_1_feature_manifest.tsv"
    write_dataframe(feature_manifest_path, feature_manifest)

    pseudotime_manifest = pd.DataFrame(pseudotime_records)
    pseudotime_manifest_path = args.metadata_dir / "driver_module6_1_pseudotime_manifest.tsv"
    write_dataframe(pseudotime_manifest_path, pseudotime_manifest)

    report = {
        "module": "6.1",
        "method": "driver analysis AnnData object preparation",
        "input_h5ad": str(args.input_h5ad),
        "output_h5ad": str(output_h5ad),
        "n_cells_input": int(adata.n_obs),
        "n_cells_driver_union": int(driver_adata.n_obs),
        "n_genes": int(driver_adata.n_vars),
        "layers": list(driver_adata.layers.keys()),
        "obsm": list(driver_adata.obsm.keys()),
        "primary_run_id": args.primary_run_id,
        "eligible_cells_by_run": {
            run_id: int(driver_cells[f"driver_{run_id}__eligible"].sum()) for run_id in PSEUDOTIME_RUNS
        },
        "candidate_sets": candidate_sets.to_dict(orient="records"),
        "pseudotime_manifest": pseudotime_manifest.to_dict(orient="records"),
        "outputs": {
            "h5ad": str(output_h5ad),
            "cells": str(cells_path),
            "run_method_long": str(long_path),
            "candidate_sets": str(candidate_sets_path),
            "feature_manifest": str(feature_manifest_path),
            "pseudotime_manifest": str(pseudotime_manifest_path),
        },
        "package_versions": {
            "anndata": version("anndata"),
            "pandas": version("pandas"),
            "numpy": version("numpy"),
        },
        "elapsed_seconds": round(time.time() - start, 3),
    }
    report_path = args.metadata_dir / "driver_module6_1_report.json"
    report["outputs"]["report"] = str(report_path)
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
