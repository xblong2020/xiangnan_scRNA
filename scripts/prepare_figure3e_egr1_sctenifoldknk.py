#!/usr/bin/env python3
"""Audit and export stress-transition subsets for EGR1 scTenifoldKnk."""

from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io, sparse

try:
    from figure3_egr1_common import (
        PROJECT_ROOT,
        PSEUDOTIME_COLUMN,
        SEED,
        TARGET_TF,
        as_bool,
        choose_stress_transition_subset,
        json_safe,
        write_json,
    )
except ModuleNotFoundError:
    from scripts.figure3_egr1_common import (
        PROJECT_ROOT,
        PSEUDOTIME_COLUMN,
        SEED,
        TARGET_TF,
        as_bool,
        choose_stress_transition_subset,
        json_safe,
        write_json,
    )


DEFAULT_H5AD = PROJECT_ROOT / "data/processed/driver/celloracle_module6_6/celloracle_module6_6_input.h5ad"
DEFAULT_DATA_DIR = PROJECT_ROOT / "data/processed/driver/figure3e_egr1_sctenifoldknk"
DEFAULT_METADATA_DIR = PROJECT_ROOT / "metadata/driver/figure3e_egr1"
SUBSET_ORDER = [
    "stressed_injured",
    "stressed_regenerative",
    "intermediate_pseudotime",
    "malignant_like",
]


def subset_masks(adata: ad.AnnData) -> dict[str, pd.Series]:
    state = adata.obs["celloracle_state"].astype(str)
    pseudotime = pd.to_numeric(adata.obs[PSEUDOTIME_COLUMN], errors="coerce")
    strict = as_bool(adata.obs["driver_main_strict__eligible"]) & pseudotime.notna()
    return {
        "stressed_injured": state.eq("stressed_injured"),
        "stressed_regenerative": state.isin(["stressed_injured", "regenerative_progenitor"]),
        "intermediate_pseudotime": strict & pseudotime.ge(0.33) & pseudotime.lt(0.67),
        "malignant_like": state.eq("malignant_or_malignant_like"),
    }


def _matrix_column(matrix, index: int) -> np.ndarray:
    column = matrix[:, index]
    return column.toarray().ravel() if sparse.issparse(column) else np.asarray(column).ravel()


def audit_subsets(adata: ad.AnnData, masks: dict[str, pd.Series]) -> pd.DataFrame:
    if TARGET_TF not in adata.var_names:
        raise ValueError(f"{TARGET_TF} is absent from the 3,000-gene network")
    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
    egr1 = _matrix_column(counts, int(adata.var_names.get_loc(TARGET_TF)))
    auc_columns = [column for column in adata.obs.columns if "egr1" in column.lower() and "auc" in column.lower()]
    rows = []
    for priority, subset in enumerate(SUBSET_ORDER, start=1):
        mask = masks[subset].to_numpy(dtype=bool)
        obs = adata.obs.loc[mask]
        datasets = obs["dataset"].astype(str)
        known_datasets = datasets.loc[~datasets.str.lower().isin({"unknown", "nan", "none", ""})]
        samples = obs["sample_id"].astype(str)
        known_samples = samples.loc[~samples.str.lower().isin({"unknown", "nan", "none", ""})]
        dataset_counts = datasets.value_counts()
        sample_counts = samples.value_counts()
        values = egr1[mask]
        rows.append(
            {
                "selection_priority": priority,
                "subset": subset,
                "n_cells": int(mask.sum()),
                "n_genes": int(adata.n_vars),
                "egr1_detection_rate": float(np.mean(values > 0)) if len(values) else 0.0,
                "egr1_mean_expression": float(np.mean(values)) if len(values) else 0.0,
                "egr1_regulon_auc_available": bool(auc_columns),
                "egr1_regulon_auc_mean": (
                    float(pd.to_numeric(obs[auc_columns[0]], errors="coerce").mean())
                    if auc_columns and len(obs)
                    else np.nan
                ),
                "n_datasets": int(known_datasets.nunique()),
                "n_datasets_including_unknown": int(datasets.nunique()),
                "n_samples_or_patients": int(known_samples.nunique()),
                "n_samples_including_unknown": int(samples.nunique()),
                "max_dataset_fraction": float(dataset_counts.iloc[0] / len(obs)) if len(obs) else 1.0,
                "max_sample_fraction": float(sample_counts.iloc[0] / len(obs)) if len(obs) else 1.0,
                "largest_dataset": str(dataset_counts.index[0]) if len(dataset_counts) else "",
                "largest_sample": str(sample_counts.index[0]) if len(sample_counts) else "",
                "dataset_composition": ";".join(f"{key}:{value}" for key, value in dataset_counts.items()),
                "sample_composition": ";".join(f"{key}:{value}" for key, value in sample_counts.items()),
            }
        )
    audit = pd.DataFrame(rows)
    selected, reason = choose_stress_transition_subset(audit)
    audit["eligible_main"] = (
        audit["n_cells"].ge(500)
        & audit["n_datasets"].ge(3)
        & audit["egr1_detection_rate"].gt(0)
        & audit["max_dataset_fraction"].lt(0.80)
    )
    audit["selected_main"] = audit["subset"].eq(selected)
    audit["selection_reason"] = np.where(audit["selected_main"], reason, "")
    return audit


def export_subset(
    adata: ad.AnnData,
    mask: pd.Series,
    subset: str,
    data_dir: Path,
) -> dict:
    selected = adata[mask.to_numpy(dtype=bool), :].copy()
    counts = selected.layers["counts"] if "counts" in selected.layers else selected.X
    counts = counts.tocsr() if sparse.issparse(counts) else sparse.csr_matrix(counts)
    genes_x_cells = counts.transpose().tocsr()
    subset_dir = data_dir / subset
    subset_dir.mkdir(parents=True, exist_ok=True)
    matrix_path = subset_dir / f"figure3e_egr1_{subset}_counts_genes_x_cells.mtx"
    genes_path = subset_dir / f"figure3e_egr1_{subset}_genes.tsv"
    cells_path = subset_dir / f"figure3e_egr1_{subset}_cells.tsv"
    metadata_path = subset_dir / f"figure3e_egr1_{subset}_cell_metadata.tsv"
    with matrix_path.open("wb") as handle:
        io.mmwrite(handle, genes_x_cells)
    pd.DataFrame({"gene": selected.var_names.astype(str)}).to_csv(genes_path, sep="\t", index=False)
    pd.DataFrame({"cell_id": selected.obs_names.astype(str)}).to_csv(cells_path, sep="\t", index=False)
    selected.obs.copy().assign(cell_id=selected.obs_names.astype(str)).to_csv(
        metadata_path, sep="\t", index=False
    )
    return {
        "subset": subset,
        "n_cells": int(selected.n_obs),
        "n_genes": int(selected.n_vars),
        "matrix": str(matrix_path.resolve()),
        "genes": str(genes_path.resolve()),
        "cells": str(cells_path.resolve()),
        "cell_metadata": str(metadata_path.resolve()),
    }


def run(h5ad_path: Path, data_dir: Path, metadata_dir: Path) -> dict:
    data_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(h5ad_path)
    masks = subset_masks(adata)
    audit = audit_subsets(adata, masks)
    selected_rows = audit.loc[audit["selected_main"]]
    if len(selected_rows) != 1:
        raise ValueError("Exactly one stress-transition main subset must be selected")
    selected_subset = str(selected_rows.iloc[0]["subset"])
    exports = [export_subset(adata, masks[subset], subset, data_dir) for subset in SUBSET_ORDER]
    audit_path = metadata_dir / "figure3e_egr1_subset_selection_audit.tsv"
    audit.to_csv(audit_path, sep="\t", index=False)
    review_risks = []
    if selected_subset != "stressed_injured":
        review_risks.append(
            {
                "flag": "stressed_injured_below_minimum_cells",
                "severity": "review_attention",
                "detail": (
                    f"stressed_injured contained {int(audit.loc[audit['subset'].eq('stressed_injured'), 'n_cells'].iloc[0])} "
                    "cells; stressed_injured + regenerative_progenitor was selected under the prespecified fallback."
                ),
            }
        )
    if not bool(audit["egr1_regulon_auc_available"].any()):
        review_risks.append(
            {
                "flag": "egr1_regulon_auc_missing",
                "severity": "review_attention",
                "detail": "No EGR1-specific regulon AUC was available for subset selection.",
            }
        )
    if int(adata.n_vars) != 3000:
        review_risks.append(
            {
                "flag": "network_gene_count_not_3000",
                "severity": "review_attention",
                "detail": f"Prepared CellOracle/scTenifoldKnk object contains {adata.n_vars} genes.",
            }
        )
    report = {
        "module": "Figure 3E input preparation",
        "target_tf": TARGET_TF,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_h5ad": str(h5ad_path.resolve()),
        "input_layer": "counts" if "counts" in adata.layers else "X",
        "network_gene_count": int(adata.n_vars),
        "subset_priority": SUBSET_ORDER,
        "selected_main_subset": selected_subset,
        "selected_main_state_definition": (
            "celloracle_state in {stressed_injured, regenerative_progenitor}"
            if selected_subset == "stressed_regenerative"
            else selected_subset
        ),
        "selection_audit": audit.to_dict(orient="records"),
        "exports": exports,
        "formal_main_parameters": {
            "nc_nNet": 10,
            "nc_nCells": 500,
            "multiple_seeds": [SEED, SEED + 1, SEED + 2],
            "network_genes": int(adata.n_vars),
        },
        "sensitivity_parameters": {
            "nc_nNet": 3,
            "nc_nCells": "min(500, subset n_cells)",
            "seed": SEED,
        },
        "review_risk_flags": review_risks,
        "outputs": {"subset_selection_audit": str(audit_path.resolve())},
        "runtime": {
            "python": sys.version,
            "executable": sys.executable,
            "platform": platform.platform(),
            "anndata": str(ad.__version__),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
    }
    report_path = metadata_dir / "figure3e_egr1_input_report.json"
    write_json(json_safe(report), report_path)
    adata.file.close() if getattr(adata, "file", None) and adata.file.is_open else None
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA_DIR)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = run(args.h5ad, args.data_dir, args.metadata_dir)
    print(
        json.dumps(
            {
                "selected_main_subset": report["selected_main_subset"],
                "network_gene_count": report["network_gene_count"],
                "report": str((args.metadata_dir / "figure3e_egr1_input_report.json").resolve()),
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

