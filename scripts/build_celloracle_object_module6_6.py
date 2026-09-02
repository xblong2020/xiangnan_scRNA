from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

try:
    from qc_celloracle_inputs_module6_5 import read_tf_list
except ModuleNotFoundError:
    from scripts.qc_celloracle_inputs_module6_5 import read_tf_list


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_H5AD = PROJECT_ROOT / "data/processed/driver/driver_cistarget_regulon_activity.module6_3c.h5ad"
DEFAULT_EMBEDDING_H5AD = PROJECT_ROOT / "data/processed/driver/driver_hepatocyte_trajectory.module6_1.h5ad"
DEFAULT_TF_LIST = PROJECT_ROOT / "metadata/driver/celloracle_input_tfs.module6_4.txt"
DEFAULT_BASE_GRN = PROJECT_ROOT / "metadata/driver/scenic_resources/celloracle_hg38_promoter_base_grn.parquet"
DEFAULT_OUT_DIR = PROJECT_ROOT / "data/processed/driver/celloracle_module6_6"
DEFAULT_REPORT = PROJECT_ROOT / "metadata/driver/celloracle_module6_6_object_report.json"
DEFAULT_GENE_TABLE = PROJECT_ROOT / "metadata/driver/celloracle_module6_6_gene_selection.tsv"
DEFAULT_CELL_TABLE = PROJECT_ROOT / "metadata/driver/celloracle_module6_6_cell_metadata.tsv.gz"

STATE_MAP = {
    "normal_reference": "normal_reference",
    "stressed_injured": "stressed_injured",
    "regenerative_progenitor": "regenerative_progenitor",
    "proliferating_candidate": "proliferating_candidate",
    "malignant_cnv_supported": "malignant_or_malignant_like",
    "malignant_proliferation_cnv_review": "malignant_or_malignant_like",
    "malignant_like_scanvi_review": "malignant_or_malignant_like",
}


def build_celloracle_state(obs: pd.DataFrame) -> pd.Series:
    if "trajectory_role" not in obs.columns:
        raise ValueError("Required obs column missing: trajectory_role")
    state = obs["trajectory_role"].astype(str).map(STATE_MAP).fillna("other_trajectory")
    return pd.Series(state.values, index=obs.index, name="celloracle_state")


def compute_gene_metrics(matrix, genes: pd.Index) -> pd.DataFrame:
    x = matrix.tocsr() if sparse.issparse(matrix) else sparse.csr_matrix(matrix)
    n_cells = x.shape[0]
    mean_counts = np.asarray(x.mean(axis=0)).ravel()
    mean_sq = np.asarray(x.multiply(x).mean(axis=0)).ravel()
    variance = np.maximum(mean_sq - mean_counts**2, 0)
    detection = np.asarray((x > 0).sum(axis=0)).ravel() / max(n_cells, 1)
    dispersion = variance / (mean_counts + 1e-8)
    return pd.DataFrame(
        {
            "gene": genes.astype(str).to_numpy(),
            "mean_counts": mean_counts,
            "variance": variance,
            "detection_rate": detection,
            "dispersion_score": dispersion,
        }
    )


def select_celloracle_genes(
    metrics: pd.DataFrame,
    input_tfs: list[str],
    max_genes: int,
    min_detection_rate: float,
) -> list[str]:
    gene_set = set(metrics["gene"].astype(str))
    forced = [tf for tf in input_tfs if tf in gene_set]
    forced_set = set(forced)

    eligible = metrics.loc[
        (metrics["detection_rate"] >= min_detection_rate)
        & (~metrics["gene"].astype(str).isin(forced_set))
    ].copy()
    eligible = eligible.sort_values(
        ["dispersion_score", "detection_rate", "mean_counts", "gene"],
        ascending=[False, False, False, True],
    )
    remaining = max(max_genes - len(forced), 0)
    selected = forced + eligible["gene"].astype(str).head(remaining).tolist()
    return selected[:max_genes]


def validate_embedding_alignment(cell_names: list[str], embedding: pd.DataFrame) -> np.ndarray:
    missing = [cell for cell in cell_names if cell not in embedding.index]
    if missing:
        raise ValueError(f"Embedding is missing {len(missing)} cells; first missing cell: {missing[0]}")
    aligned = embedding.loc[cell_names, ["UMAP1", "UMAP2"]].to_numpy(dtype=float)
    if not np.isfinite(aligned).all():
        raise ValueError("Embedding contains non-finite values after alignment")
    return aligned


def load_embedding_table(path: Path, obsm_key: str = "X_umap") -> pd.DataFrame:
    # Read only the pieces we need. This avoids incompatibilities when an older
    # anndata reader encounters newer nullable fields in unrelated /uns entries.
    with h5py.File(path, "r") as handle:
        key_path = f"obsm/{obsm_key}"
        if key_path not in handle:
            available = list(handle["obsm"].keys()) if "obsm" in handle else []
            raise ValueError(f"Embedding key {obsm_key!r} not found. Available keys: {available}")
        raw_index = handle["obs/_index"][:]
        index = [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in raw_index]
        emb = np.asarray(handle[key_path])
    table = pd.DataFrame(emb[:, :2], index=index, columns=["UMAP1", "UMAP2"])
    return table


def make_celloracle_input_anndata(
    input_h5ad: Path,
    embedding_h5ad: Path,
    tf_list: Path,
    max_genes: int,
    min_detection_rate: float,
    embedding_key: str,
) -> tuple[ad.AnnData, pd.DataFrame, pd.DataFrame]:
    input_tfs = read_tf_list(tf_list)
    source = ad.read_h5ad(input_h5ad)
    if "counts" not in source.layers:
        raise ValueError("Input h5ad must contain layers['counts']")

    counts = source.layers["counts"]
    gene_metrics = compute_gene_metrics(counts, source.var_names)
    selected_genes = select_celloracle_genes(
        gene_metrics,
        input_tfs=input_tfs,
        max_genes=max_genes,
        min_detection_rate=min_detection_rate,
    )

    selected_mask = source.var_names.astype(str).isin(selected_genes)
    prepared = source[:, selected_mask].copy()
    prepared.X = prepared.layers["counts"].copy()
    if "counts" not in prepared.layers:
        prepared.layers["counts"] = prepared.X.copy()

    celloracle_state = build_celloracle_state(prepared.obs)
    prepared.obs["celloracle_state"] = pd.Categorical(celloracle_state)
    prepared.obs["celloracle_main_strict"] = prepared.obs["cellrank_fate_prob_cnv_supported_malignant"].notna()
    prepared.obs["celloracle_perturbation_target"] = prepared.obs["driver_primary_cnv_evidence_tier"].astype(str)

    embedding_table = load_embedding_table(embedding_h5ad, obsm_key=embedding_key)
    prepared.obsm["X_celloracle_umap"] = validate_embedding_alignment(
        prepared.obs_names.astype(str).tolist(),
        embedding_table,
    )

    gene_metrics = gene_metrics.copy()
    gene_metrics["selected_for_celloracle"] = gene_metrics["gene"].astype(str).isin(set(selected_genes))
    gene_metrics["forced_input_tf"] = gene_metrics["gene"].astype(str).isin(set(input_tfs))
    gene_metrics["selection_rank"] = np.nan
    rank_map = {gene: i + 1 for i, gene in enumerate(selected_genes)}
    gene_metrics.loc[gene_metrics["selected_for_celloracle"], "selection_rank"] = (
        gene_metrics.loc[gene_metrics["selected_for_celloracle"], "gene"].map(rank_map)
    )

    cell_columns = [
        "celloracle_state",
        "celloracle_main_strict",
        "celloracle_perturbation_target",
        "dataset",
        "sample_id",
        "trajectory_role",
        "trajectory_root_end_role",
        "cellrank_fate_prob_cnv_supported_malignant",
        "driver_main_strict__pseudotime_mean",
        "driver_primary_cnv_evidence_tier",
    ]
    cell_table = prepared.obs[[col for col in cell_columns if col in prepared.obs.columns]].copy()
    cell_table.insert(0, "cell_id", prepared.obs_names.astype(str))

    return prepared, gene_metrics, cell_table


def build_and_save_oracle(
    prepared: ad.AnnData,
    base_grn: Path,
    oracle_path: Path,
    cluster_column: str = "celloracle_state",
    embedding_name: str = "X_celloracle_umap",
) -> dict:
    import celloracle as co

    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(
        adata=prepared,
        cluster_column_name=cluster_column,
        embedding_name=embedding_name,
        transform="natural_log",
    )
    oracle.import_TF_data(TF_info_matrix_path=str(base_grn))
    oracle_path.parent.mkdir(parents=True, exist_ok=True)
    oracle.to_hdf5(str(oracle_path))

    return {
        "celloracle_version": getattr(co, "__version__", None),
        "oracle_path": str(oracle_path),
        "n_oracle_cells": int(oracle.adata.shape[0]),
        "n_oracle_genes": int(oracle.adata.shape[1]),
        "n_tfdict_targets_in_expression": int(oracle.adata.var["isin_TFdict_targets"].sum()),
        "n_tfdict_regulators_in_expression": int(oracle.adata.var["isin_TFdict_regulators"].sum()),
        "n_high_var_genes": int(len(getattr(oracle, "high_var_genes", []))),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 6.6 build CellOracle Oracle object")
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT_H5AD)
    parser.add_argument("--embedding-h5ad", type=Path, default=DEFAULT_EMBEDDING_H5AD)
    parser.add_argument("--tf-list", type=Path, default=DEFAULT_TF_LIST)
    parser.add_argument("--base-grn", type=Path, default=DEFAULT_BASE_GRN)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--gene-table", type=Path, default=DEFAULT_GENE_TABLE)
    parser.add_argument("--cell-table", type=Path, default=DEFAULT_CELL_TABLE)
    parser.add_argument("--max-genes", type=int, default=3000)
    parser.add_argument("--min-detection-rate", type=float, default=0.01)
    parser.add_argument("--embedding-key", default="X_umap")
    parser.add_argument("--skip-oracle-save", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    prepared_h5ad = args.out_dir / "celloracle_module6_6_input.h5ad"
    oracle_path = args.out_dir / "celloracle_module6_6.celloracle.oracle"

    prepared, gene_metrics, cell_table = make_celloracle_input_anndata(
        input_h5ad=args.input_h5ad,
        embedding_h5ad=args.embedding_h5ad,
        tf_list=args.tf_list,
        max_genes=args.max_genes,
        min_detection_rate=args.min_detection_rate,
        embedding_key=args.embedding_key,
    )
    prepared.write_h5ad(prepared_h5ad, compression="gzip")

    args.gene_table.parent.mkdir(parents=True, exist_ok=True)
    args.cell_table.parent.mkdir(parents=True, exist_ok=True)
    gene_metrics.to_csv(args.gene_table, sep="\t", index=False)
    cell_table.to_csv(args.cell_table, sep="\t", index=False)

    oracle_summary = {}
    if not args.skip_oracle_save:
        oracle_summary = build_and_save_oracle(
            prepared=prepared,
            base_grn=args.base_grn,
            oracle_path=oracle_path,
        )

    selected_gene_table = gene_metrics.loc[gene_metrics["selected_for_celloracle"]].copy()
    input_tfs = read_tf_list(args.tf_list)
    selected_input_tfs = sorted(set(input_tfs).intersection(set(selected_gene_table["gene"].astype(str))))

    report = {
        "module": "6.6",
        "method": "Build CellOracle Oracle object from Module 6.3c driver-union cells",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_h5ad": str(args.input_h5ad),
        "embedding_h5ad": str(args.embedding_h5ad),
        "embedding_key": args.embedding_key,
        "base_grn": str(args.base_grn),
        "tf_list": str(args.tf_list),
        "input_shape": {"n_cells": 9512, "n_genes": 12000},
        "prepared_shape": {"n_cells": int(prepared.shape[0]), "n_genes": int(prepared.shape[1])},
        "max_genes": int(args.max_genes),
        "min_detection_rate": float(args.min_detection_rate),
        "n_input_tfs": len(input_tfs),
        "n_input_tfs_retained": len(selected_input_tfs),
        "input_tfs_retained": selected_input_tfs,
        "state_counts": prepared.obs["celloracle_state"].value_counts().to_dict(),
        "n_main_strict_cells": int(prepared.obs["celloracle_main_strict"].sum()),
        "outputs": {
            "prepared_h5ad": str(prepared_h5ad),
            "oracle": str(oracle_path) if not args.skip_oracle_save else None,
            "gene_table": str(args.gene_table),
            "cell_table": str(args.cell_table),
            "report": str(args.report),
        },
        "oracle_summary": oracle_summary,
        "python_runtime": {
            "executable": sys.executable,
            "version": sys.version,
            "platform": platform.platform(),
        },
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(
        {
            "prepared_h5ad": str(prepared_h5ad),
            "oracle": str(oracle_path) if not args.skip_oracle_save else None,
            "n_cells": int(prepared.shape[0]),
            "n_genes": int(prepared.shape[1]),
            "n_input_tfs_retained": len(selected_input_tfs),
            "report": str(args.report),
        },
        ensure_ascii=False,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
