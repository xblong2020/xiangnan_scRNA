from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export cluster-level logCPM matrix for SingleR.")
    parser.add_argument(
        "--input",
        type=Path,
        default=ROOT / "data/processed/scvi/scvi_integrated_counts_hvg.celltypist_major.h5ad",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=ROOT / "metadata/scvi/scvi_input_manifest.counts.tsv",
    )
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/celltype")
    parser.add_argument("--cluster-key", default="leiden_scvi")
    parser.add_argument("--exclude-column", default="excluded_doublet_cluster")
    parser.add_argument("--seed", type=int, default=20260601)
    return parser.parse_args()


def read_manifest(path: Path) -> pd.DataFrame:
    manifest = pd.read_csv(path, sep="\t")
    manifest = manifest.loc[manifest["include_in_scvi"].astype(str).str.lower().eq("true")].copy()
    manifest["study_sample"] = manifest["dataset"].astype(str) + "__" + manifest["label"].astype(str)
    return manifest


def as_csr(x) -> sparse.csr_matrix:
    if sparse.issparse(x):
        return x.tocsr()
    return sparse.csr_matrix(x)


def original_ids(global_index: pd.Index, study_sample: str) -> pd.Index:
    prefix = f"{study_sample}__"
    return pd.Index([idx[len(prefix) :] if str(idx).startswith(prefix) else idx for idx in global_index.astype(str)])


def main() -> int:
    args = parse_args()
    start = time.time()
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    matrix_path = args.metadata_dir / "singler_cluster_logcpm.tsv.gz"
    counts_path = args.metadata_dir / "singler_cluster_raw_sums.tsv.gz"
    cluster_meta_path = args.metadata_dir / "singler_cluster_metadata.tsv"
    report_path = args.metadata_dir / "singler_cluster_expression_report.json"

    manifest = read_manifest(args.manifest)
    print(f"READ OBS {args.input}", flush=True)
    backed = ad.read_h5ad(args.input, backed="r")
    obs = backed.obs.copy()
    backed.file.close()
    if args.cluster_key not in obs.columns:
        raise KeyError(f"{args.cluster_key!r} is not in obs")
    if args.exclude_column in obs.columns:
        keep = ~obs[args.exclude_column].astype(bool)
    else:
        keep = pd.Series(True, index=obs.index)
    obs = obs.loc[keep].copy()
    clusters = sorted(obs[args.cluster_key].astype(str).unique(), key=lambda x: int(x) if x.isdigit() else x)
    cluster_to_col = {cluster: i for i, cluster in enumerate(clusters)}
    print(f"CLUSTERS n={len(clusters)} cells={obs.shape[0]}", flush=True)

    gene_sets: list[pd.Index] = []
    for _, row in manifest.iterrows():
        a = ad.read_h5ad(row["output"], backed="r")
        gene_sets.append(pd.Index(a.var_names.astype(str)))
        a.file.close()
    genes = pd.Index(sorted(set().union(*[set(g) for g in gene_sets])))
    gene_to_idx = pd.Series(np.arange(len(genes)), index=genes)
    raw_sums = np.zeros((len(genes), len(clusters)), dtype=np.float64)

    for i, row in manifest.reset_index(drop=True).iterrows():
        study_sample = str(row["study_sample"])
        sample_obs = obs.loc[obs["study_sample"].astype(str).eq(study_sample)]
        if sample_obs.empty:
            continue
        print(f"SAMPLE {i + 1}/{manifest.shape[0]} {study_sample} cells={sample_obs.shape[0]}", flush=True)
        a = ad.read_h5ad(row["output"])
        wanted_original = original_ids(sample_obs.index, study_sample)
        present = wanted_original.intersection(a.obs_names)
        if len(present) != len(wanted_original):
            print(f"WARN {study_sample} missing={len(wanted_original) - len(present)}", flush=True)
        sample_obs = sample_obs.iloc[wanted_original.get_indexer(present)]
        x = as_csr(a[present, :].X).astype(np.float64)
        sample_clusters = sample_obs[args.cluster_key].astype(str).map(cluster_to_col).to_numpy()
        indicator = sparse.csr_matrix(
            (np.ones(x.shape[0], dtype=np.float64), (np.arange(x.shape[0]), sample_clusters)),
            shape=(x.shape[0], len(clusters)),
        )
        aggregate = (x.T @ indicator).toarray()
        loc = gene_to_idx.loc[pd.Index(a.var_names.astype(str))].to_numpy()
        raw_sums[loc, :] += aggregate

    total_counts = raw_sums.sum(axis=0)
    cpm = np.divide(raw_sums, total_counts.reshape(1, -1), out=np.zeros_like(raw_sums), where=total_counts.reshape(1, -1) > 0)
    logcpm = np.log1p(cpm * 1e6)
    columns = [f"cluster_{cluster}" for cluster in clusters]
    pd.DataFrame(logcpm, index=genes, columns=columns).to_csv(matrix_path, sep="\t", compression="gzip")
    pd.DataFrame(raw_sums, index=genes, columns=columns).to_csv(counts_path, sep="\t", compression="gzip")

    cluster_meta = (
        obs.groupby(args.cluster_key, observed=True)
        .agg(
            n_cells=("dataset", "size"),
            celltypist_major=("major_celltype", lambda x: x.astype(str).value_counts().idxmax()),
            celltypist_major_fraction=("major_celltype", lambda x: float(x.astype(str).value_counts().iloc[0] / len(x))),
            mean_celltypist_confidence=("celltypist_liver_confidence", "mean"),
            cycling_rate=("cell_cycle_phase", lambda x: float(x.astype(str).isin(["S", "G2M"]).mean())),
            predicted_doublet_rate=("predicted_doublet", lambda x: float(x.astype(bool).mean())),
        )
        .reset_index()
    )
    cluster_meta["singler_column"] = cluster_meta[args.cluster_key].astype(str).map(lambda x: f"cluster_{x}")
    cluster_meta.to_csv(cluster_meta_path, sep="\t", index=False)

    report = {
        "input": str(args.input.resolve()),
        "matrix_path": str(matrix_path.resolve()),
        "raw_sums_path": str(counts_path.resolve()),
        "cluster_metadata_path": str(cluster_meta_path.resolve()),
        "n_genes": int(len(genes)),
        "n_clusters": int(len(clusters)),
        "n_cells": int(obs.shape[0]),
        "cluster_key": args.cluster_key,
        "excluded_column": args.exclude_column,
        "elapsed_seconds": round(time.time() - start, 3),
    }
    with report_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)
    print(json.dumps(report, indent=2, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
