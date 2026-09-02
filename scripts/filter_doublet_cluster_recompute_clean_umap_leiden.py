from __future__ import annotations

import argparse
import json
import time
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Remove a suspected doublet Leiden cluster and recompute clean scVI neighbors/UMAP/Leiden."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=ROOT / "data/processed/scvi/scvi_integrated_counts_hvg.scvi_doublet_cell_cycle.h5ad",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "data/processed/scvi/scvi_integrated_counts_hvg.clean_no_leiden16.h5ad",
    )
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/scvi")
    parser.add_argument("--figures-dir", type=Path, default=ROOT / "figures/scvi_clean_no_leiden16")
    parser.add_argument("--exclude-cluster-key", default="leiden_scvi")
    parser.add_argument("--exclude-cluster", action="append", default=["16"])
    parser.add_argument("--cluster-key", default="leiden_scvi_clean")
    parser.add_argument("--n-neighbors", type=int, default=30)
    parser.add_argument("--min-dist", type=float, default=0.3)
    parser.add_argument("--resolution", type=float, default=1.0)
    parser.add_argument("--max-plot-cells", type=int, default=250000)
    parser.add_argument("--seed", type=int, default=20260601)
    return parser.parse_args()


def sort_cluster_frame(df: pd.DataFrame, cluster_key: str) -> pd.DataFrame:
    values = df[cluster_key].astype(str)
    if values.str.fullmatch(r"\d+").all():
        return df.assign(_cluster_sort=values.astype(int)).sort_values("_cluster_sort").drop(columns="_cluster_sort")
    return df.sort_values(cluster_key)


def export_tables(
    adata: ad.AnnData,
    excluded_obs: pd.DataFrame,
    metadata_dir: Path,
    cluster_key: str,
) -> dict[str, str]:
    metadata_dir.mkdir(parents=True, exist_ok=True)
    umap_path = metadata_dir / "scvi_clean_no_leiden16_umap.tsv.gz"
    clusters_path = metadata_dir / "scvi_clean_no_leiden16_clusters.tsv.gz"
    counts_path = metadata_dir / "scvi_clean_no_leiden16_cluster_counts.tsv"
    excluded_path = metadata_dir / "scvi_clean_excluded_leiden16_cells.tsv.gz"
    diagnostics_path = metadata_dir / "scvi_clean_no_leiden16_cluster_diagnostics.tsv"

    pd.DataFrame(
        adata.obsm["X_umap"],
        index=adata.obs_names,
        columns=["UMAP_1", "UMAP_2"],
    ).to_csv(umap_path, sep="\t", compression="gzip")

    obs_cols = [cluster_key]
    for col in [
        "leiden_scvi",
        "dataset",
        "sample_id",
        "study_sample",
        "doublet_score",
        "predicted_doublet",
        "cell_cycle_S_score",
        "cell_cycle_G2M_score",
        "cell_cycle_phase",
    ]:
        if col in adata.obs.columns and col not in obs_cols:
            obs_cols.append(col)
    adata.obs[obs_cols].to_csv(clusters_path, sep="\t", compression="gzip")
    excluded_obs.to_csv(excluded_path, sep="\t", compression="gzip")

    counts = (
        adata.obs[cluster_key]
        .astype(str)
        .value_counts()
        .rename_axis(cluster_key)
        .reset_index(name="n_cells")
    )
    counts = sort_cluster_frame(counts, cluster_key)
    counts.to_csv(counts_path, sep="\t", index=False)

    diagnostics = summarize_clean_clusters(adata.obs, cluster_key)
    diagnostics.to_csv(diagnostics_path, sep="\t", index=False)

    return {
        "umap_path": str(umap_path.resolve()),
        "clusters_path": str(clusters_path.resolve()),
        "cluster_counts_path": str(counts_path.resolve()),
        "excluded_cells_path": str(excluded_path.resolve()),
        "cluster_diagnostics_path": str(diagnostics_path.resolve()),
    }


def summarize_clean_clusters(obs: pd.DataFrame, cluster_key: str) -> pd.DataFrame:
    rows = []
    for cluster, sub in obs.groupby(cluster_key, observed=True):
        phases = sub["cell_cycle_phase"].astype(str) if "cell_cycle_phase" in sub else pd.Series([], dtype=str)
        n = int(sub.shape[0])
        predicted = sub["predicted_doublet"].astype(bool) if "predicted_doublet" in sub else pd.Series(False, index=sub.index)
        doublet_score = sub["doublet_score"].astype(float) if "doublet_score" in sub else pd.Series(np.nan, index=sub.index)
        s_score = sub["cell_cycle_S_score"].astype(float) if "cell_cycle_S_score" in sub else pd.Series(np.nan, index=sub.index)
        g2m_score = sub["cell_cycle_G2M_score"].astype(float) if "cell_cycle_G2M_score" in sub else pd.Series(np.nan, index=sub.index)
        rows.append(
            {
                cluster_key: str(cluster),
                "n_cells": n,
                "predicted_doublet_rate": float(predicted.mean()) if n else np.nan,
                "predicted_doublets": int(predicted.sum()),
                "doublet_score_mean": float(doublet_score.mean()),
                "doublet_score_p95": float(doublet_score.quantile(0.95)),
                "cycling_rate": float(phases.isin(["S", "G2M"]).mean()) if n else np.nan,
                "s_phase_rate": float((phases == "S").mean()) if n else np.nan,
                "g2m_phase_rate": float((phases == "G2M").mean()) if n else np.nan,
                "dominant_phase": str(phases.value_counts().idxmax()) if n and not phases.empty else "unknown",
                "mean_S_score": float(s_score.mean()),
                "mean_G2M_score": float(g2m_score.mean()),
            }
        )
    return sort_cluster_frame(pd.DataFrame(rows), cluster_key)


def plot_clean_umaps(adata: ad.AnnData, figures_dir: Path, max_plot_cells: int, seed: int) -> list[str]:
    figures_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    idx = np.arange(adata.n_obs)
    if adata.n_obs > max_plot_cells:
        idx = np.sort(rng.choice(idx, size=max_plot_cells, replace=False))
    obs = adata.obs.iloc[idx]
    xy = np.asarray(adata.obsm["X_umap"])[idx]
    paths: list[str] = []

    fig, ax = plt.subplots(figsize=(8, 7), constrained_layout=True)
    labels = pd.Categorical(obs["leiden_scvi_clean"].astype(str))
    cmap = plt.get_cmap("tab20", max(len(labels.categories), 1))
    for i, cat in enumerate(labels.categories):
        mask = labels == cat
        ax.scatter(xy[mask, 0], xy[mask, 1], s=0.6, c=[cmap(i % cmap.N)], linewidths=0, alpha=0.65)
    ax.set_title("Clean scVI UMAP by Leiden")
    ax.set_xticks([])
    ax.set_yticks([])
    path = figures_dir / "clean_umap_by_leiden.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    paths.append(str(path.resolve()))

    phase_colors = {"G1": "#7f7f7f", "S": "#1f77b4", "G2M": "#ff7f0e", "unknown": "#cccccc"}
    fig, ax = plt.subplots(figsize=(8, 7), constrained_layout=True)
    for phase, color in phase_colors.items():
        mask = obs["cell_cycle_phase"].astype(str).to_numpy() == phase
        if mask.any():
            ax.scatter(xy[mask, 0], xy[mask, 1], s=0.7, c=color, label=phase, linewidths=0, alpha=0.65)
    ax.legend(frameon=False, markerscale=5)
    ax.set_title("Clean scVI UMAP by cell-cycle phase")
    ax.set_xticks([])
    ax.set_yticks([])
    path = figures_dir / "clean_umap_by_cell_cycle_phase.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    paths.append(str(path.resolve()))

    fig, ax = plt.subplots(figsize=(8, 7), constrained_layout=True)
    mask = obs["predicted_doublet"].astype(bool).to_numpy()
    ax.scatter(xy[~mask, 0], xy[~mask, 1], s=0.45, c="#d0d0d0", linewidths=0, alpha=0.45)
    ax.scatter(xy[mask, 0], xy[mask, 1], s=1.1, c="#d62728", linewidths=0, alpha=0.8)
    ax.set_title("Clean scVI UMAP with residual predicted doublets")
    ax.set_xticks([])
    ax.set_yticks([])
    path = figures_dir / "clean_umap_residual_predicted_doublets.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    paths.append(str(path.resolve()))

    return paths


def plot_cluster_diagnostics(diagnostics: pd.DataFrame, cluster_key: str, figures_dir: Path) -> str:
    figures_dir.mkdir(parents=True, exist_ok=True)
    labels = diagnostics[cluster_key].astype(str).to_numpy()
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(14, 5), constrained_layout=True)
    ax.bar(x - 0.2, diagnostics["predicted_doublet_rate"], width=0.4, label="residual_doublet_rate")
    ax.bar(x + 0.2, diagnostics["cycling_rate"], width=0.4, label="cycling_rate")
    ax.axhline(0.15, color="#d62728", linestyle="--", linewidth=1)
    ax.axhline(0.50, color="#ff7f0e", linestyle="--", linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=90)
    ax.set_xlabel(cluster_key)
    ax.set_ylabel("fraction of cells")
    ax.set_ylim(0, min(1.0, max(0.1, diagnostics[["predicted_doublet_rate", "cycling_rate"]].max().max() * 1.15)))
    ax.legend(frameon=False)
    path = figures_dir / "clean_cluster_doublet_cycling_rates.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return str(path.resolve())


def main() -> int:
    args = parse_args()
    start = time.time()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    args.figures_dir.mkdir(parents=True, exist_ok=True)

    print(f"READ {args.input}", flush=True)
    adata = ad.read_h5ad(args.input)
    if "X_scVI" not in adata.obsm:
        raise KeyError("adata.obsm['X_scVI'] is required but was not found")
    if args.exclude_cluster_key not in adata.obs.columns:
        raise KeyError(f"{args.exclude_cluster_key!r} is required in adata.obs")

    exclude = set(map(str, args.exclude_cluster))
    old_clusters = adata.obs[args.exclude_cluster_key].astype(str)
    keep = ~old_clusters.isin(exclude)
    excluded_obs = adata.obs.loc[~keep].copy()
    print(
        f"FILTER {args.exclude_cluster_key} exclude={sorted(exclude)} cells {adata.n_obs}->{int(keep.sum())}",
        flush=True,
    )

    if "X_umap" in adata.obsm:
        adata.obsm["X_umap_before_doublet_filter"] = adata.obsm["X_umap"].copy()
    adata = adata[keep, :].copy()

    print(f"NEIGHBORS n_neighbors={args.n_neighbors} use_rep=X_scVI", flush=True)
    sc.pp.neighbors(
        adata,
        n_neighbors=args.n_neighbors,
        use_rep="X_scVI",
        random_state=args.seed,
    )
    print(f"UMAP min_dist={args.min_dist}", flush=True)
    sc.tl.umap(adata, min_dist=args.min_dist, random_state=args.seed)
    print(f"LEIDEN resolution={args.resolution} key_added={args.cluster_key}", flush=True)
    sc.tl.leiden(
        adata,
        resolution=args.resolution,
        key_added=args.cluster_key,
        random_state=args.seed,
        flavor="igraph",
        directed=False,
    )

    adata.uns["clean_doublet_filter"] = {
        "source_input": str(args.input.resolve()),
        "exclude_cluster_key": args.exclude_cluster_key,
        "excluded_clusters": sorted(exclude),
        "excluded_cells": int(excluded_obs.shape[0]),
        "kept_cells": int(adata.n_obs),
        "cell_cycle_columns_retained": [
            col
            for col in ["cell_cycle_S_score", "cell_cycle_G2M_score", "cell_cycle_phase"]
            if col in adata.obs.columns
        ],
    }

    print(f"WRITE {args.output}", flush=True)
    adata.write_h5ad(args.output, compression="gzip")
    exported = export_tables(adata, excluded_obs, args.metadata_dir, args.cluster_key)
    diagnostics = pd.read_csv(exported["cluster_diagnostics_path"], sep="\t")
    figures = plot_clean_umaps(adata, args.figures_dir, args.max_plot_cells, args.seed)
    figures.append(plot_cluster_diagnostics(diagnostics, args.cluster_key, args.figures_dir))

    report = {
        "input": str(args.input.resolve()),
        "output": str(args.output.resolve()),
        "n_obs_before": int(keep.shape[0]),
        "n_obs_after": int(adata.n_obs),
        "excluded_cells": int(excluded_obs.shape[0]),
        "excluded_fraction": float(excluded_obs.shape[0] / keep.shape[0]),
        "n_vars": int(adata.n_vars),
        "use_rep": "X_scVI",
        "n_neighbors": int(args.n_neighbors),
        "min_dist": float(args.min_dist),
        "resolution": float(args.resolution),
        "exclude_cluster_key": args.exclude_cluster_key,
        "excluded_clusters": sorted(exclude),
        "cluster_key": args.cluster_key,
        "n_clusters": int(adata.obs[args.cluster_key].nunique()),
        "residual_predicted_doublet_rate": float(adata.obs["predicted_doublet"].astype(bool).mean())
        if "predicted_doublet" in adata.obs
        else None,
        "global_cycling_rate": float(adata.obs["cell_cycle_phase"].astype(str).isin(["S", "G2M"]).mean())
        if "cell_cycle_phase" in adata.obs
        else None,
        "figures": figures,
        "seed": int(args.seed),
        "scanpy_version": version("scanpy"),
        "elapsed_seconds": round(time.time() - start, 3),
    }
    report.update(exported)
    report_path = args.metadata_dir / "scvi_clean_no_leiden16_report.json"
    with report_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)
    print(f"WROTE {report_path}", flush=True)
    print(json.dumps(report, indent=2, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
