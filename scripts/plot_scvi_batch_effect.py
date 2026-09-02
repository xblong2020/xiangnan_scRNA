from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot before/after batch effect for scVI integration.")
    parser.add_argument(
        "--input",
        type=Path,
        default=ROOT / "data/processed/scvi/scvi_integrated_counts_hvg.h5ad",
    )
    parser.add_argument("--out-dir", type=Path, default=ROOT / "figures/scvi_batch_effect")
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/scvi")
    parser.add_argument("--max-cells-per-sample", type=int, default=8000)
    parser.add_argument("--metric-cells", type=int, default=50000)
    parser.add_argument("--seed", type=int, default=20260601)
    return parser.parse_args()


def stratified_indices(obs: pd.DataFrame, groupby: str, max_per_group: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    selected: list[np.ndarray] = []
    for _, idx in obs.groupby(groupby, observed=True).indices.items():
        idx = np.asarray(idx)
        if len(idx) > max_per_group:
            idx = rng.choice(idx, size=max_per_group, replace=False)
        selected.append(idx)
    out = np.concatenate(selected)
    rng.shuffle(out)
    return np.sort(out)


def compute_embeddings(adata: ad.AnnData) -> tuple[ad.AnnData, np.ndarray, np.ndarray]:
    before = adata.copy()
    before.X = before.layers["counts"].copy() if "counts" in before.layers else before.X.copy()
    if sparse.issparse(before.X):
        before.X = before.X.tocsr()
    print("BEFORE normalize/log1p/PCA", flush=True)
    sc.pp.normalize_total(before, target_sum=1e4)
    sc.pp.log1p(before)
    sc.pp.pca(before, n_comps=50, random_state=20260601)
    print("BEFORE neighbors/UMAP", flush=True)
    sc.pp.neighbors(before, n_neighbors=30, n_pcs=50, random_state=20260601)
    sc.tl.umap(before, min_dist=0.3, random_state=20260601)

    after = adata.copy()
    print("AFTER scVI neighbors/UMAP", flush=True)
    sc.pp.neighbors(after, n_neighbors=30, use_rep="X_scVI", random_state=20260601)
    sc.tl.umap(after, min_dist=0.3, random_state=20260601)
    print("EMBEDDINGS complete", flush=True)

    plot_adata = ad.AnnData(
        X=sparse.csr_matrix((adata.n_obs, 1), dtype=np.float32),
        obs=adata.obs.copy(),
    )
    plot_adata.obsm["X_umap_before_logcounts"] = before.obsm["X_umap"].copy()
    plot_adata.obsm["X_umap_after_scvi"] = after.obsm["X_umap"].copy()
    plot_adata.obsm["X_pca_before_logcounts"] = before.obsm["X_pca"].copy()
    plot_adata.obsm["X_scVI"] = adata.obsm["X_scVI"].copy()
    return plot_adata, before.obsm["X_pca"].copy(), adata.obsm["X_scVI"].copy()


def color_map(values: pd.Series) -> dict[str, tuple[float, float, float, float]]:
    categories = list(pd.Categorical(values).categories)
    cmap = plt.get_cmap("tab20", max(len(categories), 1))
    return {cat: cmap(i % cmap.N) for i, cat in enumerate(categories)}


def plot_before_after(
    obs: pd.DataFrame,
    before_xy: np.ndarray,
    after_xy: np.ndarray,
    color_by: str,
    out_path: Path,
    point_size: float,
) -> None:
    colors = color_map(obs[color_by].astype(str))
    fig, axes = plt.subplots(1, 2, figsize=(13, 6), constrained_layout=True)
    for ax, xy, title in [
        (axes[0], before_xy, "Before integration: log-count PCA UMAP"),
        (axes[1], after_xy, "After integration: scVI latent UMAP"),
    ]:
        for cat, color in colors.items():
            mask = obs[color_by].astype(str).to_numpy() == cat
            ax.scatter(xy[mask, 0], xy[mask, 1], s=point_size, c=[color], label=cat, alpha=0.65, linewidths=0)
        ax.set_title(title)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines[["top", "right"]].set_visible(False)
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, markerscale=4)
    fig.suptitle(f"Batch coloring: {color_by}", y=1.02)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def normalized_batch_entropy(embedding: np.ndarray, labels: pd.Series, n_neighbors: int = 50) -> float:
    values = pd.Categorical(labels.astype(str))
    codes = values.codes
    n_categories = len(values.categories)
    if n_categories <= 1:
        return 0.0
    n_neighbors = min(n_neighbors + 1, embedding.shape[0])
    nn = NearestNeighbors(n_neighbors=n_neighbors, metric="euclidean")
    nn.fit(embedding)
    neigh = nn.kneighbors(return_distance=False)[:, 1:]
    entropies = []
    for row in neigh:
        counts = np.bincount(codes[row], minlength=n_categories).astype(float)
        probs = counts[counts > 0] / counts.sum()
        entropies.append(float(-(probs * np.log(probs)).sum() / np.log(n_categories)))
    return float(np.mean(entropies))


def compute_metrics(
    obs: pd.DataFrame,
    before_embedding: np.ndarray,
    after_embedding: np.ndarray,
    metric_cells: int,
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    n = obs.shape[0]
    metric_idx = np.arange(n)
    if n > metric_cells:
        metric_idx = rng.choice(metric_idx, size=metric_cells, replace=False)
    metric_idx = np.sort(metric_idx)
    rows = []
    for batch_key in ["dataset", "sample_id"]:
        labels = obs.iloc[metric_idx][batch_key].astype(str)
        before = before_embedding[metric_idx]
        after = after_embedding[metric_idx]
        rows.extend(
            [
                {
                    "batch_key": batch_key,
                    "embedding": "before_logcounts_pca",
                    "metric": "silhouette_score_by_batch",
                    "value": float(silhouette_score(before, labels, metric="euclidean")),
                    "interpretation": "lower_is_better_for_batch_effect",
                },
                {
                    "batch_key": batch_key,
                    "embedding": "after_scvi_latent",
                    "metric": "silhouette_score_by_batch",
                    "value": float(silhouette_score(after, labels, metric="euclidean")),
                    "interpretation": "lower_is_better_for_batch_effect",
                },
                {
                    "batch_key": batch_key,
                    "embedding": "before_logcounts_pca",
                    "metric": "mean_knn_batch_entropy",
                    "value": normalized_batch_entropy(before, labels),
                    "interpretation": "higher_is_better_for_batch_mixing",
                },
                {
                    "batch_key": batch_key,
                    "embedding": "after_scvi_latent",
                    "metric": "mean_knn_batch_entropy",
                    "value": normalized_batch_entropy(after, labels),
                    "interpretation": "higher_is_better_for_batch_mixing",
                },
            ]
        )
    return pd.DataFrame(rows)


def plot_metrics(metrics: pd.DataFrame, out_path: Path) -> None:
    labels = []
    before_values = []
    after_values = []
    for batch_key in ["dataset", "sample_id"]:
        for metric in ["silhouette_score_by_batch", "mean_knn_batch_entropy"]:
            subset = metrics[(metrics["batch_key"] == batch_key) & (metrics["metric"] == metric)]
            before = subset.loc[subset["embedding"] == "before_logcounts_pca", "value"].iloc[0]
            after = subset.loc[subset["embedding"] == "after_scvi_latent", "value"].iloc[0]
            labels.append(f"{batch_key}\n{metric}")
            before_values.append(before)
            after_values.append(after)

    x = np.arange(len(labels))
    width = 0.36
    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    ax.bar(x - width / 2, before_values, width, label="Before")
    ax.bar(x + width / 2, after_values, width, label="After scVI")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel("Metric value")
    ax.set_title("Batch effect metrics on sampled cells")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)

    backed = ad.read_h5ad(args.input, backed="r")
    idx = stratified_indices(backed.obs, "study_sample", args.max_cells_per_sample, args.seed)
    adata = backed[idx, :].to_memory()
    backed.file.close()
    print(f"SAMPLED cells={adata.n_obs} genes={adata.n_vars}", flush=True)

    plot_adata, before_embedding, after_embedding = compute_embeddings(adata)
    plot_h5ad = args.input.parent / "scvi_batch_effect_umap_sampled.h5ad"
    plot_adata.write_h5ad(plot_h5ad, compression="gzip")
    print(f"WROTE {plot_h5ad}", flush=True)

    plot_before_after(
        plot_adata.obs,
        plot_adata.obsm["X_umap_before_logcounts"],
        plot_adata.obsm["X_umap_after_scvi"],
        "dataset",
        args.out_dir / "before_after_umap_by_dataset.png",
        point_size=1.2,
    )
    plot_before_after(
        plot_adata.obs,
        plot_adata.obsm["X_umap_before_logcounts"],
        plot_adata.obsm["X_umap_after_scvi"],
        "sample_id",
        args.out_dir / "before_after_umap_by_sample.png",
        point_size=1.0,
    )

    metrics = compute_metrics(plot_adata.obs, before_embedding, after_embedding, args.metric_cells, args.seed)
    metrics_path = args.metadata_dir / "scvi_batch_effect_metrics.tsv"
    metrics.to_csv(metrics_path, sep="\t", index=False)
    plot_metrics(metrics, args.out_dir / "batch_effect_metrics.png")

    summary = {
        "input": str(args.input.resolve()),
        "sampled_h5ad": str(plot_h5ad.resolve()),
        "sampled_cells": int(plot_adata.n_obs),
        "genes": int(adata.n_vars),
        "max_cells_per_sample": int(args.max_cells_per_sample),
        "metric_cells": int(min(args.metric_cells, plot_adata.n_obs)),
        "figures": [
            str((args.out_dir / "before_after_umap_by_dataset.png").resolve()),
            str((args.out_dir / "before_after_umap_by_sample.png").resolve()),
            str((args.out_dir / "batch_effect_metrics.png").resolve()),
        ],
        "metrics": str(metrics_path.resolve()),
    }
    summary_path = args.metadata_dir / "scvi_batch_effect_plot_report.json"
    with summary_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
    print(f"WROTE {metrics_path}", flush=True)
    print(f"WROTE {summary_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
