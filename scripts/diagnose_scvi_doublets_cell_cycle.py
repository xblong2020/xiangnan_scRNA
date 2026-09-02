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
from pynndescent import NNDescent
from scipy import sparse
from sklearn.decomposition import TruncatedSVD


ROOT = Path(__file__).resolve().parents[1]


S_GENES = [
    "MCM5",
    "PCNA",
    "TYMS",
    "FEN1",
    "MCM2",
    "MCM4",
    "RRM1",
    "UNG",
    "GINS2",
    "MCM6",
    "CDCA7",
    "DTL",
    "PRIM1",
    "UHRF1",
    "MLF1IP",
    "HELLS",
    "RFC2",
    "RPA2",
    "NASP",
    "RAD51AP1",
    "GMNN",
    "WDR76",
    "SLBP",
    "CCNE2",
    "UBR7",
    "POLD3",
    "MSH2",
    "ATAD2",
    "RAD51",
    "RRM2",
    "CDC45",
    "CDC6",
    "EXO1",
    "TIPIN",
    "DSCC1",
    "BLM",
    "CASP8AP2",
    "USP1",
    "CLSPN",
    "POLA1",
    "CHAF1B",
    "BRIP1",
    "E2F8",
]

G2M_GENES = [
    "HMGB2",
    "CDK1",
    "NUSAP1",
    "UBE2C",
    "BIRC5",
    "TPX2",
    "TOP2A",
    "NDC80",
    "CKS2",
    "NUF2",
    "CKS1B",
    "MKI67",
    "TMPO",
    "CENPF",
    "TACC3",
    "FAM64A",
    "SMC4",
    "CCNB2",
    "CKAP2L",
    "CKAP2",
    "AURKB",
    "BUB1",
    "KIF11",
    "ANP32E",
    "TUBB4B",
    "GTSE1",
    "KIF20B",
    "HJURP",
    "CDCA3",
    "HN1",
    "CDC20",
    "TTK",
    "CDC25C",
    "KIF2C",
    "RANGAP1",
    "NCAPD2",
    "DLGAP5",
    "CDCA2",
    "CDCA8",
    "ECT2",
    "KIF23",
    "HMMR",
    "AURKA",
    "PSRC1",
    "ANLN",
    "LBR",
    "CKAP5",
    "CENPE",
    "CTCF",
    "NEK2",
    "G2E3",
    "GAS2L3",
    "CBX5",
    "CENPA",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run per-sample doublet detection and cell-cycle diagnostics on scVI Leiden clusters."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=ROOT / "data/processed/scvi/scvi_integrated_counts_hvg.scvi_neighbors_umap_leiden.h5ad",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "data/processed/scvi/scvi_integrated_counts_hvg.scvi_doublet_cell_cycle.h5ad",
    )
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/scvi")
    parser.add_argument("--figures-dir", type=Path, default=ROOT / "figures/scvi_doublet_cell_cycle")
    parser.add_argument("--cluster-key", default="leiden_scvi")
    parser.add_argument("--expected-doublet-rate", type=float, default=0.05)
    parser.add_argument("--sim-doublet-ratio", type=float, default=1.5)
    parser.add_argument("--max-sim-doublets", type=int, default=100000)
    parser.add_argument("--n-top-genes", type=int, default=3000)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=30)
    parser.add_argument("--doublet-cluster-rate", type=float, default=0.15)
    parser.add_argument("--doublet-cluster-fold", type=float, default=2.0)
    parser.add_argument("--cycling-cluster-rate", type=float, default=0.50)
    parser.add_argument("--min-cluster-cells", type=int, default=100)
    parser.add_argument("--max-plot-cells", type=int, default=250000)
    parser.add_argument("--seed", type=int, default=20260601)
    parser.add_argument("--only-study-sample", action="append", default=None)
    parser.add_argument("--max-samples", type=int, default=None)
    parser.add_argument("--skip-annotated-h5ad", action="store_true")
    return parser.parse_args()


def as_csr(x) -> sparse.csr_matrix:
    if sparse.issparse(x):
        return x.tocsr()
    return sparse.csr_matrix(x)


def normalize_log_csr(x: sparse.csr_matrix, target_sum: float = 1e4) -> sparse.csr_matrix:
    x = x.astype(np.float32).tocsr(copy=True)
    totals = np.asarray(x.sum(axis=1)).ravel()
    scale = np.divide(target_sum, totals, out=np.zeros_like(totals, dtype=np.float32), where=totals > 0)
    x.data *= np.repeat(scale, np.diff(x.indptr))
    np.log1p(x.data, out=x.data)
    return x


def select_hvg_indices(x: sparse.csr_matrix, genes: pd.Index, n_top_genes: int) -> np.ndarray:
    n_cells = x.shape[0]
    sums = np.asarray(x.sum(axis=0)).ravel()
    squares = np.asarray(x.multiply(x).sum(axis=0)).ravel()
    detected = np.asarray((x > 0).sum(axis=0)).ravel()
    mean = sums / max(n_cells, 1)
    var = np.maximum((squares / max(n_cells, 1)) - np.square(mean), 0)
    dispersion = var / np.maximum(mean, 1e-8)
    gene_values = genes.astype(str).to_numpy(dtype=str)
    valid = (detected >= 3) & (mean > 0) & (~np.char.startswith(gene_values, "MT-"))
    if valid.sum() == 0:
        raise ValueError("No valid genes available for doublet HVG selection")
    valid_idx = np.where(valid)[0]
    ranked = valid_idx[np.argsort(dispersion[valid_idx])[::-1]]
    return ranked[: min(n_top_genes, ranked.size)]


def simulated_doublet_knn(
    x: sparse.csr_matrix,
    genes: pd.Index,
    rng: np.random.Generator,
    expected_doublet_rate: float,
    sim_doublet_ratio: float,
    max_sim_doublets: int,
    n_top_genes: int,
    n_pcs: int,
    n_neighbors: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, float, dict[str, int | float]]:
    hvg_idx = select_hvg_indices(x, genes, n_top_genes)
    observed = x[:, hvg_idx].astype(np.float32).tocsr()
    n_obs = observed.shape[0]
    n_sim = min(int(np.ceil(n_obs * sim_doublet_ratio)), max_sim_doublets)
    left = rng.integers(0, n_obs, size=n_sim)
    right = rng.integers(0, n_obs, size=n_sim)
    simulated = (observed[left, :] + observed[right, :]).tocsr()
    combined = sparse.vstack([observed, simulated], format="csr")
    combined = normalize_log_csr(combined)

    n_components = max(2, min(n_pcs, combined.shape[0] - 1, combined.shape[1] - 1))
    embedding = TruncatedSVD(n_components=n_components, random_state=seed).fit_transform(combined)
    embedding = np.asarray(embedding, dtype=np.float32)
    k = min(n_neighbors + 1, embedding.shape[0])
    index = NNDescent(embedding, n_neighbors=k, metric="euclidean", random_state=seed, n_jobs=-1)
    neighbor_idx = index.neighbor_graph[0]

    scores = np.zeros(n_obs, dtype=np.float32)
    for i in range(n_obs):
        row = neighbor_idx[i]
        row = row[row != i][:n_neighbors]
        scores[i] = np.mean(row >= n_obs) if row.size else 0.0

    threshold = float(np.quantile(scores, max(0.0, min(1.0, 1.0 - expected_doublet_rate))))
    predicted = scores >= threshold
    stats = {
        "doublet_hvg_genes": int(len(hvg_idx)),
        "simulated_doublets": int(n_sim),
        "doublet_n_pcs": int(n_components),
        "doublet_n_neighbors": int(k - 1),
    }
    return scores, predicted, threshold, stats


def score_cell_cycle(adata: ad.AnnData) -> tuple[pd.DataFrame, dict[str, int]]:
    genes = pd.Index(adata.var_names.astype(str))
    s_present = [gene for gene in S_GENES if gene in genes]
    g2m_present = [gene for gene in G2M_GENES if gene in genes]
    if len(s_present) < 5 or len(g2m_present) < 5:
        out = pd.DataFrame(index=adata.obs_names)
        out["cell_cycle_S_score"] = np.nan
        out["cell_cycle_G2M_score"] = np.nan
        out["cell_cycle_phase"] = "unknown"
        return out, {"s_genes_used": len(s_present), "g2m_genes_used": len(g2m_present)}

    cc = ad.AnnData(X=as_csr(adata.X).copy(), obs=pd.DataFrame(index=adata.obs_names), var=adata.var.copy())
    sc.pp.normalize_total(cc, target_sum=1e4)
    sc.pp.log1p(cc)
    sc.tl.score_genes_cell_cycle(cc, s_genes=s_present, g2m_genes=g2m_present)
    out = pd.DataFrame(index=adata.obs_names)
    out["cell_cycle_S_score"] = cc.obs["S_score"].astype(float).to_numpy()
    out["cell_cycle_G2M_score"] = cc.obs["G2M_score"].astype(float).to_numpy()
    out["cell_cycle_phase"] = cc.obs["phase"].astype(str).to_numpy()
    return out, {"s_genes_used": len(s_present), "g2m_genes_used": len(g2m_present)}


def sample_table(obs: pd.DataFrame, cluster_key: str) -> pd.DataFrame:
    cols = ["dataset", "sample_id", "study_sample", "source_h5ad", cluster_key]
    missing = [col for col in cols if col not in obs.columns]
    if missing:
        raise KeyError(f"Missing required obs columns: {missing}")
    return obs[["dataset", "sample_id", "study_sample", "source_h5ad"]].drop_duplicates().sort_values(
        ["dataset", "sample_id"]
    )


def original_ids_for_sample(obs: pd.DataFrame, study_sample: str) -> pd.Index:
    prefix = f"{study_sample}__"
    idx = obs.index[obs["study_sample"].astype(str) == study_sample]
    return pd.Index([name[len(prefix) :] if name.startswith(prefix) else name for name in idx])


def process_one_sample(
    row: pd.Series,
    integrated_obs: pd.DataFrame,
    args: argparse.Namespace,
    sample_number: int,
    n_samples: int,
) -> tuple[pd.DataFrame, dict[str, object]]:
    start = time.time()
    dataset = str(row["dataset"])
    sample_id = str(row["sample_id"])
    study_sample = str(row["study_sample"])
    source = Path(str(row["source_h5ad"]))
    print(f"SAMPLE {sample_number}/{n_samples} {study_sample} READ {source}", flush=True)

    expected_ids = original_ids_for_sample(integrated_obs, study_sample)
    adata = ad.read_h5ad(source)
    present = expected_ids.intersection(adata.obs_names)
    missing_cells = int(len(expected_ids) - len(present))
    if missing_cells:
        print(f"WARN {study_sample} missing_cells={missing_cells}", flush=True)
    adata = adata[present, :].copy()
    x = as_csr(adata.X)

    rng = np.random.default_rng(args.seed + sample_number)
    scores, predicted, threshold, doublet_stats = simulated_doublet_knn(
        x=x,
        genes=pd.Index(adata.var_names.astype(str)),
        rng=rng,
        expected_doublet_rate=args.expected_doublet_rate,
        sim_doublet_ratio=args.sim_doublet_ratio,
        max_sim_doublets=args.max_sim_doublets,
        n_top_genes=args.n_top_genes,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors,
        seed=args.seed + sample_number,
    )
    cycle, cycle_stats = score_cell_cycle(adata)

    diag = pd.DataFrame(index=adata.obs_names)
    diag["cell_id"] = [f"{study_sample}__{idx}" for idx in diag.index.astype(str)]
    diag["dataset"] = dataset
    diag["sample_id"] = sample_id
    diag["study_sample"] = study_sample
    diag["doublet_score"] = scores
    diag["predicted_doublet"] = predicted
    diag["doublet_threshold"] = threshold
    diag["doublet_method"] = "simulated_doublet_knn"
    diag["cell_cycle_S_score"] = cycle["cell_cycle_S_score"].to_numpy()
    diag["cell_cycle_G2M_score"] = cycle["cell_cycle_G2M_score"].to_numpy()
    diag["cell_cycle_phase"] = cycle["cell_cycle_phase"].to_numpy()
    diag = diag.set_index("cell_id", drop=True)

    phase_counts = diag["cell_cycle_phase"].value_counts()
    summary = {
        "dataset": dataset,
        "sample_id": sample_id,
        "study_sample": study_sample,
        "source_h5ad": str(source.resolve()),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "missing_integrated_cells": missing_cells,
        "expected_doublet_rate": float(args.expected_doublet_rate),
        "predicted_doublet_rate": float(predicted.mean()) if len(predicted) else np.nan,
        "doublet_score_mean": float(np.mean(scores)) if len(scores) else np.nan,
        "doublet_score_median": float(np.median(scores)) if len(scores) else np.nan,
        "doublet_score_p95": float(np.quantile(scores, 0.95)) if len(scores) else np.nan,
        "doublet_threshold": float(threshold),
        "g1_rate": float(phase_counts.get("G1", 0) / max(adata.n_obs, 1)),
        "s_rate": float(phase_counts.get("S", 0) / max(adata.n_obs, 1)),
        "g2m_rate": float(phase_counts.get("G2M", 0) / max(adata.n_obs, 1)),
        "cycling_rate": float((phase_counts.get("S", 0) + phase_counts.get("G2M", 0)) / max(adata.n_obs, 1)),
        "elapsed_seconds": round(time.time() - start, 3),
    }
    summary.update(doublet_stats)
    summary.update(cycle_stats)
    print(
        f"DONE {study_sample} cells={adata.n_obs} doublet_rate={summary['predicted_doublet_rate']:.4f} "
        f"cycling_rate={summary['cycling_rate']:.4f} elapsed={summary['elapsed_seconds']}",
        flush=True,
    )
    return diag, summary


def sort_cluster_frame(df: pd.DataFrame, cluster_key: str) -> pd.DataFrame:
    values = df[cluster_key].astype(str)
    if values.str.fullmatch(r"\d+").all():
        return df.assign(_cluster_sort=values.astype(int)).sort_values("_cluster_sort").drop(columns="_cluster_sort")
    return df.sort_values(cluster_key)


def summarize_clusters(obs: pd.DataFrame, cluster_key: str, args: argparse.Namespace) -> pd.DataFrame:
    global_doublet_rate = float(obs["predicted_doublet"].mean())
    rows = []
    for cluster, sub in obs.groupby(cluster_key, observed=True):
        phases = sub["cell_cycle_phase"].astype(str)
        n = int(sub.shape[0])
        doublet_rate = float(sub["predicted_doublet"].mean())
        s_rate = float((phases == "S").mean())
        g2m_rate = float((phases == "G2M").mean())
        cycling_rate = s_rate + g2m_rate
        dominant_phase = str(phases.value_counts().idxmax()) if n else "unknown"
        rows.append(
            {
                cluster_key: str(cluster),
                "n_cells": n,
                "doublet_rate": doublet_rate,
                "predicted_doublets": int(sub["predicted_doublet"].sum()),
                "doublet_score_mean": float(sub["doublet_score"].mean()),
                "doublet_score_median": float(sub["doublet_score"].median()),
                "doublet_score_p95": float(sub["doublet_score"].quantile(0.95)),
                "cycling_rate": float(cycling_rate),
                "s_phase_rate": s_rate,
                "g2m_phase_rate": g2m_rate,
                "dominant_phase": dominant_phase,
                "mean_S_score": float(sub["cell_cycle_S_score"].mean()),
                "mean_G2M_score": float(sub["cell_cycle_G2M_score"].mean()),
            }
        )
    out = pd.DataFrame(rows)
    doublet_cutoff = max(args.doublet_cluster_rate, global_doublet_rate * args.doublet_cluster_fold)
    out["flag_doublet_cluster"] = (out["n_cells"] >= args.min_cluster_cells) & (out["doublet_rate"] >= doublet_cutoff)
    out["flag_cycling_cluster"] = (out["n_cells"] >= args.min_cluster_cells) & (
        out["cycling_rate"] >= args.cycling_cluster_rate
    )
    out["doublet_cluster_cutoff"] = doublet_cutoff
    out["cycling_cluster_cutoff"] = args.cycling_cluster_rate
    return sort_cluster_frame(out, cluster_key)


def plot_umap(
    obs: pd.DataFrame,
    umap: np.ndarray,
    out_dir: Path,
    max_plot_cells: int,
    seed: int,
) -> list[str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    idx = np.arange(obs.shape[0])
    if obs.shape[0] > max_plot_cells:
        idx = np.sort(rng.choice(idx, size=max_plot_cells, replace=False))
    obs_plot = obs.iloc[idx]
    xy = umap[idx]
    paths = []

    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    sca = ax.scatter(xy[:, 0], xy[:, 1], c=obs_plot["doublet_score"], s=0.8, cmap="viridis", linewidths=0)
    ax.set_title("Doublet score on scVI UMAP")
    ax.set_xticks([])
    ax.set_yticks([])
    fig.colorbar(sca, ax=ax, label="doublet_score")
    path = out_dir / "umap_doublet_score.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    paths.append(str(path.resolve()))

    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    mask = obs_plot["predicted_doublet"].to_numpy(dtype=bool)
    ax.scatter(xy[~mask, 0], xy[~mask, 1], s=0.5, c="#d0d0d0", linewidths=0, alpha=0.45)
    ax.scatter(xy[mask, 0], xy[mask, 1], s=1.2, c="#d62728", linewidths=0, alpha=0.8)
    ax.set_title("Predicted doublets on scVI UMAP")
    ax.set_xticks([])
    ax.set_yticks([])
    path = out_dir / "umap_predicted_doublets.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    paths.append(str(path.resolve()))

    phase_colors = {"G1": "#7f7f7f", "S": "#1f77b4", "G2M": "#ff7f0e", "unknown": "#cccccc"}
    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    for phase, color in phase_colors.items():
        mask = obs_plot["cell_cycle_phase"].astype(str).to_numpy() == phase
        if mask.any():
            ax.scatter(xy[mask, 0], xy[mask, 1], s=0.8, c=color, label=phase, linewidths=0, alpha=0.65)
    ax.legend(frameon=False, markerscale=5)
    ax.set_title("Cell-cycle phase on scVI UMAP")
    ax.set_xticks([])
    ax.set_yticks([])
    path = out_dir / "umap_cell_cycle_phase.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    paths.append(str(path.resolve()))

    return paths


def plot_cluster_rates(cluster_summary: pd.DataFrame, cluster_key: str, out_dir: Path) -> str:
    out_dir.mkdir(parents=True, exist_ok=True)
    labels = cluster_summary[cluster_key].astype(str).to_numpy()
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(14, 5), constrained_layout=True)
    ax.bar(x - 0.2, cluster_summary["doublet_rate"], width=0.4, label="doublet_rate")
    ax.bar(x + 0.2, cluster_summary["cycling_rate"], width=0.4, label="cycling_rate")
    ax.axhline(cluster_summary["doublet_cluster_cutoff"].iloc[0], color="#d62728", linestyle="--", linewidth=1)
    ax.axhline(cluster_summary["cycling_cluster_cutoff"].iloc[0], color="#ff7f0e", linestyle="--", linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=90)
    ax.set_ylim(0, min(1.0, max(0.1, cluster_summary[["doublet_rate", "cycling_rate"]].max().max() * 1.15)))
    ax.set_xlabel(cluster_key)
    ax.set_ylabel("fraction of cells")
    ax.legend(frameon=False)
    path = out_dir / "cluster_doublet_cycling_rates.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return str(path.resolve())


def main() -> int:
    args = parse_args()
    start = time.time()
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    args.figures_dir.mkdir(parents=True, exist_ok=True)

    backed = ad.read_h5ad(args.input, backed="r")
    integrated_obs = backed.obs.copy()
    if args.cluster_key not in integrated_obs.columns:
        raise KeyError(f"{args.cluster_key!r} is not present in input obs")
    samples = sample_table(integrated_obs, args.cluster_key)
    if args.only_study_sample:
        keep = samples["study_sample"].astype(str).isin(set(args.only_study_sample))
        samples = samples.loc[keep].copy()
    if args.max_samples is not None:
        samples = samples.head(args.max_samples).copy()
    if samples.empty:
        raise ValueError("No samples selected for diagnostics")
    backed.file.close()
    print(f"INPUT cells={integrated_obs.shape[0]} samples={samples.shape[0]}", flush=True)

    diagnostics = []
    sample_summaries = []
    for i, (_, row) in enumerate(samples.iterrows(), start=1):
        diag, summary = process_one_sample(row, integrated_obs, args, i, samples.shape[0])
        diagnostics.append(diag)
        sample_summaries.append(summary)

    per_cell = pd.concat(diagnostics, axis=0)
    per_cell = per_cell.reindex(integrated_obs.index)
    missing_diagnostics = int(per_cell["doublet_score"].isna().sum())
    if missing_diagnostics:
        print(f"WARN missing_diagnostics={missing_diagnostics}", flush=True)

    obs_diag = integrated_obs.join(
        per_cell[
            [
                "doublet_score",
                "predicted_doublet",
                "doublet_threshold",
                "doublet_method",
                "cell_cycle_S_score",
                "cell_cycle_G2M_score",
                "cell_cycle_phase",
            ]
        ]
    )
    obs_diag["predicted_doublet"] = obs_diag["predicted_doublet"].where(
        obs_diag["predicted_doublet"].notna(), False
    ).astype(bool)
    cluster_summary = summarize_clusters(obs_diag, args.cluster_key, args)
    flagged = cluster_summary[
        cluster_summary["flag_doublet_cluster"] | cluster_summary["flag_cycling_cluster"]
    ].copy()
    flagged = flagged.sort_values(["flag_doublet_cluster", "flag_cycling_cluster", "doublet_rate", "cycling_rate"], ascending=False)

    per_cell_path = args.metadata_dir / "scvi_doublet_cell_cycle_by_cell.tsv.gz"
    sample_summary_path = args.metadata_dir / "scvi_doublet_cell_cycle_by_sample.tsv"
    cluster_summary_path = args.metadata_dir / "scvi_doublet_cell_cycle_by_leiden.tsv"
    flagged_path = args.metadata_dir / "scvi_flagged_doublet_cycling_clusters.tsv"
    per_cell.to_csv(per_cell_path, sep="\t", compression="gzip")
    pd.DataFrame(sample_summaries).to_csv(sample_summary_path, sep="\t", index=False)
    cluster_summary.to_csv(cluster_summary_path, sep="\t", index=False)
    flagged.to_csv(flagged_path, sep="\t", index=False)

    figure_paths: list[str] = []
    if not args.skip_annotated_h5ad:
        print(f"WRITE annotated h5ad {args.output}", flush=True)
        adata = ad.read_h5ad(args.input)
        for col in [
            "doublet_score",
            "predicted_doublet",
            "doublet_threshold",
            "doublet_method",
            "cell_cycle_S_score",
            "cell_cycle_G2M_score",
            "cell_cycle_phase",
        ]:
            adata.obs[col] = obs_diag[col].to_numpy()
        args.output.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(args.output, compression="gzip")
        if "X_umap" in adata.obsm:
            figure_paths.extend(plot_umap(adata.obs, np.asarray(adata.obsm["X_umap"]), args.figures_dir, args.max_plot_cells, args.seed))
    else:
        backed = ad.read_h5ad(args.input, backed="r")
        if "X_umap" in backed.obsm:
            figure_paths.extend(plot_umap(obs_diag, np.asarray(backed.obsm["X_umap"]), args.figures_dir, args.max_plot_cells, args.seed))
        backed.file.close()
    figure_paths.append(plot_cluster_rates(cluster_summary, args.cluster_key, args.figures_dir))

    report = {
        "input": str(args.input.resolve()),
        "output": None if args.skip_annotated_h5ad else str(args.output.resolve()),
        "n_cells": int(obs_diag.shape[0]),
        "n_samples": int(samples.shape[0]),
        "cluster_key": args.cluster_key,
        "expected_doublet_rate": float(args.expected_doublet_rate),
        "global_predicted_doublet_rate": float(obs_diag["predicted_doublet"].mean()),
        "global_cycling_rate": float(obs_diag["cell_cycle_phase"].astype(str).isin(["S", "G2M"]).mean()),
        "flagged_doublet_clusters": cluster_summary.loc[
            cluster_summary["flag_doublet_cluster"], args.cluster_key
        ].astype(str).to_list(),
        "flagged_cycling_clusters": cluster_summary.loc[
            cluster_summary["flag_cycling_cluster"], args.cluster_key
        ].astype(str).to_list(),
        "missing_diagnostics": missing_diagnostics,
        "per_cell_path": str(per_cell_path.resolve()),
        "sample_summary_path": str(sample_summary_path.resolve()),
        "cluster_summary_path": str(cluster_summary_path.resolve()),
        "flagged_clusters_path": str(flagged_path.resolve()),
        "figures": figure_paths,
        "doublet_method": "simulated_doublet_knn",
        "scanpy_version": version("scanpy"),
        "elapsed_seconds": round(time.time() - start, 3),
    }
    report_path = args.metadata_dir / "scvi_doublet_cell_cycle_report.json"
    with report_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)
    print(f"WROTE {report_path}", flush=True)
    print(json.dumps(report, indent=2, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
