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
from ctxcore.genesig import GeneSignature
from pyscenic.aucell import aucell
from scipy import sparse
from scipy.stats import pearsonr, spearmanr


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 6.3: pySCENIC-style regulon activity with AUCell.")
    parser.add_argument(
        "--input-h5ad",
        type=Path,
        default=ROOT / "data/processed/driver/driver_cellrank_main_strict.module6_2.h5ad",
    )
    parser.add_argument(
        "--cellrank-drivers",
        type=Path,
        default=ROOT / "metadata/driver/driver_module6_2_cellrank_lineage_drivers.tsv.gz",
    )
    parser.add_argument(
        "--tf-list",
        type=Path,
        default=ROOT / "metadata/driver/scenic_resources/allTFs_hg38.txt",
    )
    parser.add_argument("--processed-dir", type=Path, default=ROOT / "data/processed/driver")
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/driver")
    parser.add_argument("--figures-dir", type=Path, default=ROOT / "figures/driver")
    parser.add_argument("--output-name", default="driver_pyscenic_regulon_activity.module6_3.h5ad")
    parser.add_argument("--fate-key", default="cellrank_fate_prob_cnv_supported_malignant")
    parser.add_argument("--time-key", default="driver_main_strict__pseudotime_median")
    parser.add_argument("--phase-key", default="driver_main_strict__pseudotime_phase")
    parser.add_argument("--min-target-corr", type=float, default=0.05)
    parser.add_argument("--min-targets", type=int, default=10)
    parser.add_argument("--max-targets", type=int, default=50)
    parser.add_argument("--max-regulons", type=int, default=150)
    parser.add_argument("--aucell-threshold", type=float, default=0.05)
    parser.add_argument("--num-workers", type=int, default=1)
    parser.add_argument("--top-n", type=int, default=30)
    parser.add_argument("--seed", type=int, default=20260606)
    parser.add_argument("--compression", default="gzip")
    return parser.parse_args()


def configure_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica"],
            "font.size": 8,
            "axes.labelsize": 8,
            "axes.titlesize": 9,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def write_dataframe(path: Path, df: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = "gzip" if path.suffix == ".gz" else None
    df.to_csv(path, sep="\t", index=False, compression=compression)


def bh_qvalues(pvalues: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvalues, errors="coerce")
    q = pd.Series(np.nan, index=p.index, dtype=float)
    valid = p.dropna()
    if valid.empty:
        return q
    order = valid.sort_values().index
    ranked = valid.loc[order].to_numpy(dtype=float)
    m = float(len(ranked))
    adjusted = ranked * m / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    q.loc[order] = np.clip(adjusted, 0.0, 1.0)
    return q


def dense_expression(adata: ad.AnnData) -> np.ndarray:
    x = adata.X
    if sparse.issparse(x):
        x = x.toarray()
    return np.asarray(x, dtype=np.float32)


def zscore_columns(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean = matrix.mean(axis=0, keepdims=True)
    std = matrix.std(axis=0, ddof=1, keepdims=True)
    std = np.where(std <= 0, np.nan, std)
    z = (matrix - mean) / std
    z = np.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)
    return z, mean.ravel(), std.ravel()


def read_tf_list(path: Path, genes: pd.Index) -> list[str]:
    if not path.exists():
        raise FileNotFoundError(f"Missing TF list: {path}")
    gene_set = set(genes.astype(str))
    tfs = []
    for line in path.read_text(encoding="utf-8").splitlines():
        tf = line.strip()
        if tf and tf in gene_set and tf not in tfs:
            tfs.append(tf)
    if not tfs:
        raise ValueError("No TFs from the TF list are present in the expression matrix.")
    return tfs


def read_gene_fate_corr(path: Path, lineage: str = "cnv_supported_malignant") -> pd.Series:
    drivers = pd.read_csv(path, sep="\t")
    if "lineage" in drivers.columns and "corr" in drivers.columns:
        sub = drivers.loc[drivers["lineage"].astype(str).eq(lineage), ["gene", "corr"]].copy()
        return pd.Series(pd.to_numeric(sub["corr"], errors="coerce").to_numpy(), index=sub["gene"].astype(str)).dropna()
    corr_col = f"{lineage}_corr"
    if corr_col not in drivers.columns:
        return pd.Series(dtype=float)
    genes = drivers["gene"].astype(str) if "gene" in drivers.columns else drivers.index.astype(str)
    return pd.Series(pd.to_numeric(drivers[corr_col], errors="coerce").to_numpy(), index=genes).dropna()


def infer_coexpression_regulons(
    x_z: np.ndarray,
    genes: pd.Index,
    tf_names: list[str],
    gene_fate_corr: pd.Series,
    min_corr: float,
    min_targets: int,
    max_targets: int,
    max_regulons: int,
) -> tuple[list[GeneSignature], pd.DataFrame, pd.DataFrame]:
    gene_list = genes.astype(str).tolist()
    gene_index = {gene: idx for idx, gene in enumerate(gene_list)}
    tf_indices = [gene_index[tf] for tf in tf_names]
    corr = (x_z[:, tf_indices].T @ x_z) / max(1, x_z.shape[0] - 1)
    fate_corr = gene_fate_corr.reindex(gene_list).fillna(0.0)

    regulon_rows = []
    edge_rows = []
    signatures: list[GeneSignature] = []
    for row_idx, tf in enumerate(tf_names):
        values = corr[row_idx, :].astype(float)
        values[tf_indices[row_idx]] = -np.inf
        candidate_idx = np.flatnonzero(values >= min_corr)
        if candidate_idx.size == 0:
            continue
        ranked_idx = candidate_idx[np.argsort(values[candidate_idx])[::-1]][:max_targets]
        targets = [gene_list[idx] for idx in ranked_idx if np.isfinite(values[idx])]
        if len(targets) < min_targets:
            continue
        weights = {gene: float(max(values[gene_index[gene]], 0.0)) for gene in targets}
        regulon_name = f"{tf}(+)"
        target_fate = fate_corr.reindex(targets).dropna()
        tf_fate = float(fate_corr.get(tf, np.nan))
        mean_target_fate = float(target_fate.mean()) if not target_fate.empty else np.nan
        positive_driver_targets = int((target_fate > 0).sum()) if not target_fate.empty else 0
        priority = np.nan_to_num(tf_fate, nan=0.0) + np.nan_to_num(mean_target_fate, nan=0.0)
        regulon_rows.append(
            {
                "regulon": regulon_name,
                "tf": tf,
                "n_targets": len(targets),
                "mean_tf_target_corr": float(np.mean([weights[g] for g in targets])),
                "max_tf_target_corr": float(np.max([weights[g] for g in targets])),
                "tf_cellrank_cnv_fate_corr": tf_fate,
                "mean_target_cellrank_cnv_fate_corr": mean_target_fate,
                "positive_driver_target_fraction": positive_driver_targets / float(len(targets)),
                "priority_score": float(priority),
                "motif_pruned": False,
                "regulon_source": "tf_target_coexpression_no_motif_pruning",
            }
        )
        for rank, gene in enumerate(targets, start=1):
            edge_rows.append(
                {
                    "regulon": regulon_name,
                    "tf": tf,
                    "target": gene,
                    "target_rank": rank,
                    "tf_target_corr": weights[gene],
                    "target_cellrank_cnv_fate_corr": float(fate_corr.get(gene, np.nan)),
                }
            )
        signatures.append(GeneSignature(regulon_name, weights))

    regulons = pd.DataFrame(regulon_rows)
    edges = pd.DataFrame(edge_rows)
    if regulons.empty:
        raise ValueError("No co-expression regulons passed filtering.")
    regulons = regulons.sort_values(["priority_score", "mean_tf_target_corr", "n_targets"], ascending=[False, False, False]).reset_index(drop=True)
    keep_names = set(regulons.head(max_regulons)["regulon"])
    regulons = regulons.loc[regulons["regulon"].isin(keep_names)].copy()
    edges = edges.loc[edges["regulon"].isin(keep_names)].copy()
    signatures = [signature for signature in signatures if signature.name in keep_names]
    return signatures, regulons.reset_index(drop=True), edges.reset_index(drop=True)


def correlate_regulon_activity(auc: pd.DataFrame, cells: pd.DataFrame, fate_key: str, time_key: str) -> pd.DataFrame:
    rows = []
    fate = pd.to_numeric(cells[fate_key], errors="coerce")
    pseudotime = pd.to_numeric(cells[time_key], errors="coerce") if time_key in cells.columns else pd.Series(np.nan, index=cells.index)
    for regulon in auc.columns.astype(str):
        values = pd.to_numeric(auc[regulon], errors="coerce")
        row = {"regulon": regulon}
        for label, target in [("cnv_fate", fate), ("pseudotime", pseudotime)]:
            mask = values.notna() & target.notna() & np.isfinite(values) & np.isfinite(target)
            if mask.sum() < 3 or values.loc[mask].nunique() < 2 or target.loc[mask].nunique() < 2:
                row[f"{label}_pearson_r"] = np.nan
                row[f"{label}_pearson_p"] = np.nan
                row[f"{label}_spearman_rho"] = np.nan
                row[f"{label}_spearman_p"] = np.nan
                continue
            pr, pp = pearsonr(values.loc[mask], target.loc[mask])
            sr, sp = spearmanr(values.loc[mask], target.loc[mask])
            row[f"{label}_pearson_r"] = float(pr)
            row[f"{label}_pearson_p"] = float(pp)
            row[f"{label}_spearman_rho"] = float(sr)
            row[f"{label}_spearman_p"] = float(sp)
        rows.append(row)
    out = pd.DataFrame(rows)
    out["cnv_fate_pearson_q"] = bh_qvalues(out["cnv_fate_pearson_p"])
    out["cnv_fate_spearman_q"] = bh_qvalues(out["cnv_fate_spearman_p"])
    return out.sort_values(["cnv_fate_pearson_r", "cnv_fate_spearman_rho"], ascending=[False, False]).reset_index(drop=True)


def summarize_activity_by_group(auc: pd.DataFrame, cells: pd.DataFrame, top_regulons: list[str], group_cols: list[str]) -> pd.DataFrame:
    rows = []
    for group_col in group_cols:
        if group_col not in cells.columns:
            continue
        for group, idx in cells.groupby(group_col, observed=True, sort=True).groups.items():
            for regulon in top_regulons:
                values = pd.to_numeric(auc.loc[idx, regulon], errors="coerce").dropna()
                rows.append(
                    {
                        "group_type": group_col,
                        "group": str(group),
                        "regulon": regulon,
                        "n_cells": int(len(idx)),
                        "mean_auc": float(values.mean()) if not values.empty else np.nan,
                        "median_auc": float(values.median()) if not values.empty else np.nan,
                    }
                )
    return pd.DataFrame(rows)


def plot_top_regulons(corr: pd.DataFrame, top_n: int, path_base: Path) -> list[str]:
    sub = corr.sort_values("cnv_fate_pearson_r", ascending=False).head(top_n).iloc[::-1]
    if sub.empty:
        return []
    fig, ax = plt.subplots(figsize=(5.4, max(3.0, 0.22 * sub.shape[0])))
    ax.barh(sub["regulon"], sub["cnv_fate_pearson_r"], color="#0072B2", height=0.75)
    ax.set_xlabel("Pearson r with CNV fate probability")
    ax.set_ylabel("")
    ax.set_title(f"Top {sub.shape[0]} AUCell regulons")
    fig.tight_layout()
    outputs = []
    for suffix in [".png", ".pdf"]:
        out = path_base.with_suffix(suffix)
        fig.savefig(out, bbox_inches="tight")
        outputs.append(str(out))
    plt.close(fig)
    return outputs


def plot_regulon_heatmap(auc: pd.DataFrame, cells: pd.DataFrame, regulons: list[str], phase_key: str, path_base: Path) -> list[str]:
    if not regulons:
        return []
    phase_order = ["early", "middle", "late"]
    rows = []
    for phase in phase_order:
        idx = cells.index[cells[phase_key].astype(str).eq(phase)] if phase_key in cells.columns else []
        if len(idx) == 0:
            continue
        rows.append(pd.to_numeric(auc.loc[idx, regulons].mean(axis=0), errors="coerce").rename(phase))
    if not rows:
        return []
    matrix = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(max(5.5, 0.28 * len(regulons)), 2.4))
    im = ax.imshow(matrix.to_numpy(dtype=float), aspect="auto", cmap="magma")
    ax.set_yticks(range(matrix.shape[0]))
    ax.set_yticklabels(matrix.index)
    ax.set_xticks(range(len(regulons)))
    ax.set_xticklabels(regulons, rotation=90)
    ax.set_title("Mean regulon AUC by pseudotime phase")
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Mean AUC")
    fig.tight_layout()
    outputs = []
    for suffix in [".png", ".pdf"]:
        out = path_base.with_suffix(suffix)
        fig.savefig(out, bbox_inches="tight")
        outputs.append(str(out))
    plt.close(fig)
    return outputs


def plot_umap_regulon(adata: ad.AnnData, auc: pd.DataFrame, regulon: str, path_base: Path) -> list[str]:
    if "X_umap_global" in adata.obsm:
        xy = np.asarray(adata.obsm["X_umap_global"])
    elif "X_umap" in adata.obsm:
        xy = np.asarray(adata.obsm["X_umap"])
    else:
        return []
    values = pd.to_numeric(auc[regulon], errors="coerce").to_numpy(dtype=float)
    order = np.argsort(np.nan_to_num(values, nan=-1.0))
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    sca = ax.scatter(xy[order, 0], xy[order, 1], c=values[order], s=5, cmap="viridis", linewidths=0, rasterized=True)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(f"{regulon} AUCell activity")
    cbar = fig.colorbar(sca, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("AUC")
    fig.tight_layout()
    outputs = []
    for suffix in [".png", ".pdf"]:
        out = path_base.with_suffix(suffix)
        fig.savefig(out, bbox_inches="tight")
        outputs.append(str(out))
    plt.close(fig)
    return outputs


def main() -> None:
    start = time.time()
    args = parse_args()
    configure_plot_style()
    args.processed_dir.mkdir(parents=True, exist_ok=True)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    args.figures_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.input_h5ad)
    if args.fate_key not in adata.obs.columns:
        raise KeyError(f"Missing fate key in adata.obs: {args.fate_key}")
    x = dense_expression(adata)
    genes = pd.Index(adata.var_names.astype(str))
    tf_names = read_tf_list(args.tf_list, genes)
    gene_fate_corr = read_gene_fate_corr(args.cellrank_drivers)
    x_z, _, _ = zscore_columns(x)

    signatures, regulons, edges = infer_coexpression_regulons(
        x_z=x_z,
        genes=genes,
        tf_names=tf_names,
        gene_fate_corr=gene_fate_corr,
        min_corr=args.min_target_corr,
        min_targets=args.min_targets,
        max_targets=args.max_targets,
        max_regulons=args.max_regulons,
    )

    expr_df = pd.DataFrame(x, index=adata.obs_names.astype(str), columns=genes)
    auc = aucell(
        expr_df,
        signatures,
        auc_threshold=args.aucell_threshold,
        normalize=True,
        seed=args.seed,
        num_workers=args.num_workers,
    )
    auc = auc.reindex(expr_df.index)
    auc.columns = auc.columns.astype(str)

    cells = adata.obs.copy()
    cells.index = adata.obs_names.astype(str)
    corr = correlate_regulon_activity(auc, cells, args.fate_key, args.time_key)
    regulons = regulons.merge(corr, on="regulon", how="left")
    regulons = regulons.sort_values(["cnv_fate_pearson_r", "priority_score"], ascending=[False, False]).reset_index(drop=True)

    top_regulons = regulons.head(args.top_n)["regulon"].tolist()
    group_summary = summarize_activity_by_group(
        auc,
        cells,
        top_regulons=top_regulons,
        group_cols=[args.phase_key, "cell_disease_stage", "trajectory_root_end_role", "dataset"],
    )

    adata.obsm["module6_3_regulon_auc"] = auc.to_numpy(dtype=np.float32)
    adata.uns["module6_3_regulon_auc_names"] = auc.columns.tolist()
    for regulon in top_regulons[:10]:
        safe = regulon.replace("(", "_").replace(")", "").replace("+", "plus").replace("-", "minus")
        adata.obs[f"module6_3_auc_{safe}"] = pd.to_numeric(auc[regulon], errors="coerce").to_numpy(dtype=float)
    adata.uns["module6_3_pyscenic_regulon_activity"] = {
        "module": "6.3",
        "method": "TF-list constrained co-expression regulons scored with pySCENIC AUCell",
        "tf_list": str(args.tf_list),
        "n_tfs_in_matrix": int(len(tf_names)),
        "n_regulons": int(len(signatures)),
        "motif_pruning_status": "not_run_missing_cistarget_motif_ranking_database",
        "grnboost2_status": "not_run_windows_dask_timeout_in_smoke_test",
        "auc_threshold": float(args.aucell_threshold),
    }

    output_h5ad = args.processed_dir / args.output_name
    adata.write_h5ad(output_h5ad, compression=args.compression)

    auc_cells = auc.copy()
    auc_cells.insert(0, "cell_id", auc_cells.index.astype(str))
    auc_path = args.metadata_dir / "driver_module6_3_pyscenic_regulon_auc.tsv.gz"
    regulons_path = args.metadata_dir / "driver_module6_3_pyscenic_regulons.tsv"
    edges_path = args.metadata_dir / "driver_module6_3_pyscenic_regulon_edges.tsv.gz"
    top_path = args.metadata_dir / "driver_module6_3_top_cnv_fate_regulons.tsv"
    summary_path = args.metadata_dir / "driver_module6_3_regulon_activity_group_summary.tsv"

    write_dataframe(auc_path, auc_cells.reset_index(drop=True))
    write_dataframe(regulons_path, regulons)
    write_dataframe(edges_path, edges)
    write_dataframe(top_path, regulons.head(args.top_n))
    write_dataframe(summary_path, group_summary)

    figure_outputs: list[str] = []
    figure_outputs += plot_top_regulons(regulons, args.top_n, args.figures_dir / "driver_module6_3_top_cnv_fate_regulons")
    figure_outputs += plot_regulon_heatmap(
        auc,
        cells,
        regulons=top_regulons[:25],
        phase_key=args.phase_key,
        path_base=args.figures_dir / "driver_module6_3_top_regulon_phase_heatmap",
    )
    if top_regulons:
        figure_outputs += plot_umap_regulon(
            adata,
            auc,
            regulon=top_regulons[0],
            path_base=args.figures_dir / f"driver_module6_3_top_regulon_umap__{top_regulons[0].replace('(+)', '')}",
        )

    report = {
        "module": "6.3",
        "method": "SCENIC/pySCENIC regulon activity with AUCell",
        "input_h5ad": str(args.input_h5ad),
        "output_h5ad": str(output_h5ad),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_tfs_in_matrix": int(len(tf_names)),
        "n_regulons_scored": int(auc.shape[1]),
        "n_regulon_edges": int(edges.shape[0]),
        "motif_pruning_status": "not_run_missing_cistarget_motif_ranking_database",
        "grnboost2_status": "not_run_windows_dask_timeout_in_smoke_test",
        "tf_list_source": "https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt",
        "top_cnv_fate_regulons": regulons.head(15).to_dict(orient="records"),
        "outputs": {
            "h5ad": str(output_h5ad),
            "auc": str(auc_path),
            "regulons": str(regulons_path),
            "edges": str(edges_path),
            "top_regulons": str(top_path),
            "group_summary": str(summary_path),
            "figures": figure_outputs,
        },
        "package_versions": {
            "pyscenic": version("pyscenic"),
            "ctxcore": version("ctxcore"),
            "anndata": version("anndata"),
            "pandas": version("pandas"),
            "numpy": version("numpy"),
        },
        "elapsed_seconds": round(time.time() - start, 3),
    }
    report_path = args.metadata_dir / "driver_module6_3_pyscenic_regulon_activity_report.json"
    report["outputs"]["report"] = str(report_path)
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
