from __future__ import annotations

import argparse
import json
import time
from collections import OrderedDict
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]


STATE_PANELS: OrderedDict[str, list[str]] = OrderedDict(
    [
        (
            "Mature_Hepatocyte",
            ["ALB", "APOA1", "APOA2", "TTR", "HPD", "ASGR1", "CYP3A4", "CYP2E1", "CYP2C9", "HNF4A", "CPS1", "ASS1"],
        ),
        (
            "Stressed_Injured",
            ["HSPA1A", "HSPA1B", "HSP90AA1", "DNAJB1", "FOS", "JUN", "JUNB", "ATF3", "DDIT3", "SAA1", "SAA2", "MT1G", "MT2A", "IER3"],
        ),
        (
            "Regenerative_Progenitor",
            ["KRT19", "EPCAM", "SOX9", "KRT7", "TACSTD2", "CD24", "PROM1", "ANXA4"],
        ),
        ("HCC_Malignant_Associated", ["AFP", "GPC3", "SPP1", "MDK", "IGF2BP3", "MUC1", "CEACAM5"]),
        ("Proliferation", ["MKI67", "TOP2A", "STMN1", "TYMS", "UBE2C", "PCNA", "MCM2", "MCM5", "HMGB2"]),
        ("Cholangiocyte", ["KRT19", "KRT7", "EPCAM", "SOX9", "TACSTD2", "ANXA4", "CLDN4", "MUC1"]),
        ("Immune", ["PTPRC", "LST1", "C1QA", "C1QB", "CD68", "LYZ", "CD3D", "NKG7", "MS4A1", "JCHAIN"]),
        ("Endothelial", ["PECAM1", "VWF", "KDR", "EMCN", "RAMP2", "ESAM"]),
        ("Stromal_HSC_Pericyte", ["COL1A1", "COL1A2", "COL3A1", "DCN", "ACTA2", "RGS5", "PDGFRB", "TAGLN"]),
    ]
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 2 hepatocyte lineage secondary annotation.")
    parser.add_argument(
        "--input",
        type=Path,
        default=ROOT / "data/processed/scvi/scvi_integrated_counts_hvg.celltypist_major.h5ad",
    )
    parser.add_argument(
        "--seed",
        type=Path,
        default=ROOT / "metadata/celltype/scanvi_seed_labels_by_cell.tsv.gz",
    )
    parser.add_argument(
        "--cluster-labels",
        type=Path,
        default=ROOT / "metadata/celltype/manual_major_label_by_cluster.tsv",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=ROOT / "metadata/scvi/scvi_input_manifest.counts.tsv",
    )
    parser.add_argument("--output-dir", type=Path, default=ROOT / "data/processed/hepatocyte")
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/hepatocyte")
    parser.add_argument("--figures-dir", type=Path, default=ROOT / "figures/hepatocyte")
    parser.add_argument("--cluster-key", default="leiden_scvi")
    parser.add_argument("--use-rep", default="X_scVI")
    parser.add_argument("--n-neighbors", type=int, default=30)
    parser.add_argument("--min-dist", type=float, default=0.25)
    parser.add_argument("--resolution", type=float, default=0.8)
    parser.add_argument("--extra-resolution", type=float, action="append", default=[0.4, 1.2])
    parser.add_argument("--seed-value", type=int, default=20260601)
    parser.add_argument("--min-panel-z", type=float, default=0.35)
    parser.add_argument("--min-panel-pct", type=float, default=0.03)
    parser.add_argument("--skip-h5ad", action="store_true")
    parser.add_argument("--reuse-existing", action="store_true")
    return parser.parse_args()


def resolve_path(value: object) -> Path:
    path = Path(str(value))
    if path.exists():
        return path
    text = str(value).replace("\\", "/")
    for anchor in ("data/processed/", "metadata/"):
        idx = text.find(anchor)
        if idx >= 0:
            candidate = ROOT / Path(text[idx:])
            if candidate.exists():
                return candidate
    return path


def read_manifest(path: Path) -> pd.DataFrame:
    manifest = pd.read_csv(path, sep="\t")
    manifest = manifest.loc[manifest["include_in_scvi"].astype(str).str.lower().eq("true")].copy()
    manifest["study_sample"] = manifest["dataset"].astype(str) + "__" + manifest["label"].astype(str)
    return manifest


def original_ids(global_index: pd.Index, study_sample: str) -> pd.Index:
    prefix = f"{study_sample}__"
    return pd.Index([idx[len(prefix) :] if str(idx).startswith(prefix) else idx for idx in global_index.astype(str)])


def as_csr(x) -> sparse.csr_matrix:
    if sparse.issparse(x):
        return x.tocsr()
    return sparse.csr_matrix(x)


def sort_key(value: object) -> tuple[int, str]:
    text = str(value)
    return (0, f"{int(text):05d}") if text.isdigit() else (1, text)


def res_key(resolution: float) -> str:
    return str(resolution).replace(".", "_")


def top_counts(values: pd.Series, n: int = 5) -> str:
    counts = values.astype(str).value_counts(dropna=False).head(n)
    total = int(counts.sum())
    if total == 0:
        return ""
    return "; ".join([f"{idx}:{count}({count / total:.1%})" for idx, count in counts.items()])


def expression_support(group: str, mean_log1p_cpm: float, mean_pct_expr: float, panel_score_z: float, args: argparse.Namespace) -> bool:
    cutoffs = {
        "Mature_Hepatocyte": (4.5, 0.20),
        "Stressed_Injured": (3.5, 0.15),
        "Regenerative_Progenitor": (2.5, 0.03),
        "HCC_Malignant_Associated": (3.0, 0.05),
        "Proliferation": (3.0, 0.05),
        "Cholangiocyte": (2.5, 0.03),
        "Immune": (2.0, 0.04),
        "Endothelial": (2.0, 0.04),
        "Stromal_HSC_Pericyte": (2.0, 0.04),
    }
    min_log, min_pct = cutoffs.get(group, (2.0, 0.03))
    z_support = panel_score_z >= args.min_panel_z and mean_pct_expr >= min_pct
    expr_support = mean_log1p_cpm >= min_log and mean_pct_expr >= min_pct
    return bool(z_support or expr_support)


def assign_hepatocyte_state(panel_rows: pd.DataFrame, args: argparse.Namespace) -> dict[str, object]:
    rows = {row["marker_group"]: row for _, row in panel_rows.iterrows()}

    def z(name: str) -> float:
        return float(rows[name]["panel_score_z"])

    def pct(name: str) -> float:
        return float(rows[name]["mean_pct_expr"])

    def logcpm(name: str) -> float:
        return float(rows[name]["mean_log1p_cpm"])

    mature = logcpm("Mature_Hepatocyte") >= 4.5 and pct("Mature_Hepatocyte") >= 0.20
    stress = z("Stressed_Injured") >= 0.50 and pct("Stressed_Injured") >= 0.15
    progenitor = (z("Regenerative_Progenitor") >= 0.50 and pct("Regenerative_Progenitor") >= 0.03) or (
        z("Cholangiocyte") >= 0.50 and pct("Cholangiocyte") >= 0.03
    )
    malignant = (z("HCC_Malignant_Associated") >= 0.80 and pct("HCC_Malignant_Associated") >= 0.05) or (
        logcpm("HCC_Malignant_Associated") >= 3.5 and pct("HCC_Malignant_Associated") >= 0.10
    )
    proliferating = z("Proliferation") >= 0.80 and pct("Proliferation") >= 0.05
    immune = z("Immune") >= 0.80 and pct("Immune") >= 0.05
    endothelial = z("Endothelial") >= 0.80 and pct("Endothelial") >= 0.05
    stromal = z("Stromal_HSC_Pericyte") >= 0.80 and pct("Stromal_HSC_Pericyte") >= 0.05
    contaminant_scores = {
        "likely_immune_contaminant": rows["Immune"]["panel_score_z"],
        "likely_endothelial_contaminant": rows["Endothelial"]["panel_score_z"],
        "likely_stromal_hsc_pericyte_contaminant": rows["Stromal_HSC_Pericyte"]["panel_score_z"],
    }
    top_contaminant = max(contaminant_scores, key=contaminant_scores.get)
    has_strong_contaminant = (immune or endothelial or stromal) and not (mature or progenitor or malignant)

    if has_strong_contaminant:
        state = top_contaminant
        confidence = "conflict"
        seed = "Unknown"
    elif malignant and (mature or progenitor or proliferating):
        state = "malignant_hepatocyte_candidate_needs_cnv"
        confidence = "needs_cnv"
        seed = "Unknown"
    elif proliferating and (mature or progenitor):
        state = "proliferating_hepatocyte_candidate"
        confidence = "needs_cnv"
        seed = "Unknown"
    elif progenitor and mature:
        state = "regenerative_progenitor_like_hepatocyte"
        confidence = "high_conf"
        seed = state
    elif progenitor:
        state = "progenitor_cholangiocyte_like_epithelial"
        confidence = "high_conf"
        seed = state
    elif stress and mature:
        state = "stressed_injured_hepatocyte"
        confidence = "high_conf"
        seed = state
    elif mature:
        state = "normal_hepatocyte_like"
        confidence = "high_conf"
        seed = state
    else:
        state = "ambiguous_epithelial_or_mixed"
        confidence = "ambiguous"
        seed = "Unknown"

    ranked = panel_rows.sort_values(["panel_score_z", "mean_log1p_cpm", "mean_pct_expr"], ascending=False)
    top = ranked.iloc[0]
    second = ranked.iloc[1]
    evidence = [
        f"top={top['marker_group']}:{top['panel_score_z']:.3f}",
        f"second={second['marker_group']}:{second['panel_score_z']:.3f}",
        f"mature={mature}",
        f"stress={stress}",
        f"progenitor={progenitor}",
        f"malignant_assoc={malignant}",
        f"proliferating={proliferating}",
        f"immune={immune}",
        f"endothelial={endothelial}",
        f"stromal={stromal}",
    ]
    return {
        "hepatocyte_state_label": state,
        "hepatocyte_state_confidence": confidence,
        "hepatocyte_state_seed_label": seed,
        "state_top_panel": str(top["marker_group"]),
        "state_top_panel_score_z": float(top["panel_score_z"]),
        "state_second_panel": str(second["marker_group"]),
        "state_second_panel_score_z": float(second["panel_score_z"]),
        "state_evidence": "; ".join(evidence),
    }


def aggregate_marker_panels(
    obs: pd.DataFrame,
    group_key: str,
    manifest: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    groups = sorted(obs[group_key].astype(str).unique(), key=sort_key)
    group_to_col = {group: i for i, group in enumerate(groups)}
    marker_genes: list[str] = []
    marker_to_groups: dict[str, list[str]] = {}
    for group, genes in STATE_PANELS.items():
        for gene in genes:
            if gene not in marker_to_groups:
                marker_genes.append(gene)
                marker_to_groups[gene] = []
            marker_to_groups[gene].append(group)
    gene_to_idx = {gene: i for i, gene in enumerate(marker_genes)}
    n_genes = len(marker_genes)
    n_groups = len(groups)
    raw_sums = np.zeros((n_genes, n_groups), dtype=np.float64)
    detected = np.zeros((n_genes, n_groups), dtype=np.float64)
    cell_denoms = np.zeros((n_genes, n_groups), dtype=np.float64)
    count_denoms = np.zeros((n_genes, n_groups), dtype=np.float64)
    sample_presence = np.zeros((n_genes, n_groups), dtype=np.float64)

    for i, row in manifest.reset_index(drop=True).iterrows():
        study_sample = str(row["study_sample"])
        sample_obs = obs.loc[obs["study_sample"].astype(str).eq(study_sample)]
        if sample_obs.empty:
            continue
        path = resolve_path(row["output"])
        print(f"MARKERS {i + 1}/{manifest.shape[0]} {study_sample} cells={sample_obs.shape[0]}", flush=True)
        a = ad.read_h5ad(path, backed="r")
        var_names = pd.Index(a.var_names.astype(str))
        present_genes = [gene for gene in marker_genes if gene in var_names]
        if not present_genes:
            a.file.close()
            continue
        wanted_original = original_ids(sample_obs.index, study_sample)
        present_cells = wanted_original.intersection(a.obs_names.astype(str))
        if len(present_cells) == 0:
            a.file.close()
            continue
        sample_obs = sample_obs.iloc[wanted_original.get_indexer(present_cells)]
        sample_group_cols = sample_obs[group_key].astype(str).map(group_to_col).to_numpy()
        indicator = sparse.csr_matrix(
            (np.ones(len(present_cells), dtype=np.float64), (np.arange(len(present_cells)), sample_group_cols)),
            shape=(len(present_cells), n_groups),
        )
        x = as_csr(a[present_cells, present_genes].X).astype(np.float64)
        raw = (x.T @ indicator).toarray()
        det = ((x > 0).astype(np.float64).T @ indicator).toarray()
        gene_locs = np.array([gene_to_idx[gene] for gene in present_genes], dtype=int)
        raw_sums[gene_locs, :] += raw
        detected[gene_locs, :] += det
        sample_group_cells = np.bincount(sample_group_cols, minlength=n_groups).astype(np.float64)
        sample_group_counts = (
            sample_obs.groupby(sample_obs[group_key].astype(str), observed=True)["qc_total_counts"].sum()
            .reindex(groups)
            .fillna(0.0)
            .to_numpy(dtype=np.float64)
        )
        for loc in gene_locs:
            cell_denoms[loc, :] += sample_group_cells
            count_denoms[loc, :] += sample_group_counts
            sample_presence[loc, :] += (sample_group_cells > 0).astype(np.float64)
        a.file.close()

    cpm = np.divide(raw_sums, count_denoms, out=np.zeros_like(raw_sums), where=count_denoms > 0) * 1e6
    log1p_cpm = np.log1p(cpm)
    pct_expr = np.divide(detected, cell_denoms, out=np.zeros_like(detected), where=cell_denoms > 0)
    gene_z = np.zeros_like(log1p_cpm)
    for i in range(n_genes):
        valid = cell_denoms[i, :] > 0
        values = log1p_cpm[i, valid]
        if values.size > 1 and float(values.std()) > 0:
            gene_z[i, valid] = (values - values.mean()) / values.std()

    marker_rows = []
    for gene in marker_genes:
        loc = gene_to_idx[gene]
        for group in groups:
            col = group_to_col[group]
            for marker_group in marker_to_groups[gene]:
                marker_rows.append(
                    {
                        group_key: group,
                        "marker_group": marker_group,
                        "gene": gene,
                        "raw_sum": raw_sums[loc, col],
                        "detected_cells": detected[loc, col],
                        "n_cells_gene_present": cell_denoms[loc, col],
                        "pct_expr": pct_expr[loc, col],
                        "log1p_cpm": log1p_cpm[loc, col],
                        "gene_z": gene_z[loc, col],
                        "gene_present_sample_count": sample_presence[loc, col],
                    }
                )
    marker_df = pd.DataFrame(marker_rows)

    panel_rows = []
    for group in groups:
        col = group_to_col[group]
        for marker_group, genes in STATE_PANELS.items():
            locs = [gene_to_idx[gene] for gene in genes if gene in gene_to_idx and cell_denoms[gene_to_idx[gene], col] > 0]
            if locs:
                panel_rows.append(
                    {
                        group_key: group,
                        "marker_group": marker_group,
                        "panel_score_z": float(np.mean(gene_z[locs, col])),
                        "mean_log1p_cpm": float(np.mean(log1p_cpm[locs, col])),
                        "mean_pct_expr": float(np.mean(pct_expr[locs, col])),
                        "n_genes_available": int(len(locs)),
                        "genes_available": ";".join([marker_genes[loc] for loc in locs]),
                    }
                )
            else:
                panel_rows.append(
                    {
                        group_key: group,
                        "marker_group": marker_group,
                        "panel_score_z": 0.0,
                        "mean_log1p_cpm": 0.0,
                        "mean_pct_expr": 0.0,
                        "n_genes_available": 0,
                        "genes_available": "",
                    }
                )
    panel_df = pd.DataFrame(panel_rows)
    return marker_df, panel_df


def recompute_panel_scores(marker_df: pd.DataFrame, group_key: str) -> pd.DataFrame:
    groups = sorted(marker_df[group_key].astype(str).unique(), key=sort_key)
    gene_df = marker_df.drop_duplicates([group_key, "gene"]).copy()
    panel_rows = []
    for group in groups:
        sub = gene_df.loc[gene_df[group_key].astype(str).eq(str(group))]
        for marker_group, genes in STATE_PANELS.items():
            panel_gene_df = sub.loc[sub["gene"].isin(genes)]
            if panel_gene_df.empty:
                panel_rows.append(
                    {
                        group_key: group,
                        "marker_group": marker_group,
                        "panel_score_z": 0.0,
                        "mean_log1p_cpm": 0.0,
                        "mean_pct_expr": 0.0,
                        "n_genes_available": 0,
                        "genes_available": "",
                    }
                )
            else:
                panel_rows.append(
                    {
                        group_key: group,
                        "marker_group": marker_group,
                        "panel_score_z": float(panel_gene_df["gene_z"].mean()),
                        "mean_log1p_cpm": float(panel_gene_df["log1p_cpm"].mean()),
                        "mean_pct_expr": float(panel_gene_df["pct_expr"].mean()),
                        "n_genes_available": int(panel_gene_df["gene"].nunique()),
                        "genes_available": ";".join([gene for gene in genes if gene in set(panel_gene_df["gene"])]),
                    }
                )
    return pd.DataFrame(panel_rows)


def plot_panel_heatmap(panel_df: pd.DataFrame, group_key: str, figures_dir: Path) -> None:
    groups = sorted(panel_df[group_key].astype(str).unique(), key=sort_key)
    panels = list(STATE_PANELS.keys())
    matrix = (
        panel_df.pivot(index=group_key, columns="marker_group", values="panel_score_z")
        .reindex(index=groups, columns=panels)
        .fillna(0.0)
    )
    fig, ax = plt.subplots(figsize=(max(9, len(panels) * 0.75), max(8, len(groups) * 0.25)))
    image = ax.imshow(matrix.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax.set_xticks(range(len(panels)))
    ax.set_xticklabels(panels, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels(groups, fontsize=7)
    ax.set_xlabel("State marker panel")
    ax.set_ylabel("Hepatocyte subcluster")
    ax.set_title("Hepatocyte state panel z-score by subcluster")
    cbar = fig.colorbar(image, ax=ax, pad=0.01)
    cbar.set_label("Mean gene z-score")
    fig.tight_layout()
    figures_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figures_dir / "hepatocyte_state_panel_heatmap.png", dpi=300)
    fig.savefig(figures_dir / "hepatocyte_state_panel_heatmap.pdf")
    plt.close(fig)


def plot_marker_dotplot(marker_df: pd.DataFrame, group_key: str, figures_dir: Path) -> None:
    groups = sorted(marker_df[group_key].astype(str).unique(), key=sort_key)
    genes: list[str] = []
    for panel_genes in STATE_PANELS.values():
        for gene in panel_genes:
            if gene not in genes:
                genes.append(gene)
    data = marker_df.loc[marker_df["gene"].isin(genes)].copy()
    data[group_key] = pd.Categorical(data[group_key].astype(str), categories=groups, ordered=True)
    data["gene"] = pd.Categorical(data["gene"], categories=genes, ordered=True)
    data = data.drop_duplicates([group_key, "gene"]).sort_values([group_key, "gene"])
    x = data["gene"].cat.codes.to_numpy()
    y = data[group_key].cat.codes.to_numpy()
    size = np.clip(data["pct_expr"].fillna(0).to_numpy(), 0, 1) * 180 + 2
    color = data["log1p_cpm"].fillna(0).to_numpy()
    fig, ax = plt.subplots(figsize=(max(18, len(genes) * 0.28), max(8, len(groups) * 0.25)))
    scatter = ax.scatter(x, y, s=size, c=color, cmap="viridis", edgecolors="none")
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=90, fontsize=6)
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels(groups, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Marker gene")
    ax.set_ylabel("Hepatocyte subcluster")
    ax.set_title("Hepatocyte lineage state markers")
    cbar = fig.colorbar(scatter, ax=ax, pad=0.01)
    cbar.set_label("log1p(CPM)")
    fig.tight_layout()
    figures_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figures_dir / "hepatocyte_state_marker_dotplot.png", dpi=300)
    fig.savefig(figures_dir / "hepatocyte_state_marker_dotplot.pdf")
    plt.close(fig)


def plot_umaps(adata: ad.AnnData, args: argparse.Namespace) -> None:
    sc.settings.figdir = str(args.figures_dir)
    sc.settings.set_figure_params(dpi=120, facecolor="white", frameon=False, figsize=(7, 6))
    args.figures_dir.mkdir(parents=True, exist_ok=True)
    for color, name in [
        ("leiden_hep", "hepatocyte_lineage_umap_subcluster.png"),
        ("hepatocyte_state_label", "hepatocyte_lineage_umap_state.png"),
        ("manual_major_label_cluster", "hepatocyte_lineage_umap_module1_label.png"),
        ("dataset", "hepatocyte_lineage_umap_dataset.png"),
    ]:
        fig = sc.pl.umap(adata, color=color, show=False, return_fig=True, frameon=False)
        fig.savefig(args.figures_dir / name, dpi=300, bbox_inches="tight")
        plt.close(fig)


def main() -> int:
    args = parse_args()
    start = time.time()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    args.figures_dir.mkdir(parents=True, exist_ok=True)
    marker_path = args.metadata_dir / "hepatocyte_state_marker_scores.tsv"
    panel_path = args.metadata_dir / "hepatocyte_state_panel_scores.tsv"
    summary_path = args.metadata_dir / "hepatocyte_subcluster_summary.tsv"
    cells_path = args.metadata_dir / "hepatocyte_lineage_cells.tsv.gz"
    seed_path = args.metadata_dir / "hepatocyte_state_seed_by_cell.tsv.gz"
    cnv_candidates_path = args.metadata_dir / "hepatocyte_cnv_candidate_cells.tsv.gz"
    cnv_by_sample_path = args.metadata_dir / "hepatocyte_cnv_candidate_by_sample.tsv"
    report_path = args.metadata_dir / "hepatocyte_lineage_module2_report.json"
    h5ad_path = args.output_dir / "hepatocyte_lineage.global_scvi_subcluster.h5ad"

    seed = pd.read_csv(args.seed, sep="\t").set_index("cell_id")
    cluster_labels = pd.read_csv(args.cluster_labels, sep="\t")
    cluster_labels["leiden_scvi"] = cluster_labels["leiden_scvi"].astype(str)
    cluster_label_map = cluster_labels.set_index("leiden_scvi").to_dict(orient="index")

    if args.reuse_existing:
        print("REUSE existing hepatocyte h5ad and marker scores", flush=True)
        adata = sc.read_h5ad(h5ad_path)
        if not marker_path.exists():
            raise FileNotFoundError(marker_path)
        marker_df = pd.read_csv(marker_path, sep="\t")
        panel_df = recompute_panel_scores(marker_df, "leiden_hep")
    else:
        print("READ backed input", args.input, flush=True)
        backed = sc.read_h5ad(args.input, backed="r")
        obs = backed.obs.copy()
        seed = seed.reindex(obs.index.astype(str))
        candidate = seed["hepatocyte_lineage_candidate"].fillna(False).astype(bool)
        if "excluded_doublet_cluster" in obs.columns:
            candidate &= ~obs["excluded_doublet_cluster"].astype(bool).to_numpy()
        print(f"HEP CANDIDATE cells={int(candidate.sum())} of {obs.shape[0]}", flush=True)
        adata = backed[candidate.to_numpy(), :].to_memory()
        backed.file.close()

        aligned_seed = seed.loc[adata.obs_names.astype(str)]
        adata.obs["module1_scanvi_seed_label_major"] = aligned_seed["scanvi_seed_label_major"].fillna("Unknown").astype(str).to_numpy()
        adata.obs["module1_hepatocyte_lineage_candidate"] = aligned_seed["hepatocyte_lineage_candidate"].fillna(False).astype(bool).to_numpy()
        adata.obs["manual_major_label_cluster"] = adata.obs[args.cluster_key].astype(str).map(
            lambda x: cluster_label_map.get(x, {}).get("manual_major_label", "Unknown")
        )
        adata.obs["module1_confidence_status"] = adata.obs[args.cluster_key].astype(str).map(
            lambda x: cluster_label_map.get(x, {}).get("confidence_status", "unknown")
        )

        if args.use_rep not in adata.obsm:
            raise KeyError(f"{args.use_rep!r} is not in adata.obsm")
        print(f"NEIGHBORS use_rep={args.use_rep} n_neighbors={args.n_neighbors}", flush=True)
        sc.pp.neighbors(adata, use_rep=args.use_rep, n_neighbors=args.n_neighbors, random_state=args.seed_value, key_added="neighbors_hep")
        print(f"UMAP min_dist={args.min_dist}", flush=True)
        sc.tl.umap(adata, neighbors_key="neighbors_hep", min_dist=args.min_dist, random_state=args.seed_value)
        resolutions = sorted(set([args.resolution] + args.extra_resolution))
        for resolution in resolutions:
            key = f"leiden_hep_r{res_key(resolution)}"
            print(f"LEIDEN {key}", flush=True)
            sc.tl.leiden(adata, resolution=resolution, key_added=key, neighbors_key="neighbors_hep", random_state=args.seed_value)
        primary_key = f"leiden_hep_r{res_key(args.resolution)}"
        adata.obs["leiden_hep"] = adata.obs[primary_key].astype(str)

        manifest = read_manifest(args.manifest)
        marker_df, panel_df = aggregate_marker_panels(adata.obs.copy(), "leiden_hep", manifest)
    state_rows = []
    for subcluster, sub_panel in panel_df.groupby("leiden_hep", observed=True):
        state = assign_hepatocyte_state(sub_panel, args)
        sub_obs = adata.obs.loc[adata.obs["leiden_hep"].astype(str).eq(str(subcluster))]
        state.update(
            {
                "leiden_hep": str(subcluster),
                "n_cells": int(sub_obs.shape[0]),
                "top_original_leiden": top_counts(sub_obs[args.cluster_key], 5),
                "top_dataset": top_counts(sub_obs["dataset"], 5),
                "top_study_sample": top_counts(sub_obs["study_sample"], 5),
                "top_module1_label": top_counts(sub_obs["manual_major_label_cluster"], 5),
                "mean_celltypist_confidence": float(pd.to_numeric(sub_obs["celltypist_liver_confidence"], errors="coerce").mean()),
                "predicted_doublet_rate": float(sub_obs["predicted_doublet"].astype(bool).mean()) if "predicted_doublet" in sub_obs else np.nan,
                "cycling_rate": float(sub_obs["cell_cycle_phase"].astype(str).isin(["S", "G2M"]).mean()) if "cell_cycle_phase" in sub_obs else np.nan,
            }
        )
        for _, row in sub_panel.iterrows():
            prefix = str(row["marker_group"]).lower()
            state[f"{prefix}_score_z"] = float(row["panel_score_z"])
            state[f"{prefix}_mean_log1p_cpm"] = float(row["mean_log1p_cpm"])
            state[f"{prefix}_mean_pct_expr"] = float(row["mean_pct_expr"])
        state_rows.append(state)
    summary = pd.DataFrame(state_rows).sort_values("leiden_hep", key=lambda s: s.map(sort_key))

    state_map = summary.set_index("leiden_hep").to_dict(orient="index")
    adata.obs["hepatocyte_state_label"] = adata.obs["leiden_hep"].astype(str).map(
        lambda x: state_map.get(x, {}).get("hepatocyte_state_label", "ambiguous_epithelial_or_mixed")
    )
    adata.obs["hepatocyte_state_confidence"] = adata.obs["leiden_hep"].astype(str).map(
        lambda x: state_map.get(x, {}).get("hepatocyte_state_confidence", "ambiguous")
    )
    adata.obs["hepatocyte_state_seed_label"] = adata.obs["leiden_hep"].astype(str).map(
        lambda x: state_map.get(x, {}).get("hepatocyte_state_seed_label", "Unknown")
    )

    cells = adata.obs[
        [
            "dataset",
            "study_sample",
            args.cluster_key,
            "leiden_hep",
            "manual_major_label_cluster",
            "module1_confidence_status",
            "module1_scanvi_seed_label_major",
            "hepatocyte_state_label",
            "hepatocyte_state_confidence",
            "hepatocyte_state_seed_label",
            "predicted_doublet",
            "cell_cycle_phase",
        ]
    ].copy()
    cells.insert(0, "cell_id", adata.obs_names.astype(str))
    if "X_umap" in adata.obsm:
        cells["umap_hep_1"] = adata.obsm["X_umap"][:, 0]
        cells["umap_hep_2"] = adata.obsm["X_umap"][:, 1]

    marker_df.to_csv(marker_path, sep="\t", index=False)
    panel_df.to_csv(panel_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)
    cells.to_csv(cells_path, sep="\t", index=False, compression="gzip")
    cells[["cell_id", "leiden_hep", "hepatocyte_state_label", "hepatocyte_state_confidence", "hepatocyte_state_seed_label"]].to_csv(
        seed_path, sep="\t", index=False, compression="gzip"
    )
    cnv_candidates = cells.loc[
        cells["hepatocyte_state_confidence"].isin(["needs_cnv", "ambiguous"])
        | cells["hepatocyte_state_label"].astype(str).str.contains("needs_cnv", regex=False)
    ].copy()
    cnv_candidates.to_csv(cnv_candidates_path, sep="\t", index=False, compression="gzip")
    cnv_by_sample = (
        cnv_candidates.groupby(["dataset", "study_sample", "hepatocyte_state_label"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["dataset", "study_sample", "hepatocyte_state_label"])
    )
    cnv_by_sample.to_csv(cnv_by_sample_path, sep="\t", index=False)

    plot_panel_heatmap(panel_df, "leiden_hep", args.figures_dir)
    plot_marker_dotplot(marker_df, "leiden_hep", args.figures_dir)
    plot_umaps(adata, args)

    if not args.skip_h5ad:
        adata.write_h5ad(h5ad_path, compression="gzip")

    report = {
        "input": str(args.input.resolve()),
        "candidate_seed": str(args.seed.resolve()),
        "n_cells": int(adata.n_obs),
        "n_vars_hvg": int(adata.n_vars),
        "use_rep": args.use_rep,
        "n_neighbors": args.n_neighbors,
        "min_dist": args.min_dist,
        "resolution": args.resolution,
        "n_subclusters": int(summary.shape[0]),
        "state_counts_by_subcluster": summary["hepatocyte_state_label"].value_counts().to_dict(),
        "state_counts_by_cell": cells["hepatocyte_state_label"].value_counts().to_dict(),
        "state_seed_counts_by_cell": cells["hepatocyte_state_seed_label"].value_counts().to_dict(),
        "cnv_candidate_cells": int(cnv_candidates.shape[0]),
        "outputs": {
            "h5ad": str(h5ad_path.resolve()) if not args.skip_h5ad else None,
            "cells": str(cells_path.resolve()),
            "state_seed_by_cell": str(seed_path.resolve()),
            "cnv_candidate_cells": str(cnv_candidates_path.resolve()),
            "cnv_candidate_by_sample": str(cnv_by_sample_path.resolve()),
            "subcluster_summary": str(summary_path.resolve()),
            "marker_scores": str(marker_path.resolve()),
            "panel_scores": str(panel_path.resolve()),
            "umap_state": str((args.figures_dir / "hepatocyte_lineage_umap_state.png").resolve()),
            "umap_subcluster": str((args.figures_dir / "hepatocyte_lineage_umap_subcluster.png").resolve()),
            "marker_dotplot": str((args.figures_dir / "hepatocyte_state_marker_dotplot.png").resolve()),
            "panel_heatmap": str((args.figures_dir / "hepatocyte_state_panel_heatmap.png").resolve()),
        },
        "elapsed_seconds": round(time.time() - start, 3),
    }
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
