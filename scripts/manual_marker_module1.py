from __future__ import annotations

import argparse
import json
import math
import time
from collections import OrderedDict
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]


MARKER_PANELS: OrderedDict[str, list[str]] = OrderedDict(
    [
        ("Hepatocyte", ["ALB", "APOA1", "APOA2", "TTR", "HPD", "ASGR1", "CYP3A4", "CYP2E1"]),
        ("Cholangiocyte_Progenitor", ["KRT19", "EPCAM", "SOX9", "KRT7", "TACSTD2"]),
        ("Endothelial", ["PECAM1", "VWF", "KDR", "EMCN", "RAMP2"]),
        ("Fibroblast_HSC_Pericyte", ["COL1A1", "COL3A1", "ACTA2", "RGS5", "PDGFRB", "DCN", "COL1A2"]),
        ("Myeloid", ["PTPRC", "LST1", "C1QA", "C1QB", "CD68", "LYZ", "FCGR3A", "CD14"]),
        ("T_NK", ["CD3D", "CD3E", "NKG7", "GNLY", "KLRD1", "TRAC"]),
        ("B_cell", ["MS4A1", "CD79A", "CD79B", "BANK1"]),
        ("Plasma_cell", ["JCHAIN", "MZB1", "XBP1", "IGHG1", "IGHM"]),
        ("Proliferation", ["MKI67", "TOP2A", "STMN1", "TYMS", "UBE2C"]),
        ("HCC_Malignant_Associated", ["AFP", "GPC3", "SPP1", "MDK", "KRT19", "EPCAM"]),
    ]
)

PRIMARY_PANELS = {
    "Hepatocyte",
    "Cholangiocyte_Progenitor",
    "Endothelial",
    "Fibroblast_HSC_Pericyte",
    "Myeloid",
    "T_NK",
    "B_cell",
    "Plasma_cell",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 1 marker-based cluster review.")
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
    parser.add_argument(
        "--singler",
        type=Path,
        default=ROOT / "metadata/celltype/singler_combined_by_cluster.tsv",
    )
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/celltype")
    parser.add_argument("--figures-dir", type=Path, default=ROOT / "figures/celltype")
    parser.add_argument("--cluster-key", default="leiden_scvi")
    parser.add_argument("--exclude-column", default="excluded_doublet_cluster")
    parser.add_argument("--min-panel-z", type=float, default=0.35)
    parser.add_argument("--min-panel-pct", type=float, default=0.05)
    parser.add_argument("--external-conflict-confidence", type=float, default=0.60)
    parser.add_argument("--seed", type=int, default=20260601)
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


def as_csr(x) -> sparse.csr_matrix:
    if sparse.issparse(x):
        return x.tocsr()
    return sparse.csr_matrix(x)


def original_ids(global_index: pd.Index, study_sample: str) -> pd.Index:
    prefix = f"{study_sample}__"
    return pd.Index([idx[len(prefix) :] if str(idx).startswith(prefix) else idx for idx in global_index.astype(str)])


def sort_cluster_key(value: object) -> tuple[int, str]:
    text = str(value)
    return (0, f"{int(text):05d}") if text.isdigit() else (1, text)


def normalize_celltypist(value: object) -> str:
    text = str(value)
    mapping = {
        "Hepatocyte": "Hepatocyte",
        "Cholangiocyte": "Cholangiocyte_Progenitor",
        "Endothelial cell": "Endothelial",
        "Fibroblast/Stromal cell": "Fibroblast_HSC_Pericyte",
        "Myeloid cell": "Myeloid",
        "Granulocyte/Basophil": "Myeloid",
        "T/NK cell": "T_NK",
        "B cell": "B_cell",
        "Plasma/B cell": "Plasma_cell",
        "Doublet_suspect": "Unknown",
    }
    return mapping.get(text, "Unknown")


def normalize_singler(value: object, reference: str) -> str:
    text = str(value).lower()
    if text in {"", "na", "nan", "none"}:
        return "Unknown"
    if "hepatocyte" in text:
        return "Hepatocyte"
    if "epithelial" in text or "cholangi" in text:
        return "Cholangiocyte_Progenitor"
    if "endothelial" in text:
        return "Endothelial"
    if any(key in text for key in ["fibro", "smooth_muscle", "smooth muscle", "pericyte", "msc", "mesangial"]):
        return "Fibroblast_HSC_Pericyte"
    if any(key in text for key in ["mono", "macrophage", "dc", "neutrophil", "myelo", "eosinophil"]):
        return "Myeloid"
    if "nk" in text or "t_cell" in text or "t-cells" in text or "t.cells" in text:
        return "T_NK"
    if "plasma" in text:
        return "Plasma_cell"
    if "b_cell" in text or "b-cells" in text or "b.cells" in text:
        return "B_cell"
    if reference == "blueprint" and "adipocyte" in text:
        return "Noninformative"
    return "Unknown"


def equivalent(left: str, right: str) -> bool:
    if left == right:
        return True
    b_plasma = {"B_cell", "Plasma_cell"}
    if left in b_plasma and right in b_plasma:
        return True
    epithelial = {"Hepatocyte", "Cholangiocyte_Progenitor"}
    return left in epithelial and right in epithelial


def make_label_safe(label: str) -> str:
    return label.replace("/", "_").replace(" ", "_")


def expression_support(group: str, mean_log1p_cpm: float, mean_pct_expr: float) -> bool:
    cutoffs = {
        "Hepatocyte": (4.5, 0.25),
        "Cholangiocyte_Progenitor": (1.2, 0.04),
        "Endothelial": (2.5, 0.05),
        "Fibroblast_HSC_Pericyte": (2.5, 0.05),
        "Myeloid": (2.5, 0.05),
        "T_NK": (2.0, 0.05),
        "B_cell": (2.0, 0.05),
        "Plasma_cell": (2.0, 0.05),
        "HCC_Malignant_Associated": (1.5, 0.03),
    }
    min_log, min_pct = cutoffs.get(group, (2.5, 0.05))
    return mean_log1p_cpm >= min_log and mean_pct_expr >= min_pct


def panel_supported(panel_row: pd.Series, args: argparse.Namespace) -> bool:
    z_support = panel_row["panel_score_z"] >= args.min_panel_z and panel_row["mean_pct_expr"] >= args.min_panel_pct
    expr_support = expression_support(
        str(panel_row["marker_group"]),
        float(panel_row["mean_log1p_cpm"]),
        float(panel_row["mean_pct_expr"]),
    )
    return bool(z_support or expr_support)


def classify_cluster(row: pd.Series, args: argparse.Namespace) -> dict[str, object]:
    panel_rows = row["panel_rows"]
    primary = panel_rows.loc[panel_rows["marker_group"].isin(PRIMARY_PANELS)].copy()
    primary["is_supported"] = primary.apply(lambda x: panel_supported(x, args), axis=1)
    supported = primary.loc[primary["is_supported"]].copy()
    if supported.empty:
        primary = primary.sort_values(["panel_score_z", "mean_pct_expr"], ascending=False)
        top = primary.iloc[0]
        second = primary.iloc[1]
        marker_label = "Unknown"
        marker_reliable = False
    else:
        supported = supported.sort_values(["panel_score_z", "mean_log1p_cpm", "mean_pct_expr"], ascending=False)
        top = supported.iloc[0]
        if supported.shape[0] > 1:
            second = supported.iloc[1]
        else:
            fallback = primary.loc[primary["marker_group"].ne(top["marker_group"])].sort_values(
                ["panel_score_z", "mean_log1p_cpm", "mean_pct_expr"], ascending=False
            )
            second = fallback.iloc[0]
        marker_label = str(top["marker_group"])
        marker_reliable = True

    celltypist = normalize_celltypist(row.get("celltypist_major", "Unknown"))
    hpca = normalize_singler(row.get("singler_hpca_main_label", "Unknown"), "hpca")
    blueprint = normalize_singler(row.get("singler_blueprint_main_label", "Unknown"), "blueprint")

    informative = []
    for name, label in [("CellTypist", celltypist), ("SingleR_HPCA", hpca), ("SingleR_Blueprint", blueprint)]:
        if label not in {"Unknown", "Noninformative"}:
            informative.append((name, label))

    agree = [name for name, label in informative if marker_label != "Unknown" and equivalent(marker_label, label)]
    conflict = []
    celltypist_conf = float(row.get("mean_celltypist_confidence", 0) or 0)
    celltypist_frac = float(row.get("celltypist_major_fraction", 0) or 0)
    hpca_delta = float(row.get("singler_hpca_main_delta_next", 0) or 0)
    bp_delta = float(row.get("singler_blueprint_main_delta_next", 0) or 0)
    if marker_label != "Unknown":
        if celltypist not in {"Unknown", "Noninformative"} and not equivalent(marker_label, celltypist):
            if celltypist_conf >= args.external_conflict_confidence or celltypist_frac >= 0.75:
                conflict.append("CellTypist")
        if hpca not in {"Unknown", "Noninformative"} and not equivalent(marker_label, hpca) and hpca_delta >= 0.03:
            conflict.append("SingleR_HPCA")
        if blueprint not in {"Unknown", "Noninformative"} and not equivalent(marker_label, blueprint) and bp_delta >= 0.03:
            conflict.append("SingleR_Blueprint")
    hpca_discordant_only = False
    if conflict == ["SingleR_HPCA"] and agree and celltypist_conf >= 0.70:
        conflict = []
        hpca_discordant_only = True

    mixed_marker = bool(marker_reliable and panel_supported(second, args) and (top["panel_score_z"] - second["panel_score_z"]) < 0.25)
    if marker_label == "Unknown":
        status = "unknown"
        manual_label = "Unknown"
    elif conflict or mixed_marker:
        status = "conflict"
        manual_label = marker_label
    elif agree:
        status = "high_conf"
        manual_label = marker_label
    else:
        status = "unknown"
        manual_label = marker_label

    hep_row = panel_rows.loc[panel_rows["marker_group"].eq("Hepatocyte")].iloc[0]
    chol_row = panel_rows.loc[panel_rows["marker_group"].eq("Cholangiocyte_Progenitor")].iloc[0]
    hcc_row = panel_rows.loc[panel_rows["marker_group"].eq("HCC_Malignant_Associated")].iloc[0]
    hep_score = float(hep_row["panel_score_z"])
    chol_score = float(chol_row["panel_score_z"])
    hcc_score = float(hcc_row["panel_score_z"])
    epithelial_marker_support = bool(panel_supported(hep_row, args) or panel_supported(chol_row, args) or panel_supported(hcc_row, args))
    epithelial_labels = {"Hepatocyte", "Cholangiocyte_Progenitor"}
    epithelial_external_count = sum([celltypist in epithelial_labels, hpca in epithelial_labels, blueprint in epithelial_labels])
    if manual_label in epithelial_labels:
        hep_candidate = True
    elif status != "high_conf" and epithelial_external_count >= 2 and epithelial_marker_support:
        hep_candidate = True
    else:
        hep_candidate = False

    seed = manual_label if status == "high_conf" else "Unknown"
    reasons = [
        f"marker={marker_label}",
        f"top_z={top['panel_score_z']:.3f}",
        f"top_pct={top['mean_pct_expr']:.3f}",
        f"second={second['marker_group']}:{second['panel_score_z']:.3f}",
        f"celltypist={celltypist}",
        f"hpca={hpca}",
        f"blueprint={blueprint}",
    ]
    if agree:
        reasons.append("agree=" + ",".join(agree))
    if conflict:
        reasons.append("conflict=" + ",".join(conflict))
    if hpca_discordant_only:
        reasons.append("hpca_discordant_only=true")
    if mixed_marker:
        reasons.append("mixed_marker=true")

    return {
        "manual_marker_top_label": marker_label,
        "manual_major_label": manual_label,
        "confidence_status": status,
        "scanvi_seed_label_major": seed,
        "hepatocyte_lineage_candidate": hep_candidate,
        "marker_top_score_z": float(top["panel_score_z"]),
        "marker_second_label": str(second["marker_group"]),
        "marker_second_score_z": float(second["panel_score_z"]),
        "marker_top_mean_pct_expr": float(top["mean_pct_expr"]),
        "normalized_celltypist_label": celltypist,
        "normalized_hpca_label": hpca,
        "normalized_blueprint_label": blueprint,
        "manual_evidence": "; ".join(reasons),
    }


def plot_marker_dotplot(marker_df: pd.DataFrame, clusters: list[str], figures_dir: Path) -> None:
    genes = []
    gene_to_group = {}
    for group, group_genes in MARKER_PANELS.items():
        for gene in group_genes:
            if gene not in genes:
                genes.append(gene)
                gene_to_group[gene] = group
    data = marker_df.loc[marker_df["gene"].isin(genes)].copy()
    data["cluster"] = data["cluster"].astype(str)
    data["gene"] = pd.Categorical(data["gene"], categories=genes, ordered=True)
    data["cluster"] = pd.Categorical(data["cluster"], categories=clusters, ordered=True)
    data = data.sort_values(["cluster", "gene"])

    x = data["gene"].cat.codes.to_numpy()
    y = data["cluster"].cat.codes.to_numpy()
    size = np.clip(data["pct_expr"].fillna(0).to_numpy(), 0, 1) * 180 + 2
    color = data["log1p_cpm"].fillna(0).to_numpy()

    width = max(16, len(genes) * 0.33)
    height = max(12, len(clusters) * 0.22)
    fig, ax = plt.subplots(figsize=(width, height))
    scatter = ax.scatter(x, y, s=size, c=color, cmap="viridis", edgecolors="none")
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=90, fontsize=7)
    ax.set_yticks(range(len(clusters)))
    ax.set_yticklabels(clusters, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Marker gene")
    ax.set_ylabel("Leiden cluster")
    ax.set_title("Classic marker expression by cluster")
    cbar = fig.colorbar(scatter, ax=ax, pad=0.01)
    cbar.set_label("log1p(CPM)")

    start = 0
    for group, group_genes in MARKER_PANELS.items():
        present = [gene for gene in group_genes if gene in genes]
        if not present:
            continue
        end = start + len(present)
        ax.axvline(end - 0.5, color="#dddddd", linewidth=0.6)
        ax.text((start + end - 1) / 2, -2.0, group, ha="center", va="bottom", rotation=45, fontsize=7)
        start = end

    figures_dir.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(figures_dir / "manual_marker_dotplot.png", dpi=300)
    fig.savefig(figures_dir / "manual_marker_dotplot.pdf")
    plt.close(fig)


def plot_panel_heatmap(panel_df: pd.DataFrame, clusters: list[str], figures_dir: Path) -> None:
    groups = list(MARKER_PANELS.keys())
    matrix = (
        panel_df.pivot(index="cluster", columns="marker_group", values="panel_score_z")
        .reindex(index=clusters, columns=groups)
        .fillna(0)
    )
    fig, ax = plt.subplots(figsize=(max(8, len(groups) * 0.7), max(12, len(clusters) * 0.22)))
    image = ax.imshow(matrix.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(clusters)))
    ax.set_yticklabels(clusters, fontsize=7)
    ax.set_xlabel("Marker panel")
    ax.set_ylabel("Leiden cluster")
    ax.set_title("Panel score z-score by cluster")
    cbar = fig.colorbar(image, ax=ax, pad=0.01)
    cbar.set_label("Mean gene z-score")
    fig.tight_layout()
    fig.savefig(figures_dir / "manual_marker_panel_score_heatmap.png", dpi=300)
    fig.savefig(figures_dir / "manual_marker_panel_score_heatmap.pdf")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    start = time.time()
    np.random.seed(args.seed)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    args.figures_dir.mkdir(parents=True, exist_ok=True)

    singler = pd.read_csv(args.singler, sep="\t")
    singler["leiden_scvi"] = singler["leiden_scvi"].astype(str)
    target_clusters = sorted(singler["leiden_scvi"].unique(), key=sort_cluster_key)
    cluster_to_col = {cluster: i for i, cluster in enumerate(target_clusters)}

    backed = ad.read_h5ad(args.input, backed="r")
    obs = backed.obs.copy()
    backed.file.close()
    obs[args.cluster_key] = obs[args.cluster_key].astype(str)
    obs_target = obs.loc[obs[args.cluster_key].isin(target_clusters)].copy()
    obs_target["qc_total_counts"] = pd.to_numeric(obs_target["qc_total_counts"], errors="coerce").fillna(0.0)
    print(f"CLUSTERS {len(target_clusters)} CELLS {obs_target.shape[0]}", flush=True)

    marker_genes = []
    marker_to_groups: dict[str, list[str]] = {}
    for group, genes in MARKER_PANELS.items():
        for gene in genes:
            if gene not in marker_to_groups:
                marker_genes.append(gene)
                marker_to_groups[gene] = []
            marker_to_groups[gene].append(group)
    gene_to_idx = {gene: i for i, gene in enumerate(marker_genes)}
    n_genes = len(marker_genes)
    n_clusters = len(target_clusters)

    raw_sums = np.zeros((n_genes, n_clusters), dtype=np.float64)
    detected = np.zeros((n_genes, n_clusters), dtype=np.float64)
    cell_denoms = np.zeros((n_genes, n_clusters), dtype=np.float64)
    count_denoms = np.zeros((n_genes, n_clusters), dtype=np.float64)
    sample_presence = np.zeros((n_genes, n_clusters), dtype=np.float64)

    manifest = read_manifest(args.manifest)
    for i, row in manifest.reset_index(drop=True).iterrows():
        study_sample = str(row["study_sample"])
        sample_obs = obs_target.loc[obs_target["study_sample"].astype(str).eq(study_sample)]
        if sample_obs.empty:
            continue
        path = resolve_path(row["output"])
        print(f"SAMPLE {i + 1}/{manifest.shape[0]} {study_sample} cells={sample_obs.shape[0]}", flush=True)
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
        sample_cluster_cols = sample_obs[args.cluster_key].astype(str).map(cluster_to_col).to_numpy()
        indicator = sparse.csr_matrix(
            (np.ones(len(present_cells), dtype=np.float64), (np.arange(len(present_cells)), sample_cluster_cols)),
            shape=(len(present_cells), n_clusters),
        )

        x = as_csr(a[present_cells, present_genes].X).astype(np.float64)
        raw = (x.T @ indicator).toarray()
        det = ((x > 0).astype(np.float64).T @ indicator).toarray()
        gene_locs = np.array([gene_to_idx[gene] for gene in present_genes], dtype=int)
        raw_sums[gene_locs, :] += raw
        detected[gene_locs, :] += det

        sample_cluster_cells = np.bincount(sample_cluster_cols, minlength=n_clusters).astype(np.float64)
        sample_cluster_counts = (
            sample_obs.groupby(sample_obs[args.cluster_key].astype(str), observed=True)["qc_total_counts"].sum()
            .reindex(target_clusters)
            .fillna(0.0)
            .to_numpy(dtype=np.float64)
        )
        for loc in gene_locs:
            cell_denoms[loc, :] += sample_cluster_cells
            count_denoms[loc, :] += sample_cluster_counts
            sample_presence[loc, :] += (sample_cluster_cells > 0).astype(np.float64)
        a.file.close()

    cpm = np.divide(raw_sums, count_denoms, out=np.zeros_like(raw_sums), where=count_denoms > 0) * 1e6
    log1p_cpm = np.log1p(cpm)
    pct_expr = np.divide(detected, cell_denoms, out=np.zeros_like(detected), where=cell_denoms > 0)

    marker_rows = []
    for gene in marker_genes:
        loc = gene_to_idx[gene]
        for cluster in target_clusters:
            col = cluster_to_col[cluster]
            for group in marker_to_groups[gene]:
                marker_rows.append(
                    {
                        "cluster": cluster,
                        "marker_group": group,
                        "gene": gene,
                        "raw_sum": raw_sums[loc, col],
                        "detected_cells": detected[loc, col],
                        "n_cells_gene_present": cell_denoms[loc, col],
                        "pct_expr": pct_expr[loc, col],
                        "log1p_cpm": log1p_cpm[loc, col],
                        "gene_present_sample_count": sample_presence[loc, col],
                    }
                )
    marker_df = pd.DataFrame(marker_rows)

    gene_z = np.zeros_like(log1p_cpm)
    for i in range(n_genes):
        valid = cell_denoms[i, :] > 0
        values = log1p_cpm[i, valid]
        if values.size > 1 and float(values.std()) > 0:
            gene_z[i, valid] = (values - values.mean()) / values.std()

    panel_rows = []
    for cluster in target_clusters:
        col = cluster_to_col[cluster]
        for group, genes in MARKER_PANELS.items():
            locs = [gene_to_idx[gene] for gene in genes if gene in gene_to_idx and cell_denoms[gene_to_idx[gene], col] > 0]
            if locs:
                panel_rows.append(
                    {
                        "cluster": cluster,
                        "marker_group": group,
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
                        "cluster": cluster,
                        "marker_group": group,
                        "panel_score_z": 0.0,
                        "mean_log1p_cpm": 0.0,
                        "mean_pct_expr": 0.0,
                        "n_genes_available": 0,
                        "genes_available": "",
                    }
                )
    panel_df = pd.DataFrame(panel_rows)

    merged = singler.copy()
    merged["leiden_scvi"] = merged["leiden_scvi"].astype(str)
    cluster_meta = []
    for _, row in merged.iterrows():
        cluster = str(row["leiden_scvi"])
        row_dict = row.to_dict()
        row_dict["panel_rows"] = panel_df.loc[panel_df["cluster"].astype(str).eq(cluster)]
        row_dict.update(classify_cluster(pd.Series(row_dict), args))
        cluster_meta.append(row_dict)
    label_df = pd.DataFrame(cluster_meta).drop(columns=["panel_rows"])
    label_df = label_df.sort_values("leiden_scvi", key=lambda s: s.map(sort_cluster_key))

    marker_path = args.metadata_dir / "manual_marker_cluster_scores.tsv"
    panel_path = args.metadata_dir / "manual_marker_panel_scores.tsv"
    labels_path = args.metadata_dir / "manual_major_label_by_cluster.tsv"
    seed_path = args.metadata_dir / "scanvi_seed_labels_by_cell.tsv.gz"
    hep_path = args.metadata_dir / "hepatocyte_lineage_candidate_clusters.tsv"
    report_path = args.metadata_dir / "manual_marker_module1_report.json"

    marker_df.to_csv(marker_path, sep="\t", index=False)
    panel_df.to_csv(panel_path, sep="\t", index=False)
    label_df.to_csv(labels_path, sep="\t", index=False)
    label_df.loc[label_df["hepatocyte_lineage_candidate"].astype(bool)].to_csv(hep_path, sep="\t", index=False)

    label_map = label_df.set_index("leiden_scvi").to_dict(orient="index")
    seed_rows = pd.DataFrame(index=obs.index)
    seed_rows["cell_id"] = obs.index.astype(str)
    seed_rows["leiden_scvi"] = obs[args.cluster_key].astype(str).to_numpy()
    seed_rows["manual_major_label_cluster"] = seed_rows["leiden_scvi"].map(
        lambda x: label_map.get(x, {}).get("manual_major_label", "Doublet_suspect" if x == "16" else "Unknown")
    )
    seed_rows["manual_confidence_status"] = seed_rows["leiden_scvi"].map(
        lambda x: label_map.get(x, {}).get("confidence_status", "excluded_or_unscored" if x == "16" else "unknown")
    )
    seed_rows["scanvi_seed_label_major"] = seed_rows["leiden_scvi"].map(
        lambda x: label_map.get(x, {}).get("scanvi_seed_label_major", "Unknown")
    )
    seed_rows["hepatocyte_lineage_candidate"] = seed_rows["leiden_scvi"].map(
        lambda x: bool(label_map.get(x, {}).get("hepatocyte_lineage_candidate", False))
    )
    if args.exclude_column in obs.columns:
        seed_rows["excluded_doublet_cluster"] = obs[args.exclude_column].astype(bool).to_numpy()
    seed_rows.to_csv(seed_path, sep="\t", index=False, compression="gzip")

    plot_marker_dotplot(marker_df, target_clusters, args.figures_dir)
    plot_panel_heatmap(panel_df, target_clusters, args.figures_dir)

    report = {
        "input": str(args.input.resolve()),
        "manifest": str(args.manifest.resolve()),
        "singler": str(args.singler.resolve()),
        "n_clusters_scored": int(len(target_clusters)),
        "n_cells_scored": int(obs_target.shape[0]),
        "n_marker_genes": int(len(marker_genes)),
        "confidence_status_counts": label_df["confidence_status"].value_counts().to_dict(),
        "manual_major_label_counts": label_df["manual_major_label"].value_counts().to_dict(),
        "scanvi_seed_label_major_counts_by_cluster": label_df["scanvi_seed_label_major"].value_counts().to_dict(),
        "hepatocyte_lineage_candidate_clusters": label_df.loc[
            label_df["hepatocyte_lineage_candidate"].astype(bool), "leiden_scvi"
        ].astype(str).tolist(),
        "thresholds": {
            "min_panel_z": args.min_panel_z,
            "min_panel_pct": args.min_panel_pct,
            "external_conflict_confidence": args.external_conflict_confidence,
        },
        "outputs": {
            "marker_scores": str(marker_path.resolve()),
            "panel_scores": str(panel_path.resolve()),
            "cluster_labels": str(labels_path.resolve()),
            "per_cell_scanvi_seed": str(seed_path.resolve()),
            "hepatocyte_lineage_candidates": str(hep_path.resolve()),
            "dotplot_png": str((args.figures_dir / "manual_marker_dotplot.png").resolve()),
            "panel_heatmap_png": str((args.figures_dir / "manual_marker_panel_score_heatmap.png").resolve()),
        },
        "elapsed_seconds": round(time.time() - start, 3),
    }
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
