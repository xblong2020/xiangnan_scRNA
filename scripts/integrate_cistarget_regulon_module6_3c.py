from __future__ import annotations

import argparse
import ast
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
from scipy.stats import pearsonr, spearmanr


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 6.3c: integrate cisTarget-pruned regulons with CellRank fate.")
    parser.add_argument(
        "--input-h5ad",
        type=Path,
        default=ROOT / "data/processed/driver/driver_union_full_expression.module6_3b.h5ad",
    )
    parser.add_argument(
        "--motifs",
        type=Path,
        default=ROOT / "metadata/driver/driver_module6_3c_cistarget_motifs.tsv",
    )
    parser.add_argument(
        "--auc",
        type=Path,
        default=ROOT / "metadata/driver/driver_module6_3c_cistarget_regulon_auc.csv",
    )
    parser.add_argument(
        "--ranking-db",
        type=Path,
        default=ROOT
        / "metadata/driver/scenic_resources/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather",
    )
    parser.add_argument("--processed-dir", type=Path, default=ROOT / "data/processed/driver")
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/driver")
    parser.add_argument("--figures-dir", type=Path, default=ROOT / "figures/driver")
    parser.add_argument("--output-name", default="driver_cistarget_regulon_activity.module6_3c.h5ad")
    parser.add_argument("--fate-key", default="cellrank_fate_prob_cnv_supported_malignant")
    parser.add_argument("--time-key", default="driver_main_strict__pseudotime_median")
    parser.add_argument("--phase-key", default="driver_main_strict__pseudotime_phase")
    parser.add_argument("--top-n", type=int, default=30)
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


def flatten_motif_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if isinstance(out.columns, pd.MultiIndex):
        flattened = []
        for top, bottom in out.columns.to_list():
            bottom_text = str(bottom)
            top_text = str(top)
            flattened.append(bottom_text if not bottom_text.startswith("Unnamed") else top_text)
        out.columns = flattened
    return out


def read_motifs(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    try:
        motifs = pd.read_csv(path, sep="\t", header=[0, 1], index_col=[0, 1])
        motifs = flatten_motif_columns(motifs)
    except Exception:
        motifs = pd.read_csv(path, sep="\t")
    if motifs.empty:
        raise ValueError(f"cisTarget motif table is empty: {path}")
    return motifs


def read_auc_matrix(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    auc = pd.read_csv(path)
    if auc.shape[1] == 1:
        auc = pd.read_csv(path, sep="\t")
    if auc.empty:
        raise ValueError(f"AUCell matrix is empty: {path}")
    return auc


def align_auc_to_cells(auc: pd.DataFrame, cell_ids: pd.Index) -> pd.DataFrame:
    work = auc.copy()
    cell_col = None
    for candidate in ["CellID", "cell_id", "Cell", "cells"]:
        if candidate in work.columns:
            cell_col = candidate
            break
    if cell_col is None and str(work.columns[0]).startswith("Unnamed"):
        cell_col = work.columns[0]
    if cell_col is not None:
        work.index = work.pop(cell_col).astype(str)
    else:
        work.index = work.index.astype(str)
    work.columns = work.columns.astype(str)
    missing = pd.Index(cell_ids.astype(str)).difference(pd.Index(work.index.astype(str)))
    if len(missing) > 0:
        raise ValueError(f"{len(missing)} cells from h5ad are missing in AUCell matrix.")
    aligned = work.loc[cell_ids.astype(str)].apply(pd.to_numeric, errors="coerce")
    variable = aligned.nunique(dropna=True) > 1
    aligned = aligned.loc[:, variable]
    if aligned.empty:
        raise ValueError("No variable regulon AUC columns remain after alignment.")
    return aligned


def target_genes_from_value(value: object) -> set[str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return set()
    if isinstance(value, (list, tuple, set, frozenset)):
        parsed = value
    else:
        text = str(value)
        try:
            parsed = ast.literal_eval(text)
        except (ValueError, SyntaxError):
            return {part.strip() for part in text.replace(";", ",").split(",") if part.strip()}
    genes = set()
    for item in parsed:
        if isinstance(item, (list, tuple)) and item:
            genes.add(str(item[0]))
        elif isinstance(item, str):
            genes.add(item)
    return genes


def first_existing(columns: pd.Index, candidates: list[str]) -> str:
    for candidate in candidates:
        if candidate in columns:
            return candidate
    raise KeyError(f"Missing expected columns; tried {candidates}")


def summarize_pruned_regulons(motifs: pd.DataFrame) -> pd.DataFrame:
    work = flatten_motif_columns(motifs)
    if not any(candidate in work.columns for candidate in ["TF", "TranscriptionFactor", "tf"]):
        work = work.reset_index()
    tf_col = first_existing(work.columns, ["TF", "TranscriptionFactor", "tf"])
    nes_col = first_existing(work.columns, ["NES", "nes"])
    target_col = first_existing(work.columns, ["TargetGenes", "target_genes", "targets"])
    motif_col = next((col for col in ["MotifID", "motif", "motif_id", "Feature"] if col in work.columns), None)
    annotation_col = next((col for col in ["Annotation", "annotation"] if col in work.columns), None)

    rows = []
    for tf, sub in work.groupby(tf_col, sort=True):
        targets: set[str] = set()
        for value in sub[target_col]:
            targets.update(target_genes_from_value(value))
        nes = pd.to_numeric(sub[nes_col], errors="coerce")
        best_idx = nes.idxmax() if nes.notna().any() else sub.index[0]
        annotations = sub[annotation_col].astype(str) if annotation_col else pd.Series([], dtype=str)
        direct = annotations.str.contains("direct", case=False, na=False).sum() if not annotations.empty else 0
        rows.append(
            {
                "regulon": f"{tf}(+)",
                "tf": str(tf),
                "n_motifs": int(sub.shape[0]),
                "n_targets": int(len(targets)),
                "best_nes": float(nes.loc[best_idx]) if pd.notna(nes.loc[best_idx]) else np.nan,
                "top_motif_id": str(sub.loc[best_idx, motif_col]) if motif_col else "",
                "top_annotation": str(sub.loc[best_idx, annotation_col]) if annotation_col else "",
                "direct_annotation_fraction": float(direct / sub.shape[0]) if sub.shape[0] > 0 else np.nan,
                "target_genes": ";".join(sorted(targets)),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError("No pruned regulons could be summarized from motif table.")
    return out.sort_values(["best_nes", "n_targets"], ascending=[False, False]).reset_index(drop=True)


def correlate_auc_with_fate(auc: pd.DataFrame, cells: pd.DataFrame, fate_key: str, time_key: str) -> pd.DataFrame:
    if fate_key not in cells.columns:
        raise KeyError(f"Missing fate key: {fate_key}")
    fate = pd.to_numeric(cells[fate_key], errors="coerce")
    pseudotime = pd.to_numeric(cells[time_key], errors="coerce") if time_key in cells.columns else pd.Series(np.nan, index=cells.index)
    rows = []
    for regulon in auc.columns.astype(str):
        values = pd.to_numeric(auc[regulon], errors="coerce")
        row = {"regulon": regulon}
        for label, target in [("cnv_fate", fate), ("pseudotime", pseudotime)]:
            mask = values.notna() & target.notna() & np.isfinite(values) & np.isfinite(target)
            row[f"n_{label}_cells"] = int(mask.sum())
            if label == "cnv_fate":
                row["n_fate_cells"] = int(mask.sum())
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
            valid_idx = pd.Index(idx).intersection(auc.index)
            for regulon in top_regulons:
                values = pd.to_numeric(auc.loc[valid_idx, regulon], errors="coerce").dropna()
                rows.append(
                    {
                        "group_type": group_col,
                        "group": str(group),
                        "regulon": regulon,
                        "n_cells": int(len(valid_idx)),
                        "mean_auc": float(values.mean()) if not values.empty else np.nan,
                        "median_auc": float(values.median()) if not values.empty else np.nan,
                    }
                )
    return pd.DataFrame(rows)


def plot_top_regulons(regulons: pd.DataFrame, top_n: int, path_base: Path) -> list[str]:
    sub = regulons.sort_values("cnv_fate_pearson_r", ascending=False).head(top_n).iloc[::-1]
    if sub.empty:
        return []
    fig, ax = plt.subplots(figsize=(5.6, max(3.0, 0.22 * sub.shape[0])))
    ax.barh(sub["regulon"], sub["cnv_fate_pearson_r"], color="#0072B2", height=0.75)
    ax.set_xlabel("Pearson r with CellRank CNV fate probability")
    ax.set_ylabel("")
    ax.set_title(f"Top {sub.shape[0]} cisTarget-pruned regulons")
    fig.tight_layout()
    outputs = []
    for suffix in [".png", ".pdf"]:
        out = path_base.with_suffix(suffix)
        fig.savefig(out, bbox_inches="tight")
        outputs.append(str(out))
    plt.close(fig)
    return outputs


def plot_phase_heatmap(auc: pd.DataFrame, cells: pd.DataFrame, regulons: list[str], phase_key: str, path_base: Path) -> list[str]:
    if not regulons or phase_key not in cells.columns:
        return []
    rows = []
    for phase in ["early", "middle", "late"]:
        idx = cells.index[cells[phase_key].astype(str).eq(phase)]
        idx = pd.Index(idx).intersection(auc.index)
        if len(idx) > 0:
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
    ax.set_title("Mean cisTarget regulon AUC by pseudotime phase")
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


def plot_group_activity(summary: pd.DataFrame, group_type: str, regulons: list[str], path_base: Path) -> list[str]:
    sub = summary.loc[summary["group_type"].eq(group_type) & summary["regulon"].isin(regulons)].copy()
    if sub.empty:
        return []
    pivot = sub.pivot_table(index="group", columns="regulon", values="mean_auc", aggfunc="mean")
    pivot = pivot.loc[:, [reg for reg in regulons if reg in pivot.columns]]
    fig, ax = plt.subplots(figsize=(max(5.0, 0.32 * pivot.shape[1]), max(3.0, 0.24 * pivot.shape[0])))
    im = ax.imshow(pivot.to_numpy(dtype=float), aspect="auto", cmap="viridis")
    ax.set_yticks(range(pivot.shape[0]))
    ax.set_yticklabels(pivot.index)
    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, rotation=90)
    ax.set_title(f"Mean regulon AUC by {group_type}")
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


def main() -> None:
    start = time.time()
    args = parse_args()
    configure_plot_style()
    args.processed_dir.mkdir(parents=True, exist_ok=True)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    args.figures_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.input_h5ad)
    motifs = read_motifs(args.motifs)
    auc_raw = read_auc_matrix(args.auc)
    auc = align_auc_to_cells(auc_raw, pd.Index(adata.obs_names.astype(str)))
    cells = adata.obs.copy()
    cells.index = adata.obs_names.astype(str)

    regulons = summarize_pruned_regulons(motifs)
    corr = correlate_auc_with_fate(auc, cells, args.fate_key, args.time_key)
    regulons = regulons.merge(corr, on="regulon", how="inner")
    regulons = regulons.sort_values(["cnv_fate_pearson_r", "best_nes"], ascending=[False, False]).reset_index(drop=True)
    top_regulons = regulons.head(args.top_n)["regulon"].tolist()
    group_summary = summarize_activity_by_group(
        auc,
        cells,
        top_regulons=top_regulons,
        group_cols=[args.phase_key, "cell_disease_stage", "trajectory_root_end_role", "dataset", "sample_id"],
    )

    adata.obsm["module6_3c_cistarget_regulon_auc"] = auc.to_numpy(dtype=np.float32)
    adata.uns["module6_3c_cistarget_regulon_auc_names"] = auc.columns.tolist()
    for regulon in top_regulons[:10]:
        safe = regulon.replace("(", "_").replace(")", "").replace("+", "plus").replace("-", "minus")
        adata.obs[f"module6_3c_auc_{safe}"] = pd.to_numeric(auc[regulon], errors="coerce").to_numpy(dtype=float)
    adata.uns["module6_3c_cistarget_regulon_activity"] = {
        "module": "6.3c",
        "method": "cisTarget motif-pruned regulons scored with pySCENIC AUCell and correlated with CellRank CNV fate",
        "ranking_db": str(args.ranking_db),
        "n_cells": int(adata.n_obs),
        "n_regulons": int(auc.shape[1]),
        "n_motif_rows": int(motifs.shape[0]),
        "n_fate_non_null_cells": int(cells[args.fate_key].notna().sum()) if args.fate_key in cells.columns else 0,
    }

    output_h5ad = args.processed_dir / args.output_name
    adata.write_h5ad(output_h5ad, compression=args.compression)

    auc_cells = auc.copy()
    auc_cells.insert(0, "cell_id", auc_cells.index.astype(str))
    auc_path = args.metadata_dir / "driver_module6_3c_cistarget_regulon_auc.tsv.gz"
    regulons_path = args.metadata_dir / "driver_module6_3c_cistarget_regulon_summary.tsv"
    top_path = args.metadata_dir / "driver_module6_3c_top_cnv_fate_regulons.tsv"
    summary_path = args.metadata_dir / "driver_module6_3c_regulon_activity_group_summary.tsv"

    write_dataframe(auc_path, auc_cells.reset_index(drop=True))
    write_dataframe(regulons_path, regulons)
    write_dataframe(top_path, regulons.head(args.top_n))
    write_dataframe(summary_path, group_summary)

    figure_outputs: list[str] = []
    figure_outputs += plot_top_regulons(regulons, args.top_n, args.figures_dir / "driver_module6_3c_top_cnv_fate_regulons")
    figure_outputs += plot_phase_heatmap(
        auc,
        cells,
        regulons=top_regulons[:25],
        phase_key=args.phase_key,
        path_base=args.figures_dir / "driver_module6_3c_top_regulon_phase_heatmap",
    )
    figure_outputs += plot_group_activity(
        group_summary,
        group_type="dataset",
        regulons=top_regulons[:15],
        path_base=args.figures_dir / "driver_module6_3c_top_regulon_dataset_heatmap",
    )

    report = {
        "module": "6.3c",
        "method": "cisTarget-pruned SCENIC regulon integration after CellRank",
        "input_h5ad": str(args.input_h5ad),
        "motifs": str(args.motifs),
        "auc_input": str(args.auc),
        "ranking_db": str(args.ranking_db),
        "output_h5ad": str(output_h5ad),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_regulons_scored": int(auc.shape[1]),
        "n_motif_rows": int(motifs.shape[0]),
        "n_pruned_regulons": int(regulons.shape[0]),
        "n_fate_non_null_cells": int(cells[args.fate_key].notna().sum()) if args.fate_key in cells.columns else 0,
        "top_cnv_fate_regulons": regulons.head(15).to_dict(orient="records"),
        "outputs": {
            "h5ad": str(output_h5ad),
            "auc": str(auc_path),
            "regulon_summary": str(regulons_path),
            "top_regulons": str(top_path),
            "group_summary": str(summary_path),
            "figures": figure_outputs,
        },
        "package_versions": {
            "anndata": version("anndata"),
            "pandas": version("pandas"),
            "numpy": version("numpy"),
            "scipy": version("scipy"),
            "pyscenic": version("pyscenic"),
            "ctxcore": version("ctxcore"),
        },
        "elapsed_seconds": round(time.time() - start, 3),
    }
    report_path = args.metadata_dir / "driver_module6_3c_cistarget_regulon_report.json"
    report["outputs"]["report"] = str(report_path)
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
