from __future__ import annotations

import argparse
import json
import time
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import celltypist
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


ROOT = Path(__file__).resolve().parents[1]


LIVER_TO_MAJOR = {
    "B cells": "B cell",
    "Basophils": "Granulocyte/Basophil",
    "Cholangiocytes": "Cholangiocyte",
    "Circulating NK/NKT": "T/NK cell",
    "Endothelial cells": "Endothelial cell",
    "Fibroblasts": "Fibroblast/Stromal cell",
    "Hepatocytes": "Hepatocyte",
    "Macrophages": "Myeloid cell",
    "Mig.cDCs": "Myeloid cell",
    "Mono+mono derived cells": "Myeloid cell",
    "Neutrophils": "Granulocyte/Basophil",
    "Plasma cells": "Plasma/B cell",
    "Resident NK": "T/NK cell",
    "T cells": "T/NK cell",
    "cDC1s": "Myeloid cell",
    "cDC2s": "Myeloid cell",
    "pDCs": "Myeloid cell",
}


SOURCE_LABEL_CANDIDATES = [
    "cell.annotation",
    "annotation_refined",
    "annotation",
    "celltype",
    "Type",
    "Global_Cluster",
    "Release_Global_Cluster",
    "seurat_clusters",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Coarse major-cell annotation with CellTypist.")
    parser.add_argument(
        "--input",
        type=Path,
        default=ROOT / "data/processed/scvi/scvi_integrated_counts_hvg.scvi_doublet_cell_cycle.h5ad",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=ROOT / "metadata/scvi/scvi_input_manifest.counts.tsv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "data/processed/scvi/scvi_integrated_counts_hvg.celltypist_major.h5ad",
    )
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/celltype")
    parser.add_argument("--figures-dir", type=Path, default=ROOT / "figures/celltype")
    parser.add_argument("--model", default="Healthy_Human_Liver.pkl")
    parser.add_argument("--exclude-cluster-key", default="leiden_scvi")
    parser.add_argument("--exclude-cluster", action="append", default=["16"])
    parser.add_argument("--cluster-key", default="leiden_scvi")
    parser.add_argument("--max-plot-cells", type=int, default=250000)
    parser.add_argument("--seed", type=int, default=20260601)
    parser.add_argument("--only-study-sample", action="append", default=None)
    parser.add_argument("--max-samples", type=int, default=None)
    parser.add_argument("--skip-h5ad", action="store_true")
    return parser.parse_args()


def read_manifest(path: Path) -> pd.DataFrame:
    manifest = pd.read_csv(path, sep="\t")
    include = manifest["include_in_scvi"].astype(str).str.lower().eq("true")
    manifest = manifest.loc[include].copy()
    manifest["study_sample"] = manifest["dataset"].astype(str) + "__" + manifest["label"].astype(str)
    return manifest


def original_ids_for_sample(obs: pd.DataFrame, study_sample: str) -> pd.Index:
    prefix = f"{study_sample}__"
    idx = obs.index[obs["study_sample"].astype(str) == study_sample]
    return pd.Index([name[len(prefix) :] if name.startswith(prefix) else name for name in idx])


def select_source_label_column(obs: pd.DataFrame) -> str | None:
    for col in SOURCE_LABEL_CANDIDATES:
        if col not in obs.columns:
            continue
        values = obs[col].astype(str)
        nunique = values.nunique(dropna=True)
        if 1 < nunique < min(500, max(3, obs.shape[0] // 2)):
            return col
    return None


def source_label_to_major(label: object) -> str:
    text = str(label).strip()
    low = text.lower()
    if text in {"", "NA", "nan", "None"}:
        return "Unknown"
    if any(key in low for key in ["doublet", "mixed"]):
        return "Doublet_suspect"
    if any(key in low for key in ["hepatocyte", "hepato", "malignant", "tumor", "epithelial"]):
        return "Hepatocyte"
    if any(key in low for key in ["cholangiocyte", "cholangio", "biliary"]):
        return "Cholangiocyte"
    if any(key in low for key in ["endo", "vascular"]):
        return "Endothelial cell"
    if any(key in low for key in ["fibro", "stellate", "stromal", "smooth muscle", "mesench"]):
        return "Fibroblast/Stromal cell"
    if any(key in low for key in ["plasma"]):
        return "Plasma/B cell"
    if low in {"b", "b cell", "b cells"} or any(key in low for key in ["bcell", "b cell", "b cells"]):
        return "B cell"
    if any(key in low for key in ["t/nk", "tcell", "t cell", "t cells", "nk", "nkt", "lymphocyte", "lymphoid"]):
        return "T/NK cell"
    if any(key in low for key in ["myeloid", "mono", "macro", "kupffer", "dendritic", "dc", "mp", "mφ"]):
        return "Myeloid cell"
    if any(key in low for key in ["neutrophil", "basophil", "mast", "granulo"]):
        return "Granulocyte/Basophil"
    if any(key in low for key in ["eryth", "red blood", "rbc"]):
        return "Erythroid cell"
    return "Unknown"


def liver_label_to_major(label: object) -> str:
    return LIVER_TO_MAJOR.get(str(label), "Unknown")


def build_sample_list(integrated_obs: pd.DataFrame, manifest: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    samples = integrated_obs[["dataset", "sample_id", "study_sample", "source_h5ad"]].drop_duplicates()
    samples = samples.merge(
        manifest[["study_sample", "output"]],
        on="study_sample",
        how="left",
    )
    samples["source_path"] = samples["output"].fillna(samples["source_h5ad"])
    samples = samples.sort_values(["dataset", "sample_id"]).reset_index(drop=True)
    if args.only_study_sample:
        samples = samples[samples["study_sample"].astype(str).isin(set(args.only_study_sample))].copy()
    if args.max_samples is not None:
        samples = samples.head(args.max_samples).copy()
    if samples.empty:
        raise ValueError("No samples selected for CellTypist annotation")
    return samples


def annotate_one_sample(
    row: pd.Series,
    integrated_obs: pd.DataFrame,
    target_index: pd.Index,
    model: str,
    sample_no: int,
    n_samples: int,
) -> tuple[pd.DataFrame, dict[str, object]]:
    start = time.time()
    study_sample = str(row["study_sample"])
    source = Path(str(row["source_path"]))
    print(f"SAMPLE {sample_no}/{n_samples} {study_sample} READ {source}", flush=True)

    expected = original_ids_for_sample(integrated_obs.loc[target_index], study_sample)
    adata = ad.read_h5ad(source)
    present = expected.intersection(adata.obs_names)
    missing = int(len(expected) - len(present))
    if missing:
        print(f"WARN {study_sample} missing={missing}", flush=True)
    adata = adata[present, :].copy()
    source_col = select_source_label_column(adata.obs)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    pred = celltypist.annotate(adata, model=model, majority_voting=False)
    labels = pred.predicted_labels["predicted_labels"].astype(str)
    probs = pred.probability_matrix.max(axis=1).astype(float)

    out = pd.DataFrame(index=adata.obs_names)
    out["cell_id"] = [f"{study_sample}__{idx}" for idx in out.index.astype(str)]
    out["dataset"] = str(row["dataset"])
    out["sample_id"] = str(row["sample_id"])
    out["study_sample"] = study_sample
    out["celltypist_liver_label"] = labels.reindex(out.index).to_numpy()
    out["celltypist_liver_confidence"] = probs.reindex(out.index).to_numpy()
    out["major_celltype"] = out["celltypist_liver_label"].map(liver_label_to_major)
    out["source_annotation_column"] = source_col or ""
    if source_col is None:
        out["source_author_label"] = ""
        out["source_author_major"] = "Unknown"
    else:
        source_labels = adata.obs[source_col].astype(str)
        out["source_author_label"] = source_labels.reindex(out.index).to_numpy()
        out["source_author_major"] = out["source_author_label"].map(source_label_to_major)
    out = out.set_index("cell_id", drop=True)

    summary = {
        "dataset": str(row["dataset"]),
        "sample_id": str(row["sample_id"]),
        "study_sample": study_sample,
        "source_h5ad": str(source.resolve()),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "missing_integrated_cells": missing,
        "source_annotation_column": source_col or "",
        "celltypist_model": model,
        "mean_confidence": float(out["celltypist_liver_confidence"].mean()),
        "top_major_celltype": str(out["major_celltype"].value_counts().idxmax()) if out.shape[0] else "",
        "elapsed_seconds": round(time.time() - start, 3),
    }
    print(
        f"DONE {study_sample} cells={adata.n_obs} top={summary['top_major_celltype']} "
        f"mean_conf={summary['mean_confidence']:.3f} elapsed={summary['elapsed_seconds']}",
        flush=True,
    )
    return out, summary


def top_counts(values: pd.Series, n: int = 3) -> str:
    counts = values.astype(str).value_counts().head(n)
    total = max(int(values.shape[0]), 1)
    return "; ".join([f"{idx}:{count}({count / total:.1%})" for idx, count in counts.items()])


def summarize_clusters(obs: pd.DataFrame, cluster_key: str) -> pd.DataFrame:
    rows = []
    for cluster, sub in obs.groupby(cluster_key, observed=True):
        n = int(sub.shape[0])
        major_counts = sub["major_celltype"].astype(str).value_counts()
        top_major = str(major_counts.idxmax()) if n else "Unknown"
        source = sub["source_author_major"].astype(str)
        source_known = source[~source.isin(["", "Unknown"])]
        rows.append(
            {
                cluster_key: str(cluster),
                "n_cells": n,
                "major_celltype": top_major,
                "major_celltype_fraction": float(major_counts.max() / n) if n else np.nan,
                "top3_celltypist_major": top_counts(sub["major_celltype"]),
                "top3_celltypist_liver_label": top_counts(sub["celltypist_liver_label"]),
                "mean_celltypist_confidence": float(sub["celltypist_liver_confidence"].mean()),
                "source_author_major_top3": top_counts(source_known) if not source_known.empty else "",
                "source_author_coverage": float(source_known.shape[0] / n) if n else 0.0,
                "predicted_doublet_rate": float(sub["predicted_doublet"].astype(bool).mean())
                if "predicted_doublet" in sub
                else np.nan,
                "cycling_rate": float(sub["cell_cycle_phase"].astype(str).isin(["S", "G2M"]).mean())
                if "cell_cycle_phase" in sub
                else np.nan,
            }
        )
    out = pd.DataFrame(rows)
    values = out[cluster_key].astype(str)
    if values.str.fullmatch(r"\d+").all():
        out = out.assign(_cluster_sort=values.astype(int)).sort_values("_cluster_sort").drop(columns="_cluster_sort")
    return out


def plot_umap(adata: ad.AnnData, figures_dir: Path, max_plot_cells: int, seed: int) -> list[str]:
    figures_dir.mkdir(parents=True, exist_ok=True)
    if "X_umap" not in adata.obsm:
        return []
    rng = np.random.default_rng(seed)
    idx = np.arange(adata.n_obs)
    if adata.n_obs > max_plot_cells:
        idx = np.sort(rng.choice(idx, size=max_plot_cells, replace=False))
    xy = np.asarray(adata.obsm["X_umap"])[idx]
    obs = adata.obs.iloc[idx]
    paths: list[str] = []

    labels = pd.Categorical(obs["major_celltype"].astype(str))
    cmap = plt.get_cmap("tab20", max(len(labels.categories), 1))
    fig, ax = plt.subplots(figsize=(8, 7), constrained_layout=True)
    for i, cat in enumerate(labels.categories):
        mask = labels == cat
        ax.scatter(xy[mask, 0], xy[mask, 1], s=0.7, c=[cmap(i % cmap.N)], label=cat, linewidths=0, alpha=0.65)
    ax.legend(frameon=False, markerscale=4, loc="center left", bbox_to_anchor=(1.0, 0.5))
    ax.set_title("CellTypist major cell types")
    ax.set_xticks([])
    ax.set_yticks([])
    path = figures_dir / "celltypist_major_umap.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    paths.append(str(path.resolve()))

    fig, ax = plt.subplots(figsize=(8, 7), constrained_layout=True)
    sca = ax.scatter(
        xy[:, 0],
        xy[:, 1],
        c=obs["celltypist_liver_confidence"].astype(float),
        s=0.7,
        cmap="viridis",
        linewidths=0,
    )
    fig.colorbar(sca, ax=ax, label="CellTypist confidence")
    ax.set_title("CellTypist confidence")
    ax.set_xticks([])
    ax.set_yticks([])
    path = figures_dir / "celltypist_confidence_umap.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    paths.append(str(path.resolve()))
    return paths


def main() -> int:
    args = parse_args()
    start = time.time()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    args.figures_dir.mkdir(parents=True, exist_ok=True)

    manifest = read_manifest(args.manifest)
    backed = ad.read_h5ad(args.input, backed="r")
    integrated_obs = backed.obs.copy()
    if args.exclude_cluster_key not in integrated_obs.columns:
        raise KeyError(f"{args.exclude_cluster_key!r} is not present in input obs")
    excluded = integrated_obs[args.exclude_cluster_key].astype(str).isin(set(map(str, args.exclude_cluster)))
    target_index = integrated_obs.index[~excluded]
    samples = build_sample_list(integrated_obs.loc[target_index], manifest, args)
    backed.file.close()
    print(
        f"INPUT cells={integrated_obs.shape[0]} target_cells={len(target_index)} "
        f"excluded={int(excluded.sum())} samples={samples.shape[0]}",
        flush=True,
    )

    cell_tables = []
    sample_summaries = []
    for i, (_, row) in enumerate(samples.iterrows(), start=1):
        per_cell, summary = annotate_one_sample(row, integrated_obs, target_index, args.model, i, samples.shape[0])
        cell_tables.append(per_cell)
        sample_summaries.append(summary)

    annotations = pd.concat(cell_tables, axis=0).reindex(target_index)
    missing = int(annotations["major_celltype"].isna().sum())
    if missing:
        print(f"WARN missing_annotations={missing}", flush=True)
    annotations["major_celltype"] = annotations["major_celltype"].fillna("Unknown")

    full_annotation = pd.DataFrame(index=integrated_obs.index)
    for col in annotations.columns:
        full_annotation[col] = annotations[col]
    full_annotation.loc[excluded, "celltypist_liver_label"] = "Doublet_suspect"
    full_annotation.loc[excluded, "celltypist_liver_confidence"] = np.nan
    full_annotation.loc[excluded, "major_celltype"] = "Doublet_suspect"
    full_annotation.loc[excluded, "source_annotation_column"] = ""
    full_annotation.loc[excluded, "source_author_label"] = ""
    full_annotation.loc[excluded, "source_author_major"] = "Doublet_suspect"
    full_annotation["excluded_doublet_cluster"] = excluded.to_numpy()

    obs_annotated = integrated_obs.join(
        full_annotation[
            [
                "celltypist_liver_label",
                "celltypist_liver_confidence",
                "major_celltype",
                "source_annotation_column",
                "source_author_label",
                "source_author_major",
                "excluded_doublet_cluster",
            ]
        ]
    )

    per_cell_path = args.metadata_dir / "celltypist_major_by_cell.tsv.gz"
    sample_summary_path = args.metadata_dir / "celltypist_major_by_sample.tsv"
    cluster_summary_path = args.metadata_dir / "celltypist_major_by_leiden.tsv"
    major_counts_path = args.metadata_dir / "celltypist_major_counts.tsv"
    full_annotation.to_csv(per_cell_path, sep="\t", compression="gzip")
    pd.DataFrame(sample_summaries).to_csv(sample_summary_path, sep="\t", index=False)
    cluster_summary = summarize_clusters(obs_annotated, args.cluster_key)
    cluster_summary.to_csv(cluster_summary_path, sep="\t", index=False)
    obs_annotated["major_celltype"].astype(str).value_counts().rename_axis("major_celltype").reset_index(
        name="n_cells"
    ).to_csv(major_counts_path, sep="\t", index=False)

    figures: list[str] = []
    if not args.skip_h5ad:
        print(f"WRITE {args.output}", flush=True)
        adata = ad.read_h5ad(args.input)
        for col in [
            "celltypist_liver_label",
            "celltypist_liver_confidence",
            "major_celltype",
            "source_annotation_column",
            "source_author_label",
            "source_author_major",
            "excluded_doublet_cluster",
        ]:
            adata.obs[col] = obs_annotated[col].to_numpy()
        adata.uns["major_celltype_annotation"] = {
            "method": "CellTypist",
            "celltypist_model": args.model,
            "excluded_doublet_cluster_key": args.exclude_cluster_key,
            "excluded_doublet_clusters": list(map(str, args.exclude_cluster)),
            "singleR_azimuth_status": "not_run_no_Rscript_or_installable_python_backend",
        }
        adata.write_h5ad(args.output, compression="gzip")
        figures = plot_umap(adata, args.figures_dir, args.max_plot_cells, args.seed)

    report = {
        "input": str(args.input.resolve()),
        "output": None if args.skip_h5ad else str(args.output.resolve()),
        "n_cells": int(integrated_obs.shape[0]),
        "target_cells": int(len(target_index)),
        "excluded_doublet_cluster_cells": int(excluded.sum()),
        "model": args.model,
        "n_samples": int(samples.shape[0]),
        "missing_annotations": missing,
        "per_cell_path": str(per_cell_path.resolve()),
        "sample_summary_path": str(sample_summary_path.resolve()),
        "cluster_summary_path": str(cluster_summary_path.resolve()),
        "major_counts_path": str(major_counts_path.resolve()),
        "figures": figures,
        "singleR_azimuth_status": "not_run_no_Rscript_or_installable_python_backend",
        "celltypist_version": version("celltypist"),
        "scanpy_version": version("scanpy"),
        "elapsed_seconds": round(time.time() - start, 3),
    }
    report_path = args.metadata_dir / "celltypist_major_report.json"
    with report_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)
    print(f"WROTE {report_path}", flush=True)
    print(json.dumps(report, indent=2, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
