from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scvi
import torch
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run global scVI integration from a counts manifest.")
    parser.add_argument(
        "--manifest",
        type=Path,
        default=ROOT / "metadata/scvi/scvi_input_manifest.counts.tsv",
    )
    parser.add_argument("--out-dir", type=Path, default=ROOT / "data/processed/scvi")
    parser.add_argument("--n-top-genes", type=int, default=2000)
    parser.add_argument("--min-detected-cells", type=int, default=100)
    parser.add_argument("--max-epochs", type=int, default=20)
    parser.add_argument("--batch-size", type=int, default=1024)
    parser.add_argument("--n-latent", type=int, default=30)
    parser.add_argument("--n-hidden", type=int, default=128)
    parser.add_argument("--n-layers", type=int, default=2)
    parser.add_argument("--batch-key", default="dataset")
    parser.add_argument("--seed", type=int, default=20260601)
    return parser.parse_args()


def as_csr(x) -> sparse.csr_matrix:
    if sparse.issparse(x):
        return x.tocsr()
    return sparse.csr_matrix(x)


def read_manifest(path: Path) -> pd.DataFrame:
    manifest = pd.read_csv(path, sep="\t")
    include = manifest["include_in_scvi"].astype(str).str.lower().eq("true")
    selected = manifest.loc[include].copy()
    if selected.empty:
        raise ValueError("No samples marked include_in_scvi=True")
    selected["output"] = selected["output"].map(lambda p: str(Path(p)))
    return selected


def find_common_genes(manifest: pd.DataFrame) -> list[str]:
    common: set[str] | None = None
    for _, row in manifest.iterrows():
        path = Path(row["output"])
        adata = ad.read_h5ad(path, backed="r")
        genes = set(map(str, adata.var_names))
        adata.file.close()
        common = genes if common is None else common.intersection(genes)
        print(f"COMMON {row['dataset']}:{row['label']} genes={len(genes)} current_common={len(common)}", flush=True)
    if not common:
        raise ValueError("No common genes across included datasets")
    return sorted(common)


def select_variable_genes(
    manifest: pd.DataFrame,
    common_genes: list[str],
    n_top_genes: int,
    min_detected_cells: int,
    out_path: Path,
) -> list[str]:
    gene_index = pd.Index(common_genes)
    n_genes = len(common_genes)
    total_cells = 0
    sum_counts = np.zeros(n_genes, dtype=np.float64)
    sum_squares = np.zeros(n_genes, dtype=np.float64)
    detected_cells = np.zeros(n_genes, dtype=np.int64)

    for _, row in manifest.iterrows():
        path = Path(row["output"])
        adata = ad.read_h5ad(path)
        loc = adata.var_names.get_indexer(gene_index)
        if np.any(loc < 0):
            missing = gene_index[np.where(loc < 0)[0][:5]].to_list()
            raise ValueError(f"Missing common genes in {path}: {missing}")
        x = as_csr(adata.X[:, loc])
        total_cells += x.shape[0]
        sum_counts += np.asarray(x.sum(axis=0)).ravel()
        sum_squares += np.asarray(x.multiply(x).sum(axis=0)).ravel()
        detected_cells += np.asarray((x > 0).sum(axis=0)).ravel()
        print(f"STATS {row['dataset']}:{row['label']} cells={x.shape[0]}", flush=True)

    mean = sum_counts / total_cells
    variance = (sum_squares / total_cells) - np.square(mean)
    variance = np.maximum(variance, 0)
    dispersion = variance / np.maximum(mean, 1e-8)
    genes = np.asarray(common_genes, dtype=str)
    technical = np.char.startswith(genes, "MT-")
    valid = (detected_cells >= min_detected_cells) & (mean > 0) & (~technical)
    if valid.sum() < n_top_genes:
        raise ValueError(f"Only {valid.sum()} valid genes available for n_top_genes={n_top_genes}")

    order = np.argsort(dispersion[valid])[::-1]
    valid_genes = genes[valid]
    selected = valid_genes[order[:n_top_genes]].tolist()

    stats = pd.DataFrame(
        {
            "gene": genes,
            "mean": mean,
            "variance": variance,
            "dispersion": dispersion,
            "detected_cells": detected_cells,
            "technical_excluded": technical,
            "selected_hvg": np.isin(genes, selected),
        }
    ).sort_values(["selected_hvg", "dispersion"], ascending=[False, False])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    stats.to_csv(out_path, sep="\t", index=False)
    print(f"WROTE {out_path}", flush=True)
    return selected


def build_adata(manifest: pd.DataFrame, genes: list[str]) -> ad.AnnData:
    gene_index = pd.Index(genes)
    parts: list[ad.AnnData] = []
    for _, row in manifest.iterrows():
        path = Path(row["output"])
        one = ad.read_h5ad(path)
        loc = one.var_names.get_indexer(gene_index)
        x = as_csr(one.X[:, loc]).astype(np.float32)
        obs = one.obs.copy()
        obs["dataset"] = str(row["dataset"])
        obs["sample_id"] = str(row["label"])
        obs["study_sample"] = f"{row['dataset']}__{row['label']}"
        obs["source_h5ad"] = str(path)
        obs.index = [f"{row['dataset']}__{row['label']}__{idx}" for idx in obs.index.astype(str)]
        var = pd.DataFrame(index=gene_index)
        part = ad.AnnData(X=x, obs=obs, var=var)
        part.layers["counts"] = x.copy()
        parts.append(part)
        print(f"LOAD {row['dataset']}:{row['label']} cells={part.n_obs} genes={part.n_vars}", flush=True)

    combined = ad.concat(parts, axis=0, join="inner", merge="same", index_unique=None)
    combined.obs_names_make_unique()
    combined.var_names_make_unique()
    for col in ["dataset", "sample_id", "study_sample"]:
        combined.obs[col] = combined.obs[col].astype("category")
    return combined


def history_tail(model: scvi.model.SCVI) -> dict[str, float | int | str]:
    out: dict[str, float | int | str] = {}
    history = getattr(model, "history", None)
    if history is None:
        return out
    for key in history.keys():
        values = history[key]
        try:
            last = values.dropna().iloc[-1]
        except Exception:
            continue
        if hasattr(last, "iloc"):
            last = last.iloc[0]
        try:
            out[f"final_{key}"] = float(last)
        except Exception:
            out[f"final_{key}"] = str(last)
    try:
        out["epochs_recorded"] = int(max(len(v) for v in history.values()))
    except Exception:
        pass
    return out


def main() -> int:
    args = parse_args()
    start = time.time()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    model_dir = args.out_dir / "model_scvi_global_counts"
    hvg_path = ROOT / "metadata/scvi/scvi_hvg_genes.tsv"
    output_h5ad = args.out_dir / "scvi_integrated_counts_hvg.h5ad"
    latent_path = ROOT / "metadata/scvi/scvi_latent.tsv.gz"
    obs_path = ROOT / "metadata/scvi/scvi_obs.tsv.gz"
    report_path = ROOT / "metadata/scvi/scvi_training_report.json"

    scvi.settings.seed = args.seed
    if torch.cuda.is_available():
        torch.set_float32_matmul_precision("high")
    manifest = read_manifest(args.manifest)
    common_genes = find_common_genes(manifest)
    selected_genes = select_variable_genes(
        manifest,
        common_genes,
        args.n_top_genes,
        args.min_detected_cells,
        hvg_path,
    )
    adata = build_adata(manifest, selected_genes)
    print(f"COMBINED cells={adata.n_obs} genes={adata.n_vars}", flush=True)

    if args.batch_key not in adata.obs.columns:
        raise ValueError(f"batch key {args.batch_key!r} is not present in adata.obs")

    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=args.batch_key)
    model = scvi.model.SCVI(
        adata,
        n_latent=args.n_latent,
        n_hidden=args.n_hidden,
        n_layers=args.n_layers,
        gene_likelihood="nb",
    )
    accelerator = "gpu" if torch.cuda.is_available() else "cpu"
    print(f"TRAIN accelerator={accelerator} max_epochs={args.max_epochs} batch_size={args.batch_size}", flush=True)
    model.train(
        max_epochs=args.max_epochs,
        batch_size=args.batch_size,
        accelerator=accelerator,
        devices=1,
        early_stopping=True,
    )

    latent = model.get_latent_representation(batch_size=args.batch_size)
    adata.obsm["X_scVI"] = latent
    model.save(model_dir, overwrite=True)
    adata.write_h5ad(output_h5ad, compression="gzip")

    latent_df = pd.DataFrame(
        latent,
        index=adata.obs_names,
        columns=[f"SCVI_{i + 1}" for i in range(latent.shape[1])],
    )
    latent_df.to_csv(latent_path, sep="\t", compression="gzip")
    adata.obs.to_csv(obs_path, sep="\t", compression="gzip")

    report = {
        "manifest": str(args.manifest.resolve()),
        "output_h5ad": str(output_h5ad.resolve()),
        "model_dir": str(model_dir.resolve()),
        "latent_path": str(latent_path.resolve()),
        "obs_path": str(obs_path.resolve()),
        "hvg_path": str(hvg_path.resolve()),
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "n_common_genes": int(len(common_genes)),
        "n_top_genes": int(args.n_top_genes),
        "batch_key": args.batch_key,
        "max_epochs": int(args.max_epochs),
        "batch_size": int(args.batch_size),
        "accelerator": accelerator,
        "torch_version": torch.__version__,
        "scvi_version": scvi.__version__,
        "elapsed_seconds": round(time.time() - start, 3),
    }
    report.update(history_tail(model))
    with report_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, ensure_ascii=False, indent=2)
    print(f"WROTE {output_h5ad}", flush=True)
    print(f"WROTE {latent_path}", flush=True)
    print(f"WROTE {obs_path}", flush=True)
    print(f"WROTE {model_dir}", flush=True)
    print(f"WROTE {report_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
