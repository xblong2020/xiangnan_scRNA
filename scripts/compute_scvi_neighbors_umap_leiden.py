from __future__ import annotations

import argparse
import json
import time
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute neighbors, UMAP, and Leiden clusters from adata.obsm['X_scVI']."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=ROOT / "data/processed/scvi/scvi_integrated_counts_hvg.h5ad",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "data/processed/scvi/scvi_integrated_counts_hvg.scvi_neighbors_umap_leiden.h5ad",
    )
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/scvi")
    parser.add_argument("--n-neighbors", type=int, default=30)
    parser.add_argument("--min-dist", type=float, default=0.3)
    parser.add_argument("--resolution", type=float, default=1.0)
    parser.add_argument("--cluster-key", default="leiden_scvi")
    parser.add_argument("--seed", type=int, default=20260601)
    return parser.parse_args()


def export_tables(adata: ad.AnnData, metadata_dir: Path, cluster_key: str) -> dict[str, str]:
    metadata_dir.mkdir(parents=True, exist_ok=True)
    umap_path = metadata_dir / "scvi_umap.tsv.gz"
    clusters_path = metadata_dir / "scvi_leiden_clusters.tsv.gz"
    counts_path = metadata_dir / "scvi_leiden_cluster_counts.tsv"

    pd.DataFrame(
        adata.obsm["X_umap"],
        index=adata.obs_names,
        columns=["UMAP_1", "UMAP_2"],
    ).to_csv(umap_path, sep="\t", compression="gzip")

    obs_cols = [cluster_key]
    for col in ["dataset", "sample_id", "study_sample"]:
        if col in adata.obs.columns:
            obs_cols.append(col)
    adata.obs[obs_cols].to_csv(clusters_path, sep="\t", compression="gzip")

    counts = (
        adata.obs[cluster_key]
        .astype(str)
        .value_counts()
        .rename_axis(cluster_key)
        .reset_index(name="n_cells")
        .sort_values(cluster_key, key=lambda s: s.astype(int) if s.str.fullmatch(r"\d+").all() else s)
    )
    counts.to_csv(counts_path, sep="\t", index=False)

    return {
        "umap_path": str(umap_path.resolve()),
        "clusters_path": str(clusters_path.resolve()),
        "cluster_counts_path": str(counts_path.resolve()),
    }


def main() -> int:
    args = parse_args()
    start = time.time()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)

    print(f"READ {args.input}", flush=True)
    adata = ad.read_h5ad(args.input)
    if "X_scVI" not in adata.obsm:
        raise KeyError("adata.obsm['X_scVI'] is required but was not found")

    print(
        f"INPUT cells={adata.n_obs} genes={adata.n_vars} X_scVI={adata.obsm['X_scVI'].shape}",
        flush=True,
    )
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

    print(f"WRITE {args.output}", flush=True)
    adata.write_h5ad(args.output, compression="gzip")
    exported = export_tables(adata, args.metadata_dir, args.cluster_key)

    report = {
        "input": str(args.input.resolve()),
        "output": str(args.output.resolve()),
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "use_rep": "X_scVI",
        "n_neighbors": int(args.n_neighbors),
        "min_dist": float(args.min_dist),
        "resolution": float(args.resolution),
        "cluster_key": args.cluster_key,
        "n_clusters": int(adata.obs[args.cluster_key].nunique()),
        "seed": int(args.seed),
        "scanpy_version": version("scanpy"),
        "elapsed_seconds": round(time.time() - start, 3),
    }
    report.update(exported)
    report_path = args.metadata_dir / "scvi_neighbors_umap_leiden_report.json"
    with report_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)
    print(f"WROTE {report_path}", flush=True)
    print(json.dumps(report, indent=2, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
