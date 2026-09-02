from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd
from scipy.io import mmread


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert a Matrix Market counts bundle to h5ad.")
    parser.add_argument("--bundle-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    bundle = args.bundle_dir
    matrix_path = bundle / "counts.mtx"
    barcodes_path = bundle / "barcodes.tsv"
    features_path = bundle / "features.tsv"
    obs_path = bundle / "obs.tsv"

    counts = mmread(matrix_path).T.tocsr()
    barcodes = pd.read_csv(barcodes_path, header=None, sep="\t")[0].astype(str).to_list()
    features = pd.read_csv(features_path, header=None, sep="\t")[0].astype(str).to_list()
    obs = pd.read_csv(obs_path, sep="\t", index_col=0)
    obs.index = obs.index.astype(str)

    if counts.shape != (len(barcodes), len(features)):
        raise ValueError(
            f"Matrix shape {counts.shape} does not match barcodes/features "
            f"{len(barcodes)}x{len(features)}"
        )
    if list(obs.index) != barcodes:
        obs = obs.reindex(barcodes)
        if obs.isna().all(axis=None):
            raise ValueError("obs.tsv row names do not match barcodes.tsv")

    var = pd.DataFrame(index=pd.Index(features, name=None))
    var["gene_symbol"] = features
    adata = ad.AnnData(X=counts, obs=obs, var=var)
    adata.obs_names = barcodes
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output, compression="gzip")
    print(f"WROTE {args.output} cells={adata.n_obs} genes={adata.n_vars}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
