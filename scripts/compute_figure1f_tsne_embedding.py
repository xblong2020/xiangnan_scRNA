from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from sklearn.manifold import TSNE


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute t-SNE embedding for Figure1F from existing module5 Slingshot cells.")
    parser.add_argument(
        "--pseudotime",
        type=Path,
        default=ROOT / "metadata/trajectory/trajectory_module5_3_main_strict_pseudotime_merged.tsv.gz",
    )
    parser.add_argument(
        "--scanvi-embedding",
        type=Path,
        default=ROOT / "data/processed/trajectory/module5_3/main_strict/embedding_x_scanvi.tsv.gz",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "metadata/trajectory/figure1f_main_strict_tsne.tsv.gz",
    )
    parser.add_argument("--perplexity", type=float, default=50.0)
    parser.add_argument("--seed", type=int, default=20260709)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    pt = pd.read_csv(args.pseudotime, sep="\t", usecols=["cell_id"])
    emb = pd.read_csv(args.scanvi_embedding, sep="\t")
    merged = pt.merge(emb, on="cell_id", how="left", validate="one_to_one")
    feature_cols = [c for c in merged.columns if c != "cell_id"]
    X = merged[feature_cols].to_numpy(dtype=float)

    tsne = TSNE(
        n_components=2,
        perplexity=args.perplexity,
        init="pca",
        learning_rate="auto",
        random_state=args.seed,
        metric="euclidean",
    )
    coords = tsne.fit_transform(X)
    out = pd.DataFrame({"cell_id": merged["cell_id"], "TSNE_1": coords[:, 0], "TSNE_2": coords[:, 1]})
    out.to_csv(args.output, sep="\t", index=False, compression="gzip")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
