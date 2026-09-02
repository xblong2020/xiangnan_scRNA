from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
import sklearn
from sklearn.manifold import TSNE


ROOT = Path(__file__).resolve().parents[1]


def select_latent_rows(latent: pd.DataFrame, cell_ids: list[str]) -> pd.DataFrame:
    target_index = pd.Index(cell_ids, name="cell_id")
    if target_index.has_duplicates:
        raise ValueError("CytoTRACE2 score table contains duplicated cell_id values.")
    if latent.index.has_duplicates:
        raise ValueError("scVI latent table contains duplicated cell_id values.")
    missing = target_index.difference(latent.index)
    if not missing.empty:
        raise ValueError(f"Missing {len(missing)} CytoTRACE2 cells from scVI latent table.")
    return latent.loc[target_index].copy()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute a t-SNE embedding for the exact Figure 1G hepatocyte CytoTRACE2 cell set."
    )
    parser.add_argument(
        "--scores",
        type=Path,
        default=ROOT / "metadata/figure1c/figure1c_cytotrace2_scores_by_cell.hepatocyte.tsv.gz",
    )
    parser.add_argument(
        "--latent",
        type=Path,
        default=ROOT / "metadata/scvi/scvi_latent.tsv.gz",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "metadata/figure1/figure1g_hepatocyte_cytotrace2_tsne.tsv.gz",
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=ROOT / "metadata/figure1/figure1g_hepatocyte_cytotrace2_tsne_report.json",
    )
    parser.add_argument("--perplexity", type=float, default=50.0)
    parser.add_argument("--max-iter", type=int, default=500)
    parser.add_argument("--angle", type=float, default=0.7)
    parser.add_argument("--seed", type=int, default=20260710)
    parser.add_argument("--chunksize", type=int, default=50000)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.perplexity <= 0:
        raise ValueError("perplexity must be positive")
    if args.max_iter < 250:
        raise ValueError("max_iter must be at least 250 for sklearn t-SNE")
    if not 0 < args.angle <= 1:
        raise ValueError("angle must be within (0, 1]")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    scores = pd.read_csv(args.scores, sep="\t", usecols=["cell_id", "CytoTRACE2_Score"])
    scores = scores.loc[pd.to_numeric(scores["CytoTRACE2_Score"], errors="coerce").notna()].copy()
    cell_ids = scores["cell_id"].astype(str).tolist()
    selected_ids = set(cell_ids)

    pieces: list[pd.DataFrame] = []
    for chunk in pd.read_csv(args.latent, sep="\t", index_col=0, chunksize=args.chunksize):
        chunk.index = chunk.index.astype(str)
        selected_chunk = chunk.loc[chunk.index.isin(selected_ids)]
        if not selected_chunk.empty:
            pieces.append(selected_chunk)
    if not pieces:
        raise ValueError("No CytoTRACE2 cells were found in the scVI latent table.")

    latent = pd.concat(pieces, axis=0)
    selected_latent = select_latent_rows(latent, cell_ids)
    tsne = TSNE(
        n_components=2,
        perplexity=args.perplexity,
        init="pca",
        learning_rate="auto",
        metric="euclidean",
        random_state=args.seed,
        max_iter=args.max_iter,
        angle=args.angle,
        verbose=1,
    )
    coordinates = tsne.fit_transform(selected_latent.to_numpy(dtype=float))
    output = pd.DataFrame(
        {"cell_id": cell_ids, "TSNE_1": coordinates[:, 0], "TSNE_2": coordinates[:, 1]}
    )
    output.to_csv(args.output, sep="\t", index=False, compression="gzip")
    args.report.write_text(
        json.dumps(
            {
                "method": "sklearn.manifold.TSNE",
                "input_scores": str(args.scores),
                "input_latent": str(args.latent),
                "output": str(args.output),
                "n_cells": int(output.shape[0]),
                "n_latent_dimensions": int(selected_latent.shape[1]),
                "perplexity": args.perplexity,
                "max_iter": args.max_iter,
                "angle": args.angle,
                "seed": args.seed,
                "sklearn_version": sklearn.__version__,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
