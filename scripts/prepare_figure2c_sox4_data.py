#!/usr/bin/env python3
"""Prepare matched SOX4 perturbation data for the Figure 2C R plot."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_B_DATA = PROJECT_ROOT / "metadata/driver/figure2b_sox4"
DEFAULT_PERTURB = PROJECT_ROOT / "metadata/driver/celloracle_module6_8_cell_shift_summary.tsv.gz"
DEFAULT_SCORE = PROJECT_ROOT / "metadata/driver/celloracle_module6_9b_cell_level_scores.tsv.gz"
DEFAULT_OUT = PROJECT_ROOT / "metadata/driver/figure2c_sox4"


def run(b_data_dir: Path, perturbation_path: Path, score_path: Path, out_dir: Path) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    cells_path = b_data_dir / "figure2b_sox4_plot_cells.tsv.gz"
    if not cells_path.exists():
        raise FileNotFoundError(
            f"{cells_path} is missing; run prepare_figure2b_sox4_data.py and plot_figure2b_sox4_baseline.R first."
        )
    cells = pd.read_csv(cells_path, sep="\t", compression="gzip")
    perturb = pd.read_csv(perturbation_path, sep="\t", compression="gzip")
    perturb = perturb.loc[perturb["tf"].astype(str).eq("SOX4")].copy()
    required = {"cell_id", "delta_embedding_1", "delta_embedding_2", "embedding_shift_norm"}
    missing = sorted(required.difference(perturb.columns))
    if missing:
        raise ValueError(f"Missing SOX4 perturbation columns: {missing}")
    out = cells.merge(
        perturb[["cell_id", "delta_embedding_1", "delta_embedding_2", "embedding_shift_norm", "malignant_axis_projection"]],
        on="cell_id",
        how="left",
        validate="one_to_one",
    )
    if out[["delta_embedding_1", "delta_embedding_2"]].isna().any().any():
        missing_n = int(out["delta_embedding_1"].isna().sum())
        raise ValueError(f"{missing_n} Figure 2B cells lack SOX4 perturbation vectors")
    score_available = score_path.exists()
    if score_available:
        score = pd.read_csv(score_path, sep="\t", compression="gzip")
        score = score.loc[score["tf"].astype(str).eq("SOX4")].copy()
        score_cols = ["cell_id", "inner_product_score", "inner_product_cell_score", "inner_product_raw_pseudotime"]
        missing_score = sorted(set(score_cols).difference(score.columns))
        if missing_score:
            raise ValueError(f"Missing inner-product columns: {missing_score}")
        out = out.merge(score[score_cols], on="cell_id", how="left", validate="one_to_one")
        if out["inner_product_score"].isna().any():
            raise ValueError("Some matched cells lack project inner_product_score")
    out.to_csv(out_dir / "figure2c_sox4_matched_cells.tsv.gz", sep="\t", index=False, compression="gzip")
    report = {
        "module": "Figure 2C",
        "target": "SOX4 knockout perturbation vector field",
        "source_cells": str(cells_path.resolve()),
        "source_perturbation": str(perturbation_path.resolve()),
        "tf": "SOX4",
        "n_source_perturbation_cells": int(len(perturb)),
        "n_matched_cells": int(len(out)),
        "cell_filter": "the exact 5,000-cell Figure 2B strict-main subset",
        "score_source": str(score_path.resolve()) if score_available else None,
        "score_available": score_available,
        "outputs": {"matched_cells": str((out_dir / "figure2c_sox4_matched_cells.tsv.gz").resolve())},
        "note": "UMAP delta_embedding is the saved CellOracle perturbation vector. t-SNE arrows are a local UMAP-to-t-SNE projection for supplementary visualization.",
    }
    (out_dir / "figure2c_sox4_data_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figure2b-data-dir", type=Path, default=DEFAULT_B_DATA)
    parser.add_argument("--perturbation", type=Path, default=DEFAULT_PERTURB)
    parser.add_argument("--score", type=Path, default=DEFAULT_SCORE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    print(json.dumps(run(args.figure2b_data_dir, args.perturbation, args.score, args.out_dir), ensure_ascii=False))


if __name__ == "__main__":
    main()
