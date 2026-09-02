#!/usr/bin/env python3
"""Match real HNF4A CellOracle virtual-knockout displacement to Figure 2B cells."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def run(
    baseline_dir: Path,
    perturbation_path: Path,
    score_path: Path,
    perturbation_report_path: Path,
    out_dir: Path,
    target_tf: str,
) -> dict:
    cells_path = baseline_dir / "figure2b_hnf4a_plot_cells.tsv.gz"
    cells = pd.read_csv(cells_path, sep="\t")
    perturb = pd.read_csv(perturbation_path, sep="\t")
    perturb = perturb.loc[perturb["tf"].astype(str).eq(target_tf)].copy()
    required = {
        "cell_id", "delta_embedding_1", "delta_embedding_2",
        "embedding_shift_norm", "malignant_axis_projection",
    }
    missing = sorted(required.difference(perturb.columns))
    if missing:
        raise ValueError(f"Missing HNF4A perturbation columns: {missing}")
    if len(perturb) == 0:
        raise ValueError(f"No CellOracle rows found for {target_tf}")
    merged = cells.merge(
        perturb[["cell_id", "delta_embedding_1", "delta_embedding_2",
                 "embedding_shift_norm", "malignant_axis_projection"]],
        on="cell_id", how="left", validate="one_to_one",
    )
    if merged[["delta_embedding_1", "delta_embedding_2"]].isna().any().any():
        raise ValueError("Some Figure 2B cells lack HNF4A delta_embedding")
    score = pd.read_csv(score_path, sep="\t")
    score = score.loc[score["tf"].astype(str).eq(target_tf)].copy()
    score_cols = [
        "cell_id", "inner_product_score", "inner_product_cell_score",
        "inner_product_raw_pseudotime",
    ]
    if sorted(set(score_cols).difference(score.columns)):
        raise ValueError("HNF4A quantitative score table lacks required columns")
    merged = merged.merge(score[score_cols], on="cell_id", how="left", validate="one_to_one")
    if merged[score_cols[1:]].isna().any().any():
        raise ValueError("Some matched HNF4A cells lack project quantitative scores")
    if not np.isfinite(merged[["delta_embedding_1", "delta_embedding_2"]]).all().all():
        raise ValueError("HNF4A delta_embedding contains non-finite values")

    with perturbation_report_path.open(encoding="utf-8-sig") as handle:
        module6_report = json.load(handle)
    params = module6_report.get("parameters", {})
    result = module6_report.get("result", {})
    if target_tf not in result.get("perturbed_tfs", []):
        raise ValueError("Module 6.8 report does not record HNF4A as a simulated TF")
    out_dir.mkdir(parents=True, exist_ok=True)
    matched_path = out_dir / "figure2c_hnf4a_matched_cells.tsv.gz"
    merged.to_csv(matched_path, sep="\t", index=False, compression="gzip")
    report = {
        "module": "Figure 2C",
        "target_tf": target_tf,
        "analysis": "HNF4A virtual knockout predicted perturbation",
        "source_cells": str(cells_path.resolve()),
        "delta_embedding_source": str(perturbation_path.resolve()),
        "celloracle_object_path": result.get("input_oracle"),
        "source_celloracle_report": str(perturbation_report_path.resolve()),
        "n_propagation": params.get("n_propagation"),
        "seed": params.get("seed"),
        "n_source_perturbation_cells": int(len(perturb)),
        "n_cells": int(len(merged)),
        "score_range_cell_level_project": [
            float(merged["inner_product_raw_pseudotime"].min()),
            float(merged["inner_product_raw_pseudotime"].max()),
        ],
        "output": str(matched_path.resolve()),
        "umap_method": "Direct saved CellOracle delta_embedding after HNF4A = 0",
        "tsne_method": "Supplementary local UMAP-to-t-SNE Jacobian projection; not a native CellOracle t-SNE simulation",
        "caveat": "Predicted perturbation from virtual knockout; computationally inferred state shift is network perturbation evidence.",
    }
    (out_dir / "figure2c_hnf4a_data_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-tf", default="HNF4A")
    parser.add_argument("--baseline-dir", type=Path, default=PROJECT_ROOT / "metadata/driver/figure2b_hnf4a")
    parser.add_argument("--perturbation", type=Path, default=PROJECT_ROOT / "metadata/driver/celloracle_module6_8_cell_shift_summary.tsv.gz")
    parser.add_argument("--score", type=Path, default=PROJECT_ROOT / "metadata/driver/celloracle_module6_9b_cell_level_scores.tsv.gz")
    parser.add_argument("--perturbation-report", type=Path, default=PROJECT_ROOT / "metadata/driver/celloracle_module6_8_perturbation_report.json")
    parser.add_argument("--out-dir", type=Path, default=PROJECT_ROOT / "metadata/driver/figure2c_hnf4a")
    args = parser.parse_args()
    print(json.dumps(run(args.baseline_dir, args.perturbation, args.score, args.perturbation_report,
                         args.out_dir, args.target_tf), ensure_ascii=False))


if __name__ == "__main__":
    main()
