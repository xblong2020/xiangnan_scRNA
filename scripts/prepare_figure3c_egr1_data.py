#!/usr/bin/env python3
"""Match true saved EGR1 CellOracle displacement to the exact Figure 3B cells."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from figure3_egr1_common import PROJECT_ROOT, TARGET_TF, json_safe, write_json
except ModuleNotFoundError:
    from scripts.figure3_egr1_common import PROJECT_ROOT, TARGET_TF, json_safe, write_json


DEFAULT_BASELINE_DIR = PROJECT_ROOT / "metadata/driver/figure3b_egr1"
DEFAULT_PERTURBATION = PROJECT_ROOT / "metadata/driver/celloracle_module6_8_cell_shift_summary.tsv.gz"
DEFAULT_SIMULATION_REPORT = PROJECT_ROOT / "metadata/driver/celloracle_module6_8_perturbation_report.json"
DEFAULT_OUT_DIR = PROJECT_ROOT / "metadata/driver/figure3c_egr1"


def run(
    baseline_dir: Path,
    perturbation_path: Path,
    simulation_report_path: Path,
    out_dir: Path,
) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    cells_path = baseline_dir / "figure3b_egr1_plot_cells.tsv.gz"
    if not cells_path.exists():
        raise FileNotFoundError(f"{cells_path} is missing; run Figure 3B preparation first.")
    cells = pd.read_csv(cells_path, sep="\t")
    perturb = pd.read_csv(perturbation_path, sep="\t")
    perturb = perturb.loc[perturb["tf"].astype(str).eq(TARGET_TF)].copy()
    required = {
        "tf",
        "cell_id",
        "delta_embedding_1",
        "delta_embedding_2",
        "embedding_shift_norm",
        "malignant_axis_projection",
    }
    missing = sorted(required.difference(perturb.columns))
    if missing:
        raise ValueError(f"Missing EGR1 perturbation columns: {missing}")
    if perturb["cell_id"].astype(str).duplicated().any():
        raise ValueError("EGR1 perturbation source contains duplicate cell_id rows")
    matched = cells.merge(
        perturb[
            [
                "cell_id",
                "tf",
                "delta_embedding_1",
                "delta_embedding_2",
                "embedding_shift_norm",
                "malignant_axis_projection",
                "mean_abs_delta_x",
                "mean_delta_x",
            ]
        ],
        on="cell_id",
        how="left",
        validate="one_to_one",
    )
    delta = matched[["delta_embedding_1", "delta_embedding_2"]].to_numpy(dtype=float)
    if len(matched) != 5000 or not np.isfinite(delta).all():
        missing_n = int(np.sum(~np.isfinite(delta).all(axis=1)))
        raise ValueError(f"Figure 3C requires 5,000 finite EGR1 vectors; invalid rows={missing_n}")
    if not matched["tf"].astype(str).eq(TARGET_TF).all():
        raise ValueError("Matched Figure 3C data contain a non-EGR1 perturbation row")

    simulation_report = json.loads(simulation_report_path.read_text(encoding="utf-8-sig"))
    parameters = simulation_report.get("parameters", {})
    egr1_run = next(
        (
            row
            for row in simulation_report.get("result", {}).get("per_tf_reports", [])
            if row.get("tf") == TARGET_TF
        ),
        None,
    )
    if not egr1_run or egr1_run.get("condition") != {TARGET_TF: 0.0}:
        raise ValueError("Module 6.8 report does not document the expected EGR1 = 0 condition")
    matched_path = out_dir / "figure3c_egr1_matched_cells.tsv.gz"
    matched.to_csv(matched_path, sep="\t", index=False, compression="gzip")
    report = {
        "module": "Figure 3C data preparation",
        "target_tf": TARGET_TF,
        "condition": {TARGET_TF: 0.0},
        "source_cells": str(cells_path.resolve()),
        "source_perturbation": str(perturbation_path.resolve()),
        "source_simulation_report": str(simulation_report_path.resolve()),
        "source_oracle": simulation_report.get("result", {}).get("input_oracle"),
        "celloracle_version": simulation_report.get("result", {}).get("celloracle_version"),
        "simulation_parameters": parameters,
        "n_source_egr1_cells": int(len(perturb)),
        "n_matched_cells": int(len(matched)),
        "cell_filter": "exact Figure 3B strict-main 5,000-cell subset in the same order",
        "delta_embedding_finite": True,
        "delta_embedding_source": "saved CellOracle Module 6.8 EGR1 knockout simulation",
        "embedding_shift_norm_range": [
            float(matched["embedding_shift_norm"].min()),
            float(matched["embedding_shift_norm"].max()),
        ],
        "output": str(matched_path.resolve()),
        "caveat": "UMAP vectors are native saved CellOracle displacement. t-SNE displacement is generated only as a supplementary local projection.",
    }
    report_path = out_dir / "figure3c_egr1_data_report.json"
    write_json(json_safe(report), report_path)
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-dir", type=Path, default=DEFAULT_BASELINE_DIR)
    parser.add_argument("--perturbation", type=Path, default=DEFAULT_PERTURBATION)
    parser.add_argument("--simulation-report", type=Path, default=DEFAULT_SIMULATION_REPORT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = run(args.baseline_dir, args.perturbation, args.simulation_report, args.out_dir)
    print(json.dumps(report, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

