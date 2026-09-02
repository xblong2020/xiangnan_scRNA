#!/usr/bin/env python3
"""Audit an independent same-seed EGR1 scTenifoldKnk repeat."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

try:
    from figure3_egr1_common import PROJECT_ROOT, TARGET_TF, json_safe, write_json
except ModuleNotFoundError:
    from scripts.figure3_egr1_common import PROJECT_ROOT, TARGET_TF, json_safe, write_json


DEFAULT_ORIGINAL = (
    PROJECT_ROOT
    / "data/processed/driver/figure3e_egr1_sctenifoldknk/stressed_regenerative/results"
    / "figure3e_egr1_stressed_regenerative_seed15071990_perturbation_genes.tsv"
)
DEFAULT_REPEAT = (
    PROJECT_ROOT
    / "metadata/driver/figure3e_egr1_determinism"
    / "figure3e_egr1_stressed_regenerative_all_seed_perturbation_genes.tsv"
)
DEFAULT_OUT_DIR = PROJECT_ROOT / "metadata/driver/figure3e_egr1_determinism"
NUMERIC_COLUMNS = ["distance", "p.adj", "p.value", "Z", "FC"]


def run(original_path: Path, repeat_path: Path, out_dir: Path) -> dict:
    original = pd.read_csv(original_path, sep="\t")
    repeat = pd.read_csv(repeat_path, sep="\t")
    original = original.loc[original["tf"].astype(str).eq(TARGET_TF)].copy()
    repeat = repeat.loc[repeat["tf"].astype(str).eq(TARGET_TF)].copy()
    merged = original[["gene", *NUMERIC_COLUMNS]].merge(
        repeat[["gene", *NUMERIC_COLUMNS]],
        on="gene",
        how="outer",
        suffixes=("_original", "_repeat"),
        indicator=True,
    )
    rows = []
    all_close = merged["_merge"].eq("both").all()
    for column in NUMERIC_COLUMNS:
        left = pd.to_numeric(merged[f"{column}_original"], errors="coerce")
        right = pd.to_numeric(merged[f"{column}_repeat"], errors="coerce")
        finite = np.isfinite(left) & np.isfinite(right)
        difference = np.abs(left[finite] - right[finite])
        close = bool(
            finite.sum() == len(merged)
            and np.allclose(left[finite], right[finite], rtol=1e-10, atol=1e-12)
        )
        all_close = all_close and close
        rho = (
            float(stats.spearmanr(left[finite], right[finite]).statistic)
            if finite.sum() > 2
            else np.nan
        )
        rows.append(
            {
                "metric": column,
                "n_compared": int(finite.sum()),
                "allclose_rtol_1e_10_atol_1e_12": close,
                "max_absolute_difference": float(difference.max()) if len(difference) else np.nan,
                "mean_absolute_difference": float(difference.mean()) if len(difference) else np.nan,
                "spearman_rho": rho,
            }
        )
    audit = pd.DataFrame(rows)
    out_dir.mkdir(parents=True, exist_ok=True)
    table_path = out_dir / "figure3e_egr1_same_seed_determinism_audit.tsv"
    audit.to_csv(table_path, sep="\t", index=False)
    report = {
        "module": "Figure 3E same-seed determinism audit",
        "target_tf": TARGET_TF,
        "seed": 15071990,
        "independent_repeat": True,
        "parameters": {"nc_nNet": 10, "nc_nCells": 500, "nCores": 8},
        "original": str(original_path.resolve()),
        "repeat": str(repeat_path.resolve()),
        "gene_sets_identical": bool(merged["_merge"].eq("both").all()),
        "n_original": int(len(original)),
        "n_repeat": int(len(repeat)),
        "numeric_values_reproducible": bool(all_close),
        "status": "pass" if all_close else "warning",
        "review_risk_flags": []
        if all_close
        else [
            {
                "flag": "same_seed_not_numerically_identical",
                "severity": "review_attention",
                "detail": "Independent fixed-seed scTenifoldKnk output was not identical within rtol=1e-10 and atol=1e-12; see per-metric differences and rank correlations.",
            }
        ],
        "outputs": {"audit": str(table_path.resolve())},
    }
    write_json(json_safe(report), out_dir / "figure3e_egr1_same_seed_determinism_report.json")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original", type=Path, default=DEFAULT_ORIGINAL)
    parser.add_argument("--repeat", type=Path, default=DEFAULT_REPEAT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    report = run(args.original, args.repeat, args.out_dir)
    print(json.dumps(report, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
