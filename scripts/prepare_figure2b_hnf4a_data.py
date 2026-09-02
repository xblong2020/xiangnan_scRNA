#!/usr/bin/env python3
"""Reuse the validated SOX4 Figure 2B source data exactly for HNF4A Figure 2B."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer")


def copy_gzip_table(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(source, "rt", encoding="utf-8") as src, gzip.open(
        destination, "wt", encoding="utf-8", newline=""
    ) as dst:
        shutil.copyfileobj(src, dst)


def run(source_dir: Path, out_dir: Path) -> dict:
    mapping = {
        "figure2b_sox4_plot_cells.tsv.gz": "figure2b_hnf4a_plot_cells.tsv.gz",
        "figure2b_sox4_baseline_grid_umap.tsv.gz": "figure2b_hnf4a_baseline_grid_umap.tsv.gz",
        "figure2b_sox4_baseline_grid_tsne.tsv.gz": "figure2b_hnf4a_baseline_grid_tsne.tsv.gz",
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    audit_rows = []
    for source_name, destination_name in mapping.items():
        source = source_dir / source_name
        destination = out_dir / destination_name
        if not source.exists():
            raise FileNotFoundError(source)
        copy_gzip_table(source, destination)
        left = read_tsv(source)
        right = read_tsv(destination)
        equal = left.equals(right)
        max_abs = 0.0
        numeric = list(left.select_dtypes(include=[np.number]).columns)
        if numeric:
            max_abs = float(np.nanmax(np.abs(left[numeric].to_numpy() - right[numeric].to_numpy())))
        audit_rows.append({
            "component": source_name,
            "source_path": str(source.resolve()),
            "hnf4a_path": str(destination.resolve()),
            "n_rows": len(left),
            "n_columns": len(left.columns),
            "cell_id_set_equal": (
                set(left["cell_id"].astype(str)) == set(right["cell_id"].astype(str))
                if "cell_id" in left.columns else np.nan
            ),
            "values_exactly_equal": equal,
            "max_absolute_numeric_difference": max_abs,
            "source_sha256": sha256(source),
            "hnf4a_sha256": sha256(destination),
        })
    audit = pd.DataFrame(audit_rows)
    audit_path = out_dir / "figure2b_baseline_equivalence_audit.tsv"
    audit.to_csv(audit_path, sep="\t", index=False)
    cells = read_tsv(out_dir / "figure2b_hnf4a_plot_cells.tsv.gz")
    report = {
        "module": "Figure 2B",
        "target_tf": "HNF4A",
        "analysis": "Baseline developmental field",
        "reuse_policy": "Exact value-level reuse of validated SOX4 baseline cell, UMAP, t-SNE, pseudotime and grid tables",
        "n_cells": int(len(cells)),
        "cell_id_unique": bool(cells["cell_id"].is_unique),
        "all_components_exactly_equal": bool(audit["values_exactly_equal"].all()),
        "equivalence_audit": str(audit_path.resolve()),
        "note": "Figure 2B is TF-independent; reuse prevents coordinate drift from re-running t-SNE.",
    }
    (out_dir / "figure2b_hnf4a_data_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=PROJECT_ROOT / "metadata/driver/figure2b_sox4")
    parser.add_argument("--out-dir", type=Path, default=PROJECT_ROOT / "metadata/driver/figure2b_hnf4a")
    args = parser.parse_args()
    print(json.dumps(run(args.source_dir, args.out_dir), ensure_ascii=False))


if __name__ == "__main__":
    main()
