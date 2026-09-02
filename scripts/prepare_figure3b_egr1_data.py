#!/usr/bin/env python3
"""Reuse the validated common 5,000-cell baseline exactly for Figure 3B."""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from figure3_egr1_common import PROJECT_ROOT, TARGET_TF, json_safe, sha256_file, write_json
except ModuleNotFoundError:
    from scripts.figure3_egr1_common import PROJECT_ROOT, TARGET_TF, json_safe, sha256_file, write_json


DEFAULT_SOURCE_DIR = PROJECT_ROOT / "metadata/driver/figure2b_sox4"
DEFAULT_OUT_DIR = PROJECT_ROOT / "metadata/driver/figure3b_egr1"


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer")


def copy_gzip_table(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(source, "rt", encoding="utf-8") as src, gzip.open(
        destination, "wt", encoding="utf-8", newline=""
    ) as dst:
        shutil.copyfileobj(src, dst)


def numeric_equivalence(
    component: str,
    source: pd.DataFrame,
    destination: pd.DataFrame,
    columns: list[str],
    source_path: Path,
    destination_path: Path,
) -> dict:
    missing = [column for column in columns if column not in source.columns or column not in destination.columns]
    if missing:
        return {
            "component": component,
            "status": "fail",
            "n_compared": 0,
            "exactly_equal": False,
            "max_absolute_difference": np.nan,
            "details": f"Missing columns: {missing}",
            "source_path": str(source_path.resolve()),
            "figure3_path": str(destination_path.resolve()),
        }
    left = source[columns].to_numpy()
    right = destination[columns].to_numpy()
    exact = left.shape == right.shape and np.array_equal(left, right, equal_nan=True)
    difference = np.asarray(left, dtype=float) - np.asarray(right, dtype=float)
    max_abs = float(np.nanmax(np.abs(difference))) if difference.size else 0.0
    return {
        "component": component,
        "status": "pass" if exact else "fail",
        "n_compared": int(left.size),
        "exactly_equal": bool(exact),
        "max_absolute_difference": max_abs,
        "details": ";".join(columns),
        "source_path": str(source_path.resolve()),
        "figure3_path": str(destination_path.resolve()),
    }


def categorical_equivalence(
    component: str,
    source: pd.Series,
    destination: pd.Series,
    source_path: Path,
    destination_path: Path,
) -> dict:
    exact = source.astype(str).tolist() == destination.astype(str).tolist()
    return {
        "component": component,
        "status": "pass" if exact else "fail",
        "n_compared": int(len(source)),
        "exactly_equal": bool(exact),
        "max_absolute_difference": np.nan,
        "details": "order-sensitive string comparison",
        "source_path": str(source_path.resolve()),
        "figure3_path": str(destination_path.resolve()),
    }


def run(source_dir: Path, out_dir: Path) -> dict:
    mapping = {
        "figure2b_sox4_plot_cells.tsv.gz": "figure3b_egr1_plot_cells.tsv.gz",
        "figure2b_sox4_baseline_grid_umap.tsv.gz": "figure3b_egr1_baseline_grid_umap.tsv.gz",
        "figure2b_sox4_baseline_grid_tsne.tsv.gz": "figure3b_egr1_baseline_grid_tsne.tsv.gz",
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    tables: dict[str, tuple[pd.DataFrame, pd.DataFrame, Path, Path]] = {}
    copies = []
    for source_name, destination_name in mapping.items():
        source_path = source_dir / source_name
        destination_path = out_dir / destination_name
        if not source_path.exists():
            raise FileNotFoundError(source_path)
        copy_gzip_table(source_path, destination_path)
        source = read_tsv(source_path)
        destination = read_tsv(destination_path)
        if not source.equals(destination):
            raise ValueError(f"Value-level copy failed for {source_name}")
        tables[source_name] = (source, destination, source_path, destination_path)
        copies.append(
            {
                "source": str(source_path.resolve()),
                "destination": str(destination_path.resolve()),
                "source_sha256": sha256_file(source_path),
                "destination_sha256": sha256_file(destination_path),
                "n_rows": int(len(source)),
                "n_columns": int(len(source.columns)),
                "value_level_equal": True,
            }
        )

    cells_key = "figure2b_sox4_plot_cells.tsv.gz"
    umap_key = "figure2b_sox4_baseline_grid_umap.tsv.gz"
    tsne_key = "figure2b_sox4_baseline_grid_tsne.tsv.gz"
    cells, cells_copy, cells_src, cells_dst = tables[cells_key]
    umap, umap_copy, umap_src, umap_dst = tables[umap_key]
    tsne, tsne_copy, tsne_src, tsne_dst = tables[tsne_key]

    audit_rows = [
        categorical_equivalence("cell_id_order", cells["cell_id"], cells_copy["cell_id"], cells_src, cells_dst),
        categorical_equivalence(
            "celloracle_state", cells["celloracle_state"], cells_copy["celloracle_state"], cells_src, cells_dst
        ),
        numeric_equivalence("umap_point_coordinates", cells, cells_copy, ["umap_1", "umap_2"], cells_src, cells_dst),
        numeric_equivalence("tsne_point_coordinates", cells, cells_copy, ["tsne_1", "tsne_2"], cells_src, cells_dst),
        numeric_equivalence("pseudotime", cells, cells_copy, ["pseudotime"], cells_src, cells_dst),
        numeric_equivalence("umap_grid_coordinates", umap, umap_copy, ["grid_x", "grid_y"], umap_src, umap_dst),
        numeric_equivalence(
            "umap_baseline_arrow_direction",
            umap,
            umap_copy,
            ["unit_x", "unit_y", "arrow_xend", "arrow_yend"],
            umap_src,
            umap_dst,
        ),
        categorical_equivalence("umap_density_mask", umap["keep"], umap_copy["keep"], umap_src, umap_dst),
        numeric_equivalence("tsne_grid_coordinates", tsne, tsne_copy, ["grid_x", "grid_y"], tsne_src, tsne_dst),
        numeric_equivalence(
            "tsne_baseline_arrow_direction",
            tsne,
            tsne_copy,
            ["unit_x", "unit_y", "arrow_xend", "arrow_yend"],
            tsne_src,
            tsne_dst,
        ),
        categorical_equivalence("tsne_density_mask", tsne["keep"], tsne_copy["keep"], tsne_src, tsne_dst),
    ]
    audit = pd.DataFrame(audit_rows)
    audit_path = out_dir / "figure3b_baseline_equivalence_audit.tsv"
    audit.to_csv(audit_path, sep="\t", index=False)
    report = {
        "module": "Figure 3B",
        "target_tf": TARGET_TF,
        "title": "Baseline developmental field",
        "reuse_policy": "Exact value-level reuse of the validated common Figure 2B baseline; t-SNE was not recalculated.",
        "source_contract": str(source_dir.resolve()),
        "n_cells": int(len(cells)),
        "cell_id_unique": bool(cells["cell_id"].is_unique),
        "all_equivalence_checks_pass": bool(audit["status"].eq("pass").all()),
        "parameters": {
            "n_grid": 20,
            "k_neighbors": 50,
            "density_quantile": 0.70,
            "seed": 15071990,
        },
        "copies": copies,
        "equivalence_audit": str(audit_path.resolve()),
        "outputs": {
            "plot_cells": str((out_dir / "figure3b_egr1_plot_cells.tsv.gz").resolve()),
            "umap_grid": str((out_dir / "figure3b_egr1_baseline_grid_umap.tsv.gz").resolve()),
            "tsne_grid": str((out_dir / "figure3b_egr1_baseline_grid_tsne.tsv.gz").resolve()),
        },
        "caveat": "The baseline field is TF-independent; EGR1 perturbation vectors are introduced only in Figure 3C.",
    }
    report_path = out_dir / "figure3b_egr1_data_report.json"
    write_json(json_safe(report), report_path)
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = run(args.source_dir, args.out_dir)
    print(json.dumps(report, ensure_ascii=False))
    return 0 if report["all_equivalence_checks_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())

