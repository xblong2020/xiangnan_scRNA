#!/usr/bin/env python3
"""Audit Figure 2 HNF4A, Figure 3 EGR1, and Figure 4 SOX4 source contracts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from figure3_egr1_common import PROJECT_ROOT, json_safe, write_json
except ModuleNotFoundError:
    from scripts.figure3_egr1_common import PROJECT_ROOT, json_safe, write_json


DEFAULT_OUT_DIR = PROJECT_ROOT / "metadata/driver/three_axis_figure_consistency"


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer")


def add(rows: list[dict], check: str, axis: str, passed: bool, observed, expected, details: str, source: Path) -> None:
    rows.append(
        {
            "check": check,
            "axis": axis,
            "status": "pass" if passed else "fail",
            "observed": json.dumps(json_safe(observed), ensure_ascii=False)
            if isinstance(observed, (dict, list))
            else observed,
            "expected": expected,
            "details": details,
            "source": str(source.resolve()),
        }
    )


def valid_grid_scores(path: Path) -> np.ndarray:
    frame = read_table(path)
    values = pd.to_numeric(frame["inner_product_score_grid"], errors="coerce")
    valid = np.isfinite(values)
    if "keep_score" in frame:
        valid &= frame["keep_score"].astype(str).str.lower().isin({"true", "1"})
    return values.loc[valid].to_numpy(dtype=float)


def run(out_dir: Path) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    canonical_cells_path = PROJECT_ROOT / "metadata/driver/figure2b_sox4/figure2b_sox4_plot_cells.tsv.gz"
    canonical_umap_path = PROJECT_ROOT / "metadata/driver/figure2b_sox4/figure2b_sox4_baseline_grid_umap.tsv.gz"
    canonical_tsne_path = PROJECT_ROOT / "metadata/driver/figure2b_sox4/figure2b_sox4_baseline_grid_tsne.tsv.gz"
    canonical_cells = read_table(canonical_cells_path)
    canonical_umap = read_table(canonical_umap_path)
    canonical_tsne = read_table(canonical_tsne_path)
    baselines = {
        "Figure 2 HNF4A": (
            PROJECT_ROOT / "metadata/driver/figure2b_hnf4a/figure2b_hnf4a_plot_cells.tsv.gz",
            PROJECT_ROOT / "metadata/driver/figure2b_hnf4a/figure2b_hnf4a_baseline_grid_umap.tsv.gz",
            PROJECT_ROOT / "metadata/driver/figure2b_hnf4a/figure2b_hnf4a_baseline_grid_tsne.tsv.gz",
        ),
        "Figure 3 EGR1": (
            PROJECT_ROOT / "metadata/driver/figure3b_egr1/figure3b_egr1_plot_cells.tsv.gz",
            PROJECT_ROOT / "metadata/driver/figure3b_egr1/figure3b_egr1_baseline_grid_umap.tsv.gz",
            PROJECT_ROOT / "metadata/driver/figure3b_egr1/figure3b_egr1_baseline_grid_tsne.tsv.gz",
        ),
        "Figure 4 SOX4": (canonical_cells_path, canonical_umap_path, canonical_tsne_path),
    }
    for axis, (cells_path, umap_path, tsne_path) in baselines.items():
        cells = read_table(cells_path)
        umap = read_table(umap_path)
        tsne = read_table(tsne_path)
        add(
            rows,
            "baseline_cell_id_order",
            axis,
            cells["cell_id"].astype(str).tolist() == canonical_cells["cell_id"].astype(str).tolist(),
            int(len(cells)),
            "exact canonical 5,000-cell order",
            "Order-sensitive identity comparison.",
            cells_path,
        )
        for check, columns in [
            ("baseline_umap_points", ["umap_1", "umap_2"]),
            ("baseline_tsne_points", ["tsne_1", "tsne_2"]),
            ("baseline_pseudotime", ["pseudotime"]),
            ("baseline_state_annotation", ["celloracle_state"]),
        ]:
            equal = cells[columns].equals(canonical_cells[columns])
            add(rows, check, axis, equal, equal, True, ";".join(columns), cells_path)
        for space, frame, canonical, path in [
            ("umap", umap, canonical_umap, umap_path),
            ("tsne", tsne, canonical_tsne, tsne_path),
        ]:
            for check, columns in [
                ("grid_coordinates", ["grid_x", "grid_y"]),
                ("baseline_arrow_direction", ["unit_x", "unit_y", "arrow_xend", "arrow_yend"]),
                ("density_mask", ["keep"]),
            ]:
                equal = frame[columns].equals(canonical[columns])
                add(rows, f"baseline_{space}_{check}", axis, equal, equal, True, ";".join(columns), path)

    displacement_specs = {
        "Figure 2 HNF4A": (
            "HNF4A",
            PROJECT_ROOT / "metadata/driver/figure2c_hnf4a/figure2c_hnf4a_matched_cells.tsv.gz",
            PROJECT_ROOT / "metadata/driver/figure2c_hnf4a/figure2c_hnf4a_data_report.json",
        ),
        "Figure 3 EGR1": (
            "EGR1",
            PROJECT_ROOT / "metadata/driver/figure3c_egr1/figure3c_egr1_matched_cells.tsv.gz",
            PROJECT_ROOT / "metadata/driver/figure3c_egr1/figure3c_egr1_data_report.json",
        ),
        "Figure 4 SOX4": (
            "SOX4",
            PROJECT_ROOT / "metadata/driver/figure2c_sox4/figure2c_sox4_matched_cells.tsv.gz",
            PROJECT_ROOT / "metadata/driver/figure2c_sox4/figure2c_sox4_data_report.json",
        ),
    }
    vectors = {}
    for axis, (tf, table_path, report_path) in displacement_specs.items():
        table = read_table(table_path)
        report = json.loads(report_path.read_text(encoding="utf-8-sig"))
        vector = table[["delta_embedding_1", "delta_embedding_2"]].to_numpy(dtype=float)
        vectors[tf] = vector
        target_values = [
            str(report.get(key, ""))
            for key in ["target_tf", "tf", "target"]
            if report.get(key) is not None
        ]
        if "EGR1" == tf and "tf" in table:
            target_values.extend(table["tf"].astype(str).unique().tolist())
        passed = any(tf.lower() in value.lower() for value in target_values)
        add(
            rows,
            "target_tf_contract",
            axis,
            passed,
            target_values,
            tf,
            "Report/table target identity.",
            report_path,
        )
        add(
            rows,
            "finite_unique_displacement",
            axis,
            np.isfinite(vector).all() and np.any(np.linalg.norm(vector, axis=1) > 0),
            {"finite": bool(np.isfinite(vector).all()), "nonzero_cells": int(np.sum(np.linalg.norm(vector, axis=1) > 0))},
            "finite and nonzero",
            "Target-specific displacement table.",
            table_path,
        )
    for left, right in [("EGR1", "SOX4"), ("EGR1", "HNF4A"), ("HNF4A", "SOX4")]:
        difference = float(np.max(np.abs(vectors[left] - vectors[right])))
        add(
            rows,
            "displacement_not_reused",
            f"{left} vs {right}",
            difference > 0,
            difference,
            "> 0",
            "Pointwise maximum absolute difference between target-specific displacement matrices.",
            displacement_specs[f"Figure 3 EGR1" if left == "EGR1" else f"Figure 2 HNF4A"][1],
        )

    grid_paths = {
        "HNF4A_umap": PROJECT_ROOT / "metadata/driver/figure2c_hnf4a/figure2c_hnf4a_inner_product_grid_umap.tsv.gz",
        "HNF4A_tsne": PROJECT_ROOT / "metadata/driver/figure2c_hnf4a/figure2c_hnf4a_inner_product_grid_tsne.tsv.gz",
        "EGR1_umap": PROJECT_ROOT / "metadata/driver/figure3c_egr1/figure3c_egr1_inner_product_grid_umap.tsv.gz",
        "EGR1_tsne": PROJECT_ROOT / "metadata/driver/figure3c_egr1/figure3c_egr1_inner_product_grid_tsne.tsv.gz",
        "SOX4_umap": PROJECT_ROOT / "metadata/driver/figure2c_sox4/figure2c_sox4_inner_product_grid_umap.tsv.gz",
        "SOX4_tsne": PROJECT_ROOT / "metadata/driver/figure2c_sox4/figure2c_sox4_inner_product_grid_tsne.tsv.gz",
    }
    limits = {name: float(np.max(np.abs(valid_grid_scores(path)))) for name, path in grid_paths.items()}
    shared_limit = float(max(limits.values()))
    limit_table = pd.DataFrame(
        [{"source": name, "panel_specific_symmetric_limit": value, "three_axis_shared_symmetric_limit": shared_limit}
         for name, value in limits.items()]
    )
    limit_path = out_dir / "figure2_figure3_figure4_shared_colour_limits.tsv"
    limit_table.to_csv(limit_path, sep="\t", index=False)
    egr1_report_path = PROJECT_ROOT / "metadata/driver/figure3c_egr1/figure3c_egr1_inner_product_report.json"
    egr1_report = json.loads(egr1_report_path.read_text(encoding="utf-8-sig"))
    reported_shared = float(egr1_report["colour_scale"]["three_axis_shared_symmetric_limit"])
    add(
        rows,
        "shared_colour_limit",
        "three-axis",
        np.isclose(reported_shared, shared_limit, rtol=0, atol=1e-12),
        reported_shared,
        shared_limit,
        "Recomputed maximum absolute valid PS across HNF4A, EGR1, and SOX4 UMAP/t-SNE grids.",
        egr1_report_path,
    )

    protected_tokens = {
        "Figure 2 HNF4A": ("hnf4a", "egr1", "sox4"),
        "Figure 3 EGR1": ("egr1", "hnf4a", "sox4"),
        "Figure 4 SOX4": ("sox4", "egr1", "hnf4a"),
    }
    for axis, (own, wrong1, wrong2) in protected_tokens.items():
        paths = [str(path).lower().replace("\\", "/") for path in displacement_specs[axis][1:]]
        passed = all(own in path and wrong1 not in path and wrong2 not in path for path in paths)
        add(
            rows,
            "target_path_isolation",
            axis,
            passed,
            paths,
            f"contains {own}; excludes {wrong1}/{wrong2}",
            "Target-specific data/report paths do not cross namespaces.",
            displacement_specs[axis][1],
        )

    audit = pd.DataFrame(rows)
    audit_path = out_dir / "figure2_figure3_figure4_consistency_audit.tsv"
    audit.to_csv(audit_path, sep="\t", index=False)
    report = {
        "module": "Three-axis figure consistency audit",
        "axes": ["Figure 2 HNF4A", "Figure 3 EGR1", "Figure 4 SOX4"],
        "n_checks": int(len(audit)),
        "n_failed": int(audit["status"].eq("fail").sum()),
        "all_checks_pass": bool(audit["status"].eq("pass").all()),
        "three_axis_shared_symmetric_limit": shared_limit,
        "panel_specific_limits": limits,
        "outputs": {
            "audit": str(audit_path.resolve()),
            "colour_limits": str(limit_path.resolve()),
        },
        "caveat": "The SOX4 implementation is stored under legacy Figure 2 filenames but is audited here as the Figure 4 SOX4 axis requested for manuscript comparison.",
    }
    report_path = out_dir / "figure2_figure3_figure4_consistency_report.json"
    write_json(json_safe(report), report_path)
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = run(args.out_dir)
    print(json.dumps(report, ensure_ascii=False))
    return 0 if report["all_checks_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())

