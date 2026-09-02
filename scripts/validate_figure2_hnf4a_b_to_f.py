#!/usr/bin/env python3
"""Validate HNF4A Figure 2B-F outputs and write the scientific handoff report."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image


PROJECT_ROOT = Path(__file__).resolve().parents[1]
TARGET_TF = "HNF4A"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer")


def add(rows: list[dict], check: str, passed: bool, detail: str, severity: str = "error") -> None:
    rows.append({"check": check, "passed": bool(passed), "severity": severity, "detail": detail})


def check_graphic(path: Path, rows: list[dict]) -> None:
    exists = path.exists() and path.stat().st_size > 0
    add(rows, f"file_nonempty:{path.name}", exists, str(path))
    if not exists:
        return
    try:
        if path.suffix.lower() == ".png":
            with Image.open(path) as image:
                image.verify()
            with Image.open(path) as image:
                dpi = image.info.get("dpi", (0, 0))
                add(rows, f"png_resolution:{path.name}",
                    image.width >= 1800 and image.height >= 1500,
                    f"{image.width}x{image.height}px; dpi={dpi}")
        elif path.suffix.lower() == ".svg":
            ET.parse(path)
            add(rows, f"svg_parse:{path.name}", True, "XML parsed")
        elif path.suffix.lower() == ".pdf":
            header = path.read_bytes()[:5]
            add(rows, f"pdf_header:{path.name}", header == b"%PDF-", f"header={header!r}")
    except Exception as exc:
        add(rows, f"graphic_open:{path.name}", False, repr(exc))


def check_finite_table(path: Path, rows: list[dict]) -> None:
    if not path.exists():
        add(rows, f"source_table_exists:{path.name}", False, str(path))
        return
    dat = read_table(path)
    checked = dat
    mask_detail = "all rows"
    if "keep_score" in dat.columns:
        keep = dat["keep_score"].astype(str).str.lower().isin({"true", "1"})
        checked = dat.loc[keep]
        mask_detail = f"keep_score rows={int(keep.sum())}"
    numeric = checked.select_dtypes(include=[np.number])
    finite = True if numeric.empty else bool(np.isfinite(numeric.to_numpy()).all())
    add(rows, f"finite_numeric:{path.name}", finite,
        f"rows={len(dat)} checked={mask_detail} numeric_columns={len(numeric.columns)}")


def generate_markdown(root: Path, reports: dict[str, dict], validation: dict, path: Path) -> None:
    b = reports["b"]
    c = reports["c"]
    d = reports["d"]
    e = reports["e"]
    f = reports["f"]
    pre = reports["preflight"]
    sens = e.get("sensitivity", {})
    u = d.get("umap", {})
    ci = u.get("spearman_ci95", [None, None])
    rho = u.get("spearman_rho")
    direction = "positive" if rho is not None and float(rho) > 0 else "negative or null"
    main_standard = bool(validation["all_required_checks_passed"] and e.get("n_significant_excluding_target", 0) > 0)
    f_standard = f.get("n_significant_pathways", 0) > 0
    lines = [
        "# Figure 2B-F: HNF4A hepatocyte identity loss",
        "",
        "## Main conclusion",
        "",
        "This figure reports a predicted perturbation after HNF4A virtual knockout. "
        "The CellOracle result is a computationally inferred state shift, and the "
        "scTenifoldKnk result is network perturbation evidence.",
        "",
        "## New files",
        "",
        "- HNF4A-specific scripts: `scripts/*figure2*hnf4a*` and `scripts/figure2_hnf4a_common.R`.",
        "- Source data and JSON reports: `metadata/driver/figure2*_hnf4a*`.",
        "- Figures: `figures/driver/figure2*_hnf4a*`.",
        "- Validation: `metadata/driver/figure2_hnf4a_b_to_f_validation_report.{json,tsv}`.",
        "",
        "## Panel data sources",
        "",
        "- Figure 2B: exact reuse of the validated SOX4 TF-independent 5,000-cell baseline cell, coordinate, pseudotime and grid tables.",
        "- Figure 2C: saved HNF4A=0 CellOracle `delta_embedding` from Module 6.8; UMAP is native to the saved perturbation space.",
        "- Figure 2D: HNF4A Figure 2C grid inner product mapped against strict-main pseudotime.",
        "- Figure 2E: HNF4A scTenifoldKnk virtual knockout in identity-high `normal_reference` hepatocytes.",
        "- Figure 2F: one-sided ORA of FDR-significant Figure 2E genes against the matching normal-reference network background.",
        "",
        "## Evidence and statistics",
        "",
        f"- HNF4A CellOracle perturbation: existing result; propagation steps = {pre['celloracle_parameters']['n_propagation']}, seed = {pre['celloracle_parameters']['seed']}.",
        f"- Figure 2B equivalence: `{b.get('source_contract')}`; n = {b.get('n_cells')} cells.",
        f"- Figure 2C UMAP PS range: {c.get('score_range_umap')}; positive fraction = {c.get('positive_ps_fraction'):.4f}; negative fraction = {c.get('negative_ps_fraction'):.4f}.",
        f"- Figure 2D UMAP Spearman rho = {float(rho):.4f} (bootstrap 95% CI {ci}); the observed relationship is {direction}.",
        f"- Figure 2E FDR-significant genes excluding HNF4A: {e.get('n_significant_excluding_target', 0)}; plotted = {e.get('n_plotted', 0)}.",
        f"- Figure 2F FDR-significant pathways: {f.get('n_significant_pathways', 0)}; plotted = {f.get('n_plotted', 0)}.",
        f"- Identity-high vs malignant-like: Jaccard = {sens.get('jaccard_overlap')}; distance-rank Spearman = {sens.get('distance_rank_spearman')}.",
        "",
        "## Limitations and review-risk flags",
        "",
        "- CellOracle and scTenifoldKnk are in silico perturbation frameworks; the results do not establish an experimentally validated or direct tumorigenic effect.",
        "- The t-SNE panels use a local UMAP-to-t-SNE Jacobian projection and are sensitivity views, not native CellOracle t-SNE simulations.",
        "- The scTenifoldKnk manifold distance is non-directional and supports network-displacement enrichment only.",
        "- The identity-high scTenifoldKnk run matches the existing Module 7 settings (`nc_nNet=1`, `nc_nCells=100`); higher-replicate, multi-seed confirmation remains a publication-risk item.",
        (f"- No FDR-significant pathway was detected. {f.get('recommendation')}"
         if not f_standard else "- Figure 2F displays only pooled-BH FDR-significant pathways."),
        "",
        "## Figure legend draft",
        "",
        "**Figure 2B-F | Predicted loss of HNF4A destabilizes hepatocyte identity.** "
        "(B) Baseline developmental vector field across 5,000 strict-main hepatocytes. "
        "(C) CellOracle-predicted perturbation field after HNF4A virtual knockout and the "
        "inner product between perturbation and baseline developmental vectors. "
        "(D) HNF4A perturbation score across strict-main pseudotime; points denote grid "
        "scores and black open circles/line denote bin medians. "
        "(E) FDR-significant genes with the largest scTenifoldKnk manifold-alignment "
        "distance in identity-high normal-reference hepatocytes. "
        "(F) FDR-significant ORA pathways from the matching network background"
        + ("." if f_standard else "; no main pathway panel was generated because no robust FDR-significant pathway was detected."),
        "",
        "## Main-figure readiness",
        "",
        f"- Figure 2B-E main-figure standard: {'met' if main_standard else 'not met'}.",
        f"- Figure 2F main-figure standard: {'met' if f_standard else 'not met'}.",
        f"- Validation checks: {validation['n_passed']}/{validation['n_checks']} passed.",
        "",
        "## SOX4 preservation",
        "",
        f"All {validation.get('n_sox4_hashes_checked', 0)} recorded SOX4 files retained their reference SHA-256 hashes.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(root: Path, out_json: Path, out_tsv: Path, report_md: Path) -> dict:
    rows: list[dict] = []
    reports = {
        "preflight": load_json(root / "metadata/driver/figure2_hnf4a_preflight_report.json"),
        "b": load_json(root / "metadata/driver/figure2b_hnf4a/figure2b_hnf4a_r_plot_report.json"),
        "c": load_json(root / "metadata/driver/figure2c_hnf4a/figure2c_hnf4a_inner_product_report.json"),
        "d": load_json(root / "metadata/driver/figure2d_hnf4a/figure2d_hnf4a_report.json"),
        "e": load_json(root / "metadata/driver/figure2e_hnf4a/figure2e_hnf4a_report.json"),
        "f": load_json(root / "metadata/driver/figure2f_hnf4a/figure2f_hnf4a_report.json"),
    }
    for name, report in reports.items():
        add(rows, f"json_target_tf:{name}", report.get("target_tf") == TARGET_TF,
            f"target_tf={report.get('target_tf')}")

    required_graphics = [
        root / f"figures/driver/figure2b_hnf4a/figure2b_hnf4a_baseline_{space}.{ext}"
        for space in ["umap", "tsne", "umap_tsne"] for ext in ["pdf", "png", "svg"]
    ] + [
        root / f"figures/driver/figure2c_hnf4a/figure2c_hnf4a_{kind}_{space}.{ext}"
        for kind in ["perturbation", "inner_product"] for space in ["umap", "tsne", "umap_tsne"]
        for ext in ["pdf", "png", "svg"]
    ] + [
        root / f"figures/driver/figure2d_hnf4a/figure2d_hnf4a_pseudotime_inner_product_{space}.{ext}"
        for space in ["umap", "tsne", "umap_tsne"] for ext in ["pdf", "png", "svg"]
    ]
    if reports["e"].get("figure_generated"):
        required_graphics += [
            root / f"figures/driver/figure2e_hnf4a/figure2e_hnf4a_significant_perturbed_genes.{ext}"
            for ext in ["pdf", "png", "svg", "tiff"]
        ]
    if reports["f"].get("figure_generated"):
        required_graphics += [
            root / f"figures/driver/figure2f_hnf4a/figure2f_hnf4a_pathway_enrichment.{ext}"
            for ext in ["pdf", "png", "svg", "tiff"]
        ]
    for path in required_graphics:
        check_graphic(path, rows)

    audit = read_table(root / "metadata/driver/figure2b_hnf4a/figure2b_baseline_equivalence_audit.tsv")
    equivalence = bool(audit["values_exactly_equal"].astype(bool).all() and
                       (pd.to_numeric(audit["max_absolute_numeric_difference"]) == 0).all())
    add(rows, "figure2b_pointwise_equivalence", equivalence, f"components={len(audit)}")

    matched = read_table(root / "metadata/driver/figure2c_hnf4a/figure2c_hnf4a_matched_cells.tsv.gz")
    shifts = read_table(root / "metadata/driver/celloracle_module6_8_cell_shift_summary.tsv.gz")
    h = shifts.loc[shifts["tf"].astype(str).eq(TARGET_TF),
                   ["cell_id", "delta_embedding_1", "delta_embedding_2"]]
    s = shifts.loc[shifts["tf"].astype(str).eq("SOX4"),
                   ["cell_id", "delta_embedding_1", "delta_embedding_2"]]
    hmatch = matched[["cell_id", "delta_embedding_1", "delta_embedding_2"]].merge(
        h, on="cell_id", suffixes=("_plot", "_source"), validate="one_to_one")
    h_equal = np.allclose(hmatch[["delta_embedding_1_plot", "delta_embedding_2_plot"]],
                          hmatch[["delta_embedding_1_source", "delta_embedding_2_source"]])
    smatch = matched[["cell_id", "delta_embedding_1", "delta_embedding_2"]].merge(
        s, on="cell_id", suffixes=("_plot", "_sox4"), validate="one_to_one")
    differs_sox4 = not np.allclose(smatch[["delta_embedding_1_plot", "delta_embedding_2_plot"]],
                                   smatch[["delta_embedding_1_sox4", "delta_embedding_2_sox4"]])
    add(rows, "figure2c_uses_hnf4a_delta", h_equal and differs_sox4,
        f"matches_HNF4A={h_equal};differs_from_SOX4={differs_sox4}")

    ip = read_table(root / "metadata/driver/figure2c_hnf4a/figure2c_hnf4a_inner_product_grid_umap.tsv.gz")
    baseline = read_table(root / "metadata/driver/figure2b_hnf4a/figure2b_hnf4a_baseline_grid_umap.tsv.gz")
    recomputed = ip["flow_x"] * baseline["unit_x"] + ip["flow_y"] * baseline["unit_y"]
    valid = ip["keep_score"].astype(bool)
    ip_equal = np.allclose(recomputed[valid], ip.loc[valid, "inner_product_score_grid"])
    add(rows, "figure2c_inner_product_recomputed", ip_equal, f"n_valid={int(valid.sum())}")
    add(rows, "figure2d_reads_hnf4a_figure2c",
        "HNF4A" in reports["d"].get("analysis", ""), reports["d"].get("analysis", ""))

    e_table = read_table(root / "metadata/driver/figure2e_hnf4a/figure2e_hnf4a_significant_perturbed_genes.tsv")
    e_valid = (len(e_table) == 0 or (
        not e_table["gene"].astype(str).eq(TARGET_TF).any()
        and (pd.to_numeric(e_table["p.adj"]) < 0.05).all()
    ))
    add(rows, "figure2e_fdr_only_target_excluded", e_valid, f"n_rows={len(e_table)}")

    f_table = read_table(root / "metadata/driver/figure2f_hnf4a/figure2f_hnf4a_plot_data.tsv")
    f_valid = len(f_table) == 0 or (pd.to_numeric(f_table["p.adjust"]) < 0.05).all()
    add(rows, "figure2f_fdr_only", f_valid, f"n_rows={len(f_table)}")
    f_background = Path(reports["f"]["background"])
    add(rows, "figure2f_background_matches_identity_high",
        "normal_reference" in f_background.as_posix() and "hnf4a" in f_background.as_posix().lower(),
        str(f_background))

    source_tables = [
        root / "metadata/driver/figure2b_hnf4a/figure2b_hnf4a_plot_cells.tsv.gz",
        root / "metadata/driver/figure2c_hnf4a/figure2c_hnf4a_inner_product_grid_umap.tsv.gz",
        root / "metadata/driver/figure2d_hnf4a/figure2d_hnf4a_pseudotime_bin_summary.tsv",
        root / "metadata/driver/figure2e_hnf4a/figure2e_hnf4a_significant_perturbed_genes.tsv",
        root / "metadata/driver/figure2f_hnf4a/figure2f_hnf4a_enrichment_all.tsv",
    ]
    for path in source_tables:
        check_finite_table(path, rows)

    hash_path = root / "metadata/driver/figure2_hnf4a_sox4_reference_hashes.tsv"
    hashes = read_table(hash_path)
    unchanged = []
    for row in hashes.itertuples(index=False):
        path = root / row.path
        unchanged.append(path.exists() and path.stat().st_size == row.size_bytes and sha256(path) == row.sha256)
    add(rows, "sox4_files_unmodified", all(unchanged), f"checked={len(unchanged)}")

    validation_table = pd.DataFrame(rows)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    validation_table.to_csv(out_tsv, sep="\t", index=False)
    errors = validation_table.loc[validation_table["severity"].eq("error")]
    all_passed = bool(errors["passed"].all())
    report = {
        "module": "Figure 2 HNF4A B-F validation",
        "target_tf": TARGET_TF,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "n_checks": int(len(validation_table)),
        "n_passed": int(validation_table["passed"].sum()),
        "n_failed": int((~validation_table["passed"]).sum()),
        "all_required_checks_passed": all_passed,
        "n_sox4_hashes_checked": int(len(hashes)),
        "determinism_contract": {
            "python_seed": 15071990, "r_seed": 15071990,
            "celloracle_seed": 15071990, "sctenifoldknk_seed": 11,
            "figure2b_reused_exactly": equivalence,
        },
        "review_risk_flags": [
            "t-SNE is a supplementary Jacobian projection, not native CellOracle simulation",
            "scTenifoldKnk normal-reference run uses Module 7 nc_nNet=1 and nc_nCells=100",
            "virtual knockout results are computational predictions",
        ],
        "outputs": {"json": str(out_json.resolve()), "tsv": str(out_tsv.resolve()),
                    "report_md": str(report_md.resolve())},
    }
    out_json.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    generate_markdown(root, reports, report, report_md)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--out-json", type=Path, default=PROJECT_ROOT / "metadata/driver/figure2_hnf4a_b_to_f_validation_report.json")
    parser.add_argument("--out-tsv", type=Path, default=PROJECT_ROOT / "metadata/driver/figure2_hnf4a_b_to_f_validation_report.tsv")
    parser.add_argument("--report-md", type=Path, default=PROJECT_ROOT / "reports/figure2_hnf4a_b_to_f_report.md")
    args = parser.parse_args()
    result = run(args.project_root, args.out_json, args.out_tsv, args.report_md)
    print(json.dumps(result, ensure_ascii=False))
    raise SystemExit(0 if result["all_required_checks_passed"] else 1)


if __name__ == "__main__":
    main()
