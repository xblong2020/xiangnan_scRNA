#!/usr/bin/env python3
"""Validate Figure 3 EGR1 A-F outputs, provenance, images, and isolation."""

from __future__ import annotations

import argparse
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path

import fitz
import numpy as np
import pandas as pd
from PIL import Image

try:
    from figure3_egr1_common import PROJECT_ROOT, TARGET_TF, json_safe, sha256_file, write_json
except ModuleNotFoundError:
    from scripts.figure3_egr1_common import PROJECT_ROOT, TARGET_TF, json_safe, sha256_file, write_json


DEFAULT_OUT_DIR = PROJECT_ROOT / "metadata/driver/figure3_egr1_validation"


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def add(
    rows: list[dict],
    check_id: str,
    requirement: str,
    status: str,
    observed,
    expected: str,
    evidence: Path | str,
    details: str,
) -> None:
    if status not in {"pass", "warning", "fail"}:
        raise ValueError(status)
    rows.append(
        {
            "check_id": check_id,
            "requirement": requirement,
            "status": status,
            "observed": json.dumps(json_safe(observed), ensure_ascii=False)
            if isinstance(observed, (dict, list))
            else observed,
            "expected": expected,
            "evidence": str(Path(evidence).resolve()) if evidence else "",
            "details": details,
        }
    )


def inspect_image(path: Path) -> dict:
    suffix = path.suffix.lower()
    result = {"path": str(path.resolve()), "exists": path.exists(), "size_bytes": path.stat().st_size if path.exists() else 0}
    if not path.exists() or path.stat().st_size == 0:
        result["readable"] = False
        return result
    try:
        if suffix in {".png", ".tif", ".tiff"}:
            with Image.open(path) as image:
                image.verify()
            with Image.open(path) as image:
                result.update(
                    {
                        "readable": True,
                        "width_px": int(image.width),
                        "height_px": int(image.height),
                        "dpi": [float(value) for value in image.info.get("dpi", ())],
                        "format": image.format,
                    }
                )
        elif suffix == ".pdf":
            document = fitz.open(path)
            result.update({"readable": len(document) > 0, "pages": len(document)})
            document.close()
        elif suffix == ".svg":
            ET.parse(path)
            result.update({"readable": True, "format": "SVG"})
        else:
            result["readable"] = True
    except Exception as error:
        result.update({"readable": False, "error": repr(error)})
    return result


def finite_plot_rows(path: Path, filter_column: str | None = None) -> tuple[bool, dict]:
    frame = pd.read_csv(path, sep="\t", compression="infer")
    if filter_column and filter_column in frame:
        frame = frame.loc[frame[filter_column].astype(str).str.lower().isin({"true", "1"})]
    numeric = frame.select_dtypes(include=[np.number])
    finite = bool(np.isfinite(numeric.to_numpy(dtype=float)).all()) if len(numeric.columns) else True
    return finite, {"n_rows": int(len(frame)), "numeric_columns": numeric.columns.tolist()}


def expected_images(e_report: dict, f_report: dict) -> list[Path]:
    images = [
        PROJECT_ROOT / "figures/driver/figure3a_stress_transition/figure3a_stress_transition_selection.pdf",
        PROJECT_ROOT / "figures/driver/figure3a_stress_transition/figure3a_stress_transition_selection.png",
        PROJECT_ROOT / "figures/driver/figure3a_stress_transition/figure3a_stress_transition_selection.svg",
        PROJECT_ROOT / "figures/driver/figure3a_stress_transition/figure3a_stress_transition_selection.tiff",
    ]
    for panel, stems in {
        "figure3b_egr1": [
            "figure3b_egr1_baseline_umap",
            "figure3b_egr1_baseline_tsne",
            "figure3b_egr1_baseline_umap_tsne",
        ],
        "figure3c_egr1": [
            "figure3c_egr1_perturbation_umap",
            "figure3c_egr1_perturbation_tsne",
            "figure3c_egr1_perturbation_umap_tsne",
            "figure3c_egr1_inner_product_umap",
            "figure3c_egr1_inner_product_tsne",
            "figure3c_egr1_inner_product_umap_tsne",
        ],
        "figure3d_egr1": [
            "figure3d_egr1_pseudotime_inner_product_umap",
            "figure3d_egr1_pseudotime_inner_product_tsne",
            "figure3d_egr1_pseudotime_inner_product_umap_tsne",
        ],
    }.items():
        for stem in stems:
            for suffix in [".pdf", ".png", ".svg"]:
                images.append(PROJECT_ROOT / "figures/driver" / panel / f"{stem}{suffix}")
    for panel, stems in {
        "figure3c_egr1": [
            "figure3c_egr1_inner_product_umap_panel_specific",
            "figure3c_egr1_inner_product_tsne_panel_specific",
            "figure3c_egr1_inner_product_umap_tsne_panel_specific",
        ],
        "figure3d_egr1": [
            "figure3d_egr1_pseudotime_inner_product_umap_panel_specific",
            "figure3d_egr1_pseudotime_inner_product_tsne_panel_specific",
            "figure3d_egr1_pseudotime_inner_product_umap_tsne_panel_specific",
        ],
    }.items():
        for stem in stems:
            for suffix in [".pdf", ".png", ".svg"]:
                images.append(PROJECT_ROOT / "figures/driver" / panel / f"{stem}{suffix}")
    for suffix in [".pdf", ".png", ".svg", ".tiff"]:
        images.append(
            PROJECT_ROOT
            / "figures/driver/figure3e_egr1_sensitivity"
            / f"figure3e_egr1_state_sensitivity_summary{suffix}"
        )
    for suffix in [".pdf", ".png"]:
        images.append(
            PROJECT_ROOT
            / "figures/driver/figure3_egr1_preview"
            / f"figure3_egr1_a_to_f_preview{suffix}"
        )
    if e_report.get("formal_plot_generated"):
        for suffix in [".pdf", ".png", ".svg", ".tiff"]:
            images.append(
                PROJECT_ROOT
                / "figures/driver/figure3e_egr1"
                / f"figure3e_egr1_significant_perturbed_genes{suffix}"
            )
    if f_report.get("formal_plot_generated"):
        for suffix in [".pdf", ".png", ".svg", ".tiff"]:
            images.append(
                PROJECT_ROOT
                / "figures/driver/figure3f_egr1"
                / f"figure3f_egr1_pathway_enrichment{suffix}"
            )
    return images


def run(out_dir: Path) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    e_report_path = PROJECT_ROOT / "metadata/driver/figure3e_egr1/figure3e_egr1_report.json"
    f_report_path = PROJECT_ROOT / "metadata/driver/figure3f_egr1/figure3f_egr1_report.json"
    e_report = read_json(e_report_path)
    f_report = read_json(f_report_path)

    required_data = [
        PROJECT_ROOT / "metadata/driver/figure3_egr1_preflight/figure3_egr1_preflight_report.json",
        PROJECT_ROOT / "metadata/driver/figure3_egr1_preflight/figure3_egr1_preflight_report.tsv",
        PROJECT_ROOT / "metadata/driver/figure3a_stress_transition/figure3a_candidate_evidence_matrix.tsv",
        PROJECT_ROOT / "metadata/driver/figure3a_stress_transition/figure3a_candidate_selection_report.json",
        PROJECT_ROOT / "metadata/driver/figure3b_egr1/figure3b_baseline_equivalence_audit.tsv",
        PROJECT_ROOT / "metadata/driver/figure3c_egr1/figure3c_egr1_cell_level_scores.tsv.gz",
        PROJECT_ROOT / "metadata/driver/figure3c_egr1/figure3c_egr1_inner_product_grid_umap.tsv.gz",
        PROJECT_ROOT / "metadata/driver/figure3c_egr1/figure3c_egr1_inner_product_grid_tsne.tsv.gz",
        PROJECT_ROOT / "metadata/driver/figure3d_egr1/figure3d_egr1_pseudotime_bin_summary.tsv",
        PROJECT_ROOT / "metadata/driver/figure3d_egr1/figure3d_egr1_stage_comparison.tsv",
        PROJECT_ROOT / "metadata/driver/figure3d_egr1/figure3d_egr1_bootstrap_peak_summary.tsv",
        PROJECT_ROOT / "metadata/driver/figure3e_egr1/figure3e_egr1_significant_perturbed_genes.tsv",
        PROJECT_ROOT / "metadata/driver/figure3f_egr1/figure3f_egr1_enrichment_all.tsv",
        PROJECT_ROOT / "metadata/driver/figure3f_egr1/figure3f_egr1_plot_data.tsv",
        PROJECT_ROOT / "metadata/driver/figure3e_egr1_sensitivity/figure3e_egr1_sensitivity_report.json",
        PROJECT_ROOT / "metadata/driver/figure3e_egr1_sensitivity/figure3e_egr1_sensitivity_summary.tsv",
        PROJECT_ROOT / "metadata/driver/figure3e_egr1_determinism/figure3e_egr1_same_seed_determinism_report.json",
        PROJECT_ROOT / "metadata/driver/figure3e_egr1_determinism/figure3e_egr1_same_seed_determinism_audit.tsv",
        PROJECT_ROOT / "figures/driver/figure3_egr1_preview/figure3_egr1_a_to_f_preview_report.json",
        PROJECT_ROOT / "metadata/driver/three_axis_figure_consistency/figure2_figure3_figure4_consistency_audit.tsv",
    ]
    missing = [str(path) for path in required_data if not path.exists() or path.stat().st_size == 0]
    add(
        rows,
        "01_required_outputs",
        "All required outputs exist and are non-empty",
        "pass" if not missing else "fail",
        {"n_expected": len(required_data), "missing": missing},
        "no missing files",
        out_dir,
        "Formal 3E/3F images are conditional on their FDR-supported report flags.",
    )

    images = expected_images(e_report, f_report)
    image_results = [inspect_image(path) for path in images]
    unreadable = [result for result in image_results if not result.get("readable")]
    add(
        rows,
        "02_image_readability",
        "PDF, PNG, SVG, and conditional TIFF files open successfully",
        "pass" if not unreadable else "fail",
        {"n_images": len(images), "unreadable": unreadable},
        "all readable",
        out_dir,
        "Images were opened with Pillow, PyMuPDF, or XML parser.",
    )
    raster = [result for result in image_results if Path(result["path"]).suffix.lower() in {".png", ".tiff", ".tif"}]
    raster_ok = all(
        result.get("readable")
        and result.get("width_px", 0) >= 1800
        and result.get("height_px", 0) >= 1500
        for result in raster
    )
    add(
        rows,
        "03_raster_resolution",
        "PNG/TIFF exports meet publication raster resolution",
        "pass" if raster_ok else "fail",
        raster,
        "at least 1800 x 1500 pixels at configured 600 dpi",
        out_dir,
        "Pixel dimensions provide a robust cross-device validation when embedded DPI metadata is rounded.",
    )

    finite_specs = [
        (PROJECT_ROOT / "metadata/driver/figure3b_egr1/figure3b_egr1_plot_cells.tsv.gz", None),
        (PROJECT_ROOT / "metadata/driver/figure3c_egr1/figure3c_egr1_inner_product_grid_umap.tsv.gz", "show_score"),
        (PROJECT_ROOT / "metadata/driver/figure3c_egr1/figure3c_egr1_inner_product_grid_tsne.tsv.gz", "show_score"),
        (PROJECT_ROOT / "metadata/driver/figure3d_egr1/figure3d_egr1_pseudotime_inner_product_umap.tsv.gz", "valid"),
        (PROJECT_ROOT / "metadata/driver/figure3d_egr1/figure3d_egr1_pseudotime_inner_product_tsne.tsv.gz", "valid"),
    ]
    finite_results = []
    for path, column in finite_specs:
        passed, observed = finite_plot_rows(path, column)
        finite_results.append({"path": str(path), "passed": passed, **observed})
    add(
        rows,
        "04_finite_plot_data",
        "No NA/Inf enters plotted rows",
        "pass" if all(result["passed"] for result in finite_results) else "fail",
        finite_results,
        "all plotted numeric values finite",
        finite_specs[0][0],
        "Invalid low-density grid rows may remain in source data but are excluded by plot masks.",
    )

    baseline_audit_path = PROJECT_ROOT / "metadata/driver/figure3b_egr1/figure3b_baseline_equivalence_audit.tsv"
    baseline = pd.read_csv(baseline_audit_path, sep="\t")
    add(
        rows,
        "05_baseline_equivalence",
        "Figure 3B exactly matches the common baseline",
        "pass" if baseline["status"].eq("pass").all() else "fail",
        baseline[["component", "status"]].to_dict(orient="records"),
        "all pass",
        baseline_audit_path,
        "Cell IDs, points, pseudotime, grids, arrows, and masks are order-sensitive comparisons.",
    )

    c_data_report_path = PROJECT_ROOT / "metadata/driver/figure3c_egr1/figure3c_egr1_data_report.json"
    c_data_report = read_json(c_data_report_path)
    true_source = (
        c_data_report.get("target_tf") == TARGET_TF
        and c_data_report.get("condition") == {TARGET_TF: 0.0}
        and "celloracle_module6_8_cell_shift_summary" in c_data_report.get("source_perturbation", "")
    )
    add(
        rows,
        "06_true_egr1_delta",
        "Figure 3C uses true saved EGR1 CellOracle delta_embedding",
        "pass" if true_source else "fail",
        {
            "target_tf": c_data_report.get("target_tf"),
            "condition": c_data_report.get("condition"),
            "source": c_data_report.get("source_perturbation"),
        },
        "EGR1=0 Module 6.8 source",
        c_data_report_path,
        "Expression correlation is not accepted as displacement.",
    )
    c_source_lower = c_data_report.get("source_perturbation", "").lower()
    add(
        rows,
        "07_no_protected_delta_reuse",
        "Figure 3C never reads SOX4/HNF4A displacement",
        "pass" if "sox4" not in c_source_lower and "hnf4a" not in c_source_lower else "fail",
        c_source_lower,
        "no sox4 or hnf4a token",
        c_data_report_path,
        "The shared baseline is allowed; target displacement is isolated.",
    )

    d_report_path = PROJECT_ROOT / "metadata/driver/figure3d_egr1/figure3d_egr1_report.json"
    d_report = read_json(d_report_path)
    add(
        rows,
        "08_figure3d_reads_figure3c",
        "Figure 3D is derived from Figure 3C EGR1 grid data",
        "pass" if d_report.get("target_tf") == TARGET_TF and d_report.get("score_definition", "").startswith("EGR1") else "fail",
        {"target_tf": d_report.get("target_tf"), "score_definition": d_report.get("score_definition")},
        "EGR1 Figure 3C score",
        d_report_path,
        "The dedicated statistics script accepts Figure 3C as its only default score source.",
    )

    e_source = pd.read_csv(e_report["source_table"], sep="\t")
    add(
        rows,
        "09_egr1_self_excluded",
        "Figure 3E excludes EGR1 itself",
        "pass" if not e_source["gene"].astype(str).eq(TARGET_TF).any() else "fail",
        e_source["gene"].astype(str).eq(TARGET_TF).sum(),
        0,
        e_report["source_table"],
        "Source table is the exact set offered to the formal plot.",
    )
    e_fdr_ok = len(e_source) == 0 or pd.to_numeric(e_source["p.adj"], errors="coerce").lt(0.05).all()
    add(
        rows,
        "10_figure3e_fdr_only",
        "Figure 3E contains only p.adj < 0.05 genes",
        "pass" if e_fdr_ok else "fail",
        int(len(e_source)),
        "all p.adj < 0.05",
        e_report["source_table"],
        "No nominally significant rows are allowed.",
    )
    no_fill = (
        bool(e_report.get("non_significant_fill_used") is False)
        and int(e_report.get("n_plotted", 0)) == min(20, int(e_report.get("n_significant_excluding_target", 0)))
    )
    add(
        rows,
        "11_no_nonsignificant_fill",
        "Figure 3E never fills Top 20 with non-significant genes",
        "pass" if no_fill else "fail",
        {
            "n_significant": e_report.get("n_significant_excluding_target"),
            "n_plotted": e_report.get("n_plotted"),
            "fill_used": e_report.get("non_significant_fill_used"),
        },
        "n_plotted=min(20,n_significant), fill=false",
        e_report_path,
        "Zero significant genes correctly suppresses the formal plot.",
    )
    add(
        rows,
        "12_stress_transition_main_state",
        "Figure 3E main subset is a stress-transition state",
        "pass" if e_report.get("subset") in {"stressed_injured", "stressed_regenerative", "intermediate_pseudotime"} else "fail",
        e_report.get("subset"),
        "stressed, stressed+regenerative, or intermediate",
        e_report_path,
        "Direct fallback to malignant-like is not accepted without upstream eligibility failure.",
    )
    background_match = (
        f_report.get("subset") == e_report.get("subset")
        and bool(f_report.get("background_matches_figure3e_subset"))
        and e_report.get("subset", "") in f_report.get("background", "")
    )
    add(
        rows,
        "13_matching_enrichment_background",
        "Figure 3F background matches Figure 3E state network",
        "pass" if background_match else "fail",
        {"figure3e_subset": e_report.get("subset"), "figure3f_subset": f_report.get("subset"), "background": f_report.get("background")},
        "same subset and network genes",
        f_report_path,
        "Background is not borrowed from malignant-like or another target.",
    )
    add(
        rows,
        "14_no_nominal_fdr_claim",
        "Figure 3F never presents nominal results as FDR significant",
        "pass" if f_report.get("nominal_results_used_as_formal") is False else "fail",
        f_report.get("nominal_results_used_as_formal"),
        False,
        f_report_path,
        "Formal plot generation is conditioned on global BH p.adjust < 0.05.",
    )

    json_paths = [
        path
        for parent in [
            PROJECT_ROOT / "metadata/driver/figure3_egr1_preflight",
            PROJECT_ROOT / "metadata/driver/figure3a_stress_transition",
            PROJECT_ROOT / "metadata/driver/figure3b_egr1",
            PROJECT_ROOT / "metadata/driver/figure3c_egr1",
            PROJECT_ROOT / "metadata/driver/figure3d_egr1",
            PROJECT_ROOT / "metadata/driver/figure3e_egr1",
            PROJECT_ROOT / "metadata/driver/figure3e_egr1_sensitivity",
            PROJECT_ROOT / "metadata/driver/figure3f_egr1",
        ]
        if parent.exists()
        for path in parent.glob("*.json")
    ]
    wrong_targets = []
    for path in json_paths:
        report = read_json(path)
        target = report.get("target_tf")
        if target is not None and str(target).upper() != TARGET_TF:
            wrong_targets.append({"path": str(path), "target_tf": target})
    add(
        rows,
        "15_json_target_tf",
        "All Figure 3 JSON target_tf values are EGR1",
        "pass" if not wrong_targets else "fail",
        {"n_json": len(json_paths), "wrong": wrong_targets},
        "target_tf=EGR1",
        json_paths[0],
        "Reports without a target_tf key are general audit reports and are not treated as target mismatches.",
    )

    fingerprint_path = PROJECT_ROOT / "metadata/driver/figure3_egr1_preflight/protected_sox4_hnf4a_fingerprints_before.tsv"
    fingerprints = pd.read_csv(fingerprint_path, sep="\t")
    changed = []
    for row in fingerprints.itertuples(index=False):
        path = Path(row.path)
        current = sha256_file(path) if path.exists() else None
        if current != row.sha256:
            stat = path.stat() if path.exists() else None
            changed.append(
                {
                    "path": str(path),
                    "before_sha256": row.sha256,
                    "after_sha256": current,
                    "before_size_bytes": int(row.size_bytes),
                    "after_size_bytes": int(stat.st_size) if stat else None,
                    "before_mtime_ns": int(row.mtime_ns),
                    "after_mtime_ns": int(stat.st_mtime_ns) if stat else None,
                }
            )
    incident_path = out_dir / "figure3_egr1_protected_asset_change_incident.tsv"
    if changed:
        pd.DataFrame(changed).to_csv(incident_path, sep="\t", index=False)
    add(
        rows,
        "16_protected_assets_unchanged",
        "No existing SOX4/HNF4A file was modified",
        "pass" if not changed else "warning",
        {"n_protected": int(len(fingerprints)), "changed": changed},
        "zero changed hashes",
        incident_path if changed else fingerprint_path,
        (
            "SHA-256 comparison against the preflight manifest. Changed protected files are "
            "reported as an unresolved concurrent/synchronization incident and are never "
            "silently re-baselined or reverted by the Figure 3 workflow."
        ),
    )

    run_report_path = PROJECT_ROOT / "metadata/driver/figure3e_egr1/figure3e_egr1_stressed_regenerative_run_report.json"
    run_report = read_json(run_report_path)
    determinism_path = (
        PROJECT_ROOT
        / "metadata/driver/figure3e_egr1_determinism/figure3e_egr1_same_seed_determinism_report.json"
    )
    determinism = read_json(determinism_path)
    deterministic_contract = (
        run_report.get("parameters", {}).get("nc_nNet") == 10
        and run_report.get("parameters", {}).get("nc_nCells_used") == 500
        and len(run_report.get("parameters", {}).get("seeds", [])) >= 3
        and run_report.get("n_successful_seeds", 0) >= 3
        and bool(determinism.get("independent_repeat"))
        and bool(determinism.get("numeric_values_reproducible"))
    )
    add(
        rows,
        "17_fixed_seed_replication",
        "Formal EGR1 network uses fixed multiple seeds and publication parameters",
        "pass" if deterministic_contract else "warning",
        {
            "parameters": run_report.get("parameters"),
            "n_successful_seeds": run_report.get("n_successful_seeds"),
            "same_seed_repeat": determinism,
        },
        "nc_nNet=10, nc_nCells=500, at least three successful fixed seeds, and an independently identical same-seed repeat",
        determinism_path,
        "Multiple fixed seeds test stability; the independent seed-15071990 repeat directly tests determinism.",
    )

    consistency_path = PROJECT_ROOT / "metadata/driver/three_axis_figure_consistency/figure2_figure3_figure4_consistency_audit.tsv"
    consistency = pd.read_csv(consistency_path, sep="\t")
    shared_ok = consistency.loc[consistency["check"].eq("shared_colour_limit"), "status"].eq("pass").all()
    add(
        rows,
        "18_shared_colour_limits",
        "Three-axis shared symmetric limits are computed correctly",
        "pass" if shared_ok else "fail",
        consistency.loc[consistency["check"].eq("shared_colour_limit")].to_dict(orient="records"),
        "pass",
        consistency_path,
        "Main Figure 3C/D outputs use the recomputed HNF4A/EGR1/SOX4 shared limit.",
    )

    a_report_path = PROJECT_ROOT / "metadata/driver/figure3a_stress_transition/figure3a_candidate_selection_report.json"
    a_report = read_json(a_report_path)
    a_matrix_path = PROJECT_ROOT / "metadata/driver/figure3a_stress_transition/figure3a_candidate_evidence_matrix.tsv"
    a_matrix = pd.read_csv(a_matrix_path, sep="\t")
    required_metrics = [
        "celloracle_evidence",
        "sctenifoldknk_evidence",
        "cross_method_concordance",
        "transition_state_specificity",
        "temporal_positioning",
        "leave_one_dataset_out_stability",
        "proliferation_dependency",
        "generic_stress_risk",
        "literature_overlap",
    ]
    a_real = all(column in a_matrix and np.isfinite(pd.to_numeric(a_matrix[column], errors="coerce")).all() for column in required_metrics)
    add(
        rows,
        "19_figure3a_real_evidence",
        "Figure 3A evidence matrix is numeric and project-derived",
        "pass" if a_real and len(a_report.get("metric_provenance", {})) >= len(required_metrics) - 2 else "fail",
        {"matrix_rows": int(len(a_matrix)), "metric_provenance": a_report.get("metric_provenance")},
        "complete numeric metrics with provenance",
        a_matrix_path,
        "No manually assigned star ratings are used.",
    )
    architecture = str(a_report.get("architecture_language", "")).lower()
    causal_safe = "overlapping" in architecture and "partially ordered" in architecture
    add(
        rows,
        "20_noncausal_architecture",
        "Figure 3A does not draw a strict causal cascade",
        "pass" if causal_safe else "fail",
        architecture,
        "overlapping and partially ordered",
        a_report_path,
        "Figure annotation explicitly states that dashed connectors are not a proven linear cascade.",
    )

    checks = pd.DataFrame(rows)
    checks_path = out_dir / "figure3_egr1_a_to_f_validation_report.tsv"
    checks.to_csv(checks_path, sep="\t", index=False)
    n_fail = int(checks["status"].eq("fail").sum())
    n_warning = int(checks["status"].eq("warning").sum())
    aggregated_risks = []
    for path in [
        a_report_path,
        d_report_path,
        e_report_path,
        f_report_path,
        run_report_path,
        determinism_path,
        PROJECT_ROOT / "metadata/driver/figure3e_egr1_sensitivity/figure3e_egr1_sensitivity_report.json",
    ]:
        report = read_json(path)
        risks = report.get("review_risk_flags", [])
        if isinstance(risks, dict):
            risks = [risks]
        aggregated_risks.extend(risks)
    if changed:
        aggregated_risks.append(
            {
                "flag": "protected_asset_concurrent_change",
                "severity": "review_attention",
                "detail": (
                    f"{len(changed)} protected SOX4/HNF4A file(s) changed after the "
                    "preflight fingerprint. The Figure 3 workflow did not re-baseline or "
                    "revert them; see the protected-asset incident table."
                ),
            }
        )
    main_blocking = any(risk.get("severity") == "main_panel_blocking" for risk in aggregated_risks)
    standard = "meets_main_figure_standard" if n_fail == 0 and not main_blocking else "conditional_or_extended_data"
    report = {
        "module": "Figure 3 EGR1 A-F validation",
        "target_tf": TARGET_TF,
        "status": "fail" if n_fail else ("warning" if n_warning or main_blocking else "pass"),
        "n_checks": int(len(checks)),
        "n_pass": int(checks["status"].eq("pass").sum()),
        "n_warning": n_warning,
        "n_fail": n_fail,
        "sci_main_figure_assessment": standard,
        "formal_figure3e_generated": bool(e_report.get("formal_plot_generated")),
        "formal_figure3f_generated": bool(f_report.get("formal_plot_generated")),
        "review_risk_flags": aggregated_risks,
        "image_audit": image_results,
        "protected_assets_unchanged": not changed,
        "outputs": {"checks": str(checks_path.resolve())},
    }
    report_path = out_dir / "figure3_egr1_a_to_f_validation_report.json"
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
    return 0 if report["n_fail"] == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
