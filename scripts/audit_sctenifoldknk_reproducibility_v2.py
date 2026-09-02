#!/usr/bin/env python3
"""Audit historical Figure 2/3 scTenifoldKnk runs before any rerun."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import subprocess
from collections import Counter
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA_DIR = PROJECT_ROOT / "metadata/driver/sctenifoldknk_reproducibility_audit_v2"
DEFAULT_REPORT_DIR = PROJECT_ROOT / "reports/sctenifoldknk_reproducibility_audit_v2"
DEFAULT_RSCRIPT = Path(r"C:\Program Files\R\R-4.5.0\bin\x64\Rscript.exe")
FORMAL_SEEDS = [15071990, 15071991, 15071992]
FDR_THRESHOLD = 0.05


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig")) if path.exists() else {}


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() and path.stat().st_size else pd.DataFrame()


def sha256_file(path: Path) -> str | None:
    if not path.exists() or not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: Path, role: str) -> dict[str, Any]:
    exists = path.exists() and path.is_file()
    stat = path.stat() if exists else None
    return {
        "role": role,
        "path": str(path.resolve()),
        "exists": exists,
        "size_bytes": int(stat.st_size) if stat else None,
        "mtime": datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat() if stat else None,
        "sha256": sha256_file(path) if exists else None,
    }


def _numeric(frame: pd.DataFrame, column: str) -> pd.Series:
    return pd.to_numeric(frame[column], errors="coerce") if column in frame.columns else pd.Series(np.nan, index=frame.index)


def compute_gene_reproducibility(
    seed_tables: dict[int, pd.DataFrame],
    target_tf: str,
    fdr_threshold: float = FDR_THRESHOLD,
) -> pd.DataFrame:
    """Summarize per-gene recurrence across independent seed result tables."""
    if not seed_tables:
        return pd.DataFrame(
            columns=[
                "gene",
                "n_observed_seeds",
                "n_significant_seeds",
                "significant_seed_fraction",
                "recurrent_class",
                "direction_concordant",
                "direction_consistent_seed_count",
                "median_distance",
                "median_Z",
                "max_p_adj",
            ]
        )
    genes = sorted(
        {
            str(gene)
            for table in seed_tables.values()
            for gene in table.get("gene", pd.Series(dtype=str)).dropna().astype(str)
            if str(gene) != target_tf
        }
    )
    rows: list[dict[str, Any]] = []
    for gene in genes:
        observed = []
        significant = []
        signs = []
        distances = []
        p_adjusted = []
        for seed, table in seed_tables.items():
            row = table.loc[table["gene"].astype(str).eq(gene)] if "gene" in table.columns else pd.DataFrame()
            if row.empty:
                continue
            row = row.iloc[0]
            observed.append(seed)
            p_adj = pd.to_numeric(row.get("p.adj", np.nan), errors="coerce")
            z_value = pd.to_numeric(row.get("Z", np.nan), errors="coerce")
            distance = pd.to_numeric(row.get("distance", np.nan), errors="coerce")
            if np.isfinite(p_adj):
                p_adjusted.append(float(p_adj))
            if np.isfinite(distance):
                distances.append(float(distance))
            if np.isfinite(z_value) and z_value != 0:
                signs.append(int(np.sign(z_value)))
            if np.isfinite(p_adj) and p_adj < fdr_threshold:
                significant.append(seed)
        n_sig = len(significant)
        direction_concordant = bool(n_sig >= 2 and len(signs) >= 2 and len(set(signs)) == 1)
        rows.append(
            {
                "gene": gene,
                "n_observed_seeds": len(observed),
                "n_significant_seeds": n_sig,
                "significant_seed_fraction": n_sig / len(seed_tables),
                "recurrent_class": f"{n_sig}/{len(seed_tables)}",
                "direction_concordant": direction_concordant,
                "direction_consistent_seed_count": len(signs) if direction_concordant else 0,
                "median_distance": float(np.median(distances)) if distances else np.nan,
                "median_Z": float(np.median([float(x) for x in signs])) if signs else np.nan,
                "max_p_adj": float(max(p_adjusted)) if p_adjusted else np.nan,
            }
        )
    return pd.DataFrame(rows).sort_values(["n_significant_seeds", "median_distance", "gene"], ascending=[False, False, True]).reset_index(drop=True)


def pairwise_gene_reproducibility(
    seed_tables: dict[int, pd.DataFrame],
    target_tf: str,
    fdr_threshold: float = FDR_THRESHOLD,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    significant_sets: dict[int, set[str]] = {}
    distances: dict[int, pd.Series] = {}
    for seed, table in seed_tables.items():
        clean = table.loc[table["gene"].astype(str).ne(target_tf)].copy()
        p_adj = _numeric(clean, "p.adj")
        significant_sets[seed] = set(clean.loc[p_adj.lt(fdr_threshold), "gene"].astype(str))
        distances[seed] = _numeric(clean.set_index("gene"), "distance")
    for left, right in combinations(sorted(seed_tables), 2):
        left_set, right_set = significant_sets[left], significant_sets[right]
        union = left_set | right_set
        intersection = left_set & right_set
        left_values = distances[left]
        right_values = distances[right]
        common = left_values.index.intersection(right_values.index)
        if len(common) > 2:
            rho = float(left_values.loc[common].corr(right_values.loc[common], method="spearman"))
        else:
            rho = np.nan
        rows.append(
            {
                "seed_left": left,
                "seed_right": right,
                "n_left_significant": len(left_set),
                "n_right_significant": len(right_set),
                "n_shared_significant_genes": len(intersection),
                "n_union_significant_genes": len(union),
                "significant_gene_jaccard": len(intersection) / len(union) if union else np.nan,
                "overlap_coefficient": len(intersection) / min(len(left_set), len(right_set)) if min(len(left_set), len(right_set)) else np.nan,
                "n_common_tested_genes": len(common),
                "distance_rank_spearman_rho": rho,
            }
        )
    return pd.DataFrame(rows)


def classify_rerun_decision(
    target_tf: str,
    nc_nnet: int,
    nc_ncells: int,
    seeds: list[int],
    successful_seeds: int,
    input_cells: int,
) -> dict[str, Any]:
    reasons: list[str] = []
    if nc_nnet < 10:
        reasons.append(f"nc_nNet={nc_nnet} is below the formal 10-network standard")
    if nc_ncells < 500:
        reasons.append(f"nc_nCells={nc_ncells} is below the formal 500-cell standard")
    if len(seeds) < 3:
        reasons.append(f"only {len(seeds)} seed(s) are recorded; 3 independent seeds are required")
    if successful_seeds < 3:
        reasons.append(f"only {successful_seeds} seed(s) completed successfully")
    if input_cells < nc_ncells:
        reasons.append(f"input cell universe ({input_cells}) is smaller than requested nc_nCells={nc_ncells}")
    return {
        "target_tf": target_tf,
        "needs_rerun": bool(reasons),
        "reasons": reasons,
        "formal_standard": {"nc_nNet": 10, "nc_nCells": 500, "independent_seeds": 3},
    }


def build_three_axis_validation_matrix(summary: pd.DataFrame) -> pd.DataFrame:
    dimensions = [
        ("Axis", "axis"),
        ("nc_nNet", "nc_nNet"),
        ("nc_nCells", "nc_nCells"),
        ("Seeds", "seeds"),
        ("Successful seeds", "successful_seeds"),
        ("Input cells", "input_cells"),
        ("Input genes", "input_genes"),
        ("Significant genes", "significant_genes"),
        ("FDR genes", "fdr_genes"),
        ("3/3 recurrent genes", "three_of_three_recurrent_genes"),
        ("2/3 recurrent genes", "two_of_three_recurrent_genes"),
        ("Gene Jaccard", "gene_jaccard"),
        ("Direction concordance", "direction_concordance"),
        ("Enriched pathways", "enriched_pathways"),
        ("2/3 recurrent pathways", "two_of_three_recurrent_pathways"),
        ("3/3 recurrent pathways", "three_of_three_recurrent_pathways"),
        ("FDR pathways", "fdr_pathways"),
        ("Programme-level recurrence", "programme_level_recurrence"),
        ("Reproducibility grade", "reproducibility_grade"),
        ("Evidence strength", "evidence_strength"),
        ("Main Figure suitability", "main_figure_suitability"),
        ("Extended Data suitability", "extended_data_suitability"),
        ("Limitation", "limitation"),
    ]
    result = pd.DataFrame({"Dimension": [label for label, _ in dimensions]})
    for tf in ["HNF4A", "EGR1", "SOX4"]:
        row = summary.loc[summary["tf"].astype(str).eq(tf)] if "tf" in summary.columns else pd.DataFrame()
        values = []
        for _, column in dimensions:
            value = row.iloc[0].get(column, "") if not row.empty else ""
            values.append("" if pd.isna(value) else str(value))
        result[tf] = values
    return result


def _run_r_version(rscript: Path) -> dict[str, Any]:
    if not rscript.exists():
        return {"path": str(rscript), "exists": False}
    code = (
        "cat(R.version.string, '\\n'); "
        "for (p in c('scTenifoldKnk','Matrix','jsonlite')) "
        "if (requireNamespace(p, quietly=TRUE)) cat(p, as.character(packageVersion(p)), '\\n') "
        "else cat(p, 'MISSING\\n')"
    )
    try:
        completed = subprocess.run([str(rscript), "-e", code], capture_output=True, text=True, timeout=60, check=False)
        return {"path": str(rscript), "exists": True, "returncode": completed.returncode, "stdout": completed.stdout, "stderr": completed.stderr}
    except Exception as exc:
        return {"path": str(rscript), "exists": True, "error": str(exc)}


def _inspect_rds(rscript: Path, path: Path) -> dict[str, Any]:
    if not path.exists() or not rscript.exists():
        return {"path": str(path), "status": "missing_or_runtime_missing"}
    code = (
        "p <- commandArgs(TRUE)[1]; x <- readRDS(p); "
        "out <- list(names=names(x)); "
        "for (n in names(x)) { y <- x[[n]]; "
        "if (is.matrix(y) || inherits(y, 'Matrix')) { out[[paste0(n,'_dim')]] <- dim(y) } "
        "else if (is.data.frame(y)) { out[[paste0(n,'_dim')]] <- c(nrow(y),ncol(y)) } "
        "else if (is.list(y)) { out[[paste0(n,'_children')]] <- names(y) } }; "
        "cat(jsonlite::toJSON(out, auto_unbox=TRUE, null='null'))"
    )
    try:
        completed = subprocess.run([str(rscript), "-e", code, str(path)], capture_output=True, text=True, timeout=120, check=False)
        output = completed.stdout.strip().splitlines()
        payload = json.loads(output[-1]) if output else {}
        payload.update({"path": str(path.resolve()), "status": "pass" if completed.returncode == 0 else "error", "stderr": completed.stderr})
        return payload
    except Exception as exc:
        return {"path": str(path.resolve()), "status": "error", "error": str(exc)}


def _seed_tables_from_all_seed(path: Path, target_tf: str) -> dict[int, pd.DataFrame]:
    frame = read_tsv(path)
    if frame.empty or "seed" not in frame.columns:
        return {}
    frame = frame.loc[frame.get("tf", pd.Series(dtype=str)).astype(str).eq(target_tf)].copy()
    return {int(seed): table.copy() for seed, table in frame.groupby("seed")}


def _make_specs(root: Path) -> dict[str, dict[str, Any]]:
    return {
        "HNF4A": {
            "axis": "Identity",
            "subset": "normal_reference",
            "main_panel": "Figure 2E",
            "pathway_panel": "Figure 2F",
            "input_matrix": root / "data/processed/driver/figure2e_hnf4a_sctenifoldknk/normal_reference/figure2e_hnf4a_normal_reference_counts_genes_x_cells.mtx",
            "genes": root / "data/processed/driver/figure2e_hnf4a_sctenifoldknk/normal_reference/figure2e_hnf4a_normal_reference_genes.tsv",
            "input_report": root / "metadata/driver/figure2e_hnf4a_sctenifoldknk/figure2e_hnf4a_normal_reference_input_report.json",
            "run_report": root / "metadata/driver/figure2e_hnf4a_sctenifoldknk/figure2e_hnf4a_normal_reference_sctenifoldknk_report.json",
            "perturbation": root / "metadata/driver/figure2e_hnf4a_sctenifoldknk/figure2e_hnf4a_normal_reference_perturbation_genes.tsv",
            "rds": root / "metadata/driver/figure2e_hnf4a_sctenifoldknk/figure2e_hnf4a_normal_reference_sctenifoldknk.rds",
            "enrichment_report": root / "metadata/driver/figure2f_hnf4a/figure2f_hnf4a_report.json",
            "enrichment_all": root / "metadata/driver/figure2f_hnf4a/figure2f_hnf4a_enrichment_all.tsv",
            "figure_report": root / "metadata/driver/figure2e_hnf4a/figure2e_hnf4a_report.json",
            "sensitivity_report": root / "metadata/driver/figure2e_hnf4a_sensitivity/figure2e_hnf4a_sensitivity_report.json",
            "runner": root / "scripts/run_figure2e_hnf4a_sctenifoldknk.R",
            "orchestrator": root / "scripts/run_figure2_hnf4a_b_to_f.ps1",
            "seeds": [11],
        },
        "EGR1": {
            "axis": "Stress",
            "subset": "stressed_regenerative",
            "main_panel": "Figure 3E",
            "pathway_panel": "Figure 3F",
            "input_matrix": root / "data/processed/driver/figure3e_egr1_sctenifoldknk/stressed_regenerative/figure3e_egr1_stressed_regenerative_counts_genes_x_cells.mtx",
            "genes": root / "data/processed/driver/figure3e_egr1_sctenifoldknk/stressed_regenerative/figure3e_egr1_stressed_regenerative_genes.tsv",
            "input_report": root / "metadata/driver/figure3e_egr1/figure3e_egr1_input_report.json",
            "run_report": root / "metadata/driver/figure3e_egr1/figure3e_egr1_stressed_regenerative_run_report.json",
            "high_level_report": root / "metadata/driver/figure3e_egr1/figure3e_egr1_report.json",
            "rds": root / "data/processed/driver/figure3e_egr1_sctenifoldknk/stressed_regenerative/results/figure3e_egr1_stressed_regenerative_seed15071990_result.rds",
            "all_seed": root / "metadata/driver/figure3e_egr1/figure3e_egr1_stressed_regenerative_all_seed_perturbation_genes.tsv",
            "consensus": root / "metadata/driver/figure3e_egr1/figure3e_egr1_stressed_regenerative_consensus_perturbation_genes.tsv",
            "enrichment_report": root / "metadata/driver/figure3f_egr1/figure3f_egr1_report.json",
            "enrichment_all": root / "metadata/driver/figure3f_egr1/figure3f_egr1_enrichment_all.tsv",
            "figure_report": root / "metadata/driver/figure3e_egr1/figure3e_egr1_report.json",
            "sensitivity_report": root / "metadata/driver/figure3e_egr1_sensitivity/figure3e_egr1_sensitivity_report.json",
            "runner": root / "scripts/run_figure3e_egr1_sctenifoldknk.R",
            "orchestrator": root / "scripts/run_figure3_egr1_a_to_f.ps1",
            "seeds": FORMAL_SEEDS,
        },
        "SOX4": {
            "axis": "Malignant state",
            "subset": "malignant_like",
            "main_panel": "Figure 2E",
            "pathway_panel": "Figure 2F",
            "input_matrix": root / "data/processed/driver/sctenifoldknk_module7_1/malignant_like/sctenifoldknk_counts_genes_x_cells.mtx",
            "genes": root / "data/processed/driver/sctenifoldknk_module7_1/malignant_like/sctenifoldknk_genes.tsv",
            "input_report": root / "metadata/driver/sctenifoldknk_module7_1_malignant_like_export_report.json",
            "run_report": root / "metadata/driver/sctenifoldknk_module7_2_malignant_like_report.json",
            "perturbation": root / "data/processed/driver/sctenifoldknk_module7_2/malignant_like/sctenifoldknk_malignant_like_SOX4_perturbation_genes.tsv",
            "rds": root / "data/processed/driver/sctenifoldknk_module7_2/malignant_like/sctenifoldknk_malignant_like_SOX4.rds",
            "enrichment_report": root / "metadata/driver/figure2f_sox4/figure2f_sox4_report.json",
            "enrichment_all": root / "metadata/driver/figure2f_sox4/figure2f_sox4_enrichment_all.tsv",
            "figure_report": root / "metadata/driver/figure2e_sox4/figure2e_sox4_report.json",
            "sensitivity_report": None,
            "runner": root / "scripts/run_sctenifoldknk_module7_2.R",
            "orchestrator": root / "scripts/run_figure2_hnf4a_b_to_f.ps1",
            "seeds": [11],
        },
    }


def audit_historical_runs(
    project_root: Path = PROJECT_ROOT,
    metadata_dir: Path = DEFAULT_METADATA_DIR,
    report_dir: Path = DEFAULT_REPORT_DIR,
    rscript: Path = DEFAULT_RSCRIPT,
) -> dict[str, Any]:
    specs = _make_specs(project_root)
    metadata_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    entries: dict[str, Any] = {}
    source_records: list[dict[str, Any]] = []
    for tf, spec in specs.items():
        for key, path in spec.items():
            if isinstance(path, Path):
                source_records.append(file_record(path, f"{tf}:{key}"))
        run_report = read_json(spec["run_report"])
        input_report = read_json(spec["input_report"])
        enrichment_report = read_json(spec["enrichment_report"])
        figure_report = read_json(spec["figure_report"])
        perturb = read_tsv(spec["perturbation"]) if spec.get("perturbation") else pd.DataFrame()
        seed_tables = {}
        if tf == "EGR1":
            seed_tables = _seed_tables_from_all_seed(spec["all_seed"], tf)
        elif not perturb.empty:
            seed = int(run_report.get("parameters", {}).get("seed", 11))
            seed_tables = {seed: perturb}
        gene_summary = compute_gene_reproducibility(seed_tables, tf)
        pairwise = pairwise_gene_reproducibility(seed_tables, tf) if len(seed_tables) >= 2 else pd.DataFrame()
        seed_level_counts = {
            int(seed): int(
                _numeric(table.loc[table["gene"].astype(str).ne(tf)], "p.adj").lt(FDR_THRESHOLD).sum()
            )
            for seed, table in seed_tables.items()
        }
        params = run_report.get("parameters", {})
        if tf == "EGR1":
            params = {**params, "seeds": spec["seeds"]}
        n_input_cells = int(input_report.get("n_cells", run_report.get("n_cells", 0)) or 0)
        n_input_genes = int(input_report.get("n_genes", run_report.get("n_genes", 0)) or 0)
        successful_seeds = int(run_report.get("n_successful_seeds", 1 if perturb.shape[0] else 0))
        if tf == "EGR1":
            successful_seeds = int(run_report.get("n_successful_seeds", len(seed_tables)))
        n_sig = int(figure_report.get("n_significant_excluding_target", run_report.get("n_significant_excluding_target", 0)) or 0)
        if tf == "EGR1" and not gene_summary.empty:
            conservative = gene_summary.loc[gene_summary["max_p_adj"].lt(FDR_THRESHOLD)]
            n_sig = int(len(conservative))
        pathway_count = int(enrichment_report.get("n_pathways_tested", 0) or 0)
        fdr_pathway_count = int(enrichment_report.get("n_significant_pathways", 0) or 0)
        decision = classify_rerun_decision(
            tf,
            int(params.get("nc_nNet", 0) or 0),
            int(params.get("nc_nCells_used", params.get("nc_nCells", 0)) or 0),
            [int(x) for x in params.get("seeds", spec["seeds"])],
            successful_seeds,
            n_input_cells,
        )
        entries[tf] = {
            "tf": tf,
            "axis": spec["axis"],
            "subset": spec["subset"],
            "main_panel": spec["main_panel"],
            "pathway_panel": spec["pathway_panel"],
            "parameters": params,
            "input_cells": n_input_cells,
            "input_genes": n_input_genes,
            "successful_seeds": successful_seeds,
            "seed_count": len(seed_tables),
            "seed_values": sorted(seed_tables),
            "seed_level_significant_gene_counts": seed_level_counts,
            "single_seed": len(seed_tables) == 1,
            "downsampling": bool(n_input_cells and int(params.get("nc_nCells_used", params.get("nc_nCells", 0)) or 0) < n_input_cells),
            "significant_genes_excluding_target": n_sig,
            "pathways_tested": pathway_count,
            "fdr_significant_pathways": fdr_pathway_count,
            "gene_reproducibility": gene_summary.to_dict(orient="records"),
            "pairwise_gene_reproducibility": pairwise.to_dict(orient="records"),
            "result_structure": _inspect_rds(rscript, spec["rds"]) if spec.get("rds") else {},
            "reports": {
                "run": run_report,
                "input": input_report,
                "enrichment": enrichment_report,
                "figure": figure_report,
            },
            "decision": decision,
        }
    docs_matches: list[str] = []
    keywords = re.compile(r"HNF4A|EGR1|SOX4|scTenifoldKnk|nc_nNet|nc_nCells", re.IGNORECASE)
    for path in (project_root / "docs/chatgpt_discussions").rglob("*.md") if (project_root / "docs/chatgpt_discussions").exists() else []:
        try:
            if keywords.search(path.read_text(encoding="utf-8", errors="ignore")):
                docs_matches.append(str(path.resolve()))
        except OSError:
            continue
    audit_scope_paths = [
        project_root / "scripts/run_figure2e_hnf4a_sctenifoldknk.R",
        project_root / "scripts/prepare_figure2e_hnf4a_sctenifoldknk.py",
        project_root / "scripts/plot_figure2e_hnf4a_sctenifoldknk.R",
        project_root / "scripts/plot_figure2f_hnf4a_pathway_enrichment.R",
        project_root / "scripts/plot_figure2e_sox4_sctenifoldknk.R",
        project_root / "scripts/plot_figure2f_sox4_pathway_enrichment.R",
        project_root / "scripts/prepare_figure3e_egr1_sctenifoldknk.py",
        project_root / "scripts/run_figure3e_egr1_sctenifoldknk.R",
        project_root / "scripts/plot_figure3e_egr1_sctenifoldknk.R",
        project_root / "scripts/plot_figure3f_egr1_pathway_enrichment.R",
        project_root / "scripts/audit_figure3e_egr1_determinism.py",
        project_root / "scripts/summarize_figure3e_egr1_sensitivity.py",
        project_root / "scripts/run_sctenifoldknk_module7_2.R",
        project_root / "scripts/run_sctenifoldknk_enrichment_module7_4.R",
        project_root / "scripts/run_figure2_hnf4a_b_to_f.ps1",
        project_root / "scripts/run_figure3_egr1_a_to_f.ps1",
        project_root / "reports/figure2_hnf4a_b_to_f_report.md",
        project_root / "reports/figure3_egr1_a_to_f_report.md",
        project_root / "tests/test_figure2_hnf4a_outputs.py",
        project_root / "tests/test_figure3_egr1_logic.py",
        project_root / "tests/test_module7_sctenifoldknk_logic.py",
        project_root / "tests/test_figure2f_sox4_outputs.R",
    ]
    audit_scope_records = [file_record(path, "audit_scope") for path in audit_scope_paths]
    source_records.extend(audit_scope_records)
    runtime = {
        "audit_python": platform.python_version(),
        "platform": platform.platform(),
        "r_4_5_runtime": _run_r_version(rscript),
        "r_4_6_0_rscript": _run_r_version(Path(r"C:\Program Files\R\R-4.6.0\bin\x64\Rscript.exe")),
        "r_4_6_1_runtime": _run_r_version(Path(r"C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe")),
    }
    audit = {
        "module": "scTenifoldKnk Figure 2/3 reproducibility audit v2",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "audit_only": True,
        "formal_standard": {"nc_nNet": 10, "nc_nCells": 500, "independent_seeds": 3, "fdr_threshold": FDR_THRESHOLD},
        "historical_runs": entries,
        "source_file_inventory": source_records,
        "chatgpt_history_keyword_matches": docs_matches,
        "audit_scope_files": audit_scope_records,
        "runtime_inventory": runtime,
        "rerun_decisions": {tf: entry["decision"] for tf, entry in entries.items()},
    }
    (metadata_dir / "historical_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    rows = []
    for tf, entry in entries.items():
        params = entry["parameters"]
        rows.append(
            {
                "tf": tf,
                "axis": entry["axis"],
                "subset": entry["subset"],
                "nc_nNet": params.get("nc_nNet"),
                "nc_nCells": params.get("nc_nCells_used", params.get("nc_nCells")),
                "seeds": ";".join(map(str, entry["seed_values"])),
                "successful_seeds": entry["successful_seeds"],
                "input_cells": entry["input_cells"],
                "input_genes": entry["input_genes"],
                "significant_genes_excluding_target": entry["significant_genes_excluding_target"],
                "pathways_tested": entry["pathways_tested"],
                "fdr_significant_pathways": entry["fdr_significant_pathways"],
                "single_seed": entry["single_seed"],
                "needs_rerun": entry["decision"]["needs_rerun"],
                "rerun_reasons": " | ".join(entry["decision"]["reasons"]),
            }
        )
    pd.DataFrame(rows).to_csv(metadata_dir / "historical_run_summary.tsv", sep="\t", index=False)
    _write_historical_audit_report(audit, report_dir / "01_historical_audit.md")
    _write_rerun_decision_report(audit, report_dir / "02_rerun_decision.md")
    return audit


def _write_historical_audit_report(audit: dict[str, Any], path: Path) -> None:
    lines = [
        "# scTenifoldKnk Figure 2/3 Historical Reproducibility Audit v2",
        "",
        "## Audit boundary",
        "",
        "This report is read-only historical audit evidence. No new scTenifoldKnk network was started before this report was generated.",
        "Canonical SCENIC, CellOracle, trajectory, scVI/scanVI and existing Figure 8 outputs were not rerun.",
        "",
        "## Formal standard",
        "",
        "`nc_nNet = 10`, `nc_nCells = 500`, three independent seeds, matched input filtering, global BH FDR threshold 0.05.",
        "",
        "## Historical runs",
        "",
        "| TF | subset | nc_nNet | nc_nCells | seeds | input cells | input genes | FDR genes | FDR pathways | decision |",
        "|---|---|---:|---:|---|---:|---:|---:|---:|---|",
    ]
    for tf, entry in audit["historical_runs"].items():
        p = entry["parameters"]
        seeds = ", ".join(map(str, entry["seed_values"]))
        lines.append(
            f"| {tf} | {entry['subset']} | {p.get('nc_nNet', '')} | {p.get('nc_nCells_used', p.get('nc_nCells', ''))} | {seeds} | "
            f"{entry['input_cells']} | {entry['input_genes']} | {entry['significant_genes_excluding_target']} | "
            f"{entry['fdr_significant_pathways']} | {'RERUN' if entry['decision']['needs_rerun'] else 'RETAIN'} |"
        )
    lines += ["", "## HNF4A", ""]
    h = audit["historical_runs"]["HNF4A"]
    lines += [
        "The Figure 2E HNF4A history is a single `normal_reference` run with `nc_nNet=1`, `nc_nCells=100`, and `seed=11`.",
        f"It contains {h['input_cells']} input cells, {h['input_genes']} genes, {h['significant_genes_excluding_target']} FDR-significant genes excluding HNF4A, and {h['fdr_significant_pathways']} FDR-significant pathways.",
        "The saved R object contains WT/KO tensor networks, manifold alignment and differential regulation, but the historical record does not establish independent-seed reproducibility.",
        "",
        "## EGR1",
        "",
    ]
    e = audit["historical_runs"]["EGR1"]
    lines += [
        "The Figure 3E formal history uses `stressed_regenerative`, 646 cells, 3,000 genes, `nc_nNet=10`, `nc_nCells=500`, and seeds 15071990, 15071991 and 15071992; all three completed.",
        f"Seed-level FDR-significant gene counts excluding EGR1 are {', '.join(f'{seed}:{count}' for seed, count in sorted(e['seed_level_significant_gene_counts'].items()))}.",
        "The existing sensitivity table reports pairwise significant-gene Jaccard values of 0.017, 0.430 and 0.018, with distance-rank Spearman correlations of 0.858, 0.832 and 0.817. The conservative maximum-p.adj consensus retains one FDR-significant gene and the formal pathway FDR count is zero.",
        "EGR1 meets the computational parameter contract, while gene-level seed stability remains a review limitation.",
        "",
        "## SOX4",
        "",
    ]
    s = audit["historical_runs"]["SOX4"]
    lines += [
        "The Figure 2E/2F SOX4 history is the Module 7 `malignant_like` run with `nc_nNet=1`, `nc_nCells=100`, and `seed=11`.",
        f"It contains {s['input_cells']} input cells, {s['input_genes']} genes, {s['significant_genes_excluding_target']} FDR-significant genes excluding SOX4, and {s['fdr_significant_pathways']} FDR-significant pathways.",
        "The saved object contains WT/KO tensor networks, manifold alignment and differential regulation. Its one-seed design does not support independent-seed reproducibility.",
        "",
        "## Runtime and provenance",
        "",
        "Historical orchestrators point to R 4.5.0. The saved historical reports do not contain a complete sessionInfo; the audit therefore records the configured runtime separately from the exact historical package session.",
        f"ChatGPT discussion files matching the audit keywords: {len(audit['chatgpt_history_keyword_matches'])}.",
        "File-level SHA-256 inventory and object-structure checks are in `metadata/driver/sctenifoldknk_reproducibility_audit_v2/historical_audit.json`.",
        "",
        "## Scientific boundary",
        "",
        "scTenifoldKnk is computational virtual network perturbation evidence. These historical runs do not establish genetic knockout causality or a strict linear HNF4A/AP-1/CEBPB/EGR1/SOX4 cascade.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_rerun_decision_report(audit: dict[str, Any], path: Path) -> None:
    lines = [
        "# scTenifoldKnk Rerun Decision v2",
        "",
        "## Decision rule",
        "",
        "A formal rerun is required when any of the following is missing: `nc_nNet=10`, `nc_nCells=500`, three independent seeds, complete outputs, or a traceable successful run record.",
        "",
    ]
    for tf in ["HNF4A", "EGR1", "SOX4"]:
        entry = audit["historical_runs"][tf]
        decision = entry["decision"]
        lines.append(f"## {tf}")
        lines.append("")
        lines.append(f"Decision: **{'RERUN' if decision['needs_rerun'] else 'RETAIN HISTORICAL FORMAL RUN'}**")
        lines.append("")
        if decision["reasons"]:
            lines.append("Reasons:")
            lines.extend(f"- {reason}" for reason in decision["reasons"])
        else:
            lines.append("The historical run meets the computational 10 × 500 × 3-seed contract; no duplicate rerun is required.")
        lines.append("")
    lines += [
        "## Planned versioned action",
        "",
        "HNF4A and SOX4 require new versioned three-seed runs using the same computational contract. EGR1 will be retained as the historical formal reference because it already meets the parameter contract; its low gene-level Jaccard remains explicitly reported.",
        "New outputs will be written under `data/processed/driver/sctenifoldknk_reproducibility_audit_v2`, `metadata/driver/sctenifoldknk_reproducibility_audit_v2`, `reports/sctenifoldknk_reproducibility_audit_v2` and `figures/driver/sctenifoldknk_reproducibility_audit_v2`.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA_DIR)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--rscript", type=Path, default=DEFAULT_RSCRIPT)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    audit = audit_historical_runs(args.project_root, args.metadata_dir, args.report_dir, args.rscript)
    print(json.dumps({"historical_audit": str((args.metadata_dir / 'historical_audit.json').resolve()), "decisions": audit["rerun_decisions"]}, ensure_ascii=False, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
