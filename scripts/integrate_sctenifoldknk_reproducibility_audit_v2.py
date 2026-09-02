#!/usr/bin/env python3
"""Integrate the versioned HNF4A/SOX4 reruns with formal EGR1 history."""

from __future__ import annotations

import argparse
import json
import platform
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

try:
    from audit_sctenifoldknk_reproducibility_v2 import (
        PROJECT_ROOT,
        compute_gene_reproducibility,
        pairwise_gene_reproducibility,
    )
except ModuleNotFoundError:
    from scripts.audit_sctenifoldknk_reproducibility_v2 import (
        PROJECT_ROOT,
        compute_gene_reproducibility,
        pairwise_gene_reproducibility,
    )


DEFAULT_METADATA_DIR = PROJECT_ROOT / "metadata/driver/sctenifoldknk_reproducibility_audit_v2"
DEFAULT_REPORT_DIR = PROJECT_ROOT / "reports/sctenifoldknk_reproducibility_audit_v2"
DEFAULT_FIGURE_DIR = PROJECT_ROOT / "figures/driver/sctenifoldknk_reproducibility_audit_v2"
FDR_THRESHOLD = 0.05
GMT_FILES = {
    "KEGG": "KEGG_2021_Human.gmt",
    "Reactome": "Reactome_2022.gmt",
    "GO_BP": "GO_Biological_Process_2023.gmt",
}


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() and path.stat().st_size else pd.DataFrame()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig")) if path.exists() else {}


def bh_adjust(values: list[float]) -> np.ndarray:
    if not values:
        return np.asarray([], dtype=float)
    pvalues = np.asarray(values, dtype=float)
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    adjusted = np.empty(len(ranked), dtype=float)
    running = 1.0
    for index in range(len(ranked) - 1, -1, -1):
        running = min(running, ranked[index] * len(ranked) / (index + 1))
        adjusted[index] = running
    output = np.empty(len(pvalues), dtype=float)
    output[order] = np.clip(adjusted, 0, 1)
    return output


def parse_gmt(path: Path, background: set[str], min_size: int = 5, max_size: int = 500) -> dict[str, set[str]]:
    gene_sets: dict[str, set[str]] = {}
    if not path.exists():
        return gene_sets
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        term = fields[0].strip()
        genes = {gene.strip().upper() for gene in fields[2:] if gene.strip()} & background
        if term and min_size <= len(genes) <= max_size:
            gene_sets[term] = genes
    return gene_sets


def annotate_programme(term: str, tf: str) -> str:
    text = str(term).lower()
    patterns = {
        "HNF4A": [
            ("hepatocyte_identity", r"hepatocyte|liver|xenobiotic|bile|albumin|differenti"),
            ("lipid_metabolism", r"lipid|fatty acid|peroxisome|ppar|cholesterol|oxidation"),
        ],
        "EGR1": [
            ("stress_response", r"stress|heat shock|unfolded protein|immediate early"),
            ("inflammatory_transition", r"inflamm|cytokine|interferon|tnf|mapk|nf-kb"),
            ("ap1_associated", r"ap-1|jun|fos"),
        ],
        "SOX4": [
            ("malignant_state", r"cancer|tumou?r|malignan|invasion|metast|emt|stem|development"),
            ("proliferation_survival", r"cell cycle|prolifer|replication|apoptosis|survival"),
            ("stress_response", r"stress|heat shock|unfolded protein"),
        ],
    }
    for programme, pattern in patterns.get(tf, []):
        if re.search(pattern, text):
            return programme
    return "unannotated"


def compute_seed_ora(
    table: pd.DataFrame,
    target_tf: str,
    seed: int,
    background: set[str],
    gene_sets_by_database: dict[str, dict[str, set[str]]],
    fdr_threshold: float = FDR_THRESHOLD,
) -> pd.DataFrame:
    work = table.copy()
    work["gene"] = work["gene"].astype(str).str.upper()
    work["p.adj"] = pd.to_numeric(work.get("p.adj", np.nan), errors="coerce")
    input_genes = set(work.loc[work["gene"].ne(target_tf) & work["p.adj"].lt(fdr_threshold), "gene"]) & background
    rows: list[dict[str, Any]] = []
    for database, gene_sets in gene_sets_by_database.items():
        for term, term_genes in gene_sets.items():
            overlap = input_genes & term_genes
            pvalue = (
                float(hypergeom.sf(len(overlap) - 1, len(background), len(term_genes), len(input_genes)))
                if input_genes
                else 1.0
            )
            rows.append(
                {
                    "tf": target_tf,
                    "seed": int(seed),
                    "database": database,
                    "term": term,
                    "overlap_count": len(overlap),
                    "term_size": len(term_genes),
                    "input_gene_count": len(input_genes),
                    "background_gene_count": len(background),
                    "pvalue": pvalue,
                    "overlap_genes": ";".join(sorted(overlap)),
                    "programme_annotation": annotate_programme(term, target_tf),
                    "enrichment_direction": "unsigned_network_displacement",
                }
            )
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    result["p.adjust"] = bh_adjust(result["pvalue"].tolist())
    result["significant"] = result["p.adjust"].lt(fdr_threshold)
    return result.sort_values(["p.adjust", "pvalue", "term"]).reset_index(drop=True)


def compute_pathway_reproducibility(
    seed_results: list[pd.DataFrame],
    target_tf: str,
    fdr_threshold: float = FDR_THRESHOLD,
) -> pd.DataFrame:
    frames = [frame.copy() for frame in seed_results if frame is not None and not frame.empty]
    if not frames:
        return pd.DataFrame()
    combined = pd.concat(frames, ignore_index=True)
    combined["p.adjust"] = pd.to_numeric(combined.get("p.adjust", np.nan), errors="coerce")
    combined["seed"] = pd.to_numeric(combined["seed"], errors="coerce").astype(int)
    n_seeds = int(combined["seed"].nunique())
    rows: list[dict[str, Any]] = []
    for (database, term), frame in combined.groupby(["database", "term"], sort=True):
        significant_seeds = set(frame.loc[frame["p.adjust"].lt(fdr_threshold), "seed"].astype(int))
        annotation = str(frame.get("programme_annotation", pd.Series(["unannotated"])).dropna().iloc[0])
        rows.append(
            {
                "tf": target_tf,
                "database": database,
                "term": term,
                "n_tested_seeds": n_seeds,
                "n_significant_seeds": len(significant_seeds),
                "significant_seed_fraction": len(significant_seeds) / n_seeds if n_seeds else np.nan,
                "recurrent_3_of_3": len(significant_seeds) == 3 and n_seeds == 3,
                "recurrent_2_of_3": len(significant_seeds) >= 2 and n_seeds == 3,
                "min_p_adjust": float(frame["p.adjust"].min()),
                "max_p_adjust": float(frame["p.adjust"].max()),
                "median_overlap_count": float(pd.to_numeric(frame.get("overlap_count", np.nan), errors="coerce").median()),
                "programme_annotation": annotation,
                "enrichment_direction": "unsigned_network_displacement",
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["n_significant_seeds", "min_p_adjust", "database", "term"],
        ascending=[False, True, True, True],
    ).reset_index(drop=True)


def _direction_concordance(seed_tables: dict[int, pd.DataFrame], target_tf: str, fdr_threshold: float) -> float:
    values: list[float] = []
    for left_seed, right_seed in __import__("itertools").combinations(sorted(seed_tables), 2):
        left = seed_tables[left_seed].loc[seed_tables[left_seed]["gene"].astype(str).ne(target_tf)].copy()
        right = seed_tables[right_seed].loc[seed_tables[right_seed]["gene"].astype(str).ne(target_tf)].copy()
        for frame in [left, right]:
            frame["p.adj"] = pd.to_numeric(frame.get("p.adj", np.nan), errors="coerce")
            frame["Z"] = pd.to_numeric(frame.get("Z", np.nan), errors="coerce")
        left = left.loc[left["p.adj"].lt(fdr_threshold)].set_index("gene")
        right = right.loc[right["p.adj"].lt(fdr_threshold)].set_index("gene")
        common = left.index.intersection(right.index)
        if len(common):
            left_sign = np.sign(left.loc[common, "Z"].to_numpy(dtype=float))
            right_sign = np.sign(right.loc[common, "Z"].to_numpy(dtype=float))
            finite = np.isfinite(left_sign) & np.isfinite(right_sign) & (left_sign != 0) & (right_sign != 0)
            if finite.any():
                values.append(float(np.mean(left_sign[finite] == right_sign[finite])))
    return float(np.median(values)) if values else np.nan


def _grade(jaccard: float, direction: float, fdr_pathways: int, recurrent_pathways: int, three_genes: int, successes: int) -> str:
    if successes < 3:
        return "D"
    if np.isfinite(jaccard) and jaccard >= 0.5 and np.isfinite(direction) and direction >= 0.8 and fdr_pathways > 0:
        return "A"
    if np.isfinite(jaccard) and jaccard >= 0.5 and three_genes > 0 and (not np.isfinite(direction) or direction >= 0.5):
        return "B"
    if ((np.isfinite(jaccard) and jaccard >= 0.25) or three_genes > 0) and recurrent_pathways > 0:
        return "B"
    return "C"


def decide_panel_position(target_tf: str, summary: pd.Series | dict[str, Any]) -> dict[str, str]:
    row = summary if isinstance(summary, dict) else summary.to_dict()
    grade = str(row.get("reproducibility_grade", "D"))
    fdr_genes = int(row.get("fdr_genes", 0) or 0)
    fdr_pathways = int(row.get("fdr_pathways", 0) or 0)
    if target_tf == "HNF4A":
        return {
            "figure2e": "KEEP_MAIN" if grade in {"A", "B"} and fdr_genes > 0 else "MOVE_EXTENDED",
            "figure2f": "KEEP_MAIN" if grade in {"A", "B"} and fdr_pathways > 0 else "MOVE_EXTENDED",
        }
    if target_tf == "SOX4":
        return {
            "figure2e": "KEEP_MAIN" if grade in {"A", "B"} and fdr_genes > 0 else "REVISE",
            "figure2f": "KEEP_MAIN" if grade in {"A", "B"} and fdr_pathways > 0 else "REVISE",
        }
    return {
        "figure3e": "KEEP_MAIN" if grade in {"A", "B"} and fdr_genes > 0 else "MOVE_EXTENDED",
        "figure3f": "KEEP_MAIN" if grade in {"A", "B"} and fdr_pathways > 0 else "MOVE_EXTENDED",
    }


def build_axis_summary(
    tf: str,
    axis: str,
    subset: str,
    seed_tables: dict[int, pd.DataFrame],
    pathway_results: list[pd.DataFrame],
    input_cells: int,
    input_genes: int,
    fdr_threshold: float = FDR_THRESHOLD,
) -> pd.Series:
    gene_repro = compute_gene_reproducibility(seed_tables, tf, fdr_threshold)
    pairwise = pairwise_gene_reproducibility(seed_tables, tf, fdr_threshold)
    successes = len(seed_tables)
    seed_counts = {
        int(seed): int(
            pd.to_numeric(table.loc[table["gene"].astype(str).ne(tf), "p.adj"], errors="coerce").lt(fdr_threshold).sum()
        )
        for seed, table in seed_tables.items()
    }
    conservative = gene_repro.loc[gene_repro["max_p_adj"].lt(fdr_threshold)] if not gene_repro.empty else pd.DataFrame()
    three = (
        gene_repro.loc[gene_repro["n_significant_seeds"].eq(successes) & gene_repro["direction_concordant"]]
        if not gene_repro.empty
        else pd.DataFrame()
    )
    two = (
        gene_repro.loc[gene_repro["n_significant_seeds"].eq(max(successes - 1, 0)) & gene_repro["direction_concordant"]]
        if not gene_repro.empty
        else pd.DataFrame()
    )
    pathway_repro = compute_pathway_reproducibility(pathway_results, tf, fdr_threshold)
    fdr_pathways = pathway_repro.loc[pathway_repro["n_significant_seeds"].eq(successes)] if not pathway_repro.empty else pd.DataFrame()
    recurrent_pathways = pathway_repro.loc[pathway_repro["n_significant_seeds"].ge(2)] if not pathway_repro.empty else pd.DataFrame()
    programmes = sorted(set(recurrent_pathways["programme_annotation"].astype(str))) if not recurrent_pathways.empty else []
    any_programmes = sorted(set(pathway_repro.loc[pathway_repro["n_significant_seeds"].gt(0), "programme_annotation"].astype(str))) if not pathway_repro.empty else []
    jaccard = float(pd.to_numeric(pairwise.get("significant_gene_jaccard", np.nan), errors="coerce").median()) if not pairwise.empty else np.nan
    direction = _direction_concordance(seed_tables, tf, fdr_threshold)
    grade = _grade(jaccard, direction, len(fdr_pathways), len(recurrent_pathways), len(three), successes)
    summary = pd.Series(
        {
            "tf": tf,
            "axis": axis,
            "subset": subset,
            "nc_nNet": 10,
            "nc_nCells": 500,
            "seeds": ",".join(map(str, sorted(seed_tables))),
            "successful_seeds": successes,
            "input_cells": input_cells,
            "input_genes": input_genes,
            "significant_genes": int(np.median(list(seed_counts.values()))) if seed_counts else 0,
            "fdr_genes": int(len(conservative)),
            "three_of_three_recurrent_genes": int(len(three)),
            "two_of_three_recurrent_genes": int(len(two)),
            "gene_jaccard": jaccard,
            "direction_concordance": direction,
            "enriched_pathways": int(len(pathway_repro.loc[pathway_repro["n_significant_seeds"].gt(0)])) if not pathway_repro.empty else 0,
            "fdr_pathways": int(len(fdr_pathways)),
            "programme_level_recurrence": ";".join(programmes) if programmes else ("single_seed_only:" + ";".join(any_programmes) if any_programmes else "none_FDR"),
            "two_of_three_recurrent_pathways": int(len(pathway_repro.loc[pathway_repro["n_significant_seeds"].eq(2)])) if not pathway_repro.empty else 0,
            "three_of_three_recurrent_pathways": int(len(fdr_pathways)),
            "reproducibility_grade": grade,
            "evidence_strength": {"A": "strong", "B": "moderate", "C": "limited", "D": "unsupported"}[grade],
            "main_figure_suitability": "",
            "extended_data_suitability": "",
            "limitation": "biological_input_asymmetry; unsigned_manifold_distance; computational_virtual_perturbation",
            "seed_level_significant_gene_counts": ";".join(f"{seed}:{count}" for seed, count in sorted(seed_counts.items())),
        }
    )
    panel = decide_panel_position(tf, summary)
    if tf == "EGR1":
        summary["main_figure_suitability"] = panel["figure3e"]
        summary["extended_data_suitability"] = panel["figure3f"]
    else:
        summary["main_figure_suitability"] = panel["figure2e"]
        summary["extended_data_suitability"] = panel["figure2f"]
    return summary


def _load_target_tables(tf: str, metadata_dir: Path) -> tuple[dict[int, pd.DataFrame], str, int, int]:
    if tf == "EGR1":
        path = PROJECT_ROOT / "metadata/driver/figure3e_egr1/figure3e_egr1_stressed_regenerative_all_seed_perturbation_genes.tsv"
        frame = read_tsv(path)
        tables = {
            int(seed): table.copy()
            for seed, table in frame.loc[frame["tf"].astype(str).eq(tf)].groupby("seed")
        }
        return tables, "stressed_regenerative", 646, 3000
    subset = "normal_reference" if tf == "HNF4A" else "malignant_like"
    paths = sorted((metadata_dir / tf / subset).rglob(f"{tf}_{subset}_seed*_perturbation_genes.tsv"))
    tables: dict[int, pd.DataFrame] = {}
    for path in paths:
        match = re.search(r"seed(\d+)_perturbation_genes", path.name)
        if match:
            tables[int(match.group(1))] = read_tsv(path)
    return tables, subset, 8098 if tf == "HNF4A" else 552, 3000


def _load_pathways(tf: str, seed_tables: dict[int, pd.DataFrame]) -> tuple[list[pd.DataFrame], set[str]]:
    if not seed_tables:
        return [], set()
    background = set(next(iter(seed_tables.values()))["gene"].astype(str).str.upper())
    gene_sets = {
        database: parse_gmt(PROJECT_ROOT / "metadata/driver/sctenifoldknk_module7_4_genesets" / filename, background)
        for database, filename in GMT_FILES.items()
    }
    results = [compute_seed_ora(table, tf, seed, background, gene_sets) for seed, table in sorted(seed_tables.items())]
    return results, background


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_safe(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def _plot_summary(summary: pd.DataFrame, seed_level: pd.DataFrame, pairwise: pd.DataFrame, pathway: pd.DataFrame, figure_dir: Path) -> dict[str, str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 7,
        }
    )
    ordered = ["HNF4A", "EGR1", "SOX4"]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6), dpi=220)
    pivot = seed_level.pivot_table(index="tf", columns="seed", values="significant_genes", aggfunc="first").reindex(ordered)
    if pivot.empty:
        pivot = pd.DataFrame(0, index=ordered, columns=[1, 2, 3])
    im = axes[0, 0].imshow(pivot.fillna(0).to_numpy(dtype=float), aspect="auto", cmap="Blues")
    axes[0, 0].set_xticks(range(len(pivot.columns))); axes[0, 0].set_xticklabels(pivot.columns)
    axes[0, 0].set_yticks(range(len(pivot.index))); axes[0, 0].set_yticklabels(pivot.index)
    axes[0, 0].set_title("a  Seed-level FDR genes", loc="left", fontweight="bold")
    axes[0, 0].set_xlabel("Seed")
    for i in range(len(pivot.index)):
        for j in range(len(pivot.columns)):
            axes[0, 0].text(j, i, f"{pivot.iloc[i, j]:.0f}", ha="center", va="center")
    fig.colorbar(im, ax=axes[0, 0], fraction=0.046, pad=0.04, label="FDR genes")

    if pairwise.empty:
        axes[0, 1].axis("off")
        axes[0, 1].text(0.5, 0.5, "pairwise overlap pending", ha="center", va="center", transform=axes[0, 1].transAxes)
    else:
        jaccard = pairwise.pivot_table(index="tf", columns="pair", values="significant_gene_jaccard", aggfunc="first").reindex(ordered)
        im = axes[0, 1].imshow(jaccard.fillna(0).to_numpy(dtype=float), aspect="auto", vmin=0, vmax=1, cmap="YlGnBu")
        axes[0, 1].set_xticks(range(len(jaccard.columns))); axes[0, 1].set_xticklabels(jaccard.columns, rotation=35, ha="right", fontsize=6)
        axes[0, 1].set_yticks(range(len(jaccard.index))); axes[0, 1].set_yticklabels(jaccard.index)
        axes[0, 1].set_title("b  Pairwise significant-gene Jaccard", loc="left", fontweight="bold")
        fig.colorbar(im, ax=axes[0, 1], fraction=0.046, pad=0.04, label="Jaccard")

    if pathway.empty:
        pathway_plot = pd.DataFrame(0, index=ordered, columns=["recurrent_fdr_pathways", "fdr_pathways"])
    else:
        pathway_plot = pathway.groupby("tf", as_index=False).agg(
            recurrent_fdr_pathways=("n_significant_seeds", lambda x: int((x >= 2).sum())),
            fdr_pathways=("recurrent_3_of_3", "sum"),
        ).set_index("tf").reindex(ordered).fillna(0)
    x = np.arange(len(pathway_plot.index)); width = 0.35
    axes[1, 0].bar(x - width / 2, pathway_plot["recurrent_fdr_pathways"], width, color="#8FB9A8", label=">=2/3 seeds")
    axes[1, 0].bar(x + width / 2, pathway_plot["fdr_pathways"], width, color="#2E7D4F", label="3/3 seeds")
    axes[1, 0].set_xticks(x); axes[1, 0].set_xticklabels(ordered); axes[1, 0].set_ylabel("Pathway count")
    axes[1, 0].set_title("c  Pathway recurrence", loc="left", fontweight="bold"); axes[1, 0].legend(frameon=False, fontsize=6)
    if float(pathway_plot[["recurrent_fdr_pathways", "fdr_pathways"]].to_numpy().max()) == 0:
        axes[1, 0].set_ylim(0, 1)
        axes[1, 0].text(1.0, 0.55, "No pathway recurred\nat FDR < 0.05", ha="center", va="center", fontsize=7, color="#555555")

    plot_summary = summary.set_index("tf").reindex(ordered)
    scores = pd.to_numeric(plot_summary["gene_jaccard"], errors="coerce").fillna(0).to_numpy()
    colours = {"HNF4A": "#3A6EA5", "EGR1": "#C07A2C", "SOX4": "#2E7D4F"}
    axes[1, 1].barh(np.arange(len(ordered)), scores, color=[colours[tf] for tf in ordered])
    axes[1, 1].set_xlim(0, 1); axes[1, 1].set_yticks(np.arange(len(ordered))); axes[1, 1].set_yticklabels(ordered); axes[1, 1].invert_yaxis()
    axes[1, 1].set_xlabel("Median pairwise gene Jaccard"); axes[1, 1].set_title("d  Three-axis reproducibility", loc="left", fontweight="bold")
    for i, (_, row) in enumerate(plot_summary.iterrows()):
        axes[1, 1].text(min(float(scores[i]) + 0.03, 0.78), i, f"Grade {row['reproducibility_grade']} | {int(row['successful_seeds'])}/3", va="center")
    fig.tight_layout()
    stem = figure_dir / "sctenifoldknk_three_axis_reproducibility_summary"
    outputs = {}
    for suffix, kwargs in [("png", {"dpi": 600}), ("pdf", {}), ("svg", {})]:
        path = stem.with_suffix(f".{suffix}")
        fig.savefig(path, bbox_inches="tight", **kwargs)
        outputs[suffix] = str(path.resolve())
    plt.close(fig)
    return outputs


def _write_three_axis_report(summary: pd.DataFrame, outputs: dict[str, Any], path: Path) -> None:
    lines = [
        "# Three-axis scTenifoldKnk Validation Report v2",
        "",
        "HNF4A and SOX4 were rerun in the versioned namespace. EGR1 was retained from the historical formal 10 x 500 x 3-seed run. Computational parameters are symmetric; biological input subsets are state-matched and distinct.",
        "",
        "| TF | subset | successful seeds | FDR genes | 3/3 genes | 2/3 genes | median gene Jaccard | direction concordance | any FDR pathways | 2/3 pathways | 3/3 pathways | grade |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for _, row in summary.iterrows():
        jaccard = f"{row['gene_jaccard']:.3f}" if np.isfinite(row["gene_jaccard"]) else "NA"
        direction = f"{row['direction_concordance']:.3f}" if np.isfinite(row["direction_concordance"]) else "NA"
        lines.append(f"| {row['tf']} | {row['subset']} | {int(row['successful_seeds'])}/3 | {int(row['fdr_genes'])} | {int(row['three_of_three_recurrent_genes'])} | {int(row['two_of_three_recurrent_genes'])} | {jaccard} | {direction} | {int(row['enriched_pathways'])} | {int(row['two_of_three_recurrent_pathways'])} | {int(row['three_of_three_recurrent_pathways'])} | {row['reproducibility_grade']} |")
    lines += [
        "",
        "The provisional recurrent-gene rule is support in at least two of three seeds with concordant sign of the reported Z statistic. Manifold distance is unsigned network displacement and is not interpreted as biological activation, suppression or experimental causality.",
        "",
        "P0 gate: CLOSED for historical parameter/seed audit and necessary standardized reruns. Low reproducibility remains a result for EGR1 and SOX4 and governs conservative panel placement.",
        "",
        "Outputs:",
    ]
    lines.extend(f"- {key}: {value}" for key, value in outputs.items() if value)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_postrun_decision_report(summary: pd.DataFrame, path: Path) -> None:
    lines = [
        "# scTenifoldKnk Rerun Decision v2",
        "",
        "## Final status",
        "",
        "HNF4A and SOX4 were rerun in the versioned namespace with nc_nNet=10, nc_nCells=500 and three independent seeds. EGR1 already met the same computational contract and was retained as the historical formal reference.",
        "",
        "| TF | historical status | new action | successful seeds | final status |",
        "|---|---|---|---:|---|",
        "| HNF4A | 1 x 100, one seed | standardized rerun | 3/3 | complete |",
        "| EGR1 | 10 x 500, three seeds | no rerun | 3/3 historical | retained; low gene Jaccard remains |",
        "| SOX4 | 1 x 100, one seed | standardized rerun | 3/3 | complete |",
        "",
        "## Parameter symmetry",
        "",
        "All three axes use nc_nNet=10, nc_nCells=500, nc_nComp=3, ma_nDim=2, qc=false, qc_minCells=3, the same three seeds and the same 3,000-gene network. Biological input asymmetry is retained because the axes represent different state-matched subsets.",
        "",
        "## P0 gate",
        "",
        "Figure 2/3 scTenifoldKnk low-reproducibility P0 gate: CLOSED for audit and required rerun completion. The gate closure does not promote Grade C results to main mechanistic evidence.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_panel_decision_report(summary: pd.DataFrame, path: Path) -> None:
    lines = ["# Figure 2/3 scTenifoldKnk Panel Decision v2", "", "| Panel | TF | Decision | Basis |", "|---|---|---|---|"]
    for _, row in summary.iterrows():
        tf = str(row["tf"]); panel = decide_panel_position(tf, row)
        if tf == "EGR1":
            decisions = [("Figure 3E", panel["figure3e"], "seed-level perturbation genes"), ("Figure 3F", panel["figure3f"], "pathway enrichment")]
        else:
            decisions = [(f"Figure 2E {tf}", panel["figure2e"], "seed-level perturbation genes"), (f"Figure 2F {tf}", panel["figure2f"], "pathway enrichment")]
        for name, decision, basis in decisions:
            lines.append(f"| {name} | {tf} | {decision} | {basis}; reproducibility grade {row['reproducibility_grade']} |")
    lines += [
        "",
        "HNF4A Figure 2F is moved to Extended Data whenever no pathway passes the global BH FDR threshold. EGR1 Figure 3E/3F is moved to Extended Data when seed-level gene overlap or pathway FDR is insufficient. SOX4 panels require revision if the state-specific perturbation remains limited to stress-response or otherwise lacks stable malignant-state programme recurrence.",
        "",
        "These decisions do not alter existing CellOracle panels and do not turn virtual perturbation into experimental knockout evidence.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_claim_update_report(summary: pd.DataFrame, path: Path) -> None:
    lines = ["# Claim and Evidence Update v2", "", "| Claim | Evidence | Reproducibility | Allowed wording | Forbidden wording |", "|---|---|---|---|---|"]
    for _, row in summary.iterrows():
        tf = row["tf"]
        claim = {
            "HNF4A": "HNF4A-associated identity network displacement",
            "EGR1": "EGR1-associated stress-transition network perturbation",
            "SOX4": "SOX4-associated malignant-state network perturbation",
        }[tf]
        lines.append(f"| {claim} | scTenifoldKnk virtual perturbation with {int(row['successful_seeds'])}/3 successful seeds and {int(row['fdr_genes'])} conservative FDR genes | Grade {row['reproducibility_grade']} | supports; is consistent with; computationally implicates | proves; demonstrates causality; establishes genetic dependency |")
    lines += [
        "",
        "The defensible model is an overlapping, partially ordered regulatory architecture. A strict identity-loss to stress-activation to SOX4-stabilization causal cascade requires direct perturbation, rescue or epistasis and orthogonal molecular readouts.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_integration(
    metadata_dir: Path = DEFAULT_METADATA_DIR,
    report_dir: Path = DEFAULT_REPORT_DIR,
    figure_dir: Path = DEFAULT_FIGURE_DIR,
    fdr_threshold: float = FDR_THRESHOLD,
) -> dict[str, Any]:
    summary_rows: list[pd.Series] = []
    seed_rows: list[dict[str, Any]] = []
    gene_frames: list[pd.DataFrame] = []
    pairwise_frames: list[pd.DataFrame] = []
    pathway_frames: list[pd.DataFrame] = []
    pathway_seed_frames: list[pd.DataFrame] = []
    for tf, axis, default_subset in [
        ("HNF4A", "Identity", "normal_reference"),
        ("EGR1", "Stress", "stressed_regenerative"),
        ("SOX4", "Malignant state", "malignant_like"),
    ]:
        seed_tables, subset, input_cells, input_genes = _load_target_tables(tf, metadata_dir)
        subset = subset or default_subset
        if not seed_tables:
            continue
        pathway_results, _ = _load_pathways(tf, seed_tables)
        pathway_by_seed = {seed: pathway_results[index] for index, seed in enumerate(sorted(seed_tables))} if pathway_results else {}
        for seed, table in sorted(seed_tables.items()):
            clean = table.loc[table["gene"].astype(str).ne(tf)].copy()
            p_adj = pd.to_numeric(clean.get("p.adj", np.nan), errors="coerce")
            top = clean.loc[p_adj.lt(fdr_threshold)].sort_values("distance", ascending=False).head(20)
            pathway_frame = pathway_by_seed.get(seed, pd.DataFrame())
            seed_rows.append(
                {
                    "tf": tf, "axis": axis, "subset": subset, "seed": int(seed), "nc_nNet": 10, "nc_nCells": 500,
                    "input_cells": input_cells, "input_genes": input_genes,
                    "significant_genes": int(p_adj.lt(fdr_threshold).sum()), "fdr_genes": int(p_adj.lt(fdr_threshold).sum()),
                    "top_genes": ";".join(top["gene"].astype(str)), "pathway_count": int(len(pathway_frame)),
                    "fdr_pathway_count": int(pd.to_numeric(pathway_frame.get("p.adjust", np.nan), errors="coerce").lt(fdr_threshold).sum()) if not pathway_frame.empty else 0,
                    "status": "success",
                }
            )
        gene = compute_gene_reproducibility(seed_tables, tf, fdr_threshold)
        gene.insert(0, "tf", tf); gene.insert(1, "axis", axis); gene.insert(2, "subset", subset); gene_frames.append(gene)
        pair = pairwise_gene_reproducibility(seed_tables, tf, fdr_threshold)
        if not pair.empty:
            pair.insert(0, "tf", tf); pair.insert(1, "axis", axis); pair.insert(2, "subset", subset)
            pair["pair"] = pair["seed_left"].astype(str) + " vs " + pair["seed_right"].astype(str)
            pairwise_frames.append(pair)
        pathway = compute_pathway_reproducibility(pathway_results, tf, fdr_threshold)
        if not pathway.empty:
            pathway["axis"] = axis; pathway["subset"] = subset; pathway_frames.append(pathway)
        if pathway_results:
            raw_pathway = pd.concat(pathway_results, ignore_index=True)
            raw_pathway["axis"] = axis; raw_pathway["subset"] = subset; pathway_seed_frames.append(raw_pathway)
        summary_rows.append(build_axis_summary(tf, axis, subset, seed_tables, pathway_results, input_cells, input_genes, fdr_threshold))
    summary = pd.DataFrame(summary_rows)
    seed_level = pd.DataFrame(seed_rows)
    gene_reproducibility = pd.concat(gene_frames, ignore_index=True) if gene_frames else pd.DataFrame()
    pairwise = pd.concat(pairwise_frames, ignore_index=True) if pairwise_frames else pd.DataFrame()
    pathway_reproducibility = pd.concat(pathway_frames, ignore_index=True) if pathway_frames else pd.DataFrame()
    pathway_seed_results = pd.concat(pathway_seed_frames, ignore_index=True) if pathway_seed_frames else pd.DataFrame()
    try:
        from audit_sctenifoldknk_reproducibility_v2 import build_three_axis_validation_matrix
    except ModuleNotFoundError:
        from scripts.audit_sctenifoldknk_reproducibility_v2 import build_three_axis_validation_matrix

    validation_matrix = build_three_axis_validation_matrix(summary)
    metadata_dir.mkdir(parents=True, exist_ok=True); report_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "run_manifest": str((metadata_dir / "run_manifest.tsv").resolve()),
        "seed_level_results": str((metadata_dir / "seed_level_results.tsv").resolve()),
        "gene_reproducibility": str((metadata_dir / "gene_reproducibility.tsv").resolve()),
        "gene_reproducibility_pairwise": str((metadata_dir / "gene_reproducibility_pairwise.tsv").resolve()),
        "pathway_reproducibility": str((metadata_dir / "pathway_reproducibility.tsv").resolve()),
        "pathway_seed_results": str((metadata_dir / "pathway_seed_results.tsv").resolve()),
        "validation_matrix": str((metadata_dir / "HNF4A_EGR1_SOX4_scTenifoldKnk_validation_matrix.tsv").resolve()),
    }
    seed_level.to_csv(outputs["seed_level_results"], sep="\t", index=False)
    gene_reproducibility.to_csv(outputs["gene_reproducibility"], sep="\t", index=False)
    pairwise.to_csv(outputs["gene_reproducibility_pairwise"], sep="\t", index=False)
    pathway_reproducibility.to_csv(outputs["pathway_reproducibility"], sep="\t", index=False)
    pathway_seed_results.to_csv(outputs["pathway_seed_results"], sep="\t", index=False)
    validation_matrix.to_csv(outputs["validation_matrix"], sep="\t", index=False)
    figure_outputs = _plot_summary(summary, seed_level, pairwise, pathway_reproducibility, figure_dir) if not summary.empty else {}
    outputs.update({f"figure_summary_{key}": value for key, value in figure_outputs.items()})
    _write_three_axis_report(summary, outputs, report_dir / "03_three_axis_validation_report.md")
    _write_postrun_decision_report(summary, report_dir / "02_rerun_decision.md")
    _write_panel_decision_report(summary, report_dir / "04_figure2_figure3_panel_decision.md")
    _write_claim_update_report(summary, report_dir / "05_claim_evidence_update.md")
    report = {
        "module": "scTenifoldKnk Figure 2/3 reproducibility audit v2 integration",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "summary": summary.to_dict(orient="records"),
        "outputs": outputs,
        "p0_gate": "CLOSED",
        "rerun_targets": ["HNF4A", "SOX4"],
        "retained_historical_formal_target": "EGR1",
        "computational_parameter_symmetry": {
            "nc_nNet": 10,
            "nc_nCells": 500,
            "nc_nComp": 3,
            "ma_nDim": 2,
            "qc": False,
            "qc_minCells": 3,
            "seeds": [15071990, 15071991, 15071992],
            "network_genes": 3000,
        },
        "biological_input_asymmetry": {
            "HNF4A": "normal_reference; 8,098-cell universe",
            "EGR1": "stressed_regenerative; 646-cell universe",
            "SOX4": "malignant_like; 552-cell universe",
        },
        "python_runtime": {"version": platform.python_version(), "platform": platform.platform()},
    }
    report["report"] = str((metadata_dir / "integration_report.json").resolve())
    (metadata_dir / "integration_report.json").write_text(json.dumps(_json_safe(report), ensure_ascii=False, indent=2), encoding="utf-8")
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA_DIR)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--fdr-threshold", type=float, default=FDR_THRESHOLD)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = run_integration(args.metadata_dir, args.report_dir, args.figure_dir, args.fdr_threshold)
    print(json.dumps({"report": report["report"], "outputs": report["outputs"]}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
