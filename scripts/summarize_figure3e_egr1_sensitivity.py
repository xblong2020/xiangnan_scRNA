#!/usr/bin/env python3
"""Compare EGR1 scTenifoldKnk results across stress-transition sensitivity subsets."""

from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

try:
    from figure3_egr1_common import PROJECT_ROOT, TARGET_TF, json_safe, write_json
except ModuleNotFoundError:
    from scripts.figure3_egr1_common import PROJECT_ROOT, TARGET_TF, json_safe, write_json


DEFAULT_MAIN_DIR = PROJECT_ROOT / "metadata/driver/figure3e_egr1"
DEFAULT_SENSITIVITY_DIR = PROJECT_ROOT / "metadata/driver/figure3e_egr1_sensitivity"
DEFAULT_DATA_ROOT = PROJECT_ROOT / "data/processed/driver/figure3e_egr1_sctenifoldknk"
DEFAULT_GMT_DIR = PROJECT_ROOT / "metadata/driver/sctenifoldknk_module7_4_genesets"
SUBSETS = [
    "stressed_injured",
    "stressed_regenerative",
    "intermediate_pseudotime",
    "malignant_like",
]


def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    values = np.asarray(pvalues, dtype=float)
    order = np.argsort(values)
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1.0)
    output = np.empty_like(values)
    output[order] = ranked
    return output


def parse_gmt(path: Path) -> list[tuple[str, set[str]]]:
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.rstrip("\n").split("\t")
        if len(fields) >= 3:
            rows.append((fields[0], {gene.upper() for gene in fields[2:] if gene}))
    return rows


def ora(
    genes: set[str],
    background: set[str],
    gmt_dir: Path,
    min_size: int = 5,
    max_size: int = 500,
) -> pd.DataFrame:
    database_files = {
        "KEGG": "KEGG_2021_Human.gmt",
        "Reactome": "Reactome_2022.gmt",
        "GO_BP": "GO_Biological_Process_2023.gmt",
    }
    genes = {gene.upper() for gene in genes}.intersection(background)
    rows = []
    for database, filename in database_files.items():
        for term, term_genes_raw in parse_gmt(gmt_dir / filename):
            term_genes = term_genes_raw.intersection(background)
            if not min_size <= len(term_genes) <= max_size:
                continue
            overlap = genes.intersection(term_genes)
            pvalue = (
                float(
                    stats.hypergeom.sf(
                        len(overlap) - 1,
                        len(background),
                        len(term_genes),
                        len(genes),
                    )
                )
                if genes
                else 1.0
            )
            rows.append(
                {
                    "database": database,
                    "term": term,
                    "overlap_count": len(overlap),
                    "term_size": len(term_genes),
                    "input_gene_count": len(genes),
                    "background_gene_count": len(background),
                    "pvalue": pvalue,
                    "overlap_genes": ";".join(sorted(overlap)),
                }
            )
    result = pd.DataFrame(rows)
    if len(result):
        result["p.adjust"] = bh_adjust(result["pvalue"].to_numpy())
        result["significant"] = result["p.adjust"].lt(0.05)
    return result


def result_path(subset: str, main_dir: Path, sensitivity_dir: Path) -> Path:
    base = main_dir if subset == "stressed_regenerative" else sensitivity_dir
    return base / f"figure3e_egr1_{subset}_consensus_perturbation_genes.tsv"


def run(
    main_dir: Path,
    sensitivity_dir: Path,
    data_root: Path,
    gmt_dir: Path,
) -> dict:
    sensitivity_dir.mkdir(parents=True, exist_ok=True)
    audit = pd.read_csv(main_dir / "figure3e_egr1_subset_selection_audit.tsv", sep="\t")
    results: dict[str, pd.DataFrame] = {}
    significant_genes: dict[str, set[str]] = {}
    pathway_sets: dict[str, set[str]] = {}
    pathway_frames = []
    summary_rows = []
    for subset in SUBSETS:
        path = result_path(subset, main_dir, sensitivity_dir)
        if not path.exists():
            raise FileNotFoundError(path)
        table = pd.read_csv(path, sep="\t")
        table["distance"] = pd.to_numeric(table["distance"], errors="coerce")
        table["p.adj"] = pd.to_numeric(table["p.adj"], errors="coerce")
        table = table.loc[
            table["tf"].astype(str).eq(TARGET_TF)
            & ~table["gene"].astype(str).eq(TARGET_TF)
        ].copy()
        results[subset] = table
        sig = table.loc[table["p.adj"].lt(0.05) & np.isfinite(table["distance"])].copy()
        significant_genes[subset] = set(sig["gene"].astype(str))
        background_path = data_root / subset / f"figure3e_egr1_{subset}_genes.tsv"
        background = set(pd.read_csv(background_path, sep="\t").iloc[:, 0].astype(str).str.upper())
        enrichment = ora(significant_genes[subset], background, gmt_dir)
        enrichment["subset"] = subset
        pathway_frames.append(enrichment)
        sig_pathways = enrichment.loc[enrichment["significant"], ["database", "term"]]
        pathway_sets[subset] = set(
            sig_pathways["database"].astype(str) + "::" + sig_pathways["term"].astype(str)
        )
        audit_row = audit.loc[audit["subset"].eq(subset)].iloc[0]
        summary_rows.append(
            {
                "subset": subset,
                "selected_main": bool(audit_row["selected_main"]),
                "n_cells": int(audit_row["n_cells"]),
                "n_genes": int(audit_row["n_genes"]),
                "egr1_detection_rate": float(audit_row["egr1_detection_rate"]),
                "egr1_mean_expression": float(audit_row["egr1_mean_expression"]),
                "n_known_datasets": int(audit_row["n_datasets"]),
                "n_samples_or_patients": int(audit_row["n_samples_or_patients"]),
                "max_dataset_fraction": float(audit_row["max_dataset_fraction"]),
                "n_tested_excluding_egr1": int(len(table)),
                "n_significant_perturbed_genes": int(len(sig)),
                "median_distance_significant": (
                    float(sig["distance"].median()) if len(sig) else np.nan
                ),
                "max_distance_significant": (
                    float(sig["distance"].max()) if len(sig) else np.nan
                ),
                "n_fdr_significant_pathways": int(len(pathway_sets[subset])),
                "dataset_stability_evidence": (
                    "composition_eligible; network-level LODO not run"
                    if int(audit_row["n_datasets"]) >= 3 and float(audit_row["max_dataset_fraction"]) < 0.80
                    else "composition_risk"
                ),
            }
        )

    pairwise_rows = []
    for left, right in combinations(SUBSETS, 2):
        left_genes, right_genes = significant_genes[left], significant_genes[right]
        gene_union = left_genes.union(right_genes)
        gene_jaccard = len(left_genes.intersection(right_genes)) / len(gene_union) if gene_union else np.nan
        merged = results[left][["gene", "distance"]].merge(
            results[right][["gene", "distance"]],
            on="gene",
            suffixes=("_left", "_right"),
        )
        rank_rho = (
            float(stats.spearmanr(merged["distance_left"], merged["distance_right"]).statistic)
            if len(merged) > 2
            else np.nan
        )
        left_paths, right_paths = pathway_sets[left], pathway_sets[right]
        pathway_union = left_paths.union(right_paths)
        pathway_jaccard = (
            len(left_paths.intersection(right_paths)) / len(pathway_union)
            if pathway_union
            else np.nan
        )
        pairwise_rows.append(
            {
                "subset_left": left,
                "subset_right": right,
                "significant_gene_jaccard": gene_jaccard,
                "distance_rank_spearman_rho": rank_rho,
                "n_common_tested_genes": int(len(merged)),
                "fdr_pathway_jaccard": pathway_jaccard,
                "n_shared_significant_genes": int(len(left_genes.intersection(right_genes))),
                "n_shared_fdr_pathways": int(len(left_paths.intersection(right_paths))),
            }
        )

    seed_path = main_dir / "figure3e_egr1_stressed_regenerative_all_seed_perturbation_genes.tsv"
    seed_stability_rows = []
    if seed_path.exists():
        seed_data = pd.read_csv(seed_path, sep="\t")
        seed_results = {
            int(seed): frame.loc[
                frame["tf"].astype(str).eq(TARGET_TF)
                & ~frame["gene"].astype(str).eq(TARGET_TF)
            ].copy()
            for seed, frame in seed_data.groupby("seed")
        }
        for left, right in combinations(sorted(seed_results), 2):
            left_frame, right_frame = seed_results[left], seed_results[right]
            left_sig = set(left_frame.loc[pd.to_numeric(left_frame["p.adj"], errors="coerce").lt(0.05), "gene"])
            right_sig = set(right_frame.loc[pd.to_numeric(right_frame["p.adj"], errors="coerce").lt(0.05), "gene"])
            union = left_sig.union(right_sig)
            merged = left_frame[["gene", "distance"]].merge(
                right_frame[["gene", "distance"]],
                on="gene",
                suffixes=("_left", "_right"),
            )
            seed_stability_rows.append(
                {
                    "seed_left": left,
                    "seed_right": right,
                    "significant_gene_jaccard": (
                        len(left_sig.intersection(right_sig)) / len(union) if union else np.nan
                    ),
                    "distance_rank_spearman_rho": float(
                        stats.spearmanr(
                            pd.to_numeric(merged["distance_left"], errors="coerce"),
                            pd.to_numeric(merged["distance_right"], errors="coerce"),
                        ).statistic
                    ),
                    "n_shared_significant_genes": int(len(left_sig.intersection(right_sig))),
                }
            )

    summary = pd.DataFrame(summary_rows)
    pairwise = pd.DataFrame(pairwise_rows)
    pathways = pd.concat(pathway_frames, ignore_index=True)
    seed_stability = pd.DataFrame(seed_stability_rows)
    summary_path = sensitivity_dir / "figure3e_egr1_sensitivity_summary.tsv"
    pairwise_path = sensitivity_dir / "figure3e_egr1_sensitivity_pairwise.tsv"
    pathways_path = sensitivity_dir / "figure3e_egr1_sensitivity_pathway_enrichment.tsv"
    seed_path_out = sensitivity_dir / "figure3e_egr1_sensitivity_seed_stability.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    pairwise.to_csv(pairwise_path, sep="\t", index=False)
    pathways.to_csv(pathways_path, sep="\t", index=False)
    seed_stability.to_csv(seed_path_out, sep="\t", index=False)
    review_risks = [
        {
            "flag": "sensitivity_networks_lower_replication",
            "severity": "review_attention",
            "detail": "Non-main sensitivity networks used nc_nNet=3 and one fixed seed; they do not replace the formal 10-network, three-seed main analysis.",
        },
        {
            "flag": "network_level_dataset_lodo_not_run",
            "severity": "review_attention",
            "detail": "Dataset stability is documented through composition/dominance and existing CellOracle LODO evidence; dedicated scTenifoldKnk leave-one-dataset-out networks were not run.",
        },
    ]
    if len(seed_stability):
        median_seed_jaccard = float(
            pd.to_numeric(
                seed_stability["significant_gene_jaccard"], errors="coerce"
            ).median()
        )
        median_seed_rank_rho = float(
            pd.to_numeric(
                seed_stability["distance_rank_spearman_rho"], errors="coerce"
            ).median()
        )
        if np.isfinite(median_seed_jaccard) and median_seed_jaccard < 0.50:
            review_risks.append(
                {
                    "flag": "main_network_significant_gene_seed_instability",
                    "severity": "review_attention",
                    "detail": (
                        "Median pairwise Jaccard of FDR-significant genes across "
                        f"formal fixed seeds was {median_seed_jaccard:.3f}; the "
                        "maximum-p.adj consensus remains deliberately conservative."
                    ),
                }
            )
    else:
        median_seed_jaccard = np.nan
        median_seed_rank_rho = np.nan

    report = {
        "module": "Figure 3E sensitivity integration",
        "target_tf": TARGET_TF,
        "subsets": SUBSETS,
        "main_subset": "stressed_regenerative",
        "comparisons": {
            "significant_gene_count": True,
            "significant_gene_jaccard": True,
            "distance_rank_correlation": True,
            "fdr_pathway_similarity": True,
            "multiple_seed_stability": bool(len(seed_stability)),
            "dataset_composition_stability": True,
            "network_level_lodo": False,
        },
        "formal_seed_stability": {
            "median_significant_gene_jaccard": median_seed_jaccard,
            "median_distance_rank_spearman_rho": median_seed_rank_rho,
        },
        "review_risk_flags": review_risks,
        "outputs": {
            "summary": str(summary_path.resolve()),
            "pairwise": str(pairwise_path.resolve()),
            "pathway_enrichment": str(pathways_path.resolve()),
            "seed_stability": str(seed_path_out.resolve()),
        },
        "caveat": "Sensitivity comparisons use different eligible cell counts and lower replication outside the prespecified main subset.",
    }
    report_path = sensitivity_dir / "figure3e_egr1_sensitivity_report.json"
    write_json(json_safe(report), report_path)
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--main-dir", type=Path, default=DEFAULT_MAIN_DIR)
    parser.add_argument("--sensitivity-dir", type=Path, default=DEFAULT_SENSITIVITY_DIR)
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--gmt-dir", type=Path, default=DEFAULT_GMT_DIR)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = run(args.main_dir, args.sensitivity_dir, args.data_root, args.gmt_dir)
    print(json.dumps(report, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
