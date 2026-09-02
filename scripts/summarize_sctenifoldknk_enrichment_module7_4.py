from __future__ import annotations

import argparse
import json
import math
import platform
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote

import numpy as np
import pandas as pd
import requests
from scipy.stats import hypergeom


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA_DIR = PROJECT_ROOT / "metadata/driver"
DEFAULT_ENRICHMENT = DEFAULT_METADATA_DIR / "sctenifoldknk_module7_4_enrichment_all.tsv"
DEFAULT_PERTURBATION = DEFAULT_METADATA_DIR / "sctenifoldknk_module7_2_driver_union_all_perturbation_genes.tsv"
DEFAULT_BACKGROUND = PROJECT_ROOT / "data/processed/driver/sctenifoldknk_module7_1/driver_union_all/sctenifoldknk_genes.tsv"
DEFAULT_GENESET_DIR = DEFAULT_METADATA_DIR / "sctenifoldknk_module7_4_genesets"
ENRICHR_LIBRARIES = {
    "KEGG_ORA": "KEGG_2021_Human",
    "Reactome_ORA": "Reactome_2022",
    "GO_BP_ORA": "GO_Biological_Process_2023",
    "KEGG_GSEA": "KEGG_2021_Human",
    "Reactome_GSEA": "Reactome_2022",
    "GO_BP_GSEA": "GO_Biological_Process_2023",
}


def _fdr_col(df: pd.DataFrame) -> str:
    for col in ["p.adj", "p_adj", "padj", "FDR", "fdr", "qvalue", "q_value"]:
        if col in df.columns:
            return col
    raise ValueError("No FDR-adjusted p-value column found")


def _distance_col(df: pd.DataFrame) -> str:
    for col in ["distance", "Distance", "dist", "perturbation_score", "score"]:
        if col in df.columns:
            return col
    raise ValueError("No distance/perturbation score column found")


def compute_preranked_metric(perturbation_genes: pd.DataFrame) -> pd.DataFrame:
    fdr_col = _fdr_col(perturbation_genes)
    distance_col = _distance_col(perturbation_genes)
    out = perturbation_genes.copy()
    fdr = pd.to_numeric(out[fdr_col], errors="coerce").fillna(1.0).clip(lower=np.nextafter(0, 1), upper=1.0)
    distance = pd.to_numeric(out[distance_col], errors="coerce").fillna(0.0)
    out["preranked_metric"] = -np.log10(fdr) * np.sign(distance)
    out.loc[np.isclose(fdr, 1.0), "preranked_metric"] = 0.0
    sort_cols = ["preranked_metric", "gene"]
    ascending = [False, True]
    if "tf" in out.columns:
        sort_cols = ["tf"] + sort_cols
        ascending = [True] + ascending
    return out.sort_values(sort_cols, ascending=ascending).reset_index(drop=True)


def parse_gmt_lines(lines: list[str]) -> dict[str, set[str]]:
    gene_sets: dict[str, set[str]] = {}
    for line in lines:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        term = parts[0].strip()
        genes = {gene.strip() for gene in parts[2:] if gene.strip()}
        if term and genes:
            gene_sets[term] = genes
    return gene_sets


def filter_gene_sets_to_background(
    gene_sets: dict[str, set[str]],
    background: set[str],
    min_size: int = 5,
    max_size: int = 500,
) -> dict[str, set[str]]:
    out = {}
    for term, genes in gene_sets.items():
        overlap = set(genes) & set(background)
        if min_size <= len(overlap) <= max_size:
            out[term] = overlap
    return out


def _bh_adjust(pvalues: list[float]) -> list[float]:
    if not pvalues:
        return []
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    n = len(p)
    adjusted = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        value = min(prev, ranked[i] * n / (i + 1))
        adjusted[i] = value
        prev = value
    out = np.empty(n, dtype=float)
    out[order] = np.clip(adjusted, 0, 1)
    return out.tolist()


def run_ora_for_tf(
    tf: str,
    subset: str,
    significant_genes: set[str],
    background: set[str],
    gene_sets: dict[str, set[str]],
    database: str,
) -> pd.DataFrame:
    sig = set(significant_genes) & set(background)
    universe_size = len(background)
    input_size = len(sig)
    rows = []
    if universe_size == 0 or input_size == 0:
        return pd.DataFrame()
    for term, genes in gene_sets.items():
        term_genes = set(genes) & set(background)
        overlap = sig & term_genes
        if not overlap:
            continue
        pvalue = float(hypergeom.sf(len(overlap) - 1, universe_size, len(term_genes), input_size))
        rows.append(
            {
                "analysis": "ORA",
                "subset": subset,
                "tf": tf,
                "database": database,
                "term_id": term,
                "term_name": term,
                "pvalue": pvalue,
                "gene_count": int(len(overlap)),
                "term_size": int(len(term_genes)),
                "input_gene_count": int(input_size),
                "background_gene_count": int(universe_size),
                "overlap_genes": ";".join(sorted(overlap)),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["p.adjust"] = _bh_adjust(out["pvalue"].tolist())
    return out.sort_values(["p.adjust", "pvalue", "gene_count"], ascending=[True, True, False]).reset_index(drop=True)


def _enrichment_score(ranked_metric: pd.Series, genes: set[str]) -> float:
    ranked = ranked_metric.dropna().sort_values(ascending=False)
    hits = ranked.index.astype(str).isin(genes)
    n_hits = int(hits.sum())
    n_miss = int((~hits).sum())
    if n_hits == 0 or n_miss == 0:
        return 0.0
    weights = ranked.abs().to_numpy(dtype=float)
    hit_weights = np.where(hits, weights, 0.0)
    hit_norm = hit_weights.sum()
    if hit_norm == 0:
        hit_weights = np.where(hits, 1.0, 0.0)
        hit_norm = hit_weights.sum()
    running = np.cumsum(np.where(hits, hit_weights / hit_norm, -1.0 / n_miss))
    max_es = float(running.max())
    min_es = float(running.min())
    return max_es if abs(max_es) >= abs(min_es) else min_es


def compute_simple_gsea(
    tf: str,
    subset: str,
    ranked_metric: pd.Series,
    gene_sets: dict[str, set[str]],
    database: str,
    permutations: int = 0,
    seed: int = 1,
) -> pd.DataFrame:
    ranked = ranked_metric.dropna()
    ranked.index = ranked.index.astype(str)
    rows = []
    rng = np.random.default_rng(seed)
    for term, genes in gene_sets.items():
        overlap = set(genes) & set(ranked.index)
        if not overlap:
            continue
        es = _enrichment_score(ranked, overlap)
        nes = es
        pvalue = np.nan
        if permutations > 0:
            null = []
            values = ranked.to_numpy(dtype=float)
            index = ranked.index.to_numpy(dtype=str)
            for _ in range(permutations):
                shuffled = pd.Series(values, index=rng.permutation(index))
                null.append(_enrichment_score(shuffled, overlap))
            null = np.asarray(null, dtype=float)
            same_tail = np.abs(null) >= abs(es)
            pvalue = float((same_tail.sum() + 1) / (len(null) + 1))
            denom = np.mean(np.abs(null[np.sign(null) == np.sign(es)]))
            if np.isfinite(denom) and denom > 0:
                nes = float(es / denom)
        rows.append(
            {
                "analysis": "GSEA",
                "subset": subset,
                "tf": tf,
                "database": database,
                "term_id": term,
                "term_name": term,
                "pvalue": pvalue,
                "p.adjust": np.nan,
                "gene_count": int(len(overlap)),
                "term_size": int(len(genes)),
                "input_gene_count": int(len(ranked)),
                "background_gene_count": int(len(ranked)),
                "ES": float(es),
                "NES": float(nes),
                "overlap_genes": ";".join(sorted(overlap)),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    if out["pvalue"].notna().any():
        out["p.adjust"] = _bh_adjust(out["pvalue"].fillna(1.0).tolist())
        return out.sort_values(["p.adjust", "NES"], ascending=[True, False]).reset_index(drop=True)
    return out.sort_values(["NES", "gene_count"], ascending=[False, False]).reset_index(drop=True)


def download_enrichr_library(library_name: str, out_path: Path, retries: int = 4) -> Path:
    if out_path.exists() and out_path.stat().st_size > 0:
        return out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)
    url = f"https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={quote(library_name)}"
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            response = requests.get(url, timeout=120)
            response.raise_for_status()
            out_path.write_text(response.text, encoding="utf-8")
            return out_path
        except Exception as exc:  # pragma: no cover - network-dependent
            last_error = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"Failed to download Enrichr library {library_name}: {last_error}")


def summarize_mapping_stats(
    tf: str,
    database: str,
    n_background: int,
    n_input: int,
    n_mapped_background: int,
    n_mapped_input: int,
) -> dict:
    return {
        "tf": tf,
        "database": database,
        "n_background": int(n_background),
        "n_input": int(n_input),
        "n_mapped_background": int(n_mapped_background),
        "n_mapped_input": int(n_mapped_input),
        "background_mapping_rate": float(n_mapped_background / n_background) if n_background else math.nan,
        "input_mapping_rate": float(n_mapped_input / n_input) if n_input else math.nan,
    }


def build_enrichment_summary(enrichment: pd.DataFrame, top_n: int = 5) -> pd.DataFrame:
    if enrichment.empty:
        return pd.DataFrame()
    required = {"tf", "database", "term_id", "term_name"}
    missing = required - set(enrichment.columns)
    if missing:
        raise ValueError(f"Enrichment table missing required columns: {sorted(missing)}")
    work = enrichment.copy()
    if "p.adjust" not in work.columns:
        if "p_adj" in work.columns:
            work["p.adjust"] = work["p_adj"]
        else:
            raise ValueError("Enrichment table must contain p.adjust or p_adj")
    if "gene_count" not in work.columns:
        work["gene_count"] = np.nan
    work["p.adjust"] = pd.to_numeric(work["p.adjust"], errors="coerce").fillna(1.0)
    work["gene_count"] = pd.to_numeric(work["gene_count"], errors="coerce").fillna(0)
    if "analysis" not in work.columns:
        work["analysis"] = "ORA"
    if "NES" not in work.columns:
        work["NES"] = np.nan
    if "subset" not in work.columns:
        work["subset"] = "driver_union_all"
    work["_sort_p"] = work["p.adjust"].fillna(1.0)
    work["_sort_nes"] = pd.to_numeric(work["NES"], errors="coerce").abs().fillna(0.0)
    work = work.sort_values(
        ["subset", "tf", "analysis", "database", "_sort_p", "_sort_nes", "gene_count", "term_name"],
        ascending=[True, True, True, True, True, False, False, True],
    )
    return (
        work.groupby(["subset", "tf", "analysis", "database"], as_index=False, group_keys=False)
        .head(top_n)
        .drop(columns=["_sort_p", "_sort_nes"])
        .reset_index(drop=True)
    )


def run_enrichment(
    perturbation_path: Path,
    background_path: Path,
    metadata_dir: Path,
    geneset_dir: Path,
    fdr_threshold: float,
    min_size: int,
    max_size: int,
    gsea_permutations: int,
) -> dict:
    perturb = pd.read_csv(perturbation_path, sep="\t")
    background = set(pd.read_csv(background_path, sep="\t").iloc[:, 0].astype(str))
    ranked = compute_preranked_metric(perturb)
    fdr_col = _fdr_col(perturb)
    gene_sets_by_db = {}
    source_files = {}
    for database, library_name in ENRICHR_LIBRARIES.items():
        if database.endswith("_GSEA"):
            continue
        gmt_path = download_enrichr_library(library_name, geneset_dir / f"{library_name}.gmt")
        source_files[database] = str(gmt_path)
        parsed = parse_gmt_lines(gmt_path.read_text(encoding="utf-8").splitlines())
        gene_sets_by_db[database] = filter_gene_sets_to_background(parsed, background, min_size=min_size, max_size=max_size)
    gene_sets_by_db["KEGG_GSEA"] = gene_sets_by_db["KEGG_ORA"]
    gene_sets_by_db["Reactome_GSEA"] = gene_sets_by_db["Reactome_ORA"]
    gene_sets_by_db["GO_BP_GSEA"] = gene_sets_by_db["GO_BP_ORA"]

    enrichment_frames = []
    mapping_rows = []
    subset = str(perturb["subset"].iloc[0]) if "subset" in perturb.columns and len(perturb) else "driver_union_all"
    for tf, tf_df in perturb.groupby("tf", sort=False):
        tf_df = tf_df.copy()
        tf_df[fdr_col] = pd.to_numeric(tf_df[fdr_col], errors="coerce").fillna(1.0)
        significant = set(tf_df.loc[tf_df[fdr_col] <= fdr_threshold, "gene"].astype(str)) & background
        for database in ["KEGG_ORA", "Reactome_ORA", "GO_BP_ORA"]:
            gene_sets = gene_sets_by_db[database]
            enrichment_frames.append(run_ora_for_tf(str(tf), subset, significant, background, gene_sets, database))
            mapped_input = set().union(*gene_sets.values()).intersection(significant) if gene_sets else set()
            mapping_rows.append(
                summarize_mapping_stats(
                    tf=str(tf),
                    database=database,
                    n_background=len(background),
                    n_input=len(significant),
                    n_mapped_background=len(set().union(*gene_sets.values())) if gene_sets else 0,
                    n_mapped_input=len(mapped_input),
                )
            )
        tf_ranked = ranked.loc[ranked["tf"].astype(str).eq(str(tf))]
        metric = tf_ranked.drop_duplicates("gene").set_index("gene")["preranked_metric"]
        for database in ["KEGG_GSEA", "Reactome_GSEA", "GO_BP_GSEA"]:
            enrichment_frames.append(
                compute_simple_gsea(
                    str(tf),
                    subset,
                    metric,
                    gene_sets_by_db[database],
                    database,
                    permutations=gsea_permutations,
                    seed=1,
                )
            )
    enrichment = pd.concat([df for df in enrichment_frames if not df.empty], ignore_index=True) if enrichment_frames else pd.DataFrame()
    mapping = pd.DataFrame(mapping_rows)
    summary = build_enrichment_summary(enrichment, top_n=5) if not enrichment.empty else pd.DataFrame()

    outputs = {
        "enrichment_all": str(metadata_dir / "sctenifoldknk_module7_4_enrichment_all.tsv"),
        "top_enrichment_summary": str(metadata_dir / "sctenifoldknk_module7_4_top_enrichment_summary.tsv"),
        "mapping_stats": str(metadata_dir / "sctenifoldknk_module7_4_mapping_stats.tsv"),
        "preranked_metrics": str(metadata_dir / "sctenifoldknk_module7_4_preranked_metrics.tsv.gz"),
        "report": str(metadata_dir / "sctenifoldknk_module7_4_report.json"),
    }
    metadata_dir.mkdir(parents=True, exist_ok=True)
    enrichment.to_csv(outputs["enrichment_all"], sep="\t", index=False)
    summary.to_csv(outputs["top_enrichment_summary"], sep="\t", index=False)
    mapping.to_csv(outputs["mapping_stats"], sep="\t", index=False)
    ranked.to_csv(outputs["preranked_metrics"], sep="\t", index=False, compression="gzip")
    report = {
        "module": "7.4",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "method": "ORA and preranked GSEA using Enrichr gene-set libraries with CellOracle/scTenifoldKnk background",
        "inputs": {"perturbation": str(perturbation_path), "background_genes": str(background_path), "geneset_sources": source_files},
        "outputs": outputs,
        "fdr_threshold": fdr_threshold,
        "min_size": min_size,
        "max_size": max_size,
        "gsea_permutations": gsea_permutations,
        "n_background_genes": int(len(background)),
        "n_enrichment_rows": int(len(enrichment)),
        "n_summary_rows": int(len(summary)),
        "databases": sorted(enrichment["database"].dropna().unique().tolist()) if not enrichment.empty else [],
        "python_runtime": {"version": platform.python_version(), "platform": platform.platform()},
    }
    Path(outputs["report"]).write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def run_summary(enrichment_path: Path, metadata_dir: Path, top_n: int) -> dict:
    enrichment = pd.read_csv(enrichment_path, sep="\t")
    summary = build_enrichment_summary(enrichment, top_n=top_n)
    outputs = {
        "top_enrichment_summary": str(metadata_dir / "sctenifoldknk_module7_4_top_enrichment_summary.tsv"),
        "report": str(metadata_dir / "sctenifoldknk_module7_4_report.json"),
    }
    metadata_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(outputs["top_enrichment_summary"], sep="\t", index=False)
    report = {
        "module": "7.4",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {"enrichment": str(enrichment_path)},
        "outputs": outputs,
        "top_n": int(top_n),
        "n_enrichment_rows": int(len(enrichment)),
        "n_summary_rows": int(len(summary)),
        "python_runtime": {"version": platform.python_version(), "platform": platform.platform()},
    }
    Path(outputs["report"]).write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 7.4 summarize scTenifoldKnk KEGG/Reactome/GSEA enrichment")
    parser.add_argument("--perturbation", type=Path, default=DEFAULT_PERTURBATION)
    parser.add_argument("--background-genes", type=Path, default=DEFAULT_BACKGROUND)
    parser.add_argument("--geneset-dir", type=Path, default=DEFAULT_GENESET_DIR)
    parser.add_argument("--enrichment", type=Path, default=DEFAULT_ENRICHMENT)
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA_DIR)
    parser.add_argument("--top-n", type=int, default=5)
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    parser.add_argument("--min-size", type=int, default=5)
    parser.add_argument("--max-size", type=int, default=500)
    parser.add_argument("--gsea-permutations", type=int, default=0)
    parser.add_argument("--summary-only", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.summary_only:
        report = run_summary(args.enrichment, args.metadata_dir, args.top_n)
    else:
        report = run_enrichment(
            args.perturbation,
            args.background_genes,
            args.metadata_dir,
            args.geneset_dir,
            args.fdr_threshold,
            args.min_size,
            args.max_size,
            args.gsea_permutations,
        )
    print(json.dumps({"report": report["outputs"]["report"], "n_summary_rows": report["n_summary_rows"]}, indent=2))


if __name__ == "__main__":
    main()
