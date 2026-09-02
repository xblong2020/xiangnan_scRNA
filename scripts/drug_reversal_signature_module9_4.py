from __future__ import annotations

import argparse
import json
import math
import platform
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA_DIR = ROOT / "metadata/driver"
DEFAULT_DRIVER_H5AD = ROOT / "data/processed/driver/driver_hepatocyte_trajectory.module6_1.h5ad"
DEFAULT_SIGNATURE_REGISTRY = DEFAULT_METADATA_DIR / "module8_signature_registry.tsv"
DEFAULT_TF_TARGETS = DEFAULT_METADATA_DIR / "module8_tf_target_signature_genes.tsv"
DEFAULT_PATHWAY_GENES = DEFAULT_METADATA_DIR / "module8_pathway_signature_genes.tsv"
DEFAULT_TEMPORAL_CELL_SCORES = DEFAULT_METADATA_DIR / "module9_1_temporal_cell_scores.tsv.gz"

MATURE_HEPATOCYTE_MARKERS = [
    "ALB",
    "APOA1",
    "APOA2",
    "TTR",
    "HPD",
    "ASGR1",
    "CYP3A4",
    "CYP2E1",
    "CYP2C9",
    "CPS1",
    "ASS1",
]

AP1_TFS = {"JUN", "FOS", "JUND", "ATF3"}
RESCUE_TFS = {"HNF4A", "PPARA"}
CEBPB_EGR1_TFS = {"CEBPB", "EGR1"}

COMPONENT_DIRECTION_PRIORITY = {
    "mature_hepatocyte": {"up": 100},
    "hnf4a_ppara_rescue": {"up": 95},
    "c_malignant_like_fate": {"down": 90},
    "sox4_state_specific": {"down": 85},
    "ap1_stress_proliferation": {"down": 70},
    "cebpb_egr1_malignant_target": {"down": 50},
    "tier1_rescue": {"up": 60},
}

OUTPUT_STEM = "module9_4"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 9.4 drug-reversal signature preparation.")
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA_DIR)
    parser.add_argument("--driver-h5ad", type=Path, default=DEFAULT_DRIVER_H5AD)
    parser.add_argument("--signature-registry", type=Path, default=DEFAULT_SIGNATURE_REGISTRY)
    parser.add_argument("--tf-targets", type=Path, default=DEFAULT_TF_TARGETS)
    parser.add_argument("--pathway-genes", type=Path, default=DEFAULT_PATHWAY_GENES)
    parser.add_argument("--temporal-cell-scores", type=Path, default=DEFAULT_TEMPORAL_CELL_SCORES)
    parser.add_argument("--max-primary-genes-per-direction", type=int, default=150)
    parser.add_argument("--min-primary-genes-per-direction", type=int, default=20)
    parser.add_argument("--c-fate-rho-threshold", type=float, default=0.15)
    parser.add_argument("--c-fate-q-threshold", type=float, default=0.05)
    parser.add_argument("--c-fate-top-n", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260616)
    return parser.parse_args()


def package_version(name: str) -> str:
    try:
        return version(name)
    except PackageNotFoundError:
        return "not_installed"


def normalize_gene_symbol(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip().upper()


def normalize_gene_symbols(values: Sequence[object]) -> list[str]:
    return [normalize_gene_symbol(value) for value in values]


def read_tsv_or_empty(path: Path, **kwargs) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", **kwargs)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def is_housekeeping_or_qc_gene(gene: object) -> bool:
    normalized = normalize_gene_symbol(gene)
    return normalized.startswith("MT-") or normalized.startswith("RPL") or normalized.startswith("RPS") or normalized == "MALAT1"


def benjamini_hochberg(p_values: Sequence[float]) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    q = np.full(p.shape, np.nan, dtype=float)
    finite_mask = np.isfinite(p)
    finite = p[finite_mask]
    if finite.size == 0:
        return q
    order = np.argsort(finite)
    ranked = finite[order]
    n = float(finite.size)
    adjusted = ranked * n / np.arange(1, finite.size + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    out = np.empty_like(finite)
    out[order] = adjusted
    q[finite_mask] = out
    return q


def component_priority(component: str, direction: str) -> int:
    return int(COMPONENT_DIRECTION_PRIORITY.get(component, {}).get(direction, 0))


def source_priority(component: str) -> int:
    return max(COMPONENT_DIRECTION_PRIORITY.get(component, {}).values(), default=0)


def make_signature_record(
    gene: str,
    desired_direction: str,
    component: str,
    source_file: str,
    source_metric: str,
    source_rank: float,
    evidence_weight: float,
) -> dict[str, object]:
    return {
        "gene": normalize_gene_symbol(gene),
        "desired_direction": desired_direction,
        "component": component,
        "source_file": source_file,
        "source_metric": source_metric,
        "source_rank": source_rank,
        "evidence_weight": evidence_weight,
    }


def display_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def numeric_rank(values: pd.Series, default: float = 999999.0) -> pd.Series:
    out = pd.to_numeric(values, errors="coerce") if values is not None else pd.Series(dtype=float)
    return out.fillna(default)


def records_from_tf_targets(tf_targets: pd.DataFrame, source_file: Path) -> list[dict[str, object]]:
    if tf_targets.empty or not {"axis", "tf", "gene"}.issubset(tf_targets.columns):
        return []
    rows = []
    df = tf_targets.copy()
    df["axis"] = df["axis"].astype(str)
    df["tf"] = df["tf"].astype(str).str.upper()
    df["source_rank_numeric"] = numeric_rank(df["rank"] if "rank" in df.columns else pd.Series(index=df.index))
    source = display_path(source_file)
    for _, row in df.iterrows():
        axis = row["axis"]
        tf = row["tf"]
        if axis == "sox4_state_specific" and tf == "SOX4":
            rows.append(make_signature_record(row["gene"], "down", "sox4_state_specific", source, "tf_target_rank", row["source_rank_numeric"], 1.0))
        elif axis == "ap1_stress_proliferation" and tf in AP1_TFS:
            rows.append(make_signature_record(row["gene"], "down", "ap1_stress_proliferation", source, "tf_target_rank", row["source_rank_numeric"], 0.8))
        elif axis in {"ap1_stress_proliferation", "sox4_state_specific"} and tf in CEBPB_EGR1_TFS:
            rows.append(make_signature_record(row["gene"], "down", "cebpb_egr1_malignant_target", source, "tf_target_rank", row["source_rank_numeric"], 0.5))
        elif axis == "tier1_rescue" and tf in RESCUE_TFS:
            rows.append(make_signature_record(row["gene"], "up", "hnf4a_ppara_rescue", source, "tf_target_rank", row["source_rank_numeric"], 1.0))
        elif axis == "tier1_rescue":
            rows.append(make_signature_record(row["gene"], "up", "tier1_rescue", source, "tf_target_rank", row["source_rank_numeric"], 0.7))
    return rows


def records_from_pathway_genes(pathway_genes: pd.DataFrame, source_file: Path) -> list[dict[str, object]]:
    if pathway_genes.empty or not {"axis", "tf", "gene"}.issubset(pathway_genes.columns):
        return []
    rows = []
    df = pathway_genes.copy()
    df["axis"] = df["axis"].astype(str)
    df["tf"] = df["tf"].astype(str).str.upper()
    df["source_rank_numeric"] = numeric_rank(df["term_rank"] if "term_rank" in df.columns else pd.Series(index=df.index))
    source = display_path(source_file)
    for _, row in df.iterrows():
        axis = row["axis"]
        tf = row["tf"]
        if axis == "sox4_state_specific" and tf == "SOX4":
            rows.append(make_signature_record(row["gene"], "down", "sox4_state_specific", source, "pathway_term_rank", row["source_rank_numeric"], 1.0))
        elif axis == "ap1_stress_proliferation" and tf in AP1_TFS:
            rows.append(make_signature_record(row["gene"], "down", "ap1_stress_proliferation", source, "pathway_term_rank", row["source_rank_numeric"], 0.8))
        elif axis in {"ap1_stress_proliferation", "sox4_state_specific"} and tf in CEBPB_EGR1_TFS:
            rows.append(make_signature_record(row["gene"], "down", "cebpb_egr1_malignant_target", source, "pathway_term_rank", row["source_rank_numeric"], 0.5))
        elif axis == "tier1_rescue" and tf in RESCUE_TFS:
            rows.append(make_signature_record(row["gene"], "up", "hnf4a_ppara_rescue", source, "pathway_term_rank", row["source_rank_numeric"], 1.0))
        elif axis == "tier1_rescue":
            rows.append(make_signature_record(row["gene"], "up", "tier1_rescue", source, "pathway_term_rank", row["source_rank_numeric"], 0.7))
    return rows


def records_from_mature_hepatocyte_markers() -> list[dict[str, object]]:
    return [
        make_signature_record(gene, "up", "mature_hepatocyte", "fixed_mature_hepatocyte_markers", "fixed_marker", rank, 0.9)
        for rank, gene in enumerate(MATURE_HEPATOCYTE_MARKERS, start=1)
    ]


def records_from_malignant_tf_markers() -> list[dict[str, object]]:
    records = [
        make_signature_record("SOX4", "down", "sox4_state_specific", "fixed_malignant_tf_markers", "fixed_tf_marker", 1, 1.0)
    ]
    for rank, gene in enumerate(["JUN", "FOS", "JUND", "ATF3"], start=1):
        records.append(
            make_signature_record(
                gene,
                "down",
                "ap1_stress_proliferation",
                "fixed_malignant_tf_markers",
                "fixed_tf_marker",
                rank,
                0.8,
            )
        )
    return records


def compute_c_fate_correlations(
    expression: pd.DataFrame,
    fate: pd.Series,
    rho_threshold: float = 0.15,
    q_threshold: float = 0.05,
    top_n: int = 100,
) -> pd.DataFrame:
    aligned_fate = pd.to_numeric(fate.reindex(expression.index), errors="coerce")
    rows = []
    for gene in expression.columns:
        values = pd.to_numeric(expression[gene], errors="coerce")
        mask = values.notna() & aligned_fate.notna()
        if int(mask.sum()) < 4 or values.loc[mask].nunique(dropna=True) < 2 or aligned_fate.loc[mask].nunique(dropna=True) < 2:
            rows.append({"gene": gene, "spearman_rho": np.nan, "p_value": np.nan, "n_cells": int(mask.sum())})
            continue
        rho, p_value = spearmanr(values.loc[mask], aligned_fate.loc[mask])
        rows.append({"gene": gene, "spearman_rho": rho, "p_value": p_value, "n_cells": int(mask.sum())})
    result = pd.DataFrame(rows)
    if result.empty:
        return pd.DataFrame(columns=["gene", "spearman_rho", "p_value", "q_value", "n_cells", "source_rank"])
    result["q_value"] = benjamini_hochberg(result["p_value"])
    keep = (
        pd.to_numeric(result["spearman_rho"], errors="coerce").ge(rho_threshold)
        & pd.to_numeric(result["q_value"], errors="coerce").lt(q_threshold)
    )
    result = result.loc[keep].sort_values(["spearman_rho", "q_value", "gene"], ascending=[False, True, True]).head(top_n).copy()
    result["source_rank"] = np.arange(1, len(result) + 1)
    return result.reset_index(drop=True)


def read_driver_expression_for_correlation(driver_h5ad: Path) -> pd.DataFrame:
    import anndata as ad

    adata = ad.read_h5ad(driver_h5ad, backed="r")
    try:
        sub = adata.to_memory()
        x = sub.X
        if sparse.issparse(x):
            x = x.toarray()
        genes = normalize_gene_symbols(sub.var_names.astype(str))
        expr = pd.DataFrame(np.asarray(x, dtype=np.float32), index=sub.obs_names.astype(str), columns=genes)
        expr = expr.loc[:, [gene for gene in expr.columns if gene]]
        if expr.columns.has_duplicates:
            expr = expr.T.groupby(level=0).mean().T
        return expr
    finally:
        if getattr(adata, "isbacked", False):
            adata.file.close()


def read_main_strict_c_fate(temporal_cell_scores: Path) -> pd.Series:
    scores = read_tsv_or_empty(
        temporal_cell_scores,
        usecols=lambda col: col in {"cell_id", "run_id", "C_malignant_like_fate"},
    )
    if scores.empty or not {"cell_id", "C_malignant_like_fate"}.issubset(scores.columns):
        return pd.Series(dtype=float)
    if "run_id" in scores.columns:
        scores = scores.loc[scores["run_id"].astype(str).eq("main_strict")].copy()
    scores["cell_id"] = scores["cell_id"].astype(str)
    scores["C_malignant_like_fate"] = pd.to_numeric(scores["C_malignant_like_fate"], errors="coerce")
    return scores.groupby("cell_id", observed=True)["C_malignant_like_fate"].mean()


def c_fate_records_from_inputs(
    driver_h5ad: Path,
    temporal_cell_scores: Path,
    rho_threshold: float,
    q_threshold: float,
    top_n: int,
) -> tuple[list[dict[str, object]], pd.DataFrame]:
    expression = read_driver_expression_for_correlation(driver_h5ad)
    fate = read_main_strict_c_fate(temporal_cell_scores)
    common = expression.index.intersection(fate.index)
    if common.empty:
        return [], pd.DataFrame(columns=["gene", "spearman_rho", "p_value", "q_value", "n_cells", "source_rank"])
    correlations = compute_c_fate_correlations(
        expression.loc[common],
        fate.loc[common],
        rho_threshold=rho_threshold,
        q_threshold=q_threshold,
        top_n=top_n,
    )
    records = [
        make_signature_record(
            row["gene"],
            "down",
            "c_malignant_like_fate",
            display_path(temporal_cell_scores),
            "spearman_rho",
            row["source_rank"],
            1.0,
        )
        for _, row in correlations.iterrows()
    ]
    return records, correlations


def collapse_values(values: pd.Series) -> str:
    return ";".join([str(value) for value in dict.fromkeys(values.dropna().astype(str))])


def choose_direction(gene: str, group: pd.DataFrame) -> tuple[str, bool]:
    directions = sorted(group["desired_direction"].dropna().astype(str).unique())
    if len(directions) == 1:
        return directions[0], False
    if gene in AP1_TFS.union({"SOX4"}) and group["desired_direction"].astype(str).eq("down").any():
        return "down", False
    if gene in RESCUE_TFS.union(MATURE_HEPATOCYTE_MARKERS) and group["desired_direction"].astype(str).eq("up").any():
        return "up", False
    direction_scores = {}
    for direction, sub in group.groupby("desired_direction", observed=True):
        direction_scores[direction] = int(max(component_priority(component, direction) for component in sub["component"].astype(str)))
    best_direction = max(direction_scores, key=direction_scores.get)
    best_score = direction_scores[best_direction]
    second_score = max(score for direction, score in direction_scores.items() if direction != best_direction)
    if best_score >= 80 and best_score > second_score:
        return best_direction, False
    return best_direction, True


def resolve_signature_records(records: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "gene",
        "desired_direction",
        "component",
        "source_file",
        "source_metric",
        "source_rank",
        "evidence_weight",
        "conflict_flag",
        "housekeeping_or_qc_flag",
        "final_weight",
        "include_primary",
        "include_sensitivity",
    ]
    if records.empty:
        return pd.DataFrame(columns=columns)
    required = {"gene", "desired_direction", "component", "source_file", "source_metric", "source_rank", "evidence_weight"}
    missing = required.difference(records.columns)
    if missing:
        raise ValueError(f"signature records missing columns: {sorted(missing)}")
    work = records.copy()
    work["gene"] = work["gene"].map(normalize_gene_symbol)
    work = work.loc[work["gene"].ne("")].copy()
    work["source_rank"] = pd.to_numeric(work["source_rank"], errors="coerce").fillna(999999.0)
    work["evidence_weight"] = pd.to_numeric(work["evidence_weight"], errors="coerce").fillna(0.0)

    resolved_rows = []
    for gene, group in work.groupby("gene", sort=True, observed=True):
        direction, conflict = choose_direction(gene, group)
        selected = group.loc[group["desired_direction"].astype(str).eq(direction)].copy()
        best_priority = selected["component"].map(source_priority).max()
        best = selected.loc[selected["component"].map(source_priority).eq(best_priority)].copy()
        max_weight = float(selected["evidence_weight"].max())
        support_bonus = min(0.2, 0.02 * max(0, len(selected) - 1))
        housekeeping = is_housekeeping_or_qc_gene(gene)
        resolved_rows.append(
            {
                "gene": gene,
                "desired_direction": direction,
                "component": collapse_values(best.sort_values(["source_rank", "component"])["component"]),
                "source_file": collapse_values(selected["source_file"]),
                "source_metric": collapse_values(selected["source_metric"]),
                "source_rank": float(selected["source_rank"].min()) if len(selected) else math.nan,
                "evidence_weight": max_weight,
                "conflict_flag": bool(conflict),
                "housekeeping_or_qc_flag": bool(housekeeping),
                "final_weight": max_weight + support_bonus,
                "include_primary": bool(not conflict and not housekeeping),
                "include_sensitivity": True,
            }
        )
    return pd.DataFrame(resolved_rows, columns=columns).sort_values(["desired_direction", "final_weight", "source_rank", "gene"], ascending=[True, False, True, True]).reset_index(drop=True)


def rank_eligible_direction(frame: pd.DataFrame, direction: str) -> pd.DataFrame:
    sub = frame.loc[
        frame["desired_direction"].eq(direction)
        & frame["include_primary"].astype(bool)
        & ~frame["conflict_flag"].astype(bool)
        & ~frame["housekeeping_or_qc_flag"].astype(bool)
    ].copy()
    if sub.empty:
        return sub
    sub["_component_priority"] = sub["component"].astype(str).map(lambda value: max(source_priority(part) for part in value.split(";") if part))
    sub["_fixed_anchor"] = sub["source_metric"].astype(str).str.contains("fixed_tf_marker|fixed_marker", regex=True)
    return sub.sort_values(
        ["_fixed_anchor", "final_weight", "_component_priority", "source_rank", "gene"],
        ascending=[False, False, False, True, True],
    )


def build_primary_signature(
    resolved: pd.DataFrame,
    max_genes_per_direction: int = 150,
    min_genes_per_direction: int = 20,
) -> dict[str, list[str]]:
    primary = {}
    for direction in ["up", "down"]:
        ranked = rank_eligible_direction(resolved, direction)
        selected = ranked.head(max_genes_per_direction)["gene"].tolist()
        primary[direction] = selected
    return primary


def apply_primary_flags(
    resolved: pd.DataFrame,
    max_genes_per_direction: int,
    min_genes_per_direction: int,
) -> tuple[pd.DataFrame, dict[str, list[str]], pd.DataFrame]:
    out = resolved.copy()
    out["include_primary"] = False
    primary = build_primary_signature(out.assign(include_primary=resolved["include_primary"]), max_genes_per_direction, min_genes_per_direction)
    qc_rows = []
    for direction, genes in primary.items():
        out.loc[out["gene"].isin(genes) & out["desired_direction"].eq(direction), "include_primary"] = True
        eligible_n = int(len(rank_eligible_direction(resolved, direction)))
        qc_rows.append(
            {
                "direction": direction,
                "n_primary": int(len(genes)),
                "n_eligible": eligible_n,
                "min_primary_genes_per_direction": int(min_genes_per_direction),
                "max_primary_genes_per_direction": int(max_genes_per_direction),
                "status": "pass" if len(genes) >= min_genes_per_direction else "below_min_primary_gene_count",
            }
        )
    return out, primary, pd.DataFrame(qc_rows)


def write_gmt(path: Path, name: str, description: str, genes: Sequence[str]) -> None:
    path.write_text("\t".join([name, description, *genes]) + "\n", encoding="utf-8")


def write_candidate_template(path: Path) -> None:
    columns = [
        "compound_id",
        "compound_name",
        "source_database",
        "cell_line",
        "dose",
        "time",
        "disease_reversal_score",
        "rescue_mimicry_score",
        "sox4_ap1_suppression_score",
        "hnf4a_ppara_support_score",
        "hepatocyte_context_score",
        "replicate_robustness_score",
        "toxicity_or_pan_stress_penalty",
        "final_rank_score",
        "rank_tier",
        "evidence_notes",
    ]
    pd.DataFrame(columns=columns).to_csv(path, sep="\t", index=False)


def build_signature_records_from_inputs(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    tf_targets = read_tsv_or_empty(args.tf_targets)
    pathway_genes = read_tsv_or_empty(args.pathway_genes)
    records = []
    records.extend(records_from_tf_targets(tf_targets, args.tf_targets))
    records.extend(records_from_pathway_genes(pathway_genes, args.pathway_genes))
    records.extend(records_from_mature_hepatocyte_markers())
    records.extend(records_from_malignant_tf_markers())
    c_fate_records, c_fate_correlations = c_fate_records_from_inputs(
        args.driver_h5ad,
        args.temporal_cell_scores,
        rho_threshold=args.c_fate_rho_threshold,
        q_threshold=args.c_fate_q_threshold,
        top_n=args.c_fate_top_n,
    )
    records.extend(c_fate_records)
    return pd.DataFrame(records), c_fate_correlations


def build_qc_table(resolved: pd.DataFrame, primary_qc: pd.DataFrame, c_fate_correlations: pd.DataFrame) -> pd.DataFrame:
    rows = []
    rows.extend(primary_qc.to_dict(orient="records"))
    rows.extend(
        [
            {"metric": "n_total_signature_genes", "value": int(resolved["gene"].nunique()), "status": "reported"},
            {"metric": "n_conflict_genes", "value": int(resolved["conflict_flag"].sum()), "status": "reported"},
            {"metric": "n_housekeeping_or_qc_genes", "value": int(resolved["housekeeping_or_qc_flag"].sum()), "status": "reported"},
            {"metric": "n_c_fate_positive_genes", "value": int(len(c_fate_correlations)), "status": "reported"},
        ]
    )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    np.random.seed(args.seed)
    metadata_dir = args.metadata_dir
    metadata_dir.mkdir(parents=True, exist_ok=True)

    raw_records, c_fate_correlations = build_signature_records_from_inputs(args)
    resolved = resolve_signature_records(raw_records)
    resolved, primary, primary_qc = apply_primary_flags(
        resolved,
        max_genes_per_direction=args.max_primary_genes_per_direction,
        min_genes_per_direction=args.min_primary_genes_per_direction,
    )
    qc = build_qc_table(resolved, primary_qc, c_fate_correlations)

    outputs = {
        "signature_tsv": metadata_dir / f"{OUTPUT_STEM}_drug_reversal_signature.tsv",
        "up_gmt": metadata_dir / f"{OUTPUT_STEM}_drug_reversal_up.gmt",
        "down_gmt": metadata_dir / f"{OUTPUT_STEM}_drug_reversal_down.gmt",
        "signature_json": metadata_dir / f"{OUTPUT_STEM}_drug_reversal_signature.json",
        "signature_qc": metadata_dir / f"{OUTPUT_STEM}_signature_qc.tsv",
        "candidate_template": metadata_dir / f"{OUTPUT_STEM}_candidate_drug_ranking_template.tsv",
        "report": metadata_dir / f"{OUTPUT_STEM}_report.json",
        "c_fate_correlations": metadata_dir / f"{OUTPUT_STEM}_c_malignant_like_fate_correlations.tsv",
    }

    resolved.to_csv(outputs["signature_tsv"], sep="\t", index=False)
    c_fate_correlations.to_csv(outputs["c_fate_correlations"], sep="\t", index=False)
    qc.to_csv(outputs["signature_qc"], sep="\t", index=False)
    write_gmt(outputs["up_gmt"], "module9_4_drug_reversal_up", "desired rescue signature", primary["up"])
    write_gmt(outputs["down_gmt"], "module9_4_drug_reversal_down", "desired disease reversal signature", primary["down"])
    write_candidate_template(outputs["candidate_template"])

    signature_payload = {
        "module": "module9_4_drug_reversal_signature",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "primary": primary,
        "sensitivity": {
            "up": resolved.loc[resolved["desired_direction"].eq("up"), "gene"].tolist(),
            "down": resolved.loc[resolved["desired_direction"].eq("down"), "gene"].tolist(),
        },
        "parameters": {
            "max_primary_genes_per_direction": args.max_primary_genes_per_direction,
            "min_primary_genes_per_direction": args.min_primary_genes_per_direction,
            "c_fate_rho_threshold": args.c_fate_rho_threshold,
            "c_fate_q_threshold": args.c_fate_q_threshold,
            "c_fate_top_n": args.c_fate_top_n,
            "seed": args.seed,
        },
    }
    outputs["signature_json"].write_text(json.dumps(signature_payload, indent=2, sort_keys=True), encoding="utf-8")

    report = {
        "module": "module9_4_drug_reversal_signature",
        "status": "completed",
        "inputs": {
            "metadata_dir": str(metadata_dir.resolve()),
            "driver_h5ad": str(args.driver_h5ad.resolve()),
            "signature_registry": str(args.signature_registry.resolve()),
            "tf_targets": str(args.tf_targets.resolve()),
            "pathway_genes": str(args.pathway_genes.resolve()),
            "temporal_cell_scores": str(args.temporal_cell_scores.resolve()),
        },
        "outputs": {key: str(path.resolve()) for key, path in outputs.items()},
        "summary": {
            "n_raw_records": int(len(raw_records)),
            "n_resolved_genes": int(resolved["gene"].nunique()),
            "n_primary_up": int(len(primary["up"])),
            "n_primary_down": int(len(primary["down"])),
            "n_conflict_genes": int(resolved["conflict_flag"].sum()),
            "n_housekeeping_or_qc_genes": int(resolved["housekeeping_or_qc_flag"].sum()),
            "n_c_fate_positive_genes": int(len(c_fate_correlations)),
        },
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "pandas": package_version("pandas"),
            "numpy": package_version("numpy"),
            "scipy": package_version("scipy"),
            "anndata": package_version("anndata"),
        },
    }
    outputs["report"].write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(report["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
