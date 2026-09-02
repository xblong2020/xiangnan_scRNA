from __future__ import annotations

"""Export raw CellOracle KO delta-expression summaries for Figure 6.

This script is the only Figure 6 Python computation. It does not integrate,
test or plot results. R performs all aggregation, bootstrap inference and
formal plotting.
"""

import argparse
import gzip
import json
import os
import platform
import sys
import time
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ORACLE = ROOT / "data/processed/driver/celloracle_module6_7/celloracle_module6_7_fitted.celloracle.oracle"
DEFAULT_METADATA = ROOT / "metadata/driver/figure6_directional_network"
DEFAULT_AVAILABILITY = ROOT / "metadata/trajectory/trajectory_module5_4_module_gene_availability.tsv"
DEFAULT_TF_TARGETS = ROOT / "metadata/driver/module8_tf_target_signature_genes.tsv"
DEFAULT_SOX4_STATE = ROOT / "metadata/driver/sctenifoldknk_module7_3_malignant_like_state_specific_genes.tsv"
DEFAULT_TFS = [
    "HNF4A", "PPARA", "EGR1", "CEBPB", "JUN", "JUNB", "JUND", "FOS", "ATF3",
    "SOX4", "HLF", "IRF1", "MAFB", "MAFF", "MYC",
]


def normalize_gene(value: object) -> str:
    return "" if pd.isna(value) else str(value).strip().upper().split(".")[0]


def unique_ordered(values: Iterable[object]) -> list[str]:
    return [x for x in dict.fromkeys(normalize_gene(v) for v in values) if x]


def split_genes(value: object) -> list[str]:
    if pd.isna(value):
        return []
    return unique_ordered(str(value).split(";"))


def select_sox4_targets(tf_targets: pd.DataFrame, state_specific: pd.DataFrame, top_n: int = 50) -> list[str]:
    genes: list[str] = []
    if {"tf", "gene"}.issubset(tf_targets.columns):
        sub = tf_targets.loc[tf_targets["tf"].astype(str).str.upper().eq("SOX4")].copy()
        if "rank" in sub.columns:
            sub["rank"] = pd.to_numeric(sub["rank"], errors="coerce")
            sub = sub.sort_values("rank")
        genes.extend(sub["gene"].head(top_n))
    if {"tf", "gene"}.issubset(state_specific.columns):
        sub = state_specific.loc[state_specific["tf"].astype(str).str.upper().eq("SOX4")].copy()
        sort_cols = [c for c in ["malignant_like_specificity_ratio", "malignant_like_fdr"] if c in sub.columns]
        if sort_cols:
            for col in sort_cols:
                sub[col] = pd.to_numeric(sub[col], errors="coerce")
            sub = sub.sort_values(sort_cols, ascending=[c != "malignant_like_specificity_ratio" for c in sort_cols])
        genes.extend(sub["gene"].head(top_n))
    return unique_ordered(genes)


def build_frozen_programmes(
    availability: pd.DataFrame,
    tf_targets: pd.DataFrame,
    state_specific: pd.DataFrame,
) -> dict[str, list[str]]:
    main = availability.loc[availability["run_id"].astype(str).eq("main_strict")].copy()
    by_module = {
        str(row.module): split_genes(row.genes_available)
        for row in main.itertuples(index=False)
    }
    sox4_targets = select_sox4_targets(tf_targets, state_specific, top_n=50)
    return {
        "identity_program_change": by_module.get("Mature_Hepatocyte", []),
        "stress_transition_change": by_module.get("Stressed_Injured", []),
        "sox4_programme_change": unique_ordered(["SOX4", *sox4_targets]),
        "proliferation_change": by_module.get("Proliferation", []),
        "cnv_malignant_signature_change": by_module.get("HCC_Malignant_Associated", []),
    }


def score_programmes(
    matrix: np.ndarray,
    gene_names: list[str],
    programmes: dict[str, list[str]],
) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    gene_to_ix = {normalize_gene(gene): ix for ix, gene in enumerate(gene_names)}
    out: dict[str, np.ndarray] = {}
    manifest: list[dict[str, object]] = []
    for programme, requested in programmes.items():
        available = [gene for gene in requested if gene in gene_to_ix]
        missing = [gene for gene in requested if gene not in gene_to_ix]
        out[programme] = (
            np.asarray(matrix[:, [gene_to_ix[g] for g in available]].mean(axis=1)).ravel()
            if available
            else np.full(matrix.shape[0], np.nan)
        )
        manifest.append(
            {
                "programme": programme,
                "n_requested": len(requested),
                "n_available": len(available),
                "genes_requested": ";".join(requested),
                "genes_available": ";".join(available),
                "genes_missing": ";".join(missing),
                "score_definition": "unweighted mean across frozen available genes",
            }
        )
    return pd.DataFrame(out), manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Figure 6 raw CellOracle delta-expression exporter")
    parser.add_argument("--oracle", type=Path, default=DEFAULT_ORACLE)
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--availability", type=Path, default=DEFAULT_AVAILABILITY)
    parser.add_argument("--tf-targets", type=Path, default=DEFAULT_TF_TARGETS)
    parser.add_argument("--sox4-state-specific", type=Path, default=DEFAULT_SOX4_STATE)
    parser.add_argument("--tfs", nargs="*", default=DEFAULT_TFS)
    parser.add_argument("--n-propagation", type=int, default=3)
    parser.add_argument("--seed", type=int, default=15071990)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--no-resume", action="store_true", help="Ignore completed per-TF Figure 6 shards")
    parser.add_argument("--compute-only", action="store_true", help="Compute/reuse TF shards without assembling the final export")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.metadata_dir / "figure6_celloracle_programme_deltas_by_cell.tsv.gz"
    report_path = args.metadata_dir / "figure6_celloracle_programme_delta_export_report.json"
    if out_path.exists() and report_path.exists() and not args.force:
        print(json.dumps({"status": "reused_existing", "output": str(out_path)}, ensure_ascii=False))
        return 0

    import celloracle as co

    started = time.time()
    availability = pd.read_csv(args.availability, sep="\t")
    tf_targets = pd.read_csv(args.tf_targets, sep="\t")
    state_specific = pd.read_csv(args.sox4_state_specific, sep="\t")
    programmes = build_frozen_programmes(availability, tf_targets, state_specific)

    oracle = co.load_hdf5(str(args.oracle))
    genes = list(map(str, oracle.adata.var_names))
    if not hasattr(oracle, "active_regulatory_genes"):
        oracle.extract_active_gene_lists(verbose=False)
    available_tfs = set(map(str, oracle.active_regulatory_genes)).intersection(genes)
    selected_tfs = [tf for tf in args.tfs if tf in available_tfs]
    missing_tfs = [tf for tf in args.tfs if tf not in available_tfs]

    baseline_matrix = np.asarray(oracle.adata.layers["simulation_input"])
    baseline_scores, programme_manifest = score_programmes(baseline_matrix, genes, programmes)
    obs_columns = [
        "dataset", "sample_id", "study_sample", "cnv_sample", "celloracle_state",
        "celloracle_main_strict", "driver_primary_module3_cnv_supported",
        "cnv_proxy_z", "proliferation_score_z",
        "driver_main_strict__module_Stressed_Injured",
        "driver_main_strict__module_Proliferation",
        "sample_disease_stage",
    ]
    obs = oracle.adata.obs.copy()
    meta = pd.DataFrame(index=obs.index.astype(str))
    for col in obs_columns:
        meta[col] = obs[col].to_numpy() if col in obs.columns else np.nan
    meta.insert(0, "cell_id", meta.index)
    for col in baseline_scores.columns:
        meta[f"baseline_{col}"] = baseline_scores[col].to_numpy()

    shard_dir = args.metadata_dir / "celloracle_programme_delta_shards"
    shard_dir.mkdir(parents=True, exist_ok=True)
    frames: list[pd.DataFrame] = []
    per_tf: list[dict[str, object]] = []
    for tf in selected_tfs:
        tf_start = time.time()
        shard_path = shard_dir / f"figure6_celloracle_{tf.lower()}_programme_deltas.tsv.gz"
        if shard_path.exists() and not args.no_resume:
            frame = pd.read_csv(shard_path, sep="\t")
            frames.append(frame)
            per_tf.append({
                "tf": tf,
                "condition": {tf: 0.0},
                "elapsed_seconds": round(time.time() - tf_start, 3),
                "mean_abs_delta_x": None,
                "status": "reused_completed_shard",
                "shard": str(shard_path.resolve()),
            })
            try:
                print(json.dumps(per_tf[-1], ensure_ascii=False), flush=True)
            except BrokenPipeError:
                pass
            continue
        oracle.simulate_shift(
            perturb_condition={tf: 0.0},
            GRN_unit="cluster",
            n_propagation=args.n_propagation,
            ignore_warning=True,
            clip_delta_X=True,
        )
        delta_x = np.asarray(oracle.adata.layers["delta_X"])
        delta_scores, _ = score_programmes(delta_x, genes, programmes)
        frame = meta.copy()
        frame.insert(0, "perturbation_type", "knockout")
        frame.insert(0, "tf", tf)
        for col in delta_scores.columns:
            frame[col] = delta_scores[col].to_numpy()
        # Persist each completed TF before progressing so a long run is resumable.
        frame.to_csv(shard_path, sep="\t", index=False, compression="gzip")
        frames.append(frame)
        per_tf.append(
            {
                "tf": tf,
                "condition": {tf: 0.0},
                "elapsed_seconds": round(time.time() - tf_start, 3),
                "mean_abs_delta_x": float(np.mean(np.abs(delta_x))),
                "status": "computed",
                "shard": str(shard_path.resolve()),
            }
        )
        try:
            print(json.dumps(per_tf[-1], ensure_ascii=False), flush=True)
        except BrokenPipeError:
            pass

    if args.compute_only:
        try:
            print(json.dumps({"status": "compute_only_complete", "tfs": selected_tfs, "n_shards": len(frames)}, ensure_ascii=False), flush=True)
        except BrokenPipeError:
            pass
        return 0

    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(out_path, sep="\t", index=False, compression="gzip")
    manifest_path = args.metadata_dir / "figure6_frozen_programme_gene_manifest.tsv"
    pd.DataFrame(programme_manifest).to_csv(manifest_path, sep="\t", index=False)
    report = {
        "module": "Figure 6 raw CellOracle delta-expression export",
        "method_scope": "raw KO delta-expression summaries only; R performs integration, statistics and plotting",
        "input_oracle": str(args.oracle.resolve()),
        "input_oracle_size_bytes": args.oracle.stat().st_size,
        "parameters": {
            "mode": "knockout",
            "expression_value": 0.0,
            "grn_unit": "cluster",
            "n_propagation": args.n_propagation,
            "clip_delta_x": True,
            "seed": args.seed,
            "embedding_recomputed": False,
            "execution_thread_limits": {
                key: os.environ.get(key)
                for key in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS"]
            },
        },
        "selected_tfs": selected_tfs,
        "missing_tfs": missing_tfs,
        "n_cells": int(oracle.adata.n_obs),
        "n_genes": int(oracle.adata.n_vars),
        "n_rows": int(combined.shape[0]),
        "programmes": programme_manifest,
        "outputs": {"cell_level": str(out_path.resolve()), "gene_manifest": str(manifest_path.resolve())},
        "per_tf": per_tf,
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "python": sys.version,
            "platform": platform.platform(),
            "celloracle": getattr(co, "__version__", None),
        },
        "guardrails": [
            "No restore or overexpression result was generated.",
            "No perturbation embedding or formal Figure 6 graphic was generated in Python.",
            "CNV-associated output is an expression signature and is not a virtual CNV burden estimate.",
        ],
    }
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    try:
        print(json.dumps({"status": "complete", "output": str(out_path), "n_rows": len(combined)}, ensure_ascii=False), flush=True)
    except BrokenPipeError:
        pass
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
