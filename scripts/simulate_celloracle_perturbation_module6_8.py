from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

try:
    from qc_celloracle_inputs_module6_5 import read_tf_list
except ModuleNotFoundError:
    from scripts.qc_celloracle_inputs_module6_5 import read_tf_list


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ORACLE = PROJECT_ROOT / "data/processed/driver/celloracle_module6_7/celloracle_module6_7_fitted.celloracle.oracle"
DEFAULT_TF_LIST = PROJECT_ROOT / "metadata/driver/celloracle_input_tfs.module6_4.txt"
DEFAULT_OUT_DIR = PROJECT_ROOT / "data/processed/driver/celloracle_module6_8"
DEFAULT_METADATA_DIR = PROJECT_ROOT / "metadata/driver"


def select_available_tfs(
    input_tfs: Iterable[str],
    available_genes: set[str],
    active_regulatory_genes: set[str],
) -> tuple[list[str], list[str]]:
    selected = []
    missing = []
    for tf in input_tfs:
        if tf in available_genes and tf in active_regulatory_genes:
            selected.append(tf)
        else:
            missing.append(tf)
    return selected, missing


def build_perturbation_condition(
    tf: str,
    mode: str,
    expression_value: float | None,
) -> dict[str, float]:
    if mode == "knockout":
        value = 0.0
    elif mode == "fixed":
        if expression_value is None:
            raise ValueError("--expression-value is required when --mode fixed")
        value = float(expression_value)
    else:
        raise ValueError(f"Unsupported perturbation mode: {mode}")
    if value < 0:
        raise ValueError("CellOracle perturbation expression value must be non-negative")
    return {tf: value}


def compute_malignant_axis(
    embedding: np.ndarray,
    states: pd.Series,
    start_state: str,
    end_state: str,
) -> np.ndarray:
    start_mask = states.astype(str).to_numpy() == start_state
    end_mask = states.astype(str).to_numpy() == end_state
    if not start_mask.any() or not end_mask.any():
        raise ValueError(f"Cannot compute malignant axis from {start_state!r} to {end_state!r}")
    vector = embedding[end_mask].mean(axis=0) - embedding[start_mask].mean(axis=0)
    norm = np.linalg.norm(vector)
    if norm == 0:
        raise ValueError("Malignant axis has zero length")
    return vector / norm


def _neighbor_indices_from_graph(embedding_knn) -> np.ndarray:
    graph = embedding_knn.tocsr()
    row_lengths = np.diff(graph.indptr)
    if len(set(row_lengths.tolist())) != 1:
        raise ValueError("Embedding KNN graph has unequal row lengths")
    return graph.indices.reshape((graph.shape[0], int(row_lengths[0])))


def compute_sparse_embedding_shift(
    embedding: np.ndarray,
    corrcoef: np.ndarray,
    neighbor_indices: np.ndarray,
    sigma_corr: float,
) -> np.ndarray:
    shifts = np.zeros_like(embedding, dtype=float)
    for cell_ix, neigh_ixs in enumerate(neighbor_indices):
        vectors = embedding[neigh_ixs] - embedding[cell_ix]
        norms = np.linalg.norm(vectors, axis=1)
        valid = norms > 0
        if not np.any(valid):
            continue
        unit_vectors = vectors[valid] / norms[valid, None]
        weights = np.exp(corrcoef[cell_ix, neigh_ixs[valid]] / sigma_corr)
        weight_sum = weights.sum()
        if weight_sum == 0 or not np.isfinite(weight_sum):
            continue
        weighted = (unit_vectors * (weights / weight_sum)[:, None]).sum(axis=0)
        baseline = unit_vectors.mean(axis=0)
        shifts[cell_ix] = weighted - baseline
    return shifts


def summarize_cell_shifts(
    tf: str,
    obs: pd.DataFrame,
    delta_embedding: np.ndarray | None,
    delta_x: np.ndarray,
    malignant_axis: np.ndarray | None,
) -> pd.DataFrame:
    summary = pd.DataFrame(
        {
            "tf": tf,
            "cell_id": obs.index.astype(str),
            "celloracle_state": obs["celloracle_state"].astype(str).to_numpy(),
            "mean_abs_delta_x": np.mean(np.abs(delta_x), axis=1),
            "mean_delta_x": np.mean(delta_x, axis=1),
        }
    )
    if delta_embedding is not None:
        summary["delta_embedding_1"] = delta_embedding[:, 0]
        summary["delta_embedding_2"] = delta_embedding[:, 1]
        summary["embedding_shift_norm"] = np.linalg.norm(delta_embedding, axis=1)
        if malignant_axis is not None:
            summary["malignant_axis_projection"] = delta_embedding @ malignant_axis
    return summary


def summarize_state_shifts(cell_summary: pd.DataFrame) -> pd.DataFrame:
    numeric_cols = [
        col
        for col in [
            "mean_abs_delta_x",
            "mean_delta_x",
            "delta_embedding_1",
            "delta_embedding_2",
            "embedding_shift_norm",
            "malignant_axis_projection",
        ]
        if col in cell_summary.columns
    ]
    grouped = cell_summary.groupby(["tf", "celloracle_state"], observed=True)
    summary = grouped[numeric_cols].agg(["mean", "median"]).reset_index()
    summary.columns = [
        "_".join([str(part) for part in col if part])
        if isinstance(col, tuple)
        else str(col)
        for col in summary.columns
    ]
    summary.insert(2, "n_cells", grouped.size().to_numpy())
    return summary


def summarize_perturbation_ranking(state_summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    state_specs = [
        ("malignant_or_malignant_like", "malignant"),
        ("normal_reference", "normal"),
        ("regenerative_progenitor", "regenerative"),
        ("stressed_injured", "stressed"),
        ("proliferating_candidate", "proliferating"),
    ]
    for tf, df in state_summary.groupby("tf", sort=False):
        total_cells = int(df["n_cells"].sum())
        row = {
            "tf": tf,
            "n_cells_total": total_cells,
            "weighted_mean_abs_delta_x": float(
                (df["mean_abs_delta_x_mean"] * df["n_cells"]).sum() / total_cells
            ),
        }
        for state_name, prefix in state_specs:
            state_df = df.loc[df["celloracle_state"].astype(str) == state_name]
            if state_df.empty:
                continue
            row[f"{prefix}_axis_projection_mean"] = float(state_df["malignant_axis_projection_mean"].iloc[0])
            row[f"{prefix}_embedding_shift_norm_mean"] = float(state_df["embedding_shift_norm_mean"].iloc[0])
            row[f"{prefix}_mean_abs_delta_x"] = float(state_df["mean_abs_delta_x_mean"].iloc[0])
        row["anti_malignant_shift_score"] = -row.get("malignant_axis_projection_mean", np.nan)
        rows.append(row)
    ranking = pd.DataFrame(rows).sort_values(
        ["anti_malignant_shift_score", "weighted_mean_abs_delta_x"],
        ascending=[False, False],
    ).reset_index(drop=True)
    ranking.insert(0, "rank", np.arange(1, len(ranking) + 1))
    return ranking


def summarize_gene_delta_by_state(
    tf: str,
    genes: list[str],
    states: pd.Series,
    delta_x: np.ndarray,
    top_n: int,
) -> pd.DataFrame:
    rows = []
    state_values = states.astype(str).to_numpy()
    genes_array = np.asarray(genes)
    for state in sorted(pd.unique(state_values)):
        mask = state_values == state
        if not mask.any():
            continue
        mean_delta = np.asarray(delta_x[mask].mean(axis=0)).ravel()
        mean_abs_delta = np.abs(mean_delta)
        top_ixs = np.argsort(mean_abs_delta)[::-1][:top_n]
        for ix in top_ixs:
            rows.append(
                {
                    "tf": tf,
                    "celloracle_state": state,
                    "gene": str(genes_array[ix]),
                    "mean_delta_x": float(mean_delta[ix]),
                    "abs_mean_delta_x": float(mean_abs_delta[ix]),
                }
            )
    return pd.DataFrame(rows)


def summarize_target_tf_delta(
    tf: str,
    genes: list[str],
    states: pd.Series,
    delta_x: np.ndarray,
    target_tfs: list[str],
) -> pd.DataFrame:
    gene_to_ix = {gene: ix for ix, gene in enumerate(genes)}
    rows = []
    state_values = states.astype(str).to_numpy()
    for state in sorted(pd.unique(state_values)):
        mask = state_values == state
        for target_tf in target_tfs:
            if target_tf not in gene_to_ix:
                continue
            values = delta_x[mask, gene_to_ix[target_tf]]
            rows.append(
                {
                    "tf": tf,
                    "target_tf": target_tf,
                    "celloracle_state": state,
                    "mean_delta_x": float(np.mean(values)),
                    "median_delta_x": float(np.median(values)),
                    "mean_abs_delta_x": float(np.mean(np.abs(values))),
                }
            )
    return pd.DataFrame(rows)


def summarize_grid_arrows(tf: str, oracle) -> pd.DataFrame:
    if not hasattr(oracle, "flow_grid") or getattr(oracle, "flow_grid") is None:
        return pd.DataFrame()
    flow_grid = np.asarray(oracle.flow_grid)
    flow = np.asarray(oracle.flow)
    flow_norm = np.asarray(getattr(oracle, "flow_norm", flow))
    total_p_mass = np.asarray(getattr(oracle, "total_p_mass", np.full(flow_grid.shape[0], np.nan)))
    return pd.DataFrame(
        {
            "tf": tf,
            "grid_x": flow_grid[:, 0],
            "grid_y": flow_grid[:, 1],
            "flow_x": flow[:, 0],
            "flow_y": flow[:, 1],
            "flow_norm_x": flow_norm[:, 0],
            "flow_norm_y": flow_norm[:, 1],
            "flow_magnitude": np.linalg.norm(flow, axis=1),
            "total_p_mass": total_p_mass,
        }
    )


def run_perturbation_simulation(
    oracle_path: Path,
    tf_list_path: Path,
    out_dir: Path,
    metadata_dir: Path,
    selected_tfs: list[str] | None,
    mode: str,
    expression_value: float | None,
    grn_unit: str,
    n_propagation: int,
    clip_delta_x: bool,
    ignore_warning: bool,
    calculate_embedding: bool,
    n_neighbors: int,
    sampled_fraction: float,
    sigma_corr: float,
    n_jobs: int,
    threads: int | None,
    calculate_grid: bool,
    grid_steps: int,
    grid_neighbors: int,
    top_genes_per_state: int,
    save_per_tf_oracle: bool,
    seed: int,
) -> dict:
    import celloracle as co

    input_tfs = read_tf_list(tf_list_path)
    if selected_tfs:
        selected_set = set(selected_tfs)
        input_tfs = [tf for tf in input_tfs if tf in selected_set]

    out_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)

    oracle = co.load_hdf5(str(oracle_path))
    available_genes = set(map(str, oracle.adata.var.index))
    if not hasattr(oracle, "active_regulatory_genes"):
        oracle.extract_active_gene_lists(verbose=False)
    active_regulatory_genes = set(map(str, oracle.active_regulatory_genes))
    perturb_tfs, skipped_tfs = select_available_tfs(input_tfs, available_genes, active_regulatory_genes)
    if not perturb_tfs:
        raise ValueError("No perturbation TFs are available in the fitted CellOracle object")

    obs = oracle.adata.obs.copy()
    genes = list(map(str, oracle.adata.var.index))
    states = obs["celloracle_state"].astype(str)
    embedding = np.asarray(oracle.adata.obsm[oracle.embedding_name])
    malignant_axis = compute_malignant_axis(
        embedding,
        states,
        start_state="normal_reference",
        end_state="malignant_or_malignant_like",
    )

    cell_frames = []
    gene_frames = []
    tf_delta_frames = []
    grid_frames = []
    per_tf_reports = []

    for i, tf in enumerate(perturb_tfs):
        tf_started = datetime.now(timezone.utc)
        condition = build_perturbation_condition(tf, mode=mode, expression_value=expression_value)
        oracle.simulate_shift(
            perturb_condition=condition,
            GRN_unit=grn_unit,
            n_propagation=n_propagation,
            ignore_warning=ignore_warning,
            clip_delta_X=clip_delta_x,
        )
        delta_x = np.asarray(oracle.adata.layers["delta_X"])
        delta_embedding = None
        if calculate_embedding:
            oracle.estimate_transition_prob(
                n_neighbors=n_neighbors,
                knn_random=True,
                sampled_fraction=sampled_fraction,
                n_jobs=n_jobs,
                threads=threads,
                calculate_randomized=False,
                random_seed=seed + i,
            )
            neighbor_indices = _neighbor_indices_from_graph(oracle.embedding_knn)
            delta_embedding = compute_sparse_embedding_shift(
                embedding=np.asarray(oracle.embedding),
                corrcoef=np.asarray(oracle.corrcoef),
                neighbor_indices=neighbor_indices,
                sigma_corr=sigma_corr,
            )
            oracle.delta_embedding = delta_embedding
            if calculate_grid:
                oracle.calculate_grid_arrows(
                    smooth=0.5,
                    steps=(grid_steps, grid_steps),
                    n_neighbors=grid_neighbors,
                    n_jobs=n_jobs,
                )
                grid_frames.append(summarize_grid_arrows(tf, oracle))

        cell_frames.append(
            summarize_cell_shifts(
                tf=tf,
                obs=obs,
                delta_embedding=delta_embedding,
                delta_x=delta_x,
                malignant_axis=malignant_axis if delta_embedding is not None else None,
            )
        )
        gene_frames.append(
            summarize_gene_delta_by_state(
                tf=tf,
                genes=genes,
                states=states,
                delta_x=delta_x,
                top_n=top_genes_per_state,
            )
        )
        tf_delta_frames.append(
            summarize_target_tf_delta(
                tf=tf,
                genes=genes,
                states=states,
                delta_x=delta_x,
                target_tfs=perturb_tfs,
            )
        )
        per_tf_oracle_path = None
        if save_per_tf_oracle:
            per_tf_oracle_path = out_dir / f"celloracle_module6_8_{tf}_{mode}.celloracle.oracle"
            oracle.to_hdf5(str(per_tf_oracle_path))

        tf_finished = datetime.now(timezone.utc)
        per_tf_reports.append(
            {
                "tf": tf,
                "condition": condition,
                "elapsed_seconds": (tf_finished - tf_started).total_seconds(),
                "per_tf_oracle": str(per_tf_oracle_path) if per_tf_oracle_path else None,
                "mean_abs_delta_x": float(np.mean(np.abs(delta_x))),
            }
        )
        print(json.dumps(per_tf_reports[-1], ensure_ascii=False), flush=True)

    cell_summary = pd.concat(cell_frames, axis=0, ignore_index=True)
    state_summary = summarize_state_shifts(cell_summary)
    ranking = summarize_perturbation_ranking(state_summary)
    gene_summary = pd.concat(gene_frames, axis=0, ignore_index=True)
    tf_delta_summary = pd.concat(tf_delta_frames, axis=0, ignore_index=True)
    grid_summary = pd.concat(grid_frames, axis=0, ignore_index=True) if grid_frames else pd.DataFrame()

    cell_summary_path = metadata_dir / "celloracle_module6_8_cell_shift_summary.tsv.gz"
    state_summary_path = metadata_dir / "celloracle_module6_8_state_shift_summary.tsv"
    ranking_path = metadata_dir / "celloracle_module6_8_perturbation_ranking.tsv"
    gene_summary_path = metadata_dir / "celloracle_module6_8_top_gene_delta_by_state.tsv.gz"
    tf_delta_summary_path = metadata_dir / "celloracle_module6_8_tf_delta_summary.tsv"
    grid_summary_path = metadata_dir / "celloracle_module6_8_grid_arrows.tsv.gz"

    cell_summary.to_csv(cell_summary_path, sep="\t", index=False)
    state_summary.to_csv(state_summary_path, sep="\t", index=False)
    ranking.to_csv(ranking_path, sep="\t", index=False)
    gene_summary.to_csv(gene_summary_path, sep="\t", index=False)
    tf_delta_summary.to_csv(tf_delta_summary_path, sep="\t", index=False)
    if len(grid_summary):
        grid_summary.to_csv(grid_summary_path, sep="\t", index=False)

    return {
        "celloracle_version": getattr(co, "__version__", None),
        "input_oracle": str(oracle_path),
        "n_cells": int(oracle.adata.n_obs),
        "n_genes": int(oracle.adata.n_vars),
        "embedding_name": str(oracle.embedding_name),
        "n_input_tfs": int(len(input_tfs)),
        "n_perturbed_tfs": int(len(perturb_tfs)),
        "perturbed_tfs": perturb_tfs,
        "skipped_tfs": skipped_tfs,
        "cell_summary": str(cell_summary_path),
        "state_summary": str(state_summary_path),
        "perturbation_ranking": str(ranking_path),
        "top_gene_delta_by_state": str(gene_summary_path),
        "tf_delta_summary": str(tf_delta_summary_path),
        "grid_arrows": str(grid_summary_path) if len(grid_summary) else None,
        "per_tf_reports": per_tf_reports,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 6.8 CellOracle TF perturbation simulation")
    parser.add_argument("--oracle", type=Path, default=DEFAULT_ORACLE)
    parser.add_argument("--tf-list", type=Path, default=DEFAULT_TF_LIST)
    parser.add_argument("--tfs", nargs="*", default=None)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA_DIR)
    parser.add_argument("--mode", choices=["knockout", "fixed"], default="knockout")
    parser.add_argument("--expression-value", type=float, default=None)
    parser.add_argument("--grn-unit", choices=["cluster", "whole"], default="cluster")
    parser.add_argument("--n-propagation", type=int, default=3)
    parser.add_argument("--clip-delta-x", action="store_true", default=True)
    parser.add_argument("--no-clip-delta-x", dest="clip_delta_x", action="store_false")
    parser.add_argument("--ignore-warning", action="store_true", default=True)
    parser.add_argument("--strict-warning", dest="ignore_warning", action="store_false")
    parser.add_argument("--calculate-embedding", action="store_true", default=True)
    parser.add_argument("--skip-embedding", dest="calculate_embedding", action="store_false")
    parser.add_argument("--n-neighbors", type=int, default=200)
    parser.add_argument("--sampled-fraction", type=float, default=0.3)
    parser.add_argument("--sigma-corr", type=float, default=0.05)
    parser.add_argument("--n-jobs", type=int, default=8)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--calculate-grid", action="store_true", default=True)
    parser.add_argument("--skip-grid", dest="calculate_grid", action="store_false")
    parser.add_argument("--grid-steps", type=int, default=40)
    parser.add_argument("--grid-neighbors", type=int, default=100)
    parser.add_argument("--top-genes-per-state", type=int, default=50)
    parser.add_argument("--save-per-tf-oracle", action="store_true")
    parser.add_argument("--seed", type=int, default=15071990)
    parser.add_argument("--report", type=Path, default=DEFAULT_METADATA_DIR / "celloracle_module6_8_perturbation_report.json")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = datetime.now(timezone.utc)
    result = run_perturbation_simulation(
        oracle_path=args.oracle,
        tf_list_path=args.tf_list,
        out_dir=args.out_dir,
        metadata_dir=args.metadata_dir,
        selected_tfs=args.tfs,
        mode=args.mode,
        expression_value=args.expression_value,
        grn_unit=args.grn_unit,
        n_propagation=args.n_propagation,
        clip_delta_x=args.clip_delta_x,
        ignore_warning=args.ignore_warning,
        calculate_embedding=args.calculate_embedding,
        n_neighbors=args.n_neighbors,
        sampled_fraction=args.sampled_fraction,
        sigma_corr=args.sigma_corr,
        n_jobs=args.n_jobs,
        threads=args.threads,
        calculate_grid=args.calculate_grid,
        grid_steps=args.grid_steps,
        grid_neighbors=args.grid_neighbors,
        top_genes_per_state=args.top_genes_per_state,
        save_per_tf_oracle=args.save_per_tf_oracle,
        seed=args.seed,
    )
    finished = datetime.now(timezone.utc)
    report = {
        "module": "6.8",
        "method": "CellOracle TF perturbation simulation",
        "created_at_utc": finished.isoformat(),
        "elapsed_seconds": (finished - started).total_seconds(),
        "parameters": {
            "mode": args.mode,
            "expression_value": args.expression_value,
            "grn_unit": args.grn_unit,
            "n_propagation": args.n_propagation,
            "clip_delta_x": args.clip_delta_x,
            "ignore_warning": args.ignore_warning,
            "calculate_embedding": args.calculate_embedding,
            "n_neighbors": args.n_neighbors,
            "sampled_fraction": args.sampled_fraction,
            "sigma_corr": args.sigma_corr,
            "n_jobs": args.n_jobs,
            "threads": args.threads,
            "calculate_grid": args.calculate_grid,
            "grid_steps": args.grid_steps,
            "grid_neighbors": args.grid_neighbors,
            "top_genes_per_state": args.top_genes_per_state,
            "save_per_tf_oracle": args.save_per_tf_oracle,
            "seed": args.seed,
        },
        "result": result,
        "python_runtime": {
            "executable": sys.executable,
            "version": sys.version,
            "platform": platform.platform(),
        },
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(
        {
            "n_perturbed_tfs": result["n_perturbed_tfs"],
            "perturbed_tfs": result["perturbed_tfs"],
            "state_summary": result["state_summary"],
            "cell_summary": result["cell_summary"],
            "report": str(args.report),
        },
        ensure_ascii=False,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
