from __future__ import annotations

import argparse
import json
import time
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import cellrank as cr
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.sparse.csgraph import connected_components


ROOT = Path(__file__).resolve().parents[1]

LINEAGE_CNV = "cnv_supported_malignant"
LINEAGE_LATE_NON_CNV = "late_non_cnv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 6.2: CellRank fate probability and fate driver analysis.")
    parser.add_argument(
        "--input-h5ad",
        type=Path,
        default=ROOT / "data/processed/driver/driver_hepatocyte_trajectory.module6_1.h5ad",
    )
    parser.add_argument("--processed-dir", type=Path, default=ROOT / "data/processed/driver")
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/driver")
    parser.add_argument("--figures-dir", type=Path, default=ROOT / "figures/driver")
    parser.add_argument("--output-name", default="driver_cellrank_main_strict.module6_2.h5ad")
    parser.add_argument("--run-id", default="main_strict")
    parser.add_argument("--time-key", default="driver_main_strict__pseudotime_median")
    parser.add_argument("--phase-key", default="driver_main_strict__pseudotime_phase")
    parser.add_argument("--eligible-key", default="driver_main_strict__eligible")
    parser.add_argument("--terminal-key", default="trajectory_root_end_role")
    parser.add_argument("--cnv-supported-key", default="driver_primary_module3_cnv_supported")
    parser.add_argument("--use-rep", default="X_scANVI")
    parser.add_argument("--n-neighbors", type=int, default=50)
    parser.add_argument("--pseudotime-weight", type=float, default=0.8)
    parser.add_argument("--threshold-scheme", choices=["soft", "hard"], default="soft")
    parser.add_argument("--frac-to-keep", type=float, default=0.5)
    parser.add_argument("--solver", default="gmres")
    parser.add_argument("--seed", type=int, default=20260606)
    parser.add_argument("--top-n-drivers", type=int, default=50)
    parser.add_argument("--compression", default="gzip")
    return parser.parse_args()


def configure_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica"],
            "font.size": 8,
            "axes.labelsize": 8,
            "axes.titlesize": 9,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def bool_series(values: pd.Series) -> pd.Series:
    labels = values.astype("object").map(lambda value: "false" if pd.isna(value) else str(value).strip().lower())
    return labels.isin({"true", "1", "yes", "y"})


def write_dataframe(path: Path, df: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = "gzip" if path.suffix == ".gz" else None
    df.to_csv(path, sep="\t", index=False, compression=compression)


def select_driver_cells(adata: ad.AnnData, eligible_key: str) -> ad.AnnData:
    if eligible_key not in adata.obs.columns:
        raise KeyError(f"Missing eligible key in adata.obs: {eligible_key}")
    mask = bool_series(adata.obs[eligible_key]).to_numpy()
    if mask.sum() == 0:
        raise ValueError(f"No cells selected by {eligible_key}.")
    sub = adata[mask, :].copy()
    sub.uns.pop("iroot", None)
    return sub


def define_terminal_states(
    adata: ad.AnnData,
    terminal_key: str,
    phase_key: str,
    cnv_supported_key: str,
) -> pd.Series:
    for key in [terminal_key, phase_key, cnv_supported_key]:
        if key not in adata.obs.columns:
            raise KeyError(f"Missing terminal definition key in adata.obs: {key}")
    states = pd.Series(pd.NA, index=adata.obs_names, dtype="object")
    cnv_terminal = adata.obs[terminal_key].astype(str).eq("end_malignant_cnv")
    late_non_cnv = adata.obs[phase_key].astype(str).eq("late") & ~bool_series(adata.obs[cnv_supported_key])
    states.loc[cnv_terminal] = LINEAGE_CNV
    states.loc[late_non_cnv] = LINEAGE_LATE_NON_CNV
    states = pd.Series(pd.Categorical(states, categories=[LINEAGE_CNV, LINEAGE_LATE_NON_CNV]), index=adata.obs_names)
    counts = states.value_counts(dropna=False)
    if counts.get(LINEAGE_CNV, 0) == 0 or counts.get(LINEAGE_LATE_NON_CNV, 0) == 0:
        raise ValueError(f"Terminal states are incomplete: {counts.to_dict()}")
    return states


def connected_component_summary(adata: ad.AnnData, connectivities_key: str, terminal_states: pd.Series) -> pd.DataFrame:
    matrix = adata.obsp[connectivities_key]
    n_components, labels = connected_components(matrix, directed=False)
    obs = pd.DataFrame({"component": labels, "terminal_state": terminal_states.astype("object").to_numpy()}, index=adata.obs_names)
    rows = []
    for component, sub in obs.groupby("component", observed=True, sort=True):
        rows.append(
            {
                "component": int(component),
                "n_cells": int(sub.shape[0]),
                f"n_{LINEAGE_CNV}": int(sub["terminal_state"].eq(LINEAGE_CNV).sum()),
                f"n_{LINEAGE_LATE_NON_CNV}": int(sub["terminal_state"].eq(LINEAGE_LATE_NON_CNV).sum()),
            }
        )
    out = pd.DataFrame(rows).sort_values("n_cells", ascending=False).reset_index(drop=True)
    out.attrs["n_components"] = n_components
    return out


def lineage_matrix_to_df(adata: ad.AnnData, estimator: cr.estimators.GPCCA) -> pd.DataFrame:
    probs = np.asarray(estimator.fate_probabilities.X, dtype=float)
    names = list(map(str, estimator.fate_probabilities.names))
    df = pd.DataFrame(probs, columns=[f"cellrank_fate_prob_{name}" for name in names])
    df.insert(0, "cell_id", adata.obs_names.astype(str).to_numpy())
    return df


def drivers_to_long(drivers: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for lineage in [LINEAGE_CNV, LINEAGE_LATE_NON_CNV]:
        corr_col = f"{lineage}_corr"
        if corr_col not in drivers.columns:
            continue
        sub = pd.DataFrame(
            {
                "gene": drivers.index.astype(str),
                "lineage": lineage,
                "corr": pd.to_numeric(drivers.get(corr_col), errors="coerce"),
                "pval": pd.to_numeric(drivers.get(f"{lineage}_pval"), errors="coerce"),
                "qval": pd.to_numeric(drivers.get(f"{lineage}_qval"), errors="coerce"),
                "ci_low": pd.to_numeric(drivers.get(f"{lineage}_ci_low"), errors="coerce"),
                "ci_high": pd.to_numeric(drivers.get(f"{lineage}_ci_high"), errors="coerce"),
            }
        )
        sub["abs_corr"] = sub["corr"].abs()
        sub = sub.sort_values(["corr", "abs_corr"], ascending=[False, False]).reset_index(drop=True)
        sub["rank_positive_corr"] = np.arange(1, sub.shape[0] + 1)
        rows.append(sub)
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, axis=0, ignore_index=True)


def summarize_fate_by_group(cells: pd.DataFrame, group_cols: list[str], fate_col: str) -> pd.DataFrame:
    rows = []
    for group_col in group_cols:
        if group_col not in cells.columns:
            continue
        for value, sub in cells.groupby(group_col, observed=True, sort=True):
            fate = pd.to_numeric(sub[fate_col], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
            rows.append(
                {
                    "group_type": group_col,
                    "group": str(value),
                    "n_cells": int(sub.shape[0]),
                    "mean_cnv_fate_probability": float(fate.mean()) if not fate.empty else np.nan,
                    "median_cnv_fate_probability": float(fate.median()) if not fate.empty else np.nan,
                    "terminal_cnv_fraction": float(bool_series(sub.get("driver_primary_module3_cnv_supported", pd.Series(False))).mean())
                    if sub.shape[0] > 0
                    else np.nan,
                }
            )
    return pd.DataFrame(rows)


def plot_umap_fate(adata: ad.AnnData, fate_col: str, path_base: Path) -> list[str]:
    if "X_umap_global" in adata.obsm:
        xy = np.asarray(adata.obsm["X_umap_global"])
    elif "X_umap" in adata.obsm:
        xy = np.asarray(adata.obsm["X_umap"])
    else:
        return []
    values = pd.to_numeric(adata.obs[fate_col], errors="coerce").to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    order = np.argsort(np.nan_to_num(values, nan=-1.0))
    sca = ax.scatter(xy[order, 0], xy[order, 1], c=values[order], s=5, cmap="viridis", linewidths=0, rasterized=True)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title("CellRank CNV-supported malignant fate probability")
    cbar = fig.colorbar(sca, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Fate probability")
    fig.tight_layout()
    outputs = []
    for suffix in [".png", ".pdf"]:
        out = path_base.with_suffix(suffix)
        fig.savefig(out, bbox_inches="tight")
        outputs.append(str(out))
    plt.close(fig)
    return outputs


def plot_fate_by_pseudotime(cells: pd.DataFrame, time_key: str, fate_col: str, path_base: Path, n_bins: int = 20) -> list[str]:
    work = cells[["cell_id", time_key, fate_col]].copy()
    work[time_key] = pd.to_numeric(work[time_key], errors="coerce")
    work[fate_col] = pd.to_numeric(work[fate_col], errors="coerce")
    work = work.replace([np.inf, -np.inf], np.nan).dropna()
    if work.empty:
        return []
    ranked = work[time_key].rank(method="first")
    work["bin"] = pd.qcut(ranked, q=min(n_bins, work.shape[0]), labels=False, duplicates="drop")
    summary = (
        work.groupby("bin", observed=True)
        .agg(
            n_cells=(fate_col, "size"),
            mean_pseudotime=(time_key, "mean"),
            mean_fate_probability=(fate_col, "mean"),
            sem_fate_probability=(fate_col, "sem"),
        )
        .reset_index()
    )
    fig, ax = plt.subplots(figsize=(5.0, 3.4))
    ax.plot(summary["mean_pseudotime"], summary["mean_fate_probability"], marker="o", markersize=3, linewidth=1.8, color="#0072B2")
    lower = summary["mean_fate_probability"] - summary["sem_fate_probability"].fillna(0)
    upper = summary["mean_fate_probability"] + summary["sem_fate_probability"].fillna(0)
    ax.fill_between(summary["mean_pseudotime"], lower, upper, color="#56B4E9", alpha=0.25, linewidth=0)
    ax.set_xlabel("Main-strict consensus pseudotime")
    ax.set_ylabel("Mean CNV fate probability")
    ax.set_title("Fate probability increases along trajectory")
    fig.tight_layout()
    outputs = []
    for suffix in [".png", ".pdf"]:
        out = path_base.with_suffix(suffix)
        fig.savefig(out, bbox_inches="tight")
        outputs.append(str(out))
    plt.close(fig)
    return outputs


def plot_top_drivers(driver_long: pd.DataFrame, lineage: str, top_n: int, path_base: Path) -> list[str]:
    sub = driver_long.loc[driver_long["lineage"].eq(lineage)].copy()
    sub = sub.sort_values("corr", ascending=False).head(top_n)
    if sub.empty:
        return []
    sub = sub.iloc[::-1]
    fig, ax = plt.subplots(figsize=(5.2, max(3.0, 0.22 * top_n)))
    ax.barh(sub["gene"], sub["corr"], color="#D55E00", height=0.75)
    ax.set_xlabel("Correlation with CNV fate probability")
    ax.set_ylabel("")
    ax.set_title(f"Top {top_n} CellRank CNV fate drivers")
    fig.tight_layout()
    outputs = []
    for suffix in [".png", ".pdf"]:
        out = path_base.with_suffix(suffix)
        fig.savefig(out, bbox_inches="tight")
        outputs.append(str(out))
    plt.close(fig)
    return outputs


def drop_nonserializable_cellrank_slots(adata: ad.AnnData) -> list[str]:
    removed = []
    for key in list(adata.obsm.keys()):
        value = adata.obsm[key]
        if key.startswith("lineages") or "cellrank" in type(value).__module__:
            removed.append(f"obsm/{key}")
            del adata.obsm[key]
            continue
        if isinstance(value, (np.ndarray, pd.DataFrame)) or sparse.issparse(value):
            continue
        removed.append(f"obsm/{key}")
        del adata.obsm[key]
    for key in list(adata.varm.keys()):
        value = adata.varm[key]
        if key.startswith("lineages") or "cellrank" in type(value).__module__:
            removed.append(f"varm/{key}")
            del adata.varm[key]
            continue
        if isinstance(value, (np.ndarray, pd.DataFrame)) or sparse.issparse(value):
            continue
        removed.append(f"varm/{key}")
        del adata.varm[key]
    return removed


def main() -> None:
    start = time.time()
    args = parse_args()
    configure_plot_style()
    args.processed_dir.mkdir(parents=True, exist_ok=True)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    args.figures_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.input_h5ad)
    sub = select_driver_cells(adata, args.eligible_key)
    if args.time_key not in sub.obs.columns:
        raise KeyError(f"Missing pseudotime key in adata.obs: {args.time_key}")
    if args.use_rep not in sub.obsm:
        raise KeyError(f"Missing representation in adata.obsm: {args.use_rep}")

    terminal_states = define_terminal_states(sub, args.terminal_key, args.phase_key, args.cnv_supported_key)
    sub.obs["cellrank_terminal_state"] = terminal_states

    sc.pp.neighbors(sub, use_rep=args.use_rep, n_neighbors=args.n_neighbors, key_added="cellrank_neighbors")
    component_summary = connected_component_summary(sub, "cellrank_neighbors_connectivities", terminal_states)
    if component_summary.attrs["n_components"] != 1:
        valid_components = component_summary.loc[
            component_summary[f"n_{LINEAGE_CNV}"].gt(0) & component_summary[f"n_{LINEAGE_LATE_NON_CNV}"].gt(0),
            "component",
        ]
        if valid_components.empty:
            raise ValueError("No connected component contains both terminal states for CellRank fate probability.")
        keep_components = set(valid_components.astype(int))
        labels = connected_components(sub.obsp["cellrank_neighbors_connectivities"], directed=False)[1]
        keep = np.array([label in keep_components for label in labels])
        sub = sub[keep, :].copy()
        terminal_states = terminal_states.loc[sub.obs_names]
        sub.obs["cellrank_terminal_state"] = terminal_states
        sc.pp.neighbors(sub, use_rep=args.use_rep, n_neighbors=args.n_neighbors, key_added="cellrank_neighbors")
        component_summary = connected_component_summary(sub, "cellrank_neighbors_connectivities", terminal_states)

    pseudotime_kernel = cr.kernels.PseudotimeKernel(
        sub,
        time_key=args.time_key,
        conn_key="cellrank_neighbors_connectivities",
    )
    pseudotime_kernel.compute_transition_matrix(
        threshold_scheme=args.threshold_scheme,
        frac_to_keep=args.frac_to_keep,
        n_jobs=1,
        show_progress_bar=False,
    )
    connectivity_kernel = cr.kernels.ConnectivityKernel(sub, conn_key="cellrank_neighbors_connectivities")
    connectivity_kernel.compute_transition_matrix()
    kernel = args.pseudotime_weight * pseudotime_kernel + (1.0 - args.pseudotime_weight) * connectivity_kernel

    estimator = cr.estimators.GPCCA(kernel)
    estimator.set_terminal_states(terminal_states)
    estimator.compute_fate_probabilities(
        solver=args.solver,
        use_petsc=False,
        n_jobs=1,
        show_progress_bar=False,
        tol=1e-6,
    )
    fate_probabilities = lineage_matrix_to_df(sub, estimator)
    for column in fate_probabilities.columns:
        if column == "cell_id":
            continue
        sub.obs[column] = fate_probabilities.set_index("cell_id").loc[sub.obs_names, column].to_numpy(dtype=float)
    sub.obsm["cellrank_fate_probabilities"] = fate_probabilities.drop(columns=["cell_id"]).to_numpy(dtype=float)
    sub.uns["cellrank_fate_probability_names"] = [
        column.replace("cellrank_fate_prob_", "") for column in fate_probabilities.columns if column != "cell_id"
    ]

    drivers = estimator.compute_lineage_drivers(method="fisher", layer=None, use_raw=False, show_progress_bar=False)
    driver_long = drivers_to_long(drivers)

    fate_cells = pd.DataFrame({"cell_id": sub.obs_names.astype(str)})
    keep_obs = [
        "dataset",
        "sample_id",
        "study_sample",
        "cnv_sample",
        "sample_source_class",
        "cell_disease_stage",
        "trajectory_root_end_role",
        args.time_key,
        args.phase_key,
        args.cnv_supported_key,
        "driver_primary_cnv_evidence_tier",
        "cellrank_terminal_state",
    ]
    for col in keep_obs:
        if col in sub.obs.columns:
            fate_cells[col] = sub.obs[col].astype(str).to_numpy() if str(sub.obs[col].dtype) == "category" else sub.obs[col].to_numpy()
    fate_cells = fate_cells.merge(fate_probabilities, on="cell_id", how="left")
    fate_col = f"cellrank_fate_prob_{LINEAGE_CNV}"
    fate_summary = summarize_fate_by_group(
        fate_cells,
        group_cols=["cell_disease_stage", args.phase_key, "dataset", "sample_id", "trajectory_root_end_role"],
        fate_col=fate_col,
    )

    output_h5ad = args.processed_dir / args.output_name
    sub.uns["module6_2_cellrank"] = {
        "module": "6.2",
        "method": "CellRank PseudotimeKernel plus ConnectivityKernel fate probability and lineage driver analysis",
        "run_id": args.run_id,
        "time_key": args.time_key,
        "terminal_states": [LINEAGE_CNV, LINEAGE_LATE_NON_CNV],
        "pseudotime_weight": float(args.pseudotime_weight),
        "connectivity_weight": float(1.0 - args.pseudotime_weight),
        "n_neighbors": int(args.n_neighbors),
        "threshold_scheme": args.threshold_scheme,
        "frac_to_keep": float(args.frac_to_keep),
        "solver": args.solver,
    }
    removed_cellrank_slots = drop_nonserializable_cellrank_slots(sub)
    sub.write_h5ad(output_h5ad, compression=args.compression)

    fate_path = args.metadata_dir / "driver_module6_2_cellrank_fate_probabilities.tsv.gz"
    drivers_path = args.metadata_dir / "driver_module6_2_cellrank_lineage_drivers.tsv.gz"
    top_drivers_path = args.metadata_dir / "driver_module6_2_cellrank_top_cnv_fate_drivers.tsv"
    summary_path = args.metadata_dir / "driver_module6_2_cellrank_fate_summary.tsv"
    terminals_path = args.metadata_dir / "driver_module6_2_cellrank_terminal_states.tsv"
    components_path = args.metadata_dir / "driver_module6_2_cellrank_components.tsv"

    write_dataframe(fate_path, fate_cells)
    write_dataframe(drivers_path, driver_long)
    top_cnv = driver_long.loc[driver_long["lineage"].eq(LINEAGE_CNV)].sort_values("corr", ascending=False).head(args.top_n_drivers)
    write_dataframe(top_drivers_path, top_cnv)
    write_dataframe(summary_path, fate_summary)
    terminal_table = fate_cells.loc[fate_cells["cellrank_terminal_state"].astype(str).isin([LINEAGE_CNV, LINEAGE_LATE_NON_CNV])].copy()
    write_dataframe(terminals_path, terminal_table)
    write_dataframe(components_path, component_summary)

    figure_outputs: list[str] = []
    figure_outputs += plot_umap_fate(sub, fate_col, args.figures_dir / "driver_module6_2_cellrank_cnv_fate_umap")
    figure_outputs += plot_fate_by_pseudotime(
        fate_cells,
        time_key=args.time_key,
        fate_col=fate_col,
        path_base=args.figures_dir / "driver_module6_2_cellrank_cnv_fate_by_pseudotime",
    )
    figure_outputs += plot_top_drivers(
        driver_long,
        lineage=LINEAGE_CNV,
        top_n=min(25, args.top_n_drivers),
        path_base=args.figures_dir / "driver_module6_2_cellrank_top_cnv_fate_drivers",
    )

    fate_matrix = fate_probabilities.drop(columns=["cell_id"]).to_numpy(dtype=float)
    fate_sum = fate_matrix.sum(axis=1)
    report = {
        "module": "6.2",
        "method": "CellRank fate probability and fate driver analysis",
        "input_h5ad": str(args.input_h5ad),
        "output_h5ad": str(output_h5ad),
        "run_id": args.run_id,
        "n_cells": int(sub.n_obs),
        "n_genes": int(sub.n_vars),
        "terminal_state_counts": terminal_states.value_counts(dropna=False).astype(int).to_dict(),
        "component_summary": component_summary.to_dict(orient="records"),
        "fate_probability_qc": {
            "min_row_sum": float(np.nanmin(fate_sum)),
            "max_row_sum": float(np.nanmax(fate_sum)),
            "mean_cnv_fate_probability": float(fate_probabilities[fate_col].mean()),
            "median_cnv_fate_probability": float(fate_probabilities[fate_col].median()),
        },
        "removed_nonserializable_cellrank_slots": removed_cellrank_slots,
        "top_cnv_fate_drivers": top_cnv.head(15).to_dict(orient="records"),
        "outputs": {
            "h5ad": str(output_h5ad),
            "fate_probabilities": str(fate_path),
            "lineage_drivers": str(drivers_path),
            "top_cnv_fate_drivers": str(top_drivers_path),
            "fate_summary": str(summary_path),
            "terminal_states": str(terminals_path),
            "components": str(components_path),
            "figures": figure_outputs,
        },
        "package_versions": {
            "cellrank": version("cellrank"),
            "scanpy": version("scanpy"),
            "anndata": version("anndata"),
            "pandas": version("pandas"),
            "numpy": version("numpy"),
        },
        "elapsed_seconds": round(time.time() - start, 3),
    }
    report_path = args.metadata_dir / "driver_module6_2_cellrank_report.json"
    report["outputs"]["report"] = str(report_path)
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
