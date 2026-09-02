from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA_DIR = PROJECT_ROOT / "metadata/driver"
DEFAULT_FIGURE_DIR = PROJECT_ROOT / "figures/driver"
DEFAULT_H5AD = PROJECT_ROOT / "data/processed/driver/celloracle_module6_6/celloracle_module6_6_input.h5ad"

STATE_ORDER = [
    "normal_reference",
    "stressed_injured",
    "regenerative_progenitor",
    "proliferating_candidate",
    "malignant_or_malignant_like",
]

STATE_LABELS = {
    "normal_reference": "Normal ref.",
    "stressed_injured": "Stressed",
    "regenerative_progenitor": "Regenerative",
    "proliferating_candidate": "Proliferating",
    "malignant_or_malignant_like": "Malignant-like",
}

EVIDENCE_COLUMNS = [
    "anti_malignant_shift_score",
    "weighted_mean_abs_delta_x",
    "total_score",
    "motif_score",
    "fate_score",
    "cnv_fate_pearson_r",
    "tf_self_cellrank_corr",
    "malignant_grn_edges_passing_p",
    "malignant_grn_mean_coef_abs",
]


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 7,
            "axes.labelsize": 7,
            "axes.titlesize": 8,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.6,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )
    sns.set_theme(style="ticks", rc=mpl.rcParams)


def save_pub_figure(fig: plt.Figure, stem: Path, dpi: int = 600) -> list[str]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    paths = []
    for suffix in [".svg", ".pdf", ".png", ".tiff"]:
        path = stem.with_suffix(suffix)
        kwargs = {"bbox_inches": "tight"}
        if suffix in {".png", ".tiff"}:
            kwargs["dpi"] = dpi
        fig.savefig(path, **kwargs)
        paths.append(str(path))
    plt.close(fig)
    return paths


def minmax_scale(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce").astype(float)
    valid = numeric.dropna()
    if valid.empty:
        return pd.Series(np.zeros(len(numeric)), index=numeric.index)
    missing = numeric.isna()
    vmin = valid.min()
    vmax = valid.max()
    if vmax == vmin:
        scaled = pd.Series(np.full(len(numeric), 0.5), index=numeric.index)
    else:
        scaled = (numeric - vmin) / (vmax - vmin)
    scaled.loc[missing] = 0.0
    return scaled.fillna(0.0)


def build_candidate_evidence_matrix(
    ranking: pd.DataFrame,
    selection: pd.DataFrame,
    grn_tf_summary: pd.DataFrame,
) -> pd.DataFrame:
    grn_malignant = grn_tf_summary.loc[
        grn_tf_summary["celloracle_state"].astype(str) == "malignant_or_malignant_like",
        ["tf", "n_edges_passing_p", "mean_coef_abs_passing_p"],
    ].rename(
        columns={
            "n_edges_passing_p": "malignant_grn_edges_passing_p",
            "mean_coef_abs_passing_p": "malignant_grn_mean_coef_abs",
        }
    )

    selection_cols = [
        "tf",
        "role",
        "tier",
        "total_score",
        "motif_score",
        "fate_score",
        "cellrank_overlap_score",
        "robustness_score",
        "compatibility_score",
        "biology_score",
        "cnv_fate_pearson_r",
        "tf_self_cellrank_corr",
        "base_grn_target_genes",
        "selected_for_main_panel",
    ]
    available_selection_cols = [col for col in selection_cols if col in selection.columns]
    evidence = ranking.merge(selection[available_selection_cols], on="tf", how="left")
    evidence = evidence.merge(grn_malignant, on="tf", how="left")

    for col in EVIDENCE_COLUMNS:
        if col not in evidence.columns:
            evidence[col] = np.nan
        evidence[f"{col}_scaled"] = minmax_scale(evidence[col])

    scaled_cols = [f"{col}_scaled" for col in EVIDENCE_COLUMNS]
    evidence["integrated_evidence_score"] = evidence[scaled_cols].mean(axis=1)
    evidence = evidence.sort_values(
        ["integrated_evidence_score", "anti_malignant_shift_score", "rank"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    evidence.insert(0, "integrated_rank", np.arange(1, len(evidence) + 1))
    return evidence


def pivot_state_projection(state_summary: pd.DataFrame, tf_order: list[str]) -> pd.DataFrame:
    pivot = state_summary.pivot_table(
        index="tf",
        columns="celloracle_state",
        values="malignant_axis_projection_mean",
        aggfunc="mean",
    )
    ordered_states = [state for state in STATE_ORDER if state in pivot.columns]
    extra_states = [state for state in pivot.columns if state not in ordered_states]
    pivot = pivot.reindex(index=tf_order, columns=ordered_states + extra_states)
    return pivot


def select_top_tfs(evidence: pd.DataFrame, n: int) -> list[str]:
    return (
        evidence.sort_values(["integrated_evidence_score", "rank"], ascending=[False, True])
        .head(n)["tf"]
        .astype(str)
        .tolist()
    )


def plot_perturbation_ranking(evidence: pd.DataFrame, figure_dir: Path, top_n: int) -> list[str]:
    data = evidence.sort_values("anti_malignant_shift_score", ascending=True).tail(top_n)
    fig, ax = plt.subplots(figsize=(3.5, 3.8))
    colors = mpl.colormaps["viridis"](minmax_scale(data["integrated_evidence_score"]))
    ax.barh(data["tf"], data["anti_malignant_shift_score"], color=colors, edgecolor="0.2", linewidth=0.3)
    ax.axvline(0, color="0.25", linewidth=0.6)
    ax.set_xlabel("Anti-malignant shift score")
    ax.set_ylabel("")
    ax.set_title("CellOracle TF knockout ranking")
    sm = mpl.cm.ScalarMappable(cmap="viridis", norm=mpl.colors.Normalize(vmin=0, vmax=1))
    cbar = fig.colorbar(sm, ax=ax, pad=0.02, fraction=0.08)
    cbar.set_label("Integrated evidence")
    sns.despine(ax=ax)
    return save_pub_figure(fig, figure_dir / "celloracle_module6_9_perturbation_ranking")


def plot_evidence_heatmap(evidence: pd.DataFrame, figure_dir: Path, top_n: int) -> list[str]:
    rows = evidence.head(top_n).copy()
    scaled_cols = [f"{col}_scaled" for col in EVIDENCE_COLUMNS]
    labels = [
        "Anti-shift",
        "Expr. effect",
        "6.4 score",
        "Motif",
        "Fate",
        "Regulon r",
        "TF corr.",
        "GRN edges",
        "GRN coef.",
    ]
    matrix = rows.set_index("tf")[scaled_cols]
    matrix.columns = labels
    fig, ax = plt.subplots(figsize=(6.8, 0.30 * len(matrix) + 1.0))
    sns.heatmap(
        matrix,
        ax=ax,
        cmap="viridis",
        vmin=0,
        vmax=1,
        linewidths=0.3,
        linecolor="white",
        cbar_kws={"label": "Scaled evidence"},
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Candidate TF evidence matrix")
    ax.tick_params(axis="x", rotation=35, labelsize=5.5)
    ax.tick_params(axis="y", labelsize=6)
    return save_pub_figure(fig, figure_dir / "celloracle_module6_9_candidate_evidence_heatmap")


def plot_state_projection_heatmap(state_summary: pd.DataFrame, tf_order: list[str], figure_dir: Path) -> list[str]:
    pivot = pivot_state_projection(state_summary, tf_order=tf_order)
    pivot = pivot.rename(columns=STATE_LABELS)
    vmax = np.nanmax(np.abs(pivot.to_numpy()))
    fig, ax = plt.subplots(figsize=(4.8, 0.34 * len(pivot) + 1.1))
    sns.heatmap(
        pivot,
        ax=ax,
        cmap="RdBu_r",
        center=0,
        vmin=-vmax,
        vmax=vmax,
        linewidths=0.3,
        linecolor="white",
        cbar_kws={"label": "Projection to malignant axis"},
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("State-specific perturbation direction")
    ax.tick_params(axis="x", rotation=35)
    return save_pub_figure(fig, figure_dir / "celloracle_module6_9_state_projection_heatmap")


def _load_umap_background(h5ad_path: Path) -> pd.DataFrame:
    adata = ad.read_h5ad(h5ad_path)
    if "X_celloracle_umap" not in adata.obsm:
        raise ValueError("X_celloracle_umap not found in h5ad")
    coords = np.asarray(adata.obsm["X_celloracle_umap"])
    return pd.DataFrame(
        {
            "cell_id": adata.obs_names.astype(str),
            "umap_1": coords[:, 0],
            "umap_2": coords[:, 1],
            "celloracle_state": adata.obs["celloracle_state"].astype(str).to_numpy(),
        }
    )


def plot_top_tf_vector_fields(
    grid_arrows: pd.DataFrame,
    evidence: pd.DataFrame,
    h5ad_path: Path,
    figure_dir: Path,
    top_n: int,
) -> list[str] | None:
    top_tfs = select_top_tfs(evidence, top_n)
    background = _load_umap_background(h5ad_path)
    fig, axes = plt.subplots(1, len(top_tfs), figsize=(2.05 * len(top_tfs), 2.2), sharex=True, sharey=True)
    if len(top_tfs) == 1:
        axes = [axes]
    state_palette = {
        "normal_reference": "#B8B8B8",
        "stressed_injured": "#56B4E9",
        "regenerative_progenitor": "#009E73",
        "proliferating_candidate": "#E69F00",
        "malignant_or_malignant_like": "#D55E00",
    }
    for ax, tf in zip(axes, top_tfs):
        for state, color in state_palette.items():
            sub = background.loc[background["celloracle_state"] == state]
            ax.scatter(sub["umap_1"], sub["umap_2"], s=0.8, c=color, alpha=0.12, linewidths=0)
        arrows = grid_arrows.loc[grid_arrows["tf"].astype(str) == tf].copy()
        if not arrows.empty:
            threshold = arrows["flow_magnitude"].quantile(0.7)
            arrows = arrows.loc[arrows["flow_magnitude"] >= threshold]
            vectors = arrows[["flow_x", "flow_y"]].to_numpy(dtype=float)
            norms = np.linalg.norm(vectors, axis=1)
            valid = norms > 0
            arrows = arrows.loc[valid]
            vectors = vectors[valid] / norms[valid, None] * 1.25
            ax.quiver(
                arrows["grid_x"],
                arrows["grid_y"],
                vectors[:, 0],
                vectors[:, 1],
                angles="xy",
                scale_units="xy",
                scale=1,
                width=0.005,
                color="black",
                alpha=0.82,
            )
        ax.set_title(f"{tf} KO")
        ax.set_xlabel("UMAP 1")
        ax.set_aspect("equal")
    axes[0].set_ylabel("UMAP 2")
    for ax in axes[1:]:
        ax.set_ylabel("")
    fig.suptitle("CellOracle perturbation vector fields", y=1.03, fontsize=8)
    return save_pub_figure(fig, figure_dir / "celloracle_module6_9_top_tf_vector_fields")


def write_figure_contract(metadata_dir: Path, top_tfs: list[str]) -> Path:
    path = metadata_dir / "celloracle_module6_9_figure_contract.md"
    text = "\n".join(
        [
            "# Module 6.9 Figure Contract",
            "",
            "Core conclusion: CellOracle perturbation simulation prioritizes TFs whose knockout shifts malignant-like cells away from the malignant trajectory axis, with multi-layer support from CellRank, SCENIC/cisTarget, and state-specific GRNs.",
            "",
            "Archetype: quantitative grid with vector-field mechanistic inset.",
            "",
            "Evidence chain:",
            "- Ranking bar plot: magnitude of anti-malignant shift after TF knockout.",
            "- Evidence heatmap: concordance across perturbation, CellRank/SCENIC, base GRN, and state-specific GRN evidence.",
            "- State projection heatmap: state specificity of simulated shifts.",
            "- Vector field panels: spatial direction of top perturbations on the CellOracle UMAP.",
            "",
            f"Top TFs selected for vector-field panels: {', '.join(top_tfs)}.",
            "",
            "Review risks handled: raw metrics are retained in source tables; integrated score is used only for ordering; all figure panels are generated from saved Module 6.4, 6.7, and 6.8 outputs.",
        ]
    )
    path.write_text(text, encoding="utf-8")
    return path


def run_module6_9(
    metadata_dir: Path,
    figure_dir: Path,
    h5ad_path: Path,
    top_n: int,
    vector_top_n: int,
) -> dict:
    configure_matplotlib()
    metadata_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    ranking = pd.read_csv(metadata_dir / "celloracle_module6_8_perturbation_ranking.tsv", sep="\t")
    state_summary = pd.read_csv(metadata_dir / "celloracle_module6_8_state_shift_summary.tsv", sep="\t")
    grid_arrows = pd.read_csv(metadata_dir / "celloracle_module6_8_grid_arrows.tsv.gz", sep="\t")
    selection = pd.read_csv(metadata_dir / "celloracle_tf_selection.module6_4.tsv", sep="\t")
    grn_tf_summary = pd.read_csv(metadata_dir / "celloracle_module6_7_tf_network_summary.tsv", sep="\t")

    evidence = build_candidate_evidence_matrix(ranking, selection, grn_tf_summary)
    top_tfs = select_top_tfs(evidence, vector_top_n)

    evidence_path = metadata_dir / "celloracle_module6_9_candidate_evidence_matrix.tsv"
    top_source_path = metadata_dir / "celloracle_module6_9_top_tf_source_data.tsv"
    state_projection_path = metadata_dir / "celloracle_module6_9_state_projection_matrix.tsv"
    evidence.to_csv(evidence_path, sep="\t", index=False)
    evidence.head(top_n).to_csv(top_source_path, sep="\t", index=False)
    pivot_state_projection(state_summary, tf_order=evidence["tf"].astype(str).tolist()).to_csv(
        state_projection_path,
        sep="\t",
    )

    figure_paths = []
    figure_paths.extend(plot_perturbation_ranking(evidence, figure_dir, top_n=top_n))
    figure_paths.extend(plot_evidence_heatmap(evidence, figure_dir, top_n=top_n))
    figure_paths.extend(
        plot_state_projection_heatmap(
            state_summary,
            tf_order=evidence["tf"].astype(str).tolist(),
            figure_dir=figure_dir,
        )
    )
    vector_paths = plot_top_tf_vector_fields(grid_arrows, evidence, h5ad_path, figure_dir, top_n=vector_top_n)
    if vector_paths:
        figure_paths.extend(vector_paths)

    contract_path = write_figure_contract(metadata_dir, top_tfs)
    return {
        "candidate_evidence_matrix": str(evidence_path),
        "top_tf_source_data": str(top_source_path),
        "state_projection_matrix": str(state_projection_path),
        "figure_contract": str(contract_path),
        "figures": figure_paths,
        "top_anti_malignant_tfs": ranking.sort_values(
            "anti_malignant_shift_score",
            ascending=False,
        ).head(5)["tf"].astype(str).tolist(),
        "top_integrated_tfs": evidence.head(top_n)["tf"].astype(str).tolist(),
        "top_vector_tfs": top_tfs,
        "n_candidates": int(len(evidence)),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 6.9 visualize and integrate CellOracle perturbation evidence")
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--top-n", type=int, default=15)
    parser.add_argument("--vector-top-n", type=int, default=5)
    parser.add_argument("--report", type=Path, default=DEFAULT_METADATA_DIR / "celloracle_module6_9_report.json")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = datetime.now(timezone.utc)
    result = run_module6_9(
        metadata_dir=args.metadata_dir,
        figure_dir=args.figure_dir,
        h5ad_path=args.h5ad,
        top_n=args.top_n,
        vector_top_n=args.vector_top_n,
    )
    finished = datetime.now(timezone.utc)
    report = {
        "module": "6.9",
        "method": "CellOracle perturbation visualization and candidate TF evidence integration",
        "created_at_utc": finished.isoformat(),
        "elapsed_seconds": (finished - started).total_seconds(),
        "parameters": {
            "metadata_dir": str(args.metadata_dir),
            "figure_dir": str(args.figure_dir),
            "h5ad": str(args.h5ad),
            "top_n": args.top_n,
            "vector_top_n": args.vector_top_n,
        },
        "result": result,
        "figure_contract": {
            "backend": "Python/matplotlib",
            "archetype": "quantitative grid with vector-field mechanistic inset",
            "export_formats": ["svg", "pdf", "png", "tiff"],
        },
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
            "n_candidates": result["n_candidates"],
            "top_anti_malignant_tfs": result["top_anti_malignant_tfs"],
            "top_integrated_tfs": result["top_integrated_tfs"][:5],
            "top_vector_tfs": result["top_vector_tfs"],
            "n_figures": len(result["figures"]),
            "report": str(args.report),
        },
        ensure_ascii=False,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
