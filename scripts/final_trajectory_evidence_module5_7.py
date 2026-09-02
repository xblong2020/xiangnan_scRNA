from __future__ import annotations

import argparse
import json
import time
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

OKABE_ITO = {
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "bluish_green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "reddish_purple": "#CC79A7",
    "black": "#000000",
    "gray": "#7F7F7F",
}

EVIDENCE_TIER_COLORS = {
    "module3_high_conf_malignant": OKABE_ITO["vermillion"],
    "module3_cnv_supported_malignant": OKABE_ITO["orange"],
    "malignant_like_needs_review": OKABE_ITO["reddish_purple"],
    "copykat_aneuploid_without_module3_call": OKABE_ITO["yellow"],
    "cnv_proxy_aneuploid_without_module3_call": OKABE_ITO["sky_blue"],
    "no_cnv_evidence_or_reference": "#BDBDBD",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 5.7: final trajectory/CNV/malignant evidence panels and conclusion tables.")
    parser.add_argument(
        "--overlay",
        type=Path,
        default=ROOT / "metadata/trajectory/trajectory_module5_5_cnv_malignant_overlay_by_cell.tsv.gz",
    )
    parser.add_argument(
        "--bin-summary",
        type=Path,
        default=ROOT / "metadata/trajectory/trajectory_module5_5_pseudotime_bin_summary.tsv",
    )
    parser.add_argument(
        "--correlations",
        type=Path,
        default=ROOT / "metadata/trajectory/trajectory_module5_5_evidence_correlations.tsv",
    )
    parser.add_argument(
        "--robustness",
        type=Path,
        default=ROOT / "metadata/trajectory/trajectory_module5_6_robustness_summary.tsv",
    )
    parser.add_argument(
        "--batch-adjusted",
        type=Path,
        default=ROOT / "metadata/trajectory/trajectory_module5_6_batch_adjusted_correlations.tsv",
    )
    parser.add_argument(
        "--sensitivity-metrics",
        type=Path,
        default=ROOT / "metadata/trajectory/trajectory_module5_3_sensitivity_metrics.tsv",
    )
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/trajectory")
    parser.add_argument("--figures-dir", type=Path, default=ROOT / "figures/trajectory")
    parser.add_argument("--max-umap-points", type=int, default=5000)
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


def first_last_delta(df: pd.DataFrame, value_col: str, bin_col: str = "pseudotime_bin") -> dict[str, float]:
    if df.empty or value_col not in df.columns or bin_col not in df.columns:
        return {"early_value": np.nan, "late_value": np.nan, "delta": np.nan}
    work = df[[bin_col, value_col]].copy()
    work[bin_col] = pd.to_numeric(work[bin_col], errors="coerce")
    work[value_col] = pd.to_numeric(work[value_col], errors="coerce")
    work = work.replace([np.inf, -np.inf], np.nan).dropna()
    if work.empty:
        return {"early_value": np.nan, "late_value": np.nan, "delta": np.nan}
    work = work.sort_values(bin_col)
    early = float(work[value_col].iloc[0])
    late = float(work[value_col].iloc[-1])
    return {"early_value": early, "late_value": late, "delta": float(late - early)}


def minmax_scale(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    finite = numeric.replace([np.inf, -np.inf], np.nan).dropna()
    if finite.empty:
        return numeric * np.nan
    low = finite.min()
    high = finite.max()
    if high == low:
        return numeric * 0.0
    return (numeric - low) / (high - low)


def lookup_robustness(
    robustness: pd.DataFrame,
    run_id: str,
    method: str,
    group_type: str,
    feature: str,
    column: str,
    default: object = np.nan,
) -> object:
    if robustness.empty:
        return default
    sub = robustness.loc[
        robustness["run_id"].astype(str).eq(str(run_id))
        & robustness["method"].astype(str).eq(str(method))
        & robustness["group_type"].astype(str).eq(str(group_type))
        & robustness["feature"].astype(str).eq(str(feature))
    ]
    if sub.empty or column not in sub.columns:
        return default
    return sub[column].iloc[0]


def lookup_correlation(correlations: pd.DataFrame, run_id: str, method: str, feature: str) -> float:
    if correlations.empty:
        return np.nan
    sub = correlations.loc[
        correlations["run_id"].astype(str).eq(str(run_id))
        & correlations["method"].astype(str).eq(str(method))
        & correlations["feature"].astype(str).eq(str(feature))
    ]
    if sub.empty or "spearman_rho" not in sub.columns:
        return np.nan
    return float(sub["spearman_rho"].iloc[0])


def lookup_centered_correlation(adjusted: pd.DataFrame, run_id: str, method: str, batch_type: str, feature: str) -> float:
    if adjusted.empty:
        return np.nan
    sub = adjusted.loc[
        adjusted["run_id"].astype(str).eq(str(run_id))
        & adjusted["method"].astype(str).eq(str(method))
        & adjusted["batch_type"].astype(str).eq(str(batch_type))
        & adjusted["feature"].astype(str).eq(str(feature))
    ]
    if sub.empty or "centered_spearman_rho" not in sub.columns:
        return np.nan
    return float(sub["centered_spearman_rho"].iloc[0])


def classify_final_evidence(row: dict[str, object] | pd.Series) -> str:
    cnv_delta = float(row.get("cnv_supported_delta", np.nan))
    hcc_delta = float(row.get("hcc_malignant_module_delta", np.nan))
    proliferation_delta = float(row.get("proliferation_module_delta", np.nan))
    sample_label = str(row.get("sample_robustness_label", ""))
    dataset_label = str(row.get("dataset_robustness_label", ""))

    cnv_positive = np.isfinite(cnv_delta) and cnv_delta >= 0.2
    module_positive = (np.isfinite(hcc_delta) and hcc_delta > 0.1) or (np.isfinite(proliferation_delta) and proliferation_delta > 0.1)
    if cnv_positive and module_positive and sample_label == "robust_positive" and dataset_label == "robust_positive":
        return "consensus_supported"
    if cnv_positive and module_positive and dataset_label == "robust_positive":
        return "supported_with_sample_composition_caveat"
    if cnv_positive and module_positive:
        return "supported_with_group_sensitivity"
    if module_positive:
        return "malignant_module_trend_without_cnv_consensus"
    if cnv_positive:
        return "cnv_trend_without_marker_consensus"
    return "insufficient_or_mixed"


def build_conclusion_table(
    bin_summary: pd.DataFrame,
    robustness: pd.DataFrame,
    correlations: pd.DataFrame,
    adjusted: pd.DataFrame | None = None,
) -> pd.DataFrame:
    adjusted = adjusted if adjusted is not None else pd.DataFrame()
    rows = []
    for (run_id, method), sub in bin_summary.groupby(["run_id", "method"], observed=True, sort=True):
        cnv = first_last_delta(sub, "cnv_supported_fraction")
        review = first_last_delta(sub, "review_fraction")
        copykat = first_last_delta(sub, "copykat_aneuploid_fraction")
        hcc = first_last_delta(sub, "mean_hcc_malignant_module")
        proliferation = first_last_delta(sub, "mean_proliferation_module")
        cnv_burden = first_last_delta(sub, "mean_cnv_proxy_burden")

        sample_label = lookup_robustness(robustness, str(run_id), str(method), "sample_id", "module3_cnv_supported", "robustness_label", "")
        sample_pos = lookup_robustness(
            robustness, str(run_id), str(method), "sample_id", "module3_cnv_supported", "positive_group_fraction", np.nan
        )
        sample_min_loo = lookup_robustness(robustness, str(run_id), str(method), "sample_id", "module3_cnv_supported", "min_loo_delta", np.nan)
        dataset_label = lookup_robustness(
            robustness, str(run_id), str(method), "dataset", "HCC_Malignant_Associated", "robustness_label", ""
        )
        dataset_pos = lookup_robustness(
            robustness, str(run_id), str(method), "dataset", "HCC_Malignant_Associated", "positive_group_fraction", np.nan
        )
        largest_sample_fraction = lookup_robustness(
            robustness, str(run_id), str(method), "sample_id", "module3_cnv_supported", "largest_group_fraction", np.nan
        )

        row = {
            "run_id": run_id,
            "method": method,
            "cnv_supported_early": cnv["early_value"],
            "cnv_supported_late": cnv["late_value"],
            "cnv_supported_delta": cnv["delta"],
            "copykat_aneuploid_delta": copykat["delta"],
            "review_fraction_delta": review["delta"],
            "hcc_malignant_module_early": hcc["early_value"],
            "hcc_malignant_module_late": hcc["late_value"],
            "hcc_malignant_module_delta": hcc["delta"],
            "proliferation_module_delta": proliferation["delta"],
            "cnv_proxy_burden_delta": cnv_burden["delta"],
            "sample_robustness_label": sample_label,
            "sample_positive_group_fraction": sample_pos,
            "sample_min_leave_one_out_delta": sample_min_loo,
            "dataset_robustness_label": dataset_label,
            "dataset_positive_group_fraction": dataset_pos,
            "largest_sample_fraction": largest_sample_fraction,
            "module3_cnv_supported_spearman_rho": lookup_correlation(
                correlations, str(run_id), str(method), "module3_cnv_supported"
            ),
            "hcc_malignant_module_spearman_rho": lookup_correlation(
                correlations, str(run_id), str(method), "HCC_Malignant_Associated"
            ),
            "proliferation_spearman_rho": lookup_correlation(correlations, str(run_id), str(method), "Proliferation"),
            "sample_centered_cnv_spearman_rho": lookup_centered_correlation(
                adjusted, str(run_id), str(method), "sample_id", "module3_cnv_supported"
            ),
            "sample_centered_hcc_spearman_rho": lookup_centered_correlation(
                adjusted, str(run_id), str(method), "sample_id", "HCC_Malignant_Associated"
            ),
        }
        row["final_evidence_label"] = classify_final_evidence(row)
        rows.append(row)
    return pd.DataFrame(rows)


def summarize_consensus(conclusion: pd.DataFrame) -> pd.DataFrame:
    if conclusion.empty:
        return pd.DataFrame()
    rows = []
    metrics = [
        "cnv_supported_delta",
        "hcc_malignant_module_delta",
        "proliferation_module_delta",
        "module3_cnv_supported_spearman_rho",
        "sample_centered_cnv_spearman_rho",
    ]
    for run_id, sub in conclusion.groupby("run_id", observed=True, sort=True):
        row = {
            "run_id": run_id,
            "n_methods": int(sub.shape[0]),
            "n_consensus_supported": int(sub["final_evidence_label"].eq("consensus_supported").sum()),
            "n_supported_with_sample_caveat": int(sub["final_evidence_label"].eq("supported_with_sample_composition_caveat").sum()),
        }
        for metric in metrics:
            if metric in sub.columns:
                row[f"median_{metric}"] = float(pd.to_numeric(sub[metric], errors="coerce").median())
                row[f"min_{metric}"] = float(pd.to_numeric(sub[metric], errors="coerce").min())
                row[f"max_{metric}"] = float(pd.to_numeric(sub[metric], errors="coerce").max())
        rows.append(row)
    return pd.DataFrame(rows)


def make_markdown_table(conclusion: pd.DataFrame, path: Path) -> None:
    columns = [
        "run_id",
        "method",
        "final_evidence_label",
        "cnv_supported_delta",
        "hcc_malignant_module_delta",
        "proliferation_module_delta",
        "sample_robustness_label",
        "dataset_robustness_label",
        "sample_positive_group_fraction",
        "sample_centered_cnv_spearman_rho",
    ]
    out = conclusion[columns].copy()
    for column in out.columns:
        if pd.api.types.is_numeric_dtype(out[column]):
            out[column] = out[column].map(lambda value: "" if pd.isna(value) else f"{value:.3f}")
    headers = [str(column) for column in out.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for _, row in out.iterrows():
        values = [str(row[column]).replace("|", "\\|") for column in out.columns]
        lines.append("| " + " | ".join(values) + " |")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.12, 1.08, label, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top", ha="left")


def plot_umap_pseudotime(ax: plt.Axes, data: pd.DataFrame) -> None:
    scatter = ax.scatter(
        data["global_umap_1"],
        data["global_umap_2"],
        c=pd.to_numeric(data["pseudotime_norm"], errors="coerce"),
        s=4,
        cmap="viridis",
        linewidths=0,
        rasterized=True,
    )
    ax.set_title("Trajectory pseudotime")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.02, label="Pseudotime")


def plot_umap_evidence(ax: plt.Axes, data: pd.DataFrame) -> None:
    tiers = [tier for tier in EVIDENCE_TIER_COLORS if tier in set(data["cnv_evidence_tier"].astype(str))]
    for tier in tiers:
        sub = data.loc[data["cnv_evidence_tier"].astype(str).eq(tier)]
        ax.scatter(
            sub["global_umap_1"],
            sub["global_umap_2"],
            s=4,
            color=EVIDENCE_TIER_COLORS[tier],
            label=tier.replace("_", " "),
            linewidths=0,
            rasterized=True,
        )
    ax.set_title("CNV / malignant evidence tier")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.legend(frameon=False, fontsize=5, markerscale=2, loc="best")


def plot_trend_panel(ax: plt.Axes, bins: pd.DataFrame) -> None:
    bins = bins.sort_values("pseudotime_bin")
    x = bins["mean_pseudotime"]
    lines = [
        ("cnv_supported_fraction", "CNV-supported fraction", OKABE_ITO["vermillion"], "-"),
        ("review_fraction", "Review fraction", OKABE_ITO["reddish_purple"], "-"),
        ("mean_hcc_malignant_module", "HCC module, scaled", OKABE_ITO["blue"], "--"),
        ("mean_proliferation_module", "Proliferation, scaled", OKABE_ITO["orange"], "--"),
    ]
    for column, label, color, linestyle in lines:
        if column not in bins.columns:
            continue
        y = pd.to_numeric(bins[column], errors="coerce")
        if "mean_" in column:
            y = minmax_scale(y)
        ax.plot(x, y, marker="o", markersize=3, linewidth=1.4, linestyle=linestyle, color=color, label=label)
    ax.set_xlabel("Normalized pseudotime")
    ax.set_ylabel("Fraction or scaled mean")
    ax.set_ylim(-0.04, 1.04)
    ax.set_title("Evidence along trajectory")
    ax.legend(frameon=False, fontsize=6)
    ax.grid(alpha=0.2)


def plot_conclusion_panel(ax: plt.Axes, row: pd.Series) -> None:
    metrics = pd.Series(
        {
            "CNV fraction": row.get("cnv_supported_delta", np.nan),
            "HCC module": row.get("hcc_malignant_module_delta", np.nan),
            "Proliferation": row.get("proliferation_module_delta", np.nan),
            "sample centered CNV rho": row.get("sample_centered_cnv_spearman_rho", np.nan),
        }
    )
    colors = [OKABE_ITO["vermillion"], OKABE_ITO["blue"], OKABE_ITO["orange"], OKABE_ITO["bluish_green"]]
    ax.barh(metrics.index, metrics.values, color=colors, alpha=0.85)
    ax.axvline(0, color="black", linewidth=0.8)
    evidence_label = str(row.get("final_evidence_label", "")).replace("_", " ")
    ax.set_title("Final evidence summary\n" + textwrap.fill(evidence_label, width=36), fontsize=8)
    ax.set_xlabel("Delta or Spearman rho")
    ax.grid(axis="x", alpha=0.2)


def plot_final_panel(
    overlay: pd.DataFrame,
    bin_summary: pd.DataFrame,
    conclusion_row: pd.Series,
    run_id: str,
    method: str,
    figures_dir: Path,
    max_umap_points: int,
) -> list[str]:
    data = overlay.loc[overlay["run_id"].eq(run_id) & overlay["method"].eq(method)].copy()
    bins = bin_summary.loc[bin_summary["run_id"].eq(run_id) & bin_summary["method"].eq(method)].copy()
    if data.empty or bins.empty:
        return []
    if data.shape[0] > max_umap_points:
        data = data.sample(max_umap_points, random_state=57)
    figures_dir.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(7.2, 6.4))
    grid = fig.add_gridspec(2, 2, hspace=0.38, wspace=0.35)
    axes = [
        fig.add_subplot(grid[0, 0]),
        fig.add_subplot(grid[0, 1]),
        fig.add_subplot(grid[1, 0]),
        fig.add_subplot(grid[1, 1]),
    ]
    plot_umap_pseudotime(axes[0], data)
    plot_umap_evidence(axes[1], data)
    plot_trend_panel(axes[2], bins)
    plot_conclusion_panel(axes[3], conclusion_row)
    for label, ax in zip(["A", "B", "C", "D"], axes):
        panel_label(ax, label)
    fig.suptitle(f"Trajectory, CNV, and malignant evidence: {run_id} / {method}", fontsize=11, y=0.99)
    png_path = figures_dir / f"trajectory_module5_7_final_panel__{run_id}__{method}.png"
    pdf_path = figures_dir / f"trajectory_module5_7_final_panel__{run_id}__{method}.pdf"
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)
    return [str(png_path.resolve()), str(pdf_path.resolve())]


def plot_consensus_heatmap(conclusion: pd.DataFrame, figures_dir: Path) -> list[str]:
    if conclusion.empty:
        return []
    metrics = [
        "cnv_supported_delta",
        "hcc_malignant_module_delta",
        "proliferation_module_delta",
        "module3_cnv_supported_spearman_rho",
        "sample_centered_cnv_spearman_rho",
    ]
    labels = [f"{row.run_id}\n{row.method}" for row in conclusion.itertuples(index=False)]
    values = conclusion[metrics].apply(pd.to_numeric, errors="coerce").T
    scaled = values.apply(minmax_scale, axis=1)
    figures_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.6, 3.2))
    image = ax.imshow(scaled, aspect="auto", cmap="cividis", vmin=0, vmax=1)
    ax.set_yticks(range(len(metrics)))
    ax.set_yticklabels([metric.replace("_", " ") for metric in metrics])
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=35, ha="right")
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            value = values.iloc[i, j]
            text = "" if pd.isna(value) else f"{value:.2f}"
            ax.text(j, i, text, ha="center", va="center", fontsize=6, color="white" if scaled.iloc[i, j] > 0.55 else "black")
    ax.set_title("Consensus evidence across trajectory methods")
    plt.colorbar(image, ax=ax, fraction=0.025, pad=0.02, label="Row-scaled value")
    fig.tight_layout()
    png_path = figures_dir / "trajectory_module5_7_consensus_heatmap.png"
    pdf_path = figures_dir / "trajectory_module5_7_consensus_heatmap.pdf"
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)
    return [str(png_path.resolve()), str(pdf_path.resolve())]


def main() -> int:
    args = parse_args()
    start = time.time()
    configure_plot_style()
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    args.figures_dir.mkdir(parents=True, exist_ok=True)

    overlay = pd.read_csv(args.overlay, sep="\t")
    bin_summary = pd.read_csv(args.bin_summary, sep="\t")
    correlations = pd.read_csv(args.correlations, sep="\t")
    robustness = pd.read_csv(args.robustness, sep="\t")
    adjusted = pd.read_csv(args.batch_adjusted, sep="\t") if args.batch_adjusted.exists() else pd.DataFrame()

    conclusion = build_conclusion_table(bin_summary, robustness, correlations, adjusted)
    consensus = summarize_consensus(conclusion)

    conclusion_path = args.metadata_dir / "trajectory_module5_7_conclusion_table.tsv"
    conclusion_md_path = args.metadata_dir / "trajectory_module5_7_conclusion_table.md"
    consensus_path = args.metadata_dir / "trajectory_module5_7_method_consensus.tsv"
    report_path = args.metadata_dir / "trajectory_module5_7_report.json"

    conclusion.to_csv(conclusion_path, sep="\t", index=False)
    consensus.to_csv(consensus_path, sep="\t", index=False)
    make_markdown_table(conclusion, conclusion_md_path)

    figure_paths = []
    for row in conclusion.itertuples(index=False):
        run_id = str(row.run_id)
        method = str(row.method)
        print(f"FINAL_PANEL {run_id} {method}", flush=True)
        conclusion_row = conclusion.loc[conclusion["run_id"].eq(run_id) & conclusion["method"].eq(method)].iloc[0]
        figure_paths.extend(plot_final_panel(overlay, bin_summary, conclusion_row, run_id, method, args.figures_dir, args.max_umap_points))
    figure_paths.extend(plot_consensus_heatmap(conclusion, args.figures_dir))

    report = {
        "module": "5.7",
        "method": "integrated final trajectory/CNV/malignant evidence figure panels and conclusion tables",
        "n_methods": int(conclusion.shape[0]),
        "final_evidence_label_counts": conclusion["final_evidence_label"].value_counts().to_dict(),
        "outputs": {
            "conclusion_table": str(conclusion_path.resolve()),
            "conclusion_markdown": str(conclusion_md_path.resolve()),
            "method_consensus": str(consensus_path.resolve()),
            "figures": figure_paths,
        },
        "elapsed_seconds": round(time.time() - start, 3),
    }
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
