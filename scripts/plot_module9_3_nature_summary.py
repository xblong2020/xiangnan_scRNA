from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colors, patches


ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = ROOT / "metadata/driver"
FIGURE_DIR = ROOT / "figures/driver"

INDIRECT_PATH = METADATA_DIR / "module9_3_indirect_effects.tsv"
COEFFICIENT_PATH = METADATA_DIR / "module9_3_path_model_coefficients.tsv"
GRADE_PATH = METADATA_DIR / "module9_3_evidence_grade.tsv"
MODULE9_1_GRADE = METADATA_DIR / "module9_1_evidence_grade.tsv"
MODULE9_2_GRADE = METADATA_DIR / "module9_2_evidence_grade.tsv"
SOURCE_DATA_PATH = METADATA_DIR / "module9_3_nature_summary_source_data.tsv"

PALETTE = {
    "blue": "#0F4D92",
    "blue_soft": "#B4C0E4",
    "red": "#B64342",
    "red_soft": "#F6CFCB",
    "gold": "#D6A21E",
    "green": "#2E9E44",
    "green_soft": "#DDF3DE",
    "neutral_light": "#D8D8D8",
    "neutral_mid": "#767676",
    "neutral_dark": "#272727",
}


def apply_style() -> None:
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["font.size"] = 7
    plt.rcParams["axes.spines.right"] = False
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.linewidth"] = 0.75
    plt.rcParams["legend.frameon"] = False


def add_panel_label(ax: plt.Axes, label: str, x: float = -0.08, y: float = 1.04) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")


def save_figure(fig: plt.Figure, base: Path) -> list[str]:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    saved = []
    for ext, kwargs in {
        "svg": {},
        "pdf": {},
        "tiff": {"dpi": 600},
        "png": {"dpi": 300},
    }.items():
        path = base.with_suffix(f".{ext}")
        fig.savefig(path, bbox_inches="tight", **kwargs)
        saved.append(str(path))
    plt.close(fig)
    return saved


def compact_dataset_label(value: str) -> str:
    text = str(value)
    if text.startswith("scRNA_pseudobulk|"):
        parts = text.split("|")
        run = "main" if parts[1] == "main_strict" else "sens"
        method = parts[2].replace("slingshot_hepatocyte_pca", "slingshot PCA").replace("slingshot_scanvi", "slingshot scANVI")
        return f"scRNA {run} {method}"
    return text.replace("bulk_pooled", "bulk pooled")


def compact_outcome(value: str) -> str:
    return (
        str(value)
        .replace("outcome_mean_malignant_fate", "mean malignant fate")
        .replace("outcome_malignant_fraction", "malignant fraction")
        .replace("outcome_sample_stage", "sample stage")
        .replace("stage", "stage")
        .replace("grade", "grade")
        .replace("time", "OS")
    )


def compact_analysis_label(dataset_id: str, outcome: str, include_method: bool = False) -> str:
    dataset = str(dataset_id)
    out = compact_outcome(outcome)
    if dataset.startswith("scRNA_pseudobulk|"):
        parts = dataset.split("|")
        run = "main" if parts[1] == "main_strict" else "sens"
        if include_method:
            method = parts[2].replace("slingshot_hepatocyte_pca", "S-PCA").replace("slingshot_scanvi", "S-scANVI").replace("monocle3", "M3")
            return f"scRNA {run} {method}: {out}"
        return f"scRNA {run}: {out}"
    return f"{dataset.replace('bulk_pooled', 'bulk')}: {out}"


def load_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.Series, pd.DataFrame]:
    indirect = pd.read_csv(INDIRECT_PATH, sep="\t")
    coefficients = pd.read_csv(COEFFICIENT_PATH, sep="\t")
    grade = pd.read_csv(GRADE_PATH, sep="\t").iloc[0]
    evidence_rows = []
    if MODULE9_1_GRADE.exists():
        g1 = pd.read_csv(MODULE9_1_GRADE, sep="\t").iloc[0]
        evidence_rows.append({"module": "9.1 temporal", "label": g1.get("final_support_label", ""), "score": 0.0})
    if MODULE9_2_GRADE.exists():
        g2 = pd.read_csv(MODULE9_2_GRADE, sep="\t").iloc[0]
        label = str(g2.get("network_direction_label", ""))
        evidence_rows.append({"module": "9.2 network", "label": label, "score": 0.5 if "partial" in label else 1.0 if "supported" in label else 0.0})
    label = str(grade.get("final_support_label", ""))
    evidence_rows.append({"module": "9.3 mediation", "label": label, "score": 0.5 if "partial" in label else 1.0 if "supported" in label else 0.0})
    evidence = pd.DataFrame(evidence_rows)

    source = indirect.loc[indirect["indirect_path"].eq("A_to_B_to_C_to_outcome")].copy()
    source["compact_dataset"] = source["dataset_id"].map(compact_dataset_label)
    source["compact_outcome"] = source["outcome"].map(compact_outcome)
    source.to_csv(SOURCE_DATA_PATH, sep="\t", index=False)
    return indirect, coefficients, grade, evidence


def draw_path_schematic(ax: plt.Axes, grade: pd.Series) -> None:
    ax.set_axis_off()
    add_panel_label(ax, "a", x=-0.02)
    ax.text(0.02, 0.94, "Observed-variable path analysis", fontsize=9, fontweight="bold", transform=ax.transAxes)
    ax.text(
        0.02,
        0.86,
        f"{grade['final_support_label'].replace('_', ' ')}; {int(grade['n_supported_sequential_indirect_effects'])}/"
        f"{int(grade['n_tested_sequential_indirect_effects'])} sequential effects supported",
        fontsize=6.5,
        color=PALETTE["neutral_mid"],
        transform=ax.transAxes,
    )
    nodes = [
        (0.13, 0.54, "A loss\nHNF4A/PPARA", PALETTE["blue_soft"]),
        (0.44, 0.54, "B activation\nAP-1/CEBPB/EGR1", PALETTE["green_soft"]),
        (0.73, 0.54, "C axis\nSOX4", PALETTE["red_soft"]),
        (0.87, 0.23, "Outcome\nOS/stage/burden", "#F2E6D9"),
    ]
    for x, y, text, color in nodes:
        box = patches.FancyBboxPatch(
            (x - 0.11, y - 0.08),
            0.20,
            0.16,
            boxstyle="round,pad=0.012,rounding_size=0.02",
            facecolor=color,
            edgecolor=PALETTE["neutral_dark"],
            linewidth=0.8,
            transform=ax.transAxes,
        )
        ax.add_patch(box)
        ax.text(x, y, text, ha="center", va="center", fontsize=5.6, transform=ax.transAxes)

    def arr(x1, y1, x2, y2, color=PALETTE["blue"], lw=1.4, rad=0.0):
        ax.add_patch(
            patches.FancyArrowPatch(
                (x1, y1),
                (x2, y2),
                arrowstyle="-|>",
                mutation_scale=9,
                linewidth=lw,
                color=color,
                connectionstyle=f"arc3,rad={rad}",
                transform=ax.transAxes,
            )
        )

    arr(0.23, 0.52, 0.32, 0.52)
    arr(0.54, 0.52, 0.62, 0.52)
    arr(0.79, 0.45, 0.86, 0.30, PALETTE["red"], 1.4)
    arr(0.18, 0.43, 0.84, 0.24, PALETTE["gold"], 1.0, -0.2)
    ax.text(0.06, 0.08, "Bootstrap 95% CI; BH-FDR across tested indirect effects", fontsize=5.8, transform=ax.transAxes)


def draw_forest(ax: plt.Axes, indirect: pd.DataFrame, title: str, max_rows: int = 12) -> None:
    seq = indirect.loc[indirect["indirect_path"].eq("A_to_B_to_C_to_outcome")].copy()
    seq = seq.loc[seq["status"].eq("tested")].copy()
    seq["analysis_label"] = [compact_analysis_label(r.dataset_id, r.outcome, include_method=True) for r in seq.itertuples()]
    seq["abs_effect"] = seq["effect"].abs()
    seq = (
        seq.sort_values(["p.adjust.global", "abs_effect"], ascending=[True, False])
        .drop_duplicates("analysis_label")
        .head(max_rows)
    )
    y = np.arange(len(seq))[::-1]
    colors_ = [PALETTE["blue"] if v > 0 else PALETTE["red"] for v in seq["effect"]]
    for yi, (_, row), color in zip(y, seq.iterrows(), colors_):
        ax.plot([row["ci_low"], row["ci_high"]], [yi, yi], color=color, linewidth=1.2)
        ax.plot(row["effect"], yi, marker="o", color=color, markersize=4)
    ax.axvline(0, color=PALETTE["neutral_mid"], linestyle="--", linewidth=0.8)
    labels = [r.analysis_label for r in seq.itertuples()]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=5.6)
    ax.set_xlabel("Sequential indirect effect")
    ax.set_title(title, fontsize=8, loc="left")


def draw_coefficient_heatmap(ax: plt.Axes, coefficients: pd.DataFrame) -> None:
    add_panel_label(ax, "b")
    paths = ["A_to_B", "B_to_C", "C_to_outcome"]
    use = coefficients.loc[coefficients["path"].isin(paths) & coefficients["status"].eq("tested")].copy()
    use["analysis"] = [compact_analysis_label(r.dataset_id, r.outcome, include_method=False) for r in use.itertuples()]
    pivot = use.pivot_table(index="analysis", columns="path", values="coef", aggfunc="mean")
    if pivot.empty:
        ax.set_axis_off()
        ax.text(0.5, 0.5, "No tested coefficients", ha="center", va="center")
        return
    preferred = [
        "TCGA-LIHC: stage",
        "TCGA-LIHC: grade",
        "TCGA-LIHC: OS",
        "ICGC-LIRI-JP: stage",
        "ICGC-LIRI-JP: OS",
        "bulk: stage",
        "bulk: grade",
        "bulk: OS",
        "scRNA main: mean malignant fate",
        "scRNA main: malignant fraction",
        "scRNA main: sample stage",
        "scRNA sens: mean malignant fate",
        "scRNA sens: malignant fraction",
        "scRNA sens: sample stage",
    ]
    pivot = pivot.reindex([idx for idx in preferred if idx in pivot.index]).reindex(columns=paths)
    vmax = np.nanmax(np.abs(pivot.to_numpy()))
    vmax = vmax if np.isfinite(vmax) and vmax > 0 else 1.0
    im = ax.imshow(pivot.to_numpy(), cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(np.arange(len(paths)))
    ax.set_xticklabels(["A->B", "B->C", "C->Y"], rotation=30, ha="right")
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=5.8)
    ax.set_title("Path coefficient consistency", fontsize=8, loc="left")
    ax.figure.colorbar(im, ax=ax, fraction=0.045, pad=0.02, label="coefficient")


def draw_evidence_matrix(ax: plt.Axes, evidence: pd.DataFrame) -> None:
    add_panel_label(ax, "d")
    matrix = evidence[["score"]].to_numpy()
    cmap = colors.ListedColormap([PALETTE["neutral_light"], PALETTE["gold"], PALETTE["blue"]])
    norm = colors.BoundaryNorm([-0.1, 0.25, 0.75, 1.1], cmap.N)
    im = ax.imshow(matrix, cmap=cmap, norm=norm, aspect="auto")
    ax.set_xticks([0])
    ax.set_xticklabels(["support"])
    ax.set_yticks(np.arange(len(evidence)))
    ax.set_yticklabels(evidence["module"])
    for i, row in evidence.reset_index(drop=True).iterrows():
        label = str(row["label"]).replace("_", " ")
        ax.text(0, i, label, ha="center", va="center", fontsize=5.5)
    ax.set_title("Integrated Module 9 evidence", fontsize=8, loc="left")
    cbar = ax.figure.colorbar(im, ax=ax, fraction=0.05, pad=0.03, ticks=[0, 0.5, 1])
    cbar.ax.set_yticklabels(["no", "partial", "full"])


def make_composite(indirect: pd.DataFrame, coefficients: pd.DataFrame, grade: pd.Series, evidence: pd.DataFrame) -> list[str]:
    fig = plt.figure(figsize=(8.4, 5.8))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.18], height_ratios=[0.95, 1.20], wspace=0.64, hspace=0.62)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    draw_path_schematic(ax_a, grade)
    draw_coefficient_heatmap(ax_b, coefficients)
    add_panel_label(ax_c, "c")
    draw_forest(ax_c, indirect, "Sequential indirect effects", max_rows=10)
    draw_evidence_matrix(ax_d, evidence)
    return save_figure(fig, FIGURE_DIR / "module9_3_mediation_path_summary")


def make_forest(indirect: pd.DataFrame) -> list[str]:
    fig, ax = plt.subplots(figsize=(5.8, 5.2))
    draw_forest(ax, indirect, "Module 9.3 sequential indirect effect forest", max_rows=14)
    return save_figure(fig, FIGURE_DIR / "module9_3_indirect_effect_forest")


def make_heatmap(coefficients: pd.DataFrame) -> list[str]:
    fig, ax = plt.subplots(figsize=(5.0, 5.6))
    draw_coefficient_heatmap(ax, coefficients)
    return save_figure(fig, FIGURE_DIR / "module9_3_dataset_consistency_heatmap")


def make_integrated(evidence: pd.DataFrame) -> list[str]:
    fig, ax = plt.subplots(figsize=(4.2, 2.0))
    draw_evidence_matrix(ax, evidence)
    return save_figure(fig, FIGURE_DIR / "module9_3_integrated_module9_evidence_summary")


def main() -> None:
    apply_style()
    indirect, coefficients, grade, evidence = load_data()
    saved = []
    saved.extend(make_composite(indirect, coefficients, grade, evidence))
    saved.extend(make_forest(indirect))
    saved.extend(make_heatmap(coefficients))
    saved.extend(make_integrated(evidence))
    print("\n".join(saved + [str(SOURCE_DATA_PATH)]))


if __name__ == "__main__":
    main()
