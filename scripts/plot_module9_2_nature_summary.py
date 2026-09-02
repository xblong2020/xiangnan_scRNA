from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import patches


ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = ROOT / "metadata/driver"
FIGURE_DIR = ROOT / "figures/driver"

ASYMMETRY_PATH = METADATA_DIR / "module9_2_asymmetry_tests.tsv"
GRADE_PATH = METADATA_DIR / "module9_2_evidence_grade.tsv"
SOURCE_DATA_PATH = METADATA_DIR / "module9_2_nature_summary_source_data.tsv"
OUT_BASE = FIGURE_DIR / "module9_2_nature_network_direction_summary"

PALETTE = {
    "neutral_dark": "#272727",
    "neutral_mid": "#767676",
    "neutral_light": "#D8D8D8",
    "blue": "#0F4D92",
    "blue_soft": "#B4C0E4",
    "teal": "#42949E",
    "green": "#2E9E44",
    "green_soft": "#DDF3DE",
    "red": "#B64342",
    "red_soft": "#F6CFCB",
    "gold": "#D6A21E",
}


COMPARISON_LABELS = {
    "A_to_B_vs_C_to_A": "A to B vs C to A",
    "A_to_C_vs_C_to_A": "A to C vs C to A",
    "B_to_C_vs_C_to_B": "B to C vs C to B",
    "SOX4_self_vs_reverse_upstream": "SOX4 self vs reverse upstream",
}

COMPARISON_ORDER = list(COMPARISON_LABELS)


def apply_style() -> None:
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["font.size"] = 7
    plt.rcParams["axes.spines.right"] = False
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.linewidth"] = 0.7
    plt.rcParams["legend.frameon"] = False


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.08, 1.05, label, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")


def load_data() -> tuple[pd.DataFrame, pd.Series]:
    asym = pd.read_csv(ASYMMETRY_PATH, sep="\t")
    grade = pd.read_csv(GRADE_PATH, sep="\t").iloc[0]
    source = asym.loc[
        asym["subset"].isin(["consensus", "joint_consensus"]),
        [
            "comparison",
            "subset",
            "evidence_source",
            "forward_impact",
            "reverse_impact",
            "control_impact",
            "directionality_index",
            "forward_reverse_ratio",
            "forward_vs_control_ratio",
            "support_label",
        ],
    ].copy()
    source["comparison_label"] = source["comparison"].map(COMPARISON_LABELS)
    source.to_csv(SOURCE_DATA_PATH, sep="\t", index=False)
    return source, grade


def support_to_score(label: str) -> float:
    label = str(label)
    if label in {"consensus_supported", "joint_consensus_supported"}:
        return 1.0
    if label in {"partial_joint_consensus_support"}:
        return 0.5
    return 0.0


def draw_schematic(ax: plt.Axes, joint: pd.DataFrame, grade: pd.Series) -> None:
    ax.set_axis_off()
    add_panel_label(ax, "a")
    ax.text(0.0, 1.02, "KO-only network direction model", fontsize=8, fontweight="bold", transform=ax.transAxes)
    ax.text(
        0.0,
        0.94,
        str(grade["network_direction_label"]).replace("_", " "),
        fontsize=6.5,
        color=PALETTE["neutral_mid"],
        transform=ax.transAxes,
    )

    nodes = {
        "A": (0.15, 0.52, "A upstream\nHNF4A/PPARA"),
        "B": (0.51, 0.52, "B transition\nAP-1/CEBPB/EGR1"),
        "C": (0.86, 0.52, "C fate axis\nSOX4\nmalignant-like"),
    }
    node_colors = {"A": PALETTE["blue_soft"], "B": PALETTE["green_soft"], "C": PALETTE["red_soft"]}
    edge_lookup = joint.set_index("comparison")

    for key, (x, y, text) in nodes.items():
        box = patches.FancyBboxPatch(
            (x - 0.115, y - 0.095),
            0.23,
            0.18,
            boxstyle="round,pad=0.015,rounding_size=0.025",
            facecolor=node_colors[key],
            edgecolor=PALETTE["neutral_dark"],
            linewidth=0.8,
            transform=ax.transAxes,
        )
        ax.add_patch(box)
        ax.text(x, y, text, ha="center", va="center", fontsize=5.9, transform=ax.transAxes)

    def arrow(start: str, end: str, comp: str, rad: float = 0.0, yoff: float = 0.0, label_xy: tuple[float, float] | None = None) -> None:
        sx, sy, _ = nodes[start]
        ex, ey, _ = nodes[end]
        row = edge_lookup.loc[comp]
        idx = float(row["directionality_index"])
        label = "supported" if row["support_label"] == "joint_consensus_supported" else "partial"
        color = PALETTE["blue"] if label == "supported" else PALETTE["gold"]
        arr = patches.FancyArrowPatch(
            (sx + 0.14, sy + yoff),
            (ex - 0.14, ey + yoff),
            arrowstyle="-|>",
            mutation_scale=10,
            linewidth=1.6 if label == "supported" else 1.1,
            color=color,
            connectionstyle=f"arc3,rad={rad}",
            transform=ax.transAxes,
        )
        ax.add_patch(arr)
        lx, ly = label_xy if label_xy is not None else ((sx + ex) / 2, sy + yoff + 0.08 + abs(rad) * 0.04)
        ax.text(
            lx,
            ly,
            f"{label}\nDI={idx:.2f}",
            ha="center",
            va="center",
            fontsize=5.5,
            color=color,
            transform=ax.transAxes,
        )

    arrow("A", "B", "A_to_B_vs_C_to_A", label_xy=(0.33, 0.72))
    arrow("B", "C", "B_to_C_vs_C_to_B", label_xy=(0.69, 0.72))
    arrow("A", "C", "A_to_C_vs_C_to_A", rad=-0.20, yoff=-0.08, label_xy=(0.50, 0.31))

    reverse = edge_lookup.loc["SOX4_self_vs_reverse_upstream"]
    ax.text(
        0.35,
        0.17,
        f"SOX4 self impact exceeds reverse-upstream impact (DI={float(reverse['directionality_index']):.2f}).",
        fontsize=6.2,
        color=PALETTE["neutral_dark"],
        transform=ax.transAxes,
    )
    ax.text(
        0.02,
        0.05,
        "Restore/fixed-mode outputs unavailable; interpretation is KO-only.",
        fontsize=5.8,
        color=PALETTE["neutral_mid"],
        transform=ax.transAxes,
    )


def draw_directionality(ax: plt.Axes, joint: pd.DataFrame) -> None:
    add_panel_label(ax, "b")
    ordered = joint.set_index("comparison").loc[COMPARISON_ORDER].reset_index()
    y = np.arange(len(ordered))[::-1]
    vals = ordered["directionality_index"].astype(float).to_numpy()
    colors = [PALETTE["blue"] if lab == "joint_consensus_supported" else PALETTE["gold"] for lab in ordered["support_label"]]
    ax.barh(y, vals, color=colors, edgecolor="black", linewidth=0.4, height=0.62)
    ax.axvline(0, color=PALETTE["neutral_mid"], linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([COMPARISON_LABELS[c] for c in ordered["comparison"]])
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("Directionality index")
    ax.set_title("Forward impact exceeds reverse impact", fontsize=8, loc="left")
    for yi, val, ratio in zip(y, vals, ordered["forward_reverse_ratio"].astype(float)):
        ax.text(val + 0.02, yi, f"{val:.2f} ({ratio:.1f}x)", va="center", fontsize=6)


def draw_method_heatmap(ax: plt.Axes, source: pd.DataFrame) -> None:
    add_panel_label(ax, "c")
    methods = ["sctenifold_signature", "celloracle_abs_delta", "combined"]
    method_labels = ["scTenifoldKnk", "CellOracle", "Joint"]
    matrix = np.zeros((len(COMPARISON_ORDER), len(methods)))
    label_matrix = [["" for _ in methods] for _ in COMPARISON_ORDER]
    for i, comp in enumerate(COMPARISON_ORDER):
        for j, method in enumerate(methods):
            if method == "combined":
                row = source.loc[source["comparison"].eq(comp) & source["subset"].eq("joint_consensus")]
            else:
                row = source.loc[
                    source["comparison"].eq(comp)
                    & source["subset"].eq("consensus")
                    & source["evidence_source"].eq(method)
                ]
            if row.empty:
                matrix[i, j] = np.nan
                label_matrix[i][j] = "NA"
            else:
                lab = str(row.iloc[0]["support_label"])
                matrix[i, j] = support_to_score(lab)
                label_matrix[i][j] = "full" if matrix[i, j] == 1 else ("part" if matrix[i, j] == 0.5 else "no")

    cmap = matplotlib.colors.ListedColormap([PALETTE["neutral_light"], PALETTE["gold"], PALETTE["blue"]])
    norm = matplotlib.colors.BoundaryNorm([-0.1, 0.25, 0.75, 1.1], cmap.N)
    im = ax.imshow(matrix, cmap=cmap, norm=norm, aspect="auto")
    ax.set_xticks(np.arange(len(method_labels)))
    ax.set_xticklabels(method_labels, rotation=25, ha="right")
    ax.set_yticks(np.arange(len(COMPARISON_ORDER)))
    ax.set_yticklabels([COMPARISON_LABELS[c] for c in COMPARISON_ORDER])
    ax.set_title("Evidence support by method", fontsize=8, loc="left")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, label_matrix[i][j], ha="center", va="center", fontsize=6, color=PALETTE["neutral_dark"])
    cbar = ax.figure.colorbar(im, ax=ax, fraction=0.045, pad=0.03, ticks=[0, 0.5, 1])
    cbar.ax.set_yticklabels(["no", "partial", "full"])
    cbar.set_label("support", fontsize=6)


def draw_forward_reverse(ax: plt.Axes, joint: pd.DataFrame) -> None:
    add_panel_label(ax, "d")
    ordered = joint.set_index("comparison").loc[COMPARISON_ORDER].reset_index()
    x = np.arange(len(ordered))
    width = 0.36
    forward = ordered["forward_impact"].astype(float).to_numpy()
    reverse = ordered["reverse_impact"].astype(float).to_numpy()
    ax.bar(x - width / 2, forward, width, label="forward", color=PALETTE["blue_soft"], edgecolor="black", linewidth=0.4)
    ax.bar(x + width / 2, reverse, width, label="reverse", color=PALETTE["red_soft"], edgecolor="black", linewidth=0.4)
    ax.set_xticks(x)
    ax.set_xticklabels(["A-B", "A-C", "B-C", "SOX4"], rotation=0)
    ax.set_ylabel("Mean impact")
    ax.set_title("Joint forward vs reverse effect size", fontsize=8, loc="left")
    ax.legend(loc="upper left", fontsize=6, ncol=2)
    ymax = max(float(np.nanmax(forward)), float(np.nanmax(reverse))) * 1.2
    ax.set_ylim(0, ymax)


def main() -> None:
    apply_style()
    source, grade = load_data()
    joint = source.loc[source["subset"].eq("joint_consensus")].copy()

    fig = plt.figure(figsize=(7.2, 5.2), constrained_layout=False)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.18, 1.0], height_ratios=[1.0, 1.0], wspace=0.58, hspace=0.68)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    draw_schematic(ax_a, joint, grade)
    draw_directionality(ax_b, joint)
    draw_method_heatmap(ax_c, source)
    draw_forward_reverse(ax_d, joint)

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{OUT_BASE}.svg", bbox_inches="tight")
    fig.savefig(f"{OUT_BASE}.pdf", bbox_inches="tight")
    fig.savefig(f"{OUT_BASE}.tiff", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT_BASE}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(
        "\n".join(
            [
                str(OUT_BASE.with_suffix(".svg")),
                str(OUT_BASE.with_suffix(".pdf")),
                str(OUT_BASE.with_suffix(".tiff")),
                str(OUT_BASE.with_suffix(".png")),
                str(SOURCE_DATA_PATH),
            ]
        )
    )


if __name__ == "__main__":
    main()
