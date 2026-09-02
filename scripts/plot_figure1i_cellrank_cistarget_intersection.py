from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

# ggsci::pal_lancet() colors, reused here in the Python plotting workflow.
LANCET_COLORS = {
    "cellrank": "#00468B",
    "cistarget": "#42B540",
    "shared": "#ED0000",
}


def compute_intersection(
    cellrank: pd.DataFrame,
    regulons: pd.DataFrame,
    known_tfs: set[str],
    top_n: int,
    qvalue_cutoff: float,
    min_corr: float,
    min_nes: float,
    min_motifs: int,
) -> dict[str, object]:
    required_cellrank = {"gene", "corr", "qval", "rank_positive_corr"}
    required_regulons = {"tf", "regulon", "best_nes", "n_targets", "n_motifs"}
    if missing := required_cellrank.difference(cellrank.columns):
        raise ValueError(f"CellRank table is missing columns: {sorted(missing)}")
    if missing := required_regulons.difference(regulons.columns):
        raise ValueError(f"cisTarget table is missing columns: {sorted(missing)}")

    known_tfs = {str(tf).upper() for tf in known_tfs}
    drivers = cellrank.copy()
    drivers["gene"] = drivers["gene"].astype(str).str.upper()
    drivers["corr"] = pd.to_numeric(drivers["corr"], errors="coerce")
    drivers["qval"] = pd.to_numeric(drivers["qval"], errors="coerce")
    drivers["rank_positive_corr"] = pd.to_numeric(drivers["rank_positive_corr"], errors="coerce")
    drivers = drivers.loc[
        drivers["gene"].isin(known_tfs)
        & drivers["rank_positive_corr"].le(top_n)
        & drivers["corr"].gt(min_corr)
        & drivers["qval"].le(qvalue_cutoff),
        ["gene", "corr", "qval", "rank_positive_corr"],
    ].drop_duplicates("gene")

    scenic = regulons.copy()
    scenic["tf"] = scenic["tf"].astype(str).str.upper()
    scenic["best_nes"] = pd.to_numeric(scenic["best_nes"], errors="coerce")
    scenic["n_targets"] = pd.to_numeric(scenic["n_targets"], errors="coerce")
    scenic["n_motifs"] = pd.to_numeric(scenic["n_motifs"], errors="coerce")
    scenic = scenic.loc[
        scenic["tf"].isin(known_tfs)
        & scenic["best_nes"].ge(min_nes)
        & scenic["n_motifs"].ge(min_motifs),
        ["tf", "regulon", "best_nes", "n_targets", "n_motifs"],
    ].drop_duplicates("tf")

    candidates = drivers.merge(scenic, left_on="gene", right_on="tf", how="inner")
    candidates = candidates.sort_values(["corr", "best_nes"], ascending=[False, False]).reset_index(drop=True)
    return {
        "cellrank_tf_set": set(drivers["gene"]),
        "cistarget_tf_set": set(scenic["tf"]),
        "overlap_set": set(candidates["tf"]),
        "candidates": candidates,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reproduce Figure 1I by intersecting CellRank CNV-fate TF drivers with cisTarget regulon TFs."
    )
    parser.add_argument(
        "--cellrank-drivers",
        type=Path,
        default=ROOT / "metadata/driver/driver_module6_2_cellrank_lineage_drivers.tsv.gz",
    )
    parser.add_argument(
        "--cistarget-regulons",
        type=Path,
        default=ROOT / "metadata/driver/driver_module6_3c_cistarget_regulon_summary.tsv",
    )
    parser.add_argument(
        "--tf-list",
        type=Path,
        default=ROOT / "metadata/driver/scenic_resources/allTFs_hg38.txt",
    )
    parser.add_argument("--lineage", default="cnv_supported_malignant")
    parser.add_argument("--cellrank-top-n", type=int, default=500)
    parser.add_argument("--qvalue-cutoff", type=float, default=0.05)
    parser.add_argument("--min-corr", type=float, default=0.0)
    parser.add_argument("--min-nes", type=float, default=3.0)
    parser.add_argument("--min-motifs", type=int, default=3)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "figures/figure1")
    parser.add_argument("--output-prefix", default="figure1I_cellrank_cistarget_tf_intersection")
    parser.add_argument("--width", type=float, default=9.6)
    parser.add_argument("--height", type=float, default=4.6)
    parser.add_argument("--dpi", type=int, default=600)
    return parser.parse_args()


def draw_intersection_panel(ax: plt.Axes, cellrank_set: set[str], cistarget_set: set[str], overlap: set[str]) -> None:
    ax.set_aspect("equal")
    ax.axis("off")
    left_color = LANCET_COLORS["cellrank"]
    right_color = LANCET_COLORS["cistarget"]
    ax.add_patch(Circle((0.42, 0.50), 0.31, facecolor=left_color, edgecolor=left_color, alpha=0.38, linewidth=1.2))
    ax.add_patch(Circle((0.60, 0.50), 0.31, facecolor=right_color, edgecolor=right_color, alpha=0.38, linewidth=1.2))
    ax.text(0.29, 0.83, "CellRank\nCNV-fate TF drivers", ha="center", va="center", fontsize=8, fontweight="bold")
    ax.text(0.73, 0.83, "cisTarget-pruned\nregulon TFs", ha="center", va="center", fontsize=8, fontweight="bold")
    ax.text(0.28, 0.50, str(len(cellrank_set.difference(cistarget_set))), ha="center", va="center", fontsize=14, fontweight="bold", color=left_color)
    ax.text(0.74, 0.50, str(len(cistarget_set.difference(cellrank_set))), ha="center", va="center", fontsize=14, fontweight="bold", color=right_color)
    ax.text(0.51, 0.54, str(len(overlap)), ha="center", va="center", fontsize=17, fontweight="bold", color="#222222")
    overlap_label = "\n".join(sorted(overlap)) if overlap else "No overlap"
    ax.text(0.51, 0.36, overlap_label, ha="center", va="top", fontsize=9, fontweight="bold", color="#222222")
    ax.text(0.51, 0.08, "Shared TF candidates", ha="center", va="center", fontsize=8, color="#555555")
    ax.set_xlim(0.02, 0.98)
    ax.set_ylim(0.0, 1.0)


def draw_rank_panel(ax: plt.Axes, candidates: pd.DataFrame) -> None:
    if candidates.empty:
        ax.text(0.5, 0.5, "No shared TF candidates", ha="center", va="center")
        ax.axis("off")
        return
    display = candidates.iloc[::-1].copy()
    bars = ax.barh(display["tf"], display["corr"], color=LANCET_COLORS["shared"], height=0.65)
    ax.set_xlabel("CellRank correlation with CNV-malignant fate")
    ax.set_ylabel("")
    ax.set_title("Shared TF candidates", loc="left", fontweight="bold", fontsize=10)
    for bar, row in zip(bars, display.itertuples(index=False)):
        ax.text(
            bar.get_width() + 0.006,
            bar.get_y() + bar.get_height() / 2,
            f"NES {row.best_nes:.2f}",
            va="center",
            ha="left",
            fontsize=7,
        )
    ax.set_xlim(0, max(display["corr"].max() * 1.34, 0.05))
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color="#E5E5E5", linewidth=0.6)
    ax.set_axisbelow(True)


def main() -> int:
    args = parse_args()
    if args.cellrank_top_n <= 0:
        raise ValueError("cellrank-top-n must be positive")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    cellrank = pd.read_csv(args.cellrank_drivers, sep="\t")
    cellrank = cellrank.loc[cellrank["lineage"].astype(str).eq(args.lineage)].copy()
    regulons = pd.read_csv(args.cistarget_regulons, sep="\t")
    known_tfs = {line.strip() for line in args.tf_list.read_text(encoding="utf-8").splitlines() if line.strip()}
    result = compute_intersection(
        cellrank,
        regulons,
        known_tfs=known_tfs,
        top_n=args.cellrank_top_n,
        qvalue_cutoff=args.qvalue_cutoff,
        min_corr=args.min_corr,
        min_nes=args.min_nes,
        min_motifs=args.min_motifs,
    )
    candidates = result["candidates"]

    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"], "font.size": 8})
    fig, axes = plt.subplots(1, 2, figsize=(args.width, args.height), gridspec_kw={"width_ratios": [1.08, 1.0]})
    draw_intersection_panel(axes[0], result["cellrank_tf_set"], result["cistarget_tf_set"], result["overlap_set"])
    draw_rank_panel(axes[1], candidates)
    fig.suptitle("Figure 1I. Key TFs jointly supported by CellRank and cisTarget", x=0.06, y=0.98, ha="left", fontsize=12, fontweight="bold")
    fig.text(
        0.06,
        0.02,
        f"CellRank: CNV-malignant fate, positive driver rank <= {args.cellrank_top_n}, q <= {args.qvalue_cutoff}; "
        f"cisTarget: motif-pruned regulon, NES >= {args.min_nes}, motifs >= {args.min_motifs}.",
        ha="left",
        va="bottom",
        fontsize=7,
        color="#555555",
    )
    fig.tight_layout(rect=(0, 0.06, 1, 0.93))

    base = args.output_dir / args.output_prefix
    fig.savefig(base.with_suffix(".png"), dpi=args.dpi, facecolor="white", bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    plt.close(fig)
    candidates.to_csv(base.with_name(base.name + "_candidates.tsv"), sep="\t", index=False)
    base.with_name(base.name + "_report.json").write_text(
        json.dumps(
            {
                "cellrank_drivers": str(args.cellrank_drivers),
                "cistarget_regulons": str(args.cistarget_regulons),
                "tf_list": str(args.tf_list),
                "lineage": args.lineage,
                "cellrank_top_n": args.cellrank_top_n,
                "qvalue_cutoff": args.qvalue_cutoff,
                "min_corr": args.min_corr,
                "min_nes": args.min_nes,
                "min_motifs": args.min_motifs,
                "n_cellrank_tf_drivers": len(result["cellrank_tf_set"]),
                "n_cistarget_regulon_tfs": len(result["cistarget_tf_set"]),
                "n_shared_candidates": len(result["overlap_set"]),
                "shared_candidates": sorted(result["overlap_set"]),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"WROTE {base.with_suffix('.png')}")
    print(f"WROTE {base.with_suffix('.pdf')}")
    print(f"WROTE {base.with_name(base.name + '_candidates.tsv')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
