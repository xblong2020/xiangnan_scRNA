from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon


ROOT = Path(__file__).resolve().parents[1]

# ggsci::pal_lancet() colors, reused in the Python plotting workflow.
LANCET_COLORS = ["#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F", "#FDAF91"]


def calculate_rss(auc: pd.DataFrame, labels: pd.Series, min_cells_per_state: int) -> pd.DataFrame:
    labels = labels.reindex(auc.index).astype(str)
    state_counts = labels.value_counts()
    valid_states = state_counts[state_counts >= min_cells_per_state].index.tolist()
    if len(valid_states) < 2:
        raise ValueError("At least two states must have the requested minimum number of cells.")

    use = labels.isin(valid_states)
    auc = auc.loc[use].apply(pd.to_numeric, errors="coerce").fillna(0.0).clip(lower=0.0)
    labels = labels.loc[use]
    rows: list[dict[str, object]] = []
    for regulon in auc.columns:
        activity = auc[regulon].to_numpy(dtype=float)
        total_activity = activity.sum()
        if total_activity <= 0:
            continue
        activity_distribution = activity / total_activity
        for state in valid_states:
            is_target = labels.eq(state).to_numpy(dtype=float)
            target_distribution = is_target / is_target.sum()
            rss = 1.0 - float(jensenshannon(activity_distribution, target_distribution, base=2))
            target_values = activity[is_target.astype(bool)]
            other_values = activity[~is_target.astype(bool)]
            rows.append(
                {
                    "state": state,
                    "regulon": str(regulon),
                    "rss": rss,
                    "n_target_cells": int(is_target.sum()),
                    "mean_auc_target": float(target_values.mean()),
                    "mean_auc_other": float(other_values.mean()),
                    "mean_auc_delta": float(target_values.mean() - other_values.mean()),
                }
            )
    return pd.DataFrame(rows).sort_values(["state", "rss"], ascending=[True, False]).reset_index(drop=True)


def select_top_regulons(rss: pd.DataFrame, target_state: str, top_n: int) -> pd.DataFrame:
    selected = rss.loc[rss["state"].astype(str).eq(str(target_state))].copy()
    if selected.empty:
        available = sorted(rss["state"].astype(str).unique())
        raise ValueError(f"Target state '{target_state}' is unavailable. Available states: {available}")
    return selected.sort_values(["rss", "mean_auc_delta"], ascending=[False, False]).head(top_n).reset_index(drop=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Reproduce Figure 1H as an RSS ranking of target-state-specific cisTarget regulons.")
    parser.add_argument(
        "--auc",
        type=Path,
        default=ROOT / "metadata/driver/driver_module6_3c_cistarget_regulon_auc.tsv.gz",
    )
    parser.add_argument(
        "--h5ad",
        type=Path,
        default=ROOT / "data/processed/driver/driver_cistarget_regulon_activity.module6_3c.h5ad",
    )
    parser.add_argument("--state-column", default="hepatocyte_state_label")
    parser.add_argument("--target-state", default="malignant_hepatocyte_candidate_needs_cnv")
    parser.add_argument("--min-cells-per-state", type=int, default=20)
    parser.add_argument("--top-n", type=int, default=5)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "figures/figure1")
    parser.add_argument("--output-prefix", default="figure1H_malignant_hepatocyte_rss_ranking")
    parser.add_argument("--width", type=float, default=6.1)
    parser.add_argument("--height", type=float, default=4.5)
    parser.add_argument("--dpi", type=int, default=600)
    return parser.parse_args()


def read_labels(path: Path, state_column: str) -> pd.Series:
    adata = ad.read_h5ad(path, backed="r")
    try:
        if state_column not in adata.obs.columns:
            raise KeyError(f"Missing state column: {state_column}")
        labels = adata.obs[state_column].astype(str).copy()
        labels.index = adata.obs_names.astype(str)
        return labels
    finally:
        adata.file.close()


def plot_rss_ranking(top: pd.DataFrame, target_state: str, width: float, height: float, dpi: int, output_base: Path) -> None:
    display = top.iloc[::-1].copy()
    colors = [LANCET_COLORS[i % len(LANCET_COLORS)] for i in range(display.shape[0])]
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"], "font.size": 9})
    fig, ax = plt.subplots(figsize=(width, height))
    bars = ax.barh(display["regulon"], display["rss"], color=colors, height=0.67)
    for bar, row in zip(bars, display.itertuples(index=False)):
        ax.text(
            bar.get_width() + 0.006,
            bar.get_y() + bar.get_height() / 2,
            f"RSS {row.rss:.3f}",
            ha="left",
            va="center",
            fontsize=8,
        )
    ax.set_xlabel("Regulon specificity score (RSS)")
    ax.set_ylabel("")
    ax.set_xlim(0, max(float(display["rss"].max()) * 1.28, 0.10))
    ax.set_title("CNV-malignant hepatocyte candidate state", loc="left", fontsize=9, color="#555555", pad=12)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color="#E5E5E5", linewidth=0.6)
    ax.set_axisbelow(True)
    fig.suptitle("Figure 1H. State-specific regulon ranking", x=0.125, y=0.99, ha="left", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(output_base.with_suffix(".png"), dpi=dpi, facecolor="white", bbox_inches="tight")
    fig.savefig(output_base.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    if args.min_cells_per_state <= 0 or args.top_n <= 0:
        raise ValueError("min-cells-per-state and top-n must be positive")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    auc = pd.read_csv(args.auc, sep="\t").set_index("cell_id")
    auc.index = auc.index.astype(str)
    labels = read_labels(args.h5ad, args.state_column)
    shared_cells = auc.index.intersection(labels.index)
    if shared_cells.empty:
        raise ValueError("No overlapping cell IDs between AUCell and state annotations.")
    auc = auc.loc[shared_cells]
    labels = labels.loc[shared_cells]
    rss = calculate_rss(auc, labels, args.min_cells_per_state)
    top = select_top_regulons(rss, args.target_state, args.top_n)

    base = args.output_dir / args.output_prefix
    plot_rss_ranking(top, args.target_state, args.width, args.height, args.dpi, base)
    rss.to_csv(base.with_name(base.name + "_all_states.tsv.gz"), sep="\t", index=False, compression="gzip")
    top.to_csv(base.with_name(base.name + "_top_regulons.tsv"), sep="\t", index=False)
    valid_counts = labels.value_counts()
    base.with_name(base.name + "_report.json").write_text(
        json.dumps(
            {
                "method": "SCENIC RSS using Jensen-Shannon distance on cisTarget-pruned AUCell regulon activities",
                "auc": str(args.auc),
                "h5ad": str(args.h5ad),
                "state_column": args.state_column,
                "target_state": args.target_state,
                "min_cells_per_state": args.min_cells_per_state,
                "top_n": args.top_n,
                "n_shared_cells": int(len(shared_cells)),
                "state_counts": {str(k): int(v) for k, v in valid_counts.items()},
                "valid_states": sorted(rss["state"].unique().tolist()),
                "top_regulons": top[["regulon", "rss"]].to_dict(orient="records"),
                "palette": "ggsci Lancet",
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"WROTE {base.with_suffix('.png')}")
    print(f"WROTE {base.with_suffix('.pdf')}")
    print(f"WROTE {base.with_name(base.name + '_top_regulons.tsv')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
