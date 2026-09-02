from __future__ import annotations

import argparse
import json
import time
from collections import OrderedDict
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]

STATE_PANELS: OrderedDict[str, list[str]] = OrderedDict(
    [
        (
            "Mature_Hepatocyte",
            ["ALB", "APOA1", "APOA2", "TTR", "HPD", "ASGR1", "CYP3A4", "CYP2E1", "CYP2C9", "HNF4A", "CPS1", "ASS1"],
        ),
        (
            "Stressed_Injured",
            ["HSPA1A", "HSPA1B", "HSP90AA1", "DNAJB1", "FOS", "JUN", "JUNB", "ATF3", "DDIT3", "SAA1", "SAA2", "MT1G", "MT2A", "IER3"],
        ),
        ("Regenerative_Progenitor", ["KRT19", "EPCAM", "SOX9", "KRT7", "TACSTD2", "CD24", "PROM1", "ANXA4"]),
        ("HCC_Malignant_Associated", ["AFP", "GPC3", "SPP1", "MDK", "IGF2BP3", "MUC1", "CEACAM5"]),
        ("Proliferation", ["MKI67", "TOP2A", "STMN1", "TYMS", "UBE2C", "PCNA", "MCM2", "MCM5", "HMGB2"]),
        ("Cholangiocyte", ["KRT19", "KRT7", "EPCAM", "SOX9", "TACSTD2", "ANXA4", "CLDN4", "MUC1"]),
        ("Immune", ["PTPRC", "LST1", "C1QA", "C1QB", "CD68", "LYZ", "CD3D", "NKG7", "MS4A1", "JCHAIN"]),
        ("Endothelial", ["PECAM1", "VWF", "KDR", "EMCN", "RAMP2", "ESAM"]),
        ("Stromal_HSC_Pericyte", ["COL1A1", "COL1A2", "COL3A1", "DCN", "ACTA2", "RGS5", "PDGFRB", "TAGLN"]),
    ]
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 5.4: gene and module trends along trajectory pseudotime.")
    parser.add_argument(
        "--input-h5ad",
        type=Path,
        default=ROOT / "data/processed/trajectory/trajectory_hepatocyte_cnv_scanvi.stage_root_end.module5_2.h5ad",
    )
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/trajectory")
    parser.add_argument("--figures-dir", type=Path, default=ROOT / "figures/trajectory")
    parser.add_argument("--n-bins", type=int, default=10)
    parser.add_argument("--top-genes-per-method", type=int, default=30)
    parser.add_argument("--max-genes", type=int, default=0, help="0 means all genes in the h5ad.")
    return parser.parse_args()


def assign_pseudotime_bins(values: pd.Series, n_bins: int = 10) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    out = pd.Series(pd.NA, index=values.index, dtype="Int64")
    finite = numeric.replace([np.inf, -np.inf], np.nan).dropna()
    if finite.empty:
        return out
    n_effective = min(n_bins, int(finite.nunique()))
    if n_effective <= 1:
        out.loc[finite.index] = 0
        return out
    ranked = finite.rank(method="first")
    labels = list(range(n_effective))
    out.loc[finite.index] = pd.qcut(ranked, q=n_effective, labels=labels, duplicates="drop").astype("Int64")
    return out


def compute_module_scores(
    expr: pd.DataFrame,
    panels: OrderedDict[str, list[str]] | dict[str, list[str]],
) -> tuple[pd.DataFrame, dict[str, dict[str, object]]]:
    scores = pd.DataFrame(index=expr.index)
    availability: dict[str, dict[str, object]] = {}
    genes = set(expr.columns.astype(str))
    for module, module_genes in panels.items():
        available = [gene for gene in module_genes if gene in genes]
        availability[module] = {
            "module": module,
            "n_requested": len(module_genes),
            "n_available": len(available),
            "genes_available": ";".join(available),
            "genes_missing": ";".join([gene for gene in module_genes if gene not in genes]),
        }
        if available:
            scores[module] = expr[available].mean(axis=1)
        else:
            scores[module] = np.nan
    return scores, availability


def summarize_trend(
    values: pd.Series,
    pseudotime: pd.Series,
    n_bins: int = 10,
    bins: pd.Series | None = None,
) -> dict[str, object]:
    values = pd.to_numeric(values, errors="coerce")
    pseudotime = pd.to_numeric(pseudotime, errors="coerce")
    finite_mask = values.notna() & pseudotime.notna() & np.isfinite(values) & np.isfinite(pseudotime)
    if finite_mask.sum() < 3 or values.loc[finite_mask].nunique() < 2 or pseudotime.loc[finite_mask].nunique() < 2:
        return {
            "n_cells": int(finite_mask.sum()),
            "spearman_rho": np.nan,
            "spearman_pvalue": np.nan,
            "delta_last_first_bin": np.nan,
            "trend_direction": "insufficient",
        }
    rho, pvalue = spearmanr(values.loc[finite_mask], pseudotime.loc[finite_mask])
    use_bins = bins if bins is not None else assign_pseudotime_bins(pseudotime, n_bins=n_bins)
    bin_df = pd.DataFrame({"bin": use_bins, "value": values, "pseudotime": pseudotime}).loc[finite_mask].dropna()
    bin_means = bin_df.groupby("bin", observed=True)["value"].mean()
    if bin_means.empty:
        delta = np.nan
    else:
        delta = float(bin_means.iloc[-1] - bin_means.iloc[0])
    if np.isfinite(rho) and np.isfinite(delta) and rho >= 0.2 and delta > 0:
        direction = "increasing"
    elif np.isfinite(rho) and np.isfinite(delta) and rho <= -0.2 and delta < 0:
        direction = "decreasing"
    else:
        direction = "flat_or_mixed"
    return {
        "n_cells": int(finite_mask.sum()),
        "spearman_rho": float(rho),
        "spearman_pvalue": float(pvalue),
        "delta_last_first_bin": delta,
        "trend_direction": direction,
    }


def pseudotime_inputs(metadata_dir: Path) -> list[dict[str, object]]:
    return [
        {
            "run_id": "main_strict",
            "path": metadata_dir / "trajectory_module5_3_main_strict_pseudotime_merged.tsv.gz",
            "methods": {
                "monocle3": "main_strict__monocle3_norm",
                "slingshot_scanvi": "main_strict__slingshot_scanvi_norm",
                "slingshot_hepatocyte_pca": "main_strict__slingshot_hepatocyte_pca_norm",
            },
        },
        {
            "run_id": "sensitivity_include_review",
            "path": metadata_dir / "trajectory_module5_3_sensitivity_include_review_pseudotime_merged.tsv.gz",
            "methods": {
                "monocle3": "sensitivity_include_review__monocle3_norm",
                "slingshot_scanvi": "sensitivity_include_review__slingshot_scanvi_norm",
                "slingshot_hepatocyte_pca": "sensitivity_include_review__slingshot_hepatocyte_pca_norm",
            },
        },
    ]


def read_expression(input_h5ad: Path, cell_ids: pd.Index) -> pd.DataFrame:
    backed = ad.read_h5ad(input_h5ad, backed="r")
    obs_index = pd.Index(backed.obs_names.astype(str))
    positions = obs_index.get_indexer(cell_ids.astype(str))
    if (positions < 0).any():
        missing = cell_ids[positions < 0]
        backed.file.close()
        raise ValueError(f"{len(missing)} pseudotime cells are missing from {input_h5ad}")
    sub = backed[positions, :].to_memory()
    backed.file.close()
    x = sub.X
    if sparse.issparse(x):
        x = x.toarray()
    return pd.DataFrame(np.asarray(x, dtype=np.float32), index=sub.obs_names.astype(str), columns=sub.var_names.astype(str))


def bin_means_for_values(
    values: pd.Series,
    pseudotime: pd.Series,
    bins: pd.Series,
    value_name: str,
    run_id: str,
    method: str,
) -> pd.DataFrame:
    df = pd.DataFrame({"pseudotime": pseudotime, "bin": bins, "value": values}).dropna()
    if df.empty:
        return pd.DataFrame()
    out = (
        df.groupby("bin", observed=True)
        .agg(n_cells=("value", "size"), mean_pseudotime=("pseudotime", "mean"), mean_value=("value", "mean"), sem_value=("value", "sem"))
        .reset_index()
    )
    out.insert(0, "feature", value_name)
    out.insert(0, "method", method)
    out.insert(0, "run_id", run_id)
    return out


def plot_module_trends(module_bins: pd.DataFrame, run_id: str, method: str, figures_dir: Path) -> str | None:
    sub = module_bins.loc[(module_bins["run_id"].eq(run_id)) & (module_bins["method"].eq(method))].copy()
    if sub.empty:
        return None
    figures_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 5))
    for module, data in sub.groupby("feature", observed=True):
        data = data.sort_values("bin")
        ax.plot(data["mean_pseudotime"], data["mean_value"], marker="o", linewidth=1.8, markersize=3, label=module)
    ax.set_xlabel("Normalized pseudotime")
    ax.set_ylabel("Mean log-normalized module score")
    ax.set_title(f"Module trends: {run_id} / {method}")
    ax.legend(fontsize=7, ncol=2, frameon=False)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    path = figures_dir / f"trajectory_module5_4_module_trends__{run_id}__{method}.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return str(path.resolve())


def main() -> int:
    args = parse_args()
    start = time.time()
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    args.figures_dir.mkdir(parents=True, exist_ok=True)

    module_trend_rows = []
    module_bin_rows = []
    gene_trend_rows = []
    gene_bin_rows = []
    module_score_rows = []
    availability_rows = []
    figure_paths = []
    run_cache: dict[str, dict[str, object]] = {}

    for run in pseudotime_inputs(args.metadata_dir):
        run_id = str(run["run_id"])
        print(f"READ_PSEUDOTIME {run_id} {run['path']}", flush=True)
        pt = pd.read_csv(run["path"], sep="\t")
        pt = pt.drop_duplicates("cell_id").set_index("cell_id", drop=False)
        expr = read_expression(args.input_h5ad, pd.Index(pt.index.astype(str)))
        if args.max_genes and args.max_genes > 0:
            expr = expr.iloc[:, : args.max_genes]
        run_cache[run_id] = {"pt": pt, "expr": expr, "methods": run["methods"]}
        scores, availability = compute_module_scores(expr, STATE_PANELS)
        for module, row in availability.items():
            availability_rows.append({"run_id": run_id, **row})
        score_out = scores.copy()
        score_out.insert(0, "cell_id", score_out.index)
        score_out.insert(0, "run_id", run_id)
        module_score_rows.append(score_out)

        gene_to_modules: dict[str, list[str]] = {}
        for module, genes in STATE_PANELS.items():
            for gene in genes:
                gene_to_modules.setdefault(gene, []).append(module)

        for method, pt_col in run["methods"].items():
            print(f"TRENDS {run_id} {method}", flush=True)
            pseudotime = pd.to_numeric(pt[pt_col], errors="coerce")
            bins = assign_pseudotime_bins(pseudotime, args.n_bins)
            for module in scores.columns:
                summary = summarize_trend(scores[module], pseudotime, n_bins=args.n_bins, bins=bins)
                module_trend_rows.append({"run_id": run_id, "method": method, "feature_type": "module", "feature": module, **summary})
                module_bin_rows.append(bin_means_for_values(scores[module], pseudotime, bins, module, run_id, method))

            for gene in expr.columns:
                summary = summarize_trend(expr[gene], pseudotime, n_bins=args.n_bins, bins=bins)
                gene_trend_rows.append(
                    {
                        "run_id": run_id,
                        "method": method,
                        "gene": gene,
                        "marker_modules": ";".join(gene_to_modules.get(gene, [])),
                        **summary,
                    }
                )

    module_trends = pd.DataFrame(module_trend_rows)
    module_bins = pd.concat([df for df in module_bin_rows if not df.empty], ignore_index=True)
    gene_trends = pd.DataFrame(gene_trend_rows)
    gene_trends["abs_spearman_rho"] = gene_trends["spearman_rho"].abs()
    gene_trends = gene_trends.sort_values(["run_id", "method", "abs_spearman_rho"], ascending=[True, True, False])

    for run_id in module_trends["run_id"].unique():
        for method in module_trends["method"].unique():
            path = plot_module_trends(module_bins, run_id, method, args.figures_dir)
            if path:
                figure_paths.append(path)

    top_gene_keys = (
        gene_trends.groupby(["run_id", "method"], observed=True)
        .head(args.top_genes_per_method)[["run_id", "method", "gene"]]
        .drop_duplicates()
    )
    for _, row in top_gene_keys.iterrows():
        cached = run_cache[str(row["run_id"])]
        pt = cached["pt"]
        expr = cached["expr"]
        pseudotime = pd.to_numeric(pt[cached["methods"][row["method"]]], errors="coerce")
        bins = assign_pseudotime_bins(pseudotime, args.n_bins)
        gene_bin_rows.append(bin_means_for_values(expr[row["gene"]], pseudotime, bins, row["gene"], row["run_id"], row["method"]))

    module_trends_path = args.metadata_dir / "trajectory_module5_4_module_trends.tsv"
    module_bins_path = args.metadata_dir / "trajectory_module5_4_module_bin_means.tsv"
    gene_trends_path = args.metadata_dir / "trajectory_module5_4_gene_trends.tsv.gz"
    top_gene_bins_path = args.metadata_dir / "trajectory_module5_4_top_gene_bin_means.tsv.gz"
    module_scores_path = args.metadata_dir / "trajectory_module5_4_module_scores_by_cell.tsv.gz"
    availability_path = args.metadata_dir / "trajectory_module5_4_module_gene_availability.tsv"
    report_path = args.metadata_dir / "trajectory_module5_4_report.json"

    module_trends.to_csv(module_trends_path, sep="\t", index=False)
    module_bins.to_csv(module_bins_path, sep="\t", index=False)
    gene_trends.to_csv(gene_trends_path, sep="\t", index=False, compression="gzip")
    pd.concat([df for df in gene_bin_rows if not df.empty], ignore_index=True).to_csv(
        top_gene_bins_path, sep="\t", index=False, compression="gzip"
    )
    pd.concat(module_score_rows, ignore_index=True).to_csv(module_scores_path, sep="\t", index=False, compression="gzip")
    pd.DataFrame(availability_rows).to_csv(availability_path, sep="\t", index=False)

    report = {
        "module": "5.4",
        "method": "gene and marker-module trends along Monocle3 and Slingshot pseudotime",
        "input_h5ad": str(args.input_h5ad.resolve()),
        "n_bins": int(args.n_bins),
        "n_module_trends": int(module_trends.shape[0]),
        "n_gene_trends": int(gene_trends.shape[0]),
        "outputs": {
            "module_trends": str(module_trends_path.resolve()),
            "module_bin_means": str(module_bins_path.resolve()),
            "gene_trends": str(gene_trends_path.resolve()),
            "top_gene_bin_means": str(top_gene_bins_path.resolve()),
            "module_scores_by_cell": str(module_scores_path.resolve()),
            "module_gene_availability": str(availability_path.resolve()),
            "figures": figure_paths,
        },
        "elapsed_seconds": round(time.time() - start, 3),
    }
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
