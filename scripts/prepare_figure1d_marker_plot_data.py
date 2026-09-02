from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GENES = ["ALB", "KRT19", "PECAM1", "COL1A1", "LST1", "CD3D", "MS4A1", "MZB1", "S100A8"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare Figure 1D marker plot data on the sampled global UMAP cells.")
    parser.add_argument(
        "--selected-cells",
        type=Path,
        default=ROOT / "metadata/figure1c/figure1c_global_selected_cells.tsv",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=ROOT / "metadata/scvi/scvi_input_manifest.counts.tsv",
    )
    parser.add_argument(
        "--umap",
        type=Path,
        default=ROOT / "metadata/scvi/scvi_umap.tsv.gz",
    )
    parser.add_argument("--output-dir", type=Path, default=ROOT / "metadata/figure1")
    parser.add_argument("--genes", nargs="+", default=DEFAULT_GENES)
    return parser.parse_args()


def read_table(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", **kwargs)


def normalize_bool(value: object) -> bool:
    if value is None or pd.isna(value):
        return False
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def original_cell_id(cell_id: str, study_sample: str) -> str:
    prefix = f"{study_sample}__"
    return cell_id[len(prefix) :] if cell_id.startswith(prefix) else cell_id


def per_cell_log1p_cp10k(expr: sparse.spmatrix | np.ndarray, total_counts: np.ndarray) -> np.ndarray:
    total = total_counts.astype(np.float64, copy=True)
    total[total <= 0] = 1.0
    if sparse.issparse(expr):
        arr = expr.toarray().astype(np.float64, copy=False)
    else:
        arr = np.asarray(expr, dtype=np.float64)
    norm = arr / total[:, None] * 1e4
    return np.log1p(norm)


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    selected = read_table(args.selected_cells)
    manifest = read_table(args.manifest)
    include = manifest["include_in_scvi"].map(normalize_bool)
    manifest = manifest.loc[include].copy()
    if "study_sample" not in manifest.columns:
        manifest["study_sample"] = manifest["dataset"].astype(str) + "__" + manifest["label"].astype(str)
    source_map = manifest.set_index("study_sample")["output"].to_dict()

    umap = read_table(args.umap)
    if umap.columns[0] in {"", "Unnamed: 0"}:
        umap = umap.rename(columns={umap.columns[0]: "cell_id"})
    selected = selected.merge(umap[["cell_id", "UMAP_1", "UMAP_2"]], on="cell_id", how="left", validate="one_to_one")

    gene_rows: list[pd.DataFrame] = []
    summary_rows: list[dict[str, object]] = []

    for study_sample, sub in selected.groupby("study_sample", sort=True, observed=True):
        source_path = Path(str(source_map[study_sample]))
        if not source_path.exists():
            raise FileNotFoundError(f"Missing source h5ad for {study_sample}: {source_path}")
        adata = ad.read_h5ad(source_path, backed="r")
        original_ids = [original_cell_id(cell_id, study_sample) for cell_id in sub["cell_id"].astype(str)]
        original_to_global = pd.DataFrame({"original_cell_id": original_ids, "cell_id": sub["cell_id"].astype(str).to_list()})
        present = pd.Index(adata.obs_names.astype(str)).intersection(original_to_global["original_cell_id"])
        if present.empty:
            adata.file.close()
            continue
        keep_map = original_to_global[original_to_global["original_cell_id"].isin(present)].drop_duplicates("original_cell_id")
        keep_map = keep_map.set_index("original_cell_id").loc[present].reset_index().rename(columns={"index": "original_cell_id"})
        sub_meta = sub.set_index("cell_id").loc[keep_map["cell_id"]].reset_index()

        var_names = pd.Index(adata.var_names.astype(str))
        present_genes = [gene for gene in args.genes if gene in set(var_names)]
        if not present_genes:
            adata.file.close()
            continue
        gene_idx = var_names.get_indexer(present_genes)
        mem = adata[present, gene_idx].to_memory()
        totals = np.asarray(adata[present, :].X.sum(axis=1)).reshape(-1)
        adata.file.close()

        expr = per_cell_log1p_cp10k(mem.X, totals)
        expr_df = pd.DataFrame(expr, columns=present_genes)
        expr_df.insert(0, "original_cell_id", keep_map["original_cell_id"].to_numpy())
        merged = keep_map.merge(sub_meta, on="cell_id", how="left", validate="one_to_one")
        merged = merged.merge(expr_df, on="original_cell_id", how="left", validate="one_to_one")

        for gene in args.genes:
            if gene not in merged.columns:
                continue
            sub_gene = merged.loc[
                :,
                ["cell_id", "dataset", "study_sample", "major_celltype", "UMAP_1", "UMAP_2", gene],
            ].copy()
            sub_gene = sub_gene.rename(columns={gene: "expression"})
            sub_gene["gene"] = gene
            gene_rows.append(sub_gene)

            by_type = sub_gene.groupby("major_celltype", observed=True)["expression"]
            for major_celltype, vals in by_type:
                vals = vals.astype(float)
                summary_rows.append(
                    {
                        "gene": gene,
                        "major_celltype": major_celltype,
                        "n_cells": int(vals.shape[0]),
                        "pct_expr_gt0": float((vals > 0).mean()),
                        "mean_log1p_cp10k": float(vals.mean()),
                        "median_log1p_cp10k": float(vals.median()),
                    }
                )

    plot_df = pd.concat(gene_rows, axis=0, ignore_index=True)
    plot_df.to_csv(args.output_dir / "figure1D_marker_plot_data.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.DataFrame(summary_rows).to_csv(args.output_dir / "figure1D_marker_summary.tsv", sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
