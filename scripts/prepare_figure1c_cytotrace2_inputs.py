from __future__ import annotations

import argparse
import gzip
import math
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare sampled CytoTRACE2 inputs for Figure 1C global and hepatocyte-lineage plots."
    )
    parser.add_argument(
        "--global-annotations",
        type=Path,
        default=ROOT / "metadata/celltype/celltypist_major_by_cell.tsv.gz",
    )
    parser.add_argument(
        "--hepatocyte-cells",
        type=Path,
        default=ROOT / "metadata/hepatocyte/hepatocyte_lineage_cells.tsv.gz",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=ROOT / "metadata/scvi/scvi_input_manifest.counts.tsv",
    )
    parser.add_argument("--output-dir", type=Path, default=ROOT / "data/processed/figure1c_inputs")
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/figure1c")
    parser.add_argument("--global-target-cells", type=int, default=24177)
    parser.add_argument("--hep-target-cells", type=int, default=20000)
    parser.add_argument("--min-per-group", type=int, default=40)
    parser.add_argument("--seed", type=int, default=20260708)
    return parser.parse_args()


def read_tsv(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", **kwargs)


def original_cell_id(cell_id: str, study_sample: str) -> str:
    prefix = f"{study_sample}__"
    return cell_id[len(prefix) :] if cell_id.startswith(prefix) else cell_id


def normalize_bool(value: object) -> bool:
    if value is None or pd.isna(value):
        return False
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def stratified_sample(
    df: pd.DataFrame,
    group_cols: list[str],
    target_n: int,
    min_per_group: int,
    seed: int,
) -> pd.DataFrame:
    if target_n <= 0:
        raise ValueError("target_n must be positive")
    if df.empty:
        raise ValueError("Cannot sample from an empty dataframe")

    group_sizes = df.groupby(group_cols, observed=True).size().rename("n_cells").reset_index()
    if target_n >= int(group_sizes["n_cells"].sum()):
        return df.copy()

    weights = group_sizes["n_cells"] / group_sizes["n_cells"].sum()
    requested = weights * target_n
    allocated = np.floor(requested).astype(int)
    allocated = np.minimum(allocated, group_sizes["n_cells"].to_numpy())

    positive_groups = group_sizes["n_cells"].to_numpy() > 0
    allocated[(allocated == 0) & positive_groups] = np.minimum(
        min_per_group,
        group_sizes.loc[(allocated == 0) & positive_groups, "n_cells"].to_numpy(),
    )
    allocated = np.minimum(allocated, group_sizes["n_cells"].to_numpy())

    total = int(allocated.sum())
    fractional = requested - np.floor(requested)

    if total > target_n:
        order = np.argsort(fractional.to_numpy())
        for idx in order:
            removable = allocated[idx] - min(1, group_sizes.loc[idx, "n_cells"])
            if removable <= 0:
                continue
            drop = min(removable, total - target_n)
            allocated[idx] -= drop
            total -= drop
            if total == target_n:
                break
    elif total < target_n:
        room = group_sizes["n_cells"].to_numpy() - allocated
        order = np.argsort((-fractional).to_numpy())
        for idx in order:
            if room[idx] <= 0:
                continue
            add = min(room[idx], target_n - total)
            allocated[idx] += add
            total += add
            if total == target_n:
                break

    sample_plan = group_sizes.copy()
    sample_plan["target_n"] = allocated

    rng = np.random.default_rng(seed)
    pieces: list[pd.DataFrame] = []
    for row in sample_plan.itertuples(index=False):
        query = pd.Series(True, index=df.index)
        for col in group_cols:
            query &= df[col].astype(str).eq(str(getattr(row, col)))
        sub = df.loc[query].copy()
        if sub.empty or int(row.target_n) <= 0:
            continue
        if int(row.target_n) >= sub.shape[0]:
            pieces.append(sub)
            continue
        chosen = rng.choice(sub.index.to_numpy(), size=int(row.target_n), replace=False)
        pieces.append(sub.loc[np.sort(chosen)].copy())

    out = pd.concat(pieces, axis=0).drop_duplicates("cell_id").reset_index(drop=True)
    if out.shape[0] > target_n:
        out = out.sample(n=target_n, random_state=seed).sort_values("cell_id").reset_index(drop=True)
    return out


def write_expression_tsv(matrix: sparse.spmatrix | np.ndarray, gene_names: pd.Index, cell_names: list[str], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out_path, "wt", encoding="utf-8", newline="") as handle:
        handle.write("gene\t" + "\t".join(cell_names) + "\n")
        if sparse.issparse(matrix):
            mat = matrix.tocsr()
            for gene_idx, gene_name in enumerate(gene_names):
                row = mat.getrow(gene_idx).toarray().ravel()
                values = "\t".join(str(int(x)) if float(x).is_integer() else f"{float(x):.6f}" for x in row)
                handle.write(f"{gene_name}\t{values}\n")
        else:
            arr = np.asarray(matrix)
            for gene_idx, gene_name in enumerate(gene_names):
                row = arr[gene_idx]
                values = "\t".join(str(int(x)) if float(x).is_integer() else f"{float(x):.6f}" for x in row)
                handle.write(f"{gene_name}\t{values}\n")


def export_mode_inputs(
    selected: pd.DataFrame,
    manifest: pd.DataFrame,
    mode_name: str,
    output_root: Path,
) -> pd.DataFrame:
    mode_dir = output_root / mode_name
    mode_dir.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, object]] = []

    manifest_map = manifest.set_index("study_sample")["output"].to_dict()
    for study_sample, sub in selected.groupby("study_sample", sort=True, observed=True):
        source_path = Path(str(manifest_map[study_sample]))
        if not source_path.exists():
            raise FileNotFoundError(f"Source h5ad missing for {study_sample}: {source_path}")
        adata = ad.read_h5ad(source_path, backed="r")
        original_ids = [original_cell_id(cell_id, study_sample) for cell_id in sub["cell_id"].astype(str)]
        present = pd.Index(adata.obs_names.astype(str)).intersection(original_ids)
        if present.empty:
            adata.file.close()
            raise ValueError(f"No selected cells from {study_sample} found in {source_path}")
        adata_mem = adata[present, :].to_memory()
        adata.file.close()

        selected_map = pd.DataFrame(
            {
                "cell_id": sub["cell_id"].astype(str).to_list(),
                "original_cell_id": original_ids,
            }
        )
        selected_map = selected_map[selected_map["original_cell_id"].isin(present)].drop_duplicates("original_cell_id")
        selected_map = selected_map.set_index("original_cell_id").loc[adata_mem.obs_names.astype(str)].reset_index()

        X = adata_mem.X
        if sparse.issparse(X):
            gene_cell = X.T.tocsr()
        else:
            gene_cell = np.asarray(X).T

        export_dir = mode_dir / study_sample.replace("/", "_")
        expr_path = export_dir / "expression.tsv.gz"
        anno_path = export_dir / "cell_annotations.tsv"
        map_path = export_dir / "cell_map.tsv"
        write_expression_tsv(gene_cell, pd.Index(adata_mem.var_names.astype(str)), selected_map["cell_id"].astype(str).to_list(), expr_path)

        annotation_cols = [c for c in ["major_celltype", "hepatocyte_state_label"] if c in sub.columns]
        cell_annotations = selected_map.merge(sub, on="cell_id", how="left")
        if annotation_cols:
            cell_annotations.loc[:, ["cell_id", annotation_cols[0]]].to_csv(anno_path, sep="\t", index=False)
        else:
            selected_map.loc[:, ["cell_id"]].assign(annotation="selected").to_csv(anno_path, sep="\t", index=False)
        cell_annotations.to_csv(map_path, sep="\t", index=False)

        records.append(
            {
                "mode": mode_name,
                "study_sample": study_sample,
                "source_h5ad": str(source_path),
                "n_selected_cells": int(selected_map.shape[0]),
                "n_genes": int(adata_mem.n_vars),
                "expression_tsv_gz": str(expr_path),
                "annotation_tsv": str(anno_path),
                "cell_map_tsv": str(map_path),
            }
        )
    return pd.DataFrame(records)


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)

    global_ann = read_tsv(args.global_annotations)
    global_ann = global_ann.rename(columns={global_ann.columns[0]: "cell_id"})
    global_ann = global_ann.loc[~global_ann["excluded_doublet_cluster"].map(normalize_bool)].copy()
    global_selected = stratified_sample(
        global_ann,
        group_cols=["major_celltype", "study_sample"],
        target_n=args.global_target_cells,
        min_per_group=args.min_per_group,
        seed=args.seed,
    )
    global_selected.to_csv(args.metadata_dir / "figure1c_global_selected_cells.tsv", sep="\t", index=False)

    hep_cells = read_tsv(args.hepatocyte_cells)
    hep_cells = hep_cells.rename(columns={hep_cells.columns[0]: "cell_id"})
    hep_cells = hep_cells.loc[~hep_cells["predicted_doublet"].map(normalize_bool)].copy()
    hep_selected = stratified_sample(
        hep_cells,
        group_cols=["hepatocyte_state_label", "study_sample"],
        target_n=args.hep_target_cells,
        min_per_group=args.min_per_group,
        seed=args.seed + 1,
    )
    hep_selected.to_csv(args.metadata_dir / "figure1c_hepatocyte_selected_cells.tsv", sep="\t", index=False)

    manifest = read_tsv(args.manifest)
    include = manifest["include_in_scvi"].map(normalize_bool)
    manifest = manifest.loc[include].copy()
    if "study_sample" not in manifest.columns:
        manifest["study_sample"] = manifest["dataset"].astype(str) + "__" + manifest["label"].astype(str)

    export_tables = [
        export_mode_inputs(global_selected, manifest, "global", args.output_dir),
        export_mode_inputs(hep_selected, manifest, "hepatocyte", args.output_dir),
    ]
    export_manifest = pd.concat(export_tables, axis=0, ignore_index=True)
    export_manifest.to_csv(args.metadata_dir / "figure1c_cytotrace2_input_manifest.tsv", sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
