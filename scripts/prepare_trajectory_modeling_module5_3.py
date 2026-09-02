from __future__ import annotations

import argparse
import gzip
import json
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io, sparse
from sklearn.decomposition import TruncatedSVD


ROOT = Path(__file__).resolve().parents[1]
STAGE_PROGENITOR = {"stage_2_regenerative_progenitor", "stage_3_proliferating_candidate"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 5.3: prepare Monocle3 and Slingshot trajectory modeling inputs.")
    parser.add_argument(
        "--input-h5ad",
        type=Path,
        default=ROOT / "data/processed/trajectory/trajectory_hepatocyte_cnv_scanvi.stage_root_end.module5_2.h5ad",
    )
    parser.add_argument("--processed-dir", type=Path, default=ROOT / "data/processed/trajectory/module5_3")
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/trajectory")
    parser.add_argument("--max-cells-per-run", type=int, default=80000)
    parser.add_argument("--pca-components", type=int, default=30)
    parser.add_argument("--cluster-key", default="leiden_trajectory")
    parser.add_argument("--seed", type=int, default=20260601)
    parser.add_argument("--skip-counts", action="store_true")
    parser.add_argument("--manifest-name", default="trajectory_module5_3_modeling_manifest.tsv")
    parser.add_argument("--report-name", default="trajectory_module5_3_prepare_report.json")
    return parser.parse_args()


def bool_series(values: pd.Series) -> pd.Series:
    return values.fillna(False).astype(str).str.lower().isin({"true", "1", "yes", "y"})


def build_model_mask(cells: pd.DataFrame, variant: str) -> pd.Series:
    if variant == "main_strict":
        return bool_series(cells["trajectory_include_cnv_strict"])
    if variant == "include_review":
        return bool_series(cells["trajectory_include_main"])
    raise ValueError(f"Unsupported model variant: {variant}")


def select_start_end_clusters(cells: pd.DataFrame, cluster_key: str = "leiden_trajectory", min_cells: int = 50) -> dict[str, object]:
    required = {cluster_key, "trajectory_root_end_role", "cell_disease_stage"}
    missing = required.difference(cells.columns)
    if missing:
        raise KeyError(f"Missing columns for cluster selection: {sorted(missing)}")

    rows = []
    for cluster, sub in cells.groupby(cluster_key, observed=True):
        n_cells = sub.shape[0]
        if n_cells < min_cells:
            continue
        rows.append(
            {
                "cluster": str(cluster),
                "n_cells": int(n_cells),
                "root_fraction": float(sub["trajectory_root_end_role"].astype(str).eq("root_reference").mean()),
                "malignant_fraction": float(sub["trajectory_root_end_role"].astype(str).eq("end_malignant_cnv").mean()),
                "progenitor_fraction": float(sub["cell_disease_stage"].astype(str).isin(STAGE_PROGENITOR).mean()),
            }
        )
    summary = pd.DataFrame(rows)
    if summary.empty:
        raise ValueError("No clusters meet the minimum cell threshold for Slingshot start/end selection.")

    start = summary.sort_values(["root_fraction", "n_cells", "cluster"], ascending=[False, False, True]).iloc[0]
    malignant = summary.sort_values(["malignant_fraction", "n_cells", "cluster"], ascending=[False, False, True]).iloc[0]
    progenitor = summary.sort_values(["progenitor_fraction", "n_cells", "cluster"], ascending=[False, False, True]).iloc[0]
    end_clusters = []
    for cluster in [str(malignant["cluster"]), str(progenitor["cluster"])]:
        if cluster not in end_clusters:
            end_clusters.append(cluster)
    return {
        "start_cluster": str(start["cluster"]),
        "malignant_end_cluster": str(malignant["cluster"]),
        "progenitor_end_cluster": str(progenitor["cluster"]),
        "end_clusters": end_clusters,
        "cluster_summary": summary.sort_values("cluster").to_dict(orient="records"),
    }


def stratified_model_sample(
    cells: pd.DataFrame,
    max_cells: int,
    seed: int,
    strata_cols: list[str] | None = None,
) -> pd.DataFrame:
    if cells.shape[0] <= max_cells:
        return cells.copy()
    if max_cells <= 0:
        raise ValueError("max_cells must be positive.")
    strata_cols = strata_cols or ["trajectory_root_end_role", "cell_disease_stage"]
    rng = np.random.default_rng(seed)

    reserve_mask = bool_series(cells.get("trajectory_root_cell_selected", pd.Series(False, index=cells.index))) | bool_series(
        cells.get("trajectory_end_cell_selected", pd.Series(False, index=cells.index))
    )
    reserve = cells.loc[reserve_mask].copy()
    if reserve.shape[0] > max_cells:
        return reserve.sample(n=max_cells, random_state=seed).sort_index()

    remaining = cells.loc[~reserve_mask].copy()
    slots = max_cells - reserve.shape[0]
    if slots == 0 or remaining.empty:
        return reserve.sort_index()

    group_sizes = remaining.groupby(strata_cols, observed=True).size().reset_index(name="n")
    group_sizes["raw_alloc"] = slots * group_sizes["n"] / float(group_sizes["n"].sum())
    group_sizes["alloc"] = np.floor(group_sizes["raw_alloc"]).astype(int)
    group_sizes.loc[group_sizes["alloc"].eq(0) & group_sizes["n"].gt(0), "alloc"] = 1
    while group_sizes["alloc"].sum() > slots:
        candidates = group_sizes.loc[group_sizes["alloc"].gt(1)].copy()
        if candidates.empty:
            break
        idx = candidates.sort_values(["raw_alloc", "n"], ascending=[True, True]).index[0]
        group_sizes.loc[idx, "alloc"] -= 1
    while group_sizes["alloc"].sum() < slots:
        idx = group_sizes.sort_values(["raw_alloc", "n"], ascending=[False, False]).index[0]
        group_sizes.loc[idx, "alloc"] += 1

    sampled_parts = []
    for _, row in group_sizes.iterrows():
        mask = pd.Series(True, index=remaining.index)
        for col in strata_cols:
            mask &= remaining[col].astype(str).eq(str(row[col]))
        sub = remaining.loc[mask]
        n = min(int(row["alloc"]), sub.shape[0])
        if n > 0:
            sampled_parts.append(sub.sample(n=n, random_state=int(rng.integers(0, 2**31 - 1))))
    sampled = pd.concat([reserve] + sampled_parts, axis=0).drop_duplicates("cell_id")
    if sampled.shape[0] > max_cells:
        nonreserve = sampled.loc[
            ~(bool_series(sampled.get("trajectory_root_cell_selected", pd.Series(False, index=sampled.index))) | bool_series(
                sampled.get("trajectory_end_cell_selected", pd.Series(False, index=sampled.index))
            ))
        ]
        keep_extra = max_cells - reserve.shape[0]
        sampled = pd.concat([reserve, nonreserve.sample(n=keep_extra, random_state=seed)], axis=0)
    if sampled.shape[0] < max_cells:
        extra_pool = cells.loc[~cells["cell_id"].isin(set(sampled["cell_id"]))]
        extra_n = min(max_cells - sampled.shape[0], extra_pool.shape[0])
        if extra_n > 0:
            sampled = pd.concat([sampled, extra_pool.sample(n=extra_n, random_state=seed)], axis=0)
    return sampled.sort_index()


def write_dataframe(path: Path, df: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, compression="gzip" if path.suffix == ".gz" else None)


def write_embedding(path: Path, cell_ids: pd.Index, values: np.ndarray, prefix: str) -> None:
    columns = [f"{prefix}_{idx + 1}" for idx in range(values.shape[1])]
    df = pd.DataFrame(values, columns=columns)
    df.insert(0, "cell_id", cell_ids.astype(str).to_numpy())
    write_dataframe(path, df)


def write_counts_mtx(path: Path, matrix: sparse.spmatrix) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wb") as handle:
        io.mmwrite(handle, matrix)


def variant_definitions() -> list[dict[str, str]]:
    return [
        {
            "run_id": "main_strict",
            "variant": "main_strict",
            "description": "CNV-strict trajectory: normal mature hepatocytes to CNV-supported malignant cells.",
        },
        {
            "run_id": "sensitivity_include_review",
            "variant": "include_review",
            "description": "Sensitivity trajectory including malignant-like review endpoints.",
        },
    ]


def prepare_run(
    adata: ad.AnnData,
    cells: pd.DataFrame,
    run: dict[str, str],
    args: argparse.Namespace,
) -> dict[str, object]:
    mask = build_model_mask(cells, run["variant"])
    selected = cells.loc[mask].copy()
    sampled = stratified_model_sample(
        selected,
        max_cells=args.max_cells_per_run,
        seed=args.seed,
        strata_cols=["trajectory_root_end_role", "cell_disease_stage"],
    )
    sampled = sampled.copy()
    sampled["module5_3_run_id"] = run["run_id"]
    sampled["module5_3_variant"] = run["variant"]
    cluster_selection = select_start_end_clusters(sampled, cluster_key=args.cluster_key, min_cells=20)
    sampled["slingshot_cluster"] = sampled[args.cluster_key].astype(str)
    sampled["monocle3_root_cell"] = sampled["trajectory_root_end_role"].astype(str).eq("root_reference")
    sampled["slingshot_start_cluster"] = cluster_selection["start_cluster"]
    sampled["slingshot_end_clusters"] = ",".join(cluster_selection["end_clusters"])

    run_dir = args.processed_dir / run["run_id"]
    run_dir.mkdir(parents=True, exist_ok=True)
    selected_index = pd.Index(sampled.index)
    sub = adata[selected_index, :].copy()

    write_dataframe(run_dir / "cell_metadata.tsv.gz", sampled.reset_index(drop=True))
    genes = pd.DataFrame(
        {
            "gene_id": sub.var_names.astype(str),
            "gene_short_name": sub.var_names.astype(str),
        }
    )
    genes.to_csv(run_dir / "genes.tsv", sep="\t", index=False)
    if not args.skip_counts:
        counts = sub.layers["counts"] if "counts" in sub.layers else sub.X
        counts = sparse.csr_matrix(counts).transpose().tocoo()
        write_counts_mtx(run_dir / "counts_gene_by_cell.mtx.gz", counts)

    write_embedding(run_dir / "embedding_x_scanvi.tsv.gz", pd.Index(sub.obs_names), np.asarray(sub.obsm["X_scANVI"]), "SCANVI")
    if "X_umap_global" in sub.obsm:
        write_embedding(run_dir / "embedding_umap.tsv.gz", pd.Index(sub.obs_names), np.asarray(sub.obsm["X_umap_global"]), "UMAP")
    else:
        write_embedding(run_dir / "embedding_umap.tsv.gz", pd.Index(sub.obs_names), np.asarray(sub.obsm["X_umap"]), "UMAP")

    svd = TruncatedSVD(n_components=args.pca_components, random_state=args.seed)
    x = sparse.csr_matrix(sub.X)
    pca = svd.fit_transform(x)
    write_embedding(run_dir / "embedding_hepatocyte_pca.tsv.gz", pd.Index(sub.obs_names), pca, "HepPCA")

    cluster_summary_path = run_dir / "slingshot_cluster_summary.tsv"
    pd.DataFrame(cluster_selection["cluster_summary"]).to_csv(cluster_summary_path, sep="\t", index=False)
    config = {
        "run_id": run["run_id"],
        "variant": run["variant"],
        "description": run["description"],
        "n_selected_before_sampling": int(selected.shape[0]),
        "n_cells": int(sampled.shape[0]),
        "n_genes": int(sub.n_vars),
        "max_cells_per_run": int(args.max_cells_per_run),
        "cluster_key": args.cluster_key,
        "monocle3": {
            "method": "learn_graph + order_cells",
            "root_cells_column": "monocle3_root_cell",
            "root_cells_definition": "normal mature hepatocytes: trajectory_root_end_role == root_reference",
        },
        "slingshot": {
            "cluster_column": "slingshot_cluster",
            "start_cluster": cluster_selection["start_cluster"],
            "end_clusters": cluster_selection["end_clusters"],
            "reduced_dims": ["embedding_x_scanvi.tsv.gz", "embedding_hepatocyte_pca.tsv.gz"],
            "start_definition": "cluster enriched for root_reference normal mature hepatocytes",
            "end_definition": "clusters enriched for malignant CNV-supported and regenerative/proliferating progenitor cells",
        },
        "outputs": {
            "cell_metadata": str((run_dir / "cell_metadata.tsv.gz").resolve()),
            "genes": str((run_dir / "genes.tsv").resolve()),
            "counts_gene_by_cell": str((run_dir / "counts_gene_by_cell.mtx.gz").resolve()) if not args.skip_counts else None,
            "embedding_x_scanvi": str((run_dir / "embedding_x_scanvi.tsv.gz").resolve()),
            "embedding_hepatocyte_pca": str((run_dir / "embedding_hepatocyte_pca.tsv.gz").resolve()),
            "embedding_umap": str((run_dir / "embedding_umap.tsv.gz").resolve()),
            "slingshot_cluster_summary": str(cluster_summary_path.resolve()),
        },
    }
    (run_dir / "model_config.json").write_text(json.dumps(config, indent=2, ensure_ascii=False), encoding="utf-8")
    return config


def main() -> int:
    args = parse_args()
    start = time.time()
    args.processed_dir.mkdir(parents=True, exist_ok=True)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    print(f"READ {args.input_h5ad}", flush=True)
    adata = ad.read_h5ad(args.input_h5ad)
    required = [
        "trajectory_include_cnv_strict",
        "trajectory_include_main",
        "trajectory_root_end_role",
        "cell_disease_stage",
        args.cluster_key,
    ]
    missing = [col for col in required if col not in adata.obs.columns]
    if missing:
        raise KeyError(f"Input AnnData is missing module 5.2 columns: {missing}")
    if "X_scANVI" not in adata.obsm:
        raise KeyError("adata.obsm['X_scANVI'] is required for module 5.3.")

    cells = adata.obs.copy()
    cells.insert(0, "cell_id", adata.obs_names.astype(str))
    run_reports = []
    for run in variant_definitions():
        print(f"PREPARE_RUN {run['run_id']}", flush=True)
        run_reports.append(prepare_run(adata, cells, run, args))

    manifest = pd.DataFrame(
        [
            {
                "run_id": report["run_id"],
                "variant": report["variant"],
                "n_cells": report["n_cells"],
                "n_genes": report["n_genes"],
                "start_cluster": report["slingshot"]["start_cluster"],
                "end_clusters": ",".join(report["slingshot"]["end_clusters"]),
                "run_dir": str((args.processed_dir / report["run_id"]).resolve()),
            }
            for report in run_reports
        ]
    )
    manifest_path = args.metadata_dir / args.manifest_name
    manifest.to_csv(manifest_path, sep="\t", index=False)
    report = {
        "module": "5.3",
        "method": "prepare Monocle3 and Slingshot trajectory modeling inputs with sensitivity variants",
        "input_h5ad": str(args.input_h5ad.resolve()),
        "processed_dir": str(args.processed_dir.resolve()),
        "max_cells_per_run": int(args.max_cells_per_run),
        "pca_components": int(args.pca_components),
        "runs": run_reports,
        "manifest": str(manifest_path.resolve()),
        "elapsed_seconds": round(time.time() - start, 3),
    }
    report_path = args.metadata_dir / args.report_name
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
