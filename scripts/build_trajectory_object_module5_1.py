from __future__ import annotations

import argparse
import json
import re
import time
from importlib.metadata import version
from pathlib import Path
from typing import Mapping

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


ROOT = Path(__file__).resolve().parents[1]

MALIGNANT_SUPPORTED_CALLS = {
    "malignant_hcc_high_conf",
    "malignant_hcc_cnv_support",
    "malignant_hcc_probable",
}
MALIGNANT_REVIEW_CALLS = {
    "malignant_hcc_marker_proliferation_needs_cnv_review",
    "cnv_not_available",
}
ROLE_ORDER = {
    "normal_reference": 0,
    "stressed_injured": 1,
    "regenerative_progenitor": 2,
    "proliferating_candidate": 3,
    "malignant_cnv_supported": 4,
    "malignant_proliferation_cnv_review": 5,
    "malignant_like_scanvi_review": 6,
    "ambiguous_epithelial_or_mixed": 7,
}
MAIN_TRAJECTORY_ROLES = {
    "normal_reference",
    "stressed_injured",
    "regenerative_progenitor",
    "proliferating_candidate",
    "malignant_cnv_supported",
    "malignant_proliferation_cnv_review",
    "malignant_like_scanvi_review",
}
CNV_STRICT_TRAJECTORY_ROLES = {
    "normal_reference",
    "stressed_injured",
    "regenerative_progenitor",
    "proliferating_candidate",
    "malignant_cnv_supported",
}
ROOT_SOURCE_CLASSES = {"normal_adjacent", "non_hcc_liver"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 5.1: build the hepatocyte/CNV trajectory AnnData object.")
    parser.add_argument(
        "--scanvi-h5ad",
        type=Path,
        default=ROOT / "data/processed/scanvi/scanvi_unified_labels.h5ad",
    )
    parser.add_argument(
        "--hepatocyte-h5ad",
        type=Path,
        default=ROOT / "data/processed/hepatocyte/hepatocyte_lineage.global_scvi_subcluster.h5ad",
    )
    parser.add_argument(
        "--malignant-calls",
        type=Path,
        default=ROOT / "metadata/malignant/malignant_hcc_calls_by_cell.tsv.gz",
    )
    parser.add_argument("--processed-dir", type=Path, default=ROOT / "data/processed/trajectory")
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/trajectory")
    parser.add_argument("--output-name", default="trajectory_hepatocyte_cnv_scanvi.module5_1.h5ad")
    parser.add_argument("--use-rep", default="X_scANVI")
    parser.add_argument("--n-neighbors", type=int, default=30)
    parser.add_argument("--min-dist", type=float, default=0.25)
    parser.add_argument("--resolution", type=float, default=0.8)
    parser.add_argument("--seed", type=int, default=20260601)
    parser.add_argument("--skip-umap", action="store_true")
    parser.add_argument("--skip-paga", action="store_true")
    parser.add_argument("--keep-existing-x", action="store_true")
    return parser.parse_args()


def text_value(value: object, default: str = "") -> str:
    if value is None or pd.isna(value):
        return default
    text = str(value)
    if text.lower() in {"nan", "none", "<na>"}:
        return default
    return text


def sample_source_class(dataset: object, sample: object) -> str:
    dataset = text_value(dataset)
    sample = text_value(sample)
    up = sample.upper()
    if dataset == "GSE149614":
        if up.endswith("T"):
            return "tumor"
        if up.endswith("P"):
            return "pvtt_tumor"
        if up.endswith("L"):
            return "metastatic_tumor_lymphnode"
        if up.endswith("N"):
            return "normal_adjacent"
    if dataset == "GSE185477":
        if "TST" in up:
            return "tumor"
        if "NST" in up:
            return "normal_adjacent"
        if "CST" in up:
            return "cirrhotic_or_chronic_liver"
        if up.endswith("_SC") or "_SC" in up:
            return "unknown_liver"
    if dataset in {"GSE202379", "GSE174748"}:
        return "non_hcc_liver"
    if dataset in {"GSE151530", "GSE212046"}:
        return "unknown_hcc_dataset"
    return "unknown"


def original_id(cell_id: object, study_sample: object) -> str:
    cell_id_text = text_value(cell_id)
    prefix = f"{text_value(study_sample)}__"
    if cell_id_text.startswith(prefix):
        return cell_id_text[len(prefix) :]
    return cell_id_text


def derive_cnv_sample(row: Mapping[str, object]) -> str:
    explicit = text_value(row.get("cnv_sample"))
    if explicit:
        return explicit
    dataset = text_value(row.get("dataset"))
    if dataset == "GSE149614":
        match = re.match(r"^(HCC\d+[A-Z]?)_", original_id(row.get("cell_id"), row.get("study_sample")))
        if match:
            return match.group(1)
    sample_id = text_value(row.get("sample_id"))
    if sample_id:
        return sample_id
    return text_value(row.get("study_sample"), "unknown")


def assign_trajectory_role(row: Mapping[str, object]) -> str:
    call = text_value(row.get("malignant_hcc_call"), "not_module3_candidate")
    state = text_value(row.get("hepatocyte_state_label"), "ambiguous_epithelial_or_mixed")
    strict_label = text_value(row.get("scanvi_unified_final_strict_label"))

    if call in MALIGNANT_SUPPORTED_CALLS:
        return "malignant_cnv_supported"
    if call == "malignant_hcc_marker_proliferation_needs_cnv_review":
        return "malignant_proliferation_cnv_review"
    if strict_label == "malignant_like_hepatocyte_needs_review" or call == "cnv_not_available":
        return "malignant_like_scanvi_review"
    if state == "normal_hepatocyte_like":
        return "normal_reference"
    if state == "stressed_injured_hepatocyte":
        return "stressed_injured"
    if state == "regenerative_progenitor_like_hepatocyte":
        return "regenerative_progenitor"
    if state == "proliferating_hepatocyte_candidate":
        return "proliferating_candidate"
    return "ambiguous_epithelial_or_mixed"


def read_hepatocyte_obs(path: Path) -> pd.DataFrame:
    backed = ad.read_h5ad(path, backed="r")
    obs = backed.obs.copy()
    obs.insert(0, "cell_id", backed.obs_names.astype(str))
    backed.file.close()
    return obs


def read_scanvi_subset(path: Path, cell_ids: pd.Index) -> ad.AnnData:
    backed = ad.read_h5ad(path, backed="r")
    obs_names = pd.Index(backed.obs_names.astype(str))
    mask = obs_names.isin(set(cell_ids.astype(str)))
    adata = backed[mask, :].to_memory()
    backed.file.close()
    missing = cell_ids.difference(pd.Index(adata.obs_names.astype(str)))
    if len(missing) > 0:
        raise ValueError(f"{len(missing)} hepatocyte cells are missing from the scANVI object.")
    return adata[cell_ids.astype(str), :].copy()


def merge_obs_annotations(adata: ad.AnnData, hepatocyte_obs: pd.DataFrame, malignant_path: Path) -> pd.DataFrame:
    hep_cols = [
        "cell_id",
        "leiden_hep",
        "hepatocyte_state_label",
        "hepatocyte_state_confidence",
        "hepatocyte_state_seed_label",
        "module1_scanvi_seed_label_major",
        "manual_major_label_cluster",
        "module1_confidence_status",
    ]
    hep = hepatocyte_obs[[col for col in hep_cols if col in hepatocyte_obs.columns]].drop_duplicates("cell_id")

    malignant_cols = [
        "cell_id",
        "cnv_sample",
        "sample_source_class",
        "cnv_proxy_burden",
        "cnv_proxy_z",
        "cnv_proxy_high_bin_fraction",
        "cnv_proxy_max_abs_bin_log2",
        "cnv_proxy_status",
        "hcc_malignant_associated_score_z",
        "proliferation_score_z",
        "regenerative_progenitor_score_z",
        "copykat_pred",
        "copykat_status",
        "malignant_hcc_call",
        "malignant_hcc_cnv_method",
        "malignant_hcc_evidence",
        "malignant_hcc_evidence_copykat",
    ]
    malignant = pd.read_csv(malignant_path, sep="\t", usecols=lambda col: col in malignant_cols)
    malignant = malignant.drop_duplicates("cell_id")

    cells = pd.DataFrame({"cell_id": adata.obs_names.astype(str)})
    cells = cells.merge(hep, on="cell_id", how="left").merge(malignant, on="cell_id", how="left")

    for col in [
        "scanvi_unified_seed_label",
        "scanvi_unified_seed_source",
        "scanvi_unified_pred_label",
        "scanvi_unified_pred_max_prob",
        "scanvi_unified_pred_margin",
        "scanvi_unified_final_label",
        "scanvi_unified_final_source",
        "scanvi_unified_final_strict_label",
        "scanvi_unified_final_strict_source",
        "dataset",
        "study_sample",
        "sample_id",
        "leiden_scvi",
        "major_celltype",
        "predicted_doublet",
        "cell_cycle_phase",
        "qc_total_counts",
        "qc_n_genes_by_counts",
        "qc_pct_counts_mt",
    ]:
        if col in adata.obs.columns and col not in cells.columns:
            cells[col] = adata.obs[col].astype(str).to_numpy() if adata.obs[col].dtype.name == "category" else adata.obs[col].to_numpy()

    cells["malignant_hcc_call"] = cells["malignant_hcc_call"].fillna("not_module3_candidate").astype(str)
    cells["cnv_sample"] = [derive_cnv_sample(row) for row in cells.to_dict(orient="records")]
    derived_source = [sample_source_class(row["dataset"], row["cnv_sample"]) for row in cells.to_dict(orient="records")]
    cells["sample_source_class"] = cells["sample_source_class"].fillna(pd.Series(derived_source, index=cells.index)).astype(str)
    cells["trajectory_role"] = [assign_trajectory_role(row) for row in cells.to_dict(orient="records")]
    cells["trajectory_role_order"] = cells["trajectory_role"].map(ROLE_ORDER).fillna(99).astype(int)
    predicted_doublet = cells["predicted_doublet"].fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})
    cells["trajectory_include_main"] = cells["trajectory_role"].isin(MAIN_TRAJECTORY_ROLES) & ~predicted_doublet
    cells["trajectory_include_cnv_strict"] = cells["trajectory_role"].isin(CNV_STRICT_TRAJECTORY_ROLES) & ~predicted_doublet
    cells["trajectory_root_candidate"] = (
        cells["trajectory_role"].eq("normal_reference")
        & cells["sample_source_class"].isin(ROOT_SOURCE_CLASSES)
        & ~predicted_doublet
    )
    cells["trajectory_root_priority"] = np.select(
        [
            cells["trajectory_root_candidate"],
            cells["trajectory_role"].eq("normal_reference") & ~predicted_doublet,
            cells["trajectory_role"].isin(["stressed_injured", "regenerative_progenitor"]) & ~predicted_doublet,
        ],
        [0, 1, 2],
        default=9,
    ).astype(int)
    return cells


def attach_obs(adata: ad.AnnData, cells: pd.DataFrame) -> None:
    aligned = cells.set_index("cell_id").loc[adata.obs_names.astype(str)]
    for col in aligned.columns:
        values = aligned[col]
        if values.dtype == bool:
            adata.obs[col] = values.to_numpy(dtype=bool)
        elif pd.api.types.is_numeric_dtype(values):
            adata.obs[col] = values.to_numpy()
        elif col == "cell_id":
            continue
        else:
            text = values.astype(object).where(pd.notna(values), "Unknown").astype(str)
            adata.obs[col] = pd.Categorical(text.to_numpy())


def normalize_x_from_counts(adata: ad.AnnData) -> None:
    if "counts" not in adata.layers:
        raise KeyError("layers['counts'] is required unless --keep-existing-x is used.")
    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.uns["module5_1_x_basis"] = "log1p(normalize_total(layers['counts'], target_sum=1e4))"


def export_tables(adata: ad.AnnData, cells: pd.DataFrame, metadata_dir: Path) -> dict[str, str]:
    metadata_dir.mkdir(parents=True, exist_ok=True)
    cells_path = metadata_dir / "trajectory_module5_1_cells.tsv.gz"
    role_counts_path = metadata_dir / "trajectory_module5_1_role_counts.tsv"
    sample_counts_path = metadata_dir / "trajectory_module5_1_role_by_sample.tsv"
    root_candidates_path = metadata_dir / "trajectory_module5_1_root_candidates.tsv.gz"

    export = cells.copy()
    if "X_umap_trajectory" in adata.obsm:
        export["trajectory_umap_1"] = adata.obsm["X_umap_trajectory"][:, 0]
        export["trajectory_umap_2"] = adata.obsm["X_umap_trajectory"][:, 1]
    if "X_umap_global" in adata.obsm:
        export["global_umap_1"] = adata.obsm["X_umap_global"][:, 0]
        export["global_umap_2"] = adata.obsm["X_umap_global"][:, 1]
    export.to_csv(cells_path, sep="\t", index=False, compression="gzip")

    role_counts = (
        cells.groupby(
            [
                "trajectory_role",
                "trajectory_include_main",
                "trajectory_include_cnv_strict",
                "malignant_hcc_call",
            ],
            observed=True,
        )
        .size()
        .reset_index(name="n_cells")
        .sort_values(["trajectory_role", "malignant_hcc_call"])
    )
    role_counts.to_csv(role_counts_path, sep="\t", index=False)

    sample_counts = (
        cells.groupby(["dataset", "study_sample", "cnv_sample", "sample_source_class", "trajectory_role"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["dataset", "study_sample", "trajectory_role"])
    )
    sample_counts.to_csv(sample_counts_path, sep="\t", index=False)

    root_cols = [
        "cell_id",
        "dataset",
        "study_sample",
        "cnv_sample",
        "sample_source_class",
        "trajectory_role",
        "trajectory_root_priority",
    ]
    cells.loc[cells["trajectory_root_candidate"], root_cols].to_csv(
        root_candidates_path, sep="\t", index=False, compression="gzip"
    )

    return {
        "cells": str(cells_path.resolve()),
        "role_counts": str(role_counts_path.resolve()),
        "role_by_sample": str(sample_counts_path.resolve()),
        "root_candidates": str(root_candidates_path.resolve()),
    }


def main() -> int:
    args = parse_args()
    start = time.time()
    args.processed_dir.mkdir(parents=True, exist_ok=True)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    output_h5ad = args.processed_dir / args.output_name
    report_path = args.metadata_dir / "trajectory_module5_1_report.json"

    print(f"READ_HEPATOCYTE_OBS {args.hepatocyte_h5ad}", flush=True)
    hepatocyte_obs = read_hepatocyte_obs(args.hepatocyte_h5ad)
    cell_ids = pd.Index(hepatocyte_obs["cell_id"].astype(str))

    print(f"READ_SCANVI_SUBSET {args.scanvi_h5ad} cells={len(cell_ids)}", flush=True)
    adata = read_scanvi_subset(args.scanvi_h5ad, cell_ids)
    if args.use_rep not in adata.obsm:
        fallback = "X_scVI"
        if fallback not in adata.obsm:
            raise KeyError(f"Neither {args.use_rep!r} nor {fallback!r} is present in adata.obsm.")
        print(f"USE_REP_FALLBACK {args.use_rep} -> {fallback}", flush=True)
        args.use_rep = fallback

    if "X_umap" in adata.obsm:
        adata.obsm["X_umap_global"] = adata.obsm["X_umap"].copy()

    print("MERGE_MODULE2_MODULE3_MODULE4_ANNOTATIONS", flush=True)
    cells = merge_obs_annotations(adata, hepatocyte_obs, args.malignant_calls)
    attach_obs(adata, cells)

    if not args.keep_existing_x:
        print("NORMALIZE_X_FROM_COUNTS", flush=True)
        normalize_x_from_counts(adata)

    print(f"NEIGHBORS key=neighbors_trajectory use_rep={args.use_rep} n_neighbors={args.n_neighbors}", flush=True)
    sc.pp.neighbors(
        adata,
        n_neighbors=args.n_neighbors,
        use_rep=args.use_rep,
        random_state=args.seed,
        key_added="neighbors_trajectory",
    )
    if not args.skip_umap:
        print(f"UMAP key=X_umap_trajectory min_dist={args.min_dist}", flush=True)
        sc.tl.umap(
            adata,
            min_dist=args.min_dist,
            random_state=args.seed,
            neighbors_key="neighbors_trajectory",
            key_added="X_umap_trajectory",
        )
    else:
        print("SKIP_UMAP", flush=True)
    print(f"LEIDEN key=leiden_trajectory resolution={args.resolution}", flush=True)
    sc.tl.leiden(
        adata,
        resolution=args.resolution,
        key_added="leiden_trajectory",
        neighbors_key="neighbors_trajectory",
        random_state=args.seed,
        flavor="igraph",
        directed=False,
    )
    cells["leiden_trajectory"] = adata.obs["leiden_trajectory"].astype(str).to_numpy()

    if not args.skip_paga:
        print("PAGA groups=trajectory_role", flush=True)
        sc.tl.paga(adata, groups="trajectory_role", neighbors_key="neighbors_trajectory")

    print(f"WRITE {output_h5ad}", flush=True)
    adata.write_h5ad(output_h5ad, compression="gzip")
    exported = export_tables(adata, cells, args.metadata_dir)

    report = {
        "module": "5.1",
        "method": "hepatocyte-lineage trajectory object with module 3 CNV malignant calls and module 4 scANVI labels",
        "scanvi_h5ad": str(args.scanvi_h5ad.resolve()),
        "hepatocyte_h5ad": str(args.hepatocyte_h5ad.resolve()),
        "malignant_calls": str(args.malignant_calls.resolve()),
        "output_h5ad": str(output_h5ad.resolve()),
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "use_rep": args.use_rep,
        "n_neighbors": int(args.n_neighbors),
        "min_dist": float(args.min_dist),
        "resolution": float(args.resolution),
        "skip_umap": bool(args.skip_umap),
        "seed": int(args.seed),
        "role_counts": cells["trajectory_role"].value_counts().to_dict(),
        "main_include_cells": int(cells["trajectory_include_main"].sum()),
        "cnv_strict_include_cells": int(cells["trajectory_include_cnv_strict"].sum()),
        "root_candidate_cells": int(cells["trajectory_root_candidate"].sum()),
        "module3_call_counts": cells["malignant_hcc_call"].value_counts().to_dict(),
        "sample_source_class_counts": cells["sample_source_class"].value_counts().to_dict(),
        "scanpy_version": version("scanpy"),
        "anndata_version": version("anndata"),
        "elapsed_seconds": round(time.time() - start, 3),
        "outputs": exported,
    }
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
