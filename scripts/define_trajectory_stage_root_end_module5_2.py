from __future__ import annotations

import argparse
import json
import time
from importlib.metadata import version
from pathlib import Path
from typing import Mapping

import anndata as ad
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

SAMPLE_STAGE_ORDER = {
    "reference_non_hcc_liver": 0,
    "reference_adjacent_liver": 0,
    "chronic_liver": 1,
    "unknown_liver": 2,
    "hcc_dataset_unspecified_tissue": 3,
    "primary_hcc_tumor": 4,
    "pvtt_tumor": 5,
    "metastatic_lymphnode_tumor": 6,
    "unknown": 99,
}
CELL_STAGE_ORDER = {
    "stage_0_reference_hepatocyte": 0,
    "stage_1_stressed_injured": 1,
    "stage_2_regenerative_progenitor": 2,
    "stage_3_proliferating_candidate": 3,
    "stage_4_cnv_supported_malignant": 4,
    "stage_4_malignant_like_review": 4,
    "unresolved_epithelial_or_mixed": 9,
}
ROOT_END_ROLE_ORDER = {
    "root_reference": 0,
    "intermediate_trajectory": 1,
    "end_malignant_cnv": 2,
    "end_malignant_review": 3,
    "excluded_from_main_trajectory": 9,
}
ROOT_SOURCE_CLASSES = {"non_hcc_liver", "normal_adjacent"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 5.2: define disease stage and trajectory root/end cells.")
    parser.add_argument(
        "--input-h5ad",
        type=Path,
        default=ROOT / "data/processed/trajectory/trajectory_hepatocyte_cnv_scanvi.module5_1.h5ad",
    )
    parser.add_argument("--processed-dir", type=Path, default=ROOT / "data/processed/trajectory")
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/trajectory")
    parser.add_argument("--output-name", default="trajectory_hepatocyte_cnv_scanvi.stage_root_end.module5_2.h5ad")
    parser.add_argument("--use-rep", default="X_scANVI")
    parser.add_argument("--no-write-h5ad", action="store_true")
    return parser.parse_args()


def text_value(value: object, default: str = "") -> str:
    if value is None or pd.isna(value):
        return default
    text = str(value)
    if text.lower() in {"nan", "none", "<na>"}:
        return default
    return text


def bool_value(value: object) -> bool:
    if value is None or pd.isna(value):
        return False
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def assign_sample_disease_stage(row: Mapping[str, object]) -> str:
    source = text_value(row.get("sample_source_class"), "unknown")
    return {
        "non_hcc_liver": "reference_non_hcc_liver",
        "normal_adjacent": "reference_adjacent_liver",
        "cirrhotic_or_chronic_liver": "chronic_liver",
        "unknown_liver": "unknown_liver",
        "unknown_hcc_dataset": "hcc_dataset_unspecified_tissue",
        "tumor": "primary_hcc_tumor",
        "pvtt_tumor": "pvtt_tumor",
        "metastatic_tumor_lymphnode": "metastatic_lymphnode_tumor",
    }.get(source, "unknown")


def assign_cell_disease_stage(row: Mapping[str, object]) -> str:
    role = text_value(row.get("trajectory_role"), "ambiguous_epithelial_or_mixed")
    if role == "malignant_cnv_supported":
        return "stage_4_cnv_supported_malignant"
    if role in {"malignant_proliferation_cnv_review", "malignant_like_scanvi_review"}:
        return "stage_4_malignant_like_review"
    if role == "proliferating_candidate":
        return "stage_3_proliferating_candidate"
    if role == "regenerative_progenitor":
        return "stage_2_regenerative_progenitor"
    if role == "stressed_injured":
        return "stage_1_stressed_injured"
    if role == "normal_reference":
        return "stage_0_reference_hepatocyte"
    return "unresolved_epithelial_or_mixed"


def assign_root_end_role(row: Mapping[str, object]) -> str:
    role = text_value(row.get("trajectory_role"), "ambiguous_epithelial_or_mixed")
    source = text_value(row.get("sample_source_class"), "unknown")
    include_main = bool_value(row.get("trajectory_include_main"))
    include_strict = bool_value(row.get("trajectory_include_cnv_strict"))

    if role == "normal_reference" and source in ROOT_SOURCE_CLASSES and include_main:
        return "root_reference"
    if role == "malignant_cnv_supported" and include_strict:
        return "end_malignant_cnv"
    if role in {"malignant_proliferation_cnv_review", "malignant_like_scanvi_review"} and include_main:
        return "end_malignant_review"
    if include_main:
        return "intermediate_trajectory"
    return "excluded_from_main_trajectory"


def build_stage_table(obs: pd.DataFrame) -> pd.DataFrame:
    cells = obs.copy()
    cells.insert(0, "cell_id", cells.index.astype(str))
    rows = cells.to_dict(orient="records")
    cells["sample_disease_stage"] = [assign_sample_disease_stage(row) for row in rows]
    cells["sample_disease_stage_order"] = cells["sample_disease_stage"].map(SAMPLE_STAGE_ORDER).fillna(99).astype(int)
    cells["cell_disease_stage"] = [assign_cell_disease_stage(row) for row in rows]
    cells["cell_disease_stage_order"] = cells["cell_disease_stage"].map(CELL_STAGE_ORDER).fillna(99).astype(int)
    cells["trajectory_root_end_role"] = [assign_root_end_role(row) for row in cells.to_dict(orient="records")]
    cells["trajectory_root_end_role_order"] = cells["trajectory_root_end_role"].map(ROOT_END_ROLE_ORDER).fillna(99).astype(int)
    cells["trajectory_root_cell_candidate"] = cells["trajectory_root_end_role"].eq("root_reference")
    cells["trajectory_end_malignant_cnv_candidate"] = cells["trajectory_root_end_role"].eq("end_malignant_cnv")
    cells["trajectory_end_malignant_review_candidate"] = cells["trajectory_root_end_role"].eq("end_malignant_review")
    cells["trajectory_terminal_candidate"] = cells[
        ["trajectory_end_malignant_cnv_candidate", "trajectory_end_malignant_review_candidate"]
    ].any(axis=1)
    return cells


def choose_medoid_index(rep: np.ndarray, mask: np.ndarray) -> int:
    positions = np.flatnonzero(mask)
    if positions.size == 0:
        raise ValueError("No candidate cells available for medoid selection.")
    sub = np.asarray(rep[positions], dtype=np.float32)
    centroid = sub.mean(axis=0)
    distances = np.square(sub - centroid).sum(axis=1)
    return int(positions[int(np.argmin(distances))])


def choose_farthest_index(rep: np.ndarray, mask: np.ndarray, reference_vector: np.ndarray) -> int:
    positions = np.flatnonzero(mask)
    if positions.size == 0:
        raise ValueError("No candidate cells available for endpoint selection.")
    sub = np.asarray(rep[positions], dtype=np.float32)
    distances = np.square(sub - reference_vector).sum(axis=1)
    return int(positions[int(np.argmax(distances))])


def choose_root_end(adata: ad.AnnData, cells: pd.DataFrame, use_rep: str) -> dict[str, object]:
    if use_rep not in adata.obsm:
        fallback = "X_scVI"
        if fallback not in adata.obsm:
            raise KeyError(f"Neither {use_rep!r} nor {fallback!r} is present in adata.obsm.")
        use_rep = fallback
    rep = np.asarray(adata.obsm[use_rep], dtype=np.float32)
    root_mask = cells["trajectory_root_cell_candidate"].to_numpy(dtype=bool)
    end_mask = cells["trajectory_end_malignant_cnv_candidate"].to_numpy(dtype=bool)
    review_mask = cells["trajectory_end_malignant_review_candidate"].to_numpy(dtype=bool)
    if not end_mask.any() and review_mask.any():
        end_mask = review_mask

    root_index = choose_medoid_index(rep, root_mask)
    root_vector = rep[root_index]
    end_index = choose_farthest_index(rep, end_mask, root_vector)

    root_centroid = rep[root_mask].mean(axis=0)
    end_centroid = rep[end_mask].mean(axis=0)
    cells["trajectory_distance_to_root_centroid"] = np.sqrt(np.square(rep - root_centroid).sum(axis=1))
    cells["trajectory_distance_to_cnv_end_centroid"] = np.sqrt(np.square(rep - end_centroid).sum(axis=1))
    cells["trajectory_root_cell_selected"] = False
    cells["trajectory_end_cell_selected"] = False
    cells.loc[cells.index[root_index], "trajectory_root_cell_selected"] = True
    cells.loc[cells.index[end_index], "trajectory_end_cell_selected"] = True

    return {
        "use_rep": use_rep,
        "iroot": int(root_index),
        "root_cell_id": str(adata.obs_names[root_index]),
        "root_cell_stage": str(cells.iloc[root_index]["cell_disease_stage"]),
        "root_sample_stage": str(cells.iloc[root_index]["sample_disease_stage"]),
        "root_candidate_cells": int(root_mask.sum()),
        "end_index": int(end_index),
        "end_cell_id": str(adata.obs_names[end_index]),
        "end_cell_stage": str(cells.iloc[end_index]["cell_disease_stage"]),
        "end_sample_stage": str(cells.iloc[end_index]["sample_disease_stage"]),
        "end_cnv_candidate_cells": int(cells["trajectory_end_malignant_cnv_candidate"].sum()),
        "end_review_candidate_cells": int(cells["trajectory_end_malignant_review_candidate"].sum()),
    }


def attach_stage_obs(adata: ad.AnnData, cells: pd.DataFrame) -> None:
    aligned = cells.set_index("cell_id").loc[adata.obs_names.astype(str)]
    for col in [
        "sample_disease_stage",
        "cell_disease_stage",
        "trajectory_root_end_role",
    ]:
        adata.obs[col] = pd.Categorical(aligned[col].astype(str).to_numpy())
    for col in [
        "sample_disease_stage_order",
        "cell_disease_stage_order",
        "trajectory_root_end_role_order",
    ]:
        adata.obs[col] = aligned[col].to_numpy(dtype=int)
    for col in [
        "trajectory_root_cell_candidate",
        "trajectory_end_malignant_cnv_candidate",
        "trajectory_end_malignant_review_candidate",
        "trajectory_terminal_candidate",
        "trajectory_root_cell_selected",
        "trajectory_end_cell_selected",
    ]:
        adata.obs[col] = aligned[col].to_numpy(dtype=bool)
    for col in [
        "trajectory_distance_to_root_centroid",
        "trajectory_distance_to_cnv_end_centroid",
    ]:
        adata.obs[col] = aligned[col].to_numpy(dtype=float)


def export_tables(adata: ad.AnnData, cells: pd.DataFrame, metadata_dir: Path) -> dict[str, str]:
    metadata_dir.mkdir(parents=True, exist_ok=True)
    cells_path = metadata_dir / "trajectory_module5_2_stage_root_end_cells.tsv.gz"
    counts_path = metadata_dir / "trajectory_module5_2_stage_counts.tsv"
    sample_counts_path = metadata_dir / "trajectory_module5_2_stage_by_sample.tsv"
    root_path = metadata_dir / "trajectory_module5_2_root_cells.tsv.gz"
    end_path = metadata_dir / "trajectory_module5_2_end_cells.tsv.gz"

    export = cells.copy()
    if "X_umap_global" in adata.obsm:
        export["global_umap_1"] = adata.obsm["X_umap_global"][:, 0]
        export["global_umap_2"] = adata.obsm["X_umap_global"][:, 1]
    export.to_csv(cells_path, sep="\t", index=False, compression="gzip")

    counts = (
        cells.groupby(["sample_disease_stage", "cell_disease_stage", "trajectory_root_end_role"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["sample_disease_stage", "cell_disease_stage", "trajectory_root_end_role"])
    )
    counts.to_csv(counts_path, sep="\t", index=False)

    sample_counts = (
        cells.groupby(
            ["dataset", "study_sample", "cnv_sample", "sample_source_class", "sample_disease_stage", "cell_disease_stage"],
            observed=True,
        )
        .size()
        .reset_index(name="n_cells")
        .sort_values(["dataset", "study_sample", "sample_disease_stage", "cell_disease_stage"])
    )
    sample_counts.to_csv(sample_counts_path, sep="\t", index=False)

    keep_cols = [
        "cell_id",
        "dataset",
        "study_sample",
        "cnv_sample",
        "sample_source_class",
        "sample_disease_stage",
        "cell_disease_stage",
        "trajectory_root_end_role",
        "trajectory_distance_to_root_centroid",
        "trajectory_distance_to_cnv_end_centroid",
        "trajectory_root_cell_selected",
        "trajectory_end_cell_selected",
    ]
    cells.loc[cells["trajectory_root_cell_candidate"], keep_cols].to_csv(
        root_path, sep="\t", index=False, compression="gzip"
    )
    cells.loc[cells["trajectory_terminal_candidate"], keep_cols].to_csv(
        end_path, sep="\t", index=False, compression="gzip"
    )

    return {
        "cells": str(cells_path.resolve()),
        "stage_counts": str(counts_path.resolve()),
        "stage_by_sample": str(sample_counts_path.resolve()),
        "root_cells": str(root_path.resolve()),
        "end_cells": str(end_path.resolve()),
    }


def main() -> int:
    args = parse_args()
    start = time.time()
    args.processed_dir.mkdir(parents=True, exist_ok=True)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    output_h5ad = args.processed_dir / args.output_name
    report_path = args.metadata_dir / "trajectory_module5_2_report.json"

    print(f"READ {args.input_h5ad}", flush=True)
    adata = ad.read_h5ad(args.input_h5ad)
    required = ["trajectory_role", "trajectory_include_main", "trajectory_include_cnv_strict", "sample_source_class"]
    missing = [col for col in required if col not in adata.obs.columns]
    if missing:
        raise KeyError(f"Input AnnData is missing required module 5.1 obs columns: {missing}")

    print("DEFINE_STAGE_ROOT_END", flush=True)
    cells = build_stage_table(adata.obs)
    selection = choose_root_end(adata, cells, args.use_rep)
    attach_stage_obs(adata, cells)
    adata.uns["iroot"] = selection["iroot"]
    adata.uns["module5_2_stage_root_end"] = selection
    adata.uns["module5_2_stage_root_end"]["stage_note"] = (
        "Disease stage is derived from sample_source_class and trajectory_role; it is not a clinical TNM/BCLC stage."
    )

    exported = export_tables(adata, cells, args.metadata_dir)
    if not args.no_write_h5ad:
        print(f"WRITE {output_h5ad}", flush=True)
        adata.write_h5ad(output_h5ad, compression="gzip")

    report = {
        "module": "5.2",
        "method": "define sample/cell disease stage plus trajectory root/end candidates",
        "input_h5ad": str(args.input_h5ad.resolve()),
        "output_h5ad": str(output_h5ad.resolve()) if not args.no_write_h5ad else None,
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "sample_disease_stage_counts": cells["sample_disease_stage"].value_counts().to_dict(),
        "cell_disease_stage_counts": cells["cell_disease_stage"].value_counts().to_dict(),
        "root_end_role_counts": cells["trajectory_root_end_role"].value_counts().to_dict(),
        "selection": selection,
        "stage_note": adata.uns["module5_2_stage_root_end"]["stage_note"],
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
