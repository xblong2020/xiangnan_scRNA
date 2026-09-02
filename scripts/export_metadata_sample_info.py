from __future__ import annotations

import gzip
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
QC_ROOT = ROOT / "data" / "processed" / "qc_h5ad"
OUT_ROOT = ROOT / "metadata" / "sample_info"
CELL_META = OUT_ROOT / "cell_metadata"
SAMPLE_SUMMARY = OUT_ROOT / "sample_summary"
RAW_META = OUT_ROOT / "raw_metadata"


SAMPLE_COL_PRIORITY = [
    "sample",
    "Sample",
    "orig.ident",
    "S_ID",
    "PatientID",
    "Patient.ID",
    "patient",
    "Cell",
]
PATIENT_COL_PRIORITY = ["PatientID", "Patient.ID", "patient", "S_ID"]
GROUP_COL_CANDIDATES = [
    "condition",
    "Disease.status",
    "Tissue",
    "site",
    "stage",
    "Stage",
    "virus",
    "Virus",
    "celltype",
    "annotation",
    "annotation_refined",
    "cell.annotation",
    "Global_Cluster",
    "Sub_Cluster",
    "Release_Global_Cluster",
    "Release_Cluster",
    "Type",
]


def choose_col(columns: pd.Index, candidates: list[str]) -> str | None:
    for col in candidates:
        if col in columns:
            return col
    return None


def safe_str_series(series: pd.Series) -> pd.Series:
    return series.astype("string").fillna("NA").astype(str)


def dataset_obs(dataset_dir: Path) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for path in sorted(dataset_dir.glob("*.h5ad")):
        obj = ad.read_h5ad(path, backed="r")
        obs = obj.obs.copy()
        obs.insert(0, "cell_id", obs.index.astype(str))
        obs.insert(1, "source_file", path.name)
        obs.insert(2, "metadata_dataset", dataset_dir.name)
        frames.append(obs.reset_index(drop=True))
        obj.file.close()
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True, sort=False)


def summarize_values(frame: pd.DataFrame, group_col: str) -> pd.DataFrame:
    grouped = frame.groupby(group_col, dropna=False)
    summary = grouped.size().rename("n_cells").reset_index()
    used_columns = {c.lower() for c in summary.columns}

    for col in ["qc_total_counts", "qc_n_genes_by_counts", "qc_pct_counts_mt", "nCount_RNA", "nFeature_RNA"]:
        if col in frame.columns:
            values = pd.to_numeric(frame[col], errors="coerce")
            stat = (
                frame.assign(_value=values)
                .groupby(group_col, dropna=False)["_value"]
                .agg(["median", "mean", "min", "max"])
                .reset_index()
                .rename(
                    columns={
                        "median": f"{col}_median",
                        "mean": f"{col}_mean",
                        "min": f"{col}_min",
                        "max": f"{col}_max",
                    }
                )
            )
            summary = summary.merge(stat, on=group_col, how="left")

    for col in GROUP_COL_CANDIDATES:
        if col == group_col or col not in frame.columns:
            continue
        nunique = frame[col].nunique(dropna=True)
        if nunique <= 50:
            value_col = f"{col}_values"
            suffix = 2
            while value_col.lower() in used_columns:
                value_col = f"{col}_values_{suffix}"
                suffix += 1
            used_columns.add(value_col.lower())
            collapsed = (
                frame.groupby(group_col, dropna=False)[col]
                .apply(lambda x: ";".join(sorted(map(str, pd.unique(x.dropna()))))[:2000])
                .rename(value_col)
                .reset_index()
            )
            summary = summary.merge(collapsed, on=group_col, how="left")
    return summary


def composition_table(frame: pd.DataFrame, sample_col: str, dataset: str) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for col in ["celltype", "annotation_refined", "annotation", "cell.annotation", "Global_Cluster", "Release_Global_Cluster", "Type"]:
        if col not in frame.columns:
            continue
        table = (
            frame.groupby([sample_col, col], dropna=False)
            .size()
            .rename("n_cells")
            .reset_index()
            .rename(columns={sample_col: "sample", col: "category"})
        )
        table.insert(0, "dataset", dataset)
        table.insert(2, "category_type", col)
        rows.append(table)
    if rows:
        return pd.concat(rows, ignore_index=True, sort=False)
    return pd.DataFrame(columns=["dataset", "sample", "category_type", "category", "n_cells"])


def write_dataset_outputs(dataset_dir: Path) -> dict[str, object]:
    dataset = dataset_dir.name
    frame = dataset_obs(dataset_dir)
    if frame.empty:
        return {"dataset": dataset, "status": "empty"}

    cell_path = CELL_META / f"{dataset}_cell_metadata.tsv.gz"
    frame.to_csv(cell_path, sep="\t", index=False, compression="gzip")

    sample_col = choose_col(frame.columns, SAMPLE_COL_PRIORITY)
    patient_col = choose_col(frame.columns, PATIENT_COL_PRIORITY)
    if sample_col is None:
        frame["_sample_id"] = frame["qc_file"].astype(str)
        sample_col = "_sample_id"
    frame[sample_col] = safe_str_series(frame[sample_col])

    sample_summary = summarize_values(frame, sample_col).rename(columns={sample_col: "sample"})
    sample_summary.insert(0, "dataset", dataset)
    sample_path = SAMPLE_SUMMARY / f"{dataset}_sample_summary.tsv"
    sample_summary.to_csv(sample_path, sep="\t", index=False)

    patient_path = ""
    if patient_col is not None:
        frame[patient_col] = safe_str_series(frame[patient_col])
        patient_summary = summarize_values(frame, patient_col).rename(columns={patient_col: "patient"})
        patient_summary.insert(0, "dataset", dataset)
        patient_path = SAMPLE_SUMMARY / f"{dataset}_patient_summary.tsv"
        patient_summary.to_csv(patient_path, sep="\t", index=False)

    comp = composition_table(frame, sample_col, dataset)
    comp_path = SAMPLE_SUMMARY / f"{dataset}_sample_composition.tsv"
    comp.to_csv(comp_path, sep="\t", index=False)

    return {
        "dataset": dataset,
        "status": "complete",
        "n_cells": len(frame),
        "n_columns": frame.shape[1],
        "sample_col": sample_col,
        "patient_col": patient_col or "",
        "cell_metadata": str(cell_path),
        "sample_summary": str(sample_path),
        "patient_summary": str(patient_path) if patient_path else "",
        "sample_composition": str(comp_path),
    }


def copy_raw_metadata() -> pd.DataFrame:
    files = [
        ("HCC_atlas", ROOT / "data/public/figshare/HCC_atlas/HCC_atlas_metadata_batch_effect.csv"),
        ("GSE149614", ROOT / "data/public/geo/GSE149614/GSE149614_HCC.metadata.updated.txt.gz"),
        ("GSE151530", ROOT / "data/public/geo/GSE151530/GSE151530_Info.txt.gz"),
        ("GSE185477", ROOT / "data/public/geo/GSE185477/GSE185477_Final_Metadata.txt.gz"),
        ("GSE202379", ROOT / "data/public/geo/GSE202379/filelist.txt"),
        ("GSE174748", ROOT / "data/public/geo/GSE174748/filelist.txt"),
        ("GSE185477", ROOT / "data/public/geo/GSE185477/filelist.txt"),
        ("GSE212046", ROOT / "data/public/geo/GSE212046/filelist.txt"),
    ]
    rows = []
    for dataset, src in files:
        if not src.exists():
            rows.append({"dataset": dataset, "source": str(src), "output": "", "status": "missing"})
            continue
        suffix = "".join(src.suffixes)
        out = RAW_META / f"{dataset}_{src.name}"
        shutil.copy2(src, out)
        rows.append({"dataset": dataset, "source": str(src), "output": str(out), "status": "copied"})
    return pd.DataFrame(rows)


def dedupe_columns_case_insensitive(columns: list[str]) -> list[str]:
    used: set[str] = set()
    out: list[str] = []
    for col in columns:
        candidate = str(col)
        suffix = 2
        while candidate.lower() in used:
            candidate = f"{col}_{suffix}"
            suffix += 1
        used.add(candidate.lower())
        out.append(candidate)
    return out


def main() -> int:
    for folder in [CELL_META, SAMPLE_SUMMARY, RAW_META]:
        folder.mkdir(parents=True, exist_ok=True)

    rows = []
    for dataset_dir in sorted([p for p in QC_ROOT.iterdir() if p.is_dir()]):
        print(f"EXPORT {dataset_dir.name}", flush=True)
        rows.append(write_dataset_outputs(dataset_dir))

    index = pd.DataFrame(rows)
    raw_index = copy_raw_metadata()

    index_path = OUT_ROOT / "metadata_export_index.tsv"
    raw_index_path = RAW_META / "raw_metadata_index.tsv"
    index.to_csv(index_path, sep="\t", index=False)
    raw_index.to_csv(raw_index_path, sep="\t", index=False)

    combined_sample = []
    for p in sorted(SAMPLE_SUMMARY.glob("*_sample_summary.tsv")):
        combined_sample.append(pd.read_csv(p, sep="\t"))
    if combined_sample:
        combined = pd.concat(combined_sample, ignore_index=True, sort=False)
        combined.columns = dedupe_columns_case_insensitive(list(combined.columns))
        combined.to_csv(SAMPLE_SUMMARY / "all_datasets_sample_summary.tsv", sep="\t", index=False)

    print(f"WROTE {index_path}", flush=True)
    print(f"WROTE {raw_index_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
