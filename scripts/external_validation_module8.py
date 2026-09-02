from __future__ import annotations

import argparse
import json
import math
import platform
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA_DIR = PROJECT_ROOT / "metadata/driver"
DEFAULT_FIGURE_DIR = PROJECT_ROOT / "figures/driver"
DEFAULT_LOCAL_SCRNA_ROOT = Path(r"G:\wanyi_HCC_scRNA\HCCscRNA")
DEFAULT_TCGA_EXPRESSION_PATH = Path(r"F:\Charley\2024\134panImmune\03.download\expression\TCGA-LIHC.htseq_fpkm.tsv.gz")
DEFAULT_TCGA_SURVIVAL_PATH = Path(r"F:\Charley\2024\134panImmune\03.download\survival\TCGA-LIHC.survival.tsv.gz")
DEFAULT_CLINICAL_ROOT = Path(r"G:\万亿肝癌")

MODULE8_REQUIRED_MANIFEST_COLUMNS = [
    "dataset_id",
    "source",
    "modality",
    "cancer_type",
    "species",
    "n_samples",
    "n_cells",
    "tissue_type",
    "access_url",
    "raw_available",
    "processed_available",
    "included",
    "exclusion_reason",
]

MODULE8_AXES = [
    "tier1_rescue",
    "ap1_stress_proliferation",
    "sox4_state_specific",
    "control_calibration",
]


def read_tsv_or_empty(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def normalize_gene_symbols(genes: Iterable[object]) -> list[str]:
    normalized = []
    for gene in genes:
        if pd.isna(gene):
            normalized.append("")
            continue
        value = str(gene).strip().upper()
        value = re.sub(r"\.\d+$", "", value)
        normalized.append(value)
    return normalized


def strip_ensembl_version(gene_id: object) -> str:
    if pd.isna(gene_id):
        return ""
    return re.sub(r"\.\d+$", "", str(gene_id).strip())


def infer_tcga_sample_type(samples: Iterable[object]) -> list[str]:
    sample_types = []
    for sample in samples:
        parts = str(sample).split("-")
        if len(parts) >= 4:
            code = parts[3][:2]
            if code in {"01", "02", "03", "05"}:
                sample_types.append("tumor")
            elif code in {"10", "11", "12", "13", "14"}:
                sample_types.append("normal")
            else:
                sample_types.append("unknown")
        else:
            sample_types.append("unknown")
    return sample_types


def infer_icgc_sample_type(samples: Iterable[object]) -> list[str]:
    sample_types = []
    for sample in samples:
        value = str(sample)
        if value.endswith("-T"):
            sample_types.append("tumor")
        elif value.endswith("-N"):
            sample_types.append("normal")
        else:
            sample_types.append("unknown")
    return sample_types


def sample_to_patient_id(sample: object) -> str:
    value = str(sample)
    if value.startswith("TCGA-"):
        return "-".join(value.split("-")[:3])
    match = re.search(r"(DO\d+)", value)
    if match:
        return match.group(1)
    return value


def load_tcga_clinical_table(path: Path) -> pd.DataFrame:
    """Read TCGA clinical tables saved either as tab text or Excel."""
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t")
    except (UnicodeDecodeError, pd.errors.ParserError):
        pass
    try:
        return pd.read_excel(path)
    except Exception:
        return pd.read_csv(path, sep="\t", engine="python")


def _ordinal_from_text(value: object, prefix: str | None = None) -> float:
    if pd.isna(value):
        return np.nan
    text = str(value).strip().upper()
    if text in {"", "NA", "NAN", "NONE", "NOT REPORTED", "UNKNOWN"}:
        return np.nan
    if prefix:
        match = re.search(rf"{re.escape(prefix.upper())}\s*(\d+)", text)
        if match:
            return float(match.group(1))
    match = re.search(r"\b([IVX]+)[A-C]?\b", text)
    if match:
        roman = match.group(1)
        roman_map = {"I": 1, "II": 2, "III": 3, "IV": 4, "V": 5}
        if roman in roman_map:
            return float(roman_map[roman])
    match = re.search(r"(\d+)", text)
    if match:
        return float(match.group(1))
    return np.nan


def prepare_tcga_clinical_covariates(clinical: pd.DataFrame) -> pd.DataFrame:
    if clinical.empty:
        return clinical.copy()
    prepared = clinical.copy()
    if "Id" in prepared.columns:
        prepared["Id"] = prepared["Id"].astype(str).str.slice(0, 12)
    numeric_columns = ["age", "futime", "fustat"]
    for column in numeric_columns:
        if column in prepared.columns:
            prepared[column] = pd.to_numeric(prepared[column], errors="coerce")
    if "gender" in prepared.columns:
        gender_text = prepared["gender"].astype(str).str.strip().str.upper()
        mapped_gender = gender_text.map({"MALE": 1, "M": 1, "1": 1, "FEMALE": 0, "F": 0, "0": 0})
        prepared["gender"] = pd.to_numeric(mapped_gender, errors="coerce")
    if "grade" in prepared.columns:
        prepared["grade"] = prepared["grade"].map(lambda value: _ordinal_from_text(value, prefix="G"))
    if "stage" in prepared.columns:
        prepared["stage"] = prepared["stage"].map(_ordinal_from_text)
    for column in ["T", "M", "N"]:
        if column in prepared.columns:
            prepared[column] = prepared[column].map(lambda value, c=column: _ordinal_from_text(value, prefix=c))
    return prepared


def discover_local_scrna_sources(local_root: str | Path | None) -> pd.DataFrame:
    columns = [
        "dataset_id",
        "source",
        "modality",
        "input_format",
        "expression_path",
        "metadata_path",
        "sample_id",
        "status",
    ]
    if local_root is None:
        return pd.DataFrame(columns=columns)
    root = Path(local_root)
    if not root.exists():
        return pd.DataFrame(columns=columns)
    rows = []
    gse156625_h5ad = root / "GSE156625-HCC/python-data/GSE156625_HCCscanpyobj.h5ad"
    if gse156625_h5ad.exists():
        rows.append(
            {
                "dataset_id": "GSE156625",
                "source": "local:G:/wanyi_HCC_scRNA/HCCscRNA",
                "modality": "scRNA-seq",
                "input_format": "h5ad",
                "expression_path": str(gse156625_h5ad),
                "metadata_path": str(gse156625_h5ad),
                "sample_id": "GSE156625_all",
                "status": "ready",
            }
        )
    cnp_expression = root / "CNP0000650/HCC_log_tpm_expression_matrix.txt.gz"
    cnp_metadata = root / "CNP0000650/HCC_cell_metadata.txt"
    if cnp_expression.exists() and cnp_metadata.exists():
        rows.append(
            {
                "dataset_id": "CNP0000650",
                "source": "local:G:/wanyi_HCC_scRNA/HCCscRNA",
                "modality": "scRNA-seq",
                "input_format": "gene_by_cell_tsv_gz",
                "expression_path": str(cnp_expression),
                "metadata_path": str(cnp_metadata),
                "sample_id": "CNP0000650_all",
                "status": "ready",
            }
        )
    for set_id in ["Set1", "Set2"]:
        matrix = root / f"GSE125449/GSE125449_{set_id}_matrix.mtx.gz"
        genes = root / f"GSE125449/GSE125449_{set_id}_genes.tsv.gz"
        barcodes = root / f"GSE125449/GSE125449_{set_id}_barcodes.tsv.gz"
        samples = root / f"GSE125449/GSE125449_{set_id}_samples.txt.gz"
        if matrix.exists() and genes.exists() and barcodes.exists() and samples.exists():
            rows.append(
                {
                    "dataset_id": f"GSE125449_{set_id}",
                    "source": "local:G:/wanyi_HCC_scRNA/HCCscRNA",
                    "modality": "scRNA-seq",
                    "input_format": "10x_mtx",
                    "expression_path": str(matrix),
                    "metadata_path": str(samples),
                    "sample_id": f"GSE125449_{set_id}",
                    "status": "ready",
                }
            )
    return pd.DataFrame(rows, columns=columns)


def infer_comparison_group(metadata: pd.DataFrame) -> pd.Series:
    groups = pd.Series("reference", index=metadata.index, dtype=object)
    if "NormalvsTumor" in metadata.columns:
        values = metadata["NormalvsTumor"].astype(object).where(metadata["NormalvsTumor"].notna(), "").astype(str).str.upper().str.strip()
        groups.loc[values.eq("T")] = "malignant"
        groups.loc[values.eq("N")] = "reference"
        decided = values.isin(["T", "N"])
    else:
        decided = pd.Series(False, index=metadata.index)

    undecided = ~decided
    if "Type" in metadata.columns:
        values = metadata["Type"].astype(object).where(metadata["Type"].notna(), "").astype(str).str.lower()
        groups.loc[undecided & values.str.contains("malignant|tumou?r|hcc")] = "malignant"
        reference_labels = values.isin({"t cell", "b cell", "caf", "tam", "tec", "nk", "unclassified", "hpc-like"})
        groups.loc[undecided & reference_labels] = "reference"
        decided = decided | values.str.contains("malignant|tumou?r|hcc") | reference_labels

    undecided = ~decided
    if "cell_type" in metadata.columns:
        values = metadata["cell_type"].astype(object).where(metadata["cell_type"].notna(), "").astype(str).str.lower()
        groups.loc[undecided & values.str.contains("tumou?r|malignant")] = "malignant"
        groups.loc[undecided & values.str.contains("tcell|bcell|nk|mye|hsc|caf|tec|endo|macro")] = "reference"
        decided = decided | values.str.contains("tumou?r|malignant|tcell|bcell|nk|mye|hsc|caf|tec|endo|macro")

    undecided = ~decided
    if "tissue_source" in metadata.columns:
        values = metadata["tissue_source"].astype(object).where(metadata["tissue_source"].notna(), "").astype(str).str.lower()
        groups.loc[undecided & values.str.contains("tumou?r")] = "malignant"
        groups.loc[undecided & values.str.contains("adjacent|normal|liver")] = "reference"
    return groups


def expression_matrix_from_h5ad(path: str | Path, signatures: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    import anndata as ad
    from scipy import sparse

    adata = ad.read_h5ad(path, backed="r")
    try:
        var_names = pd.Index(normalize_gene_symbols(adata.var_names))
        wanted = set(normalize_gene_symbols(signatures["gene"]))
        mask = var_names.isin(wanted)
        if not mask.any():
            return pd.DataFrame(), pd.DataFrame()
        sub = adata[:, mask].to_memory()
        genes = var_names[mask].tolist()
        x = sub.X
        if sparse.issparse(x):
            x = x.toarray()
        expr = pd.DataFrame(np.asarray(x).T, index=genes, columns=sub.obs_names.astype(str))
        metadata = sub.obs.copy()
        metadata["cell_id"] = sub.obs_names.astype(str)
        metadata["comparison_group"] = infer_comparison_group(metadata)
        return expr, metadata[["cell_id", "comparison_group"]]
    finally:
        if getattr(adata, "isbacked", False):
            adata.file.close()


def expression_matrix_from_gene_by_cell_tsv(
    expression_path: str | Path,
    metadata_path: str | Path,
    signatures: pd.DataFrame,
    chunksize: int = 500,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    wanted = set(normalize_gene_symbols(signatures["gene"]))
    selected = []
    for chunk in pd.read_csv(expression_path, sep="\t", chunksize=chunksize):
        gene_col = chunk.columns[0]
        chunk["_gene_symbol"] = normalize_gene_symbols(chunk[gene_col])
        keep = chunk["_gene_symbol"].isin(wanted)
        if keep.any():
            kept = chunk.loc[keep].drop(columns=[gene_col]).set_index("_gene_symbol")
            selected.append(kept)
    expr = pd.concat(selected, axis=0) if selected else pd.DataFrame()
    if not expr.empty:
        expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0)
        expr = expr.groupby(expr.index).mean()
    metadata = pd.read_csv(metadata_path, sep="\t", skiprows=[1])
    metadata["cell_id"] = metadata["name"].astype(str)
    metadata["comparison_group"] = infer_comparison_group(metadata)
    common = [cell for cell in expr.columns.astype(str) if cell in set(metadata["cell_id"])] if not expr.empty else []
    if common:
        expr = expr.loc[:, common]
        metadata = metadata.loc[metadata["cell_id"].isin(common), ["cell_id", "comparison_group"]]
    return expr, metadata


def expression_matrix_from_10x_source(source_row: pd.Series, signatures: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    from scipy import sparse
    from scipy.io import mmread

    matrix_path = Path(source_row["expression_path"])
    base_dir = matrix_path.parent
    dataset_prefix = matrix_path.name.replace("_matrix.mtx.gz", "")
    genes_path = base_dir / f"{dataset_prefix}_genes.tsv.gz"
    barcodes_path = base_dir / f"{dataset_prefix}_barcodes.tsv.gz"
    metadata_path = Path(source_row["metadata_path"])
    genes = pd.read_csv(genes_path, sep="\t", header=None)
    gene_symbols = pd.Index(normalize_gene_symbols(genes.iloc[:, 1] if genes.shape[1] > 1 else genes.iloc[:, 0]))
    wanted = set(normalize_gene_symbols(signatures["gene"]))
    mask = gene_symbols.isin(wanted)
    if not mask.any():
        return pd.DataFrame(), pd.DataFrame()
    barcodes = pd.read_csv(barcodes_path, sep="\t", header=None).iloc[:, 0].astype(str).tolist()
    mat = mmread(matrix_path).tocsr()
    sub = mat[np.asarray(mask), :]
    if sparse.issparse(sub):
        x = sub.toarray()
    else:
        x = np.asarray(sub)
    expr = pd.DataFrame(x, index=gene_symbols[mask].tolist(), columns=barcodes).groupby(level=0).mean()
    metadata = pd.read_csv(metadata_path, sep="\t")
    metadata["cell_id"] = metadata["Cell Barcode"].astype(str)
    metadata["comparison_group"] = infer_comparison_group(metadata)
    common = [cell for cell in expr.columns.astype(str) if cell in set(metadata["cell_id"])]
    expr = expr.loc[:, common]
    metadata = metadata.loc[metadata["cell_id"].isin(common), ["cell_id", "comparison_group"]]
    return expr, metadata


def score_local_scrna_sources(sources: pd.DataFrame, signatures: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    score_summaries = []
    recurrence_rows = []
    detections = []
    if sources.empty or signatures.empty:
        return empty_scrna_scores_schema(), pd.DataFrame(), pd.DataFrame()
    for _, source in sources.iterrows():
        try:
            if source["input_format"] == "h5ad":
                expr, metadata = expression_matrix_from_h5ad(source["expression_path"], signatures)
            elif source["input_format"] == "gene_by_cell_tsv_gz":
                expr, metadata = expression_matrix_from_gene_by_cell_tsv(
                    source["expression_path"], source["metadata_path"], signatures
                )
            elif source["input_format"] == "10x_mtx":
                expr, metadata = expression_matrix_from_10x_source(source, signatures)
            else:
                continue
            scores, detection = compute_signature_scores(expr, signatures, dataset_id=source["dataset_id"])
            if scores.empty:
                continue
            axis_scores = collapse_scores_to_axis(scores)
            recurrence = compute_group_recurrence(
                axis_scores,
                metadata,
                dataset_id=source["dataset_id"],
                modality="scRNA-seq",
                positive_group="malignant",
                reference_group="reference",
            )
            summary = summarize_scrna_dataset_scores(scores, metadata, detection, recurrence)
            score_summaries.append(summary)
            recurrence_rows.append(recurrence)
            detections.append(detection)
        except Exception as exc:
            score_summaries.append(
                pd.DataFrame(
                    [
                        {
                            "dataset_id": source["dataset_id"],
                            "modality": "scRNA-seq",
                            "axis": "",
                            "tf": "",
                            "cell_state": "",
                            "n_cells": 0,
                            "mean_signature_score": np.nan,
                            "gene_detection_rate": np.nan,
                            "effect_size": np.nan,
                            "pvalue": 1.0,
                            "p.adjust": 1.0,
                            "direction": "not_testable",
                            "missingness_rate": 1.0,
                            "status": f"read_error:{type(exc).__name__}",
                        }
                    ]
                )
            )
    summary = pd.concat(score_summaries, ignore_index=True) if score_summaries else empty_scrna_scores_schema()
    recurrence = pd.concat(recurrence_rows, ignore_index=True) if recurrence_rows else pd.DataFrame()
    detection = pd.concat(detections, ignore_index=True) if detections else pd.DataFrame()
    return summary, recurrence, detection


def parse_gtf_ensembl_to_symbol(gtf_path: str | Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    path = Path(gtf_path)
    if not path.exists():
        return mapping
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = fields[8]
            gene_id_match = re.search(r'gene_id "([^"]+)"', attrs)
            gene_name_match = re.search(r'gene_name "([^"]+)"', attrs)
            if gene_id_match and gene_name_match:
                mapping[strip_ensembl_version(gene_id_match.group(1))] = normalize_gene_symbols([gene_name_match.group(1)])[0]
    return mapping


def load_tcga_expression_signature_genes(
    expression_path: str | Path,
    signatures: pd.DataFrame,
    ensembl_to_symbol: dict[str, str],
    chunksize: int = 1000,
) -> pd.DataFrame:
    wanted = set(normalize_gene_symbols(signatures["gene"]))
    selected = []
    for chunk in pd.read_csv(expression_path, sep="\t", chunksize=chunksize):
        gene_col = chunk.columns[0]
        chunk["_ensembl_base"] = chunk[gene_col].map(strip_ensembl_version)
        chunk["_gene_symbol"] = chunk["_ensembl_base"].map(ensembl_to_symbol).fillna("")
        keep = chunk["_gene_symbol"].isin(wanted)
        if keep.any():
            kept = chunk.loc[keep].drop(columns=[gene_col, "_ensembl_base"]).set_index("_gene_symbol")
            selected.append(kept)
    if not selected:
        return pd.DataFrame()
    expr = pd.concat(selected, axis=0)
    expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0)
    return expr.groupby(expr.index).mean()


def load_symbol_expression_signature_genes(
    expression_path: str | Path,
    signatures: pd.DataFrame,
    chunksize: int = 1000,
) -> pd.DataFrame:
    wanted = set(normalize_gene_symbols(signatures["gene"]))
    selected = []
    for chunk in pd.read_csv(expression_path, sep="\t", chunksize=chunksize):
        gene_col = chunk.columns[0]
        chunk["_gene_symbol"] = normalize_gene_symbols(chunk[gene_col])
        keep = chunk["_gene_symbol"].isin(wanted)
        if keep.any():
            kept = chunk.loc[keep].drop(columns=[gene_col]).set_index("_gene_symbol")
            selected.append(kept)
    if not selected:
        return pd.DataFrame()
    expr = pd.concat(selected, axis=0)
    expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0)
    return expr.groupby(expr.index).mean()


def extract_discovery_dataset_ids(metadata_dir: Path) -> set[str]:
    dataset_ids: set[str] = set()
    for path in metadata_dir.glob("celloracle_module6*.tsv"):
        text = path.read_text(encoding="utf-8", errors="ignore")
        dataset_ids.update(re.findall(r"GSE\d+", text))
    for path in metadata_dir.glob("sctenifoldknk_module7*.tsv"):
        text = path.read_text(encoding="utf-8", errors="ignore")
        dataset_ids.update(re.findall(r"GSE\d+", text))
    return dataset_ids


def build_external_dataset_manifest(discovery_dataset_ids: set[str] | None = None) -> pd.DataFrame:
    discovery_dataset_ids = discovery_dataset_ids or set()
    rows = [
        {
            "dataset_id": "GSE156625",
            "source": "NCBI GEO",
            "modality": "scRNA-seq",
            "cancer_type": "hepatocellular carcinoma",
            "species": "Homo sapiens",
            "n_samples": np.nan,
            "n_cells": np.nan,
            "tissue_type": "tumor/adjacent liver candidate",
            "access_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156625",
            "raw_available": True,
            "processed_available": True,
            "included": True,
            "exclusion_reason": "",
        },
        {
            "dataset_id": "GSE98638",
            "source": "NCBI GEO",
            "modality": "scRNA-seq",
            "cancer_type": "liver cancer / hepatocellular carcinoma",
            "species": "Homo sapiens",
            "n_samples": np.nan,
            "n_cells": np.nan,
            "tissue_type": "tumor immune microenvironment candidate",
            "access_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98638",
            "raw_available": True,
            "processed_available": True,
            "included": True,
            "exclusion_reason": "",
        },
        {
            "dataset_id": "GSE151530",
            "source": "NCBI GEO",
            "modality": "scRNA-seq",
            "cancer_type": "hepatocellular carcinoma",
            "species": "Homo sapiens",
            "n_samples": np.nan,
            "n_cells": np.nan,
            "tissue_type": "discovery leakage audit control",
            "access_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151530",
            "raw_available": True,
            "processed_available": True,
            "included": True,
            "exclusion_reason": "",
        },
        {
            "dataset_id": "GSE149614",
            "source": "NCBI GEO",
            "modality": "scRNA-seq",
            "cancer_type": "hepatocellular carcinoma",
            "species": "Homo sapiens",
            "n_samples": np.nan,
            "n_cells": np.nan,
            "tissue_type": "discovery leakage audit control",
            "access_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614",
            "raw_available": True,
            "processed_available": True,
            "included": True,
            "exclusion_reason": "",
        },
        {
            "dataset_id": "TCGA-LIHC",
            "source": "GDC/TCGA",
            "modality": "bulk_RNA-seq",
            "cancer_type": "liver hepatocellular carcinoma",
            "species": "Homo sapiens",
            "n_samples": np.nan,
            "n_cells": np.nan,
            "tissue_type": "tumor/solid tissue normal",
            "access_url": "https://portal.gdc.cancer.gov/projects/TCGA-LIHC",
            "raw_available": True,
            "processed_available": True,
            "included": True,
            "exclusion_reason": "",
        },
        {
            "dataset_id": "LIRI-JP",
            "source": "ICGC",
            "modality": "bulk_RNA-seq",
            "cancer_type": "liver cancer",
            "species": "Homo sapiens",
            "n_samples": np.nan,
            "n_cells": np.nan,
            "tissue_type": "tumor cohort",
            "access_url": "https://dcc.icgc.org/projects/LIRI-JP",
            "raw_available": True,
            "processed_available": True,
            "included": True,
            "exclusion_reason": "",
        },
        {
            "dataset_id": "PUBLIC_HCC_SPATIAL_TBD",
            "source": "public spatial transcriptomics repository",
            "modality": "spatial_transcriptomics",
            "cancer_type": "hepatocellular carcinoma",
            "species": "Homo sapiens",
            "n_samples": np.nan,
            "n_cells": np.nan,
            "tissue_type": "tumor spatial candidate",
            "access_url": "",
            "raw_available": False,
            "processed_available": False,
            "included": False,
            "exclusion_reason": "pending_public_dataset_selection",
        },
    ]
    manifest = pd.DataFrame(rows, columns=MODULE8_REQUIRED_MANIFEST_COLUMNS)
    return audit_dataset_leakage(manifest, discovery_dataset_ids)


def audit_dataset_leakage(manifest: pd.DataFrame, discovery_dataset_ids: set[str]) -> pd.DataFrame:
    audited = manifest.copy()
    audited["dataset_id"] = audited["dataset_id"].astype(str)
    leaked = audited["dataset_id"].isin(discovery_dataset_ids)
    audited.loc[leaked, "included"] = False
    audited.loc[leaked, "exclusion_reason"] = audited.loc[leaked, "exclusion_reason"].apply(
        lambda value: "discovery_or_lodo_dataset"
        if not str(value).strip() or str(value).lower() == "nan"
        else str(value)
        if "discovery_or_lodo_dataset" in str(value)
        else f"{value};discovery_or_lodo_dataset"
    )
    audited["dataset_leakage_status"] = np.where(leaked, "leaked_excluded", "clean")
    return audited


def axis_from_display_group(display_group: str) -> tuple[str, str]:
    group = str(display_group)
    if group == "main_tier1":
        return "tier1_rescue", "main"
    if group == "main_ap1_axis":
        return "ap1_stress_proliferation", "main"
    if group == "main_state_specific":
        return "sox4_state_specific", "main"
    return "control_calibration", "control"


def build_signature_registry(tf_matrix: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in tf_matrix.iterrows():
        axis, axis_group = axis_from_display_group(row.get("display_group", ""))
        rows.append(
            {
                "tf": str(row["tf"]),
                "axis": axis,
                "axis_group": axis_group,
                "signature_class": "main" if axis_group == "main" else "control",
                "display_group": row.get("display_group", ""),
                "candidate_tier": row.get("candidate_tier", ""),
                "module7_integrated_rank": row.get("module7_integrated_rank", np.nan),
                "module7_integrated_score": row.get("integrated_module7_score", np.nan),
                "signature_direction": expected_axis_direction(axis),
            }
        )
    return pd.DataFrame(rows).sort_values(["axis_group", "axis", "module7_integrated_rank", "tf"]).reset_index(drop=True)


def expected_axis_direction(axis: str) -> str:
    if axis == "tier1_rescue":
        return "reference_or_differentiated_high"
    if axis == "ap1_stress_proliferation":
        return "malignant_or_transition_high"
    if axis == "sox4_state_specific":
        return "malignant_like_state_high"
    return "calibration_no_systematic_outperformance"


def build_tf_target_signature_genes(registry: pd.DataFrame, perturbation: pd.DataFrame, top_n: int = 50) -> pd.DataFrame:
    if registry.empty or perturbation.empty:
        return pd.DataFrame(
            columns=["axis", "tf", "gene", "signature_source", "rank", "distance", "p.adj", "signature_class"]
        )
    reg = registry[["tf", "axis", "signature_class"]].drop_duplicates()
    pert = perturbation.copy()
    pert["tf"] = pert["tf"].astype(str)
    pert["gene"] = normalize_gene_symbols(pert["gene"])
    pert["distance_abs"] = pd.to_numeric(pert.get("distance", 0), errors="coerce").abs().fillna(0)
    pert["p.adj"] = pd.to_numeric(pert.get("p.adj", 1), errors="coerce").fillna(1)
    rows = []
    for _, tf_row in reg.iterrows():
        tf_pert = pert.loc[pert["tf"].eq(tf_row["tf"])].sort_values(["p.adj", "distance_abs"], ascending=[True, False])
        significant = tf_pert.loc[tf_pert["p.adj"].le(0.05)]
        selected = significant.head(top_n) if not significant.empty else tf_pert.head(top_n)
        for rank, (_, gene_row) in enumerate(selected.iterrows(), start=1):
            rows.append(
                {
                    "axis": tf_row["axis"],
                    "tf": tf_row["tf"],
                    "gene": gene_row["gene"],
                    "signature_source": "sctenifold_top_disturbed",
                    "rank": rank,
                    "distance": gene_row.get("distance", np.nan),
                    "p.adj": gene_row.get("p.adj", np.nan),
                    "signature_class": tf_row["signature_class"],
                }
            )
    return pd.DataFrame(rows)


def build_pathway_signature_genes(
    registry: pd.DataFrame, pathway_matrix: pd.DataFrame, top_terms_per_tf_database: int = 3
) -> pd.DataFrame:
    columns = ["axis", "tf", "database", "term_name", "gene", "signature_source", "term_rank"]
    if registry.empty or pathway_matrix.empty or "overlap_genes" not in pathway_matrix.columns:
        return pd.DataFrame(columns=columns)
    reg = registry[["tf", "axis"]].drop_duplicates()
    pathway = pathway_matrix.merge(reg, on="tf", how="inner")
    pathway["p.adjust_sort"] = pd.to_numeric(pathway.get("p.adjust", 1), errors="coerce").fillna(1)
    pathway["term_rank"] = pathway.sort_values(["p.adjust_sort", "tf"]).groupby(["tf", "database"]).cumcount() + 1
    pathway = pathway.loc[pathway["term_rank"].le(top_terms_per_tf_database)]
    rows = []
    for _, row in pathway.iterrows():
        genes = [gene for gene in str(row.get("overlap_genes", "")).split(";") if gene and gene.lower() != "nan"]
        for gene in normalize_gene_symbols(genes):
            rows.append(
                {
                    "axis": row["axis"],
                    "tf": row["tf"],
                    "database": row.get("database", ""),
                    "term_name": row.get("term_name", ""),
                    "gene": gene,
                    "signature_source": "module7_top_pathway_overlap",
                    "term_rank": row["term_rank"],
                }
            )
    return pd.DataFrame(rows, columns=columns).drop_duplicates()


def _prepare_expression_by_gene(expression: pd.DataFrame) -> pd.DataFrame:
    expr = expression.copy()
    expr.index = normalize_gene_symbols(expr.index)
    expr = expr.loc[expr.index.astype(bool)]
    return expr.groupby(expr.index).mean()


def compute_signature_scores(
    expression: pd.DataFrame, signatures: pd.DataFrame, dataset_id: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = _prepare_expression_by_gene(expression)
    if expr.empty or signatures.empty:
        return pd.DataFrame(), pd.DataFrame()
    gene_mean = expr.mean(axis=1)
    gene_std = expr.std(axis=1, ddof=0).replace(0, np.nan)
    z = expr.sub(gene_mean, axis=0).div(gene_std, axis=0).fillna(0)
    rows = []
    detection_rows = []
    for (axis, tf), sig in signatures.groupby(["axis", "tf"], dropna=False):
        genes = [gene for gene in normalize_gene_symbols(sig["gene"]) if gene]
        present = sorted(set(genes).intersection(z.index))
        detection_rate = len(present) / len(set(genes)) if genes else 0.0
        detection_rows.append(
            {
                "dataset_id": dataset_id,
                "axis": axis,
                "tf": tf,
                "n_signature_genes": len(set(genes)),
                "n_detected_genes": len(present),
                "gene_detection_rate": detection_rate,
            }
        )
        if not present:
            continue
        scores = z.loc[present].mean(axis=0)
        for cell_id, score in scores.items():
            rows.append(
                {
                    "dataset_id": dataset_id,
                    "axis": axis,
                    "tf": tf,
                    "cell_id": cell_id,
                    "signature_score": float(score),
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(detection_rows)


def compute_bulk_signature_scores(
    expression: pd.DataFrame, signatures: pd.DataFrame, dataset_id: str, sample_type_mode: str = "tcga"
) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = _prepare_expression_by_gene(expression)
    if expr.empty or signatures.empty:
        return pd.DataFrame(), pd.DataFrame()
    gene_mean = expr.mean(axis=1)
    gene_std = expr.std(axis=1, ddof=0).replace(0, np.nan)
    z = expr.sub(gene_mean, axis=0).div(gene_std, axis=0).fillna(0)
    axis_signatures = signatures[["axis", "gene"]].drop_duplicates()
    rows = []
    detection_rows = []
    if sample_type_mode == "icgc":
        inferred_types = infer_icgc_sample_type(expr.columns)
    elif sample_type_mode == "tcga":
        inferred_types = infer_tcga_sample_type(expr.columns)
    else:
        inferred_types = ["unknown"] * len(expr.columns)
    sample_types = dict(zip(expr.columns.astype(str), inferred_types))
    for axis, sig in axis_signatures.groupby("axis", dropna=False):
        genes = [gene for gene in normalize_gene_symbols(sig["gene"]) if gene]
        present = sorted(set(genes).intersection(z.index))
        detection_rate = len(present) / len(set(genes)) if genes else 0.0
        detection_rows.append(
            {
                "dataset_id": dataset_id,
                "axis": axis,
                "n_signature_genes": len(set(genes)),
                "n_detected_genes": len(present),
                "gene_detection_rate": detection_rate,
            }
        )
        if not present:
            continue
        scores = z.loc[present].mean(axis=0)
        for sample, score in scores.items():
            rows.append(
                {
                    "dataset_id": dataset_id,
                    "axis": axis,
                    "sample": str(sample),
                    "signature_score": float(score),
                    "sample_type": sample_types.get(str(sample), "unknown"),
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(detection_rows)


def compute_bulk_tumor_normal_association(scores: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "dataset_id",
        "modality",
        "axis",
        "clinical_variable",
        "n_samples",
        "effect_size",
        "hazard_ratio",
        "pvalue",
        "p.adjust",
        "direction",
        "missingness_rate",
        "status",
    ]
    if scores.empty:
        return pd.DataFrame(columns=columns)
    rows = []
    missingness = float(scores["sample_type"].eq("unknown").mean()) if "sample_type" in scores.columns else 1.0
    for (dataset_id, axis), axis_df in scores.groupby(["dataset_id", "axis"], dropna=False):
        tumor = axis_df.loc[axis_df["sample_type"].eq("tumor"), "signature_score"].dropna().to_numpy()
        normal = axis_df.loc[axis_df["sample_type"].eq("normal"), "signature_score"].dropna().to_numpy()
        if len(tumor) and len(normal):
            effect = float(np.mean(tumor) - np.mean(normal))
            pvalue = mann_whitney_pvalue(tumor, normal)
            direction = "positive" if effect > 0 else "negative" if effect < 0 else "no_effect"
            status = "tested"
        else:
            effect = np.nan
            pvalue = 1.0
            direction = "not_testable"
            status = "insufficient_tumor_or_normal_samples"
        rows.append(
            {
                "dataset_id": dataset_id,
                "modality": "bulk_RNA-seq",
                "axis": axis,
                "clinical_variable": "tumor_vs_normal",
                "n_samples": int(len(tumor) + len(normal)),
                "effect_size": effect,
                "hazard_ratio": np.nan,
                "pvalue": pvalue,
                "direction": direction,
                "missingness_rate": missingness,
                "status": status,
            }
        )
    out = pd.DataFrame(rows)
    out["p.adjust"] = benjamini_hochberg(out["pvalue"])
    return out[columns]


def compute_exploratory_survival_association(scores: pd.DataFrame, survival: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "dataset_id",
        "modality",
        "axis",
        "clinical_variable",
        "n_samples",
        "effect_size",
        "hazard_ratio",
        "pvalue",
        "p.adjust",
        "direction",
        "missingness_rate",
        "status",
    ]
    if scores.empty or survival.empty or not {"sample", "OS", "OS.time"}.issubset(survival.columns):
        return pd.DataFrame(columns=columns)
    surv = survival[["sample", "OS", "OS.time"]].copy()
    surv["sample"] = surv["sample"].astype(str)
    rows = []
    for (dataset_id, axis), axis_df in scores.groupby(["dataset_id", "axis"], dropna=False):
        merged = axis_df[["sample", "signature_score"]].merge(surv, on="sample", how="inner")
        merged = merged.dropna(subset=["signature_score", "OS", "OS.time"])
        if len(merged) >= 6 and merged["OS"].nunique() > 1:
            event = pd.to_numeric(merged["OS"], errors="coerce").fillna(0)
            score = pd.to_numeric(merged["signature_score"], errors="coerce").fillna(0)
            time = pd.to_numeric(merged["OS.time"], errors="coerce").replace(0, np.nan)
            event_scores = score.loc[event.astype(bool)]
            censored_scores = score.loc[~event.astype(bool)]
            effect = float(event_scores.mean() - censored_scores.mean()) if len(event_scores) and len(censored_scores) else np.nan
            pvalue = (
                mann_whitney_pvalue(event_scores.to_numpy(), censored_scores.to_numpy())
                if len(event_scores) and len(censored_scores)
                else 1.0
            )
            hazard_ratio = float(np.exp(effect)) if np.isfinite(effect) else np.nan
            direction = "positive" if effect > 0 else "negative" if effect < 0 else "no_effect"
            missingness = float(time.isna().mean())
        else:
            effect = np.nan
            pvalue = 1.0
            hazard_ratio = np.nan
            direction = "not_testable"
            missingness = 1.0
        rows.append(
            {
                "dataset_id": dataset_id,
                "modality": "bulk_RNA-seq",
                "axis": axis,
                "clinical_variable": "OS_univariate_exploratory",
                "n_samples": int(len(merged)),
                "effect_size": effect,
                "hazard_ratio": hazard_ratio,
                "pvalue": pvalue,
                "direction": direction,
                "missingness_rate": missingness,
                "status": "missing_covariates_exploratory",
            }
        )
    out = pd.DataFrame(rows)
    out["p.adjust"] = benjamini_hochberg(out["pvalue"])
    return out[columns]


def survival_summary_columns() -> list[str]:
    return [
        "dataset_id",
        "modality",
        "axis",
        "clinical_variable",
        "n_samples",
        "effect_size",
        "hazard_ratio",
        "pvalue",
        "p.adjust",
        "direction",
        "missingness_rate",
        "status",
        "covariates_used",
        "model_type",
    ]


def compute_adjusted_survival_association(
    scores: pd.DataFrame,
    survival: pd.DataFrame,
    clinical: pd.DataFrame,
    dataset_id: str,
    survival_id_col: str,
    event_col: str,
    time_col: str,
    clinical_id_col: str,
    covariates: list[str],
) -> pd.DataFrame:
    columns = survival_summary_columns()
    if scores.empty or survival.empty or clinical.empty:
        return pd.DataFrame(columns=columns)
    surv = survival[[survival_id_col, event_col, time_col]].copy()
    surv = surv.rename(columns={survival_id_col: "patient_id", event_col: "event", time_col: "time"})
    clin_cols = [clinical_id_col] + [cov for cov in covariates if cov in clinical.columns]
    clin = clinical[clin_cols].copy().rename(columns={clinical_id_col: "patient_id"})
    usable_covariates = [cov for cov in covariates if cov in clin.columns]
    rows = []
    for axis, axis_df in scores.groupby("axis", dropna=False):
        analysis_df = axis_df
        if "sample_type" in axis_df.columns and axis_df["sample_type"].astype(str).eq("tumor").any():
            analysis_df = axis_df.loc[axis_df["sample_type"].astype(str).eq("tumor")]
        merged = analysis_df[["sample", "signature_score"]].copy()
        merged["patient_id"] = merged["sample"].map(sample_to_patient_id)
        merged = merged.merge(surv, on="patient_id", how="inner").merge(clin, on="patient_id", how="left")
        merged = merged.dropna(subset=["signature_score", "event", "time"])
        n_samples = int(len(merged))
        covariates_present = [
            cov
            for cov in usable_covariates
            if cov in merged.columns and pd.to_numeric(merged[cov], errors="coerce").notna().sum() >= max(4, n_samples // 4)
        ]
        missingness = float(merged[covariates_present].isna().mean().mean()) if covariates_present else 1.0
        if n_samples >= 8 and merged["event"].nunique() > 1 and covariates_present:
            try:
                from statsmodels.duration.hazard_regression import PHReg

                model_df = merged[["time", "event", "signature_score"] + covariates_present].copy()
                for col in ["time", "event", "signature_score"] + covariates_present:
                    model_df[col] = pd.to_numeric(model_df[col], errors="coerce")
                model_df = model_df.dropna()
                exog = model_df[["signature_score"] + covariates_present]
                result = PHReg(model_df["time"], exog, status=model_df["event"]).fit(disp=False)
                coef = float(result.params[0])
                pvalue = float(result.pvalues[0])
                hazard_ratio = float(np.exp(coef))
                direction = "positive" if coef > 0 else "negative" if coef < 0 else "no_effect"
                status = "adjusted_cox"
                effect = coef
                model_type = "cox_phreg"
                n_samples = int(len(model_df))
            except Exception:
                event_scores = pd.to_numeric(merged.loc[merged["event"].astype(bool), "signature_score"], errors="coerce").dropna()
                censored_scores = pd.to_numeric(merged.loc[~merged["event"].astype(bool), "signature_score"], errors="coerce").dropna()
                effect = float(event_scores.mean() - censored_scores.mean()) if len(event_scores) and len(censored_scores) else np.nan
                pvalue = mann_whitney_pvalue(event_scores.to_numpy(), censored_scores.to_numpy()) if len(event_scores) and len(censored_scores) else 1.0
                hazard_ratio = float(np.exp(effect)) if np.isfinite(effect) else np.nan
                direction = "positive" if effect > 0 else "negative" if effect < 0 else "not_testable"
                status = "adjusted_cox_failed_exploratory"
                model_type = "event_score_mannwhitney"
        else:
            event_scores = pd.to_numeric(merged.loc[merged["event"].astype(bool), "signature_score"], errors="coerce").dropna()
            censored_scores = pd.to_numeric(merged.loc[~merged["event"].astype(bool), "signature_score"], errors="coerce").dropna()
            effect = float(event_scores.mean() - censored_scores.mean()) if len(event_scores) and len(censored_scores) else np.nan
            pvalue = mann_whitney_pvalue(event_scores.to_numpy(), censored_scores.to_numpy()) if len(event_scores) and len(censored_scores) else 1.0
            hazard_ratio = float(np.exp(effect)) if np.isfinite(effect) else np.nan
            direction = "positive" if effect > 0 else "negative" if effect < 0 else "not_testable"
            status = "missing_covariates_exploratory"
            model_type = "event_score_mannwhitney"
        rows.append(
            {
                "dataset_id": dataset_id,
                "modality": "bulk_RNA-seq",
                "axis": axis,
                "clinical_variable": "OS_adjusted",
                "n_samples": n_samples,
                "effect_size": effect,
                "hazard_ratio": hazard_ratio,
                "pvalue": pvalue,
                "direction": direction,
                "missingness_rate": missingness,
                "status": status,
                "covariates_used": ";".join(covariates_present),
                "model_type": model_type,
            }
        )
    out = pd.DataFrame(rows)
    out["p.adjust"] = benjamini_hochberg(out["pvalue"])
    return out[columns]


def compute_bulk_clinical_variable_association(
    scores: pd.DataFrame,
    clinical: pd.DataFrame,
    dataset_id: str,
    variables: list[str],
    clinical_id_col: str = "Id",
) -> pd.DataFrame:
    columns = [
        "dataset_id",
        "modality",
        "axis",
        "clinical_variable",
        "n_samples",
        "effect_size",
        "hazard_ratio",
        "pvalue",
        "p.adjust",
        "direction",
        "missingness_rate",
        "status",
    ]
    if scores.empty or clinical.empty:
        return pd.DataFrame(columns=columns)
    clin = clinical.copy().rename(columns={clinical_id_col: "patient_id"})
    rows = []
    for variable in variables:
        if variable not in clin.columns:
            continue
        for (score_dataset, axis), axis_df in scores.groupby(["dataset_id", "axis"], dropna=False):
            analysis_df = axis_df
            if "sample_type" in axis_df.columns and axis_df["sample_type"].astype(str).eq("tumor").any():
                analysis_df = axis_df.loc[axis_df["sample_type"].astype(str).eq("tumor")]
            merged = analysis_df[["sample", "signature_score"]].copy()
            merged["patient_id"] = merged["sample"].map(sample_to_patient_id)
            merged = merged.merge(clin[["patient_id", variable]], on="patient_id", how="inner")
            merged[variable] = pd.to_numeric(merged[variable], errors="coerce")
            missingness = float(merged[variable].isna().mean()) if len(merged) else 1.0
            merged = merged.dropna(subset=["signature_score", variable])
            if merged.empty:
                continue
            unique_values = sorted(float(value) for value in merged[variable].dropna().unique())
            if len(unique_values) == 2 and set(unique_values).issubset({0.0, 1.0}):
                high = merged.loc[merged[variable].eq(1), "signature_score"].dropna().to_numpy()
                low = merged.loc[merged[variable].eq(0), "signature_score"].dropna().to_numpy()
                comparison_label = f"{variable}_1_vs_0"
            else:
                threshold = merged[variable].median()
                high = merged.loc[merged[variable] > threshold, "signature_score"].dropna().to_numpy()
                low = merged.loc[merged[variable] <= threshold, "signature_score"].dropna().to_numpy()
                comparison_label = f"{variable}_high_vs_low"
            if len(high) and len(low):
                effect = float(np.mean(high) - np.mean(low))
                pvalue = mann_whitney_pvalue(high, low)
                direction = "positive" if effect > 0 else "negative" if effect < 0 else "no_effect"
                status = "tested"
            else:
                effect = np.nan
                pvalue = 1.0
                direction = "not_testable"
                status = "insufficient_high_or_low_samples"
            rows.append(
                {
                    "dataset_id": dataset_id or score_dataset,
                    "modality": "bulk_RNA-seq",
                    "axis": axis,
                    "clinical_variable": comparison_label,
                    "n_samples": int(len(merged)),
                    "effect_size": effect,
                    "hazard_ratio": np.nan,
                    "pvalue": pvalue,
                    "direction": direction,
                    "missingness_rate": missingness,
                    "status": status,
                }
            )
    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(columns=columns)
    out["p.adjust"] = benjamini_hochberg(out["pvalue"])
    return out[columns]


def collapse_scores_to_axis(scores: pd.DataFrame) -> pd.DataFrame:
    if scores.empty:
        return pd.DataFrame(columns=["dataset_id", "axis", "cell_id", "signature_score"])
    return (
        scores.groupby(["dataset_id", "axis", "cell_id"], as_index=False)["signature_score"]
        .mean()
        .sort_values(["dataset_id", "axis", "cell_id"])
        .reset_index(drop=True)
    )


def summarize_scrna_dataset_scores(
    scores: pd.DataFrame,
    metadata: pd.DataFrame,
    detection: pd.DataFrame,
    recurrence: pd.DataFrame,
    modality: str = "scRNA-seq",
) -> pd.DataFrame:
    columns = [
        "dataset_id",
        "modality",
        "axis",
        "tf",
        "cell_state",
        "n_cells",
        "mean_signature_score",
        "gene_detection_rate",
        "effect_size",
        "pvalue",
        "p.adjust",
        "direction",
        "missingness_rate",
        "status",
    ]
    if scores.empty:
        return pd.DataFrame(columns=columns)
    meta = metadata.copy()
    if "cell_id" not in meta.columns:
        meta = meta.reset_index().rename(columns={meta.index.name or "index": "cell_id"})
    if "comparison_group" not in meta.columns:
        meta["comparison_group"] = infer_comparison_group(meta)
    merged = scores.merge(meta[["cell_id", "comparison_group"]], on="cell_id", how="left")
    summary = (
        merged.groupby(["dataset_id", "axis", "tf", "comparison_group"], dropna=False)
        .agg(n_cells=("cell_id", "nunique"), mean_signature_score=("signature_score", "mean"))
        .reset_index()
        .rename(columns={"comparison_group": "cell_state"})
    )
    summary["modality"] = modality
    if not detection.empty:
        summary = summary.merge(
            detection[["dataset_id", "axis", "tf", "gene_detection_rate"]],
            on=["dataset_id", "axis", "tf"],
            how="left",
        )
    else:
        summary["gene_detection_rate"] = np.nan
    recurrence_cols = ["axis", "effect_size", "pvalue", "p.adjust", "direction", "missingness_rate"]
    if not recurrence.empty and set(recurrence_cols).issubset(recurrence.columns):
        summary = summary.merge(recurrence[recurrence_cols].drop_duplicates("axis"), on="axis", how="left")
    else:
        for col in recurrence_cols[1:]:
            summary[col] = np.nan
    summary["status"] = "scored"
    return summary[columns]


def benjamini_hochberg(pvalues: Iterable[float]) -> list[float]:
    pvals = np.asarray([1.0 if pd.isna(p) else float(p) for p in pvalues], dtype=float)
    n = len(pvals)
    if n == 0:
        return []
    order = np.argsort(pvals)
    ranks = np.empty(n, dtype=float)
    ranks[order] = np.arange(1, n + 1)
    adjusted = pvals * n / ranks
    adjusted_ordered = np.minimum.accumulate(adjusted[order][::-1])[::-1]
    result = np.empty(n, dtype=float)
    result[order] = np.clip(adjusted_ordered, 0, 1)
    return result.tolist()


def mann_whitney_pvalue(a: np.ndarray, b: np.ndarray) -> float:
    try:
        from scipy.stats import mannwhitneyu

        return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)
    except Exception:
        return 1.0


def compute_group_recurrence(
    scores: pd.DataFrame,
    metadata: pd.DataFrame,
    dataset_id: str,
    modality: str,
    positive_group: str = "malignant",
    reference_group: str = "reference",
) -> pd.DataFrame:
    columns = [
        "dataset_id",
        "modality",
        "axis",
        "comparison",
        "n_positive",
        "n_reference",
        "effect_size",
        "pvalue",
        "p.adjust",
        "direction",
        "missingness_rate",
        "included",
        "dataset_leakage_status",
    ]
    if scores.empty or metadata.empty:
        return pd.DataFrame(columns=columns)
    merged = scores.merge(metadata[["cell_id", "comparison_group"]], on="cell_id", how="left")
    missingness_rate = float(merged["comparison_group"].isna().mean())
    rows = []
    for axis, axis_df in merged.groupby("axis", dropna=False):
        positive = axis_df.loc[axis_df["comparison_group"].eq(positive_group), "signature_score"].dropna().to_numpy()
        reference = axis_df.loc[axis_df["comparison_group"].eq(reference_group), "signature_score"].dropna().to_numpy()
        if len(positive) == 0 or len(reference) == 0:
            effect_size = np.nan
            pvalue = 1.0
            direction = "not_testable"
        else:
            effect_size = float(np.mean(positive) - np.mean(reference))
            pvalue = mann_whitney_pvalue(positive, reference)
            direction = "positive" if effect_size > 0 else "negative" if effect_size < 0 else "no_effect"
        rows.append(
            {
                "dataset_id": dataset_id,
                "modality": modality,
                "axis": axis,
                "comparison": f"{positive_group}_vs_{reference_group}",
                "n_positive": int(len(positive)),
                "n_reference": int(len(reference)),
                "effect_size": effect_size,
                "pvalue": pvalue,
                "direction": direction,
                "missingness_rate": missingness_rate,
                "included": True,
                "dataset_leakage_status": "clean",
            }
        )
    recurrence = pd.DataFrame(rows, columns=[c for c in columns if c != "p.adjust"])
    recurrence["p.adjust"] = benjamini_hochberg(recurrence["pvalue"])
    return recurrence[columns]


def fisher_meta_pvalue(pvalues: Iterable[float]) -> float:
    vals = [max(min(float(p), 1.0), 1e-300) for p in pvalues if not pd.isna(p)]
    if not vals:
        return 1.0
    try:
        from scipy.stats import chi2

        statistic = -2 * sum(math.log(p) for p in vals)
        return float(chi2.sf(statistic, 2 * len(vals)))
    except Exception:
        return min(vals)


def build_axis_level_evidence_grade(recurrence: pd.DataFrame, known_axes: list[str] | None = None) -> pd.DataFrame:
    axes = sorted(set(known_axes or MODULE8_AXES).union(set(recurrence["axis"]) if "axis" in recurrence.columns else set()))
    rows = []
    for axis in axes:
        axis_df = recurrence.loc[recurrence["axis"].eq(axis)].copy() if not recurrence.empty else pd.DataFrame()
        included = axis_df.loc[axis_df.get("included", True).astype(bool)] if not axis_df.empty else axis_df
        significant = included.loc[
            pd.to_numeric(included.get("p.adjust", 1), errors="coerce").fillna(1).le(0.05)
            & included.get("direction", "").astype(str).isin(["positive", "negative"])
        ] if not included.empty else included
        directions = significant["direction"].astype(str).tolist() if not significant.empty else []
        if directions:
            direction_consistency = max(directions.count("positive"), directions.count("negative")) / len(directions)
            consensus_direction = "positive" if directions.count("positive") >= directions.count("negative") else "negative"
        else:
            direction_consistency = 0.0
            consensus_direction = "not_replicated"
        modalities = set(significant["modality"].astype(str)) if not significant.empty else set()
        has_scrna = "scRNA-seq" in modalities
        has_bulk = "bulk_RNA-seq" in modalities
        has_spatial_or_clinical = bool(modalities.intersection({"spatial_transcriptomics", "clinical"}))
        if has_scrna and has_bulk and has_spatial_or_clinical:
            grade = "A"
        elif has_scrna and has_bulk:
            grade = "B"
        elif len(modalities) == 1:
            grade = "C"
        else:
            grade = "D"
        median_effect = (
            float(pd.to_numeric(significant["effect_size"], errors="coerce").median()) if not significant.empty else 0.0
        )
        n_replicated = int(significant["dataset_id"].nunique()) if not significant.empty else 0
        fdr_meta = fisher_meta_pvalue(significant["p.adjust"]) if not significant.empty else 1.0
        rows.append(
            {
                "axis": axis,
                "axis_group": "control" if axis == "control_calibration" else "main",
                "direction_consistency": direction_consistency,
                "consensus_direction": consensus_direction,
                "effect_size_median": median_effect,
                "fdr_meta": fdr_meta,
                "n_replicated_cohorts": n_replicated,
                "dataset_leakage_status": "clean" if not included.empty else "pending_external_data",
                "evidence_grade": grade,
                "recurrence_score": n_replicated * direction_consistency * abs(median_effect),
            }
        )
    return pd.DataFrame(rows)


def flag_external_control_outperformance(grades: pd.DataFrame) -> pd.DataFrame:
    if grades.empty or not {"axis_group", "recurrence_score"}.issubset(grades.columns):
        return pd.DataFrame(columns=["risk_type", "axis", "risk_detail", "severity"])
    main = grades.loc[grades["axis_group"].eq("main"), "recurrence_score"]
    controls = grades.loc[grades["axis_group"].eq("control")]
    if main.empty or controls.empty:
        return pd.DataFrame(columns=["risk_type", "axis", "risk_detail", "severity"])
    threshold = float(main.max())
    rows = []
    for _, row in controls.iterrows():
        if float(row["recurrence_score"]) > threshold:
            rows.append(
                {
                    "risk_type": "review_risk_external_control_outperformance",
                    "axis": row["axis"],
                    "risk_detail": f"Control recurrence score {row['recurrence_score']:.3f} exceeds strongest main axis {threshold:.3f}.",
                    "severity": "review_attention",
                }
            )
    return pd.DataFrame(rows, columns=["risk_type", "axis", "risk_detail", "severity"])


def empty_scrna_scores_schema() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "dataset_id",
            "modality",
            "axis",
            "tf",
            "cell_state",
            "n_cells",
            "mean_signature_score",
            "gene_detection_rate",
            "effect_size",
            "pvalue",
            "p.adjust",
            "direction",
            "missingness_rate",
            "status",
        ]
    )


def empty_bulk_schema() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "dataset_id",
            "modality",
            "axis",
            "clinical_variable",
            "n_samples",
            "effect_size",
            "hazard_ratio",
            "pvalue",
            "p.adjust",
            "direction",
            "missingness_rate",
            "status",
        ]
    )


def empty_spatial_schema() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "dataset_id",
            "modality",
            "axis",
            "region_label",
            "n_spots",
            "mean_signature_score",
            "effect_size",
            "pvalue",
            "p.adjust",
            "direction",
            "missingness_rate",
            "status",
        ]
    )


def _set_nature_matplotlib_style() -> None:
    import matplotlib as mpl

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 6.5,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "legend.frameon": False,
        }
    )


def _ordered_module8_axes(axis_grade: pd.DataFrame) -> list[str]:
    preferred = ["tier1_rescue", "sox4_state_specific", "ap1_stress_proliferation", "control_calibration"]
    if axis_grade.empty or "axis" not in axis_grade.columns:
        return preferred
    observed = axis_grade["axis"].astype(str).tolist()
    return [axis for axis in preferred if axis in observed] + [axis for axis in observed if axis not in preferred]


def _axis_label(axis: object) -> str:
    labels = {
        "tier1_rescue": "HNF4A/PPARA\nCEBPB/EGR1",
        "sox4_state_specific": "SOX4\nstate-specific",
        "ap1_stress_proliferation": "AP-1 stress\nJUN/FOS/ATF3",
        "control_calibration": "Control /\ncalibration",
    }
    return labels.get(str(axis), str(axis).replace("_", " "))


def _panel_label(ax, label: str) -> None:
    ax.text(-0.08, 1.06, label, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom", ha="left")


def _draw_dataset_flow(
    ax,
    manifest: pd.DataFrame,
    scrna_recurrence: pd.DataFrame | None = None,
    bulk_assoc: pd.DataFrame | None = None,
) -> None:
    ax.axis("off")
    if manifest.empty:
        counts = {"candidate": 0, "included": 0, "scrna": 0, "bulk": 0, "excluded": 0}
    else:
        included = manifest.loc[manifest["included"].astype(bool)] if "included" in manifest.columns else manifest
        excluded = manifest.loc[~manifest["included"].astype(bool)] if "included" in manifest.columns else manifest.iloc[0:0]
        modality = included.get("modality", pd.Series(dtype=str)).astype(str).str.lower()
        scored_scrna = (
            int(scrna_recurrence["dataset_id"].nunique())
            if scrna_recurrence is not None and not scrna_recurrence.empty and "dataset_id" in scrna_recurrence.columns
            else int(modality.str.contains("scrna|single", regex=True).sum())
        )
        scored_bulk = (
            int(bulk_assoc["dataset_id"].nunique())
            if bulk_assoc is not None and not bulk_assoc.empty and "dataset_id" in bulk_assoc.columns
            else int(modality.str.contains("bulk").sum())
        )
        scored_total = scored_scrna + scored_bulk
        counts = {
            "candidate": int(len(manifest)),
            "included": scored_total,
            "scrna": scored_scrna,
            "bulk": scored_bulk,
            "excluded": int(len(excluded)),
        }
    boxes = [
        (0.04, 0.58, "Public candidates", counts["candidate"], "#E8EEF7"),
        (0.36, 0.58, "Leakage-audited\nscored", counts["included"], "#DDEFE5"),
        (0.68, 0.72, "scRNA recurrence", counts["scrna"], "#E9F1E3"),
        (0.68, 0.44, "bulk clinical\nvalidation", counts["bulk"], "#F3E8D8"),
        (0.36, 0.18, "excluded / pending", counts["excluded"], "#F1E4E4"),
    ]
    for x, y, title, count, color in boxes:
        ax.add_patch(
            plt.Rectangle((x, y), 0.24, 0.17, transform=ax.transAxes, facecolor=color, edgecolor="#333333", linewidth=0.6)
        )
        ax.text(x + 0.12, y + 0.105, title, transform=ax.transAxes, ha="center", va="center", fontsize=6.5)
        ax.text(x + 0.12, y + 0.035, f"n = {count}", transform=ax.transAxes, ha="center", va="center", fontsize=7, fontweight="bold")
    arrowprops = dict(arrowstyle="-|>", color="#555555", lw=0.8, shrinkA=2, shrinkB=2)
    ax.annotate("", xy=(0.36, 0.665), xytext=(0.28, 0.665), xycoords=ax.transAxes, arrowprops=arrowprops)
    ax.annotate("", xy=(0.68, 0.795), xytext=(0.60, 0.665), xycoords=ax.transAxes, arrowprops=arrowprops)
    ax.annotate("", xy=(0.68, 0.525), xytext=(0.60, 0.665), xycoords=ax.transAxes, arrowprops=arrowprops)
    ax.annotate("", xy=(0.48, 0.35), xytext=(0.48, 0.58), xycoords=ax.transAxes, arrowprops=arrowprops)
    ax.text(0.04, 0.05, "Scored cohorts after leakage audit; discovery/LODO excluded.", transform=ax.transAxes, fontsize=5.8, color="#555555")
    ax.set_title("External validation cohort map", loc="left", fontsize=7.5, fontweight="bold")


def _draw_scrna_heatmap(ax, scrna_recurrence: pd.DataFrame | None, axes_order: list[str]):
    if scrna_recurrence is None or scrna_recurrence.empty:
        heat = pd.DataFrame(0.0, index=["pending"], columns=axes_order)
    else:
        heat = (
            scrna_recurrence.loc[scrna_recurrence["direction"].astype(str).isin(["positive", "negative"])]
            .pivot_table(index="dataset_id", columns="axis", values="effect_size", aggfunc="median")
            .reindex(columns=axes_order)
            .fillna(0)
        )
    values = heat.to_numpy(dtype=float)
    max_abs = max(float(np.nanmax(np.abs(values))) if values.size else 0.0, 0.1)
    im = ax.imshow(values, aspect="auto", cmap="RdBu_r", vmin=-max_abs, vmax=max_abs)
    ax.set_xticks(range(len(heat.columns)))
    ax.set_xticklabels([_axis_label(axis) for axis in heat.columns], rotation=35, ha="right", fontsize=5.8)
    ax.set_yticks(range(len(heat.index)))
    ax.set_yticklabels(heat.index.astype(str), fontsize=5.8)
    ax.set_title("scRNA signature recurrence", loc="left", fontsize=7.5, fontweight="bold")
    ax.tick_params(length=0)
    return im


def _draw_bulk_forest(ax, survival_summary: pd.DataFrame | None, bulk_assoc: pd.DataFrame | None, axes_order: list[str]) -> None:
    plot_df = pd.DataFrame()
    if survival_summary is not None and not survival_summary.empty and "hazard_ratio" in survival_summary.columns:
        plot_df = survival_summary.loc[
            survival_summary["status"].astype(str).eq("adjusted_cox") & survival_summary["axis"].astype(str).isin(axes_order)
        ].copy()
        plot_df["metric"] = plot_df["hazard_ratio"]
        plot_df["label"] = plot_df["dataset_id"].astype(str) + " | " + plot_df["axis"].map(_axis_label).str.replace("\n", " ", regex=False)
        xlabel = "Adjusted OS hazard ratio"
        ref = 1.0
        xscale = "log"
    if plot_df.empty and bulk_assoc is not None and not bulk_assoc.empty:
        plot_df = bulk_assoc.loc[
            bulk_assoc["clinical_variable"].astype(str).eq("tumor_vs_normal") & bulk_assoc["axis"].astype(str).isin(axes_order)
        ].copy()
        plot_df["metric"] = plot_df["effect_size"]
        plot_df["label"] = plot_df["dataset_id"].astype(str) + " | " + plot_df["axis"].map(_axis_label).str.replace("\n", " ", regex=False)
        xlabel = "Tumor-vs-normal effect size"
        ref = 0.0
        xscale = "linear"
    if plot_df.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, "bulk validation pending", ha="center", va="center", transform=ax.transAxes)
        return
    plot_df = plot_df.sort_values(["dataset_id", "axis"]).reset_index(drop=True)
    y = np.arange(len(plot_df))
    pvals = pd.to_numeric(plot_df.get("p.adjust", pd.Series(np.nan, index=plot_df.index)), errors="coerce").fillna(1.0)
    sizes = np.clip(-np.log10(pvals + 1e-300) * 14, 18, 80)
    colors = plot_df["dataset_id"].astype("category").cat.codes.map({0: "#3A6EA5", 1: "#C07A2C"}).fillna("#555555")
    ax.axvline(ref, color="#777777", lw=0.7, ls="--")
    ax.scatter(plot_df["metric"], y, s=sizes, color=colors, edgecolor="white", linewidth=0.4, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["label"], fontsize=5.5)
    ax.invert_yaxis()
    ax.set_xscale(xscale)
    ax.set_xlabel(xlabel, fontsize=6.5)
    ax.set_title("TCGA/ICGC clinical forest plot", loc="left", fontsize=7.5, fontweight="bold")
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.5)


def _draw_clinical_heatmap(ax, bulk_assoc: pd.DataFrame | None, axes_order: list[str]):
    if bulk_assoc is None or bulk_assoc.empty:
        heat = pd.DataFrame(0.0, index=["pending"], columns=axes_order)
    else:
        clinical = bulk_assoc.loc[
            ~bulk_assoc["clinical_variable"].astype(str).eq("tumor_vs_normal")
            & bulk_assoc["status"].astype(str).eq("tested")
            & bulk_assoc["axis"].astype(str).isin(axes_order)
        ].copy()
        if clinical.empty:
            heat = pd.DataFrame(0.0, index=["pending"], columns=axes_order)
        else:
            clinical["row"] = clinical["dataset_id"].astype(str) + "\n" + clinical["clinical_variable"].astype(str)
            heat = clinical.pivot_table(index="row", columns="axis", values="effect_size", aggfunc="median").reindex(columns=axes_order).fillna(0)
    values = heat.to_numpy(dtype=float)
    max_abs = max(float(np.nanmax(np.abs(values))) if values.size else 0.0, 0.1)
    im = ax.imshow(values, aspect="auto", cmap="PuOr_r", vmin=-max_abs, vmax=max_abs)
    ax.set_xticks(range(len(heat.columns)))
    ax.set_xticklabels([_axis_label(axis) for axis in heat.columns], rotation=35, ha="right", fontsize=5.4)
    ax.set_yticks(range(len(heat.index)))
    ax.set_yticklabels(heat.index.astype(str), fontsize=5.2)
    ax.tick_params(length=0)
    ax.set_title("Clinical association heatmap", loc="left", fontsize=7.5, fontweight="bold")
    return im


def _draw_axis_summary(ax, axis_grade: pd.DataFrame, axes_order: list[str]) -> None:
    if axis_grade.empty:
        df = pd.DataFrame({"axis": axes_order, "recurrence_score": 0, "evidence_grade": "D", "n_replicated_cohorts": 0})
    else:
        df = axis_grade.set_index("axis").reindex(axes_order).reset_index()
    scores = pd.to_numeric(df.get("recurrence_score", pd.Series(0, index=df.index)), errors="coerce").fillna(0)
    grades = df.get("evidence_grade", pd.Series("D", index=df.index)).fillna("D").astype(str)
    cohorts = pd.to_numeric(df.get("n_replicated_cohorts", pd.Series(0, index=df.index)), errors="coerce").fillna(0).astype(int)
    color_map = {"A": "#2E7D4F", "B": "#3A6EA5", "C": "#C07A2C", "D": "#888888"}
    y = np.arange(len(df))
    ax.barh(y, scores, color=[color_map.get(grade, "#888888") for grade in grades], height=0.62)
    for idx, (score, grade, n) in enumerate(zip(scores, grades, cohorts)):
        ax.text(max(float(score) * 0.5, 0.05), idx, f"{grade} | {n} cohorts", ha="center", va="center", color="white", fontsize=5.8, fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels([_axis_label(axis).replace("\n", " ") for axis in df["axis"]], fontsize=5.8)
    ax.invert_yaxis()
    ax.set_xlabel("Recurrence score", fontsize=6.5)
    ax.set_title("Axis-level recurrence summary", loc="left", fontsize=7.5, fontweight="bold")
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.5)


def write_module8_nature_figures(
    manifest: pd.DataFrame,
    axis_grade: pd.DataFrame,
    figure_dir: Path,
    scrna_recurrence: pd.DataFrame | None = None,
    bulk_assoc: pd.DataFrame | None = None,
    survival_summary: pd.DataFrame | None = None,
) -> dict:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    globals()["plt"] = plt
    _set_nature_matplotlib_style()
    figure_dir.mkdir(parents=True, exist_ok=True)
    outputs = {}

    flow_stem = figure_dir / "module8_external_dataset_flow"
    heatmap_stem = figure_dir / "module8_external_scrna_signature_heatmap"
    forest_stem = figure_dir / "module8_tcga_icgc_signature_forestplot"
    clinical_heatmap_stem = figure_dir / "module8_tcga_icgc_clinical_association_heatmap"
    recurrence_stem = figure_dir / "module8_external_axis_recurrence_summary"
    multipanel_stem = figure_dir / "module8_external_validation_figure8"

    axes = _ordered_module8_axes(axis_grade)

    fig, ax = plt.subplots(figsize=(3.6, 2.25), dpi=220)
    _draw_dataset_flow(ax, manifest, scrna_recurrence=scrna_recurrence, bulk_assoc=bulk_assoc)
    fig.tight_layout()
    outputs.update(save_figure_all_formats(fig, flow_stem, "figure8_dataset_flow"))
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.0, 3.0), dpi=220)
    im = _draw_scrna_heatmap(ax, scrna_recurrence, axes)
    fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02, label="Effect size")
    fig.tight_layout()
    outputs.update(save_figure_all_formats(fig, heatmap_stem, "scrna_heatmap"))
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.4, 3.0), dpi=220)
    _draw_bulk_forest(ax, survival_summary, bulk_assoc, axes)
    fig.tight_layout()
    outputs.update(save_figure_all_formats(fig, forest_stem, "bulk_forestplot"))
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.4, 3.4), dpi=220)
    im = _draw_clinical_heatmap(ax, bulk_assoc, axes)
    fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02, label="Effect size")
    fig.tight_layout()
    outputs.update(save_figure_all_formats(fig, clinical_heatmap_stem, "clinical_heatmap"))
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.2, 2.6), dpi=220)
    _draw_axis_summary(ax, axis_grade, axes)
    fig.tight_layout()
    outputs.update(save_figure_all_formats(fig, recurrence_stem, "axis_recurrence_summary"))
    plt.close(fig)

    fig = plt.figure(figsize=(7.2, 7.6), dpi=220)
    gs = fig.add_gridspec(3, 2, height_ratios=[0.95, 1.05, 1.0], width_ratios=[1.0, 1.25], hspace=0.72, wspace=0.55)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c1 = fig.add_subplot(gs[1, 0])
    ax_c2 = fig.add_subplot(gs[1, 1])
    ax_d = fig.add_subplot(gs[2, :])
    _draw_dataset_flow(ax_a, manifest, scrna_recurrence=scrna_recurrence, bulk_assoc=bulk_assoc)
    im_b = _draw_scrna_heatmap(ax_b, scrna_recurrence, axes)
    _draw_bulk_forest(ax_c1, survival_summary, bulk_assoc, axes)
    im_c = _draw_clinical_heatmap(ax_c2, bulk_assoc, axes)
    _draw_axis_summary(ax_d, axis_grade, axes)
    _panel_label(ax_a, "a")
    _panel_label(ax_b, "b")
    _panel_label(ax_c1, "c")
    _panel_label(ax_d, "d")
    fig.colorbar(im_b, ax=ax_b, fraction=0.046, pad=0.02, label="Effect")
    fig.colorbar(im_c, ax=ax_c2, fraction=0.046, pad=0.02, label="Effect")
    outputs.update(save_figure_all_formats(fig, multipanel_stem, "figure8_multipanel"))
    plt.close(fig)

    return outputs


def write_placeholder_figures(
    axis_grade: pd.DataFrame,
    figure_dir: Path,
    scrna_recurrence: pd.DataFrame | None = None,
    bulk_assoc: pd.DataFrame | None = None,
) -> dict:
    return write_module8_nature_figures(
        manifest=pd.DataFrame(),
        axis_grade=axis_grade,
        figure_dir=figure_dir,
        scrna_recurrence=scrna_recurrence,
        bulk_assoc=bulk_assoc,
        survival_summary=None,
    )


def save_figure_all_formats(fig, stem: Path, key_prefix: str) -> dict:
    outputs = {}
    for suffix in [".png", ".pdf", ".svg"]:
        path = stem.with_suffix(suffix)
        save_kwargs = {"bbox_inches": "tight"}
        if suffix == ".png":
            save_kwargs["dpi"] = 600
        fig.savefig(path, **save_kwargs)
        outputs[f"{key_prefix}_{suffix.lstrip('.')}"] = str(path)
    return outputs


def build_module8_report_payload(
    inputs: dict,
    outputs: dict,
    manifest: pd.DataFrame,
    registry: pd.DataFrame,
    axis_grade: pd.DataFrame,
    risks: pd.DataFrame,
    recurrence: pd.DataFrame,
) -> dict:
    n_scored = (
        int(recurrence.loc[recurrence["included"].astype(bool), "dataset_id"].nunique())
        if not recurrence.empty and {"dataset_id", "included"}.issubset(recurrence.columns)
        else 0
    )
    return {
        "module": "8",
        "method": "External public-cohort validation framework for Module 7 TF axes",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": inputs,
        "outputs": outputs,
        "n_manifest_datasets": int(len(manifest)),
        "n_included_public_datasets": int(pd.Series(manifest.get("included", [])).astype(bool).sum()) if not manifest.empty else 0,
        "n_signature_tfs": int(registry["tf"].nunique()) if "tf" in registry.columns else 0,
        "n_axes": int(axis_grade["axis"].nunique()) if "axis" in axis_grade.columns else 0,
        "n_scored_external_cohorts": n_scored,
        "n_review_risk_flags": int(len(risks)),
        "external_result_status": "scored_external_cohorts" if n_scored else "pending_external_expression_matrices",
        "python_runtime": {"version": platform.python_version(), "platform": platform.platform()},
    }


def build_main_conclusions(
    axis_grade: pd.DataFrame,
    manifest: pd.DataFrame,
    recurrence: pd.DataFrame,
    bulk_assoc: pd.DataFrame | None = None,
    survival_summary: pd.DataFrame | None = None,
) -> str:
    lines = ["# Module 8 External Validation Conclusions", ""]
    lines.append("## Current Status")
    n_scored = (
        int(recurrence.loc[recurrence["included"].astype(bool), "dataset_id"].nunique())
        if not recurrence.empty and {"dataset_id", "included"}.issubset(recurrence.columns)
        else 0
    )
    if n_scored:
        lines.append(
            f"Module 8 scored {n_scored} external cohort(s) using frozen Module 7 TF-axis signatures. "
            "Evidence grades are based on current public scRNA/bulk inputs and should be updated when additional cohorts are added."
        )
    else:
        lines.append(
            "Module 8 has frozen the Module 7 TF-axis signatures and prepared a public-cohort validation framework. "
            "No external expression matrix has been scored yet, so recurrence grades remain provisional."
        )
    lines.append("")
    if not manifest.empty:
        included = manifest.loc[manifest["included"].astype(bool), "dataset_id"].astype(str).tolist()
        excluded = manifest.loc[~manifest["included"].astype(bool), "dataset_id"].astype(str).tolist()
        lines.append("## Dataset Audit")
        lines.append(f"Included public candidates: {', '.join(included) if included else 'none'}")
        lines.append(f"Excluded or pending candidates: {', '.join(excluded) if excluded else 'none'}")
        lines.append("")
    if not axis_grade.empty:
        lines.append("## Axis Evidence Grades")
        for _, row in axis_grade.iterrows():
            lines.append(f"- {row['axis']}: grade {row['evidence_grade']}, replicated cohorts {row['n_replicated_cohorts']}")
        lines.append("")
    if bulk_assoc is not None and not bulk_assoc.empty:
        lines.append("## Bulk Clinical Validation")
        tumor_normal = bulk_assoc.loc[
            bulk_assoc["clinical_variable"].astype(str).eq("tumor_vs_normal")
            & bulk_assoc["status"].astype(str).eq("tested")
        ]
        if not tumor_normal.empty:
            for dataset_id, dataset_df in tumor_normal.groupby("dataset_id", sort=True):
                n_positive = int(dataset_df["direction"].astype(str).eq("positive").sum())
                median_effect = float(dataset_df["effect_size"].median())
                lines.append(
                    f"- {dataset_id}: tumor-vs-normal recurrence was positive for {n_positive}/{len(dataset_df)} axes "
                    f"(median effect {median_effect:.3f})."
                )
        clinical_hits = bulk_assoc.loc[
            ~bulk_assoc["clinical_variable"].astype(str).eq("tumor_vs_normal")
            & bulk_assoc["status"].astype(str).eq("tested")
            & (bulk_assoc["p.adjust"] < 0.05)
        ].copy()
        if not clinical_hits.empty:
            for dataset_id, dataset_df in clinical_hits.groupby("dataset_id", sort=True):
                variables = sorted(dataset_df["clinical_variable"].astype(str).unique())
                lines.append(f"- {dataset_id}: significant clinical associations were detected for {', '.join(variables)}.")
        lines.append("")
    if survival_summary is not None and not survival_summary.empty:
        lines.append("## Adjusted Survival Models")
        for dataset_id, dataset_df in survival_summary.groupby("dataset_id", sort=True):
            adjusted = dataset_df.loc[dataset_df["status"].astype(str).eq("adjusted_cox")]
            if adjusted.empty:
                lines.append(f"- {dataset_id}: adjusted Cox models were not available.")
                continue
            best = adjusted.sort_values("pvalue").iloc[0]
            lines.append(
                f"- {dataset_id}: adjusted Cox models used {best['covariates_used']}; strongest axis was "
                f"{best['axis']} (HR {float(best['hazard_ratio']):.2f}, FDR {float(best['p.adjust']):.3g})."
            )
        lines.append("")
    return "\n".join(lines)


def add_local_sources_to_manifest(manifest: pd.DataFrame, sources: pd.DataFrame, discovery_ids: set[str]) -> pd.DataFrame:
    if sources.empty:
        return manifest
    existing = set(manifest["dataset_id"].astype(str))
    rows = []
    for _, source in sources.iterrows():
        if str(source["dataset_id"]) in existing:
            continue
        rows.append(
            {
                "dataset_id": source["dataset_id"],
                "source": source["source"],
                "modality": source["modality"],
                "cancer_type": "hepatocellular carcinoma",
                "species": "Homo sapiens",
                "n_samples": np.nan,
                "n_cells": np.nan,
                "tissue_type": "local external validation candidate",
                "access_url": source["expression_path"],
                "raw_available": True,
                "processed_available": True,
                "included": True,
                "exclusion_reason": "",
            }
        )
    if rows:
        manifest = pd.concat([manifest, pd.DataFrame(rows)], ignore_index=True)
    return audit_dataset_leakage(manifest, discovery_ids)


def bulk_assoc_to_recurrence(assoc: pd.DataFrame) -> pd.DataFrame:
    if assoc.empty:
        return pd.DataFrame(
            columns=[
                "dataset_id",
                "modality",
                "axis",
                "comparison",
                "n_positive",
                "n_reference",
                "effect_size",
                "pvalue",
                "p.adjust",
                "direction",
                "missingness_rate",
                "included",
                "dataset_leakage_status",
            ]
        )
    out = pd.DataFrame(
        {
            "dataset_id": assoc["dataset_id"],
            "modality": assoc["modality"],
            "axis": assoc["axis"],
            "comparison": assoc["clinical_variable"],
            "n_positive": assoc["n_samples"],
            "n_reference": 0,
            "effect_size": assoc["effect_size"],
            "pvalue": assoc["pvalue"],
            "p.adjust": assoc["p.adjust"],
            "direction": assoc["direction"],
            "missingness_rate": assoc["missingness_rate"],
            "included": assoc["status"].astype(str).eq("tested"),
            "dataset_leakage_status": "clean",
        }
    )
    return out


def run_module8(
    metadata_dir: Path,
    figure_dir: Path,
    top_n_targets: int = 50,
    local_scrna_root: Path | None = DEFAULT_LOCAL_SCRNA_ROOT,
    tcga_expression_path: Path | None = DEFAULT_TCGA_EXPRESSION_PATH,
    tcga_survival_path: Path | None = DEFAULT_TCGA_SURVIVAL_PATH,
    clinical_root: Path | None = DEFAULT_CLINICAL_ROOT,
) -> dict:
    tf_matrix_path = metadata_dir / "sctenifoldknk_module7_5_tf_level_replication_matrix.tsv"
    perturbation_path = metadata_dir / "sctenifoldknk_module7_2_driver_union_all_perturbation_genes.tsv"
    pathway_path = metadata_dir / "sctenifoldknk_module7_5_pathway_level_enrichment_matrix.tsv"

    tf_matrix = read_tsv_or_empty(tf_matrix_path)
    perturbation = read_tsv_or_empty(perturbation_path)
    pathway_matrix = read_tsv_or_empty(pathway_path)
    discovery_ids = extract_discovery_dataset_ids(metadata_dir)

    registry = build_signature_registry(tf_matrix)
    target_genes = build_tf_target_signature_genes(registry, perturbation, top_n=top_n_targets)
    pathway_genes = build_pathway_signature_genes(registry, pathway_matrix)
    combined_signature_genes = pd.concat(
        [
            target_genes[["axis", "tf", "gene"]].drop_duplicates() if not target_genes.empty else pd.DataFrame(),
            pathway_genes[["axis", "tf", "gene"]].drop_duplicates() if not pathway_genes.empty else pd.DataFrame(),
        ],
        ignore_index=True,
    ).drop_duplicates()

    local_sources = discover_local_scrna_sources(local_scrna_root)
    manifest = add_local_sources_to_manifest(build_external_dataset_manifest(discovery_ids), local_sources, discovery_ids)
    local_sources = local_sources.loc[~local_sources["dataset_id"].isin(discovery_ids)].reset_index(drop=True)
    scrna_scores, scrna_recurrence, scrna_detection = score_local_scrna_sources(local_sources, combined_signature_genes)

    bulk_assoc_frames = []
    survival_frames = []
    bulk_recurrence_frames = []
    if tcga_expression_path and Path(tcga_expression_path).exists() and not combined_signature_genes.empty:
        gtf_path = Path(r"G:\wanyi_HCC_scRNA\HCCscRNA\GSE156625-HCC\cellranger\hg38\refdata-gex-GRCh38-2020-A\genes\genes.gtf")
        ensembl_map = parse_gtf_ensembl_to_symbol(gtf_path)
        tcga_expr = load_tcga_expression_signature_genes(tcga_expression_path, combined_signature_genes, ensembl_map)
        bulk_scores, bulk_detection = compute_bulk_signature_scores(tcga_expr, combined_signature_genes, "TCGA-LIHC")
        tcga_assoc = compute_bulk_tumor_normal_association(bulk_scores)
        tcga_clinical = pd.DataFrame()
        if clinical_root:
            tcga_complete_clinical_path = Path(clinical_root) / "clinical.xls"
            tcga_binary_clinical_path = Path(clinical_root) / "tcgaClinical.txt"
            if tcga_complete_clinical_path.exists():
                tcga_clinical = prepare_tcga_clinical_covariates(load_tcga_clinical_table(tcga_complete_clinical_path))
            elif tcga_binary_clinical_path.exists():
                tcga_clinical = pd.read_csv(tcga_binary_clinical_path, sep="\t")
            if not tcga_clinical.empty:
                tcga_variables = [column for column in ["age", "gender", "grade", "stage", "T", "M", "N"] if column in tcga_clinical.columns]
                tcga_assoc = pd.concat(
                    [
                        tcga_assoc,
                        compute_bulk_clinical_variable_association(
                            bulk_scores,
                            tcga_clinical,
                            dataset_id="TCGA-LIHC",
                            variables=tcga_variables,
                            clinical_id_col="Id",
                        ),
                    ],
                    ignore_index=True,
                )
        if not tcga_clinical.empty and {"Id", "fustat", "futime"}.issubset(tcga_clinical.columns):
            survival_frames.append(
                compute_adjusted_survival_association(
                    bulk_scores,
                    tcga_clinical,
                    tcga_clinical,
                    dataset_id="TCGA-LIHC",
                    survival_id_col="Id",
                    event_col="fustat",
                    time_col="futime",
                    clinical_id_col="Id",
                    covariates=["age", "gender", "stage"],
                )
            )
        elif tcga_survival_path and Path(tcga_survival_path).exists():
            survival = pd.read_csv(tcga_survival_path, sep="\t")
            if not tcga_clinical.empty:
                survival_frames.append(
                    compute_adjusted_survival_association(
                        bulk_scores,
                        survival,
                        tcga_clinical,
                        dataset_id="TCGA-LIHC",
                        survival_id_col="_PATIENT",
                        event_col="OS",
                        time_col="OS.time",
                        clinical_id_col="Id",
                        covariates=["age", "gender", "stage"],
                    )
                )
            else:
                survival_frames.append(compute_exploratory_survival_association(bulk_scores, survival))
        bulk_assoc_frames.append(tcga_assoc)
        bulk_recurrence_frames.append(bulk_assoc_to_recurrence(tcga_assoc.loc[tcga_assoc["clinical_variable"].eq("tumor_vs_normal")]))

    if clinical_root:
        clinical_root = Path(clinical_root)
        icgc_expression_path = clinical_root / "ICGCsymbol.txt"
        icgc_survival_path = clinical_root / "ICGCtime.txt"
        icgc_clinical_path = clinical_root / "icgcClinical.txt"
        if icgc_expression_path.exists() and not combined_signature_genes.empty:
            icgc_expr = load_symbol_expression_signature_genes(icgc_expression_path, combined_signature_genes)
            icgc_scores, icgc_detection = compute_bulk_signature_scores(
                icgc_expr, combined_signature_genes, "ICGC-LIRI-JP", sample_type_mode="icgc"
            )
            icgc_assoc = compute_bulk_tumor_normal_association(icgc_scores)
            if icgc_clinical_path.exists():
                icgc_clinical = pd.read_csv(icgc_clinical_path, sep="\t")
                icgc_assoc = pd.concat(
                    [
                        icgc_assoc,
                        compute_bulk_clinical_variable_association(
                            icgc_scores,
                            icgc_clinical,
                            dataset_id="ICGC-LIRI-JP",
                            variables=["Age", "Gender", "Stage"],
                            clinical_id_col="Id",
                        ),
                    ],
                    ignore_index=True,
                )
            if icgc_survival_path.exists() and icgc_clinical_path.exists():
                icgc_survival = pd.read_csv(icgc_survival_path, sep="\t")
                icgc_clinical = pd.read_csv(icgc_clinical_path, sep="\t")
                survival_frames.append(
                    compute_adjusted_survival_association(
                        icgc_scores,
                        icgc_survival,
                        icgc_clinical,
                        dataset_id="ICGC-LIRI-JP",
                        survival_id_col="id",
                        event_col="fustat",
                        time_col="futime",
                        clinical_id_col="Id",
                        covariates=["Age", "Gender", "Stage"],
                    )
                )
            bulk_assoc_frames.append(icgc_assoc)
            bulk_recurrence_frames.append(bulk_assoc_to_recurrence(icgc_assoc.loc[icgc_assoc["clinical_variable"].eq("tumor_vs_normal")]))

    bulk_assoc = pd.concat(bulk_assoc_frames, ignore_index=True) if bulk_assoc_frames else empty_bulk_schema()
    survival_summary = pd.concat(survival_frames, ignore_index=True) if survival_frames else pd.DataFrame(columns=survival_summary_columns())
    bulk_recurrence = pd.concat(bulk_recurrence_frames, ignore_index=True) if bulk_recurrence_frames else bulk_assoc_to_recurrence(empty_bulk_schema())

    spatial_qc = manifest.loc[manifest["modality"].eq("spatial_transcriptomics")].copy()
    spatial_scores = empty_spatial_schema()
    recurrence = pd.concat([scrna_recurrence, bulk_recurrence], ignore_index=True)
    axis_grade = build_axis_level_evidence_grade(recurrence, known_axes=MODULE8_AXES)
    risks = flag_external_control_outperformance(axis_grade)
    figure_outputs = write_module8_nature_figures(
        manifest=manifest,
        axis_grade=axis_grade,
        figure_dir=figure_dir,
        scrna_recurrence=scrna_recurrence,
        bulk_assoc=bulk_assoc,
        survival_summary=survival_summary,
    )

    outputs = {
        "external_dataset_manifest": str(metadata_dir / "module8_external_dataset_manifest.tsv"),
        "signature_registry": str(metadata_dir / "module8_signature_registry.tsv"),
        "tf_target_signature_genes": str(metadata_dir / "module8_tf_target_signature_genes.tsv"),
        "pathway_signature_genes": str(metadata_dir / "module8_pathway_signature_genes.tsv"),
        "scrna_dataset_level_scores": str(metadata_dir / "module8_scrna_dataset_level_scores.tsv"),
        "scrna_cell_state_recurrence": str(metadata_dir / "module8_scrna_cell_state_recurrence.tsv"),
        "bulk_signature_clinical_association": str(metadata_dir / "module8_bulk_signature_clinical_association.tsv"),
        "bulk_survival_model_summary": str(metadata_dir / "module8_bulk_survival_model_summary.tsv"),
        "spatial_dataset_qc": str(metadata_dir / "module8_spatial_dataset_qc.tsv"),
        "spatial_region_signature_scores": str(metadata_dir / "module8_spatial_region_signature_scores.tsv"),
        "scrna_gene_detection": str(metadata_dir / "module8_scrna_gene_detection.tsv"),
        "external_recurrence_matrix": str(metadata_dir / "module8_external_recurrence_matrix.tsv"),
        "axis_level_evidence_grade": str(metadata_dir / "module8_axis_level_evidence_grade.tsv"),
        "review_risk_flags": str(metadata_dir / "module8_review_risk_flags.tsv"),
        "report": str(metadata_dir / "module8_external_validation_report.json"),
        "main_conclusions": str(metadata_dir / "module8_external_validation_main_conclusions.md"),
        "supplementary_index": str(metadata_dir / "module8_external_validation_supplementary_index.tsv"),
        **figure_outputs,
    }
    inputs = {
        "tf_level_replication_matrix": str(tf_matrix_path),
        "driver_union_perturbation_genes": str(perturbation_path),
        "pathway_level_enrichment_matrix": str(pathway_path),
        "local_scrna_root": str(local_scrna_root) if local_scrna_root else "",
        "tcga_expression": str(tcga_expression_path) if tcga_expression_path else "",
        "tcga_survival": str(tcga_survival_path) if tcga_survival_path else "",
        "clinical_root": str(clinical_root) if clinical_root else "",
        "tcga_complete_clinical": str(Path(clinical_root) / "clinical.xls") if clinical_root else "",
        "tcga_binary_clinical_fallback": str(Path(clinical_root) / "tcgaClinical.txt") if clinical_root else "",
        "icgc_clinical": str(Path(clinical_root) / "icgcClinical.txt") if clinical_root else "",
    }

    metadata_dir.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(outputs["external_dataset_manifest"], sep="\t", index=False)
    registry.to_csv(outputs["signature_registry"], sep="\t", index=False)
    target_genes.to_csv(outputs["tf_target_signature_genes"], sep="\t", index=False)
    pathway_genes.to_csv(outputs["pathway_signature_genes"], sep="\t", index=False)
    scrna_scores.to_csv(outputs["scrna_dataset_level_scores"], sep="\t", index=False)
    scrna_recurrence.to_csv(outputs["scrna_cell_state_recurrence"], sep="\t", index=False)
    bulk_assoc.to_csv(outputs["bulk_signature_clinical_association"], sep="\t", index=False)
    survival_summary.to_csv(outputs["bulk_survival_model_summary"], sep="\t", index=False)
    spatial_qc.to_csv(outputs["spatial_dataset_qc"], sep="\t", index=False)
    spatial_scores.to_csv(outputs["spatial_region_signature_scores"], sep="\t", index=False)
    scrna_detection.to_csv(outputs["scrna_gene_detection"], sep="\t", index=False)
    recurrence.to_csv(outputs["external_recurrence_matrix"], sep="\t", index=False)
    axis_grade.to_csv(outputs["axis_level_evidence_grade"], sep="\t", index=False)
    risks.to_csv(outputs["review_risk_flags"], sep="\t", index=False)

    supplementary = pd.DataFrame(
        [
            {
                "table_id": "Supplementary Table 8.1",
                "primary_source_file": Path(outputs["external_dataset_manifest"]).name,
                "description": "Public external cohort manifest and leakage audit.",
                "intended_manuscript_use": "External validation dataset audit.",
            },
            {
                "table_id": "Supplementary Table 8.2",
                "primary_source_file": Path(outputs["signature_registry"]).name,
                "description": "Frozen Module 7 TF-axis signature registry.",
                "intended_manuscript_use": "External validation methods and reproducibility.",
            },
            {
                "table_id": "Supplementary Table 8.3",
                "primary_source_file": Path(outputs["axis_level_evidence_grade"]).name,
                "description": "Axis-level external recurrence score and provisional evidence grade.",
                "intended_manuscript_use": "Cross-cohort recurrence summary.",
            },
            {
                "table_id": "Figure 8A source",
                "primary_source_file": Path(outputs["external_dataset_manifest"]).name,
                "description": "External cohort map and leakage-audited inclusion flow.",
                "intended_manuscript_use": "Main Figure 8A.",
            },
            {
                "table_id": "Figure 8B source",
                "primary_source_file": Path(outputs["external_recurrence_matrix"]).name,
                "description": "scRNA signature recurrence effect sizes across external cohorts.",
                "intended_manuscript_use": "Main Figure 8B.",
            },
            {
                "table_id": "Figure 8C source",
                "primary_source_file": Path(outputs["bulk_signature_clinical_association"]).name,
                "description": "TCGA/ICGC tumor-normal and clinical variable associations.",
                "intended_manuscript_use": "Main Figure 8C clinical heatmap.",
            },
            {
                "table_id": "Figure 8C survival source",
                "primary_source_file": Path(outputs["bulk_survival_model_summary"]).name,
                "description": "TCGA/ICGC adjusted Cox models for TF-axis signature scores.",
                "intended_manuscript_use": "Main Figure 8C forest plot.",
            },
            {
                "table_id": "Figure 8D source",
                "primary_source_file": Path(outputs["axis_level_evidence_grade"]).name,
                "description": "Axis-level recurrence score, evidence grade and replicated cohort count.",
                "intended_manuscript_use": "Main Figure 8D.",
            },
        ]
    )
    supplementary.to_csv(outputs["supplementary_index"], sep="\t", index=False)

    report = build_module8_report_payload(inputs, outputs, manifest, registry, axis_grade, risks, recurrence)
    Path(outputs["report"]).write_text(json.dumps(report, indent=2), encoding="utf-8")
    Path(outputs["main_conclusions"]).write_text(
        build_main_conclusions(axis_grade, manifest, recurrence, bulk_assoc, survival_summary),
        encoding="utf-8",
    )
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 8 external public-cohort validation framework")
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--top-n-targets", type=int, default=50)
    parser.add_argument("--local-scrna-root", type=Path, default=DEFAULT_LOCAL_SCRNA_ROOT)
    parser.add_argument("--tcga-expression", type=Path, default=DEFAULT_TCGA_EXPRESSION_PATH)
    parser.add_argument("--tcga-survival", type=Path, default=DEFAULT_TCGA_SURVIVAL_PATH)
    parser.add_argument("--clinical-root", type=Path, default=DEFAULT_CLINICAL_ROOT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    report = run_module8(
        args.metadata_dir,
        args.figure_dir,
        top_n_targets=args.top_n_targets,
        local_scrna_root=args.local_scrna_root,
        tcga_expression_path=args.tcga_expression,
        tcga_survival_path=args.tcga_survival,
        clinical_root=args.clinical_root,
    )
    print(
        json.dumps(
            {
                "report": report["outputs"]["report"],
                "n_manifest_datasets": report["n_manifest_datasets"],
                "n_signature_tfs": report["n_signature_tfs"],
                "external_result_status": report["external_result_status"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
