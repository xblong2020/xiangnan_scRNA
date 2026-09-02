from __future__ import annotations

import argparse
import re
import json
import math
import time
from io import StringIO
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import requests
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
AUTOSOMES = [str(i) for i in range(1, 23)]
HGNC_COMPLETE_SET_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
NON_EPITHELIAL_SEEDS = {
    "T_NK",
    "Myeloid",
    "Endothelial",
    "Fibroblast_HSC_Pericyte",
    "B_cell",
    "Plasma_cell",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 3 malignant HCC evidence integration.")
    parser.add_argument(
        "--integrated",
        type=Path,
        default=ROOT / "data/processed/scvi/scvi_integrated_counts_hvg.celltypist_major.h5ad",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=ROOT / "metadata/scvi/scvi_input_manifest.counts.tsv",
    )
    parser.add_argument(
        "--major-seed",
        type=Path,
        default=ROOT / "metadata/celltype/scanvi_seed_labels_by_cell.tsv.gz",
    )
    parser.add_argument(
        "--hepatocyte-cells",
        type=Path,
        default=ROOT / "metadata/hepatocyte/hepatocyte_lineage_cells.tsv.gz",
    )
    parser.add_argument(
        "--cnv-candidate-cells",
        type=Path,
        default=ROOT / "metadata/hepatocyte/hepatocyte_cnv_candidate_cells.tsv.gz",
    )
    parser.add_argument(
        "--hepatocyte-subclusters",
        type=Path,
        default=ROOT / "metadata/hepatocyte/hepatocyte_subcluster_summary.tsv",
    )
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/malignant")
    parser.add_argument("--gene-map", type=Path, default=ROOT / "metadata/malignant/gene_chromosome_map_hgnc.tsv")
    parser.add_argument("--gene-bin-size", type=int, default=100)
    parser.add_argument("--min-candidate-cells", type=int, default=20)
    parser.add_argument("--min-reference-cells", type=int, default=100)
    parser.add_argument("--max-reference-cells", type=int, default=8000)
    parser.add_argument("--cnv-z-threshold", type=float, default=3.0)
    parser.add_argument("--cnv-high-bin-fraction", type=float, default=0.08)
    parser.add_argument("--seed", type=int, default=20260601)
    return parser.parse_args()


def read_manifest(path: Path) -> pd.DataFrame:
    manifest = pd.read_csv(path, sep="\t")
    manifest = manifest.loc[manifest["include_in_scvi"].astype(str).str.lower().eq("true")].copy()
    manifest["study_sample"] = manifest["dataset"].astype(str) + "__" + manifest["label"].astype(str)
    return manifest


def resolve_path(value: object) -> Path:
    path = Path(str(value))
    if path.exists():
        return path
    text = str(value).replace("\\", "/")
    for anchor in ("data/processed/", "metadata/"):
        idx = text.find(anchor)
        if idx >= 0:
            candidate = ROOT / Path(text[idx:])
            if candidate.exists():
                return candidate
    return path


def as_csr(x) -> sparse.csr_matrix:
    if sparse.issparse(x):
        return x.tocsr()
    return sparse.csr_matrix(x)


def original_id(cell_id: str, study_sample: str) -> str:
    prefix = f"{study_sample}__"
    return cell_id[len(prefix) :] if cell_id.startswith(prefix) else cell_id


def cnv_sample_id(cell_id: str, dataset: str, study_sample: str, sample_id: str) -> str:
    orig = original_id(cell_id, study_sample)
    if dataset == "GSE149614" and "_" in orig:
        return orig.split("_", 1)[0]
    return str(sample_id) if str(sample_id) and str(sample_id) != "nan" else str(study_sample)


def sample_source_class(dataset: str, sample: str) -> str:
    sample = str(sample)
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
    if dataset == "GSE202379":
        return "non_hcc_liver"
    if dataset == "GSE174748":
        return "non_hcc_liver"
    if dataset in {"GSE151530", "GSE212046"}:
        return "unknown_hcc_dataset"
    return "unknown"


def hgnc_location_chromosome(location: object) -> str | None:
    match = re.match(r"^(?:CHR)?([0-9]{1,2}|X|Y)", str(location).upper())
    if not match:
        return None
    chrom = match.group(1)
    return chrom if chrom in AUTOSOMES else None


def split_hgnc_symbol_list(value: object) -> list[str]:
    if pd.isna(value) or str(value).strip() == "":
        return []
    text = str(value).replace('"', "")
    return [part.strip() for part in text.split("|") if part.strip()]


def download_hgnc_gene_map(symbols: list[str]) -> pd.DataFrame:
    response = requests.get(HGNC_COMPLETE_SET_URL, timeout=120)
    response.raise_for_status()
    hgnc = pd.read_csv(StringIO(response.text), sep="\t", dtype=str).fillna("")
    hgnc = hgnc.loc[hgnc["status"].eq("Approved")].copy()
    hgnc["chromosome"] = hgnc["location"].map(hgnc_location_chromosome)
    hgnc = hgnc.loc[hgnc["chromosome"].isin(AUTOSOMES)].copy()
    hgnc["chrom_order"] = hgnc["chromosome"].astype(int)
    hgnc = hgnc.sort_values(["chrom_order", "location_sortable", "symbol"], kind="mergesort")
    hgnc["approx_start"] = hgnc.groupby("chromosome", sort=False).cumcount() + 1

    wanted = set(map(str, symbols))
    rows = []
    seen = set()
    for _, row in hgnc.iterrows():
        candidates = [row["symbol"]]
        candidates.extend(split_hgnc_symbol_list(row.get("prev_symbol", "")))
        candidates.extend(split_hgnc_symbol_list(row.get("alias_symbol", "")))
        for gene in candidates:
            if gene not in wanted or gene in seen:
                continue
            rows.append(
                {
                    "gene": gene,
                    "symbol": row["symbol"],
                    "chromosome": row["chromosome"],
                    "start": int(row["approx_start"]),
                    "end": int(row["approx_start"]),
                    "location": row["location"],
                    "location_sortable": row["location_sortable"],
                    "source": "HGNC complete set cytoband location",
                }
            )
            seen.add(gene)
    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(columns=["gene", "symbol", "chromosome", "start", "end", "location", "location_sortable", "source"])
    out["chrom_order"] = out["chromosome"].astype(int)
    out = out.sort_values(["chrom_order", "start", "gene"], kind="mergesort").drop(columns=["chrom_order"])
    return out.drop_duplicates("gene")


def load_or_create_gene_map(manifest: pd.DataFrame, gene_map_path: Path) -> pd.DataFrame:
    gene_map_path.parent.mkdir(parents=True, exist_ok=True)
    if gene_map_path.exists():
        gene_map = pd.read_csv(gene_map_path, sep="\t")
        required = {"gene", "chromosome", "start", "end"}
        if required.issubset(gene_map.columns):
            return gene_map
    genes = set()
    for _, row in manifest.iterrows():
        a = ad.read_h5ad(resolve_path(row["output"]), backed="r")
        genes.update(map(str, a.var_names))
        a.file.close()
    print(f"DOWNLOAD HGNC gene map symbols={len(genes)}", flush=True)
    gene_map = download_hgnc_gene_map(sorted(genes))
    gene_map.to_csv(gene_map_path, sep="\t", index=False)
    return gene_map


def make_gene_bins(gene_map: pd.DataFrame, genes: pd.Index, bin_size: int) -> tuple[pd.DataFrame, sparse.csr_matrix]:
    present = gene_map.loc[gene_map["gene"].isin(set(genes.astype(str)))].copy()
    present["gene_idx"] = present["gene"].map(pd.Series(np.arange(len(genes)), index=genes.astype(str))).astype(int)
    present["chromosome"] = present["chromosome"].astype(str)
    present["chrom_order"] = present["chromosome"].astype(int)
    present["start"] = pd.to_numeric(present["start"], errors="coerce").fillna(0).astype(int)
    present["end"] = pd.to_numeric(present["end"], errors="coerce").fillna(present["start"]).astype(int)
    present = present.sort_values(["chrom_order", "start", "gene"], kind="mergesort")
    bin_rows = []
    gene_bin_pairs = []
    bin_idx = 0
    for chrom in AUTOSOMES:
        sub = present.loc[present["chromosome"].eq(chrom)].reset_index(drop=True)
        for start in range(0, sub.shape[0], bin_size):
            chunk = sub.iloc[start : start + bin_size]
            if chunk.shape[0] < max(20, bin_size // 4):
                continue
            bin_rows.append(
                {
                    "bin_id": f"chr{chrom}_bin{start // bin_size:03d}",
                    "chromosome": chrom,
                    "start": int(chunk["start"].min()),
                    "end": int(chunk["end"].max()),
                    "n_genes": int(chunk.shape[0]),
                    "genes": ";".join(chunk["gene"].astype(str).tolist()),
                }
            )
            for gene_idx in chunk["gene_idx"].astype(int).tolist():
                gene_bin_pairs.append((gene_idx, bin_idx))
            bin_idx += 1
    bins = pd.DataFrame(bin_rows)
    indicator = sparse.csr_matrix(
        (
            np.ones(len(gene_bin_pairs), dtype=np.float64),
            ([pair[0] for pair in gene_bin_pairs], [pair[1] for pair in gene_bin_pairs]),
        ),
        shape=(len(genes), bins.shape[0]),
    )
    return bins, indicator


def robust_mad(values: np.ndarray) -> float:
    med = np.median(values)
    mad = np.median(np.abs(values - med)) * 1.4826
    return float(mad) if mad > 1e-8 else float(np.std(values) + 1e-8)


def compute_cnv_proxy_for_sample(
    a: ad.AnnData,
    sample_cells: pd.DataFrame,
    ref_cells: pd.DataFrame,
    gene_map: pd.DataFrame,
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    all_cells = pd.concat([sample_cells, ref_cells], axis=0)
    original_ids = [original_id(cell_id, ss) for cell_id, ss in zip(all_cells["cell_id"], all_cells["study_sample"])]
    present_ids = pd.Index(original_ids).intersection(a.obs_names.astype(str))
    if len(present_ids) == 0:
        raise ValueError("No cells matched source h5ad")
    keep = pd.Series(np.arange(len(original_ids)), index=original_ids).loc[present_ids].to_numpy()
    all_cells = all_cells.iloc[keep].copy()
    is_candidate = all_cells["cnv_role"].eq("candidate").to_numpy()
    is_ref = all_cells["cnv_role"].eq("reference").to_numpy()
    x = as_csr(a[present_ids, :].X).astype(np.float64)
    bins, indicator = make_gene_bins(gene_map, pd.Index(a.var_names.astype(str)), args.gene_bin_size)
    if bins.empty:
        raise ValueError("No gene bins available")
    bin_counts = (x @ indicator).astype(np.float64)
    totals = np.asarray(x.sum(axis=1)).ravel()
    fractions = np.divide(
        bin_counts.toarray(),
        totals[:, None],
        out=np.zeros((x.shape[0], bins.shape[0]), dtype=np.float64),
        where=totals[:, None] > 0,
    )
    baseline = np.median(fractions[is_ref, :], axis=0)
    eps = max(float(np.median(baseline[baseline > 0])) * 0.05, 1e-9) if np.any(baseline > 0) else 1e-9
    log2ratio = np.log2((fractions + eps) / (baseline[None, :] + eps))
    burden = np.mean(np.abs(log2ratio), axis=1)
    high_frac = np.mean(np.abs(log2ratio) >= 0.25, axis=1)
    max_abs = np.max(np.abs(log2ratio), axis=1)
    ref_burden = burden[is_ref]
    ref_med = float(np.median(ref_burden))
    ref_mad = robust_mad(ref_burden)
    z = (burden - ref_med) / ref_mad
    all_cells["cnv_proxy_burden"] = burden
    all_cells["cnv_proxy_z"] = z
    all_cells["cnv_proxy_high_bin_fraction"] = high_frac
    all_cells["cnv_proxy_max_abs_bin_log2"] = max_abs
    all_cells["cnv_proxy_status"] = np.where(
        (z >= args.cnv_z_threshold) & (high_frac >= args.cnv_high_bin_fraction),
        "aneuploid_proxy",
        np.where((z >= 2.0) | (high_frac >= args.cnv_high_bin_fraction), "borderline_proxy", "diploid_like_proxy"),
    )
    bin_stats = bins.copy()
    candidate_log2 = log2ratio[is_candidate, :]
    ref_log2 = log2ratio[is_ref, :]
    bin_stats["candidate_mean_log2"] = candidate_log2.mean(axis=0)
    bin_stats["candidate_median_log2"] = np.median(candidate_log2, axis=0)
    bin_stats["candidate_abs_mean_log2"] = np.mean(np.abs(candidate_log2), axis=0)
    bin_stats["reference_abs_mean_log2"] = np.mean(np.abs(ref_log2), axis=0)
    summary = pd.DataFrame(
        [
            {
                "cnv_sample": sample_cells["cnv_sample"].iloc[0],
                "study_sample": sample_cells["study_sample"].iloc[0],
                "dataset": sample_cells["dataset"].iloc[0],
                "n_candidate": int(is_candidate.sum()),
                "n_reference": int(is_ref.sum()),
                "n_bins": int(bins.shape[0]),
                "reference_burden_median": ref_med,
                "reference_burden_mad": ref_mad,
                "candidate_burden_median": float(np.median(burden[is_candidate])),
                "candidate_burden_mean": float(np.mean(burden[is_candidate])),
                "candidate_z_median": float(np.median(z[is_candidate])),
                "aneuploid_proxy_fraction": float(np.mean(all_cells.loc[is_candidate, "cnv_proxy_status"].eq("aneuploid_proxy"))),
                "borderline_proxy_fraction": float(np.mean(all_cells.loc[is_candidate, "cnv_proxy_status"].eq("borderline_proxy"))),
            }
        ]
    )
    return all_cells.loc[is_candidate].copy(), summary, bin_stats


def build_module3_inputs(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cnv_candidates = pd.read_csv(args.cnv_candidate_cells, sep="\t")
    hep = pd.read_csv(args.hepatocyte_cells, sep="\t")
    major = pd.read_csv(args.major_seed, sep="\t", usecols=["cell_id", "scanvi_seed_label_major", "excluded_doublet_cluster"])
    integrated = ad.read_h5ad(args.integrated, backed="r")
    obs = integrated.obs[["dataset", "study_sample", "sample_id", "major_celltype", "celltypist_liver_confidence"]].copy()
    obs.insert(0, "cell_id", integrated.obs_names.astype(str))
    integrated.file.close()
    major = major.merge(obs, on="cell_id", how="left", suffixes=("", "_obs"))
    cnv_candidates = cnv_candidates.merge(obs[["cell_id", "sample_id", "major_celltype"]], on="cell_id", how="left")
    hep = hep.merge(obs[["cell_id", "sample_id", "major_celltype"]], on="cell_id", how="left")

    for df in (cnv_candidates, hep, major):
        df["cnv_sample"] = [
            cnv_sample_id(c, d, s, sid)
            for c, d, s, sid in zip(df["cell_id"], df["dataset"], df["study_sample"], df["sample_id"])
        ]
        df["sample_source_class"] = [sample_source_class(d, sample) for d, sample in zip(df["dataset"], df["cnv_sample"])]

    refs = major.loc[
        major["scanvi_seed_label_major"].isin(NON_EPITHELIAL_SEEDS)
        & ~major["excluded_doublet_cluster"].fillna(False).astype(bool)
    ].copy()
    return cnv_candidates, hep, refs


def integrate_malignant_calls(
    cnv_calls: pd.DataFrame,
    cnv_candidates: pd.DataFrame,
    subclusters: pd.DataFrame,
    copykat_available: bool,
    infercnv_available: bool,
) -> pd.DataFrame:
    sub = subclusters.copy()
    sub["leiden_hep"] = sub["leiden_hep"].astype(str)
    evidence_cols = [
        "leiden_hep",
        "hcc_malignant_associated_score_z",
        "hcc_malignant_associated_mean_log1p_cpm",
        "hcc_malignant_associated_mean_pct_expr",
        "proliferation_score_z",
        "proliferation_mean_pct_expr",
        "regenerative_progenitor_score_z",
        "state_evidence",
    ]
    cnv_evidence_cols = [
        "cell_id",
        "cnv_proxy_burden",
        "cnv_proxy_z",
        "cnv_proxy_high_bin_fraction",
        "cnv_proxy_max_abs_bin_log2",
        "cnv_proxy_status",
    ]
    cnv_for_merge = cnv_calls[[col for col in cnv_evidence_cols if col in cnv_calls.columns]].drop_duplicates("cell_id")
    out = cnv_candidates.merge(cnv_for_merge, on="cell_id", how="left")
    out["leiden_hep"] = out["leiden_hep"].astype(str)
    out = out.merge(sub[evidence_cols], on="leiden_hep", how="left")
    out["copykat_status"] = "not_run_package_unavailable" if not copykat_available else "not_run"
    out["infercnv_status"] = "not_run_package_unavailable" if not infercnv_available else "not_run"
    marker_high = (
        (pd.to_numeric(out["hcc_malignant_associated_score_z"], errors="coerce") >= 0.8)
        | (pd.to_numeric(out["hcc_malignant_associated_mean_log1p_cpm"], errors="coerce") >= 3.5)
        | out["hepatocyte_state_label"].astype(str).str.contains("malignant_hepatocyte_candidate", regex=False)
    )
    prolif_high = (
        (pd.to_numeric(out["proliferation_score_z"], errors="coerce") >= 0.8)
        | out["hepatocyte_state_label"].astype(str).str.contains("proliferating", regex=False)
    )
    tumor_source = out["sample_source_class"].isin(
        ["tumor", "pvtt_tumor", "metastatic_tumor_lymphnode", "unknown_hcc_dataset"]
    )
    cnv_high = out["cnv_proxy_status"].eq("aneuploid_proxy")
    cnv_border = out["cnv_proxy_status"].eq("borderline_proxy")
    out["malignant_hcc_call"] = "non_malignant_or_unresolved"
    out.loc[cnv_high & tumor_source & marker_high, "malignant_hcc_call"] = "malignant_hcc_high_conf"
    out.loc[cnv_high & tumor_source & ~marker_high, "malignant_hcc_call"] = "malignant_hcc_cnv_support"
    out.loc[cnv_border & tumor_source & marker_high, "malignant_hcc_call"] = "malignant_hcc_probable"
    out.loc[~cnv_high & tumor_source & marker_high & prolif_high, "malignant_hcc_call"] = "malignant_hcc_marker_proliferation_needs_cnv_review"
    out.loc[out["sample_source_class"].isin(["normal_adjacent", "non_hcc_liver"]) & ~cnv_high, "malignant_hcc_call"] = "not_malignant_source_or_cnv"
    out.loc[out["cnv_proxy_status"].isna(), "malignant_hcc_call"] = "cnv_not_available"
    out["malignant_hcc_evidence"] = (
        "cnv_proxy="
        + out["cnv_proxy_status"].fillna("NA").astype(str)
        + "; source="
        + out["sample_source_class"].fillna("NA").astype(str)
        + "; hcc_score_z="
        + pd.to_numeric(out["hcc_malignant_associated_score_z"], errors="coerce").round(3).astype(str)
        + "; prolif_score_z="
        + pd.to_numeric(out["proliferation_score_z"], errors="coerce").round(3).astype(str)
        + "; state="
        + out["hepatocyte_state_label"].astype(str)
    )
    return out


def main() -> int:
    args = parse_args()
    start = time.time()
    np.random.seed(args.seed)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    cnv_candidates, hepatocyte_cells, refs = build_module3_inputs(args)
    manifest = read_manifest(args.manifest)
    gene_map = load_or_create_gene_map(manifest, args.gene_map)
    subclusters = pd.read_csv(args.hepatocyte_subclusters, sep="\t")

    copykat_available = False
    infercnv_available = False
    try:
        import subprocess

        r_cmd = [
            r"C:\Program Files\R\R-4.6.0\bin\Rscript.exe",
            "-e",
            "cat(requireNamespace('copykat', quietly=TRUE), requireNamespace('infercnv', quietly=TRUE), sep='\\t')",
        ]
        result = subprocess.run(r_cmd, capture_output=True, text=True, timeout=60)
        parts = result.stdout.strip().split("\t")
        copykat_available = parts[0].upper() == "TRUE" if parts else False
        infercnv_available = parts[1].upper() == "TRUE" if len(parts) > 1 else False
    except Exception:
        pass

    candidate_rows = []
    sample_summaries = []
    bin_rows = []
    skipped = []
    manifest_by_study = manifest.set_index("study_sample").to_dict(orient="index")
    for (study_sample, cnv_sample), sample_candidates in cnv_candidates.groupby(["study_sample", "cnv_sample"], observed=True):
        sample_candidates = sample_candidates.copy()
        if sample_candidates.shape[0] < args.min_candidate_cells:
            skipped.append({"study_sample": study_sample, "cnv_sample": cnv_sample, "reason": "candidate_cells_below_threshold", "n_candidate": sample_candidates.shape[0]})
            continue
        same_sample_refs = refs.loc[refs["study_sample"].eq(study_sample) & refs["cnv_sample"].eq(cnv_sample)].copy()
        if same_sample_refs.shape[0] < args.min_reference_cells:
            same_sample_refs = refs.loc[refs["study_sample"].eq(study_sample)].copy()
        if same_sample_refs.shape[0] < args.min_reference_cells:
            same_sample_refs = refs.loc[refs["dataset"].eq(sample_candidates["dataset"].iloc[0])].copy()
        if same_sample_refs.shape[0] < args.min_reference_cells:
            skipped.append({"study_sample": study_sample, "cnv_sample": cnv_sample, "reason": "reference_cells_below_threshold", "n_candidate": sample_candidates.shape[0], "n_reference": same_sample_refs.shape[0]})
            continue
        if same_sample_refs.shape[0] > args.max_reference_cells:
            same_sample_refs = same_sample_refs.sample(args.max_reference_cells, random_state=args.seed)
        row = manifest_by_study.get(study_sample)
        if row is None:
            skipped.append({"study_sample": study_sample, "cnv_sample": cnv_sample, "reason": "missing_manifest", "n_candidate": sample_candidates.shape[0]})
            continue
        sample_candidates["cnv_role"] = "candidate"
        same_sample_refs["cnv_role"] = "reference"
        print(f"CNV_PROXY {cnv_sample} candidates={sample_candidates.shape[0]} refs={same_sample_refs.shape[0]}", flush=True)
        a = ad.read_h5ad(resolve_path(row["output"]), backed="r")
        try:
            calls, summary, bins = compute_cnv_proxy_for_sample(a, sample_candidates, same_sample_refs, gene_map, args)
        finally:
            a.file.close()
        candidate_rows.append(calls)
        sample_summaries.append(summary)
        bins["cnv_sample"] = cnv_sample
        bins["study_sample"] = study_sample
        bins["dataset"] = sample_candidates["dataset"].iloc[0]
        bin_rows.append(bins)

    cnv_calls = pd.concat(candidate_rows, ignore_index=True) if candidate_rows else pd.DataFrame(columns=["cell_id"])
    sample_summary = pd.concat(sample_summaries, ignore_index=True) if sample_summaries else pd.DataFrame()
    bin_summary = pd.concat(bin_rows, ignore_index=True) if bin_rows else pd.DataFrame()
    skipped_df = pd.DataFrame(skipped)
    final = integrate_malignant_calls(cnv_calls, cnv_candidates, subclusters, copykat_available, infercnv_available)

    cnv_calls_path = args.metadata_dir / "cnv_proxy_calls_by_cell.tsv.gz"
    sample_summary_path = args.metadata_dir / "cnv_proxy_by_sample.tsv"
    bin_summary_path = args.metadata_dir / "cnv_proxy_bin_summary.tsv.gz"
    skipped_path = args.metadata_dir / "cnv_proxy_skipped_samples.tsv"
    final_path = args.metadata_dir / "malignant_hcc_calls_by_cell.tsv.gz"
    final_sample_path = args.metadata_dir / "malignant_hcc_calls_by_sample.tsv"
    report_path = args.metadata_dir / "malignant_hcc_module3_report.json"
    cnv_calls.to_csv(cnv_calls_path, sep="\t", index=False, compression="gzip")
    sample_summary.to_csv(sample_summary_path, sep="\t", index=False)
    bin_summary.to_csv(bin_summary_path, sep="\t", index=False, compression="gzip")
    skipped_df.to_csv(skipped_path, sep="\t", index=False)
    final.to_csv(final_path, sep="\t", index=False, compression="gzip")
    final_sample = (
        final.groupby(["dataset", "study_sample", "cnv_sample", "sample_source_class", "malignant_hcc_call"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["dataset", "cnv_sample", "malignant_hcc_call"])
    )
    final_sample.to_csv(final_sample_path, sep="\t", index=False)
    report = {
        "method": "hgnc_cytoband_binned_expression_cnv_proxy",
        "method_note": "Fallback expression CNV proxy based on HGNC cytoband chromosome bins; CopyKAT/inferCNV were not available in this R environment.",
        "sample_source_note": "GSE149614 T=tumor, N=normal_adjacent, P=PVTT tumor, L=lymph-node metastatic tumor based on project sample metadata.",
        "copykat_available": copykat_available,
        "infercnv_available": infercnv_available,
        "copykat_status": "not_run_package_unavailable" if not copykat_available else "available_not_run_by_this_script",
        "infercnv_status": "not_run_package_unavailable" if not infercnv_available else "available_not_run_by_this_script",
        "n_cnv_candidate_input": int(cnv_candidates.shape[0]),
        "n_cnv_proxy_called": int(cnv_calls.shape[0]),
        "n_skipped_candidate": int(cnv_candidates.shape[0] - cnv_calls.shape[0]),
        "n_gene_map": int(gene_map.shape[0]),
        "gene_bin_size": args.gene_bin_size,
        "cnv_thresholds": {
            "cnv_z_threshold": args.cnv_z_threshold,
            "cnv_high_bin_fraction": args.cnv_high_bin_fraction,
        },
        "cnv_proxy_status_counts": cnv_calls["cnv_proxy_status"].value_counts(dropna=False).to_dict() if not cnv_calls.empty else {},
        "malignant_hcc_call_counts": final["malignant_hcc_call"].value_counts(dropna=False).to_dict() if not final.empty else {},
        "outputs": {
            "cnv_proxy_calls_by_cell": str(cnv_calls_path.resolve()),
            "cnv_proxy_by_sample": str(sample_summary_path.resolve()),
            "cnv_proxy_bin_summary": str(bin_summary_path.resolve()),
            "skipped_samples": str(skipped_path.resolve()),
            "malignant_hcc_calls_by_cell": str(final_path.resolve()),
            "malignant_hcc_calls_by_sample": str(final_sample_path.resolve()),
            "gene_map": str(args.gene_map.resolve()),
        },
        "elapsed_seconds": round(time.time() - start, 3),
    }
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
