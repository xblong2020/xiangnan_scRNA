from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io, sparse


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import malignant_hcc_module3 as m3  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export module 3 candidate/reference matrices for CopyKAT.")
    parser.add_argument("--integrated", type=Path, default=ROOT / "data/processed/scvi/scvi_integrated_counts_hvg.celltypist_major.h5ad")
    parser.add_argument("--manifest", type=Path, default=ROOT / "metadata/scvi/scvi_input_manifest.counts.tsv")
    parser.add_argument("--major-seed", type=Path, default=ROOT / "metadata/celltype/scanvi_seed_labels_by_cell.tsv.gz")
    parser.add_argument("--cnv-candidate-cells", type=Path, default=ROOT / "metadata/hepatocyte/hepatocyte_cnv_candidate_cells.tsv.gz")
    parser.add_argument("--hepatocyte-cells", type=Path, default=ROOT / "metadata/hepatocyte/hepatocyte_lineage_cells.tsv.gz")
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/copykat_module3")
    parser.add_argument("--out-dir", type=Path, default=ROOT / "data/processed/copykat_module3")
    parser.add_argument("--min-candidate-cells", type=int, default=20)
    parser.add_argument("--min-reference-cells", type=int, default=100)
    parser.add_argument("--max-reference-cells", type=int, default=1000)
    parser.add_argument("--max-candidate-cells-per-run", type=int, default=4000)
    parser.add_argument("--samples", default="", help="Optional comma-separated cnv_sample allowlist.")
    parser.add_argument("--seed", type=int, default=20260601)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def safe_name(value: object) -> str:
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value))
    return text.strip("_") or "sample"


def as_csr(x) -> sparse.csr_matrix:
    if sparse.issparse(x):
        return x.tocsr()
    return sparse.csr_matrix(x)


def source_positions(a: ad.AnnData, original_ids: list[str]) -> tuple[np.ndarray, list[str]]:
    pos = pd.Series(np.arange(a.n_obs), index=a.obs_names.astype(str))
    present = [cell for cell in original_ids if cell in pos.index]
    return pos.loc[present].to_numpy(), present


def write_run(
    a: ad.AnnData,
    rows: pd.DataFrame,
    run_dir: Path,
    run_id: str,
    study_sample: str,
    overwrite: bool,
) -> dict[str, object]:
    matrix_path = run_dir / "raw_counts_gene_by_cell.mtx"
    map_path = run_dir / "cell_map.tsv"
    gene_path = run_dir / "genes.tsv"
    norm_path = run_dir / "normal_cell_keys.txt"
    if matrix_path.exists() and map_path.exists() and gene_path.exists() and norm_path.exists() and not overwrite:
        return {
            "run_id": run_id,
            "run_dir": str(run_dir),
            "status": "exists",
            "n_cells": int(pd.read_csv(map_path, sep="\t").shape[0]),
        }

    run_dir.mkdir(parents=True, exist_ok=True)
    rows = rows.copy().reset_index(drop=True)
    rows["original_cell_id"] = [m3.original_id(c, study_sample) for c in rows["cell_id"].astype(str)]
    pos, present = source_positions(a, rows["original_cell_id"].tolist())
    keep = pd.Index(rows["original_cell_id"]).get_indexer(present)
    rows = rows.iloc[keep].copy().reset_index(drop=True)
    rows["cell_key"] = [f"C{i + 1:06d}" for i in range(rows.shape[0])]
    if rows.empty:
        raise ValueError(f"{run_id}: no cells present in source AnnData")

    x = as_csr(a[pos, :].X).astype(np.float32)
    with matrix_path.open("wb") as handle:
        io.mmwrite(handle, x.T.tocoo())
    pd.DataFrame({"gene": a.var_names.astype(str)}).to_csv(gene_path, sep="\t", index=False)
    rows.to_csv(map_path, sep="\t", index=False)
    norm_keys = rows.loc[rows["cnv_role"].eq("reference"), "cell_key"].astype(str)
    norm_path.write_text("\n".join(norm_keys.tolist()) + "\n", encoding="utf-8")
    return {
        "run_id": run_id,
        "run_dir": str(run_dir),
        "status": "exported",
        "n_cells": int(rows.shape[0]),
        "n_candidate": int(rows["cnv_role"].eq("candidate").sum()),
        "n_reference": int(rows["cnv_role"].eq("reference").sum()),
        "n_genes": int(a.n_vars),
        "matrix_path": str(matrix_path),
        "cell_map": str(map_path),
        "genes": str(gene_path),
        "normal_cell_keys": str(norm_path),
    }


def main() -> int:
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    allow = {x.strip() for x in args.samples.split(",") if x.strip()}

    module_args = argparse.Namespace(
        integrated=args.integrated,
        major_seed=args.major_seed,
        cnv_candidate_cells=args.cnv_candidate_cells,
        hepatocyte_cells=args.hepatocyte_cells,
    )
    candidates, _, refs = m3.build_module3_inputs(module_args)
    manifest = m3.read_manifest(args.manifest)
    manifest_by_study = manifest.set_index("study_sample").to_dict(orient="index")

    manifest_rows = []
    skipped = []
    for (study_sample, cnv_sample), sample_candidates in candidates.groupby(["study_sample", "cnv_sample"], observed=True):
        if allow and str(cnv_sample) not in allow:
            continue
        sample_candidates = sample_candidates.copy()
        if sample_candidates.shape[0] < args.min_candidate_cells:
            skipped.append({"study_sample": study_sample, "cnv_sample": cnv_sample, "reason": "candidate_cells_below_threshold", "n_candidate": sample_candidates.shape[0]})
            continue
        same_refs = refs.loc[refs["study_sample"].eq(study_sample) & refs["cnv_sample"].eq(cnv_sample)].copy()
        if same_refs.shape[0] < args.min_reference_cells:
            same_refs = refs.loc[refs["study_sample"].eq(study_sample)].copy()
        if same_refs.shape[0] < args.min_reference_cells:
            same_refs = refs.loc[refs["dataset"].eq(sample_candidates["dataset"].iloc[0])].copy()
        if same_refs.shape[0] < args.min_reference_cells:
            skipped.append({"study_sample": study_sample, "cnv_sample": cnv_sample, "reason": "reference_cells_below_threshold", "n_candidate": sample_candidates.shape[0], "n_reference": same_refs.shape[0]})
            continue
        if same_refs.shape[0] > args.max_reference_cells:
            same_refs = same_refs.sample(args.max_reference_cells, random_state=args.seed)

        row = manifest_by_study.get(study_sample)
        if row is None:
            skipped.append({"study_sample": study_sample, "cnv_sample": cnv_sample, "reason": "missing_manifest", "n_candidate": sample_candidates.shape[0]})
            continue
        source_path = m3.resolve_path(row["output"])
        a = ad.read_h5ad(source_path, backed="r")
        try:
            cand_indices = np.arange(sample_candidates.shape[0])
            chunks = np.array_split(cand_indices, max(1, int(np.ceil(sample_candidates.shape[0] / args.max_candidate_cells_per_run))))
            for chunk_idx, chunk in enumerate(chunks, start=1):
                chunk_candidates = sample_candidates.iloc[chunk].copy()
                chunk_candidates["cnv_role"] = "candidate"
                ref_rows = same_refs.copy()
                ref_rows["cnv_role"] = "reference"
                run_rows = pd.concat([chunk_candidates, ref_rows], axis=0, ignore_index=True, sort=False)
                run_id = safe_name(f"{cnv_sample}__chunk{chunk_idx:02d}")
                run_dir = args.out_dir / run_id
                info = write_run(a, run_rows, run_dir, run_id, study_sample, args.overwrite)
                info.update(
                    {
                        "dataset": sample_candidates["dataset"].iloc[0],
                        "study_sample": study_sample,
                        "cnv_sample": cnv_sample,
                        "chunk_index": chunk_idx,
                        "source_h5ad": str(source_path),
                    }
                )
                manifest_rows.append(info)
                print(f"EXPORT {run_id} candidates={info.get('n_candidate')} refs={info.get('n_reference')} cells={info.get('n_cells')}", flush=True)
        finally:
            a.file.close()

    out_manifest = pd.DataFrame(manifest_rows)
    out_manifest.to_csv(args.metadata_dir / "copykat_input_manifest.tsv", sep="\t", index=False)
    pd.DataFrame(skipped).to_csv(args.metadata_dir / "copykat_input_skipped.tsv", sep="\t", index=False)
    print(f"WROTE {args.metadata_dir / 'copykat_input_manifest.tsv'} rows={out_manifest.shape[0]}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
