from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
QC_ROOT = ROOT / "data" / "processed" / "qc_h5ad"
OUT_DIR = ROOT / "metadata" / "scvi"

MIN_CELLS_EXCLUDE = 1000
MIN_GENES_EXCLUDE = 5000
LOW_CELL_REVIEW = 3000
MAX_NON_INTEGER_RATE = 0.001


def matrix_check(path: Path) -> dict[str, object]:
    adata = ad.read_h5ad(path)
    x = adata.X
    gene_names = pd.Index(adata.var_names.astype(str))
    numeric_gene_rate = float(gene_names.str.fullmatch(r"\d+").mean()) if len(gene_names) else 0.0
    values = x.data if sparse.issparse(x) else np.ravel(np.asarray(x))
    sample = values[: min(values.size, 100_000)]
    non_integer_rate = 0.0
    min_value = 0.0
    max_value = 0.0
    if sample.size:
        non_integer_rate = float(np.mean(np.abs(sample - np.round(sample)) > 1e-6))
        min_value = float(sample.min())
        max_value = float(sample.max())
    return {
        "kept_cells": int(adata.n_obs),
        "kept_genes": int(adata.n_vars),
        "x_dtype": str(x.dtype),
        "x_non_integer_rate_sample": non_integer_rate,
        "x_min_sample": min_value,
        "x_max_sample": max_value,
        "layers": ",".join(adata.layers.keys()),
        "numeric_gene_name_rate": numeric_gene_rate,
    }


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for path in sorted(QC_ROOT.glob("*/*.qc.h5ad")):
        dataset = path.parent.name
        label = path.name.removesuffix(".qc.h5ad")
        row = {
            "dataset": dataset,
            "label": label,
            "output": str(path.resolve()),
            "bytes": path.stat().st_size,
        }
        row.update(matrix_check(path))
        rows.append(row)

    manifest = pd.DataFrame(rows)
    exclusion_reasons: list[str] = []
    review_flags: list[str] = []
    for _, row in manifest.iterrows():
        reasons: list[str] = []
        flags: list[str] = []
        if row["kept_cells"] < MIN_CELLS_EXCLUDE:
            reasons.append(f"kept_cells<{MIN_CELLS_EXCLUDE}")
        if row["kept_genes"] < MIN_GENES_EXCLUDE:
            reasons.append(f"kept_genes<{MIN_GENES_EXCLUDE}")
        if row["x_non_integer_rate_sample"] > MAX_NON_INTEGER_RATE:
            reasons.append("not_raw_counts_matrix")
        if row["numeric_gene_name_rate"] > 0.5:
            reasons.append("gene_names_are_numeric_indices")
        if not reasons and row["kept_cells"] < LOW_CELL_REVIEW:
            flags.append(f"low_cells_review<{LOW_CELL_REVIEW}")
        exclusion_reasons.append(";".join(reasons))
        review_flags.append(";".join(flags))

    manifest["exclude_from_scvi"] = [bool(reason) for reason in exclusion_reasons]
    manifest["exclude_reason"] = exclusion_reasons
    manifest["review_flag"] = review_flags
    manifest["include_in_scvi"] = ~manifest["exclude_from_scvi"]
    manifest = manifest.sort_values(["include_in_scvi", "dataset", "label"], ascending=[False, True, True])

    matrix_check_out = OUT_DIR / "scvi_counts_matrix_check.tsv"
    manifest_out = OUT_DIR / "scvi_input_manifest.counts.tsv"
    excluded_out = OUT_DIR / "excluded_samples.counts.tsv"
    review_out = OUT_DIR / "review_samples.counts.tsv"

    manifest.to_csv(matrix_check_out, sep="\t", index=False)
    manifest.to_csv(manifest_out, sep="\t", index=False)
    manifest.loc[manifest["exclude_from_scvi"]].to_csv(excluded_out, sep="\t", index=False)
    manifest.loc[(~manifest["exclude_from_scvi"]) & (manifest["review_flag"] != "")].to_csv(
        review_out, sep="\t", index=False
    )

    print(f"WROTE {matrix_check_out}")
    print(f"WROTE {manifest_out}")
    print(f"WROTE {excluded_out}")
    print(f"WROTE {review_out}")
    print("Included samples:", int(manifest["include_in_scvi"].sum()))
    print("Excluded samples:", int((~manifest["include_in_scvi"]).sum()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
