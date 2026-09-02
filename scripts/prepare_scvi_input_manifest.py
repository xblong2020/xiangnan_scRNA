from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
QC_SUMMARY = ROOT / "metadata" / "qc" / "qc_summary.tsv"
QC_INVENTORY = ROOT / "metadata" / "qc" / "qc_final_inventory.tsv"
OUT_DIR = ROOT / "metadata" / "scvi"


MIN_CELLS_EXCLUDE = 1000
MIN_GENES_EXCLUDE = 5000
LOW_CELL_REVIEW = 3000


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary = pd.read_csv(QC_SUMMARY, sep="\t")
    inventory = pd.read_csv(QC_INVENTORY, sep="\t")

    for col in ["input_cells", "kept_cells", "input_genes", "kept_genes"]:
        if col in summary.columns:
            summary[col] = pd.to_numeric(summary[col], errors="coerce")
    summary["retention_rate"] = summary["kept_cells"] / summary["input_cells"]

    exclusion_reasons: list[str] = []
    review_flags: list[str] = []
    for _, row in summary.iterrows():
        reasons: list[str] = []
        flags: list[str] = []
        if row["kept_cells"] < MIN_CELLS_EXCLUDE:
            reasons.append(f"kept_cells<{MIN_CELLS_EXCLUDE}")
        if row["kept_genes"] < MIN_GENES_EXCLUDE:
            reasons.append(f"kept_genes<{MIN_GENES_EXCLUDE}")
        if not reasons and row["kept_cells"] < LOW_CELL_REVIEW:
            flags.append(f"low_cells_review<{LOW_CELL_REVIEW}")
        exclusion_reasons.append(";".join(reasons))
        review_flags.append(";".join(flags))

    summary["exclude_from_scvi"] = [bool(x) for x in exclusion_reasons]
    summary["exclude_reason"] = exclusion_reasons
    summary["review_flag"] = review_flags

    exclusions = summary.loc[summary["exclude_from_scvi"]].copy()
    review = summary.loc[(~summary["exclude_from_scvi"]) & (summary["review_flag"] != "")].copy()

    curated = inventory.merge(
        summary[["dataset", "label", "exclude_from_scvi", "exclude_reason", "review_flag"]],
        on=["dataset", "label"],
        how="left",
    )
    curated["include_in_scvi"] = ~curated["exclude_from_scvi"].fillna(False)
    curated = curated.sort_values(["include_in_scvi", "dataset", "label"], ascending=[False, True, True])

    summary.to_csv(OUT_DIR / "sample_qc_flags.tsv", sep="\t", index=False)
    exclusions.to_csv(OUT_DIR / "excluded_samples.tsv", sep="\t", index=False)
    review.to_csv(OUT_DIR / "review_samples.tsv", sep="\t", index=False)
    curated.to_csv(OUT_DIR / "scvi_input_manifest.tsv", sep="\t", index=False)

    print(f"WROTE {OUT_DIR / 'sample_qc_flags.tsv'}")
    print(f"WROTE {OUT_DIR / 'excluded_samples.tsv'}")
    print(f"WROTE {OUT_DIR / 'review_samples.tsv'}")
    print(f"WROTE {OUT_DIR / 'scvi_input_manifest.tsv'}")
    print("Excluded samples:")
    if exclusions.empty:
        print("  none")
    else:
        for _, row in exclusions.iterrows():
            print(f"  {row['dataset']}:{row['label']} - {row['exclude_reason']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
