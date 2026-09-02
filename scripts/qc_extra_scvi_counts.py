from __future__ import annotations

from pathlib import Path

import pandas as pd

from qc_unify_h5ad import META_ROOT, QCTask, qc_one


ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    tasks = [
        QCTask(
            "GSE202379",
            "GSE202379_SeuratObject_AllCells_counts_fixed",
            ROOT / "data/processed/h5ad_from_seurat/GSE202379/GSE202379_SeuratObject_AllCells.h5ad",
            use_raw=False,
        ),
        QCTask(
            "GSE174748",
            "GSE174748_hl_nuclei_counts",
            ROOT / "data/processed/h5ad_counts_from_seurat/GSE174748/GSE174748_hl_nuclei_counts.h5ad",
            use_raw=False,
        )
    ]
    rows = [qc_one(task) for task in tasks]
    report = pd.DataFrame(rows)
    META_ROOT.mkdir(parents=True, exist_ok=True)
    out = META_ROOT / "qc_extra_scvi_counts.tsv"
    report.to_csv(out, sep="\t", index=False)
    print(f"WROTE {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
