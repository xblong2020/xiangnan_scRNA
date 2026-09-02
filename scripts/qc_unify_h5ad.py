from __future__ import annotations

import gc
import re
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
OUT_ROOT = ROOT / "data" / "processed" / "qc_h5ad"
META_ROOT = ROOT / "metadata" / "qc"


@dataclass(frozen=True)
class QCTask:
    dataset: str
    label: str
    path: Path
    use_raw: bool = False
    min_genes: int = 200
    min_counts: int = 500
    max_mito_pct: float = 25.0
    min_cells_per_gene: int = 3
    upper_quantile: float = 0.995


def tasks() -> list[QCTask]:
    return [
        QCTask(
            "GSE202379",
            "GSE202379_SeuratObject_AllCells",
            ROOT / "data/processed/h5ad_from_seurat/GSE202379/GSE202379_SeuratObject_AllCells.h5ad",
            use_raw=True,
        ),
        QCTask(
            "GSE174748",
            "GSE174748_hl_nuclei",
            ROOT / "data/processed/h5ad_from_seurat/GSE174748/GSE174748_hl_nuclei.h5ad",
            use_raw=True,
        ),
        QCTask(
            "HCC_atlas",
            "HCC_atlas_all_release",
            ROOT / "data/processed/h5ad_from_seurat/HCC_atlas/HCC_atlas_all_release.h5ad",
            use_raw=True,
        ),
        QCTask(
            "GSE149614",
            "GSE149614_raw_counts",
            ROOT / "data/processed/h5ad/GSE149614/GSE149614_raw_counts.h5ad",
        ),
        QCTask(
            "GSE151530",
            "GSE151530",
            ROOT / "data/processed/h5ad/GSE151530/GSE151530.h5ad",
        ),
        *[
            QCTask("GSE185477", p.stem, p)
            for p in sorted((ROOT / "data/processed/h5ad/GSE185477").glob("*.h5ad"))
        ],
        *[
            QCTask("GSE212046", p.stem, p)
            for p in sorted((ROOT / "data/processed/h5ad/GSE212046").glob("*.h5ad"))
        ],
    ]


def clean_gene_symbol(value: object) -> str:
    symbol = str(value).strip()
    symbol = re.sub(r"\s+", "_", symbol)
    return symbol.upper()


def make_obs_metrics(adata: ad.AnnData, x: sparse.spmatrix, gene_symbols: np.ndarray) -> pd.DataFrame:
    obs = adata.obs.copy()
    total_counts = np.asarray(x.sum(axis=1)).ravel()
    n_genes = np.asarray((x > 0).sum(axis=1)).ravel()
    mt_mask = np.char.startswith(gene_symbols.astype(str), "MT-")
    if mt_mask.any():
        mt_counts = np.asarray(x[:, mt_mask].sum(axis=1)).ravel()
        pct_mt = np.divide(mt_counts, total_counts, out=np.zeros_like(mt_counts, dtype=float), where=total_counts > 0) * 100
    else:
        pct_mt = np.zeros(x.shape[0], dtype=float)

    obs["qc_total_counts"] = total_counts
    obs["qc_n_genes_by_counts"] = n_genes
    obs["qc_pct_counts_mt"] = pct_mt
    return obs


def aggregate_duplicate_genes(x: sparse.spmatrix, gene_symbols: np.ndarray) -> tuple[sparse.csr_matrix, np.ndarray]:
    codes, uniques = pd.factorize(gene_symbols, sort=True)
    mapper = sparse.csr_matrix(
        (np.ones(len(codes), dtype=np.float32), (np.arange(len(codes)), codes)),
        shape=(len(codes), len(uniques)),
    )
    return (x @ mapper).tocsr(), uniques.astype(str)


def qc_one(task: QCTask) -> dict[str, object]:
    out_dir = OUT_ROOT / task.dataset
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{task.label}.qc.h5ad"
    if out_path.exists():
        print(f"SKIP {out_path}", flush=True)
        a = ad.read_h5ad(out_path, backed="r")
        result = {
            "dataset": task.dataset,
            "label": task.label,
            "source": str(task.path),
            "output": str(out_path),
            "input_cells": "",
            "input_genes": "",
            "kept_cells": a.n_obs,
            "kept_genes": a.n_vars,
            "status": "skipped_existing",
        }
        a.file.close()
        return result

    print(f"READ {task.dataset}/{task.label}", flush=True)
    source = ad.read_h5ad(task.path)
    if task.use_raw and source.raw is not None:
        adata = source.raw.to_adata()
        adata.obs = source.obs.copy()
    else:
        adata = source

    x = adata.X
    if not sparse.issparse(x):
        x = sparse.csr_matrix(x)
    else:
        x = x.tocsr()

    gene_symbols = np.array([clean_gene_symbol(v) for v in adata.var_names], dtype=str)
    valid_gene = np.array([g not in {"", "NAN", "NONE"} for g in gene_symbols], dtype=bool)
    x = x[:, valid_gene]
    gene_symbols = gene_symbols[valid_gene]

    x, unified_genes = aggregate_duplicate_genes(x, gene_symbols)
    obs = make_obs_metrics(adata, x, unified_genes)

    total = obs["qc_total_counts"].to_numpy(dtype=float)
    n_genes = obs["qc_n_genes_by_counts"].to_numpy(dtype=float)
    pct_mt = obs["qc_pct_counts_mt"].to_numpy(dtype=float)
    max_counts = np.quantile(total[total > 0], task.upper_quantile) if np.any(total > 0) else 0
    max_genes = np.quantile(n_genes[n_genes > 0], task.upper_quantile) if np.any(n_genes > 0) else 0

    keep_cells = (
        (total >= task.min_counts)
        & (n_genes >= task.min_genes)
        & (total <= max_counts)
        & (n_genes <= max_genes)
        & (pct_mt <= task.max_mito_pct)
    )
    x = x[keep_cells, :].tocsr()
    obs = obs.loc[keep_cells].copy()
    gene_cells = np.asarray((x > 0).sum(axis=0)).ravel()
    keep_genes = gene_cells >= task.min_cells_per_gene
    x = x[:, keep_genes].tocsr()
    kept_genes = unified_genes[keep_genes]

    var = pd.DataFrame(index=kept_genes)
    var["gene_symbol"] = kept_genes
    qc = ad.AnnData(X=x, obs=obs, var=var)
    qc.obs_names_make_unique()
    qc.var_names_make_unique()
    qc.uns["qc_thresholds"] = {
        "min_genes": task.min_genes,
        "min_counts": task.min_counts,
        "max_mito_pct": task.max_mito_pct,
        "min_cells_per_gene": task.min_cells_per_gene,
        "upper_quantile": task.upper_quantile,
        "max_counts": float(max_counts),
        "max_genes": float(max_genes),
        "gene_symbol_rule": "strip whitespace, replace internal whitespace with underscore, uppercase, aggregate duplicates by sum",
    }
    qc.write_h5ad(out_path, compression="gzip")
    print(
        f"WROTE {out_path} cells {adata.n_obs}->{qc.n_obs}, genes {adata.n_vars}->{qc.n_vars}",
        flush=True,
    )
    result = {
        "dataset": task.dataset,
        "label": task.label,
        "source": str(task.path),
        "output": str(out_path),
        "input_cells": adata.n_obs,
        "input_genes": adata.n_vars,
        "kept_cells": qc.n_obs,
        "kept_genes": qc.n_vars,
        "min_counts": task.min_counts,
        "min_genes": task.min_genes,
        "max_counts": float(max_counts),
        "max_genes": float(max_genes),
        "max_mito_pct": task.max_mito_pct,
        "status": "complete",
    }
    del source, adata, qc, x
    gc.collect()
    return result


def main() -> int:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    META_ROOT.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for task in tasks():
        if not task.path.exists():
            rows.append(
                {
                    "dataset": task.dataset,
                    "label": task.label,
                    "source": str(task.path),
                    "status": "missing_source",
                }
            )
            continue
        rows.append(qc_one(task))
    report = pd.DataFrame(rows)
    report.to_csv(META_ROOT / "qc_summary.tsv", sep="\t", index=False)
    print(f"WROTE {META_ROOT / 'qc_summary.tsv'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
