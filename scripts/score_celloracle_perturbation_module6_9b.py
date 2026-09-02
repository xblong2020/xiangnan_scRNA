from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_H5AD = PROJECT_ROOT / "data/processed/driver/celloracle_module6_6/celloracle_module6_6_input.h5ad"
DEFAULT_METADATA_DIR = PROJECT_ROOT / "metadata/driver"

FATE_COL = "cellrank_fate_prob_cnv_supported_malignant"
PSEUDOTIME_COL = "driver_main_strict__pseudotime_mean"
HCC_MODULE_COL = "driver_main_strict__module_HCC_Malignant_Associated"
PROLIFERATION_MODULE_COL = "driver_main_strict__module_Proliferation"
MATURE_MODULE_COL = "driver_main_strict__module_Mature_Hepatocyte"

SCORE_COLUMNS = [
    "malignant_fate_direction_score",
    "inner_product_score",
    "cnv_fate_probability_association_score",
    "module_rescue_score",
    "state_specificity_score",
]


def minmax_scale(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce").astype(float)
    valid = numeric.dropna()
    if valid.empty:
        return pd.Series(np.zeros(len(numeric)), index=numeric.index)
    missing = numeric.isna()
    vmin = valid.min()
    vmax = valid.max()
    if vmax == vmin:
        scaled = pd.Series(np.full(len(numeric), 0.5), index=numeric.index)
    else:
        scaled = (numeric - vmin) / (vmax - vmin)
    scaled.loc[missing] = 0.0
    return scaled.fillna(0.0)


def _safe_numeric(values: pd.Series) -> np.ndarray:
    numeric = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    if np.isnan(numeric).all():
        return np.zeros_like(numeric)
    median = np.nanmedian(numeric)
    numeric[np.isnan(numeric)] = median
    return numeric


def compute_local_gradients(embedding: np.ndarray, values: np.ndarray, k: int = 50) -> np.ndarray:
    if embedding.ndim != 2 or embedding.shape[1] != 2:
        raise ValueError("Embedding must be a two-column matrix")
    if len(values) != embedding.shape[0]:
        raise ValueError("Values and embedding must have the same number of cells")
    k_eff = min(max(3, k), embedding.shape[0] - 1)
    nn = NearestNeighbors(n_neighbors=k_eff + 1)
    nn.fit(embedding)
    _, indices = nn.kneighbors(embedding)
    neighbor_indices = indices[:, 1:]
    gradients = np.zeros_like(embedding, dtype=float)
    values = np.asarray(values, dtype=float)
    for cell_ix, neigh_ixs in enumerate(neighbor_indices):
        xdiff = embedding[neigh_ixs] - embedding[cell_ix]
        ydiff = values[neigh_ixs] - values[cell_ix]
        gradients[cell_ix] = np.linalg.lstsq(xdiff, ydiff, rcond=None)[0]
    return gradients


def load_cell_metadata_and_gradients(h5ad_path: Path, gradient_k: int) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    adata = ad.read_h5ad(h5ad_path)
    if "X_celloracle_umap" not in adata.obsm:
        raise ValueError("X_celloracle_umap is required for Module 6.9b scoring")
    required = [FATE_COL, PSEUDOTIME_COL, HCC_MODULE_COL, PROLIFERATION_MODULE_COL, MATURE_MODULE_COL, "celloracle_state"]
    missing = [col for col in required if col not in adata.obs.columns]
    if missing:
        raise ValueError(f"Missing required h5ad obs columns: {missing}")

    embedding = np.asarray(adata.obsm["X_celloracle_umap"], dtype=float)
    metadata = pd.DataFrame(
        {
            "cell_id": adata.obs_names.astype(str),
            "celloracle_state": adata.obs["celloracle_state"].astype(str).to_numpy(),
            "umap_1": embedding[:, 0],
            "umap_2": embedding[:, 1],
            FATE_COL: _safe_numeric(adata.obs[FATE_COL]),
            PSEUDOTIME_COL: _safe_numeric(adata.obs[PSEUDOTIME_COL]),
            HCC_MODULE_COL: _safe_numeric(adata.obs[HCC_MODULE_COL]),
            PROLIFERATION_MODULE_COL: _safe_numeric(adata.obs[PROLIFERATION_MODULE_COL]),
            MATURE_MODULE_COL: _safe_numeric(adata.obs[MATURE_MODULE_COL]),
        }
    )
    gradient_specs = {
        "cnv_fate": FATE_COL,
        "pseudotime": PSEUDOTIME_COL,
        "hcc_module": HCC_MODULE_COL,
        "proliferation_module": PROLIFERATION_MODULE_COL,
        "mature_module": MATURE_MODULE_COL,
    }
    gradients = {
        name: compute_local_gradients(embedding, metadata[col].to_numpy(dtype=float), k=gradient_k)
        for name, col in gradient_specs.items()
    }
    return metadata, gradients


def score_cell_level_vectors(
    cell_shift: pd.DataFrame,
    metadata: pd.DataFrame,
    gradients: dict[str, np.ndarray],
) -> pd.DataFrame:
    metadata = metadata.reset_index(drop=True).copy()
    metadata["cell_order_ix"] = np.arange(len(metadata))
    merge_keys = ["cell_id", "celloracle_state"] if "celloracle_state" in metadata.columns else ["cell_id"]
    merged = cell_shift.merge(metadata, on=merge_keys, how="left", validate="many_to_one")
    if merged["cell_order_ix"].isna().any():
        missing = int(merged["cell_order_ix"].isna().sum())
        raise ValueError(f"{missing} perturbation rows could not be aligned to h5ad cell metadata")
    cell_ix = merged["cell_order_ix"].to_numpy(dtype=int)
    delta = merged[["delta_embedding_1", "delta_embedding_2"]].to_numpy(dtype=float)

    merged["malignant_fate_direction_cell_score"] = -pd.to_numeric(
        merged["malignant_axis_projection"],
        errors="coerce",
    ).astype(float)
    merged["inner_product_raw_pseudotime"] = np.einsum("ij,ij->i", delta, gradients["pseudotime"][cell_ix])
    merged["inner_product_score"] = -merged["inner_product_raw_pseudotime"]
    merged["inner_product_cell_score"] = merged["inner_product_score"]
    merged["cnv_fate_inner_product_raw"] = np.einsum("ij,ij->i", delta, gradients["cnv_fate"][cell_ix])
    merged["cnv_fate_probability_association_cell_score"] = -merged["cnv_fate_inner_product_raw"]
    merged["hcc_module_inner_product_raw"] = np.einsum("ij,ij->i", delta, gradients["hcc_module"][cell_ix])
    merged["proliferation_module_inner_product_raw"] = np.einsum(
        "ij,ij->i",
        delta,
        gradients["proliferation_module"][cell_ix],
    )
    merged["mature_module_inner_product_raw"] = np.einsum("ij,ij->i", delta, gradients["mature_module"][cell_ix])
    merged["hcc_malignant_module_rescue_cell_score"] = -merged["hcc_module_inner_product_raw"]
    merged["proliferation_module_rescue_cell_score"] = -merged["proliferation_module_inner_product_raw"]
    merged["mature_hepatocyte_module_rescue_cell_score"] = merged["mature_module_inner_product_raw"]
    merged["module_rescue_cell_score"] = (
        merged["hcc_malignant_module_rescue_cell_score"]
        + merged["proliferation_module_rescue_cell_score"]
        + merged["mature_hepatocyte_module_rescue_cell_score"]
    ) / 3.0
    return merged.drop(columns=["cell_order_ix"])


def _weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    values_arr = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    weights_arr = pd.to_numeric(weights, errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(values_arr) & np.isfinite(weights_arr)
    if not valid.any():
        return float("nan")
    weights_arr = np.clip(weights_arr[valid], 0, None)
    if weights_arr.sum() == 0:
        return float(np.nanmean(values_arr[valid]))
    return float(np.average(values_arr[valid], weights=weights_arr))


def compute_state_specificity(cell_scores: pd.DataFrame, eps: float = 1e-9) -> pd.DataFrame:
    rows = []
    for tf, df in cell_scores.groupby("tf", sort=False):
        malignant = df.loc[df["celloracle_state"].astype(str) == "malignant_or_malignant_like"]
        other = df.loc[df["celloracle_state"].astype(str) != "malignant_or_malignant_like"]
        malignant_effect = float(malignant["malignant_fate_direction_cell_score"].mean()) if len(malignant) else 0.0
        malignant_abs = float(malignant["embedding_shift_norm"].mean()) if len(malignant) else 0.0
        other_abs = float(other["embedding_shift_norm"].mean()) if len(other) else 0.0
        specificity_ratio = malignant_abs / (other_abs + eps)
        malignant_specific_rescue = malignant_effect * specificity_ratio
        rows.append(
            {
                "tf": tf,
                "malignant_state_direction_mean": malignant_effect,
                "malignant_state_embedding_shift_mean": malignant_abs,
                "non_malignant_embedding_shift_mean": other_abs,
                "state_specificity_ratio": specificity_ratio,
                "malignant_specific_rescue_score": malignant_specific_rescue,
                "state_specificity_score": malignant_specific_rescue,
            }
        )
    return pd.DataFrame(rows)


def aggregate_tf_scores(cell_scores: pd.DataFrame, fate_high_threshold: float) -> pd.DataFrame:
    specificity = compute_state_specificity(cell_scores)
    rows = []
    for tf, df in cell_scores.groupby("tf", sort=False):
        high_fate = df.loc[pd.to_numeric(df[FATE_COL], errors="coerce") >= fate_high_threshold]
        malignant = df.loc[df["celloracle_state"].astype(str) == "malignant_or_malignant_like"]
        focus = high_fate if len(high_fate) else malignant
        if not len(focus):
            focus = df
        weights = focus[FATE_COL].fillna(0.0) + 0.05
        rows.append(
            {
                "tf": tf,
                "n_cells_total": int(df["cell_id"].nunique()) if "cell_id" in df.columns else int(len(df)),
                "n_high_cnv_fate_cells": int(high_fate["cell_id"].nunique()) if "cell_id" in df.columns else int(len(high_fate)),
                "n_malignant_like_cells": int(malignant["cell_id"].nunique()) if "cell_id" in df.columns else int(len(malignant)),
                "malignant_fate_direction_score": _weighted_mean(
                    focus["malignant_fate_direction_cell_score"],
                    weights,
                ),
                "inner_product_score": _weighted_mean(focus["inner_product_cell_score"], weights),
                "cnv_fate_probability_association_score": _weighted_mean(
                    focus["cnv_fate_probability_association_cell_score"],
                    weights,
                ),
                "hcc_malignant_module_rescue_score": _weighted_mean(
                    focus["hcc_malignant_module_rescue_cell_score"],
                    weights,
                ),
                "proliferation_module_rescue_score": _weighted_mean(
                    focus["proliferation_module_rescue_cell_score"],
                    weights,
                ),
                "mature_hepatocyte_module_rescue_score": _weighted_mean(
                    focus["mature_hepatocyte_module_rescue_cell_score"],
                    weights,
                ),
                "module_rescue_score": _weighted_mean(focus["module_rescue_cell_score"], weights),
            }
        )
    summary = pd.DataFrame(rows).merge(specificity, on="tf", how="left")
    for col in SCORE_COLUMNS:
        summary[f"{col}_scaled"] = minmax_scale(summary[col])
    summary["quantitative_perturbation_score"] = summary[[f"{col}_scaled" for col in SCORE_COLUMNS]].mean(axis=1)
    summary = summary.sort_values(
        ["quantitative_perturbation_score", "malignant_fate_direction_score"],
        ascending=[False, False],
    ).reset_index(drop=True)
    summary.insert(0, "quantitative_rank", np.arange(1, len(summary) + 1))
    return summary


def run_module6_9b(
    h5ad_path: Path,
    metadata_dir: Path,
    gradient_k: int,
    fate_high_threshold: float,
) -> dict:
    metadata_dir.mkdir(parents=True, exist_ok=True)
    cell_shift = pd.read_csv(metadata_dir / "celloracle_module6_8_cell_shift_summary.tsv.gz", sep="\t")
    metadata, gradients = load_cell_metadata_and_gradients(h5ad_path, gradient_k=gradient_k)
    cell_scores = score_cell_level_vectors(cell_shift, metadata, gradients)
    tf_scores = aggregate_tf_scores(cell_scores, fate_high_threshold=fate_high_threshold)

    existing_evidence_path = metadata_dir / "celloracle_module6_9_candidate_evidence_matrix.tsv"
    updated_evidence_path = metadata_dir / "celloracle_module6_9b_candidate_evidence_matrix_with_quant_scores.tsv"
    if existing_evidence_path.exists():
        evidence = pd.read_csv(existing_evidence_path, sep="\t")
        merged_evidence = evidence.merge(tf_scores, on="tf", how="left", suffixes=("", "_6_9b"))
        merged_evidence.to_csv(updated_evidence_path, sep="\t", index=False)
    else:
        updated_evidence_path = None

    cell_scores_path = metadata_dir / "celloracle_module6_9b_cell_level_scores.tsv.gz"
    tf_scores_path = metadata_dir / "celloracle_module6_9b_quantitative_tf_scores.tsv"
    definitions_path = metadata_dir / "celloracle_module6_9b_score_definitions.md"
    cell_scores.to_csv(cell_scores_path, sep="\t", index=False)
    tf_scores.to_csv(tf_scores_path, sep="\t", index=False)
    definitions_path.write_text(
        "\n".join(
            [
                "# Module 6.9b Quantitative Perturbation Scores",
                "",
                "All scores are oriented so that larger values indicate stronger predicted rescue away from CNV-supported malignant fate.",
                "",
                "- malignant_fate_direction_score: negative projection of perturbation vectors onto the global normal-to-malignant CellOracle UMAP axis, evaluated in high CNV-fate or malignant-like cells.",
                "- inner_product_score: negative inner product between perturbation vectors and the local pseudotime gradient.",
                "- cnv_fate_probability_association_score: negative inner product between perturbation vectors and the local CNV malignant fate probability gradient.",
                "- module_rescue_score: mean of decreased HCC malignant module direction, decreased proliferation module direction, and increased mature hepatocyte module direction.",
                "- state_specificity_ratio: malignant-like embedding shift magnitude divided by non-malignant embedding shift magnitude.",
                "- state_specificity_score: malignant-like anti-malignant direction score multiplied by state_specificity_ratio.",
                "",
                "Local gradients are estimated by least-squares regression over UMAP nearest neighbors.",
            ]
        ),
        encoding="utf-8",
    )

    return {
        "cell_level_scores": str(cell_scores_path),
        "quantitative_tf_scores": str(tf_scores_path),
        "updated_candidate_evidence_matrix": str(updated_evidence_path) if updated_evidence_path else None,
        "score_definitions": str(definitions_path),
        "n_tfs": int(tf_scores["tf"].nunique()),
        "n_cell_score_rows": int(len(cell_scores)),
        "top_quantitative_tfs": tf_scores.head(5)["tf"].astype(str).tolist(),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 6.9b quantitative CellOracle perturbation scoring")
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA_DIR)
    parser.add_argument("--gradient-k", type=int, default=50)
    parser.add_argument("--fate-high-threshold", type=float, default=0.5)
    parser.add_argument("--report", type=Path, default=DEFAULT_METADATA_DIR / "celloracle_module6_9b_report.json")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = datetime.now(timezone.utc)
    result = run_module6_9b(
        h5ad_path=args.h5ad,
        metadata_dir=args.metadata_dir,
        gradient_k=args.gradient_k,
        fate_high_threshold=args.fate_high_threshold,
    )
    finished = datetime.now(timezone.utc)
    report = {
        "module": "6.9b",
        "method": "Quantitative CellOracle perturbation effect scoring",
        "created_at_utc": finished.isoformat(),
        "elapsed_seconds": (finished - started).total_seconds(),
        "parameters": {
            "h5ad": str(args.h5ad),
            "metadata_dir": str(args.metadata_dir),
            "gradient_k": args.gradient_k,
            "fate_high_threshold": args.fate_high_threshold,
        },
        "result": result,
        "python_runtime": {
            "executable": sys.executable,
            "version": sys.version,
            "platform": platform.platform(),
        },
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(
        {
            "n_tfs": result["n_tfs"],
            "n_cell_score_rows": result["n_cell_score_rows"],
            "top_quantitative_tfs": result["top_quantitative_tfs"],
            "report": str(args.report),
        },
        ensure_ascii=False,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
