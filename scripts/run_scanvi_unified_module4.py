from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scvi
import torch


ROOT = Path(__file__).resolve().parents[1]
UNLABELED = "Unknown"
MAJOR_SEED_LABELS = {
    "T_NK",
    "Myeloid",
    "Endothelial",
    "Fibroblast_HSC_Pericyte",
    "B_cell",
    "Plasma_cell",
    "Cholangiocyte_Progenitor",
}
HEPATOCYTE_STATE_SEED_LABELS = {
    "normal_hepatocyte_like",
    "stressed_injured_hepatocyte",
    "regenerative_progenitor_like_hepatocyte",
}
MALIGNANT_SEED_CALL = "malignant_hcc_high_conf"
MALIGNANT_LABEL = "malignant_hepatocyte"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 4: scANVI unified label propagation.")
    parser.add_argument(
        "--input-h5ad",
        type=Path,
        default=ROOT / "data/processed/scvi/scvi_integrated_counts_hvg.celltypist_major.h5ad",
    )
    parser.add_argument(
        "--scvi-model-dir",
        type=Path,
        default=ROOT / "data/processed/scvi/model_scvi_global_counts",
    )
    parser.add_argument(
        "--major-seed",
        type=Path,
        default=ROOT / "metadata/celltype/scanvi_seed_labels_by_cell.tsv.gz",
    )
    parser.add_argument(
        "--hepatocyte-seed",
        type=Path,
        default=ROOT / "metadata/hepatocyte/hepatocyte_state_seed_by_cell.tsv.gz",
    )
    parser.add_argument(
        "--malignant-calls",
        type=Path,
        default=ROOT / "metadata/malignant/malignant_hcc_calls_by_cell.tsv.gz",
    )
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/scanvi")
    parser.add_argument("--processed-dir", type=Path, default=ROOT / "data/processed/scanvi")
    parser.add_argument("--max-epochs", type=int, default=20)
    parser.add_argument("--batch-size", type=int, default=1024)
    parser.add_argument("--n-samples-per-label", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260601)
    return parser.parse_args()


def bool_series(values: pd.Series) -> pd.Series:
    return values.fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})


def build_seed_table(args: argparse.Namespace, obs_names: pd.Index) -> pd.DataFrame:
    major = pd.read_csv(args.major_seed, sep="\t")
    hepatocyte = pd.read_csv(args.hepatocyte_seed, sep="\t")
    malignant = pd.read_csv(
        args.malignant_calls,
        sep="\t",
        usecols=[
            "cell_id",
            "malignant_hcc_call",
            "cnv_proxy_status",
            "sample_source_class",
            "malignant_hcc_evidence",
        ],
    )

    seed = pd.DataFrame({"cell_id": obs_names.astype(str)})
    seed = seed.merge(
        major[
            [
                "cell_id",
                "leiden_scvi",
                "manual_major_label_cluster",
                "manual_confidence_status",
                "scanvi_seed_label_major",
                "hepatocyte_lineage_candidate",
                "excluded_doublet_cluster",
            ]
        ],
        on="cell_id",
        how="left",
    )
    seed["scanvi_unified_seed_label"] = UNLABELED
    seed["scanvi_unified_seed_source"] = "unlabeled"

    major_mask = (
        seed["manual_confidence_status"].eq("high_conf")
        & seed["scanvi_seed_label_major"].isin(MAJOR_SEED_LABELS)
        & ~bool_series(seed["excluded_doublet_cluster"])
    )
    seed.loc[major_mask, "scanvi_unified_seed_label"] = seed.loc[major_mask, "scanvi_seed_label_major"].astype(str)
    seed.loc[major_mask, "scanvi_unified_seed_source"] = "module1_major_high_conf"

    hepatocyte_labeled = hepatocyte.loc[
        hepatocyte["hepatocyte_state_seed_label"].isin(HEPATOCYTE_STATE_SEED_LABELS),
        ["cell_id", "hepatocyte_state_label", "hepatocyte_state_confidence", "hepatocyte_state_seed_label"],
    ].copy()
    seed = seed.merge(hepatocyte_labeled, on="cell_id", how="left")
    hep_mask = seed["hepatocyte_state_seed_label"].isin(HEPATOCYTE_STATE_SEED_LABELS)
    seed.loc[hep_mask, "scanvi_unified_seed_label"] = seed.loc[hep_mask, "hepatocyte_state_seed_label"].astype(str)
    seed.loc[hep_mask, "scanvi_unified_seed_source"] = "module2_hepatocyte_state_high_conf"

    malignant = malignant.drop_duplicates("cell_id")
    seed = seed.merge(malignant, on="cell_id", how="left")
    malignant_mask = seed["malignant_hcc_call"].eq(MALIGNANT_SEED_CALL)
    seed.loc[malignant_mask, "scanvi_unified_seed_label"] = MALIGNANT_LABEL
    seed.loc[malignant_mask, "scanvi_unified_seed_source"] = "module3_malignant_high_conf"

    seed["scanvi_unified_seed_label"] = seed["scanvi_unified_seed_label"].fillna(UNLABELED).astype(str)
    seed["scanvi_unified_seed_source"] = seed["scanvi_unified_seed_source"].fillna("unlabeled").astype(str)
    return seed


def history_tail(model: scvi.model.SCANVI) -> dict[str, float | int | str]:
    out: dict[str, float | int | str] = {}
    history = getattr(model, "history", None)
    if history is None:
        return out
    for key, values in history.items():
        try:
            last = values.dropna().iloc[-1]
        except Exception:
            continue
        if hasattr(last, "iloc"):
            last = last.iloc[0]
        try:
            out[f"final_{key}"] = float(last)
        except Exception:
            out[f"final_{key}"] = str(last)
    try:
        out["epochs_recorded"] = int(max(len(v) for v in history.values()))
    except Exception:
        pass
    return out


def confidence_bin(prob: pd.Series) -> pd.Series:
    return pd.cut(
        prob,
        bins=[-np.inf, 0.5, 0.7, 0.85, 0.95, np.inf],
        labels=["lt_0.50", "0.50_0.70", "0.70_0.85", "0.85_0.95", "gte_0.95"],
    ).astype(str)


def main() -> int:
    args = parse_args()
    start = time.time()
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    args.processed_dir.mkdir(parents=True, exist_ok=True)

    model_dir = args.processed_dir / "model_scanvi_unified"
    output_h5ad = args.processed_dir / "scanvi_unified_labels.h5ad"
    seed_path = args.metadata_dir / "scanvi_unified_seed_labels_by_cell.tsv.gz"
    pred_path = args.metadata_dir / "scanvi_unified_predictions_by_cell.tsv.gz"
    prob_path = args.metadata_dir / "scanvi_unified_probability_matrix.tsv.gz"
    final_path = args.metadata_dir / "scanvi_unified_final_labels_by_cell.tsv.gz"
    strict_path = args.metadata_dir / "scanvi_unified_final_strict_labels_by_cell.tsv.gz"
    counts_path = args.metadata_dir / "scanvi_unified_label_counts.tsv"
    final_counts_path = args.metadata_dir / "scanvi_unified_final_label_counts.tsv"
    strict_counts_path = args.metadata_dir / "scanvi_unified_final_strict_label_counts.tsv"
    sample_path = args.metadata_dir / "scanvi_unified_by_sample.tsv"
    final_sample_path = args.metadata_dir / "scanvi_unified_final_by_sample.tsv"
    cluster_path = args.metadata_dir / "scanvi_unified_by_cluster.tsv"
    final_cluster_path = args.metadata_dir / "scanvi_unified_final_by_cluster.tsv"
    seed_pred_path = args.metadata_dir / "scanvi_unified_seed_prediction_crosstab.tsv"
    final_vs_malignant_path = args.metadata_dir / "scanvi_unified_final_vs_malignant_module3.tsv"
    strict_vs_malignant_path = args.metadata_dir / "scanvi_unified_final_strict_vs_malignant_module3.tsv"
    report_path = args.metadata_dir / "scanvi_unified_module4_report.json"

    scvi.settings.seed = args.seed
    if torch.cuda.is_available():
        torch.set_float32_matmul_precision("high")
    accelerator = "gpu" if torch.cuda.is_available() else "cpu"

    print(f"LOAD {args.input_h5ad}", flush=True)
    adata = ad.read_h5ad(args.input_h5ad)
    if "counts" not in adata.layers:
        raise ValueError("Input AnnData must contain layers['counts'] for scANVI.")
    for col in ["dataset", "sample_id", "study_sample"]:
        if col in adata.obs:
            adata.obs[col] = adata.obs[col].astype("category")

    seed = build_seed_table(args, pd.Index(adata.obs_names))
    if seed["cell_id"].duplicated().any():
        raise ValueError("Duplicated cell_id values in module 4 seed table.")
    seed = seed.set_index("cell_id").loc[adata.obs_names].rename_axis("cell_id").reset_index()
    labeled_count = int(seed["scanvi_unified_seed_label"].ne(UNLABELED).sum())
    print(f"SEEDS labeled={labeled_count} unlabeled={adata.n_obs - labeled_count}", flush=True)
    seed.to_csv(seed_path, sep="\t", index=False, compression="gzip")

    adata.obs["scanvi_unified_seed_label"] = pd.Categorical(seed["scanvi_unified_seed_label"].to_numpy())
    adata.obs["scanvi_unified_seed_source"] = pd.Categorical(seed["scanvi_unified_seed_source"].to_numpy())

    print(f"LOAD_SCVI_MODEL {args.scvi_model_dir}", flush=True)
    scvi_model = scvi.model.SCVI.load(
        str(args.scvi_model_dir),
        adata=adata,
        accelerator=accelerator,
        device=1 if accelerator == "gpu" else "auto",
    )
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        labels_key="scanvi_unified_seed_label",
        unlabeled_category=UNLABELED,
    )
    print(
        f"TRAIN_SCANVI accelerator={accelerator} max_epochs={args.max_epochs} "
        f"batch_size={args.batch_size} n_samples_per_label={args.n_samples_per_label}",
        flush=True,
    )
    scanvi_model.train(
        max_epochs=args.max_epochs,
        n_samples_per_label=args.n_samples_per_label,
        batch_size=args.batch_size,
        accelerator=accelerator,
        devices=1,
        early_stopping=True,
    )

    print("PREDICT", flush=True)
    pred = scanvi_model.predict(adata, batch_size=args.batch_size)
    probs = scanvi_model.predict(adata, soft=True, batch_size=args.batch_size)
    if not isinstance(probs, pd.DataFrame):
        raise TypeError("Expected scANVI soft predictions to return a DataFrame.")
    probs.index = adata.obs_names
    max_prob = probs.max(axis=1)
    second_prob = probs.mask(probs.eq(max_prob, axis=0)).max(axis=1).fillna(0.0)
    margin = max_prob - second_prob
    pred = pd.Series(np.asarray(pred).astype(str), index=adata.obs_names, name="scanvi_unified_pred_label")

    latent = scanvi_model.get_latent_representation(adata=adata, batch_size=args.batch_size)
    adata.obsm["X_scANVI"] = latent
    adata.obs["scanvi_unified_pred_label"] = pd.Categorical(pred.loc[adata.obs_names].to_numpy())
    adata.obs["scanvi_unified_pred_max_prob"] = max_prob.loc[adata.obs_names].to_numpy()
    adata.obs["scanvi_unified_pred_second_prob"] = second_prob.loc[adata.obs_names].to_numpy()
    adata.obs["scanvi_unified_pred_margin"] = margin.loc[adata.obs_names].to_numpy()
    adata.obs["scanvi_unified_pred_confidence_bin"] = pd.Categorical(confidence_bin(max_prob).loc[adata.obs_names].to_numpy())

    pred_df = pd.DataFrame(
        {
            "cell_id": adata.obs_names,
            "scanvi_unified_seed_label": adata.obs["scanvi_unified_seed_label"].astype(str).to_numpy(),
            "scanvi_unified_seed_source": adata.obs["scanvi_unified_seed_source"].astype(str).to_numpy(),
            "scanvi_unified_pred_label": pred.loc[adata.obs_names].to_numpy(),
            "scanvi_unified_pred_max_prob": max_prob.loc[adata.obs_names].to_numpy(),
            "scanvi_unified_pred_second_prob": second_prob.loc[adata.obs_names].to_numpy(),
            "scanvi_unified_pred_margin": margin.loc[adata.obs_names].to_numpy(),
            "scanvi_unified_pred_confidence_bin": confidence_bin(max_prob).loc[adata.obs_names].to_numpy(),
            "dataset": adata.obs["dataset"].astype(str).to_numpy() if "dataset" in adata.obs else "",
            "study_sample": adata.obs["study_sample"].astype(str).to_numpy() if "study_sample" in adata.obs else "",
            "sample_id": adata.obs["sample_id"].astype(str).to_numpy() if "sample_id" in adata.obs else "",
            "leiden_scvi": adata.obs["leiden_scvi"].astype(str).to_numpy() if "leiden_scvi" in adata.obs else "",
            "major_celltype": adata.obs["major_celltype"].astype(str).to_numpy() if "major_celltype" in adata.obs else "",
            "module3_malignant_hcc_call": seed["malignant_hcc_call"].fillna("not_module3_candidate").astype(str).to_numpy(),
        }
    )
    seeded_mask = pred_df["scanvi_unified_seed_label"].ne(UNLABELED)
    pred_df["scanvi_unified_final_label"] = pred_df["scanvi_unified_pred_label"]
    pred_df.loc[seeded_mask, "scanvi_unified_final_label"] = pred_df.loc[seeded_mask, "scanvi_unified_seed_label"]
    pred_df["scanvi_unified_final_source"] = np.where(
        seeded_mask,
        pred_df["scanvi_unified_seed_source"],
        "scanvi_propagated_from_unknown",
    )
    malignant_supported = pred_df["module3_malignant_hcc_call"].isin(
        ["malignant_hcc_high_conf", "malignant_hcc_cnv_support", "malignant_hcc_probable"]
    )
    unsupported_malignant = pred_df["scanvi_unified_final_label"].eq(MALIGNANT_LABEL) & ~malignant_supported
    pred_df["scanvi_unified_final_strict_label"] = pred_df["scanvi_unified_final_label"]
    pred_df.loc[unsupported_malignant, "scanvi_unified_final_strict_label"] = "malignant_like_hepatocyte_needs_review"
    pred_df["scanvi_unified_final_strict_source"] = pred_df["scanvi_unified_final_source"]
    pred_df.loc[
        pred_df["scanvi_unified_final_label"].eq(MALIGNANT_LABEL) & malignant_supported,
        "scanvi_unified_final_strict_source",
    ] = "module3_supported_malignant"
    pred_df.loc[unsupported_malignant, "scanvi_unified_final_strict_source"] = "scanvi_malignant_like_needs_module3_review"
    adata.obs["scanvi_unified_final_label"] = pd.Categorical(pred_df["scanvi_unified_final_label"].to_numpy())
    adata.obs["scanvi_unified_final_source"] = pd.Categorical(pred_df["scanvi_unified_final_source"].to_numpy())
    adata.obs["scanvi_unified_final_strict_label"] = pd.Categorical(pred_df["scanvi_unified_final_strict_label"].to_numpy())
    adata.obs["scanvi_unified_final_strict_source"] = pd.Categorical(pred_df["scanvi_unified_final_strict_source"].to_numpy())
    pred_df.to_csv(pred_path, sep="\t", index=False, compression="gzip")
    pred_df.to_csv(final_path, sep="\t", index=False, compression="gzip")
    pred_df.to_csv(strict_path, sep="\t", index=False, compression="gzip")
    probs.insert(0, "cell_id", probs.index.astype(str))
    probs.to_csv(prob_path, sep="\t", index=False, compression="gzip")

    counts = (
        pred_df.groupby(["scanvi_unified_pred_label", "scanvi_unified_pred_confidence_bin"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["scanvi_unified_pred_label", "scanvi_unified_pred_confidence_bin"])
    )
    counts.to_csv(counts_path, sep="\t", index=False)
    final_counts = (
        pred_df.groupby(["scanvi_unified_final_label", "scanvi_unified_final_source"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["scanvi_unified_final_label", "scanvi_unified_final_source"])
    )
    final_counts.to_csv(final_counts_path, sep="\t", index=False)
    strict_counts = (
        pred_df.groupby(["scanvi_unified_final_strict_label", "scanvi_unified_final_strict_source"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["scanvi_unified_final_strict_label", "scanvi_unified_final_strict_source"])
    )
    strict_counts.to_csv(strict_counts_path, sep="\t", index=False)

    by_sample = (
        pred_df.groupby(["dataset", "study_sample", "sample_id", "scanvi_unified_pred_label"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["dataset", "study_sample", "scanvi_unified_pred_label"])
    )
    by_sample.to_csv(sample_path, sep="\t", index=False)
    final_by_sample = (
        pred_df.groupby(["dataset", "study_sample", "sample_id", "scanvi_unified_final_label"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["dataset", "study_sample", "scanvi_unified_final_label"])
    )
    final_by_sample.to_csv(final_sample_path, sep="\t", index=False)

    by_cluster = (
        pred_df.groupby(["leiden_scvi", "scanvi_unified_pred_label"], observed=True)
        .agg(n_cells=("cell_id", "size"), mean_max_prob=("scanvi_unified_pred_max_prob", "mean"))
        .reset_index()
        .sort_values(["leiden_scvi", "n_cells"], ascending=[True, False])
    )
    by_cluster.to_csv(cluster_path, sep="\t", index=False)
    final_by_cluster = (
        pred_df.groupby(["leiden_scvi", "scanvi_unified_final_label"], observed=True)
        .agg(n_cells=("cell_id", "size"), mean_max_prob=("scanvi_unified_pred_max_prob", "mean"))
        .reset_index()
        .sort_values(["leiden_scvi", "n_cells"], ascending=[True, False])
    )
    final_by_cluster.to_csv(final_cluster_path, sep="\t", index=False)

    labeled_pred = pred_df.loc[pred_df["scanvi_unified_seed_label"].ne(UNLABELED)].copy()
    seed_pred = pd.crosstab(labeled_pred["scanvi_unified_seed_label"], labeled_pred["scanvi_unified_pred_label"])
    seed_pred.to_csv(seed_pred_path, sep="\t")
    final_vs_malignant = pd.crosstab(
        pred_df["module3_malignant_hcc_call"],
        pred_df["scanvi_unified_final_label"],
    )
    strict_vs_malignant = pd.crosstab(
        pred_df["module3_malignant_hcc_call"],
        pred_df["scanvi_unified_final_strict_label"],
    )
    final_vs_malignant.to_csv(final_vs_malignant_path, sep="\t")
    strict_vs_malignant.to_csv(strict_vs_malignant_path, sep="\t")

    scanvi_model.save(model_dir, overwrite=True)
    adata.write_h5ad(output_h5ad, compression="gzip")

    report = {
        "method": "scANVI semi-supervised unified label propagation",
        "input_h5ad": str(args.input_h5ad.resolve()),
        "scvi_model_dir": str(args.scvi_model_dir.resolve()),
        "model_dir": str(model_dir.resolve()),
        "output_h5ad": str(output_h5ad.resolve()),
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "unlabeled_category": UNLABELED,
        "seed_label_counts": seed["scanvi_unified_seed_label"].value_counts(dropna=False).to_dict(),
        "seed_source_counts": seed["scanvi_unified_seed_source"].value_counts(dropna=False).to_dict(),
        "pred_label_counts": pred_df["scanvi_unified_pred_label"].value_counts(dropna=False).to_dict(),
        "final_label_counts": pred_df["scanvi_unified_final_label"].value_counts(dropna=False).to_dict(),
        "final_strict_label_counts": pred_df["scanvi_unified_final_strict_label"].value_counts(dropna=False).to_dict(),
        "final_strict_label_rule": "Preserve seeds and propagated labels, but keep malignant_hepatocyte only when module 3 call is high_conf, cnv_support, or probable.",
        "mean_prediction_max_prob": float(pred_df["scanvi_unified_pred_max_prob"].mean()),
        "median_prediction_max_prob": float(pred_df["scanvi_unified_pred_max_prob"].median()),
        "max_epochs": int(args.max_epochs),
        "batch_size": int(args.batch_size),
        "n_samples_per_label": int(args.n_samples_per_label),
        "accelerator": accelerator,
        "torch_version": torch.__version__,
        "scvi_version": scvi.__version__,
        "outputs": {
            "seed_labels_by_cell": str(seed_path.resolve()),
            "predictions_by_cell": str(pred_path.resolve()),
            "probability_matrix": str(prob_path.resolve()),
            "final_labels_by_cell": str(final_path.resolve()),
            "final_strict_labels_by_cell": str(strict_path.resolve()),
            "label_counts": str(counts_path.resolve()),
            "final_label_counts": str(final_counts_path.resolve()),
            "final_strict_label_counts": str(strict_counts_path.resolve()),
            "by_sample": str(sample_path.resolve()),
            "final_by_sample": str(final_sample_path.resolve()),
            "by_cluster": str(cluster_path.resolve()),
            "final_by_cluster": str(final_cluster_path.resolve()),
            "seed_prediction_crosstab": str(seed_pred_path.resolve()),
            "final_vs_malignant_module3": str(final_vs_malignant_path.resolve()),
            "final_strict_vs_malignant_module3": str(strict_vs_malignant_path.resolve()),
        },
        "elapsed_seconds": round(time.time() - start, 3),
    }
    report.update(history_tail(scanvi_model))
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
