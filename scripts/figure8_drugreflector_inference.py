from __future__ import annotations

"""Figure 8-only DrugReflector checkpoint inference.

This script performs model loading/inference only. Signature construction,
statistics, integration, and plotting remain in R.
"""

import argparse
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.special import softmax


ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "metadata" / "driver" / "figure8_transcriptomic_reversal"
CHECKPOINTS = ROOT / "metadata" / "driver" / "drugreflector_checkpoints"
SOURCE = ROOT / "tmp" / "drugreflector-main-from-zip" / "drugreflector-main"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Figure 8 DrugReflector inference only")
    p.add_argument("--mode", choices=["export", "variants", "toxicity", "random"], required=True)
    p.add_argument("--input", type=Path)
    p.add_argument("--watchlist", type=Path)
    p.add_argument("--metadata-dir", type=Path, default=META)
    p.add_argument("--checkpoint-dir", type=Path, default=CHECKPOINTS)
    p.add_argument("--source-dir", type=Path, default=SOURCE)
    p.add_argument("--top-n", type=int, default=200)
    p.add_argument("--batch-size", type=int, default=25)
    p.add_argument("--seed", type=int, default=20260805)
    return p.parse_args()


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_model(args: argparse.Namespace):
    if args.source_dir.exists():
        sys.path.insert(0, str(args.source_dir.resolve()))
    from drugreflector import DrugReflector

    paths = [args.checkpoint_dir / f"model_fold_{idx}.pt" for idx in range(3)]
    missing = [str(p) for p in paths if not p.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing checkpoints: {missing}")
    return DrugReflector(checkpoint_paths=[str(p) for p in paths]), paths


def read_wide(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", compression="infer")
    if "signature_id" not in frame.columns:
        raise ValueError("signature wide input must contain signature_id")
    frame = frame.set_index("signature_id")
    frame.columns = frame.columns.astype(str).str.strip().str.upper()
    frame = frame.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    if frame.index.duplicated().any():
        raise ValueError("signature_id values must be unique")
    return frame


def prepare_for_folds(model, frame: pd.DataFrame):
    from drugreflector.utils import clip_rescale_rows

    adata = model._prepare_vscores(frame)
    clip_rescale_rows(adata.X, clip=2, target_std=1)
    adata.var_names = adata.var_names.str.upper()
    adata.var_names = adata.var_names.str.replace(r"^ENSG\d+\.", "", regex=True)
    adata.var_names = adata.var_names.str.replace(r"_AT$", "", regex=True)
    adata.var_names = adata.var_names.str.replace(r"\..*$", "", regex=True)
    adata.var_names = adata.var_names.str.replace(r"[^A-Z0-9\-]", "", regex=True)
    adata.var_names_make_unique()
    formatted = model.model.format_vscores(adata)
    fold_scores = model.model.get_predictions(formatted, average=False)
    return fold_scores


def infer_batch(model, frame: pd.DataFrame):
    fold_scores = prepare_for_folds(model, frame)
    ensemble_scores = fold_scores.mean(axis=0)
    ensemble_ranks = stats.rankdata(-ensemble_scores, axis=1, method="average").astype(np.int32)
    fold_ranks = stats.rankdata(-fold_scores, axis=2, method="average").astype(np.int32)
    probabilities = softmax(ensemble_scores, axis=1)
    return ensemble_scores, ensemble_ranks, probabilities, fold_scores, fold_ranks


def long_predictions(names, compounds, scores, ranks, probs) -> pd.DataFrame:
    rows = []
    for idx, signature_id in enumerate(names):
        rows.append(
            pd.DataFrame(
                {
                    "signature_id": signature_id,
                    "compound": compounds,
                    "rank_1based": ranks[idx],
                    "logit": scores[idx],
                    "prob": probs[idx],
                }
            )
        )
    return pd.concat(rows, ignore_index=True)


def run_export(model, args: argparse.Namespace, checkpoints: list[Path]) -> dict:
    genes = list(map(str, model.model.dimensions["var_names"][0]))
    if any(list(map(str, x)) != genes for x in model.model.dimensions["var_names"][1:]):
        raise ValueError("DrugReflector folds have different input gene order")
    out = args.metadata_dir / "figure8_drugreflector_model_genes.tsv"
    pd.DataFrame({"model_gene_order": np.arange(1, len(genes) + 1), "gene": genes}).to_csv(out, sep="\t", index=False)
    return {"model_genes": str(out.resolve()), "n_model_genes": len(genes), "n_compounds": int(model.n_compounds)}


def run_variants(model, args: argparse.Namespace) -> dict:
    if args.input is None:
        raise ValueError("--input is required for variants")
    wide = read_wide(args.input)
    scores, ranks, probs, _, fold_ranks = infer_batch(model, wide)
    compounds = np.asarray(model.compound_names, dtype=str)
    pred = long_predictions(wide.index.astype(str), compounds, scores, ranks, probs)
    pred_path = args.metadata_dir / "figure8_drugreflector_variant_predictions.tsv.gz"
    pred.to_csv(pred_path, sep="\t", index=False, compression="gzip")

    fold_rows = []
    for fold in range(fold_ranks.shape[0]):
        for idx, signature_id in enumerate(wide.index.astype(str)):
            fold_rows.append(
                pd.DataFrame(
                    {
                        "signature_id": signature_id,
                        "fold": fold,
                        "compound": compounds,
                        "fold_rank_1based": fold_ranks[fold, idx],
                    }
                )
            )
    fold = pd.concat(fold_rows, ignore_index=True)
    fold_path = args.metadata_dir / "figure8_drugreflector_fold_predictions.tsv.gz"
    fold.to_csv(fold_path, sep="\t", index=False, compression="gzip")
    return {
        "variant_predictions": str(pred_path.resolve()),
        "fold_predictions": str(fold_path.resolve()),
        "n_signatures": int(wide.shape[0]),
        "n_compounds": int(len(compounds)),
    }


def run_toxicity(model, args: argparse.Namespace) -> dict:
    if args.input is None:
        raise ValueError("--input is required for toxicity")
    wide = read_wide(args.input)
    scores, ranks, probs, _, _ = infer_batch(model, wide)
    compounds = np.asarray(model.compound_names, dtype=str)
    pred = long_predictions(wide.index.astype(str), compounds, scores, ranks, probs)
    pred_path = args.metadata_dir / "figure8_drugreflector_toxicity_control_predictions.tsv.gz"
    pred.to_csv(pred_path, sep="\t", index=False, compression="gzip")
    return {
        "toxicity_control_predictions": str(pred_path.resolve()),
        "n_signatures": int(wide.shape[0]),
        "n_compounds": int(len(compounds)),
    }


def run_random(model, args: argparse.Namespace) -> dict:
    if args.input is None:
        raise ValueError("--input is required for random")
    wide = read_wide(args.input)
    watch = []
    if args.watchlist and args.watchlist.is_file():
        watch_frame = pd.read_csv(args.watchlist, sep="\t")
        if "compound" in watch_frame.columns:
            watch = watch_frame["compound"].dropna().astype(str).drop_duplicates().tolist()
    compounds = np.asarray(model.compound_names, dtype=str)
    compound_index = {name: idx for idx, name in enumerate(compounds)}
    watch = [name for name in watch if name in compound_index]

    summaries, tops, watched = [], [], []
    for start in range(0, wide.shape[0], args.batch_size):
        batch = wide.iloc[start : start + args.batch_size]
        scores, ranks, probs, _, _ = infer_batch(model, batch)
        for row_idx, signature_id in enumerate(batch.index.astype(str)):
            order = np.argsort(ranks[row_idx])[: args.top_n]
            summaries.append(
                {
                    "signature_id": signature_id,
                    "top_compound": compounds[order[0]],
                    "max_probability": float(probs[row_idx, order[0]]),
                    "max_logit": float(scores[row_idx, order[0]]),
                    "n_nonzero_input_genes": int(np.count_nonzero(batch.iloc[row_idx].to_numpy())),
                }
            )
            tops.append(
                pd.DataFrame(
                    {
                        "signature_id": signature_id,
                        "compound": compounds[order],
                        "rank_1based": ranks[row_idx, order],
                        "prob": probs[row_idx, order],
                    }
                )
            )
            if watch:
                idxs = np.asarray([compound_index[name] for name in watch])
                watched.append(
                    pd.DataFrame(
                        {
                            "signature_id": signature_id,
                            "compound": np.asarray(watch),
                            "rank_1based": ranks[row_idx, idxs],
                            "prob": probs[row_idx, idxs],
                        }
                    )
                )
    summary_path = args.metadata_dir / "figure8_random_signature_inference_summary.tsv.gz"
    top_path = args.metadata_dir / "figure8_random_signature_top_predictions.tsv.gz"
    watch_path = args.metadata_dir / "figure8_random_signature_watchlist_predictions.tsv.gz"
    pd.DataFrame(summaries).to_csv(summary_path, sep="\t", index=False, compression="gzip")
    pd.concat(tops, ignore_index=True).to_csv(top_path, sep="\t", index=False, compression="gzip")
    (pd.concat(watched, ignore_index=True) if watched else pd.DataFrame(columns=["signature_id", "compound", "rank_1based", "prob"])).to_csv(
        watch_path, sep="\t", index=False, compression="gzip"
    )
    return {
        "random_summary": str(summary_path.resolve()),
        "random_top_predictions": str(top_path.resolve()),
        "random_watchlist_predictions": str(watch_path.resolve()),
        "n_signatures": int(wide.shape[0]),
        "n_watchlist_compounds": len(watch),
        "n_compounds": int(len(compounds)),
    }


def main() -> None:
    args = parse_args()
    np.random.seed(args.seed)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    model, checkpoints = load_model(args)
    if args.mode == "export":
        outputs = run_export(model, args, checkpoints)
    elif args.mode == "variants":
        outputs = run_variants(model, args)
    elif args.mode == "toxicity":
        outputs = run_toxicity(model, args)
    else:
        outputs = run_random(model, args)
    report = {
        "module": "figure8_drugreflector_inference",
        "mode": args.mode,
        "status": "completed",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "seed": args.seed,
        "input": str(args.input.resolve()) if args.input else None,
        "checkpoint_files": [
            {"path": str(p.resolve()), "size_bytes": p.stat().st_size, "md5": md5(p)} for p in checkpoints
        ],
        "source_dir": str(args.source_dir.resolve()),
        "outputs": outputs,
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "interpretation_boundary": "Model scores are computational prioritization outputs and do not establish efficacy, safety, or direct expression rescue.",
    }
    report_path = args.metadata_dir / f"figure8_drugreflector_{args.mode}_inference_report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(outputs, sort_keys=True))


if __name__ == "__main__":
    main()
