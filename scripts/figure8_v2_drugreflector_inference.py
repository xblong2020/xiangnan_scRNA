from __future__ import annotations

"""Frozen DrugReflector inference adapter for Figure 8 v2.

This module performs model inference only. Signature construction, statistics,
integration, and every formal figure remain in R.
"""

import argparse
import hashlib
import importlib.util
import json
import platform
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import scipy.stats as stats


ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "metadata" / "driver" / "figure8_transcriptomic_reversal_v2_mainfigure"
DATA = ROOT / "data" / "processed" / "driver" / "figure8_transcriptomic_reversal_v2_mainfigure"
CHECKPOINTS = ROOT / "metadata" / "driver" / "drugreflector_checkpoints"
SOURCE = ROOT / "tmp" / "drugreflector-main-from-zip" / "drugreflector-main"
MODEL_GENES = ROOT / "metadata" / "driver" / "figure8_transcriptomic_reversal" / "figure8_drugreflector_model_genes.tsv"

EXPECTED_CHECKPOINT_MD5 = [
    "0a27e253713c37f4874318b5ba0c27a9",
    "0e785196fd046d946f84e4480c81ff53",
    "d8e36f6a8f9fa7a22feda7acdd0bee86",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Figure 8 v2 frozen DrugReflector inference only")
    parser.add_argument("--mode", choices=["variants", "random"], required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--watchlist", type=Path)
    parser.add_argument("--metadata-dir", type=Path, default=META)
    parser.add_argument("--checkpoint-dir", type=Path, default=CHECKPOINTS)
    parser.add_argument("--source-dir", type=Path, default=SOURCE)
    parser.add_argument("--top-n", type=int, default=200)
    parser.add_argument("--batch-size", type=int, default=25)
    parser.add_argument("--seed", type=int, default=20260805)
    return parser.parse_args()


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_checkpoints(checkpoint_dir: Path) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for fold, expected in enumerate(EXPECTED_CHECKPOINT_MD5):
        path = checkpoint_dir / f"model_fold_{fold}.pt"
        if not path.is_file():
            raise FileNotFoundError(f"Missing frozen checkpoint: {path}")
        observed = md5(path)
        if observed != expected:
            raise ValueError(f"Checkpoint MD5 mismatch for fold {fold}: {observed} != {expected}")
        records.append({"fold": fold, "path": str(path.resolve()), "size_bytes": path.stat().st_size, "md5": observed})
    return records


def validate_model_order(frame: pd.DataFrame, model_genes: list[str]) -> None:
    observed = [str(column).strip().upper() for column in frame.columns]
    expected = [str(gene).strip().upper() for gene in model_genes]
    if len(expected) != 978 or observed != expected:
        raise ValueError("Input columns must match the exact frozen 978-gene order")


def rank_descending(scores: np.ndarray) -> np.ndarray:
    return stats.rankdata(-np.asarray(scores), axis=1, method="min").astype(np.int32)


def load_v1_adapter():
    path = ROOT / "scripts" / "figure8_drugreflector_inference.py"
    spec = importlib.util.spec_from_file_location("figure8_v1_inference", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import frozen v1 inference adapter: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_model(args: argparse.Namespace):
    adapter = load_v1_adapter()
    model_args = SimpleNamespace(source_dir=args.source_dir, checkpoint_dir=args.checkpoint_dir)
    model, checkpoint_paths = adapter.load_model(model_args)
    checkpoint_records = validate_checkpoints(args.checkpoint_dir)
    genes = [str(x).upper() for x in model.model.dimensions["var_names"][0]]
    if len(genes) != 978 or any([str(x).upper() for x in fold] != genes for fold in model.model.dimensions["var_names"]):
        raise ValueError("DrugReflector folds do not share the frozen 978-gene input order")
    if [str(p.resolve()) for p in checkpoint_paths] != [row["path"] for row in checkpoint_records]:
        raise ValueError("Loaded checkpoint paths differ from the validated paths")
    return adapter, model, genes, checkpoint_records


def read_input(path: Path, model_genes: list[str]) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", compression="infer")
    if "signature_id" not in frame.columns:
        raise ValueError("Input matrix must contain signature_id")
    frame = frame.set_index("signature_id")
    if frame.index.duplicated().any():
        raise ValueError("signature_id values must be unique")
    validate_model_order(frame, model_genes)
    frame.columns = [str(column).strip().upper() for column in frame.columns]
    return frame.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def infer(adapter, model, frame: pd.DataFrame):
    scores, _, probabilities, fold_scores, _ = adapter.infer_batch(model, frame)
    ensemble_ranks = rank_descending(scores)
    fold_ranks = np.stack([rank_descending(fold_scores[fold]) for fold in range(fold_scores.shape[0])], axis=0)
    return scores, ensemble_ranks, probabilities, fold_scores, fold_ranks


def run_variants(adapter, model, frame: pd.DataFrame, metadata_dir: Path) -> dict[str, object]:
    scores, ranks, probabilities, fold_scores, fold_ranks = infer(adapter, model, frame)
    compounds = np.asarray(model.compound_names, dtype=str)
    prediction_rows = []
    fold_rows = []
    for row_idx, signature_id in enumerate(frame.index.astype(str)):
        prediction_rows.append(
            pd.DataFrame(
                {
                    "signature_id": signature_id,
                    "compound": compounds,
                    "rank_1based": ranks[row_idx],
                    "logit": scores[row_idx],
                    "probability": probabilities[row_idx],
                }
            )
        )
        for fold in range(3):
            fold_rows.append(
                pd.DataFrame(
                    {
                        "signature_id": signature_id,
                        "fold": fold,
                        "compound": compounds,
                        "fold_rank_1based": fold_ranks[fold, row_idx],
                        "fold_logit": fold_scores[fold, row_idx],
                    }
                )
            )
    prediction_path = metadata_dir / "figure8_v2_drugreflector_variant_predictions.tsv.gz"
    fold_path = metadata_dir / "figure8_v2_drugreflector_fold_predictions.tsv.gz"
    pd.concat(prediction_rows, ignore_index=True).to_csv(prediction_path, sep="\t", index=False, compression="gzip")
    pd.concat(fold_rows, ignore_index=True).to_csv(fold_path, sep="\t", index=False, compression="gzip")
    return {
        "variant_predictions": str(prediction_path.resolve()),
        "fold_predictions": str(fold_path.resolve()),
        "n_signatures": int(frame.shape[0]),
        "n_compounds": int(len(compounds)),
    }


def read_watchlist(path: Path | None, compounds: np.ndarray) -> list[str]:
    if path is None or not path.is_file():
        return []
    frame = pd.read_csv(path, sep="\t")
    if "compound" not in frame.columns:
        raise ValueError("Watchlist must contain compound")
    universe = set(compounds)
    return [x for x in frame["compound"].dropna().astype(str).drop_duplicates() if x in universe]


def run_random(
    adapter,
    model,
    frame: pd.DataFrame,
    metadata_dir: Path,
    watchlist_path: Path | None,
    top_n: int,
    batch_size: int,
) -> dict[str, object]:
    compounds = np.asarray(model.compound_names, dtype=str)
    compound_index = {name: idx for idx, name in enumerate(compounds)}
    watchlist = read_watchlist(watchlist_path, compounds)
    summaries: list[dict[str, object]] = []
    top_rows: list[pd.DataFrame] = []
    watch_rows: list[pd.DataFrame] = []
    for start in range(0, frame.shape[0], batch_size):
        batch = frame.iloc[start : start + batch_size]
        scores, ranks, probabilities, _, fold_ranks = infer(adapter, model, batch)
        for row_idx, signature_id in enumerate(batch.index.astype(str)):
            order = np.argsort(ranks[row_idx], kind="stable")[:top_n]
            top_idx = order[0]
            fold_top = fold_ranks[:, row_idx, top_idx]
            summaries.append(
                {
                    "signature_id": signature_id,
                    "top_compound": compounds[top_idx],
                    "max_probability": float(probabilities[row_idx, top_idx]),
                    "max_logit": float(scores[row_idx, top_idx]),
                    "top10_probability_concentration": float(probabilities[row_idx, order[:10]].sum()),
                    "top100_probability_concentration": float(probabilities[row_idx, order[:100]].sum()),
                    "top200_probability_concentration": float(probabilities[row_idx, order[:200]].sum()),
                    "top_candidate_fold_rank_min": int(fold_top.min()),
                    "top_candidate_fold_rank_max": int(fold_top.max()),
                    "top_candidate_model_agreement": float(max(0, 1 - (fold_top.max() - fold_top.min()) / max(1, len(compounds) - 1))),
                    "n_nonzero_input_genes": int(np.count_nonzero(batch.iloc[row_idx].to_numpy())),
                }
            )
            top_rows.append(
                pd.DataFrame(
                    {
                        "signature_id": signature_id,
                        "compound": compounds[order],
                        "rank_1based": ranks[row_idx, order],
                        "probability": probabilities[row_idx, order],
                    }
                )
            )
            if watchlist:
                indices = np.asarray([compound_index[x] for x in watchlist], dtype=int)
                watch_rows.append(
                    pd.DataFrame(
                        {
                            "signature_id": signature_id,
                            "compound": np.asarray(watchlist),
                            "rank_1based": ranks[row_idx, indices],
                            "probability": probabilities[row_idx, indices],
                        }
                    )
                )
    summary_path = metadata_dir / "figure8_v2_matched_random_inference_summary.tsv.gz"
    top_path = metadata_dir / "figure8_v2_matched_random_top_predictions.tsv.gz"
    watch_path = metadata_dir / "figure8_v2_matched_random_watchlist_predictions.tsv.gz"
    pd.DataFrame(summaries).to_csv(summary_path, sep="\t", index=False, compression="gzip")
    pd.concat(top_rows, ignore_index=True).to_csv(top_path, sep="\t", index=False, compression="gzip")
    (pd.concat(watch_rows, ignore_index=True) if watch_rows else pd.DataFrame(columns=["signature_id", "compound", "rank_1based", "probability"])).to_csv(
        watch_path, sep="\t", index=False, compression="gzip"
    )
    return {
        "random_summary": str(summary_path.resolve()),
        "random_top_predictions": str(top_path.resolve()),
        "random_watchlist_predictions": str(watch_path.resolve()),
        "n_signatures": int(frame.shape[0]),
        "n_compounds": int(len(compounds)),
        "n_watchlist": len(watchlist),
    }


def main() -> None:
    args = parse_args()
    np.random.seed(args.seed)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    adapter, model, model_genes, checkpoints = load_model(args)
    frame = read_input(args.input, model_genes)
    if args.mode == "variants":
        outputs = run_variants(adapter, model, frame, args.metadata_dir)
    else:
        outputs = run_random(adapter, model, frame, args.metadata_dir, args.watchlist, args.top_n, args.batch_size)
    report = {
        "module": "figure8_v2_drugreflector_inference",
        "status": "completed",
        "mode": args.mode,
        "seed": args.seed,
        "model_version": "DrugReflector V3.5 frozen three-fold ensemble",
        "input": str(args.input.resolve()),
        "input_md5": md5(args.input),
        "checkpoint_files": checkpoints,
        "model_gene_count": len(model_genes),
        "compound_count": int(model.n_compounds),
        "source_dir": str(args.source_dir.resolve()),
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__},
        "outputs": outputs,
        "interpretation_boundary": "Computational prioritization only; no efficacy, normal-cell safety, or phenotypic rescue claim.",
    }
    report_path = args.metadata_dir / f"figure8_v2_drugreflector_{args.mode}_inference_report.json"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(outputs, ensure_ascii=False))


if __name__ == "__main__":
    main()

