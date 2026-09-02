from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA_DIR = ROOT / "metadata/driver"
DEFAULT_PRIMARY_VSCORE = DEFAULT_METADATA_DIR / "module9_6_drugreflector_vscore_primary.tsv"
DEFAULT_SENSITIVITY_VSCORE = DEFAULT_METADATA_DIR / "module9_6_drugreflector_vscore_sensitivity.tsv"
DEFAULT_CHECKPOINT_DIR = DEFAULT_METADATA_DIR / "drugreflector_checkpoints"
DEFAULT_DRUGREFLECTOR_SOURCE = ROOT / "tmp/drugreflector-main-from-zip/drugreflector-main"
OUTPUT_STEM = "module9_7_drugreflector"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 9.7 DrugReflector inference.")
    parser.add_argument("--primary-vscore", type=Path, default=DEFAULT_PRIMARY_VSCORE)
    parser.add_argument("--sensitivity-vscore", type=Path, default=DEFAULT_SENSITIVITY_VSCORE)
    parser.add_argument("--checkpoint-dir", type=Path, default=DEFAULT_CHECKPOINT_DIR)
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA_DIR)
    parser.add_argument("--drugreflector-source-dir", type=Path, default=DEFAULT_DRUGREFLECTOR_SOURCE)
    parser.add_argument("--top-n", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260616)
    parser.add_argument("--skip-sensitivity", action="store_true")
    return parser.parse_args()


def package_version(name: str) -> str:
    try:
        return version(name)
    except PackageNotFoundError:
        return "not_installed"


def configure_drugreflector_source(source_dir: Path) -> None:
    if source_dir.exists():
        sys.path.insert(0, str(source_dir.resolve()))


def checkpoint_paths(checkpoint_dir: Path) -> list[Path]:
    paths = [checkpoint_dir / f"model_fold_{idx}.pt" for idx in range(3)]
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing DrugReflector checkpoints: {missing}")
    return paths


def read_vscore_series(path: Path, label: str) -> pd.Series:
    frame = pd.read_csv(path, sep="\t")
    required = {"gene", "v_score"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"{path} missing required columns: {sorted(missing)}")
    frame = frame.loc[:, ["gene", "v_score"]].copy()
    frame["gene"] = frame["gene"].astype(str).str.strip().str.upper()
    frame["v_score"] = pd.to_numeric(frame["v_score"], errors="coerce")
    frame = frame.loc[frame["gene"].ne("") & frame["v_score"].notna()].copy()
    if frame["gene"].duplicated().any():
        frame = frame.groupby("gene", as_index=False)["v_score"].mean()
    series = pd.Series(frame["v_score"].to_numpy(dtype=float), index=frame["gene"], name=label)
    if series.empty:
        raise ValueError(f"{path} produced an empty v-score series")
    return series


def model_gene_sets(model: object) -> list[set[str]]:
    var_names = model.model.dimensions["var_names"]
    return [set(str(gene).upper() for gene in fold_names) for fold_names in var_names]


def build_gene_coverage_rows(label: str, series: pd.Series, fold_gene_sets: Sequence[set[str]]) -> list[dict[str, object]]:
    query_genes = set(str(gene).upper() for gene in series.index)
    rows: list[dict[str, object]] = []
    for fold_idx, genes in enumerate(fold_gene_sets):
        overlap = sorted(query_genes.intersection(genes))
        rows.append(
            {
                "signature": label,
                "fold": fold_idx,
                "n_query_genes": len(query_genes),
                "n_model_landmark_genes": len(genes),
                "n_overlap_genes": len(overlap),
                "query_gene_coverage_fraction": len(overlap) / len(query_genes) if query_genes else 0.0,
                "model_gene_coverage_fraction": len(overlap) / len(genes) if genes else 0.0,
                "overlap_genes": ",".join(overlap),
            }
        )
    union_genes = set().union(*fold_gene_sets) if fold_gene_sets else set()
    union_overlap = sorted(query_genes.intersection(union_genes))
    rows.append(
        {
            "signature": label,
            "fold": "union",
            "n_query_genes": len(query_genes),
            "n_model_landmark_genes": len(union_genes),
            "n_overlap_genes": len(union_overlap),
            "query_gene_coverage_fraction": len(union_overlap) / len(query_genes) if query_genes else 0.0,
            "model_gene_coverage_fraction": len(union_overlap) / len(union_genes) if union_genes else 0.0,
            "overlap_genes": ",".join(union_overlap),
        }
    )
    return rows


def flatten_top_compounds(results: dict[str, pd.DataFrame], source_label: str) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for signature, frame in results.items():
        out = frame.copy()
        out.insert(0, "signature", signature)
        out.insert(1, "source_label", source_label)
        out["rank_0based"] = pd.to_numeric(out["rank"], errors="coerce").astype(int)
        out["rank_1based"] = out["rank_0based"] + 1
        out = out.drop(columns=["rank"])
        out = out[["signature", "source_label", "compound", "rank_0based", "rank_1based", "logit", "prob"]]
        frames.append(out)
    if not frames:
        return pd.DataFrame(columns=["signature", "source_label", "compound", "rank_0based", "rank_1based", "logit", "prob"])
    return pd.concat(frames, ignore_index=True).sort_values(["signature", "rank_0based", "compound"]).reset_index(drop=True)


def build_consensus_table(primary: pd.DataFrame, sensitivity: pd.DataFrame | None = None) -> pd.DataFrame:
    primary_cols = primary.rename(
        columns={
            "rank_0based": "primary_rank_0based",
            "rank_1based": "primary_rank_1based",
            "logit": "primary_logit",
            "prob": "primary_prob",
        }
    )[["compound", "primary_rank_0based", "primary_rank_1based", "primary_logit", "primary_prob"]]
    if sensitivity is None or sensitivity.empty:
        out = primary_cols.copy()
        out["sensitivity_rank_0based"] = np.nan
        out["sensitivity_rank_1based"] = np.nan
        out["sensitivity_logit"] = np.nan
        out["sensitivity_prob"] = np.nan
    else:
        sensitivity_cols = sensitivity.rename(
            columns={
                "rank_0based": "sensitivity_rank_0based",
                "rank_1based": "sensitivity_rank_1based",
                "logit": "sensitivity_logit",
                "prob": "sensitivity_prob",
            }
        )[["compound", "sensitivity_rank_0based", "sensitivity_rank_1based", "sensitivity_logit", "sensitivity_prob"]]
        out = primary_cols.merge(sensitivity_cols, on="compound", how="outer")

    out["in_primary_top"] = out["primary_rank_0based"].notna()
    out["in_sensitivity_top"] = out["sensitivity_rank_0based"].notna()
    out["in_both_top_lists"] = out["in_primary_top"] & out["in_sensitivity_top"]
    rank_cols = ["primary_rank_1based", "sensitivity_rank_1based"]
    out["mean_rank_1based"] = out[rank_cols].mean(axis=1, skipna=True)
    out["best_rank_1based"] = out[rank_cols].min(axis=1, skipna=True)
    out["rank_delta_abs"] = (out["primary_rank_1based"] - out["sensitivity_rank_1based"]).abs()
    out = out.sort_values(
        ["in_both_top_lists", "mean_rank_1based", "best_rank_1based", "compound"],
        ascending=[False, True, True, True],
    ).reset_index(drop=True)
    return out


def run_prediction(model: object, series: pd.Series, label: str, top_n: int) -> pd.DataFrame:
    top = model.get_top_compounds(series, n_top=top_n)
    return flatten_top_compounds(top, label)


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def main() -> None:
    args = parse_args()
    np.random.seed(args.seed)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    configure_drugreflector_source(args.drugreflector_source_dir)

    from drugreflector import DrugReflector

    ckpts = checkpoint_paths(args.checkpoint_dir)
    model = DrugReflector(checkpoint_paths=[str(path) for path in ckpts])
    fold_gene_sets = model_gene_sets(model)

    primary_series = read_vscore_series(args.primary_vscore, "module9_7_primary")
    prediction_frames = {
        "primary": run_prediction(model, primary_series, "primary", args.top_n),
    }
    coverage_rows = build_gene_coverage_rows("primary", primary_series, fold_gene_sets)

    if not args.skip_sensitivity:
        sensitivity_series = read_vscore_series(args.sensitivity_vscore, "module9_7_sensitivity")
        prediction_frames["sensitivity"] = run_prediction(model, sensitivity_series, "sensitivity", args.top_n)
        coverage_rows.extend(build_gene_coverage_rows("sensitivity", sensitivity_series, fold_gene_sets))

    outputs = {
        "primary_predictions": args.metadata_dir / f"{OUTPUT_STEM}_primary_predictions.tsv",
        "sensitivity_predictions": args.metadata_dir / f"{OUTPUT_STEM}_sensitivity_predictions.tsv",
        "consensus_predictions": args.metadata_dir / f"{OUTPUT_STEM}_consensus_predictions.tsv",
        "gene_coverage": args.metadata_dir / f"{OUTPUT_STEM}_gene_coverage.tsv",
        "report": args.metadata_dir / f"{OUTPUT_STEM}_report.json",
    }
    prediction_frames["primary"].to_csv(outputs["primary_predictions"], sep="\t", index=False)
    if "sensitivity" in prediction_frames:
        prediction_frames["sensitivity"].to_csv(outputs["sensitivity_predictions"], sep="\t", index=False)
    else:
        pd.DataFrame().to_csv(outputs["sensitivity_predictions"], sep="\t", index=False)

    consensus = build_consensus_table(prediction_frames["primary"], prediction_frames.get("sensitivity"))
    consensus.to_csv(outputs["consensus_predictions"], sep="\t", index=False)

    coverage = pd.DataFrame(coverage_rows)
    coverage.to_csv(outputs["gene_coverage"], sep="\t", index=False)

    report = {
        "module": "module9_7_drugreflector_inference",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "completed_drugreflector_inference",
        "inputs": {
            "primary_vscore": str(args.primary_vscore.resolve()),
            "sensitivity_vscore": str(args.sensitivity_vscore.resolve()),
            "checkpoint_dir": str(args.checkpoint_dir.resolve()),
            "drugreflector_source_dir": str(args.drugreflector_source_dir.resolve()),
            "top_n": args.top_n,
            "seed": args.seed,
            "skip_sensitivity": bool(args.skip_sensitivity),
        },
        "outputs": {key: str(path.resolve()) for key, path in outputs.items()},
        "summary": {
            "n_compounds_total": int(model.n_compounds),
            "n_primary_predictions": int(len(prediction_frames["primary"])),
            "n_sensitivity_predictions": int(len(prediction_frames.get("sensitivity", pd.DataFrame()))),
            "n_consensus_predictions": int(len(consensus)),
            "n_compounds_in_both_top_lists": int(consensus["in_both_top_lists"].sum()),
            "primary_union_overlap_genes": int(
                coverage.loc[(coverage["signature"].eq("primary")) & (coverage["fold"].astype(str).eq("union")), "n_overlap_genes"].iloc[0]
            ),
            "sensitivity_union_overlap_genes": int(
                coverage.loc[(coverage["signature"].eq("sensitivity")) & (coverage["fold"].astype(str).eq("union")), "n_overlap_genes"].iloc[0]
            )
            if "sensitivity" in prediction_frames
            else None,
        },
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": package_version("numpy"),
            "pandas": package_version("pandas"),
            "torch": package_version("torch"),
            "drugreflector": package_version("drugreflector"),
        },
    }
    write_json(outputs["report"], report)
    print(json.dumps(report["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
