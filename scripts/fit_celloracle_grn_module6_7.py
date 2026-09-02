from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

try:
    from qc_celloracle_inputs_module6_5 import read_tf_list
except ModuleNotFoundError:
    from scripts.qc_celloracle_inputs_module6_5 import read_tf_list


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ORACLE = PROJECT_ROOT / "data/processed/driver/celloracle_module6_6/celloracle_module6_6.celloracle.oracle"
DEFAULT_TF_LIST = PROJECT_ROOT / "metadata/driver/celloracle_input_tfs.module6_4.txt"
DEFAULT_OUT_DIR = PROJECT_ROOT / "data/processed/driver/celloracle_module6_7"
DEFAULT_METADATA_DIR = PROJECT_ROOT / "metadata/driver"


def merge_links_dict(links_dict: dict[str, pd.DataFrame]) -> pd.DataFrame:
    frames = []
    for state, table in links_dict.items():
        df = table.copy()
        df.insert(0, "celloracle_state", state)
        frames.append(df)
    if not frames:
        return pd.DataFrame(
            columns=["celloracle_state", "source", "target", "coef_mean", "coef_abs", "p", "-logp"]
        )
    merged = pd.concat(frames, axis=0, ignore_index=True)
    expected = ["celloracle_state", "source", "target", "coef_mean", "coef_abs", "p", "-logp"]
    missing = [col for col in expected if col not in merged.columns]
    if missing:
        raise ValueError(f"Links table missing required columns: {missing}")
    return merged


def summarize_grn_links(links: pd.DataFrame, p_threshold: float) -> pd.DataFrame:
    rows = []
    for state, df in links.groupby("celloracle_state", observed=True):
        passing = df["p"] <= p_threshold
        rows.append(
            {
                "celloracle_state": state,
                "n_edges_total": int(len(df)),
                "n_edges_passing_p": int(passing.sum()),
                "n_source_tfs_total": int(df["source"].nunique()),
                "n_target_genes_total": int(df["target"].nunique()),
                "n_source_tfs_passing_p": int(df.loc[passing, "source"].nunique()),
                "n_target_genes_passing_p": int(df.loc[passing, "target"].nunique()),
                "median_coef_abs": float(df["coef_abs"].median()) if len(df) else 0.0,
                "median_coef_abs_passing_p": float(df.loc[passing, "coef_abs"].median()) if passing.any() else 0.0,
            }
        )
    return pd.DataFrame(rows).sort_values("celloracle_state").reset_index(drop=True)


def summarize_tf_network(links: pd.DataFrame, input_tfs: list[str], p_threshold: float) -> pd.DataFrame:
    states = sorted(links["celloracle_state"].astype(str).unique())
    rows = []
    for state in states:
        state_df = links.loc[links["celloracle_state"].astype(str) == state]
        for tf in input_tfs:
            tf_df = state_df.loc[state_df["source"].astype(str) == tf]
            passing = tf_df["p"] <= p_threshold if len(tf_df) else pd.Series(dtype=bool)
            passing_df = tf_df.loc[passing] if len(tf_df) else tf_df
            rows.append(
                {
                    "celloracle_state": state,
                    "tf": tf,
                    "n_edges_total": int(len(tf_df)),
                    "n_edges_passing_p": int(len(passing_df)),
                    "n_positive_edges_passing_p": int((passing_df["coef_mean"] > 0).sum()) if len(passing_df) else 0,
                    "n_negative_edges_passing_p": int((passing_df["coef_mean"] < 0).sum()) if len(passing_df) else 0,
                    "n_target_genes_passing_p": int(passing_df["target"].nunique()) if len(passing_df) else 0,
                    "mean_coef_abs_passing_p": float(passing_df["coef_abs"].mean()) if len(passing_df) else 0.0,
                    "max_coef_abs_passing_p": float(passing_df["coef_abs"].max()) if len(passing_df) else 0.0,
                }
            )
    return pd.DataFrame(rows)


def export_filtered_links(links, merged_links: pd.DataFrame, p_threshold: float, threshold_number: int) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    links.filter_links(p=p_threshold, weight="coef_abs", threshold_number=threshold_number)
    filtered = merge_links_dict(links.filtered_links)
    network_score = None
    try:
        links.get_network_score()
        network_score = links.merged_score.copy()
    except Exception:
        network_score = None
    return filtered, network_score


def run_celloracle_grn(
    oracle_path: Path,
    tf_list_path: Path,
    out_dir: Path,
    metadata_dir: Path,
    alpha_links: float,
    alpha_simulation: float,
    bagging_number: int,
    n_pca_dims: int,
    k: int,
    n_jobs: int,
    p_threshold: float,
    threshold_number: int,
    test_mode: bool,
    skip_simulation_fit: bool,
) -> dict:
    import celloracle as co

    input_tfs = read_tf_list(tf_list_path)
    out_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)

    oracle = co.load_hdf5(str(oracle_path))
    oracle.perform_PCA(n_components=n_pca_dims)
    oracle.knn_imputation(
        n_pca_dims=n_pca_dims,
        k=k,
        balanced=True,
        b_sight=3000,
        b_maxl=1500,
        n_jobs=n_jobs,
    )
    preprocessed_oracle_path = out_dir / "celloracle_module6_7_preprocessed.celloracle.oracle"
    oracle.to_hdf5(str(preprocessed_oracle_path))

    links = oracle.get_links(
        cluster_name_for_GRN_unit="celloracle_state",
        alpha=alpha_links,
        bagging_number=bagging_number,
        verbose_level=1,
        test_mode=test_mode,
        model_method="bagging_ridge",
        n_jobs=n_jobs,
    )
    links_path = out_dir / "celloracle_module6_7_state_specific_grn.celloracle.links"
    links.to_hdf5(str(links_path))

    merged_links = merge_links_dict(links.links_dict)
    raw_links_path = metadata_dir / "celloracle_module6_7_grn_links_raw.tsv.gz"
    merged_links.to_csv(raw_links_path, sep="\t", index=False)

    filtered_links, network_score = export_filtered_links(
        links,
        merged_links=merged_links,
        p_threshold=p_threshold,
        threshold_number=threshold_number,
    )
    filtered_links_path = metadata_dir / "celloracle_module6_7_grn_links_filtered.tsv.gz"
    filtered_links.to_csv(filtered_links_path, sep="\t", index=False)

    grn_summary = summarize_grn_links(merged_links, p_threshold=p_threshold)
    grn_summary_path = metadata_dir / "celloracle_module6_7_grn_state_summary.tsv"
    grn_summary.to_csv(grn_summary_path, sep="\t", index=False)

    tf_summary = summarize_tf_network(merged_links, input_tfs=input_tfs, p_threshold=p_threshold)
    tf_summary_path = metadata_dir / "celloracle_module6_7_tf_network_summary.tsv"
    tf_summary.to_csv(tf_summary_path, sep="\t", index=False)

    network_score_path = None
    if network_score is not None:
        network_score_path = metadata_dir / "celloracle_module6_7_network_scores.tsv.gz"
        network_score.to_csv(network_score_path, sep="\t", index=False)

    fitted_oracle_path = None
    if not skip_simulation_fit:
        oracle.fit_GRN_for_simulation(
            GRN_unit="cluster",
            alpha=alpha_simulation,
            verbose_level=1,
        )
        fitted_oracle_path = out_dir / "celloracle_module6_7_fitted.celloracle.oracle"
        oracle.to_hdf5(str(fitted_oracle_path))

    return {
        "celloracle_version": getattr(co, "__version__", None),
        "input_oracle": str(oracle_path),
        "preprocessed_oracle": str(preprocessed_oracle_path),
        "links_object": str(links_path),
        "fitted_oracle": str(fitted_oracle_path) if fitted_oracle_path else None,
        "raw_links": str(raw_links_path),
        "filtered_links": str(filtered_links_path),
        "grn_state_summary": str(grn_summary_path),
        "tf_network_summary": str(tf_summary_path),
        "network_scores": str(network_score_path) if network_score_path else None,
        "n_states": int(len(links.cluster)),
        "states": list(map(str, links.cluster)),
        "n_edges_raw": int(len(merged_links)),
        "n_edges_filtered": int(len(filtered_links)),
        "n_input_tfs": int(len(input_tfs)),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 6.7 fit CellOracle state-specific GRN")
    parser.add_argument("--oracle", type=Path, default=DEFAULT_ORACLE)
    parser.add_argument("--tf-list", type=Path, default=DEFAULT_TF_LIST)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA_DIR)
    parser.add_argument("--alpha-links", type=float, default=10.0)
    parser.add_argument("--alpha-simulation", type=float, default=1.0)
    parser.add_argument("--bagging-number", type=int, default=20)
    parser.add_argument("--n-pca-dims", type=int, default=30)
    parser.add_argument("--k", type=int, default=30)
    parser.add_argument("--n-jobs", type=int, default=8)
    parser.add_argument("--p-threshold", type=float, default=0.001)
    parser.add_argument("--threshold-number", type=int, default=10000)
    parser.add_argument("--test-mode", action="store_true")
    parser.add_argument("--skip-simulation-fit", action="store_true")
    parser.add_argument("--report", type=Path, default=DEFAULT_METADATA_DIR / "celloracle_module6_7_grn_report.json")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = datetime.now(timezone.utc)
    result = run_celloracle_grn(
        oracle_path=args.oracle,
        tf_list_path=args.tf_list,
        out_dir=args.out_dir,
        metadata_dir=args.metadata_dir,
        alpha_links=args.alpha_links,
        alpha_simulation=args.alpha_simulation,
        bagging_number=args.bagging_number,
        n_pca_dims=args.n_pca_dims,
        k=args.k,
        n_jobs=args.n_jobs,
        p_threshold=args.p_threshold,
        threshold_number=args.threshold_number,
        test_mode=args.test_mode,
        skip_simulation_fit=args.skip_simulation_fit,
    )
    finished = datetime.now(timezone.utc)
    report = {
        "module": "6.7",
        "method": "CellOracle state-specific GRN fitting",
        "created_at_utc": finished.isoformat(),
        "elapsed_seconds": (finished - started).total_seconds(),
        "parameters": {
            "alpha_links": args.alpha_links,
            "alpha_simulation": args.alpha_simulation,
            "bagging_number": args.bagging_number,
            "n_pca_dims": args.n_pca_dims,
            "k": args.k,
            "n_jobs": args.n_jobs,
            "p_threshold": args.p_threshold,
            "threshold_number": args.threshold_number,
            "test_mode": args.test_mode,
            "skip_simulation_fit": args.skip_simulation_fit,
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
            "n_states": result["n_states"],
            "states": result["states"],
            "n_edges_raw": result["n_edges_raw"],
            "n_edges_filtered": result["n_edges_filtered"],
            "links_object": result["links_object"],
            "fitted_oracle": result["fitted_oracle"],
            "report": str(args.report),
        },
        ensure_ascii=False,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
