from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
TRAJ_DIR = ROOT / "data/processed/trajectory"
META_DIR = ROOT / "metadata/trajectory"


def finite_spearman(left: pd.Series, right: pd.Series) -> tuple[float, int]:
    df = pd.DataFrame({"left": left, "right": right}).replace([np.inf, -np.inf], np.nan).dropna()
    if df.shape[0] < 3 or df["left"].nunique() < 2 or df["right"].nunique() < 2:
        return float("nan"), int(df.shape[0])
    return float(spearmanr(df["left"], df["right"]).correlation), int(df.shape[0])


def normalize_series(values: pd.Series) -> pd.Series:
    values = pd.to_numeric(values, errors="coerce")
    finite = values[np.isfinite(values)]
    if finite.empty:
        return values * np.nan
    low = finite.min()
    high = finite.max()
    if high == low:
        return values * 0
    return (values - low) / (high - low)


def read_monocle(run_id: str) -> pd.DataFrame:
    path = TRAJ_DIR / "module5_3_quick5k_counts" / run_id / "monocle3" / "monocle3_pseudotime.tsv"
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={"monocle3_pseudotime": f"{run_id}__monocle3"})
    df[f"{run_id}__monocle3_norm"] = normalize_series(df[f"{run_id}__monocle3"])
    return df


def read_slingshot(run_id: str, embedding: str) -> pd.DataFrame:
    path = TRAJ_DIR / "module5_3_quick5k" / run_id / f"slingshot_{embedding}" / "slingshot_pseudotime.tsv"
    df = pd.read_csv(path, sep="\t")
    lineage_cols = [col for col in df.columns if col != "cell_id"]
    values = df[lineage_cols].apply(pd.to_numeric, errors="coerce")
    df[f"{run_id}__slingshot_{embedding}"] = values.mean(axis=1, skipna=True)
    df[f"{run_id}__slingshot_{embedding}_finite_lineages"] = values.notna().sum(axis=1)
    df[f"{run_id}__slingshot_{embedding}_norm"] = normalize_series(df[f"{run_id}__slingshot_{embedding}"])
    return df[["cell_id", f"{run_id}__slingshot_{embedding}", f"{run_id}__slingshot_{embedding}_norm", f"{run_id}__slingshot_{embedding}_finite_lineages"]]


def stage_monotonicity(df: pd.DataFrame, method_col: str) -> tuple[float, int]:
    stage_order = {
        "stage_0_reference_hepatocyte": 0,
        "stage_1_stressed_injured": 1,
        "stage_2_regenerative_progenitor": 2,
        "stage_3_proliferating_candidate": 3,
        "stage_4_cnv_supported_malignant": 4,
        "stage_4_malignant_like_review": 4,
        "unresolved_epithelial_or_mixed": np.nan,
    }
    sub = df[["cell_disease_stage", method_col]].copy()
    sub["stage_order"] = sub["cell_disease_stage"].map(stage_order)
    grouped = sub.dropna().groupby("stage_order", observed=True)[method_col].mean().reset_index()
    rho, n = finite_spearman(grouped["stage_order"], grouped[method_col])
    return rho, n


def method_status_rows() -> list[dict[str, object]]:
    rows = []
    for base, method, status_file in [
        ("module5_3_quick5k_counts", "monocle3", "monocle3/monocle3_status.json"),
        ("module5_3_quick5k", "slingshot_scanvi", "slingshot_scanvi/slingshot_status.json"),
        ("module5_3_quick5k", "slingshot_hepatocyte_pca", "slingshot_hepatocyte_pca/slingshot_status.json"),
    ]:
        for run_id in ["main_strict", "sensitivity_include_review"]:
            path = TRAJ_DIR / base / run_id / status_file
            if path.exists():
                status = json.loads(path.read_text(encoding="utf-8"))
                rows.append({"run_id": run_id, "method": method, **status})
            else:
                rows.append({"run_id": run_id, "method": method, "status": "missing_output", "message": str(path)})
    return rows


def main() -> int:
    META_DIR.mkdir(parents=True, exist_ok=True)
    metrics = []
    merged_outputs = {}

    for run_id in ["main_strict", "sensitivity_include_review"]:
        mono = read_monocle(run_id)
        scanvi = read_slingshot(run_id, "scanvi")
        pca = read_slingshot(run_id, "hepatocyte_pca")
        merged = mono.merge(scanvi, on="cell_id", how="outer").merge(pca, on="cell_id", how="outer")
        merged_outputs[run_id] = merged
        merged.to_csv(META_DIR / f"trajectory_module5_3_{run_id}_pseudotime_merged.tsv.gz", sep="\t", index=False, compression="gzip")

        pairs = [
            ("monocle3_vs_slingshot_scanvi", f"{run_id}__monocle3_norm", f"{run_id}__slingshot_scanvi_norm"),
            ("monocle3_vs_slingshot_hepatocyte_pca", f"{run_id}__monocle3_norm", f"{run_id}__slingshot_hepatocyte_pca_norm"),
            ("slingshot_scanvi_vs_hepatocyte_pca", f"{run_id}__slingshot_scanvi_norm", f"{run_id}__slingshot_hepatocyte_pca_norm"),
        ]
        for comparison, left, right in pairs:
            rho, n = finite_spearman(merged[left], merged[right])
            metrics.append({"run_id": run_id, "metric": comparison, "spearman_rho": rho, "n_common_cells": n})

        for method_col in [
            f"{run_id}__monocle3_norm",
            f"{run_id}__slingshot_scanvi_norm",
            f"{run_id}__slingshot_hepatocyte_pca_norm",
        ]:
            rho, n = stage_monotonicity(merged, method_col)
            metrics.append({"run_id": run_id, "metric": f"stage_monotonicity__{method_col}", "spearman_rho": rho, "n_common_cells": n})

    common = merged_outputs["main_strict"].merge(
        merged_outputs["sensitivity_include_review"],
        on="cell_id",
        how="inner",
        suffixes=("", "_review"),
    )
    for method in ["monocle3", "slingshot_scanvi", "slingshot_hepatocyte_pca"]:
        rho, n = finite_spearman(common[f"main_strict__{method}_norm"], common[f"sensitivity_include_review__{method}_norm"])
        metrics.append({"run_id": "main_vs_include_review", "metric": f"variant_stability__{method}", "spearman_rho": rho, "n_common_cells": n})

    metrics_df = pd.DataFrame(metrics)
    metrics_path = META_DIR / "trajectory_module5_3_sensitivity_metrics.tsv"
    metrics_df.to_csv(metrics_path, sep="\t", index=False)

    status_df = pd.DataFrame(method_status_rows())
    status_path = META_DIR / "trajectory_module5_3_method_status.tsv"
    status_df.to_csv(status_path, sep="\t", index=False)

    report = {
        "module": "5.3",
        "method": "Monocle3 and Slingshot trajectory modeling with sensitivity validation",
        "metrics": str(metrics_path.resolve()),
        "method_status": str(status_path.resolve()),
        "merged_outputs": {
            run_id: str((META_DIR / f"trajectory_module5_3_{run_id}_pseudotime_merged.tsv.gz").resolve())
            for run_id in merged_outputs
        },
        "status_counts": status_df["status"].value_counts().to_dict(),
        "n_metrics": int(metrics_df.shape[0]),
    }
    report_path = META_DIR / "trajectory_module5_3_summary_report.json"
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
