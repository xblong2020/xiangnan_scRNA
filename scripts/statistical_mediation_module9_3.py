from __future__ import annotations

import argparse
import json
import math
import platform
import sys
import time
from collections.abc import Iterable, Sequence
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts import external_validation_module8 as module8
DEFAULT_METADATA_DIR = ROOT / "metadata/driver"
DEFAULT_FIGURE_DIR = ROOT / "figures/driver"
DEFAULT_MODULE8_REPORT = DEFAULT_METADATA_DIR / "module8_external_validation_report.json"
DEFAULT_MODULE9_1_CELL_SCORES = DEFAULT_METADATA_DIR / "module9_1_temporal_cell_scores.tsv.gz"
DEFAULT_STAGE_BY_SAMPLE = ROOT / "metadata/trajectory/trajectory_module5_2_stage_by_sample.tsv"
DEFAULT_GTF_PATH = Path(
    r"G:\wanyi_HCC_scRNA\HCCscRNA\GSE156625-HCC\cellranger\hg38\refdata-gex-GRCh38-2020-A\genes\genes.gtf"
)

AXIS_MAP = {
    "tier1_rescue": "A_hnf4a_ppara_loss",
    "ap1_stress_proliferation": "B_transition_activation",
    "sox4_state_specific": "C_sox4_axis",
}
PATH_COLUMNS = ["A_hnf4a_ppara_loss", "B_transition_activation", "C_sox4_axis"]
INDIRECT_PATHS = [
    "A_to_B_to_outcome",
    "A_to_C_to_outcome",
    "B_to_C_to_outcome",
    "A_to_B_to_C_to_outcome",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 9.3 statistical mediation/path-analysis evidence.")
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--module8-report", type=Path, default=DEFAULT_MODULE8_REPORT)
    parser.add_argument("--module9-1-cell-scores", type=Path, default=DEFAULT_MODULE9_1_CELL_SCORES)
    parser.add_argument("--stage-by-sample", type=Path, default=DEFAULT_STAGE_BY_SAMPLE)
    parser.add_argument("--gtf-path", type=Path, default=DEFAULT_GTF_PATH)
    parser.add_argument("--n-bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260615)
    return parser.parse_args()


def package_version(name: str) -> str:
    try:
        return version(name)
    except PackageNotFoundError:
        return "not_installed"


def read_tsv_or_empty(path: Path, **kwargs) -> pd.DataFrame:
    if path is None or not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", **kwargs)


def safe_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def zscore(series: pd.Series) -> pd.Series:
    values = safe_numeric(series)
    sd = values.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.zeros(len(values)), index=series.index, dtype=float)
    return (values - values.mean()) / sd


def benjamini_hochberg(pvalues: Iterable[float]) -> list[float]:
    p = np.asarray([1.0 if pd.isna(v) else float(v) for v in pvalues], dtype=float)
    n = len(p)
    if n == 0:
        return []
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.empty(n, dtype=float)
    running = 1.0
    for i in range(n - 1, -1, -1):
        running = min(running, ranked[i] * n / (i + 1))
        adjusted[order[i]] = running
    return adjusted.clip(0, 1).tolist()


def build_bulk_axis_scores(scores: pd.DataFrame) -> pd.DataFrame:
    if scores.empty:
        return pd.DataFrame()
    required = {"dataset_id", "sample", "axis", "signature_score"}
    missing = required.difference(scores.columns)
    if missing:
        raise ValueError(f"bulk signature score table missing columns: {sorted(missing)}")

    table = scores.loc[scores["axis"].isin(AXIS_MAP)].copy()
    if "sample_type" in table.columns:
        table = table.loc[table["sample_type"].astype(str).eq("tumor")].copy()
    wide = table.pivot_table(
        index=["dataset_id", "sample"],
        columns="axis",
        values="signature_score",
        aggfunc="mean",
    ).reset_index()
    if "sample_type" in table.columns:
        types = table.groupby(["dataset_id", "sample"], as_index=False)["sample_type"].first()
        wide = wide.merge(types, on=["dataset_id", "sample"], how="left")
    for axis, out_col in AXIS_MAP.items():
        if axis not in wide.columns:
            wide[axis] = np.nan
    rows = []
    for dataset_id, group in wide.groupby("dataset_id", dropna=False):
        out = group[["dataset_id", "sample"]].copy()
        out["sample_type"] = group["sample_type"].values if "sample_type" in group.columns else "tumor"
        out["A_hnf4a_ppara_loss"] = -zscore(group["tier1_rescue"]).values
        out["B_transition_activation"] = zscore(group["ap1_stress_proliferation"]).values
        out["C_sox4_axis"] = zscore(group["sox4_state_specific"]).values
        rows.append(out)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def stage_to_ordinal(value: object) -> float:
    if pd.isna(value):
        return np.nan
    text = str(value).strip().lower()
    if text in {"", "na", "nan", "unknown"}:
        return np.nan
    if "metastatic" in text:
        return 3.0
    if "pvtt" in text:
        return 2.0
    if "primary" in text or "hcc_dataset" in text:
        return 1.0
    if "chronic" in text or "reference" in text or "adjacent" in text or "non_hcc" in text:
        return 0.0
    roman = {"i": 1, "ii": 2, "iii": 3, "iv": 4}
    for key, value_ in roman.items():
        if f"stage {key}" in text or f"stage_{key}" in text:
            return float(value_)
    digits = "".join(ch for ch in text if ch.isdigit())
    return float(digits[0]) if digits else np.nan


def build_scrna_pseudobulk_scores(cell_scores: pd.DataFrame, stage_by_sample: pd.DataFrame) -> pd.DataFrame:
    if cell_scores.empty:
        return pd.DataFrame()
    required = {"run_id", "method", "dataset", "cnv_sample", *PATH_COLUMNS, "C_malignant_like_fate", "cell_disease_stage"}
    missing = required.difference(cell_scores.columns)
    if missing:
        raise ValueError(f"Module 9.1 cell score table missing columns: {sorted(missing)}")
    table = cell_scores.copy()
    table["is_malignant_cell"] = table["cell_disease_stage"].astype(str).str.contains("stage_4|malignant", case=False, regex=True)
    group_cols = ["run_id", "method", "dataset", "cnv_sample"]
    aggregated = (
        table.groupby(group_cols, dropna=False)
        .agg(
            A_hnf4a_ppara_loss=("A_hnf4a_ppara_loss", "mean"),
            B_transition_activation=("B_transition_activation", "mean"),
            C_sox4_axis=("C_sox4_axis", "mean"),
            mean_C_malignant_like_fate=("C_malignant_like_fate", "mean"),
            malignant_fraction=("is_malignant_cell", "mean"),
            n_cells=("cell_disease_stage", "size"),
        )
        .reset_index()
    )
    if not stage_by_sample.empty and {"dataset", "cnv_sample", "sample_disease_stage"}.issubset(stage_by_sample.columns):
        sample_stage = (
            stage_by_sample[["dataset", "cnv_sample", "sample_disease_stage"]]
            .dropna(subset=["cnv_sample"])
            .drop_duplicates(["dataset", "cnv_sample"])
        )
        aggregated = aggregated.merge(sample_stage, on=["dataset", "cnv_sample"], how="left")
    else:
        aggregated["sample_disease_stage"] = np.nan
    aggregated["sample_stage_ordinal"] = aggregated["sample_disease_stage"].map(stage_to_ordinal)
    aggregated["analysis_unit"] = (
        aggregated["run_id"].astype(str)
        + "|"
        + aggregated["method"].astype(str)
        + "|"
        + aggregated["dataset"].astype(str)
        + "|"
        + aggregated["cnv_sample"].astype(str)
    )
    return aggregated


def _fit_ols(y: pd.Series, x: pd.DataFrame) -> tuple[pd.Series, pd.Series, int, str]:
    model_df = pd.concat([y.rename("y"), x], axis=1).dropna()
    if len(model_df) < max(8, x.shape[1] + 3) or model_df["y"].nunique() < 2:
        return pd.Series(dtype=float), pd.Series(dtype=float), int(len(model_df)), "not_testable_insufficient_variation"
    exog = sm.add_constant(model_df.drop(columns=["y"]), has_constant="add")
    try:
        result = sm.OLS(model_df["y"], exog).fit()
        return result.params, result.pvalues, int(len(model_df)), "tested"
    except Exception:
        return pd.Series(dtype=float), pd.Series(dtype=float), int(len(model_df)), "model_failed"


def compute_path_coefficients(
    data: pd.DataFrame,
    outcome_col: str,
    outcome_type: str = "continuous",
    covariates: Sequence[str] | None = None,
    dataset_id: str | None = None,
) -> pd.DataFrame:
    covariates = [cov for cov in (covariates or []) if cov in data.columns]
    if dataset_id is None:
        dataset_id = str(data["dataset_id"].iloc[0]) if "dataset_id" in data.columns and len(data) else "unknown"
    rows: list[dict[str, object]] = []
    base_cols = [*PATH_COLUMNS, outcome_col, *covariates]
    if outcome_type == "cox" and "event" in data.columns:
        base_cols.append("event")
    df = data[[col for col in base_cols if col in data.columns]].copy()
    for col in df.columns:
        df[col] = safe_numeric(df[col])

    params, pvalues, n, status = _fit_ols(df["B_transition_activation"], df[["A_hnf4a_ppara_loss"] + covariates])
    rows.append(
        {
            "dataset_id": dataset_id,
            "outcome": outcome_col,
            "outcome_type": outcome_type,
            "path": "A_to_B",
            "predictor": "A_hnf4a_ppara_loss",
            "response": "B_transition_activation",
            "coef": float(params.get("A_hnf4a_ppara_loss", np.nan)),
            "pvalue": float(pvalues.get("A_hnf4a_ppara_loss", np.nan)),
            "n_samples": n,
            "model_type": "ols",
            "status": status,
        }
    )

    params, pvalues, n, status = _fit_ols(
        df["C_sox4_axis"], df[["A_hnf4a_ppara_loss", "B_transition_activation"] + covariates]
    )
    for path, predictor in [("A_to_C", "A_hnf4a_ppara_loss"), ("B_to_C", "B_transition_activation")]:
        rows.append(
            {
                "dataset_id": dataset_id,
                "outcome": outcome_col,
                "outcome_type": outcome_type,
                "path": path,
                "predictor": predictor,
                "response": "C_sox4_axis",
                "coef": float(params.get(predictor, np.nan)),
                "pvalue": float(pvalues.get(predictor, np.nan)),
                "n_samples": n,
                "model_type": "ols",
                "status": status,
            }
        )

    if outcome_type == "cox":
        cox = fit_cox_path_model(df.rename(columns={outcome_col: "time"}), time_col="time", event_col="event", covariates=covariates)
        for _, row in cox.iterrows():
            path = {
                "A_hnf4a_ppara_loss": "A_to_outcome",
                "B_transition_activation": "B_to_outcome",
                "C_sox4_axis": "C_to_outcome",
            }.get(row["predictor"], f"{row['predictor']}_to_outcome")
            rows.append(
                {
                    "dataset_id": dataset_id,
                    "outcome": outcome_col,
                    "outcome_type": outcome_type,
                    "path": path,
                    "predictor": row["predictor"],
                    "response": outcome_col,
                    "coef": row["coef"],
                    "pvalue": row["pvalue"],
                    "n_samples": row["n_samples"],
                    "model_type": row["model_type"],
                    "status": row["status"],
                }
            )
    else:
        params, pvalues, n, status = _fit_ols(
            df[outcome_col], df[["A_hnf4a_ppara_loss", "B_transition_activation", "C_sox4_axis"] + covariates]
        )
        for path, predictor in [
            ("A_to_outcome", "A_hnf4a_ppara_loss"),
            ("B_to_outcome", "B_transition_activation"),
            ("C_to_outcome", "C_sox4_axis"),
        ]:
            rows.append(
                {
                    "dataset_id": dataset_id,
                    "outcome": outcome_col,
                    "outcome_type": outcome_type,
                    "path": path,
                    "predictor": predictor,
                    "response": outcome_col,
                    "coef": float(params.get(predictor, np.nan)),
                    "pvalue": float(pvalues.get(predictor, np.nan)),
                    "n_samples": n,
                    "model_type": "ols",
                    "status": status,
                }
            )
    return pd.DataFrame(rows)


def fit_cox_path_model(
    data: pd.DataFrame,
    time_col: str,
    event_col: str,
    covariates: Sequence[str] | None = None,
) -> pd.DataFrame:
    covariates = [cov for cov in (covariates or []) if cov in data.columns]
    predictors = [*PATH_COLUMNS, *covariates]
    columns = ["predictor", "coef", "hazard_ratio", "pvalue", "n_samples", "n_events", "model_type", "status"]
    if time_col not in data.columns or event_col not in data.columns:
        return pd.DataFrame([{col: np.nan for col in columns} | {"status": "not_testable_missing_time_or_event"}])
    df = data[[time_col, event_col] + predictors].copy()
    for col in df.columns:
        df[col] = safe_numeric(df[col])
    df = df.dropna()
    n = int(len(df))
    n_events = int(df[event_col].sum()) if n else 0
    if n < max(12, len(predictors) + 5) or n_events < 5 or df[event_col].nunique() < 2:
        return pd.DataFrame(
            [
                {
                    "predictor": PATH_COLUMNS[0],
                    "coef": np.nan,
                    "hazard_ratio": np.nan,
                    "pvalue": np.nan,
                    "n_samples": n,
                    "n_events": n_events,
                    "model_type": "cox_phreg",
                    "status": "not_testable_insufficient_events",
                }
            ]
        )
    try:
        from statsmodels.duration.hazard_regression import PHReg

        exog = df[predictors]
        keep = [col for col in exog.columns if exog[col].std(ddof=0) > 0]
        exog = exog[keep]
        result = PHReg(df[time_col], exog, status=df[event_col]).fit(disp=False)
        rows = []
        for idx, predictor in enumerate(keep):
            rows.append(
                {
                    "predictor": predictor,
                    "coef": float(result.params[idx]),
                    "hazard_ratio": float(np.exp(result.params[idx])),
                    "pvalue": float(result.pvalues[idx]),
                    "n_samples": n,
                    "n_events": n_events,
                    "model_type": "cox_phreg",
                    "status": "tested",
                }
            )
        return pd.DataFrame(rows, columns=columns)
    except Exception:
        return pd.DataFrame(
            [
                {
                    "predictor": PATH_COLUMNS[0],
                    "coef": np.nan,
                    "hazard_ratio": np.nan,
                    "pvalue": np.nan,
                    "n_samples": n,
                    "n_events": n_events,
                    "model_type": "cox_phreg",
                    "status": "model_failed",
                }
            ]
        )


def _coefficient_lookup(coefficients: pd.DataFrame) -> dict[str, float]:
    lookup = coefficients.set_index("path")["coef"].to_dict() if not coefficients.empty else {}
    return {key: float(value) for key, value in lookup.items() if pd.notna(value)}


def _indirect_from_coefficients(coefficients: pd.DataFrame) -> dict[str, float]:
    coef = _coefficient_lookup(coefficients)
    return {
        "A_to_B_to_outcome": coef.get("A_to_B", np.nan) * coef.get("B_to_outcome", np.nan),
        "A_to_C_to_outcome": coef.get("A_to_C", np.nan) * coef.get("C_to_outcome", np.nan),
        "B_to_C_to_outcome": coef.get("B_to_C", np.nan) * coef.get("C_to_outcome", np.nan),
        "A_to_B_to_C_to_outcome": coef.get("A_to_B", np.nan) * coef.get("B_to_C", np.nan) * coef.get("C_to_outcome", np.nan),
    }


def bootstrap_indirect_effects(
    data: pd.DataFrame,
    outcome_col: str,
    outcome_type: str = "continuous",
    covariates: Sequence[str] | None = None,
    n_bootstrap: int = 1000,
    random_state: int = 0,
    dataset_id: str | None = None,
) -> pd.DataFrame:
    dataset_id = dataset_id or (str(data["dataset_id"].iloc[0]) if "dataset_id" in data.columns and len(data) else "unknown")
    point_coefficients = compute_path_coefficients(data, outcome_col, outcome_type, covariates, dataset_id)
    point = _indirect_from_coefficients(point_coefficients)
    rng = np.random.default_rng(random_state)
    boot: dict[str, list[float]] = {path: [] for path in INDIRECT_PATHS}
    usable = data.dropna(subset=[col for col in [*PATH_COLUMNS, outcome_col, *(covariates or [])] if col in data.columns]).copy()
    if len(usable) < 8:
        return pd.DataFrame(
            [
                {
                    "dataset_id": dataset_id,
                    "outcome": outcome_col,
                    "outcome_type": outcome_type,
                    "indirect_path": path,
                    "effect": point.get(path, np.nan),
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "pvalue": np.nan,
                    "p.adjust": np.nan,
                    "n_bootstrap": 0,
                    "n_successful_bootstrap": 0,
                    "status": "not_testable_insufficient_samples",
                }
                for path in INDIRECT_PATHS
            ]
        )
    for _ in range(max(0, int(n_bootstrap))):
        sample = usable.iloc[rng.integers(0, len(usable), len(usable))].copy()
        try:
            coefs = compute_path_coefficients(sample, outcome_col, outcome_type, covariates, dataset_id)
            effects = _indirect_from_coefficients(coefs)
        except Exception:
            continue
        for path, value in effects.items():
            if np.isfinite(value):
                boot[path].append(float(value))
    rows = []
    pvalues = []
    for path in INDIRECT_PATHS:
        values = np.asarray(boot[path], dtype=float)
        effect = float(point.get(path, np.nan))
        if len(values):
            ci_low, ci_high = np.percentile(values, [2.5, 97.5])
            pvalue = float(2 * min(np.mean(values <= 0), np.mean(values >= 0)))
            pvalue = max(pvalue, 1.0 / max(len(values), 1)) if pvalue == 0 else pvalue
            status = "tested"
        else:
            ci_low = ci_high = pvalue = np.nan
            status = "bootstrap_failed"
        pvalues.append(pvalue)
        rows.append(
            {
                "dataset_id": dataset_id,
                "outcome": outcome_col,
                "outcome_type": outcome_type,
                "indirect_path": path,
                "effect": effect,
                "ci_low": float(ci_low) if np.isfinite(ci_low) else np.nan,
                "ci_high": float(ci_high) if np.isfinite(ci_high) else np.nan,
                "pvalue": pvalue,
                "n_bootstrap": int(n_bootstrap),
                "n_successful_bootstrap": int(len(values)),
                "status": status,
            }
        )
    adjusted = benjamini_hochberg(pvalues)
    for row, padj in zip(rows, adjusted):
        row["p.adjust"] = padj
    return pd.DataFrame(rows)


def assess_outcome_availability(data: pd.DataFrame, outcomes: Sequence[str], min_n: int = 3) -> pd.DataFrame:
    rows = []
    for outcome in outcomes:
        if outcome not in data.columns:
            status = "not_testable_missing_column"
            n_observed = 0
            n_unique = 0
        else:
            values = safe_numeric(data[outcome]).dropna()
            n_observed = int(len(values))
            n_unique = int(values.nunique())
            if n_observed == 0:
                status = "not_testable_no_observed_values"
            elif n_observed < min_n:
                status = "not_testable_too_few_samples"
            elif n_unique < 2:
                status = "not_testable_no_variation"
            else:
                status = "available"
        rows.append(
            {
                "outcome": outcome,
                "n_observed": n_observed,
                "n_unique": n_unique,
                "missingness_rate": 1.0 - n_observed / len(data) if len(data) else 1.0,
                "status": status,
            }
        )
    return pd.DataFrame(rows)


def load_module8_inputs(module8_report: Path) -> dict[str, str]:
    with open(module8_report, encoding="utf-8") as handle:
        report = json.load(handle)
    return report.get("inputs", {})


def load_combined_signature_genes(metadata_dir: Path) -> pd.DataFrame:
    target = read_tsv_or_empty(metadata_dir / "module8_tf_target_signature_genes.tsv")
    pathway = read_tsv_or_empty(metadata_dir / "module8_pathway_signature_genes.tsv")
    frames = []
    for df in [target, pathway]:
        if not df.empty and {"axis", "gene"}.issubset(df.columns):
            keep = df[["axis", "gene"]].copy()
            keep["tf"] = df["tf"].values if "tf" in df.columns else ""
            frames.append(keep)
    return pd.concat(frames, ignore_index=True).drop_duplicates() if frames else pd.DataFrame(columns=["axis", "tf", "gene"])


def build_real_bulk_scores(metadata_dir: Path, module8_inputs: dict[str, str], gtf_path: Path) -> pd.DataFrame:
    signatures = load_combined_signature_genes(metadata_dir)
    frames = []
    tcga_expression = Path(module8_inputs.get("tcga_expression", ""))
    if tcga_expression.exists() and not signatures.empty:
        ensembl_map = module8.parse_gtf_ensembl_to_symbol(gtf_path) if gtf_path.exists() else {}
        tcga_expr = module8.load_tcga_expression_signature_genes(tcga_expression, signatures, ensembl_map)
        tcga_scores, _ = module8.compute_bulk_signature_scores(tcga_expr, signatures, "TCGA-LIHC", sample_type_mode="tcga")
        frames.append(build_bulk_axis_scores(tcga_scores))
    clinical_root = Path(module8_inputs.get("clinical_root", ""))
    icgc_expression = clinical_root / "ICGCsymbol.txt"
    if icgc_expression.exists() and not signatures.empty:
        icgc_expr = module8.load_symbol_expression_signature_genes(icgc_expression, signatures)
        icgc_scores, _ = module8.compute_bulk_signature_scores(icgc_expr, signatures, "ICGC-LIRI-JP", sample_type_mode="icgc")
        frames.append(build_bulk_axis_scores(icgc_scores))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def prepare_tcga_outcomes(module8_inputs: dict[str, str]) -> pd.DataFrame:
    clinical_root = Path(module8_inputs.get("clinical_root", ""))
    complete = clinical_root / "clinical.xls"
    fallback = clinical_root / "tcgaClinical.txt"
    clinical = pd.DataFrame()
    if complete.exists():
        clinical = module8.prepare_tcga_clinical_covariates(module8.load_tcga_clinical_table(complete))
    if clinical.empty and fallback.exists():
        clinical = pd.read_csv(fallback, sep="\t")
    if clinical.empty:
        return pd.DataFrame()
    out = clinical.copy()
    if "Id" not in out.columns and "sample" in out.columns:
        out = out.rename(columns={"sample": "Id"})
    rename = {"fustat": "OS_event", "futime": "OS_time"}
    out = out.rename(columns=rename)
    keep = [col for col in ["Id", "age", "gender", "grade", "stage", "OS_event", "OS_time"] if col in out.columns]
    out = out[keep].copy()
    out["dataset_id"] = "TCGA-LIHC"
    out = out.rename(columns={"Id": "patient_id"})
    return out


def prepare_icgc_outcomes(module8_inputs: dict[str, str]) -> pd.DataFrame:
    clinical_root = Path(module8_inputs.get("clinical_root", ""))
    clinical_path = clinical_root / "icgcClinical.txt"
    survival_path = clinical_root / "ICGCtime.txt"
    if not clinical_path.exists():
        return pd.DataFrame()
    clinical = pd.read_csv(clinical_path, sep="\t").rename(
        columns={"Id": "patient_id", "Age": "age", "Gender": "gender", "Stage": "stage"}
    )
    if survival_path.exists():
        survival = pd.read_csv(survival_path, sep="\t").rename(
            columns={"id": "patient_id", "fustat": "OS_event", "futime": "OS_time"}
        )
        clinical = clinical.merge(survival, on="patient_id", how="left")
    clinical["dataset_id"] = "ICGC-LIRI-JP"
    keep = [col for col in ["dataset_id", "patient_id", "age", "gender", "stage", "OS_event", "OS_time"] if col in clinical.columns]
    return clinical[keep].copy()


def merge_bulk_outcomes(axis_scores: pd.DataFrame, module8_inputs: dict[str, str]) -> pd.DataFrame:
    if axis_scores.empty:
        return axis_scores.copy()
    scores = axis_scores.copy()
    scores["patient_id"] = scores["sample"].map(module8.sample_to_patient_id)
    outcomes = pd.concat([prepare_tcga_outcomes(module8_inputs), prepare_icgc_outcomes(module8_inputs)], ignore_index=True)
    if outcomes.empty:
        return scores
    merged = scores.merge(outcomes, on=["dataset_id", "patient_id"], how="left")
    for col in ["age", "gender", "stage", "grade", "OS_event", "OS_time"]:
        if col in merged.columns:
            merged[col] = safe_numeric(merged[col])
    return merged


def select_covariates(data: pd.DataFrame, outcome: str, outcome_type: str) -> list[str]:
    candidates = ["age", "gender"]
    if outcome_type == "cox" and outcome != "stage":
        candidates.append("stage")
    if outcome == "grade":
        candidates.append("stage")
    covariates = []
    for col in candidates:
        if col in data.columns and safe_numeric(data[col]).notna().sum() >= max(8, len(data) // 5) and safe_numeric(data[col]).nunique() > 1:
            covariates.append(col)
    return covariates


def run_path_analyses(
    bulk: pd.DataFrame,
    scrna: pd.DataFrame,
    n_bootstrap: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    coef_frames = []
    indirect_frames = []
    availability_frames = []

    analyses: list[tuple[str, pd.DataFrame, str, str]] = []
    for dataset_id, group in bulk.groupby("dataset_id", dropna=False):
        analyses.extend(
            [
                (str(dataset_id), group, "stage", "continuous"),
                (str(dataset_id), group, "grade", "continuous"),
                (str(dataset_id), group.rename(columns={"OS_time": "OS_time_for_model"}), "OS_time_for_model", "cox"),
            ]
        )
    if not bulk.empty:
        pooled = bulk.copy()
        pooled["dataset_code"] = (pooled["dataset_id"].astype("category").cat.codes).astype(float)
        analyses.extend(
            [
                ("bulk_pooled", pooled, "stage", "continuous"),
                ("bulk_pooled", pooled, "grade", "continuous"),
                ("bulk_pooled", pooled.rename(columns={"OS_time": "OS_time_for_model"}), "OS_time_for_model", "cox"),
            ]
        )
    if not scrna.empty:
        for (run_id, method), group in scrna.groupby(["run_id", "method"], dropna=False):
            dataset_id = f"scRNA_pseudobulk|{run_id}|{method}"
            analyses.extend(
                [
                    (dataset_id, group.rename(columns={"malignant_fraction": "outcome_malignant_fraction"}), "outcome_malignant_fraction", "continuous"),
                    (dataset_id, group.rename(columns={"mean_C_malignant_like_fate": "outcome_mean_malignant_fate"}), "outcome_mean_malignant_fate", "continuous"),
                    (dataset_id, group.rename(columns={"sample_stage_ordinal": "outcome_sample_stage"}), "outcome_sample_stage", "continuous"),
                ]
            )

    for idx, (dataset_id, data, outcome, outcome_type) in enumerate(analyses):
        outcome_check_col = "OS_event" if outcome_type == "cox" else outcome
        availability = assess_outcome_availability(data, [outcome_check_col], min_n=8)
        availability["dataset_id"] = dataset_id
        availability["outcome_model_col"] = outcome
        availability["outcome_type"] = outcome_type
        availability_frames.append(availability)
        if availability.loc[0, "status"] != "available":
            continue
        model_data = data.copy()
        if outcome_type == "cox":
            if "OS_event" not in model_data.columns or outcome not in model_data.columns:
                continue
            model_data = model_data.rename(columns={outcome: "time", "OS_event": "event"})
            covariates = select_covariates(model_data, "time", "cox")
            coefs = compute_path_coefficients(model_data, "time", "cox", covariates, dataset_id)
            effects = bootstrap_indirect_effects(
                model_data, "time", "cox", covariates, n_bootstrap=max(100, min(n_bootstrap, 250)), random_state=seed + idx, dataset_id=dataset_id
            )
        else:
            covariates = select_covariates(model_data, outcome, "continuous")
            coefs = compute_path_coefficients(model_data, outcome, "continuous", covariates, dataset_id)
            effects = bootstrap_indirect_effects(
                model_data, outcome, "continuous", covariates, n_bootstrap=n_bootstrap, random_state=seed + idx, dataset_id=dataset_id
            )
        coef_frames.append(coefs)
        indirect_frames.append(effects)

    coefficients = pd.concat(coef_frames, ignore_index=True) if coef_frames else pd.DataFrame()
    indirect = pd.concat(indirect_frames, ignore_index=True) if indirect_frames else pd.DataFrame()
    if not indirect.empty:
        indirect["p.adjust.global"] = benjamini_hochberg(indirect["pvalue"])
    availability = pd.concat(availability_frames, ignore_index=True) if availability_frames else pd.DataFrame()
    return coefficients, indirect, availability


def build_evidence_grade(indirect: pd.DataFrame, coefficients: pd.DataFrame, availability: pd.DataFrame) -> pd.DataFrame:
    if indirect.empty:
        label = "not_testable"
        n_supported = 0
        n_partial = 0
        n_tested = 0
    else:
        tested = indirect.loc[indirect["status"].eq("tested") & indirect["indirect_path"].eq("A_to_B_to_C_to_outcome")].copy()
        n_tested = int(len(tested))
        supported_mask = (tested["effect"] > 0) & (tested["ci_low"] > 0) & (tested["p.adjust.global"] <= 0.10)
        partial_mask = (tested["effect"] > 0) & (tested["ci_low"] > 0)
        n_supported = int(supported_mask.sum())
        n_partial = int((partial_mask & ~supported_mask).sum())
        has_bulk_pooled = bool((tested["dataset_id"].eq("bulk_pooled") & supported_mask).any())
        has_independent_bulk = bool((tested["dataset_id"].isin(["TCGA-LIHC", "ICGC-LIRI-JP"]) & supported_mask).any())
        if has_bulk_pooled and has_independent_bulk:
            label = "mediation_supported"
        elif n_supported > 0 or n_partial > 0:
            label = "partial_mediation_support"
        else:
            label = "mediation_not_supported"
    return pd.DataFrame(
        [
            {
                "evidence_domain": "Module9.3_statistical_mediation",
                "final_support_label": label,
                "n_tested_sequential_indirect_effects": n_tested,
                "n_supported_sequential_indirect_effects": n_supported,
                "n_partial_sequential_indirect_effects": n_partial,
                "n_available_outcomes": int(availability["status"].eq("available").sum()) if not availability.empty else 0,
                "n_path_coefficients": int(len(coefficients)),
            }
        ]
    )


def write_conclusions(path: Path, grade: pd.DataFrame, indirect: pd.DataFrame) -> None:
    row = grade.iloc[0].to_dict() if not grade.empty else {}
    lines = [
        "# Module 9.3 Statistical Mediation Evidence",
        "",
        f"- Evidence label: `{row.get('final_support_label', 'not_testable')}`.",
        f"- Tested sequential indirect effects: {row.get('n_tested_sequential_indirect_effects', 0)}.",
        f"- Supported sequential indirect effects: {row.get('n_supported_sequential_indirect_effects', 0)}.",
        "",
        "## Primary sequential path",
        "",
    ]
    if indirect.empty:
        lines.append("- No indirect effects were testable.")
    else:
        seq = indirect.loc[indirect["indirect_path"].eq("A_to_B_to_C_to_outcome")].copy()
        seq = seq.sort_values(["p.adjust.global", "dataset_id"], na_position="last")
        for _, r in seq.head(12).iterrows():
            lines.append(
                "- `{dataset}` / `{outcome}`: effect={effect:.3g}, 95% CI [{lo:.3g}, {hi:.3g}], q={q:.3g}, status=`{status}`.".format(
                    dataset=r["dataset_id"],
                    outcome=r["outcome"],
                    effect=float(r["effect"]) if pd.notna(r["effect"]) else math.nan,
                    lo=float(r["ci_low"]) if pd.notna(r["ci_low"]) else math.nan,
                    hi=float(r["ci_high"]) if pd.notna(r["ci_high"]) else math.nan,
                    q=float(r["p.adjust.global"]) if pd.notna(r["p.adjust.global"]) else math.nan,
                    status=r["status"],
                )
            )
    lines.extend(
        [
            "",
            "## Interpretation guardrails",
            "",
            "- SEM is implemented as observed-variable path analysis with bootstrap indirect effects.",
            "- Bulk A axis is the inverse of the tier1 rescue signature, approximating HNF4A/PPARA loss.",
            "- scRNA pseudo-bulk evaluates malignant burden/stage-like outcomes, not OS.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_table(df: pd.DataFrame, path: Path, compress: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, compression="gzip" if compress else None)


def run_module(args: argparse.Namespace) -> dict[str, object]:
    start = time.time()
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    inputs = load_module8_inputs(args.module8_report)

    bulk_axis = build_real_bulk_scores(args.metadata_dir, inputs, args.gtf_path)
    bulk_axis = merge_bulk_outcomes(bulk_axis, inputs)
    cell_scores = read_tsv_or_empty(args.module9_1_cell_scores)
    stage_by_sample = read_tsv_or_empty(args.stage_by_sample)
    scrna_pseudo = build_scrna_pseudobulk_scores(cell_scores, stage_by_sample)
    coefficients, indirect, availability = run_path_analyses(bulk_axis, scrna_pseudo, args.n_bootstrap, args.seed)
    grade = build_evidence_grade(indirect, coefficients, availability)

    outputs = {
        "bulk_sample_axis_scores": args.metadata_dir / "module9_3_bulk_sample_axis_scores.tsv.gz",
        "scrna_pseudobulk_axis_scores": args.metadata_dir / "module9_3_scrna_pseudobulk_axis_scores.tsv.gz",
        "path_model_coefficients": args.metadata_dir / "module9_3_path_model_coefficients.tsv",
        "indirect_effects": args.metadata_dir / "module9_3_indirect_effects.tsv",
        "outcome_availability": args.metadata_dir / "module9_3_outcome_availability.tsv",
        "evidence_grade": args.metadata_dir / "module9_3_evidence_grade.tsv",
        "main_conclusions": args.metadata_dir / "module9_3_main_conclusions.md",
        "report": args.metadata_dir / "module9_3_report.json",
    }
    write_table(bulk_axis, outputs["bulk_sample_axis_scores"], compress=True)
    write_table(scrna_pseudo, outputs["scrna_pseudobulk_axis_scores"], compress=True)
    write_table(coefficients, outputs["path_model_coefficients"])
    write_table(indirect, outputs["indirect_effects"])
    write_table(availability, outputs["outcome_availability"])
    write_table(grade, outputs["evidence_grade"])
    write_conclusions(outputs["main_conclusions"], grade, indirect)

    report = {
        "module": "9.3_statistical_mediation_evidence",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "runtime_seconds": round(time.time() - start, 3),
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "package_versions": {
            "numpy": package_version("numpy"),
            "pandas": package_version("pandas"),
            "scipy": package_version("scipy"),
            "statsmodels": package_version("statsmodels"),
            "matplotlib": package_version("matplotlib"),
        },
        "inputs": {
            "module8_report": str(args.module8_report),
            "module9_1_cell_scores": str(args.module9_1_cell_scores),
            "stage_by_sample": str(args.stage_by_sample),
            "gtf_path": str(args.gtf_path),
            "module8_inputs": inputs,
        },
        "parameters": {"n_bootstrap": args.n_bootstrap, "seed": args.seed},
        "counts": {
            "n_bulk_samples": int(len(bulk_axis)),
            "n_scrna_pseudobulk_samples": int(len(scrna_pseudo)),
            "n_path_coefficients": int(len(coefficients)),
            "n_indirect_effect_rows": int(len(indirect)),
            "n_outcome_availability_rows": int(len(availability)),
        },
        "evidence_grade": grade.iloc[0].to_dict() if not grade.empty else {},
        "outputs": {key: str(value) for key, value in outputs.items()},
    }
    outputs["report"].write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    return report


def main() -> None:
    args = parse_args()
    report = run_module(args)
    print(json.dumps({"module": report["module"], "evidence_grade": report["evidence_grade"], "outputs": report["outputs"]}, indent=2))


if __name__ == "__main__":
    main()
