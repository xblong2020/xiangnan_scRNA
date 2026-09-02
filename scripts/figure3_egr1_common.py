from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
TARGET_TF = "EGR1"
SEED = 15071990
PSEUDOTIME_COLUMN = "driver_main_strict__pseudotime_rank"
STRICT_COLUMN = "driver_main_strict__eligible"
STATE_COLUMN = "celloracle_state"

STATE_ORDER = [
    "normal_reference",
    "stressed_injured",
    "regenerative_progenitor",
    "proliferating_candidate",
    "malignant_or_malignant_like",
]
STATE_LABELS = {
    "normal_reference": "Normal/reference",
    "stressed_injured": "Stressed/injured",
    "regenerative_progenitor": "Regenerative/progenitor",
    "proliferating_candidate": "Proliferating candidate",
    "malignant_or_malignant_like": "Malignant/malignant-like",
}
STATE_PALETTE = {
    "normal_reference": "#B8B8B8",
    "stressed_injured": "#56B4E9",
    "regenerative_progenitor": "#009E73",
    "proliferating_candidate": "#E69F00",
    "malignant_or_malignant_like": "#D55E00",
}
DIVERGING_PALETTE = {"low": "#3B4CC0", "mid": "#F7F7F7", "high": "#B40426"}

SELECTION_WEIGHTS = {
    "celloracle_robustness": 0.20,
    "sctenifoldknk_robustness": 0.20,
    "transition_state_specificity": 0.20,
    "temporal_positioning": 0.15,
    "cross_dataset_stability": 0.15,
    "pathway_interpretability": 0.10,
}
PENALTY_WEIGHTS = {
    "generic_stress_penalty": 0.10,
    "proliferation_dependency_penalty": 0.10,
    "literature_overlap_penalty": 0.05,
}


def as_bool(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False)
    return values.astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y"})


def minmax_scale(values: pd.Series, neutral: float = 0.5) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce").astype(float)
    valid = numeric[np.isfinite(numeric)]
    if valid.empty:
        return pd.Series(np.full(len(numeric), neutral), index=numeric.index, dtype=float)
    lo, hi = float(valid.min()), float(valid.max())
    if math.isclose(lo, hi):
        scaled = pd.Series(np.full(len(numeric), neutral), index=numeric.index, dtype=float)
    else:
        scaled = (numeric - lo) / (hi - lo)
    return scaled.fillna(neutral).clip(0.0, 1.0)


def robust_inverse_scale(values: pd.Series) -> pd.Series:
    return 1.0 - minmax_scale(values)


def assign_fixed_stage(pseudotime: pd.Series) -> pd.Series:
    pt = pd.to_numeric(pseudotime, errors="coerce")
    stage = pd.Series(pd.NA, index=pt.index, dtype="object")
    stage.loc[(pt >= 0.0) & (pt < 0.33)] = "early"
    stage.loc[(pt >= 0.33) & (pt < 0.67)] = "intermediate"
    stage.loc[(pt >= 0.67) & (pt <= 1.0)] = "late"
    return stage


def choose_stress_transition_subset(audits: pd.DataFrame) -> tuple[str | None, str]:
    priority = [
        "stressed_injured",
        "stressed_regenerative",
        "intermediate_pseudotime",
        "malignant_like",
    ]
    required = {
        "subset",
        "n_cells",
        "n_datasets",
        "egr1_detection_rate",
        "max_dataset_fraction",
    }
    missing = required.difference(audits.columns)
    if missing:
        raise ValueError(f"Missing subset-audit columns: {sorted(missing)}")
    indexed = audits.set_index("subset", drop=False)
    for subset in priority:
        if subset not in indexed.index:
            continue
        row = indexed.loc[subset]
        eligible = (
            int(row["n_cells"]) >= 500
            and int(row["n_datasets"]) >= 3
            and float(row["egr1_detection_rate"]) > 0
            and float(row["max_dataset_fraction"]) < 0.80
        )
        if eligible:
            return subset, "first eligible subset under the prespecified priority"
    return None, "no subset met the prespecified publication-oriented eligibility thresholds"


def select_significant_genes(
    table: pd.DataFrame,
    target_tf: str = TARGET_TF,
    fdr_cutoff: float = 0.05,
    top_n: int = 20,
) -> tuple[pd.DataFrame, int]:
    required = {"tf", "gene", "distance", "p.adj"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"Missing perturbation columns: {sorted(missing)}")
    data = table.loc[
        table["tf"].astype(str).eq(target_tf)
        & ~table["gene"].astype(str).eq(target_tf)
    ].copy()
    data["distance"] = pd.to_numeric(data["distance"], errors="coerce")
    data["p.adj"] = pd.to_numeric(data["p.adj"], errors="coerce")
    significant = data.loc[
        np.isfinite(data["distance"])
        & np.isfinite(data["p.adj"])
        & data["p.adj"].lt(fdr_cutoff)
    ].sort_values(["distance", "p.adj", "gene"], ascending=[False, True, True])
    return significant.head(top_n).reset_index(drop=True), int(len(significant))


def compute_selection_score(frame: pd.DataFrame) -> pd.DataFrame:
    required = set(SELECTION_WEIGHTS).union(PENALTY_WEIGHTS)
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Missing selection-score columns: {sorted(missing)}")
    out = frame.copy()
    positive = sum(out[column].astype(float) * weight for column, weight in SELECTION_WEIGHTS.items())
    penalties = sum(out[column].astype(float) * weight for column, weight in PENALTY_WEIGHTS.items())
    out["positive_evidence_score"] = positive
    out["total_penalty"] = penalties
    out["selection_score"] = positive - penalties
    out = out.sort_values(["selection_score", "candidate"], ascending=[False, True]).reset_index(drop=True)
    out["selection_rank"] = np.arange(1, len(out) + 1)
    return out


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def file_fingerprint(path: Path) -> dict:
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
        "sha256": sha256_file(path),
    }


def fingerprint_paths(paths: Iterable[Path]) -> pd.DataFrame:
    rows = [file_fingerprint(path) for path in sorted(set(paths)) if path.is_file()]
    return pd.DataFrame(rows, columns=["path", "size_bytes", "mtime_ns", "sha256"])


def assert_figure3_path_isolated(path: Path) -> None:
    lower = str(path).replace("\\", "/").lower()
    if "sox4" in lower or "hnf4a" in lower:
        raise ValueError(f"Figure 3 output path crosses protected TF namespace: {path}")
    if not any(token in lower for token in ("figure3", "egr1", "three_axis_figure_consistency")):
        raise ValueError(f"Figure 3 output path lacks a dedicated namespace: {path}")


def write_json(value: Mapping | list, path: Path) -> None:
    assert_figure3_path_isolated(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")


def json_safe(value):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, np.ndarray):
        return [json_safe(x) for x in value.tolist()]
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(x) for x in value]
    if pd.isna(value):
        return None
    return value

