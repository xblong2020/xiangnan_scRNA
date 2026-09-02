from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA_DIR = ROOT / "metadata/driver"
DEFAULT_SIGNATURE = DEFAULT_METADATA_DIR / "module9_4_drug_reversal_signature.tsv"
DEFAULT_CHECKPOINT_DIR = DEFAULT_METADATA_DIR / "drugreflector_checkpoints"
DEFAULT_DRUGREFLECTOR_SOURCE = ROOT / "tmp/drugreflector-main-from-zip/drugreflector-main"
DEFAULT_MODULE98 = DEFAULT_METADATA_DIR / "module9_8_drugreflector_metadata_crossvalidation.tsv"
OUTPUT_STEM = "module9_9_landmark_decomposition"
MALIGNANT_COMPONENTS = {
    "sox4_state_specific",
    "ap1_stress_proliferation",
    "c_malignant_like_fate",
}
RESCUE_COMPONENTS = {
    "hnf4a_ppara_rescue",
    "mature_hepatocyte",
    "tier1_rescue",
}
DRUGREFLECTOR_PROFILES = [
    "malignant_only",
    "rescue_only",
    "combined_balanced",
    "combined_original_landmark",
    "primary_landmark_baseline",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Module 9.9 landmark-signature sensitivity decomposition."
    )
    parser.add_argument("--signature", type=Path, default=DEFAULT_SIGNATURE)
    parser.add_argument("--checkpoint-dir", type=Path, default=DEFAULT_CHECKPOINT_DIR)
    parser.add_argument("--drugreflector-source-dir", type=Path, default=DEFAULT_DRUGREFLECTOR_SOURCE)
    parser.add_argument("--module9-8-table", type=Path, default=DEFAULT_MODULE98)
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA_DIR)
    parser.add_argument("--top-n", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20260616)
    return parser.parse_args()


def package_version(name: str) -> str:
    try:
        return version(name)
    except PackageNotFoundError:
        return "not_installed"


def bool_series(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})


def select_landmark_signature(
    signature: pd.DataFrame,
    landmark_genes: set[str],
    sensitivity: bool = True,
) -> pd.DataFrame:
    required = {
        "gene",
        "desired_direction",
        "component",
        "final_weight",
        "include_primary",
        "include_sensitivity",
        "conflict_flag",
        "housekeeping_or_qc_flag",
    }
    missing = required.difference(signature.columns)
    if missing:
        raise ValueError(f"signature missing required columns: {sorted(missing)}")
    frame = signature.copy()
    frame["gene"] = frame["gene"].astype(str).str.strip().str.upper()
    frame["final_weight"] = pd.to_numeric(frame["final_weight"], errors="coerce")
    for column in [
        "include_primary",
        "include_sensitivity",
        "conflict_flag",
        "housekeeping_or_qc_flag",
    ]:
        frame[column] = bool_series(frame[column])
    include_column = "include_sensitivity" if sensitivity else "include_primary"
    frame = frame.loc[
        frame[include_column]
        & ~frame["conflict_flag"]
        & ~frame["housekeeping_or_qc_flag"]
        & frame["gene"].isin({str(gene).upper() for gene in landmark_genes})
        & frame["final_weight"].notna()
        & frame["desired_direction"].isin(["up", "down"])
    ].copy()
    return frame.sort_values(
        ["desired_direction", "component", "final_weight", "gene"],
        ascending=[True, True, False, True],
    ).reset_index(drop=True)


def signed_weights(frame: pd.DataFrame) -> pd.Series:
    signs = frame["desired_direction"].map({"up": 1.0, "down": -1.0})
    return signs * pd.to_numeric(frame["final_weight"], errors="coerce").fillna(0.0)


def normalize_absolute_mass(values: pd.Series) -> pd.Series:
    mass = float(values.abs().sum())
    if mass <= 0:
        return values.copy()
    return values / mass


def profile_rows(
    frame: pd.DataFrame,
    profile: str,
    normalize_mass: bool,
) -> pd.DataFrame:
    out = frame.copy()
    out["v_score"] = signed_weights(out)
    if normalize_mass:
        out["v_score"] = normalize_absolute_mass(out["v_score"])
    out["profile"] = profile
    return out


def build_decomposed_profiles(selected: pd.DataFrame) -> pd.DataFrame:
    malignant = selected.loc[selected["component"].isin(MALIGNANT_COMPONENTS)].copy()
    rescue = selected.loc[selected["component"].isin(RESCUE_COMPONENTS)].copy()
    combined = selected.copy()
    combined["v_score"] = signed_weights(combined)
    combined["v_score"] = combined.groupby("desired_direction")["v_score"].transform(
        normalize_absolute_mass
    )
    combined["profile"] = "combined_balanced"

    original = selected.copy()
    original["v_score"] = signed_weights(original)
    original["profile"] = "combined_original_landmark"

    primary = selected.loc[selected["include_primary"]].copy()
    primary["v_score"] = signed_weights(primary)
    primary["profile"] = "primary_landmark_baseline"

    profiles = pd.concat(
        [
            profile_rows(malignant, "malignant_only", normalize_mass=True),
            profile_rows(rescue, "rescue_only", normalize_mass=True),
            combined,
            original,
            primary,
        ],
        ignore_index=True,
        sort=False,
    )
    return profiles.sort_values(["profile", "v_score", "gene"], ascending=[True, False, True]).reset_index(drop=True)


def profiles_to_wide(profiles: pd.DataFrame) -> pd.DataFrame:
    wide = profiles.pivot_table(
        index="profile",
        columns="gene",
        values="v_score",
        aggfunc="sum",
        fill_value=0.0,
    )
    return wide.reindex(DRUGREFLECTOR_PROFILES).fillna(0.0)


def flatten_predictions(predictions: pd.DataFrame) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    profiles = list(dict.fromkeys(column[1] for column in predictions.columns))
    for profile in profiles:
        frame = pd.DataFrame(
            {
                "compound": predictions.index.astype(str),
                "profile": profile,
                "rank_0based": pd.to_numeric(predictions[("rank", profile)], errors="coerce"),
                "logit": pd.to_numeric(predictions[("logit", profile)], errors="coerce"),
                "prob": pd.to_numeric(predictions[("prob", profile)], errors="coerce"),
            }
        )
        frame["rank_1based"] = frame["rank_0based"] + 1
        rows.append(frame)
    return pd.concat(rows, ignore_index=True).sort_values(
        ["profile", "rank_0based", "compound"]
    ).reset_index(drop=True)


def rank_score(rank: pd.Series, n_compounds: int) -> pd.Series:
    denominator = max(n_compounds - 1, 1)
    return 1.0 - (pd.to_numeric(rank, errors="coerce") - 1.0) / denominator


def build_final_priority_table(
    predictions: pd.DataFrame,
    n_compounds: int,
) -> pd.DataFrame:
    ranks = predictions.pivot_table(
        index="compound",
        columns="profile",
        values="rank_1based",
        aggfunc="min",
    )
    for profile in DRUGREFLECTOR_PROFILES:
        if profile not in ranks.columns:
            ranks[profile] = np.nan
    output = ranks.reset_index()
    for profile in DRUGREFLECTOR_PROFILES:
        output[f"{profile}_rank_score"] = rank_score(output[profile], n_compounds)

    malignant_score = output["malignant_only_rank_score"]
    rescue_score = output["rescue_only_rank_score"]
    combined_score = output["combined_balanced_rank_score"]
    branch_min = pd.concat([malignant_score, rescue_score], axis=1).min(axis=1)
    branch_mean = pd.concat([malignant_score, rescue_score], axis=1).mean(axis=1)
    output["branch_balance_score"] = branch_min
    output["branch_rank_gap"] = (
        pd.to_numeric(output["malignant_only"], errors="coerce")
        - pd.to_numeric(output["rescue_only"], errors="coerce")
    ).abs()
    output["decomposition_score"] = (
        0.5 * branch_min + 0.3 * combined_score + 0.2 * branch_mean
    )
    rank_columns = DRUGREFLECTOR_PROFILES
    output["n_profiles_top_50"] = output[rank_columns].le(50).sum(axis=1)
    output["n_profiles_top_200"] = output[rank_columns].le(200).sum(axis=1)
    output["n_biological_branches_top_200"] = output[
        ["malignant_only", "rescue_only"]
    ].le(200).sum(axis=1)
    output["both_biological_branches_top_200"] = output[
        "n_biological_branches_top_200"
    ].eq(2)
    return output.sort_values(
        [
            "both_biological_branches_top_200",
            "decomposition_score",
            "n_profiles_top_200",
            "combined_balanced",
            "compound",
        ],
        ascending=[False, False, False, True, True],
    ).reset_index(drop=True)


def profile_qc(profiles: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for profile, group in profiles.groupby("profile", sort=False):
        rows.append(
            {
                "profile": profile,
                "n_genes": int(group["gene"].nunique()),
                "n_up_genes": int(group["v_score"].gt(0).sum()),
                "n_down_genes": int(group["v_score"].lt(0).sum()),
                "up_absolute_mass": float(group.loc[group["v_score"].gt(0), "v_score"].abs().sum()),
                "down_absolute_mass": float(group.loc[group["v_score"].lt(0), "v_score"].abs().sum()),
                "max_absolute_score": float(group["v_score"].abs().max()),
            }
        )
    return pd.DataFrame(rows)


def configure_source(source_dir: Path) -> None:
    if source_dir.exists():
        sys.path.insert(0, str(source_dir.resolve()))


def checkpoint_paths(checkpoint_dir: Path) -> list[Path]:
    paths = [checkpoint_dir / f"model_fold_{idx}.pt" for idx in range(3)]
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing checkpoints: {missing}")
    return paths


def main() -> None:
    args = parse_args()
    np.random.seed(args.seed)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    configure_source(args.drugreflector_source_dir)
    from drugreflector import DrugReflector

    model = DrugReflector(
        checkpoint_paths=[str(path) for path in checkpoint_paths(args.checkpoint_dir)]
    )
    landmark_sets = [
        set(str(gene).upper() for gene in genes)
        for genes in model.model.dimensions["var_names"]
    ]
    if not all(genes == landmark_sets[0] for genes in landmark_sets[1:]):
        raise ValueError("DrugReflector folds have different landmark gene sets")
    landmarks = landmark_sets[0]

    signature = pd.read_csv(args.signature, sep="\t")
    selected = select_landmark_signature(signature, landmarks, sensitivity=True)
    profiles = build_decomposed_profiles(selected)
    wide = profiles_to_wide(profiles)
    predictions = flatten_predictions(model.predict(wide))
    priority = build_final_priority_table(predictions, model.n_compounds)

    if args.module9_8_table.is_file():
        module98 = pd.read_csv(args.module9_8_table, sep="\t")
        external_columns = [
            column
            for column in [
                "compound",
                "pert_iname",
                "pert_type",
                "canonical_smiles",
                "pubchem_cid",
                "l1000_similar_best_rank",
                "l1000_opposite_best_rank",
                "clue_tau",
                "clue_hepg2_tau",
                "clue_ha1e_tau",
                "clue_hcc515_tau",
                "clue_match_type",
            ]
            if column in module98.columns
        ]
        priority = priority.merge(
            module98[external_columns].drop_duplicates("compound"),
            on="compound",
            how="left",
            validate="one_to_one",
        )

    outputs = {
        "signature_long": args.metadata_dir / f"{OUTPUT_STEM}_signature_long.tsv",
        "signature_wide": args.metadata_dir / f"{OUTPUT_STEM}_signature_wide.tsv",
        "profile_qc": args.metadata_dir / f"{OUTPUT_STEM}_profile_qc.tsv",
        "predictions_full": args.metadata_dir / f"{OUTPUT_STEM}_drugreflector_predictions.tsv.gz",
        "predictions_top": args.metadata_dir / f"{OUTPUT_STEM}_drugreflector_top.tsv",
        "priority": args.metadata_dir / f"{OUTPUT_STEM}_priority.tsv",
        "report": args.metadata_dir / f"{OUTPUT_STEM}_report.json",
    }
    profiles.to_csv(outputs["signature_long"], sep="\t", index=False)
    wide.to_csv(outputs["signature_wide"], sep="\t")
    qc = profile_qc(profiles)
    qc.to_csv(outputs["profile_qc"], sep="\t", index=False)
    predictions.to_csv(outputs["predictions_full"], sep="\t", index=False, compression="gzip")
    predictions.loc[predictions["rank_1based"].le(args.top_n)].to_csv(
        outputs["predictions_top"], sep="\t", index=False
    )
    priority.to_csv(outputs["priority"], sep="\t", index=False)

    report = {
        "module": "module9_9_landmark_signature_decomposition",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "completed_drugreflector_decomposition",
        "inputs": {
            "signature": str(args.signature.resolve()),
            "checkpoint_dir": str(args.checkpoint_dir.resolve()),
            "drugreflector_source_dir": str(args.drugreflector_source_dir.resolve()),
            "module9_8_table": str(args.module9_8_table.resolve()),
            "seed": args.seed,
        },
        "outputs": {key: str(path.resolve()) for key, path in outputs.items()},
        "summary": {
            "n_model_landmark_genes": len(landmarks),
            "n_selected_landmark_genes": int(selected["gene"].nunique()),
            "n_malignant_landmark_genes": int(
                selected["component"].isin(MALIGNANT_COMPONENTS).sum()
            ),
            "n_rescue_landmark_genes": int(
                selected["component"].isin(RESCUE_COMPONENTS).sum()
            ),
            "n_profiles": int(wide.shape[0]),
            "n_compounds": int(model.n_compounds),
            "n_compounds_both_branches_top_200": int(
                priority["both_biological_branches_top_200"].sum()
            ),
        },
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": package_version("numpy"),
            "pandas": package_version("pandas"),
            "torch": package_version("torch"),
        },
        "interpretation_boundary": (
            "Decomposition scores prioritize compounds that rank well in both malignant "
            "suppression and hepatocyte rescue profiles. They are comparative model scores, "
            "not evidence of clinical efficacy."
        ),
    }
    outputs["report"].write_text(
        json.dumps(report, indent=2, sort_keys=True), encoding="utf-8"
    )
    print(json.dumps(report["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
