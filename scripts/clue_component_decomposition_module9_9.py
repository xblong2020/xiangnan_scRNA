from __future__ import annotations

import argparse
import io
import json
import re
import tarfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

import h5py
import numpy as np
import pandas as pd

try:
    from scripts.drugreflector_metadata_crossvalidation_module9_8 import (
        collapse_pert_metadata,
        normalized_name_series,
        read_pert_info,
    )
except ModuleNotFoundError:
    from drugreflector_metadata_crossvalidation_module9_8 import (
        collapse_pert_metadata,
        normalized_name_series,
        read_pert_info,
    )


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA_DIR = ROOT / "metadata/driver"
DEFAULT_ARCHIVE = DEFAULT_METADATA_DIR / "module9_8_clue_results.tar.gz"
DEFAULT_PRIORITY = DEFAULT_METADATA_DIR / "module9_9_landmark_decomposition_priority.tsv"
DEFAULT_PHASE1 = DEFAULT_METADATA_DIR / "GSE92742_Broad_LINCS_pert_info.txt.gz"
DEFAULT_PHASE2 = DEFAULT_METADATA_DIR / "GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz"
DEFAULT_L1000 = DEFAULT_METADATA_DIR / "module9_8_drugreflector_metadata_crossvalidation_l1000fwd_mapped.tsv"
OUTPUT_STEM = "module9_9_landmark_decomposition"
LIVER_CELL_LINES = ["HEPG2", "HA1E", "HCC515"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 9.9 CLUE component decomposition.")
    parser.add_argument("--clue-archive", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument("--drugreflector-priority", type=Path, default=DEFAULT_PRIORITY)
    parser.add_argument("--phase1-pert-info", type=Path, default=DEFAULT_PHASE1)
    parser.add_argument("--phase2-pert-info", type=Path, default=DEFAULT_PHASE2)
    parser.add_argument("--l1000-mapped", type=Path, default=DEFAULT_L1000)
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA_DIR)
    parser.add_argument("--seed", type=int, default=20260616)
    return parser.parse_args()


def extract_gctx(archive_path: Path, suffix: str) -> tuple[np.ndarray, list[str]]:
    with tarfile.open(archive_path, "r:gz") as archive:
        members = [
            member
            for member in archive.getmembers()
            if member.isfile() and member.name.endswith(f"/{suffix}")
        ]
        if len(members) != 1:
            raise ValueError(f"expected one {suffix}, found {len(members)}")
        handle = archive.extractfile(members[0])
        if handle is None:
            raise ValueError(f"could not extract {suffix}")
        data = handle.read()
    with h5py.File(io.BytesIO(data), "r") as h5:
        matrix = h5["0/DATA/0/matrix"][:].reshape(-1)
        ids = [value.decode("utf-8") for value in h5["0/META/ROW/id"][:]]
    return matrix, ids


def parse_signature_ids(signature_ids: Sequence[str]) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    pattern = re.compile(
        r"^[^_]+_(?P<cell>[^_]+)_[^:]+:(?P<compound>BRD-[A-Z]\d+)"
    )
    for sig_id in signature_ids:
        match = pattern.search(str(sig_id))
        rows.append(
            {
                "sig_id": str(sig_id),
                "compound": match.group("compound") if match else "",
                "cell_line": match.group("cell") if match else "",
            }
        )
    return pd.DataFrame(rows)


def load_signature_components(archive_path: Path) -> pd.DataFrame:
    up, ids = extract_gctx(archive_path, "cs_up_n1x476251.gctx")
    down, down_ids = extract_gctx(archive_path, "cs_dn_n1x476251.gctx")
    combined, combined_ids = extract_gctx(archive_path, "cs_n1x476251.gctx")
    if ids != down_ids or ids != combined_ids:
        raise ValueError("CLUE component matrices have different signature ordering")
    frame = parse_signature_ids(ids)
    frame["cs_up"] = up
    frame["cs_down"] = down
    frame["cs_combined"] = combined
    return frame.loc[frame["compound"].ne("")].reset_index(drop=True)


def aggregate_clue_components(signatures: pd.DataFrame) -> pd.DataFrame:
    required = {"compound", "cell_line", "cs_up", "cs_down", "cs_combined"}
    missing = required.difference(signatures.columns)
    if missing:
        raise ValueError(f"signature components missing columns: {sorted(missing)}")
    frame = signatures.copy()
    frame["malignant_suppression"] = -pd.to_numeric(
        frame["cs_down"], errors="coerce"
    )
    frame["rescue_activation"] = pd.to_numeric(frame["cs_up"], errors="coerce")
    frame["combined_connectivity"] = pd.to_numeric(
        frame["cs_combined"], errors="coerce"
    )
    overall = (
        frame.groupby("compound")[
            ["rescue_activation", "malignant_suppression", "combined_connectivity"]
        ]
        .mean()
        .rename(
            columns={
                "rescue_activation": "clue_rescue_component",
                "malignant_suppression": "clue_malignant_suppression_component",
                "combined_connectivity": "clue_combined_component",
            }
        )
    )
    counts = frame.groupby("compound").size().rename("clue_signature_count")
    output = overall.join(counts)
    for cell_line in LIVER_CELL_LINES:
        subset = frame.loc[frame["cell_line"].eq(cell_line)]
        cell = (
            subset.groupby("compound")[
                ["rescue_activation", "malignant_suppression", "combined_connectivity"]
            ]
            .mean()
            .rename(
                columns={
                    "rescue_activation": f"clue_{cell_line.lower()}_rescue_component",
                    "malignant_suppression": f"clue_{cell_line.lower()}_malignant_suppression_component",
                    "combined_connectivity": f"clue_{cell_line.lower()}_combined_component",
                }
            )
        )
        output = output.join(cell, how="left")
    return output.reset_index()


def percentile(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").rank(pct=True, method="average")


def add_component_percentiles(aggregate: pd.DataFrame) -> pd.DataFrame:
    output = aggregate.copy()
    output["clue_rescue_percentile"] = percentile(output["clue_rescue_component"])
    output["clue_malignant_suppression_percentile"] = percentile(
        output["clue_malignant_suppression_component"]
    )
    output["clue_combined_percentile"] = percentile(output["clue_combined_component"])
    output["clue_branch_balance_raw"] = output[
        ["clue_rescue_component", "clue_malignant_suppression_component"]
    ].min(axis=1)
    output["clue_branch_balance_percentile"] = percentile(
        output["clue_branch_balance_raw"]
    )

    liver_balance_columns: list[str] = []
    for cell_line in LIVER_CELL_LINES:
        prefix = f"clue_{cell_line.lower()}"
        rescue_column = f"{prefix}_rescue_component"
        malignant_column = f"{prefix}_malignant_suppression_component"
        if rescue_column not in output or malignant_column not in output:
            continue
        rescue_pct = f"{prefix}_rescue_percentile"
        malignant_pct = f"{prefix}_malignant_suppression_percentile"
        balance_pct = f"{prefix}_branch_balance_percentile"
        output[rescue_pct] = percentile(output[rescue_column])
        output[malignant_pct] = percentile(output[malignant_column])
        balance_raw = output[[rescue_column, malignant_column]].min(axis=1)
        output[balance_pct] = percentile(balance_raw)
        liver_balance_columns.append(balance_pct)
    output["clue_liver_context_percentile"] = (
        output[liver_balance_columns].mean(axis=1, skipna=True)
        if liver_balance_columns
        else np.nan
    )
    return output


def map_components_to_compounds(
    compounds: pd.DataFrame,
    clue: pd.DataFrame,
    metadata: pd.DataFrame,
) -> pd.DataFrame:
    drug_meta = compounds[["compound"]].merge(
        metadata, on="compound", how="left", validate="one_to_one"
    )
    clue_meta = clue.merge(metadata, on="compound", how="left", validate="one_to_one")
    drug_meta["normalized_name"] = normalized_name_series(drug_meta["pert_iname"])
    clue_meta["normalized_name"] = normalized_name_series(clue_meta["pert_iname"])

    exact = {str(row["compound"]): row for _, row in clue_meta.iterrows()}
    by_inchi = {
        str(key): group
        for key, group in clue_meta.dropna(subset=["inchi_key"]).groupby("inchi_key")
    }
    by_name = {
        str(key): group
        for key, group in clue_meta.loc[clue_meta["normalized_name"].ne("")].groupby(
            "normalized_name"
        )
    }
    score_columns = [
        column
        for column in clue.columns
        if column != "compound"
    ]
    rows: list[dict[str, object]] = []
    for _, row in drug_meta.iterrows():
        compound = str(row["compound"])
        if compound in exact:
            match_type = "exact_id"
            matches = pd.DataFrame([exact[compound]])
        elif pd.notna(row.get("inchi_key")) and str(row["inchi_key"]) in by_inchi:
            match_type = "inchi_key"
            matches = by_inchi[str(row["inchi_key"])]
        elif row["normalized_name"] and row["normalized_name"] in by_name:
            match_type = "normalized_name"
            matches = by_name[row["normalized_name"]]
        else:
            match_type = "unmapped"
            matches = clue_meta.iloc[0:0]
        output: dict[str, object] = {
            "compound": compound,
            "clue_component_match_type": match_type,
            "clue_component_matched_ids": ",".join(
                sorted(set(matches["compound"].astype(str)))
            ),
        }
        for column in score_columns:
            values = pd.to_numeric(matches[column], errors="coerce").dropna()
            output[column] = float(values.mean()) if not values.empty else np.nan
        rows.append(output)
    return pd.DataFrame(rows)


def enrich_compound_metadata(
    compounds: pd.DataFrame,
    metadata: pd.DataFrame,
) -> pd.DataFrame:
    metadata_columns = [
        column
        for column in [
            "compound",
            "pert_iname",
            "pert_type",
            "is_touchstone",
            "inchi_key",
            "canonical_smiles",
            "pubchem_cid",
        ]
        if column in metadata.columns
    ]
    output = compounds.merge(
        metadata[metadata_columns].drop_duplicates("compound"),
        on="compound",
        how="left",
        suffixes=("", "_metadata"),
        validate="one_to_one",
    )
    for column in metadata_columns:
        if column == "compound":
            continue
        metadata_column = f"{column}_metadata"
        if metadata_column not in output:
            continue
        if column in output:
            output[column] = output[column].where(
                output[column].notna(), output[metadata_column]
            )
            output = output.drop(columns=metadata_column)
        else:
            output = output.rename(columns={metadata_column: column})
    return output


def evidence_tier(row: pd.Series) -> str:
    both = bool(row.get("both_biological_branches_top_200", False))
    clue_balance = float(row.get("clue_branch_balance_percentile", 0.0) or 0.0)
    clue_combined = float(row.get("clue_combined_percentile", 0.0) or 0.0)
    if both and clue_balance >= 0.9 and clue_combined >= 0.75:
        return "A_cross_platform_balanced"
    if both and clue_balance >= 0.75:
        return "B_cross_platform_support"
    if both:
        return "C_drugreflector_balanced"
    if clue_balance >= 0.9:
        return "D_clue_balanced_only"
    return "E_exploratory"


def build_evidence_adjusted_priority(
    drugreflector: pd.DataFrame,
    clue: pd.DataFrame,
    l1000: pd.DataFrame | None = None,
) -> pd.DataFrame:
    output = drugreflector.merge(clue, on="compound", how="left", validate="one_to_one")
    output["drugreflector_decomposition_percentile"] = percentile(
        output["decomposition_score"]
    )
    if l1000 is not None and not l1000.empty:
        keep = [
            column
            for column in [
                "compound",
                "l1000_similar_best_rank",
                "l1000_opposite_best_rank",
                "l1000_support_score",
            ]
            if column in l1000.columns
        ]
        output = output.merge(
            l1000[keep].drop_duplicates("compound"),
            on="compound",
            how="left",
            suffixes=("", "_external"),
            validate="one_to_one",
        )
        for column in keep:
            if column == "compound":
                continue
            external_column = f"{column}_external"
            if external_column not in output:
                continue
            if column in output:
                output[column] = output[column].where(
                    output[column].notna(), output[external_column]
                )
                output = output.drop(columns=external_column)
            else:
                output = output.rename(columns={external_column: column})
    l1000_support = pd.to_numeric(
        output.get("l1000_support_score", pd.Series(0.0, index=output.index)),
        errors="coerce",
    ).fillna(0.0)
    l1000_positive = l1000_support.clip(lower=0.0, upper=1.0)
    clue_balance = pd.to_numeric(
        output.get(
            "clue_branch_balance_percentile", pd.Series(0.0, index=output.index)
        ),
        errors="coerce",
    ).fillna(0.0)
    clue_combined = pd.to_numeric(
        output.get("clue_combined_percentile", pd.Series(0.0, index=output.index)),
        errors="coerce",
    ).fillna(0.0)
    clue_liver = pd.to_numeric(
        output.get(
            "clue_liver_context_percentile", pd.Series(0.0, index=output.index)
        ),
        errors="coerce",
    ).fillna(0.0)
    output["evidence_adjusted_score"] = (
        0.55 * output["drugreflector_decomposition_percentile"]
        + 0.20 * clue_balance
        + 0.10 * clue_combined
        + 0.10 * clue_liver
        + 0.05 * l1000_positive
    )
    output["external_component_available"] = output.get(
        "clue_component_match_type", pd.Series("unmapped", index=output.index)
    ).ne("unmapped")
    output["evidence_tier"] = output.apply(evidence_tier, axis=1)
    tier_order = {
        "A_cross_platform_balanced": 0,
        "B_cross_platform_support": 1,
        "C_drugreflector_balanced": 2,
        "D_clue_balanced_only": 3,
        "E_exploratory": 4,
    }
    output["_tier_order"] = output["evidence_tier"].map(tier_order)
    output = output.sort_values(
        ["_tier_order", "evidence_adjusted_score", "decomposition_score", "compound"],
        ascending=[True, False, False, True],
    ).drop(columns="_tier_order")
    output["final_priority_rank"] = np.arange(1, len(output) + 1)
    return output.reset_index(drop=True)


def main() -> None:
    args = parse_args()
    np.random.seed(args.seed)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)

    signatures = load_signature_components(args.clue_archive)
    aggregate = add_component_percentiles(aggregate_clue_components(signatures))
    drugreflector = pd.read_csv(args.drugreflector_priority, sep="\t")
    phase1 = read_pert_info(args.phase1_pert_info, "GSE92742_phase1", 1)
    phase2 = read_pert_info(args.phase2_pert_info, "GSE70138_phase2", 2)
    metadata = collapse_pert_metadata([phase1, phase2])
    drugreflector = enrich_compound_metadata(drugreflector, metadata)
    mapped = map_components_to_compounds(drugreflector, aggregate, metadata)
    l1000 = (
        pd.read_csv(args.l1000_mapped, sep="\t")
        if args.l1000_mapped.is_file()
        else pd.DataFrame()
    )
    final = build_evidence_adjusted_priority(drugreflector, mapped, l1000)

    outputs = {
        "signature_components": args.metadata_dir / f"{OUTPUT_STEM}_clue_signature_components.tsv.gz",
        "compound_components": args.metadata_dir / f"{OUTPUT_STEM}_clue_compound_components.tsv",
        "component_crosswalk": args.metadata_dir / f"{OUTPUT_STEM}_clue_component_crosswalk.tsv",
        "final_priority": args.metadata_dir / f"{OUTPUT_STEM}_final_priority.tsv",
        "report": args.metadata_dir / f"{OUTPUT_STEM}_final_report.json",
    }
    signatures.to_csv(
        outputs["signature_components"], sep="\t", index=False, compression="gzip"
    )
    aggregate.to_csv(outputs["compound_components"], sep="\t", index=False)
    mapped.to_csv(outputs["component_crosswalk"], sep="\t", index=False)
    final.to_csv(outputs["final_priority"], sep="\t", index=False)

    report = {
        "module": "module9_9_clue_component_decomposition",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "completed",
        "inputs": {
            "clue_archive": str(args.clue_archive.resolve()),
            "drugreflector_priority": str(args.drugreflector_priority.resolve()),
            "seed": args.seed,
        },
        "outputs": {key: str(path.resolve()) for key, path in outputs.items()},
        "summary": {
            "n_clue_signatures": int(len(signatures)),
            "n_clue_compounds": int(len(aggregate)),
            "n_drugreflector_compounds": int(len(drugreflector)),
            "n_exact_component_matches": int(
                mapped["clue_component_match_type"].eq("exact_id").sum()
            ),
            "n_alias_component_matches": int(
                mapped["clue_component_match_type"]
                .isin(["inchi_key", "normalized_name"])
                .sum()
            ),
            "n_unmapped_components": int(
                mapped["clue_component_match_type"].eq("unmapped").sum()
            ),
            "evidence_tier_counts": {
                str(key): int(value)
                for key, value in final["evidence_tier"].value_counts().items()
            },
        },
        "scoring": {
            "drugreflector_decomposition_percentile": 0.55,
            "clue_branch_balance_percentile": 0.20,
            "clue_combined_percentile": 0.10,
            "clue_liver_context_percentile": 0.10,
            "positive_l1000_support": 0.05,
        },
        "interpretation_boundary": (
            "CLUE rescue and malignant components are decomposed from cs_up and -cs_down "
            "within the completed paired query. Evidence tiers prioritize cross-platform "
            "balance and do not imply clinical efficacy."
        ),
    }
    outputs["report"].write_text(
        json.dumps(report, indent=2, sort_keys=True), encoding="utf-8"
    )
    print(json.dumps(report["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
