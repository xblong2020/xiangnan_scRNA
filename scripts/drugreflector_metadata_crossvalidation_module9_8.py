from __future__ import annotations

import argparse
import io
import json
import platform
import tarfile
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

try:
    from scripts.drug_query_module9_5 import run_l1000fwd_query
except ModuleNotFoundError:
    from drug_query_module9_5 import run_l1000fwd_query


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA_DIR = ROOT / "metadata/driver"
DEFAULT_CONSENSUS = DEFAULT_METADATA_DIR / "module9_7_drugreflector_consensus_predictions.tsv"
DEFAULT_PHASE1_PERT_INFO = DEFAULT_METADATA_DIR / "GSE92742_Broad_LINCS_pert_info.txt.gz"
DEFAULT_PHASE2_PERT_INFO = DEFAULT_METADATA_DIR / "GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz"
DEFAULT_L1000_PAYLOAD = DEFAULT_METADATA_DIR / "module9_5_l1000fwd_query_payload.json"
DEFAULT_L1000_CANDIDATES = DEFAULT_METADATA_DIR / "module9_5_l1000fwd_candidate_ranking.tsv"
DEFAULT_CLUE_STATUS = DEFAULT_METADATA_DIR / "module9_5_clue_query_status.json"
DEFAULT_CLUE_RESULTS_ARCHIVE = DEFAULT_METADATA_DIR / "module9_8_clue_results.tar.gz"
DEFAULT_L1000_BASE_URL = "https://maayanlab.cloud/l1000fwd"
OUTPUT_STEM = "module9_8_drugreflector_metadata_crossvalidation"
MISSING_VALUES = {"", "-666", "-666.0", "nan", "none", "null"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Module 9.8 BRD metadata mapping and cross-validation.")
    parser.add_argument("--consensus", type=Path, default=DEFAULT_CONSENSUS)
    parser.add_argument("--phase1-pert-info", type=Path, default=DEFAULT_PHASE1_PERT_INFO)
    parser.add_argument("--phase2-pert-info", type=Path, default=DEFAULT_PHASE2_PERT_INFO)
    parser.add_argument("--l1000-payload", type=Path, default=DEFAULT_L1000_PAYLOAD)
    parser.add_argument("--l1000-candidates", type=Path, default=DEFAULT_L1000_CANDIDATES)
    parser.add_argument("--clue-status", type=Path, default=DEFAULT_CLUE_STATUS)
    parser.add_argument("--clue-results-archive", type=Path, default=DEFAULT_CLUE_RESULTS_ARCHIVE)
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA_DIR)
    parser.add_argument("--l1000fwd-base-url", default=DEFAULT_L1000_BASE_URL)
    parser.add_argument("--refresh-l1000fwd", action="store_true")
    parser.add_argument("--timeout-seconds", type=int, default=120)
    parser.add_argument("--seed", type=int, default=20260616)
    return parser.parse_args()


def package_version(name: str) -> str:
    try:
        return version(name)
    except PackageNotFoundError:
        return "not_installed"


def clean_value(value: object) -> object:
    if pd.isna(value):
        return pd.NA
    text = str(value).strip()
    if text.lower() in MISSING_VALUES:
        return pd.NA
    return text


def read_pert_info(path: Path, source_dataset: str, source_priority: int) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", compression="infer", low_memory=False)
    if "pert_id" not in frame.columns:
        raise ValueError(f"{path} missing pert_id")
    expected = [
        "pert_id",
        "pert_iname",
        "pert_type",
        "is_touchstone",
        "inchi_key_prefix",
        "inchi_key",
        "canonical_smiles",
        "pubchem_cid",
    ]
    for column in expected:
        if column not in frame.columns:
            frame[column] = pd.NA
        frame[column] = frame[column].map(clean_value)
    frame["source_dataset"] = source_dataset
    frame["source_priority"] = source_priority
    return frame[expected + ["source_dataset", "source_priority"]].copy()


def collapse_pert_metadata(frames: Sequence[pd.DataFrame]) -> pd.DataFrame:
    combined = pd.concat(frames, ignore_index=True, sort=False)
    combined = combined.loc[combined["pert_id"].notna()].copy()
    combined = combined.sort_values(["pert_id", "source_priority"], ascending=[True, False])
    value_columns = [
        "pert_iname",
        "pert_type",
        "is_touchstone",
        "inchi_key_prefix",
        "inchi_key",
        "canonical_smiles",
        "pubchem_cid",
    ]
    rows: list[dict[str, object]] = []
    for pert_id, group in combined.groupby("pert_id", sort=True):
        row: dict[str, object] = {
            "compound": str(pert_id),
            "metadata_source_datasets": ",".join(dict.fromkeys(group["source_dataset"].astype(str))),
            "metadata_source_row_count": int(len(group)),
        }
        conflict_columns: list[str] = []
        for column in value_columns:
            values = [value for value in group[column].tolist() if not pd.isna(value)]
            unique_values = list(dict.fromkeys(str(value) for value in values))
            row[column] = unique_values[0] if unique_values else pd.NA
            if len(unique_values) > 1:
                conflict_columns.append(column)
        row["metadata_conflict_flag"] = bool(conflict_columns)
        row["metadata_conflict_columns"] = ",".join(conflict_columns)
        rows.append(row)
    return pd.DataFrame(rows)


def aggregate_l1000fwd_candidates(candidates: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "compound",
        "l1000_similar_best_rank",
        "l1000_opposite_best_rank",
        "l1000_similar_signature_count",
        "l1000_opposite_signature_count",
        "l1000_similar_best_score",
        "l1000_opposite_best_score",
        "l1000_support_score",
        "l1000_cell_lines",
        "l1000_sig_ids",
    ]
    if candidates.empty or "compound_id" not in candidates.columns:
        return pd.DataFrame(columns=columns)
    frame = candidates.loc[candidates["compound_id"].notna()].copy()
    if frame.empty:
        return pd.DataFrame(columns=columns)
    frame["compound_id"] = frame["compound_id"].astype(str)
    rows: list[dict[str, object]] = []
    for compound, group in frame.groupby("compound_id", sort=True):
        similar = group.loc[group["result_group"].eq("similar")]
        opposite = group.loc[group["result_group"].eq("opposite")]
        similar_rank = pd.to_numeric(similar["rank_within_group"], errors="coerce").min()
        opposite_rank = pd.to_numeric(opposite["rank_within_group"], errors="coerce").min()
        similar_rr = 0.0 if pd.isna(similar_rank) else 1.0 / float(similar_rank)
        opposite_rr = 0.0 if pd.isna(opposite_rank) else 1.0 / float(opposite_rank)
        rows.append(
            {
                "compound": compound,
                "l1000_similar_best_rank": similar_rank,
                "l1000_opposite_best_rank": opposite_rank,
                "l1000_similar_signature_count": int(len(similar)),
                "l1000_opposite_signature_count": int(len(opposite)),
                "l1000_similar_best_score": pd.to_numeric(similar["raw_score"], errors="coerce").max(),
                "l1000_opposite_best_score": pd.to_numeric(opposite["raw_score"], errors="coerce").min(),
                "l1000_support_score": similar_rr - opposite_rr,
                "l1000_cell_lines": ",".join(sorted(set(group["cell_line"].dropna().astype(str)))),
                "l1000_sig_ids": ",".join(group["sig_id"].dropna().astype(str)),
            }
        )
    return pd.DataFrame(rows, columns=columns)


def classify_crossvalidation(row: pd.Series) -> str:
    similar = int(row.get("l1000_similar_signature_count", 0) or 0) > 0
    opposite = int(row.get("l1000_opposite_signature_count", 0) or 0) > 0
    if similar and opposite:
        return "mixed_l1000fwd_evidence"
    if similar:
        return "supported_by_l1000fwd_similar"
    if opposite:
        return "discordant_l1000fwd_opposite"
    return "drugreflector_only"


def read_gct_from_tar(archive_path: Path, filename: str) -> pd.DataFrame:
    with tarfile.open(archive_path, "r:gz") as archive:
        members = [
            member
            for member in archive.getmembers()
            if member.isfile() and member.name.endswith(f"/{filename}")
        ]
        if len(members) != 1:
            raise ValueError(f"expected one {filename} in {archive_path}, found {len(members)}")
        handle = archive.extractfile(members[0])
        if handle is None:
            raise ValueError(f"could not extract {filename} from {archive_path}")
        text = handle.read().decode("utf-8")
    frame = pd.read_csv(io.StringIO(text), sep="\t", skiprows=2)
    if frame.shape[1] != 2:
        raise ValueError(f"{filename} expected two columns, found {frame.shape[1]}")
    frame.columns = ["id", "clue_tau"]
    frame["clue_tau"] = pd.to_numeric(frame["clue_tau"], errors="coerce")
    return frame


def load_clue_results(archive_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary = read_gct_from_tar(archive_path, "pert_id_summary.gct").rename(
        columns={"id": "clue_compound"}
    )
    cell = read_gct_from_tar(archive_path, "pert_id_cell.gct").rename(
        columns={"id": "compound_cell"}
    )
    split = cell["compound_cell"].astype(str).str.rsplit(":", n=1, expand=True)
    cell["clue_compound"] = split[0]
    cell["cell_line"] = split[1]
    return summary, cell[["clue_compound", "cell_line", "clue_tau"]]


def build_clue_crosswalk(
    consensus_metadata: pd.DataFrame,
    clue_summary: pd.DataFrame,
    clue_cell: pd.DataFrame,
    all_metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    clue_meta = clue_summary.merge(
        all_metadata.rename(columns={"compound": "clue_compound"}),
        on="clue_compound",
        how="left",
        validate="one_to_one",
    )
    clue_meta["normalized_pert_name"] = normalized_name_series(clue_meta["pert_iname"])

    exact_groups = {
        str(compound): group
        for compound, group in clue_meta.groupby("clue_compound", sort=False)
    }
    inchi_groups = {
        str(inchi): group
        for inchi, group in clue_meta.dropna(subset=["inchi_key"]).groupby("inchi_key", sort=False)
    }
    name_groups = {
        str(name): group
        for name, group in clue_meta.loc[clue_meta["normalized_pert_name"].ne("")].groupby(
            "normalized_pert_name", sort=False
        )
    }
    cell_groups = {
        (str(compound), str(cell_line)): group["clue_tau"]
        for (compound, cell_line), group in clue_cell.groupby(
            ["clue_compound", "cell_line"], sort=False
        )
    }

    rows: list[dict[str, object]] = []
    for _, row in consensus_metadata.iterrows():
        compound = str(row["compound"])
        inchi_key = "" if pd.isna(row.get("inchi_key")) else str(row.get("inchi_key"))
        normalized_name = normalized_name_series(
            pd.Series([row.get("pert_iname", pd.NA)])
        ).iloc[0]
        if compound in exact_groups:
            match_type = "exact_id"
            matches = exact_groups[compound]
        elif inchi_key and inchi_key in inchi_groups:
            match_type = "inchi_key"
            matches = inchi_groups[inchi_key]
        elif normalized_name and normalized_name in name_groups:
            match_type = "normalized_name"
            matches = name_groups[normalized_name]
        else:
            match_type = "unmapped"
            matches = clue_meta.iloc[0:0]

        matched_ids = sorted(set(matches["clue_compound"].astype(str)))
        tau_values = pd.to_numeric(matches["clue_tau"], errors="coerce").dropna()
        clue_tau = float(tau_values.median()) if not tau_values.empty else np.nan
        output: dict[str, object] = {
            "compound": compound,
            "clue_match_type": match_type,
            "clue_matched_ids": ",".join(matched_ids),
            "clue_alias_count": len(matched_ids),
            "clue_tau": clue_tau,
            "clue_positive": bool(pd.notna(clue_tau) and clue_tau > 0),
            "clue_negative": bool(pd.notna(clue_tau) and clue_tau < 0),
            "clue_strong_support": bool(pd.notna(clue_tau) and clue_tau >= 90),
            "clue_strong_opposition": bool(pd.notna(clue_tau) and clue_tau <= -90),
        }
        for cell_line in ["HEPG2", "HA1E", "HCC515"]:
            values = [
                value
                for matched_id in matched_ids
                for value in pd.to_numeric(
                    cell_groups.get((matched_id, cell_line), pd.Series(dtype=float)),
                    errors="coerce",
                ).dropna()
            ]
            output[f"clue_{cell_line.lower()}_tau"] = (
                float(np.median(values)) if values else np.nan
            )
        rows.append(output)

    crosswalk = pd.DataFrame(rows)
    crosswalk["clue_rank_desc"] = crosswalk["clue_tau"].rank(
        method="min", ascending=False, na_option="bottom"
    )
    return crosswalk, clue_meta


def classify_integrated_support(row: pd.Series) -> str:
    l1000_similar = int(row.get("l1000_similar_signature_count", 0) or 0) > 0
    l1000_opposite = int(row.get("l1000_opposite_signature_count", 0) or 0) > 0
    clue_strong_positive = bool(row.get("clue_strong_support", False))
    clue_strong_negative = bool(row.get("clue_strong_opposition", False))
    clue_positive = bool(row.get("clue_positive", False))
    clue_negative = bool(row.get("clue_negative", False))
    if l1000_similar and clue_strong_positive:
        return "three_method_strong_support"
    if l1000_opposite or clue_strong_negative:
        return "independent_opposition_flag"
    if l1000_similar:
        return "drugreflector_l1000fwd_support"
    if clue_strong_positive:
        return "drugreflector_clue_strong_support"
    if clue_positive:
        return "drugreflector_clue_weak_positive"
    if clue_negative:
        return "drugreflector_clue_weak_negative"
    return "drugreflector_only_or_clue_neutral"


def build_crossvalidation_table(
    consensus: pd.DataFrame,
    metadata: pd.DataFrame,
    l1000_aggregate: pd.DataFrame,
    clue_status: str,
    clue_crosswalk: pd.DataFrame | None = None,
) -> pd.DataFrame:
    table = consensus.merge(metadata, on="compound", how="left", validate="one_to_one")
    table = table.merge(l1000_aggregate, on="compound", how="left", validate="one_to_one")
    if clue_crosswalk is not None:
        table = table.merge(clue_crosswalk, on="compound", how="left", validate="one_to_one")
    for column in ["l1000_similar_signature_count", "l1000_opposite_signature_count"]:
        table[column] = pd.to_numeric(table[column], errors="coerce").fillna(0).astype(int)
    for column in ["l1000_similar_best_rank", "l1000_opposite_best_rank", "l1000_support_score"]:
        table[column] = pd.to_numeric(table[column], errors="coerce")
    primary_rr = 1.0 / pd.to_numeric(table["primary_rank_1based"], errors="coerce")
    sensitivity_rr = 1.0 / pd.to_numeric(table["sensitivity_rank_1based"], errors="coerce")
    table["drugreflector_reciprocal_rank_score"] = (
        primary_rr.fillna(0.0) + sensitivity_rr.fillna(0.0)
    ) / 2.0
    table["crossvalidation_status"] = table.apply(classify_crossvalidation, axis=1)
    table["clue_query_status"] = clue_status
    table["clue_exact_match_available"] = table.get(
        "clue_match_type", pd.Series("", index=table.index)
    ).eq("exact_id")
    clue_component = pd.to_numeric(
        table.get("clue_tau", pd.Series(np.nan, index=table.index)),
        errors="coerce",
    ).fillna(0.0) / 100.0
    table["cross_method_score"] = (
        table["drugreflector_reciprocal_rank_score"]
        + table["l1000_support_score"].fillna(0.0)
        + clue_component
    )
    table["integrated_support_status"] = table.apply(
        classify_integrated_support, axis=1
    )
    table["metadata_mapped"] = table["pert_iname"].notna()
    return table.sort_values(
        ["cross_method_score", "in_both_top_lists", "mean_rank_1based", "compound"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)


def normalized_name_series(series: pd.Series) -> pd.Series:
    return (
        series.fillna("")
        .astype(str)
        .str.lower()
        .str.replace(r"[^a-z0-9]+", "", regex=True)
    )


def compute_entity_overlap(drugreflector_metadata: pd.DataFrame, l1000_metadata: pd.DataFrame) -> dict[str, int]:
    dr_ids = set(drugreflector_metadata["compound"].dropna().astype(str))
    l1000_ids = set(l1000_metadata["compound"].dropna().astype(str))
    dr_names = set(normalized_name_series(drugreflector_metadata["pert_iname"])) - {""}
    l1000_names = set(normalized_name_series(l1000_metadata["pert_iname"])) - {""}
    dr_inchi = set(drugreflector_metadata["inchi_key"].dropna().astype(str)) - {""}
    l1000_inchi = set(l1000_metadata["inchi_key"].dropna().astype(str)) - {""}
    return {
        "n_exact_brd_id_overlap": len(dr_ids.intersection(l1000_ids)),
        "n_normalized_name_overlap": len(dr_names.intersection(l1000_names)),
        "n_inchi_key_overlap": len(dr_inchi.intersection(l1000_inchi)),
    }


def read_json(path: Path, default: Any) -> Any:
    if not path.exists() or path.stat().st_size == 0:
        return default
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    args = parse_args()
    np.random.seed(args.seed)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)

    consensus = pd.read_csv(args.consensus, sep="\t")
    phase1 = read_pert_info(args.phase1_pert_info, "GSE92742_phase1", 1)
    phase2 = read_pert_info(args.phase2_pert_info, "GSE70138_phase2", 2)
    metadata = collapse_pert_metadata([phase1, phase2])

    l1000_status: dict[str, Any] = {"status": "not_refreshed"}
    l1000_topn: dict[str, Any] = {}
    if args.refresh_l1000fwd:
        payload = read_json(args.l1000_payload, {})
        try:
            search, l1000_topn, candidates = run_l1000fwd_query(
                args.l1000fwd_base_url,
                payload,
                args.timeout_seconds,
            )
            l1000_status = {"status": "completed", **search}
        except Exception as exc:
            l1000_status = {
                "status": "failed",
                "error_type": type(exc).__name__,
                "error": str(exc),
            }
            candidates = pd.read_csv(args.l1000_candidates, sep="\t")
    else:
        candidates = pd.read_csv(args.l1000_candidates, sep="\t")

    l1000_aggregate = aggregate_l1000fwd_candidates(candidates)
    clue_payload = read_json(args.clue_status, {"status": "missing_status_file"})
    clue_query_status = str(clue_payload.get("status", "unknown"))
    mapped_metadata = metadata.loc[metadata["compound"].isin(set(consensus["compound"]))].copy()
    clue_crosswalk = None
    clue_summary_mapped = pd.DataFrame()
    clue_cell = pd.DataFrame()
    if args.clue_results_archive.is_file():
        clue_summary, clue_cell = load_clue_results(args.clue_results_archive)
        clue_crosswalk, clue_summary_mapped = build_clue_crosswalk(
            mapped_metadata,
            clue_summary,
            clue_cell,
            metadata,
        )
        clue_query_status = "completed"
    crossvalidation = build_crossvalidation_table(
        consensus,
        metadata,
        l1000_aggregate,
        clue_query_status,
        clue_crosswalk=clue_crosswalk,
    )
    l1000_mapped = l1000_aggregate.merge(metadata, on="compound", how="left", validate="one_to_one")
    entity_overlap = compute_entity_overlap(mapped_metadata, l1000_mapped)

    outputs = {
        "perturbagen_metadata": args.metadata_dir / f"{OUTPUT_STEM}_perturbagen_metadata.tsv",
        "l1000fwd_candidates": args.metadata_dir / f"{OUTPUT_STEM}_l1000fwd_candidates.tsv",
        "l1000fwd_aggregate": args.metadata_dir / f"{OUTPUT_STEM}_l1000fwd_aggregate.tsv",
        "l1000fwd_mapped": args.metadata_dir / f"{OUTPUT_STEM}_l1000fwd_mapped.tsv",
        "l1000fwd_search": args.metadata_dir / f"{OUTPUT_STEM}_l1000fwd_search.json",
        "l1000fwd_topn": args.metadata_dir / f"{OUTPUT_STEM}_l1000fwd_topn.json",
        "clue_crosswalk": args.metadata_dir / f"{OUTPUT_STEM}_clue_crosswalk.tsv",
        "clue_perturbagen_summary": args.metadata_dir / f"{OUTPUT_STEM}_clue_perturbagen_summary.tsv",
        "clue_cell_scores": args.metadata_dir / f"{OUTPUT_STEM}_clue_cell_scores.tsv.gz",
        "clue_job_status": args.metadata_dir / "module9_8_clue_job_status.json",
        "clue_results_archive": args.clue_results_archive,
        "clue_up_entrez_gmt": args.metadata_dir / "module9_8_clue_up_entrez.gmt",
        "clue_down_entrez_gmt": args.metadata_dir / "module9_8_clue_down_entrez.gmt",
        "crossvalidation": args.metadata_dir / f"{OUTPUT_STEM}.tsv",
        "report": args.metadata_dir / f"{OUTPUT_STEM}_report.json",
    }
    mapped_metadata.to_csv(outputs["perturbagen_metadata"], sep="\t", index=False)
    candidates.to_csv(outputs["l1000fwd_candidates"], sep="\t", index=False)
    l1000_aggregate.to_csv(outputs["l1000fwd_aggregate"], sep="\t", index=False)
    l1000_mapped.to_csv(outputs["l1000fwd_mapped"], sep="\t", index=False)
    outputs["l1000fwd_search"].write_text(json.dumps(l1000_status, indent=2, sort_keys=True), encoding="utf-8")
    outputs["l1000fwd_topn"].write_text(json.dumps(l1000_topn, indent=2, sort_keys=True), encoding="utf-8")
    (clue_crosswalk if clue_crosswalk is not None else pd.DataFrame()).to_csv(
        outputs["clue_crosswalk"], sep="\t", index=False
    )
    clue_summary_mapped.to_csv(
        outputs["clue_perturbagen_summary"], sep="\t", index=False
    )
    clue_cell.to_csv(
        outputs["clue_cell_scores"], sep="\t", index=False, compression="gzip"
    )
    crossvalidation.to_csv(outputs["crossvalidation"], sep="\t", index=False)

    report = {
        "module": "module9_8_drugreflector_metadata_crossvalidation",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "completed",
        "inputs": {
            "consensus": str(args.consensus.resolve()),
            "phase1_pert_info": str(args.phase1_pert_info.resolve()),
            "phase2_pert_info": str(args.phase2_pert_info.resolve()),
            "l1000_payload": str(args.l1000_payload.resolve()),
            "clue_status": str(args.clue_status.resolve()),
            "clue_results_archive": str(args.clue_results_archive.resolve()),
            "refresh_l1000fwd": bool(args.refresh_l1000fwd),
            "seed": args.seed,
        },
        "outputs": {key: str(path.resolve()) for key, path in outputs.items()},
        "summary": {
            "n_drugreflector_consensus_compounds": int(len(consensus)),
            "n_metadata_mapped": int(crossvalidation["metadata_mapped"].sum()),
            "n_metadata_unmapped": int((~crossvalidation["metadata_mapped"]).sum()),
            "n_metadata_conflicts": int(crossvalidation["metadata_conflict_flag"].fillna(False).sum()),
            "n_l1000fwd_signature_hits": int(len(candidates)),
            "n_l1000fwd_unique_compounds": int(len(l1000_aggregate)),
            "n_l1000fwd_metadata_mapped": int(l1000_mapped["pert_iname"].notna().sum()),
            "n_exact_cross_method_compounds": int(
                crossvalidation["crossvalidation_status"].ne("drugreflector_only").sum()
            ),
            "n_l1000fwd_similar_supported": int(
                crossvalidation["crossvalidation_status"].isin(
                    ["supported_by_l1000fwd_similar", "mixed_l1000fwd_evidence"]
                ).sum()
            ),
            "n_l1000fwd_opposite_flagged": int(
                crossvalidation["crossvalidation_status"].isin(
                    ["discordant_l1000fwd_opposite", "mixed_l1000fwd_evidence"]
                ).sum()
            ),
            "l1000fwd_status": l1000_status.get("status"),
            "clue_status": clue_query_status,
            "n_clue_summary_compounds": int(len(clue_summary_mapped)),
            "n_clue_exact_id_matches": int(
                crossvalidation["clue_match_type"].eq("exact_id").sum()
            )
            if "clue_match_type" in crossvalidation
            else 0,
            "n_clue_inchi_matches": int(
                crossvalidation["clue_match_type"].eq("inchi_key").sum()
            )
            if "clue_match_type" in crossvalidation
            else 0,
            "n_clue_name_matches": int(
                crossvalidation["clue_match_type"].eq("normalized_name").sum()
            )
            if "clue_match_type" in crossvalidation
            else 0,
            "n_clue_unmapped": int(
                crossvalidation["clue_match_type"].eq("unmapped").sum()
            )
            if "clue_match_type" in crossvalidation
            else int(len(crossvalidation)),
            "n_clue_strong_support": int(
                crossvalidation["clue_strong_support"].fillna(False).sum()
            )
            if "clue_strong_support" in crossvalidation
            else 0,
            "n_clue_strong_opposition": int(
                crossvalidation["clue_strong_opposition"].fillna(False).sum()
            )
            if "clue_strong_opposition" in crossvalidation
            else 0,
            "n_clue_weak_positive": int(
                crossvalidation["clue_positive"].fillna(False).sum()
            )
            if "clue_positive" in crossvalidation
            else 0,
            "n_clue_weak_negative": int(
                crossvalidation["clue_negative"].fillna(False).sum()
            )
            if "clue_negative" in crossvalidation
            else 0,
            "n_three_method_strong_support": int(
                crossvalidation["integrated_support_status"]
                .eq("three_method_strong_support")
                .sum()
            )
            if "integrated_support_status" in crossvalidation
            else 0,
            **entity_overlap,
        },
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "pandas": package_version("pandas"),
            "numpy": package_version("numpy"),
        },
        "data_sources": [
            {
                "name": "LINCS Phase I perturbagen metadata",
                "accession": "GSE92742",
                "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_pert_info.txt.gz",
            },
            {
                "name": "LINCS Phase II perturbagen metadata",
                "accession": "GSE70138",
                "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz",
            },
            {
                "name": "L1000FWD signature search",
                "url": args.l1000fwd_base_url,
            },
            {
                "name": "CLUE Connectivity Map sig_gutc_tool",
                "url": "https://api.clue.io/api/jobs",
                "tool_version": "1.1.1.2",
            },
        ],
    }
    outputs["report"].write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(report["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
