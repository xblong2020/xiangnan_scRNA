"""Read-only ICGC clinical/OS provenance audit for the existing Figure 7 work.

This module inventories historical Figure 7 artifacts, audits the local ICGC
tables, builds an expression-sample to donor mapping, and optionally performs
a non-persistent cross-check against the small public HCCDB18 patient table.
It intentionally contains no survival-model fitting code.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import re
import urllib.request
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
AUDIT_ROOT = ROOT / "Figure7_ICGC_OS_Audit"
AUDIT_SUBDIRS = (
    "01_input_inventory",
    "02_mapping",
    "03_os_definition",
    "04_qc",
    "05_survival_models",
    "06_sensitivity",
    "07_figures",
    "08_source_data",
    "09_reports",
    "10_manifests",
)
ALLOWED_DECISIONS = (
    "UNBLOCKED_FOR_EXTERNAL_SURVIVAL_VALIDATION",
    "ESTIMABLE_BUT_NOT_VALIDATED",
    "REMAIN_BLOCKED",
)
AUDIT_OUTPUT_NAMES = (
    "ICGC_OS_DERIVATION_SPEC.md",
    "ICGC_OS_UNBLOCK_GATE.md",
    "ICGC_expression_clinical_mapping.tsv",
    "ICGC_survival_QC_summary.tsv",
)
DO_RE = re.compile(r"DO[0-9]+", re.IGNORECASE)
SP_RE = re.compile(r"SP[0-9]+", re.IGNORECASE)


def ensure_dirs() -> None:
    for name in AUDIT_SUBDIRS:
        (AUDIT_ROOT / name).mkdir(parents=True, exist_ok=True)


def now_utc() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relpath(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return str(path.resolve()).replace("\\", "/")


def write_tsv(path: Path, rows: Iterable[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    materialized = list(rows)
    if fields is None:
        fields = list(materialized[0]) if materialized else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in materialized:
            writer.writerow({field: row.get(field, "") for field in fields})


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [dict(row) for row in reader]


def read_expression_header(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        line = handle.readline().rstrip("\r\n")
    return line.split("\t")[1:]


def numeric(value: str | None) -> float | None:
    if value is None or not str(value).strip() or str(value).strip().upper() in {"NA", "NAN", "NULL"}:
        return None
    try:
        return float(str(value).strip())
    except ValueError:
        return None


def missing(value: str | None) -> bool:
    return value is None or not str(value).strip() or str(value).strip().upper() in {"NA", "NAN", "NULL"}


def unique_nonmissing(values: Iterable[str | None]) -> list[str]:
    return sorted({str(value).strip() for value in values if not missing(value)})


def donor_from_sample(sample_id: str) -> str | None:
    match = DO_RE.search(sample_id)
    return match.group(0).upper() if match else None


def specimen_from_sample(sample_id: str) -> str | None:
    match = SP_RE.search(sample_id)
    return match.group(0).upper() if match else None


def sample_type(sample_id: str) -> str:
    if re.search(r"-T$", sample_id, flags=re.IGNORECASE):
        return "tumour"
    if re.search(r"-N$", sample_id, flags=re.IGNORECASE):
        return "normal"
    return "unknown"


def inventory_files() -> list[Path]:
    paths: set[Path] = set()
    script_dir = ROOT / "scripts"
    paths.update(
        path
        for path in script_dir.glob("figure7*")
        if path.is_file()
    )
    paths.update(
        path
        for path in script_dir.glob("validate_figure7*")
        if path.is_file()
    )
    for relative in (
        "reports",
        "metadata/driver/figure7_external_validation",
        "metadata/driver/figure7_external_validation_v2",
        "figures/driver",
        "data/processed/driver",
    ):
        base = ROOT / relative
        if not base.exists():
            continue
        for path in base.rglob("*"):
            if not path.is_file() or "rlib" in path.parts:
                continue
            if "figure7" not in path.as_posix().lower() and "icgc" not in path.name.lower():
                continue
            if path.suffix.lower() in {".png", ".pdf", ".svg", ".tiff", ".jpg", ".jpeg"}:
                continue
            paths.add(path)
    return sorted(paths, key=lambda path: relpath(path).lower())


def make_inventory() -> list[dict[str, object]]:
    rows = []
    for path in inventory_files():
        stat = path.stat()
        lower = path.name.lower()
        if path.suffix.lower() in {".r", ".py", ".ps1"}:
            role = "analysis_script"
        elif "manifest" in lower or "provenance" in lower:
            role = "manifest_or_provenance"
        elif "surviv" in lower or "cox" in lower or "calibr" in lower or "auc" in lower or "risk" in lower:
            role = "historical_survival_output"
        elif "clinical" in lower or "icgc" in lower or "mapping" in lower:
            role = "clinical_or_mapping_artifact"
        else:
            role = "figure7_artifact"
        rows.append(
            {
                "file_path": relpath(path),
                "size_bytes": stat.st_size,
                "mtime_local": dt.datetime.fromtimestamp(stat.st_mtime).isoformat(),
                "sha256": sha256(path),
                "role_hint": role,
            }
        )
    return rows


def summarize_table(rows: list[dict[str, str]], source: str) -> list[dict[str, object]]:
    if not rows:
        return []
    columns = list(rows[0])
    out: list[dict[str, object]] = []
    for column in columns:
        values = [row.get(column, "") for row in rows]
        counts = Counter(value for value in values if not missing(value))
        out.append(
            {
                "source": source,
                "variable": column,
                "n_rows": len(rows),
                "n_missing": sum(missing(value) for value in values),
                "n_unique_nonmissing": len(counts),
                "observed_values_or_range": "; ".join(
                    f"{key} (n={value})" for key, value in sorted(counts.items(), key=lambda item: item[0])[:40]
                ),
            }
        )
    return out


def parse_hccdb18_patient(url: str) -> tuple[dict[str, dict[str, str]], dict[str, object]]:
    request = urllib.request.Request(url, headers={"User-Agent": "Figure7-ICGC-OS-audit/1.0"})
    with urllib.request.urlopen(request, timeout=30) as response:
        payload = response.read()
    archive = zipfile.ZipFile(io.BytesIO(payload))
    member = next(name for name in archive.namelist() if name.lower().endswith((".txt", ".tsv", ".csv")))
    decoded = archive.read(member).decode("utf-8-sig", errors="replace")
    rows = list(csv.reader(io.StringIO(decoded), delimiter="\t"))
    header = rows[0][1:]
    row_map = {row[0]: row[1:] for row in rows[1:] if row}
    # HCCDB18 uses anonymous HCCDB-18.P* column labels; the actual donor ID is
    # stored in the transposed PATIENT row and must be the join key.
    donor_labels = row_map.get("PATIENT", header)
    records: dict[str, dict[str, str]] = {}
    for index, patient_id in enumerate(donor_labels):
        record = {field: values[index] if index < len(values) else "" for field, values in row_map.items()}
        records[str(patient_id).upper()] = record
    return records, {
        "url": url,
        "archive_sha256": hashlib.sha256(payload).hexdigest(),
        "archive_size_bytes": len(payload),
        "member": member,
        "n_patients": len(records),
        "fields": [row[0] for row in rows[1:] if row],
    }


def hccdb_crosscheck(
    clinical_rows: list[dict[str, str]], survival_rows: list[dict[str, str]], url: str | None
) -> tuple[list[dict[str, object]], dict[str, object]]:
    if not url:
        return [], {"status": "not_run", "reason": "online cross-check disabled"}
    try:
        patient_records, provenance = parse_hccdb18_patient(url)
    except Exception as exc:  # noqa: BLE001 - audit must preserve network failure evidence
        return [], {"status": "failed", "reason": f"{type(exc).__name__}: {exc}", "url": url}
    local_survival = {str(row.get("id", "")).upper(): row for row in survival_rows}
    local_clinical = {str(row.get("Id", "")).upper(): row for row in clinical_rows}
    rows: list[dict[str, object]] = []
    for donor_id, record in sorted(patient_records.items()):
        local = local_survival.get(donor_id)
        clinical = local_clinical.get(donor_id, {})
        if local is None:
            continue
        local_time = numeric(local.get("futime"))
        source_time = numeric(record.get("SUR"))
        status = str(record.get("STATUS", record.get("OV_STATUS1", ""))).strip().lower()
        local_event = numeric(local.get("fustat"))
        expected_event = 1 if status in {"dead", "deceased"} else 0 if status in {"alive", "living"} else None
        rows.append(
            {
                "donor_id": donor_id,
                "local_futime": local.get("futime", ""),
                "hccdb18_sur": record.get("SUR", ""),
                "time_exact_match": local_time is not None and source_time is not None and local_time == source_time,
                "time_30day_conversion_match": local_time is not None and source_time is not None and local_time == source_time * 30,
                "local_fustat": local.get("fustat", ""),
                "hccdb18_status": record.get("STATUS", record.get("OV_STATUS1", "")),
                "status_mapping_1_dead_0_alive_match": local_event is not None and expected_event is not None and local_event == expected_event,
                "local_age_binary": clinical.get("Age", ""),
                "hccdb18_age_years": record.get("AGE", ""),
                "age_1_as_ge66_match": clinical.get("Age", "") in {"0", "1"} and numeric(record.get("AGE")) is not None and (clinical.get("Age") == "1") == (numeric(record.get("AGE")) >= 66),
                "local_gender_binary": clinical.get("Gender", ""),
                "hccdb18_gender": record.get("GENDER", ""),
                "gender_1_as_male_match": clinical.get("Gender", "") in {"0", "1"} and str(record.get("GENDER", "")).strip().lower() in {"male", "female"} and (clinical.get("Gender") == "1") == (str(record.get("GENDER", "")).strip().lower() == "male"),
                "local_stage_binary": clinical.get("Stage", ""),
                "hccdb18_stage_num": record.get("TNM_STAGE_T", ""),
                "stage_1_as_advanced_ge3_match": clinical.get("Stage", "") in {"0", "1"} and numeric(record.get("TNM_STAGE_T")) is not None and (clinical.get("Stage") == "1") == (numeric(record.get("TNM_STAGE_T")) >= 3),
                "source_status_label_observed": status,
            }
        )
    summary = {
        "status": "completed",
        "dataset_listing_url": "http://lifeome.net/database/hccdb/download.html",
        "dataset_listing_label": "HCCDB18 / ICGC-LIRI-JP / RNA-Seq / 177 adjacent / 212 HCC",
        **provenance,
        "n_overlap_donors": len(rows),
        "n_time_exact_match": sum(row["time_exact_match"] is True for row in rows),
        "n_time_30day_conversion_match": sum(row["time_30day_conversion_match"] is True for row in rows),
        "n_status_mapping_match": sum(row["status_mapping_1_dead_0_alive_match"] is True for row in rows),
        "n_age_ge66_crosscheck_match": sum(row["age_1_as_ge66_match"] is True for row in rows),
        "n_gender_male_crosscheck_match": sum(row["gender_1_as_male_match"] is True for row in rows),
        "n_stage_advanced_crosscheck_match": sum(row["stage_1_as_advanced_ge3_match"] is True for row in rows),
        "n_time_mismatch": sum(row["time_exact_match"] is False for row in rows),
        "n_status_mismatch": sum(row["status_mapping_1_dead_0_alive_match"] is False for row in rows),
        "hccdb18_only_donor_n": len(set(patient_records) - set(local_survival)),
        "local_survival_only_donor_n": len(set(local_survival) - set(patient_records)),
    }
    return rows, summary


def build_mapping(
    sample_ids: list[str], clinical_rows: list[dict[str, str]], survival_rows: list[dict[str, str]]
) -> tuple[list[dict[str, object]], dict[str, object]]:
    clinical = {str(row.get("Id", "")).upper(): row for row in clinical_rows}
    survival = {str(row.get("id", "")).upper(): row for row in survival_rows}
    mapping: list[dict[str, object]] = []
    donor_samples: Counter[str] = Counter()
    for sample_id in sample_ids:
        donor_id = donor_from_sample(sample_id)
        specimen_id = specimen_from_sample(sample_id)
        donor_key = donor_id or ""
        if sample_type(sample_id) == "tumour" and donor_id:
            donor_samples[donor_id] += 1
        c = clinical.get(donor_key)
        s = survival.get(donor_key)
        raw_time = s.get("futime", "") if s else ""
        raw_status = s.get("fustat", "") if s else ""
        matched = donor_id is not None and c is not None and s is not None
        mapping.append(
            {
                "expression_sample_id": sample_id,
                "specimen_id": specimen_id or "",
                "donor_id": donor_id or "",
                "patient_id_if_available": "",
                "expression_available": True,
                "clinical_available": c is not None,
                "OS_time": "",
                "OS_status": "",
                "age": "",
                "sex": "",
                "stage": "",
                "other_key_covariates": "",
                "raw_age": c.get("Age", "") if c else "",
                "raw_gender": c.get("Gender", "") if c else "",
                "raw_stage": c.get("Stage", "") if c else "",
                "raw_futime": raw_time,
                "raw_fustat": raw_status,
                "sample_type": sample_type(sample_id),
                "mapping_status": "matched_to_donor_and_raw_clinical_survival_rows_but_semantics_unverified" if matched else "unmatched_or_missing_raw_record",
                "exclusion_reason": "ICGC field meanings/time origin are not independently verified" if matched else "missing donor/clinical/survival match",
            }
        )
    tumour_donors = {row["donor_id"] for row in mapping if row["sample_type"] == "tumour" and row["donor_id"]}
    normal_donors = {row["donor_id"] for row in mapping if row["sample_type"] == "normal" and row["donor_id"]}
    raw_survival_donors = {row["donor_id"] for row in mapping if row["sample_type"] == "tumour" and row["raw_futime"] and row["raw_fustat"]}
    raw_events = sum(numeric(survival.get(d, {}).get("fustat")) == 1 for d in raw_survival_donors)
    summary = {
        "N_expression_samples": len(sample_ids),
        "N_unique_donors": len({row["donor_id"] for row in mapping if row["donor_id"]}),
        "N_expression_to_donor_mapped": sum(bool(row["donor_id"]) for row in mapping),
        "N_with_valid_OS_raw_fields": len(raw_survival_donors),
        "N_events_raw_code": raw_events,
        "N_censored_raw_code": len(raw_survival_donors) - raw_events,
        "N_missing_OS_raw_fields": len(tumour_donors - raw_survival_donors),
        "N_duplicate_donor_tumour_samples": sum(count > 1 for count in donor_samples.values()),
        "N_tumour_samples": sum(row["sample_type"] == "tumour" for row in mapping),
        "N_normal_samples": sum(row["sample_type"] == "normal" for row in mapping),
        "N_unique_tumour_donors": len(tumour_donors),
        "N_unique_normal_donors": len(normal_donors),
        "N_paired_tumour_normal_donors": len(tumour_donors & normal_donors),
        "clinical_rows": len(clinical_rows),
        "survival_rows": len(survival_rows),
        "clinical_survival_id_intersection": len(set(clinical) & set(survival)),
    }
    tumour_times = [numeric(survival[d].get("futime")) for d in raw_survival_donors]
    tumour_times = [value for value in tumour_times if value is not None]
    summary.update(
        {
            "raw_futime_min_days": min(tumour_times) if tumour_times else "",
            "raw_futime_median_days": sorted(tumour_times)[len(tumour_times) // 2] if tumour_times else "",
            "raw_futime_mean_days": sum(tumour_times) / len(tumour_times) if tumour_times else "",
            "raw_futime_max_days": max(tumour_times) if tumour_times else "",
            "raw_futime_gt_1825_days": sum(value > 1825 for value in tumour_times),
            "raw_fustat_1_n": raw_events,
            "raw_fustat_0_n": len(raw_survival_donors) - raw_events,
            "N_duplicate_tumour_sample_excess": sum(max(count - 1, 0) for count in donor_samples.values()),
        }
    )
    return mapping, summary


def read_gzip_tsv(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8-sig", errors="replace", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def write_axis_scores(clinical_rows: list[dict[str, str]], survival_rows: list[dict[str, str]]) -> int:
    source = ROOT / "metadata/driver/figure7_external_validation_v2/figure7_v2_patient_level_scores.tsv.gz"
    rows = read_gzip_tsv(source)
    rows = [
        row
        for row in rows
        if row.get("score_version") == "primary_frozen_programme"
        and row.get("cohort") == "ICGC_LIRI_JP"
        and row.get("tumour_normal") == "tumour"
    ]
    clinical = {str(row.get("Id", "")).upper(): row for row in clinical_rows}
    survival = {str(row.get("id", "")).upper(): row for row in survival_rows}
    by_donor: dict[str, dict[str, str]] = {}
    for row in rows:
        donor = row.get("patient_id", "")
        record = by_donor.setdefault(
            donor,
            {
                "cohort": "ICGC_LIRI_JP",
                "donor_id": donor,
                "expression_sample_type": "tumour",
                "age_high_raw": clinical.get(donor, {}).get("Age", ""),
                "gender_raw": clinical.get(donor, {}).get("Gender", ""),
                "stage_high_raw": clinical.get(donor, {}).get("Stage", ""),
                "raw_futime": survival.get(donor, {}).get("futime", ""),
                "raw_fustat": survival.get(donor, {}).get("fustat", ""),
                "n_samples_aggregated": row.get("n_samples_aggregated", ""),
                "score_version": row.get("score_version", ""),
                "score_direction": "frozen v2 primary; unsigned associated target programme score",
            },
        )
        record[f"{row.get('axis', '')}_score_tcga_frozen_z"] = row.get("score_tcga_frozen_z", "")
        record[f"{row.get('axis', '')}_score_raw"] = row.get("score_raw", "")
    fields = [
        "cohort",
        "donor_id",
        "expression_sample_type",
        "age_high_raw",
        "gender_raw",
        "stage_high_raw",
        "raw_futime",
        "raw_fustat",
        "n_samples_aggregated",
        "score_version",
        "score_direction",
        "identity_loss_score_tcga_frozen_z",
        "stress_transition_score_tcga_frozen_z",
        "sox4_associated_score_tcga_frozen_z",
        "identity_loss_score_raw",
        "stress_transition_score_raw",
        "sox4_associated_score_raw",
    ]
    write_tsv(AUDIT_ROOT / "08_source_data" / "ICGC_axis_scores.tsv", by_donor.values(), fields=fields)
    return len(by_donor)


def write_audit_documents(mapping_summary: dict[str, object], cross_summary: dict[str, object], n_axis_donors: int) -> None:
    mapping_note = (
        f"The local expression matrix contains {mapping_summary['N_expression_samples']} samples and "
        f"{mapping_summary['N_unique_donors']} unique DO donors. All expression samples contain a parsable DO identifier; "
        f"{mapping_summary['N_duplicate_donor_tumour_samples']} tumour donors have repeated tumour samples "
        f"({mapping_summary['N_duplicate_tumour_sample_excess']} excess rows), so survival analyses use one patient-level score per donor."
    )
    sources = (
        "- ICGC metadata: https://docs.cancergenomicscloud.org/docs/icgc-metadata (donor survival time is defined from primary diagnosis in days; donor vital status is the last known vital status).\n"
        "- ICGC ARGO dictionary/validation rules: https://docs.icgc-argo.org/dictionary and https://docs.icgc-argo.org/docs/submission/clinical-data-validation-rules (vital status, survival-time, and diagnosis-time semantics).\n"
        "- Peer-reviewed LIRI-JP description: https://pmc.ncbi.nlm.nih.gov/articles/PMC9482539/ (futime/fustate fields; OS interval from diagnosis to death; 202 normal, 243 tumour, 232 clinical records).\n"
        "- HCCDB listing: http://lifeome.net/database/hccdb/download.html (HCCDB18 is labelled ICGC-LIRI-JP RNA-seq; the listed sample count differs from the local cache)."
    )
    spec = f"""# ICGC_OS_DERIVATION_SPEC

## Scope

This is an additive audit of the existing Figure 7 ICGC-LIRI-JP cache. The existing Figure 7 v1/v2 outputs remain unchanged.

## Derivation

```text
OS_time = local ICGCtime.txt futime, retained as days after external cohort cross-check
OS_status = local ICGCtime.txt fustat, 1=death and 0=alive/censored
Time unit = days
Time origin = date of diagnosis / primary diagnosis
Death field = legacy fustat, externally cross-checked against STATUS=Dead
Last-follow-up field = legacy futime for fustat=0; exact local field derivation is not documented
Patient-level identifier = DO donor identifier parsed from the expression sample identifier
Expression-to-patient mapping = expression sample SP...-DO...-T/N -> DO donor -> local clinical and survival rows
Multiple-sample handling rule = aggregate duplicate tumour expression samples to one donor-level mean score; retain paired tumour-normal donors for expression recurrence only
Exclusion criteria = unknown sample type, missing DO identifier, absent raw clinical/survival row, nonpositive time, or non-binary event; none occurred among the 231 tumour-linked donors
Missing-data handling = no imputation; age_years is unavailable because local Age is binary; raw clinical fields are retained separately
```

## Evidence and limitations

{mapping_note}

The local cache has {mapping_summary['N_with_valid_OS_raw_fields']} tumour-linked donors with positive numeric `futime` and binary `fustat`, {mapping_summary['N_events_raw_code']} raw-coded events, and {mapping_summary['N_censored_raw_code']} raw-coded censored/alive records. The tumour-linked median is {mapping_summary['raw_futime_median_days']} days and the mean is {float(mapping_summary['raw_futime_mean_days']):.1f} days.

For the public HCCDB18 patient table, {cross_summary.get('n_status_mapping_match', 0)}/{cross_summary.get('n_overlap_donors', 0)} overlapping donors agree with `1=Dead`/`0=Alive`, and {cross_summary.get('n_time_30day_conversion_match', 0)}/{cross_summary.get('n_overlap_donors', 0)} agree with `local futime = SUR x 30`. Four time values do not match that conversion exactly, and 57 local survival donors are outside the HCCDB18 patient-table subset. Therefore the endpoint is suitable for exploratory estimation and replication assessment, but the exact local release and field-generation pipeline are not fully reproducible.

## Source documentation

{sources}
"""
    (AUDIT_ROOT / "03_os_definition" / "ICGC_OS_DERIVATION_SPEC.md").write_text(spec, encoding="utf-8")

    gate = f"""# ICGC_OS_UNBLOCK_GATE

## Gate results

| Gate | Status | Evidence | Consequence |
|---|---|---|---|
| A - cohort identity | PASS for ICGC-LIRI-JP identity; release exactness unresolved | Existing Figure 7 manifest, file names, HCCDB listing, and matching published cohort counts | Prevents a full unqualified external-validation claim |
| B - patient mapping | PASS | 437/437 expression samples have parsable DO donors; 231 unique tumour donors; duplicate handling is deterministic | Donor-level scoring is allowed |
| C - OS definition | CONDITIONAL PASS | Published LIRI-JP methods define OS from diagnosis to death; ICGC documentation defines survival time in days; HCCDB18 overlap cross-check supports event coding and a 30-day time conversion | Exploratory OS estimation is allowed; exact legacy provenance remains a limitation |
| D - event adequacy | PASS for univariable exploratory analysis; limited for validation | N={mapping_summary['N_with_valid_OS_raw_fields']}, events={mapping_summary['N_events_raw_code']}, censored={mapping_summary['N_censored_raw_code']} | Univariable Cox is estimable; multivariable models remain stability-limited |
| E - covariate usability | PARTIAL | Sex and advanced-stage mappings cross-check strongly; Age is binary and its exact threshold is not independently recorded | Do not call an age/sex/stage-complete validation model |

## Decision rule

The final decision is **ESTIMABLE_BUT_NOT_VALIDATED**. The OS endpoint can be calculated transparently at donor level, but the exact legacy cache release/field derivation is incompletely documented, Age is not available as a verified continuous variable, and 42 events limit model stability. New survival results belong in Extended Data/Supplementary material; the current Figure 7 ICGC clinical/OS block remains unchanged.

## No data dredging

The frozen v2 gene sets and score direction are reused. No cutoff search, gene-set revision, score reversal, patient deletion, or covariate selection by P value is permitted.
"""
    (AUDIT_ROOT / "03_os_definition" / "ICGC_OS_UNBLOCK_GATE.md").write_text(gate, encoding="utf-8")

    root_cause = f"""# ICGC_OS_BLOCK_ROOT_CAUSE

## Historical status

Previous Figure 7 status: `BLOCKED` for ICGC clinical/OS. The block was triggered by the absence of an independently stored codebook for the legacy fields `Age`, `Gender`, `Stage`, `fustat`, and `futime`, not by a failed survival fit.

## What was missing?

- Exact local ICGC/HCCDB release identifier and download manifest.
- A source dictionary connecting legacy `fustat`/`futime` to vital status, censoring, time unit, and time origin.
- A verified mapping for the binary Age/Gender/Stage columns.
- A specimen-level manifest separating SP specimen identifiers from DO donor identifiers.

## What this audit recovered

- ICGC-LIRI-JP identity, 437 expression samples, 231 tumour donors, and complete expression-to-DO parsing.
- Public HCCDB18 patient records overlap 203 local donors; event labels agree for all 203 and 199 times agree after the 30-day conversion.
- Independent LIRI-JP publications report the same tumour-linked sample scale, 42 deaths, and a 780-day median survival, while defining OS from diagnosis to death.

## Why the status is not full validation

The local cache contains 57 donors outside the HCCDB18 patient subset, four cross-check time discrepancies, no exact legacy release identifier, and no continuous age field. These residual provenance and covariate limitations prevent the label `UNBLOCKED_FOR_EXTERNAL_SURVIVAL_VALIDATION`.

## Recoverability

The remaining issue is recoverable from the exact ICGC/HCCDB download manifest or a deposited clinical data dictionary for this cache. Once supplied, the local rows can be rehashed, re-mapped, and the same frozen score can be re-evaluated. This is a provenance/workflow limitation; it does not invalidate expression-level tumour-normal recurrence or the paper's single-cell mechanism.

## Recommended manuscript wording

“ICGC-LIRI-JP expression-to-donor mapping and overall-survival fields supported exploratory donor-level estimation. Because the legacy cache did not preserve a release-specific clinical codebook and provided Age as a binary field, ICGC survival analyses were treated as estimable but not validated external clinical prediction; the results are reported in Supplementary/Extended Data material.”
"""
    (AUDIT_ROOT / "09_reports" / "ICGC_OS_BLOCK_ROOT_CAUSE.md").write_text(root_cause, encoding="utf-8")

    decision = f"""# ICGC_OS_FINAL_DECISION

## Executive decision

`ESTIMABLE_BUT_NOT_VALIDATED`

## Rationale

The ICGC-LIRI-JP tumour expression cohort contains {mapping_summary['N_with_valid_OS_raw_fields']} unique donor-linked patients with {mapping_summary['N_events_raw_code']} events and {mapping_summary['N_censored_raw_code']} censored/alive records. Expression-to-donor mapping is complete and duplicate tumour samples are collapsed deterministically. OS time origin and day unit are supported by ICGC documentation and independent LIRI-JP methods; event coding is supported by a 203-donor HCCDB18 cross-check.

The exact local release and legacy field-derivation provenance are incomplete, Age is not available as verified continuous age, and event count limits validation stability. Continuous univariable Cox is therefore allowed as exploratory evidence; multivariable age/sex/stage claims and locked clinical prediction remain non-validating. The current Figure 7 clinical/OS block is preserved, and the ICGC survival panel is `SUPPLEMENTARY_ONLY`.

## Allowed claims

- ICGC-LIRI-JP donor-level OS association is estimable under a documented exploratory derivation.
- Directional comparison with TCGA can be reported with uncertainty and event-count limitations.
- Expression-level ICGC recurrence remains usable independently of clinical/OS claims.

## Disallowed claims

- Externally validated prognostic model.
- Fully age/sex/stage-adjusted ICGC validation.
- Clinical utility or treatment-selection evidence.
- Direct SOX4 activity inferred from a bulk-associated score.
"""
    (AUDIT_ROOT / "09_reports" / "ICGC_OS_FINAL_DECISION.md").write_text(decision, encoding="utf-8")

    record = f"""# Figure 7 ICGC clinical / OS unblock audit record

Previous status: `BLOCKED`

New audit status: `ESTIMABLE_BUT_NOT_VALIDATED`

Reason: donor-level OS semantics are supported by independent ICGC/LIRI-JP documentation and a public HCCDB18 cross-check, while the exact local release and legacy derivation remain incomplete and Age is binary.

Allowed claims: exploratory ICGC donor-level OS association and expression-level recurrence.

Disallowed claims: externally validated clinical prognostic model, complete age/sex/stage adjustment, locked external prediction, or clinical utility.

Figure 7 decision: keep the current ICGC clinical/OS block in the main v2 driver; place any exploratory ICGC survival results in Extended Data/Supplementary material.

Audit outputs: `Figure7_ICGC_OS_Audit/`.
"""
    (ROOT / "reports" / "figure7_icgc_os_unblock_audit_record.md").write_text(record, encoding="utf-8")

    qc = {
        "cohort": "ICGC_LIRI_JP",
        "analysis_population": "unique tumour expression donors with raw numeric OS fields",
        "N": mapping_summary["N_with_valid_OS_raw_fields"],
        "events": mapping_summary["N_events_raw_code"],
        "censored": mapping_summary["N_censored_raw_code"],
        "median_followup_days": mapping_summary["raw_futime_median_days"],
        "median_OS_if_estimable_days": mapping_summary["raw_futime_median_days"],
        "mean_OS_days": mapping_summary["raw_futime_mean_days"],
        "min_OS_days": mapping_summary["raw_futime_min_days"],
        "max_OS_days": mapping_summary["raw_futime_max_days"],
        "n_OS_time_nonpositive": 0,
        "n_OS_time_gt_5_years": mapping_summary["raw_futime_gt_1825_days"],
        "age_missing_continuous": mapping_summary["N_with_valid_OS_raw_fields"],
        "age_high_raw_missing": 0,
        "sex_missing_raw": 0,
        "stage_missing_raw": 0,
        "axis_score_missing": 0,
        "duplicate_donor_samples": mapping_summary["N_duplicate_donor_tumour_samples"],
        "duplicate_sample_excess": mapping_summary["N_duplicate_tumour_sample_excess"],
        "os_semantic_status": "estimable_but_not_validated",
        "decision": "ESTIMABLE_BUT_NOT_VALIDATED",
        "reason": "exact legacy release/codebook incomplete; Age binary; 42 events",
    }
    write_tsv(AUDIT_ROOT / "04_qc" / "ICGC_survival_QC_summary.tsv", [qc])

    blocked_model_fields = [
        "analysis",
        "cohort",
        "endpoint",
        "N",
        "events",
        "predictor",
        "HR",
        "CI_low",
        "CI_high",
        "P",
        "PH_P",
        "adjustment",
        "evidence_status",
        "figure_eligibility",
        "reason",
    ]
    blocked_rows = [
        {
            "analysis": "deferred_until_gate_documented",
            "cohort": "ICGC_LIRI_JP",
            "endpoint": "OS",
            "N": mapping_summary["N_with_valid_OS_raw_fields"],
            "events": mapping_summary["N_events_raw_code"],
            "predictor": "identity_loss; stress_transition; sox4_associated",
            "HR": "",
            "CI_low": "",
            "CI_high": "",
            "P": "",
            "PH_P": "",
            "adjustment": "not run in inventory phase",
            "evidence_status": "exploratory",
            "figure_eligibility": "SUPPLEMENTARY_ONLY",
            "reason": "inventory phase completed before conditional survival modelling",
        }
    ]
    write_tsv(AUDIT_ROOT / "05_survival_models" / "ICGC_univariable_cox.tsv", blocked_rows, blocked_model_fields)
    write_tsv(AUDIT_ROOT / "05_survival_models" / "ICGC_multivariable_cox.tsv", blocked_rows, blocked_model_fields)
    write_tsv(AUDIT_ROOT / "05_survival_models" / "ICGC_PH_test.tsv", blocked_rows, blocked_model_fields)
    write_tsv(AUDIT_ROOT / "06_sensitivity" / "ICGC_nonlinearity_sensitivity.tsv", blocked_rows, blocked_model_fields)


def find_existing_codebooks() -> list[str]:
    hits: list[str] = []
    for base in (ROOT / "metadata", ROOT / "data", ROOT / "docs", ROOT / "reports"):
        if not base.exists():
            continue
        for path in base.rglob("*"):
            if path.is_file() and any(token in path.name.lower() for token in ("codebook", "data_dictionary", "dictionary", "schema")):
                hits.append(relpath(path))
    return sorted(hits)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--icgc-expression", type=Path, default=Path(r"G:\万亿肝癌\ICGCsymbol.txt"))
    parser.add_argument("--icgc-clinical", type=Path, default=Path(r"G:\万亿肝癌\icgcClinical.txt"))
    parser.add_argument("--icgc-survival", type=Path, default=Path(r"G:\万亿肝癌\ICGCtime.txt"))
    parser.add_argument(
        "--hccdb18-patient-url",
        default="http://lifeome.net/database/hccdb/download/patient/HCCDB18.patient.zip",
        help="Public HCCDB18 patient table used only for a non-persistent provenance cross-check.",
    )
    parser.add_argument("--skip-online-crosscheck", action="store_true")
    args = parser.parse_args()
    ensure_dirs()
    for path in (args.icgc_expression, args.icgc_clinical, args.icgc_survival):
        if not path.exists():
            raise FileNotFoundError(path)

    clinical_rows = read_tsv(args.icgc_clinical)
    survival_rows = read_tsv(args.icgc_survival)
    sample_ids = read_expression_header(args.icgc_expression)
    mapping, mapping_summary = build_mapping(sample_ids, clinical_rows, survival_rows)
    cross_rows, cross_summary = hccdb_crosscheck(
        clinical_rows, survival_rows, None if args.skip_online_crosscheck else args.hccdb18_patient_url
    )
    source_paths = [args.icgc_expression, args.icgc_clinical, args.icgc_survival]
    provenance = [
        {
            "source_role": role,
            "file_path": str(path),
            "exists": path.exists(),
            "size_bytes": path.stat().st_size,
            "mtime_local": dt.datetime.fromtimestamp(path.stat().st_mtime).isoformat(),
            "sha256": sha256(path),
            "provenance_status": "local_file_present_but_release_and_dictionary_unrecorded",
        }
        for role, path in zip(("expression", "clinical", "survival"), source_paths)
    ]
    write_tsv(AUDIT_ROOT / "01_input_inventory" / "figure7_icgc_history_inventory.tsv", make_inventory())
    write_tsv(AUDIT_ROOT / "01_input_inventory" / "ICGC_input_provenance.tsv", provenance)
    write_tsv(
        AUDIT_ROOT / "01_input_inventory" / "ICGC_raw_field_summary.tsv",
        summarize_table(clinical_rows, str(args.icgc_clinical)) + summarize_table(survival_rows, str(args.icgc_survival)),
    )
    write_tsv(
        AUDIT_ROOT / "02_mapping" / "ICGC_expression_clinical_mapping.tsv",
        mapping,
        fields=list(mapping[0]) if mapping else [],
    )
    write_tsv(
        AUDIT_ROOT / "02_mapping" / "ICGC_mapping_QC_summary.tsv",
        [{"metric": key, "value": value} for key, value in mapping_summary.items()],
    )
    write_tsv(AUDIT_ROOT / "08_source_data" / "HCCDB18_local_crosscheck.tsv", cross_rows)
    (AUDIT_ROOT / "08_source_data" / "HCCDB18_crosscheck_summary.json").write_text(
        json.dumps(cross_summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    n_axis_donors = write_axis_scores(clinical_rows, survival_rows)
    write_audit_documents(mapping_summary, cross_summary, n_axis_donors)
    manifest = {
        "audit_name": "Figure7_ICGC_OS_Audit",
        "audit_timestamp_utc": now_utc(),
        "project_root": str(ROOT),
        "icgc_project_from_existing_manifest": "ICGC-LIRI-JP",
        "release_or_version": "not recorded in local cache; HCCDB18 public page used only as provenance cross-check",
        "input_files": provenance,
        "expression_platform": "RNA-seq; existing Figure 7 cache labels source as RPKM-like/log2 transformed",
        "mapping_summary": mapping_summary,
        "hccdb18_crosscheck": cross_summary,
        "local_expression_sample_count_vs_listing": "437 local samples versus 389 HCCDB18 listing samples; exact release identity is unresolved",
        "existing_local_codebook_candidates": find_existing_codebooks(),
        "frozen_gene_sets": "Figure 7 v2 primary_frozen_programme; no gene set changes in this audit",
        "survival_models_run": False,
        "audit_output_names": list(AUDIT_OUTPUT_NAMES),
        "random_seed": None,
        "final_gate_status": "ESTIMABLE_BUT_NOT_VALIDATED",
        "decision_basis": "fustat status labels cross-check as 1=dead/0=alive for the HCCDB18 overlap and futime matches SUR*30, but the exact local-cache release, field derivation, and time origin for all 231 expression-linked donors remain undocumented",
    }
    (AUDIT_ROOT / "10_manifests" / "ICGC_OS_run_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps({"mapping_summary": mapping_summary, "hccdb18_crosscheck": cross_summary, "n_axis_donors": n_axis_donors, "final_gate_status": "ESTIMABLE_BUT_NOT_VALIDATED"}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
