from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
TUMOR_SOURCE_CLASSES = {"tumor", "pvtt_tumor", "metastatic_tumor_lymphnode", "unknown_hcc_dataset"}
NON_HCC_SOURCE_CLASSES = {"normal_adjacent", "non_hcc_liver"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Merge CopyKAT predictions into module 3 malignant HCC calls.")
    parser.add_argument("--manifest", type=Path, default=ROOT / "metadata/copykat_module3/copykat_input_manifest.tsv")
    parser.add_argument("--run-status", type=Path, default=ROOT / "metadata/copykat_module3/copykat_run_status.tsv")
    parser.add_argument("--base-calls", type=Path, default=ROOT / "metadata/malignant/malignant_hcc_calls_by_cell.tsv.gz")
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "metadata/malignant")
    parser.add_argument("--copykat-metadata-dir", type=Path, default=ROOT / "metadata/copykat_module3")
    return parser.parse_args()


def resolve_path(value: object) -> Path:
    path = Path(str(value))
    if path.exists():
        return path
    text = str(value).replace("\\", "/")
    for anchor in ("data/processed/", "metadata/"):
        idx = text.find(anchor)
        if idx >= 0:
            candidate = ROOT / Path(text[idx:])
            if candidate.exists():
                return candidate
    return path


def read_table(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", encoding="utf-8", **kwargs)


def prediction_path_for_run(row: pd.Series, status_by_run: dict[str, pd.Series]) -> Path | None:
    run_id = str(row["run_id"])
    status = status_by_run.get(run_id)
    if status is not None and str(status.get("prediction_path", "")).strip():
        path = resolve_path(status["prediction_path"])
        if path.exists():
            return path
    run_dir = resolve_path(row["run_dir"])
    direct = run_dir / f"{run_id}_copykat_prediction.txt"
    if direct.exists():
        return direct
    matches = sorted(run_dir.glob("*copykat_prediction.txt"))
    return matches[0] if matches else None


def load_copykat_calls(manifest: pd.DataFrame, status: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    status_by_run = {str(row["run_id"]): row for _, row in status.iterrows()} if not status.empty else {}
    rows = []
    skipped = []
    for _, run in manifest.iterrows():
        run_id = str(run["run_id"])
        run_dir = resolve_path(run["run_dir"])
        map_path = run_dir / "cell_map.tsv"
        if not map_path.exists():
            skipped.append({"run_id": run_id, "reason": "missing_cell_map"})
            continue
        cell_map = read_table(map_path)
        pred_path = prediction_path_for_run(run, status_by_run)
        if pred_path is None:
            pred = pd.DataFrame(columns=["cell_key", "copykat_pred"])
            skipped.append({"run_id": run_id, "reason": "missing_prediction"})
        else:
            pred = read_table(pred_path)
            pred = pred.rename(columns={"cell.names": "cell_key", "copykat.pred": "copykat_pred"})
            pred = pred[["cell_key", "copykat_pred"]].drop_duplicates("cell_key")

        merged = cell_map.merge(pred, on="cell_key", how="left")
        candidates = merged.loc[merged["cnv_role"].eq("candidate")].copy()
        candidates["copykat_run_id"] = run_id
        candidates["copykat_prediction_path"] = "" if pred_path is None else str(pred_path.resolve())
        candidates["copykat_pred"] = candidates["copykat_pred"].fillna("not_called_by_copykat")
        rows.append(
            candidates[
                [
                    "cell_id",
                    "dataset",
                    "study_sample",
                    "cnv_sample",
                    "sample_source_class",
                    "cell_key",
                    "copykat_run_id",
                    "copykat_pred",
                    "copykat_prediction_path",
                ]
            ]
        )
    calls = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["cell_id"])
    return calls, pd.DataFrame(skipped)


def copykat_status(pred: object) -> str:
    text = str(pred).strip().lower()
    if text == "" or text == "nan":
        return "not_run"
    if "aneuploid" in text:
        return "aneuploid"
    if "diploid" in text:
        return "diploid"
    if "not_called" in text:
        return "not_called_by_copykat"
    return text.replace(" ", "_")


def integrate_calls(base: pd.DataFrame, copykat_calls: pd.DataFrame) -> pd.DataFrame:
    out = base.copy()
    out = out.rename(columns={"malignant_hcc_call": "malignant_hcc_call_proxy"})
    out["copykat_status_proxy_script"] = out.get("copykat_status", pd.Series(index=out.index, dtype=object))
    out = out.drop(columns=["copykat_status"], errors="ignore")

    keep = ["cell_id", "cell_key", "copykat_run_id", "copykat_pred", "copykat_prediction_path"]
    ck = copykat_calls[keep].drop_duplicates("cell_id") if not copykat_calls.empty else pd.DataFrame(columns=keep)
    out = out.merge(ck, on="cell_id", how="left")
    out["copykat_pred"] = out["copykat_pred"].fillna("not_run")
    out["copykat_status"] = out["copykat_pred"].map(copykat_status)

    marker_high = (
        pd.to_numeric(out["hcc_malignant_associated_score_z"], errors="coerce").ge(0.8)
        | pd.to_numeric(out["hcc_malignant_associated_mean_log1p_cpm"], errors="coerce").ge(3.5)
        | out["hepatocyte_state_label"].astype(str).str.contains("malignant_hepatocyte_candidate", regex=False)
    )
    prolif_high = (
        pd.to_numeric(out["proliferation_score_z"], errors="coerce").ge(0.8)
        | out["hepatocyte_state_label"].astype(str).str.contains("proliferating", regex=False)
    )
    tumor_source = out["sample_source_class"].isin(TUMOR_SOURCE_CLASSES)
    non_hcc_source = out["sample_source_class"].isin(NON_HCC_SOURCE_CLASSES)
    copykat_aneuploid = out["copykat_status"].eq("aneuploid")
    copykat_diploid_or_no_aneuploid = ~copykat_aneuploid & ~out["copykat_status"].eq("not_run")
    copykat_missing = out["copykat_status"].isin(["not_run", "not_called_by_copykat"])

    out["malignant_hcc_call_copykat"] = "non_malignant_or_unresolved"
    out.loc[copykat_aneuploid & tumor_source & marker_high, "malignant_hcc_call_copykat"] = "malignant_hcc_high_conf"
    out.loc[copykat_aneuploid & tumor_source & ~marker_high, "malignant_hcc_call_copykat"] = "malignant_hcc_cnv_support"
    out.loc[
        copykat_diploid_or_no_aneuploid & tumor_source & marker_high & prolif_high,
        "malignant_hcc_call_copykat",
    ] = "malignant_hcc_marker_proliferation_needs_cnv_review"
    out.loc[non_hcc_source & ~copykat_aneuploid, "malignant_hcc_call_copykat"] = "not_malignant_source_or_cnv"
    out.loc[copykat_missing, "malignant_hcc_call_copykat"] = "copykat_not_available"

    official_available = ~copykat_missing
    out["malignant_hcc_call"] = out["malignant_hcc_call_proxy"]
    out.loc[official_available, "malignant_hcc_call"] = out.loc[official_available, "malignant_hcc_call_copykat"]
    out["malignant_hcc_cnv_method"] = "copykat"
    out.loc[~official_available, "malignant_hcc_cnv_method"] = "proxy_retained_no_copykat_call"
    out["malignant_hcc_evidence_copykat"] = (
        "copykat="
        + out["copykat_status"].fillna("NA").astype(str)
        + "; source="
        + out["sample_source_class"].fillna("NA").astype(str)
        + "; hcc_score_z="
        + pd.to_numeric(out["hcc_malignant_associated_score_z"], errors="coerce").round(3).astype(str)
        + "; prolif_score_z="
        + pd.to_numeric(out["proliferation_score_z"], errors="coerce").round(3).astype(str)
        + "; state="
        + out["hepatocyte_state_label"].astype(str)
    )
    return out


def main() -> int:
    args = parse_args()
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    args.copykat_metadata_dir.mkdir(parents=True, exist_ok=True)
    manifest = read_table(args.manifest)
    status = read_table(args.run_status) if args.run_status.exists() else pd.DataFrame()
    base = read_table(args.base_calls)

    copykat_calls, skipped = load_copykat_calls(manifest, status)
    if not copykat_calls.empty:
        copykat_calls["copykat_status"] = copykat_calls["copykat_pred"].map(copykat_status)
    final = integrate_calls(base, copykat_calls)

    calls_path = args.metadata_dir / "copykat_calls_by_cell.tsv.gz"
    calls_sample_path = args.metadata_dir / "copykat_calls_by_sample.tsv"
    final_path = args.metadata_dir / "malignant_hcc_calls_by_cell.copykat.tsv.gz"
    final_sample_path = args.metadata_dir / "malignant_hcc_calls_by_sample.copykat.tsv"
    skipped_path = args.copykat_metadata_dir / "copykat_merge_skipped_runs.tsv"
    report_path = args.metadata_dir / "malignant_hcc_module3_copykat_report.json"

    copykat_calls.to_csv(calls_path, sep="\t", index=False, compression="gzip")
    (
        copykat_calls.groupby(["dataset", "study_sample", "cnv_sample", "sample_source_class", "copykat_status"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["dataset", "cnv_sample", "copykat_status"])
        .to_csv(calls_sample_path, sep="\t", index=False)
    )
    skipped.to_csv(skipped_path, sep="\t", index=False)
    final.to_csv(final_path, sep="\t", index=False, compression="gzip")
    (
        final.groupby(["dataset", "study_sample", "cnv_sample", "sample_source_class", "malignant_hcc_call"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["dataset", "cnv_sample", "malignant_hcc_call"])
        .to_csv(final_sample_path, sep="\t", index=False)
    )

    report = {
        "method": "copykat_official_with_proxy_retained_for_missing_copykat_calls",
        "n_manifest_runs": int(manifest.shape[0]),
        "n_copykat_candidate_rows": int(copykat_calls.shape[0]),
        "n_base_candidate_rows": int(base.shape[0]),
        "n_copykat_predicted_rows": int(copykat_calls.loc[~copykat_calls["copykat_status"].eq("not_called_by_copykat")].shape[0]) if not copykat_calls.empty else 0,
        "copykat_status_counts": copykat_calls["copykat_status"].value_counts(dropna=False).to_dict() if not copykat_calls.empty else {},
        "malignant_hcc_call_copykat_counts": final["malignant_hcc_call_copykat"].value_counts(dropna=False).to_dict(),
        "malignant_hcc_call_recommended_counts": final["malignant_hcc_call"].value_counts(dropna=False).to_dict(),
        "outputs": {
            "copykat_calls_by_cell": str(calls_path.resolve()),
            "copykat_calls_by_sample": str(calls_sample_path.resolve()),
            "malignant_hcc_calls_by_cell_copykat": str(final_path.resolve()),
            "malignant_hcc_calls_by_sample_copykat": str(final_sample_path.resolve()),
            "copykat_merge_skipped_runs": str(skipped_path.resolve()),
        },
    }
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
