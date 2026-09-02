from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path


# This builder lives at <project_root>/21-SCI.../过程记录/; the project root
# is therefore two parents above the script file.
ROOT = Path(__file__).resolve().parents[2]
STAGE = ROOT / "21-SCI生信研究Methods撰写器"
OUT = STAGE / "输出"
INPUT = STAGE / "输入"
LOG = STAGE / "过程记录"
QA = STAGE / "质量核查"

WATERMARK = (
    "> 版权声明：本文件由杨师兄原创“研究型论文 Skill 系统”生成。  \n"
    "> 未经书面授权，禁止复制、传播、改编、转售、商用或用于第三方交付。  \n"
    "> 授权请联系杨师兄。\n"
)


def now_local() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rel(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return str(path).replace("\\", "/")


def write_text(path: Path, text: str, watermark: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if watermark and not text.startswith(WATERMARK):
        text = WATERMARK + "\n" + text.lstrip("\n")
    path.write_text(text, encoding="utf-8", newline="\n")


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str], delimiter: str = ",") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=delimiter, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def source_inventory() -> list[Path]:
    rel_paths = [
        "08-SCI数据集纳入排除与下载记录器/输出/下载记录与版本表.csv",
        "08-SCI数据集纳入排除与下载记录器/输出/样本分组与元数据表.csv",
        "08-SCI数据集纳入排除与下载记录器/输出/数据集纳排标准表.csv",
        "08-SCI数据集纳入排除与下载记录器/输出/dataset_intake_audit.json",
        "metadata/public_data_accession_version_license.md",
        "metadata/public_dataset_download_manifest.after.tsv",
        "10-SCI生信分析流程设计器/输出/生信分析流程图.md",
        "10-SCI生信分析流程设计器/输出/分析模块与软件包表.csv",
        "10-SCI生信分析流程设计器/输出/关键参数预设表.md",
        "10-SCI生信分析流程设计器/输出/analysis_workflow_audit.json",
        "11-SCI数据预处理与批次效应处理器/输出/预处理方案.md",
        "11-SCI数据预处理与批次效应处理器/输出/预处理输入对象清单.csv",
        "11-SCI数据预处理与批次效应处理器/输出/标准化与输入对象记录.csv",
        "11-SCI数据预处理与批次效应处理器/输出/批次效应处理记录表.csv",
        "11-SCI数据预处理与批次效应处理器/输出/preprocessing_input_audit.json",
        "11-SCI数据预处理与批次效应处理器/输出/sessionInfo.txt",
        "12-SCI差异分析与功能富集分析器/输出/差异分析方案.md",
        "12-SCI差异分析与功能富集分析器/输出/富集分析方案.md",
        "12-SCI差异分析与功能富集分析器/输出/差异分析对比设计表.csv",
        "12-SCI差异分析与功能富集分析器/输出/deg_enrichment_input_audit.json",
        "13-SCI网络分析与关键基因筛选器/输出/网络分析方案.md",
        "13-SCI网络分析与关键基因筛选器/输出/关键基因筛选逻辑表.csv",
        "13-SCI网络分析与关键基因筛选器/输出/三轴候选汇总.csv",
        "13-SCI网络分析与关键基因筛选器/输出/network_analysis_audit.json",
        "14-SCI机器学习建模与验证器/输出/机器学习建模方案.md",
        "14-SCI机器学习建模与验证器/输出/machine_learning_audit.json",
        "15-SCI外部验证与实验验证规划器/输出/外部验证方案.md",
        "15-SCI外部验证与实验验证规划器/输出/验证结果汇总表.csv",
        "15-SCI外部验证与实验验证规划器/输出/external_validation_audit.json",
        "16-SCI生信结果图表制作器/输出/图件契约记录.md",
        "16-SCI生信结果图表制作器/输出/图注与统计标注草稿.md",
        "16-SCI生信结果图表制作器/输出/图表数据来源索引.csv",
        "16-SCI生信结果图表制作器/输出/figure_contract_audit.json",
        "17-SCI生信代码与可重复性核查器/输出/软件版本与环境记录.md",
        "17-SCI生信代码与可重复性核查器/输出/代码可重复性核查表.csv",
        "17-SCI生信代码与可重复性核查器/输出/图件状态风险表.csv",
        "17-SCI生信代码与可重复性核查器/输出/code_reproducibility_audit.json",
        "18-SCI生信结果主线提炼器/输出/生信结果主线.md",
        "18-SCI生信结果主线提炼器/输出/证据等级与禁用表述.md",
        "18-SCI生信结果主线提炼器/输出/Results段落蓝图.md",
        "19-SCI生信研究Results撰写器/输出/results_draft_audit.json",
        "19-SCI生信研究Results撰写器/质量核查/质量核查表.csv",
        "19-SCI生信研究Results撰写器/reproducibility_closure/19_stage_gate_reassessment.json",
        "19-SCI生信研究Results撰写器/reproducibility_closure/软件环境与精确版本审计.md",
        "20-SCI生信研究Discussion撰写器/输出/Discussion草稿.md",
        "20-SCI生信研究Discussion撰写器/过程记录/stage20_run_record.json",
        "20-SCI生信研究Discussion撰写器/质量核查/质量核查表.csv",
        "patient_level_validation_cohort_rescue/stage19_final_closure_adjudication_v1/07_stage19_gate/STAGE19_FINAL_GATE.json",
        "patient_level_validation_cohort_rescue/stage19_final_closure_adjudication_v1/08_stage20_handoff/STAGE20_HANDOFF.json",
        "patient_level_validation_cohort_rescue/stage19_final_closure_adjudication_v1/09_qa/QA_SUMMARY.json",
        "patient_level_validation_cohort_rescue/stage19_final_closure_adjudication_v1/10_reports/STAGE19_FINAL_CLOSURE_ADJUDICATION_REPORT.md",
        "patient_level_validation_cohort_rescue/stage19b_gse189175_validation_v1/11_reports/STAGE19B_GSE189175_FINAL_REPORT.md",
        "patient_level_validation_cohort_rescue/stage19b_gse189175_validation_v1/05_pseudobulk/pseudobulk_build_manifest.json",
        "patient_level_validation_cohort_rescue/stage19b_gse189175_validation_v1/05_pseudobulk/pseudobulk_qc_summary.json",
        "patient_level_validation_cohort_rescue/stage19b_gse189175_validation_v1/07_three_axis_validation/FROZEN_AXIS_DEFINITION.json",
        "patient_level_validation_cohort_rescue/stage19b_gse189175_validation_v1/07_three_axis_validation/THREE_AXIS_ADJUDICATION_RULES_PREREGISTERED.json",
        "patient_level_validation_cohort_rescue/stage19b_gse189175_validation_v1/10_logs/run_paired_edgeR.R",
        "patient_level_validation_cohort_rescue/stage19b_gse189175_validation_v1/10_logs/R_sessionInfo_paired_edgeR.txt",
        "patient_level_validation_cohort_rescue/GSE326201_tier3_validation/GSE326201_tier3_validation_report.md",
        "patient_level_validation_cohort_rescue/GSE326201_tier3_validation/GSE326201_tier3_validation_audit.json",
        "patient_level_validation_cohort_rescue/GSE326201_tier3_validation/GSE326201_three_axis_validation_audit.json",
        "patient_level_validation_cohort_rescue/GSE326201_tier3_validation/edgeR_results/GSE326201_edgeR_audit.json",
        "patient_level_validation_cohort_rescue/GSE326201_tier3_validation/software_environment_manifest.json",
        "patient_level_validation_cohort_rescue/GSE282701_rescue/GSE282701_tier3_validation_report.md",
        "patient_level_validation_cohort_rescue/GSE282701_rescue/GSE282701_tier3_validation_audit.json",
        "patient_level_validation_cohort_rescue/stage20_5_public_data_compliance_audit_v1/08_compliance_gate/PUBLIC_DATA_COMPLIANCE_GATE.json",
        "patient_level_validation_cohort_rescue/stage20_5_public_data_compliance_audit_v1/09_qa/QA_SUMMARY.json",
        "patient_level_validation_cohort_rescue/stage20_5_public_data_compliance_audit_v1/10_reports/COMPLIANCE_DECISION_SUMMARY.md",
        "patient_level_validation_cohort_rescue/stage20_5_public_data_compliance_audit_v1/07_manuscript_ready/DATA_AVAILABILITY_FINAL_DRAFT.md",
        "patient_level_validation_cohort_rescue/stage20_5_public_data_compliance_audit_v1/07_manuscript_ready/ETHICS_SECONDARY_ANALYSIS_WORDING_CANDIDATES.md",
        "reports/module6_3b_canonical_scenic_final_report.md",
        "metadata/driver/scenic_module6_3b/driver_module6_3b_final_report.json",
        "metadata/driver/scenic_module6_3b/driver_module6_3b_grnboost2_report.json",
        "metadata/driver/scenic_module6_3b/driver_module6_3b_ctx_report.json",
        "metadata/driver/scenic_module6_3b/driver_module6_3b_aucell_report.json",
        "metadata/driver/driver_module6_2_cellrank_report.json",
        "metadata/driver/celloracle_module6_7_grn_report.json",
        "metadata/driver/celloracle_module6_8_perturbation_report.json",
        "metadata/driver/module9_1_report.json",
        "metadata/driver/module9_2_report.json",
        "metadata/driver/module9_3_report.json",
        "metadata/driver/module9_4_report.json",
        "metadata/driver/figure5_temporal_positioning/figure5_total_report.json",
        "metadata/driver/figure6_directional_network/figure6f_model_report.json",
        "metadata/driver/figure6_directional_network/figure6_r_package_versions.tsv",
        "reports/figure2_hnf4a_b_to_f_report.md",
        "reports/figure3_egr1_a_to_f_report.md",
        "reports/figure5_temporal_positioning_report.md",
        "reports/figure6_directional_network_report.md",
        "reports/figure7_external_bulk_clinical_validation_v2_report.md",
        "reports/figure8_transcriptomic_reversal_v2_mainfigure/figure8_v2_mainfigure_reanalysis_report.md",
        "reports/figure8_transcriptomic_reversal_v2_mainfigure/figure8_v2_integrated_score_definition.md",
        "metadata/driver/figure8_transcriptomic_reversal_v2_mainfigure/figure8_v2_validation_report.json",
        "metadata/reproducibility/software_environment_manifest.json",
    ]
    return [ROOT / p for p in rel_paths]


def build_input_manifest() -> list[dict[str, object]]:
    rows = []
    for path in source_inventory():
        exists = path.is_file()
        rows.append(
            {
                "relative_path": rel(path),
                "exists": str(exists).upper(),
                "file_size_bytes": path.stat().st_size if exists else "",
                "modified_time": datetime.fromtimestamp(path.stat().st_mtime).astimezone().isoformat(timespec="seconds") if exists else "",
                "sha256": sha256(path) if exists and path.stat().st_size <= 100 * 1024 * 1024 else ("LARGE_FILE_HASH_IN_UPSTREAM_MANIFEST" if exists else ""),
                "role": "Stage21_frozen_input_or_audit_record",
            }
        )
    return rows


def build_dependency_outputs() -> tuple[list[dict[str, object]], list[dict[str, object]], dict]:
    manifest_path = ROOT / "metadata/reproducibility/software_environment_manifest.json"
    manifest = read_json(manifest_path)
    r_rows = []
    for row in manifest.get("r_environment", {}).get("package_rows", []):
        r_rows.append(
            {
                "package": row.get("package", ""),
                "call_types": row.get("call_types", ""),
                "file_count": row.get("file_count", ""),
                "source_files": row.get("source_files", ""),
                "package_class": row.get("package_class", ""),
                "bioconductor_signal": row.get("bioconductor_signal", ""),
                "current_status": row.get("current_status", ""),
                "current_version": row.get("current_version", ""),
                "current_lib_path": row.get("current_lib_path", ""),
                "historical_version_status": row.get("historical_version_status", ""),
                "historical_versions_seen": row.get("historical_versions_seen", ""),
                "historical_evidence_files": row.get("historical_evidence_files", ""),
                "closure_status": row.get("closure_status", ""),
            }
        )
    py_rows = []
    for row in manifest.get("python_environments", {}).get("dependency_rows", []):
        py_rows.append(
            {
                "module": row.get("module", ""),
                "distribution": row.get("distribution", ""),
                "source_file_count": row.get("source_file_count", ""),
                "source_files": row.get("source_files", ""),
                "roles": row.get("roles", ""),
                "active_python_version": row.get("active_python_version", ""),
                "project_venv_version": row.get("project_venv_version", ""),
                "active_status": row.get("active_status", ""),
                "project_venv_status": row.get("project_venv_status", ""),
                "historical_version_status": row.get("historical_version_status", ""),
                "historical_versions_seen": row.get("historical_versions_seen", ""),
                "historical_evidence_files": row.get("historical_evidence_files", ""),
                "closure_status": row.get("closure_status", ""),
            }
        )
    return r_rows, py_rows, manifest


def methods_parameter_rows() -> list[dict[str, object]]:
    return [
        {"module_or_figure": "Study design", "analysis_role": "evidence architecture", "parameter_or_item": "architecture", "recorded_value": "cross-sectional, partial-order three-axis framework; no strict biological cascade", "source_path": "10-SCI生信分析流程设计器/输出/生信分析流程图.md;18-SCI生信结果主线提炼器/输出/证据等级与禁用表述.md", "status": "CLOSED", "provenance_status": "verified_frozen_record", "notes": "Computational associations are not direct causality."},
        {"module_or_figure": "Primary data", "analysis_role": "integrated discovery", "parameter_or_item": "raw-count datasets", "recorded_value": "GSE149614; GSE151530; GSE174748; GSE185477; GSE202379; GSE212046", "source_path": "11-SCI数据预处理与批次效应处理器/输出/preprocessing_input_audit.json", "status": "CLOSED", "provenance_status": "verified_frozen_record", "notes": "17 counts objects; 420435 included cells; dataset roles remain distinct."},
        {"module_or_figure": "Reference", "analysis_role": "annotation/reference", "parameter_or_item": "HCC atlas", "recorded_value": "figshare article 22332568, version 1, DOI 10.6084/m9.figshare.22332568.v1; CC BY 4.0 per audit", "source_path": "08-SCI数据集纳入排除与下载记录器/输出/数据集纳排标准表.csv;metadata/public_data_accession_version_license.md", "status": "CLOSED_WITH_LIMITATION", "provenance_status": "verified_source_record", "notes": "Processed/non-raw object excluded from scVI count input; use for reference/annotation only."},
        {"module_or_figure": "Data provenance", "analysis_role": "repository version", "parameter_or_item": "GEO/TCGA/ICGC access dates and release", "recorded_value": "GEO local download date is session-level 2026-05-31; exact per-file download version not recorded; TCGA/ICGC cache release not recorded", "source_path": "08-SCI数据集纳入排除与下载记录器/输出/下载记录与版本表.csv;metadata/public_data_accession_version_license.md", "status": "AUTHOR ACTION REQUIRED", "provenance_status": "partial_provenance", "notes": "Do not replace missing historical release with current release."},
        {"module_or_figure": "QC", "analysis_role": "cell/gene filtering", "parameter_or_item": "QC thresholds", "recorded_value": "min_genes=200; min_counts=500; max_mito_pct=25; min_cells_per_gene=3; upper_quantile=0.995", "source_path": "10-SCI生信分析流程设计器/输出/关键参数预设表.md;11-SCI数据预处理与批次效应处理器/输出/预处理方案.md", "status": "CLOSED", "provenance_status": "verified_script_and_record", "notes": "Low-cell review threshold kept_cells<3000; no automatic equivalence to exclusion."},
        {"module_or_figure": "QC", "analysis_role": "input eligibility", "parameter_or_item": "scVI input exclusions", "recorded_value": "kept_cells<1000; kept_genes<5000; non-integer rate>0.001; numeric gene-name rate>0.5", "source_path": "10-SCI生信分析流程设计器/输出/关键参数预设表.md;11-SCI数据预处理与批次效应处理器/输出/预处理输入对象清单.csv", "status": "CLOSED", "provenance_status": "verified_frozen_record", "notes": "HCC atlas and normalized/non-count objects are excluded from scVI input."},
        {"module_or_figure": "scVI", "analysis_role": "batch integration", "parameter_or_item": "model and training", "recorded_value": "n_top_genes=2000; n_latent=30; n_hidden=128; n_layers=2; batch_key=dataset; max_epochs=10; batch_size=1024; seed=20260601; GPU", "source_path": "metadata/scvi/scvi_training_report.json;10-SCI生信分析流程设计器/输出/关键参数预设表.md", "status": "CLOSED", "provenance_status": "verified_run_record", "notes": "420435 cells x 2000 variables; 11010 common genes; scvi-tools=1.3.3; torch=2.4.1+cu121."},
        {"module_or_figure": "scanVI", "analysis_role": "semi-supervised labels", "parameter_or_item": "model and prediction", "recorded_value": "max_epochs=20; batch_size=1024; n_samples_per_label=1000; seed=20260601; high-confidence seeds plus malignant-call constraint", "source_path": "metadata/scanvi/scanvi_unified_module4_report.json", "status": "CLOSED", "provenance_status": "verified_run_record", "notes": "scvi-tools=1.3.3; torch=2.4.1+cu121; predicted labels are constrained for strict malignant calls."},
        {"module_or_figure": "Embedding", "analysis_role": "neighbors/UMAP/Leiden", "parameter_or_item": "low-dimensional graph", "recorded_value": "use_rep=X_scVI; n_neighbors=30; min_dist=0.3; Leiden resolution=1.0; seed=20260601", "source_path": "metadata/scvi/scvi_neighbors_umap_leiden_report.json", "status": "CLOSED", "provenance_status": "verified_run_record", "notes": "62 Leiden clusters."},
        {"module_or_figure": "Cell annotation", "analysis_role": "major cell types", "parameter_or_item": "annotation models", "recorded_value": "CellTypist Healthy_Human_Liver.pkl; CellTypist=1.7.1; Scanpy=1.11.5; SingleR=2.14.0; celldex=1.22.0; manual 57-marker review", "source_path": "metadata/celltype/celltypist_major_report.json;metadata/celltype/manual_marker_module1_report.json;metadata/celltype/singler_cluster_report.json", "status": "CLOSED_WITH_REVIEW", "provenance_status": "verified_run_record", "notes": "Manual marker thresholds min_panel_z=0.35; min_panel_pct=0.05; external conflict confidence=0.6."},
        {"module_or_figure": "Doublet/cell cycle", "analysis_role": "QC diagnostic", "parameter_or_item": "simulated-doublet kNN", "recorded_value": "expected doublet rate=0.05; global predicted rate=0.0559801; cluster 16 flagged; cycling clusters retained as flags", "source_path": "metadata/scvi/scvi_doublet_cell_cycle_report.json", "status": "CLOSED_WITH_REVIEW", "provenance_status": "verified_run_record", "notes": "Diagnostic labels are not a separate biological endpoint."},
        {"module_or_figure": "CNV/malignant", "analysis_role": "CopyKAT", "parameter_or_item": "CopyKAT call", "recorded_value": "id.type=S; cell.line=no; ngene.chr=5; min.gene.per.cell=200; LOW.DR=0.05; UP.DR=0.1; win.size=25; KS.cut=0.1; distance=euclidean; genome=hg20; n.cores=1", "source_path": "scripts/run_copykat_module3.R;metadata/malignant/malignant_hcc_module3_copykat_report.json", "status": "CLOSED_WITH_VERSION_LIMITATION", "provenance_status": "verified_script_parameters", "notes": "CopyKAT package/runtime exact historical version is not recoverable; proxy calls are retained separately."},
        {"module_or_figure": "CNV/malignant", "analysis_role": "expression CNV proxy", "parameter_or_item": "cytoband-binned proxy", "recorded_value": "gene_bin_size=100; cnv_z_threshold=3.0; cnv_high_bin_fraction=0.08; gene map=22411", "source_path": "scripts/malignant_hcc_module3.py;metadata/malignant/malignant_hcc_module3_report.proxy.json", "status": "CLOSED_WITH_BOUNDARY", "provenance_status": "verified_script_and_record", "notes": "Proxy is not DNA-level CNV confirmation."},
        {"module_or_figure": "CNV/malignant", "analysis_role": "marker integration", "parameter_or_item": "malignant evidence", "recorded_value": "marker_high if score_z>=0.8 or mean_log1p_CPM>=3.5 or malignant state label; proliferation_high if score_z>=0.8 or proliferating label", "source_path": "scripts/malignant_hcc_module3.py", "status": "CLOSED_WITH_BOUNDARY", "provenance_status": "verified_script_parameters", "notes": "Tumour-source classes and CNV/proliferation review categories are retained separately."},
        {"module_or_figure": "Trajectory", "analysis_role": "trajectory object", "parameter_or_item": "object and orientation", "recorded_value": "283498 cells x 2000 variables; use_rep=X_scANVI; root=reference hepatocyte; end=CNV-supported malignant; main_strict plus include-review sensitivity", "source_path": "metadata/trajectory/trajectory_module5_1_report.json;metadata/trajectory/trajectory_module5_2_report.json", "status": "CLOSED_WITH_BOUNDARY", "provenance_status": "verified_run_record", "notes": "Stage labels are derived computational/source classes, not TNM/BCLC stages."},
        {"module_or_figure": "Trajectory", "analysis_role": "trajectory methods", "parameter_or_item": "Monocle3/Slingshot", "recorded_value": "Monocle3; Slingshot on X_scANVI and hepatocyte PCA; versions not present in frozen run records", "source_path": "metadata/trajectory/trajectory_module5_7_report.json;reports/figure5_temporal_positioning_report.md", "status": "AUTHOR ACTION REQUIRED", "provenance_status": "historical_exact_version_not_recoverable", "notes": "Do not insert a current package version."},
        {"module_or_figure": "CellRank", "analysis_role": "fate probability and drivers", "parameter_or_item": "main_strict", "recorded_value": "5000 cells; 2000 genes; terminal states late_non_cnv and CNV-supported malignant; fate rows audited to approximately 1", "source_path": "metadata/driver/driver_module6_2_cellrank_report.json;10-SCI生信分析流程设计器/输出/关键参数预设表.md", "status": "CLOSED_WITH_BOUNDARY", "provenance_status": "verified_run_record", "notes": "CellRank=2.0.7; Scanpy=1.11.5; anndata=0.12.16; driver associations are not causality."},
        {"module_or_figure": "Canonical SCENIC", "analysis_role": "GRNBoost2", "parameter_or_item": "formal GRN", "recorded_value": "9512 cells; 11923 genes; 1767 TFs; 1131800 edges; seed=777; report num_workers=8", "source_path": "metadata/driver/scenic_module6_3b/driver_module6_3b_grnboost2_report.json;reports/module6_3b_canonical_scenic_final_report.md", "status": "CLOSED_WITH_SINGLE_SEED", "provenance_status": "verified_run_record", "notes": "Environment snapshot separately lists formal_grn_workers=32; reconcile worker count before final Methods."},
        {"module_or_figure": "Canonical SCENIC", "analysis_role": "cisTarget/ctx", "parameter_or_item": "motif pruning", "recorded_value": "hg38 RefSeq80; mc_v10_clust 10-kb and proximal resources; rank_threshold=5000; auc_threshold=0.05; NES>=3.0; top_n_targets=50; top_n_regulators=5,10,50; min_genes=20", "source_path": "metadata/driver/scenic_module6_3b/driver_module6_3b_ctx_report.json;reports/module6_3b_canonical_scenic_final_report.md", "status": "CLOSED_WITH_SINGLE_SEED", "provenance_status": "verified_run_record", "notes": "422 regulons; 319 primary regulons; 103 small regulons; resource hashes are in upstream audit."},
        {"module_or_figure": "Canonical SCENIC", "analysis_role": "AUCell", "parameter_or_item": "regulon activity", "recorded_value": "9512 x 422 regulon matrix; rank_threshold=5000; auc_threshold=0.05; NES threshold=3.0; seed=777", "source_path": "metadata/driver/scenic_module6_3b/driver_module6_3b_aucell_report.json", "status": "CLOSED_WITH_SINGLE_SEED", "provenance_status": "verified_run_record", "notes": "CellRank association uses the 5000-cell fate subset; 4512 cells lack fate metadata by design."},
        {"module_or_figure": "CellOracle", "analysis_role": "state-specific GRN", "parameter_or_item": "GRN fitting", "recorded_value": "CellOracle=0.20.0; Python=3.10.20 in WSL2; alpha_links=10; alpha_simulation=1; bagging_number=20; n_pca_dims=30; k=30; p_threshold=0.001; threshold_number=10000; n_jobs=8", "source_path": "metadata/driver/celloracle_module6_7_grn_report.json", "status": "CLOSED_WITH_BOUNDARY", "provenance_status": "verified_run_record", "notes": "Historical execution environment is separate from current Windows Python; no local re-execution is claimed."},
        {"module_or_figure": "CellOracle", "analysis_role": "virtual perturbation", "parameter_or_item": "KO simulation", "recorded_value": "15 TFs; knockout; n_propagation=3; clip_delta_x=true; n_neighbors=100; sampled_fraction=0.3; sigma_corr=0.05; grid_steps=20; grid_neighbors=50; top_genes_per_state=50; seed=15071990", "source_path": "metadata/driver/celloracle_module6_8_perturbation_report.json;10-SCI生信分析流程设计器/输出/关键参数预设表.md", "status": "CLOSED_WITH_BOUNDARY", "provenance_status": "verified_run_record", "notes": "KO only; controls n=5; no restore/OE; AP-1 summaries are member-wise summaries, not combined suppression."},
        {"module_or_figure": "scTenifoldKnk", "analysis_role": "primary Module7 network", "parameter_or_item": "historical primary settings", "recorded_value": "45 TF-KO runs; qc=false; nc_nNet=1; nc_nCells=100; nCores=1; seed=11; 3000-gene background", "source_path": "10-SCI生信分析流程设计器/输出/关键参数预设表.md;reports/sctenifoldknk_reproducibility_audit_v2/01_historical_audit.md", "status": "CLOSED_WITH_BOUNDARY", "provenance_status": "reported_run_record", "notes": "Unsigned manifold displacement; package exact historical version remains not fully recoverable."},
        {"module_or_figure": "scTenifoldKnk", "analysis_role": "v2 reproducibility audit", "parameter_or_item": "standardized selected TF runs", "recorded_value": "nc_nNet=10; nc_nCells=500; nc_nComp=3; ma_nDim=2; qc=false; qc_minCells=3; three seeds=15071990,15071991,15071992", "source_path": "reports/sctenifoldknk_reproducibility_audit_v2/02_rerun_decision.md;reports/sctenifoldknk_reproducibility_audit_v2/03_three_axis_validation_report.md", "status": "CLOSED_WITH_BOUNDARY", "provenance_status": "verified_audit_record", "notes": "HNF4A and SOX4 rerun; EGR1 retained from matching historical three-seed run; do not conflate with original Figure2/3 frozen panel counts."},
        {"module_or_figure": "Temporal hierarchy", "analysis_role": "relative positioning", "parameter_or_item": "bootstrap and landmarks", "recorded_value": "grid_n=201; baseline_end=0.10; search_start=0.10; slope_search_end=0.95; onset_min_run=5; crossing_min_run=3; derivative_fraction=0.10; min_effect_fraction=0.05; min_effect_baseline_sd=0.25; precedence_tolerance=0.005; bootstrap=1000; seed=20260805", "source_path": "10-SCI生信分析流程设计器/输出/关键参数预设表.md;reports/figure5_temporal_positioning_report.md", "status": "CLOSED_WITH_BOUNDARY", "provenance_status": "verified_frozen_record", "notes": "Pseudotime is relative; right edges are fades, not termination events."},
        {"module_or_figure": "Temporal hierarchy", "analysis_role": "GAM trends", "parameter_or_item": "smoothing", "recorded_value": "patient-pseudotime-bin pseudobulk; mgcv GAM with fixed dataset effects and patient random-effect smooth where estimable; 95% model-based intervals", "source_path": "reports/figure5_temporal_positioning_report.md", "status": "CLOSED_WITH_VERSION_LIMITATION", "provenance_status": "verified_method_record", "notes": "mgcv=1.9.4 in Figure5 report; exact Monocle3/Slingshot versions remain unresolved."},
        {"module_or_figure": "Network direction", "analysis_role": "CellOracle/scTenifold concordance", "parameter_or_item": "cross-method integration", "recorded_value": "8 signature sets; 48 asymmetry rows; 180 concordance rows; restore unavailable", "source_path": "metadata/driver/module9_2_report.json;reports/figure6_directional_network_report.md", "status": "CLOSED_WITH_BOUNDARY", "provenance_status": "verified_run_record", "notes": "Direction and rank evidence only; no genetic epistasis claim."},
        {"module_or_figure": "Statistical mediation", "analysis_role": "partial mediation", "parameter_or_item": "bootstrap", "recorded_value": "bootstrap=1000; seed=20260615; path models evaluated only where outcomes available", "source_path": "metadata/driver/module9_3_report.json", "status": "CLOSED_WITH_BOUNDARY", "provenance_status": "verified_run_record", "notes": "Partial computational mediation does not establish a causal chain."},
        {"module_or_figure": "Differential analysis", "analysis_role": "cell-level exploratory specification", "parameter_or_item": "FindMarkers/Wilcoxon", "recorded_value": "min.pct=0.10; abs(log2FC)>=0.25; BH-FDR<0.05", "source_path": "12-SCI差异分析与功能富集分析器/输出/差异分析方案.md;12-SCI差异分析与功能富集分析器/输出/deg_enrichment_input_audit.json", "status": "SPECIFIED_NOT_COMPLETE_AS_GLOBAL_TABLE", "provenance_status": "prespecified_record", "notes": "No complete traditional cell-level DEG table was found in the frozen inventory; do not describe as a completed global result."},
        {"module_or_figure": "Differential analysis", "analysis_role": "sample/patient pseudobulk", "parameter_or_item": "edgeR/limma", "recorded_value": "DGEList -> filterByExpr(y, design) -> calcNormFactors(TMM) -> estimateDisp(robust=TRUE) -> glmQLFit(robust=TRUE) -> glmQLFTest; BH; abs(log2FC)>=0.50; paired ~patient_id+condition or unpaired ~condition", "source_path": "19-SCI生信研究Results撰写器/sample_level_deg_closure/scripts/run_edgeR_pseudobulk.R;19-SCI生信研究Results撰写器/sample_level_deg_closure/sample_level_DEG_QA.md", "status": "CLOSED_WITH_LIMITATIONS", "provenance_status": "verified_script_and_run_record", "notes": "Cells/nuclei are not biological replicates; one QL-F-derived SE row remains NA in historical output."},
        {"module_or_figure": "GSE189175", "analysis_role": "matched patient support", "parameter_or_item": "fixed compartment and pseudobulk", "recorded_value": "58100 genes x 39995 cells; marker rule >=2 hepatocyte markers or >=2 epithelial markers; 33419 retained; 6 profiles from 3 matched patients; raw integer counts", "source_path": "patient_level_validation_cohort_rescue/stage19b_gse189175_validation_v1/11_reports/STAGE19B_GSE189175_FINAL_REPORT.md;.../05_pseudobulk/pseudobulk_build_manifest.json", "status": "CLOSED_WITH_LIMITATIONS", "provenance_status": "verified_frozen_record", "notes": "No clustering, InferCNV or malignant calling in this validation compartment."},
        {"module_or_figure": "GSE189175", "analysis_role": "paired edgeR", "parameter_or_item": "count model and QC", "recorded_value": "min compartment cells=100; min library=100000; min detected genes=1000; low-depth fraction=0.2; pair log2 library ratio review=1.5; design ~patient_id_original+condition; Tumor-Adjacent", "source_path": "patient_level_validation_cohort_rescue/stage19b_gse189175_validation_v1/05_pseudobulk/pseudobulk_qc_summary.json;.../10_logs/run_paired_edgeR.R", "status": "CLOSED_WITH_LIMITATIONS", "provenance_status": "verified_run_record", "notes": "R=4.3.1; edgeR=3.42.4; limma=3.56.2; filterByExpr retained 23995/58100 features; n=3 limits inference."},
        {"module_or_figure": "GSE326201", "analysis_role": "exploratory matched validation", "parameter_or_item": "CNV recall and pseudobulk", "recorded_value": "direct official mapping; InferCNV eligible call thresholds reference>=50 and malignant>=20; 3 eligible pairs of 10 patients; final Tier 1 exploratory", "source_path": "patient_level_validation_cohort_rescue/GSE326201_tier3_validation/GSE326201_tier3_validation_report.md;.../GSE326201_tier3_validation_audit.json", "status": "CLOSED_WITH_EXPLORATORY_BOUNDARY", "provenance_status": "verified_frozen_record", "notes": "Formal Tier 2+ validation remains unavailable."},
        {"module_or_figure": "GSE326201", "analysis_role": "InferCNV/edgeR runtime", "parameter_or_item": "software", "recorded_value": "InferCNV R=4.6.1, infercnv=1.28.0, Matrix=1.7.5; edgeR R=4.5.0, edgeR=4.8.2, limma=3.66.0", "source_path": "patient_level_validation_cohort_rescue/GSE326201_tier3_validation/software_environment_manifest.json;.../edgeR_results/GSE326201_edgeR_audit.json", "status": "CLOSED_WITH_LIMITATIONS", "provenance_status": "verified_rescue_record", "notes": "Current runtime is provenance only; historical full environment is not asserted."},
        {"module_or_figure": "Enrichment", "analysis_role": "ORA/GSEA", "parameter_or_item": "pathway tests", "recorded_value": "one-sided hypergeometric ORA for GO/KEGG/Reactome; matched background; network background=3000 genes where specified; BH within analysis family; preranked GSEA/GSVA only when recorded", "source_path": "12-SCI差异分析与功能富集分析器/输出/富集分析方案.md;12-SCI差异分析与功能富集分析器/输出/deg_enrichment_input_audit.json", "status": "CLOSED_WITH_BOUNDARY", "provenance_status": "prespecified_record", "notes": "Network displacement enrichment is separate from traditional DEG enrichment."},
        {"module_or_figure": "Figure 7", "analysis_role": "bulk expression recurrence", "parameter_or_item": "cohort-specific programme score", "recorded_value": "frozen Figure7 v2 signatures; patient-level mean aggregation; unsigned associated target programme score; Hedges g with 95% CI; exploratory random-effects k=2; seed=20260812; 1000 matched random signatures per cohort-axis", "source_path": "scripts/figure7_v2_core.R;reports/figure7_external_bulk_clinical_validation_v2_report.md;metadata/driver/figure7_external_validation_v2/figure7_v2_signature_manifest.tsv", "status": "CLOSED_WITH_LIMITATIONS", "provenance_status": "verified_frozen_record", "notes": "Signature derivation did not use bulk outcomes; expression recurrence is the main Figure7 candidate layer."},
        {"module_or_figure": "Figure 7", "analysis_role": "ICGC OS", "parameter_or_item": "exploratory Cox", "recorded_value": "continuous frozen score; HR per 1 within-ICGC tumour SD; 231 donor-linked records; 42 events; no cutoff search; PH and nonlinearity audit", "source_path": "Figure7_ICGC_OS_Audit;metadata/driver/figure7_external_validation_v2/figure7_v2_icgc_os_policy_update_20260828.md", "status": "ESTIMABLE_BUT_NOT_VALIDATED", "provenance_status": "verified_scope_with_manual_semantics", "notes": "Supplementary/Extended Data only; release/codebook and variable semantics remain limited."},
        {"module_or_figure": "Figure 8", "analysis_role": "transcriptomic reversal", "parameter_or_item": "continuous landmark-space score", "recorded_value": "978 landmarks; robust scale [-1,1]; weights 0.35/0.15/0.10/0.10/0.10/0.15/0.05; directional agreement shrink; unsupported coordinates=0", "source_path": "reports/figure8_transcriptomic_reversal_v2_mainfigure/figure8_v2_mainfigure_reanalysis_report.md;.../figure8_v2_integrated_score_definition.md", "status": "CLOSED_AS_EXTENDED_DATA_ONLY", "provenance_status": "verified_frozen_record", "notes": "839 reliable scores; no Tier A/B/C; no efficacy, safety or clinical actionability."},
        {"module_or_figure": "Figure 8", "analysis_role": "rescue candidate ranking", "parameter_or_item": "null and external resources", "recorded_value": "DrugReflector V3.5; 9597 checkpoints; 2000 matched null signatures; exact identifier hierarchy InChIKey>name>BRD; CMap/L1000FWD/CLUE; PRISM 23Q2 primary and 19Q4 secondary", "source_path": "reports/figure8_transcriptomic_reversal_v2_mainfigure/figure8_v2_mainfigure_reanalysis_report.md;metadata/driver/figure8_transcriptomic_reversal_v2_mainfigure/figure8_v2_validation_report.json", "status": "CLOSED_AS_EXTENDED_DATA_ONLY", "provenance_status": "verified_frozen_record", "notes": "Missing data remain NA; no fuzzy name matching; compound-level MoA/viability/safety not inferred."},
        {"module_or_figure": "Machine learning", "analysis_role": "conditional only", "parameter_or_item": "new model", "recorded_value": "not activated; any future task requires fixed label/features/independent cohort, nested CV and locked external validation", "source_path": "14-SCI机器学习建模与验证器/输出/机器学习建模方案.md;14-SCI机器学习建模与验证器/输出/machine_learning_audit.json", "status": "NOT_APPLICABLE_CURRENT_STUDY", "provenance_status": "verified_scope_record", "notes": "Legacy Figure7 model artifacts are not treated as a new clinical prediction result."},
        {"module_or_figure": "Reproducibility", "analysis_role": "software environment", "parameter_or_item": "current and historical version policy", "recorded_value": "primary Windows Python=3.11.5; local recorded R=4.5.0; historical exact versions marked historical_exact_version_not_recoverable when evidence absent", "source_path": "metadata/reproducibility/software_environment_manifest.json;17-SCI生信代码与可重复性核查器/输出/软件版本与环境记录.md", "status": "PARTIALLY_CLOSED", "provenance_status": "verified_audit_record", "notes": "Current environment drift is not substituted for historical run versions."},
    ]


METHODS_TEXT = r"""# Methods

## 1. Study design and evidence architecture

This study used a cross-sectional, multi-cohort computational design to organize hepatocyte states from reference and injury-associated contexts to a CNV-supported malignant-like endpoint in hepatocellular carcinoma (HCC). The analytical framework was specified as a partially ordered three-axis architecture comprising HNF4A/PPARA-associated hepatocyte identity loss, AP-1/CEBPB/EGR1-associated stress transition, and a later SOX4-associated malignant-state stabilization programme. Arrows in the computational workflow represent data and analysis dependencies; they do not specify an obligatory biological cascade. Cell states, CNV-supported malignant-like labels, trajectory positions, regulatory activities and virtual perturbations were treated as distinct evidence layers.

The discovery layer used integrated single-cell and single-nucleus expression objects. The reference layer supplied healthy or injury-related hepatocyte contexts and cross-study annotation support. CNV, trajectory and fate analyses were used to define and position computational states. Canonical SCENIC/cisTarget, CellOracle and scTenifoldKnk were used for regulatory association, candidate prioritization and virtual perturbation. TCGA-LIHC, ICGC-LIRI-JP and independent single-cell resources were kept as external corroboration layers. Patient-level matched analyses were reported as supportive or exploratory according to their prespecified eligibility and sample size. No new human participants, biospecimens or adult-HCC cell experiments were included in this computational study.

## 2. Public datasets, accessions and analytical roles

The primary integrated count layer comprised six public human datasets: GSE149614 (HCC malignant endpoint and hepatocyte-state discovery), GSE151530 (conditional HCC/ICC tumour-evolution reference), GSE174748 (healthy-to-NAFLD injury reference), GSE185477 (healthy hepatocyte and modality reference), GSE202379 (healthy-to-NASH/cirrhosis continuum), and GSE212046 (paired cirrhotic-background/HCC context). The audited counts input contained 420,435 cells after object-level eligibility filtering. Local QC summaries recorded 71,305 cells for GSE149614, 52,107 for GSE151530, 11,909 for the counts-fixed GSE174748 object, 123,437 across the GSE185477 objects, 99,191 for GSE202379 and 62,730 for GSE212046. These local counts describe the project objects and are not substituted for the repository's study-level sample descriptions.

The HCC atlas was obtained from figshare article 22332568, version 1 (DOI: 10.6084/m9.figshare.22332568.v1), and was used as a cross-study HCC tumour-microenvironment and annotation reference. The available atlas object was not treated as a raw-count matrix and was excluded from the scVI count input. The associated audit records the deposited files as CC BY 4.0; the final manuscript will preserve the dataset attribution and associated publication citation.

External single-cell corroboration used the HCC subsets of GSE156625, CNP0000650 and GSE125449 Set 1/Set 2 as available in the frozen Module 8 cache. The CNP0000650 matrix was recorded as log-TPM rather than established raw counts. These resources were used for axis-level expression recurrence and were not treated as new discovery data or as direct patient-level causal evidence. GSE189175 was analysed separately as a small matched patient-level supportive cohort. GSE326201 was analysed as an exploratory Tier 1 matched cohort; it did not satisfy the formal Tier 2+ validation requirement. TCGA-LIHC and ICGC-LIRI-JP were used as external bulk expression cohorts. GSE282701 was retained as a design-level rescue candidate but was not used for formal validation because the official processing record documents equal-depth normalization and feature-barcode matrix recomputation, leaving unnormalized raw-count provenance unresolved. The duplicate GSE149614 external cache was excluded as discovery leakage.

Repository revision labels, accessions, local download records, source URLs and file inventories were retained in the project intake and public-data provenance records. The Stage20.5 compliance audit classified citation, licence/access and repository acknowledgement as `PASS_WITH_LIMITATION`, ethics/secondary-use wording as `MANUAL_CONFIRMATION_REQUIRED`, and Data Availability as `MANUSCRIPT_READY_WITH_PLACEHOLDERS`. The exact repository release/access date for the TCGA/ICGC caches and the final journal-specific wording remain `AUTHOR ACTION REQUIRED` and must be completed before submission.

## 3. Data provenance, ethics and secondary use

All analysed human data were obtained from public or previously cached repository records within the project scope. The study did not recruit participants, contact participants or collect new biospecimens. The current candidate wording is that the analyses constitute secondary analyses of de-identified public datasets and that original consent and ethics procedures were those of the data-generating studies. This statement requires institutional confirmation before submission and is not an inferred ethics approval or waiver. Public repository access was not interpreted as a blanket waiver of submitter, third-party, privacy or data-use restrictions.

Original source metadata, accession identifiers and repository links were retained in derived records. Raw human matrices and non-public raw reads are not redistributed by this project. The GSE149614 record links raw-data access to EGA managed access, and the GSE212046 record documents that human raw data were not uploaded under the source study's IRB condition. Dataset-specific GEO licence terms were not explicitly identified in the audit and therefore are not represented as CC0 or unrestricted reuse. The HCC atlas file-level attribution follows its figshare record. Final institutional secondary-use wording and any dataset-specific licence/waiver language remain `MANUAL_CONFIRMATION_REQUIRED`.

## 4. Matrix preparation and input quality control

GEO and figshare files were converted into analysis objects while preserving the original downloaded files. The primary single-cell/single-nucleus input was required to contain non-negative count-like values. The frozen QC rules were `min_genes=200`, `min_counts=500`, `max_mito_pct=25`, `min_cells_per_gene=3` and an upper quantile of 0.995 where recorded by the preprocessing scripts. Objects with fewer than 1,000 retained cells, fewer than 5,000 retained genes, a sampled non-integer rate greater than 0.001 or a numeric gene-name rate greater than 0.5 were excluded from the scVI count input. Objects with fewer than 3,000 retained cells were retained as review flags unless a separate hard exclusion rule was met.

Raw counts were retained for model inputs and count-based analyses. A log1p representation was used only for visualization and PCA diagnostics. HCC atlas normalized/non-count objects, the original normalized GSE174748 object, low-cell GSE185477 objects that failed the hard object gate, and the GSE202379 object with unresolved numeric gene indexing were not silently promoted to the primary count input. The resulting counts manifest contained 17 included objects and 5 excluded objects. QC, object paths and exclusion flags are recorded in the preprocessing input manifest.

Doublet and cell-cycle diagnostics were treated as quality annotations. A simulated-doublet k-nearest-neighbour diagnostic used an expected doublet rate of 0.05; the recorded global predicted rate was 0.0559801 and cluster 16 was flagged for review. Cycling clusters were retained as interpretive flags. These diagnostics did not create a separate biological endpoint.

## 5. Integration, dimensionality reduction and cell annotation

For the global integration, scVI was fitted to the raw-count input using 2,000 highly variable genes, a 30-dimensional latent representation, two hidden layers of 128 units, `batch_key=dataset`, 10 epochs, batch size 1,024 and seed 20260601. The recorded run contained 420,435 cells, 2,000 retained variables and 11,010 common genes and used scvi-tools 1.3.3 with PyTorch 2.4.1+cu121 on a GPU. Neighbours, UMAP and Leiden clustering were computed from `X_scVI` with 30 neighbours, `min_dist=0.3`, Leiden resolution 1.0 and seed 20260601; 62 Leiden clusters were recorded. Batch diagnostics were calculated at both dataset and `sample_id` levels. Improvement in dataset mixing was not interpreted as complete removal of sample-level residual structure.

Major cell-type labels were assigned using the CellTypist `Healthy_Human_Liver.pkl` model (CellTypist 1.7.1; Scanpy 1.11.5), existing SingleR/celldex reference outputs (SingleR 2.14.0; celldex 1.22.0) and a manual marker review. The manual review used 57 markers, with recorded thresholds `min_panel_z=0.35`, `min_panel_pct=0.05` and external-conflict confidence 0.6. High-confidence labels were used as seeds; ambiguous, mixed and conflict labels were retained. scANVI was then used for semi-supervised label propagation with 20 epochs, batch size 1,024, 1,000 samples per label and seed 20260601. Strict malignant labels were retained only when supported by the module 3 malignant call categories; review labels were not converted into definitive malignant calls.

## 6. CNV inference and computational malignant-like state definition

CNV evidence was derived from the project CopyKAT workflow and an explicitly separated expression-based proxy. The recorded CopyKAT invocation used `id.type=S`, `cell.line=no`, `ngene.chr=5`, `min.gene.per.cell=200`, `LOW.DR=0.05`, `UP.DR=0.1`, `win.size=25`, `KS.cut=0.1`, Euclidean distance, `genome=hg20`, `output.seg=FALSE`, `plot.genes=FALSE` and one core. Normal reference cells were selected within the corresponding sample and analyses were chunked according to the input manifest. CopyKAT outputs were retained as `aneuploid`, `diploid` or `not.defined` where available.

For missing or unavailable CopyKAT calls, a cytoband-binned expression proxy used an HGNC gene map, bins of 100 genes, a burden z-score threshold of 3.0 and a high-bin fraction threshold of 0.08. The proxy compared candidate cells with within-sample reference cells and retained borderline calls separately from aneuploid-proxy calls. Malignant-associated expression evidence was flagged when the recorded marker score z-score was at least 0.8, the mean log1p CPM was at least 3.5 or a malignant-hepatocyte state label was present. Proliferation review used a score z-score threshold of 0.8 or a proliferating state label. Tumour-source classes, normal/adjacent classes and unknown HCC-source classes were kept explicit.

The final label `CNV-supported malignant-like` denotes a computational state assembled from expression, source context and CNV-associated evidence. It is not equivalent to DNA-sequencing-confirmed clonality, a clinical stage or a proven malignant transformation event. CopyKAT's exact historical package version and runtime are not recoverable from the frozen record and are therefore flagged in the parameter table rather than reconstructed.

## 7. Trajectory, relative state positioning and fate analysis

The trajectory object contained 283,498 cells and 2,000 variables. The main trajectory was oriented from reference hepatocytes to CNV-supported malignant-like cells using the scANVI representation; malignant review cells and a `main_strict`/`include_review` sensitivity distinction were retained. Disease-stage labels were derived from source/sample class and computational trajectory roles and were not mapped to TNM or BCLC clinical stages.

Monocle3, Slingshot on the scANVI representation, and Slingshot on hepatocyte PCA were used as available trajectory methods. The exact historical Monocle3 and Slingshot versions are not present in the recoverable run records and remain `historical_exact_version_not_recoverable`. Relative programme trends used the recorded 10-bin trajectory summaries. The corrected temporal-positioning implementation used a 201-point grid, baseline/search boundaries of 0.10, a slope-search end of 0.95, minimum onset and crossing runs of 5 and 3, minimum effect fraction 0.05, minimum effect of 0.25 baseline standard deviations, derivative fraction 0.10, precedence tolerance 0.005, 1,000 bootstrap replicates and seed 20260805. Activity bands represent relative programme prominence; their right edges are prominence-weighted fades and do not estimate discrete termination events. DPT and CellRank pseudotime were unavailable in the Figure 5 positioning record.

CellRank fate analysis used the `main_strict` object with 5,000 cells and 2,000 genes. Terminal-state annotations included late non-CNV and CNV-supported malignant states; fate-probability row sums were audited to approximately one. CellRank 2.0.7, Scanpy 1.11.5, anndata 0.12.16, NumPy 1.26.4 and pandas 2.2.2 were recorded for this run. Fate-driver associations were used for prioritization and were not interpreted as causal drivers.

## 8. Frozen three-axis programme definitions and scoring

The three-axis definitions were frozen before external validation and were not selected using bulk outcomes. The identity axis was defined around HNF4A/PPARA-associated identity loss, the stress-transition axis around AP-1/CEBPB/EGR1-associated transition genes, and the malignant-state axis around a SOX4-associated programme. The provisional core transcription factors were HNF4A and PPARA for identity, CEBPB, EGR1, JUN, FOS and JUND for stress transition, and SOX4 for the malignant-state axis. The direct SOX4 candidate count remained below the target range; unsupported transcription factors were not added to satisfy a quota.

Different analysis modules used their own frozen gene-set version and these versions were not silently merged. The Figure 5 state-positioning report recorded identity, stress and SOX4 programme sizes of 35, 153 and 41 genes, respectively, with dataset-wise median/MAD robust z-scores clipped to [-5,5] and UCell used for sensitivity scoring. The Figure 7 v2/GSE189175 patient-level validation definition recorded 15 identity genes (core HNF4A and PPARA), 95 stress-transition genes (core CEBPB, EGR1, JUN, JUNB, FOS, JUND and ATF3), and 31 SOX4-associated genes (core SOX4). For this validation layer, scores were calculated from `log2(CPM+1)`, gene-wise z-scores across the six patient-by-condition profiles, an unweighted mean across genes and zero for zero-variance genes; identity retention was sign-inverted to obtain the identity-loss score. Alternative scoring definitions were sensitivity analyses only.

Candidate prioritization integrated six evidence dimensions: CellRank/trajectory association, canonical SCENIC/cisTarget support, CellOracle support, scTenifoldKnk support, state/CNV/trajectory association and axis-level external recurrence. A transparent support fraction was calculated as the number of supported dimensions divided by six. Degree, MCC or a single network-connectivity statistic was not sufficient for candidate selection. External recurrence was used at the axis level and was not converted into TF-level clinical or prognostic evidence.

## 9. Regulatory-network inference and virtual perturbation

### Canonical SCENIC/cisTarget

The formal canonical SCENIC branch used 9,512 full-expression driver-union cells and 11,923 genes after all-zero cleanup. GRNBoost2 used 1,767 expressed TFs, seed 777 and the run-recorded eight workers, generating 1,131,800 edges. cisTarget motif pruning used hg38/RefSeq 80 resources, mc_v10_clust 10-kb and proximal promoter ranking families and the HGNC motif annotation table. The recorded context settings were a rank threshold of 5,000, an AUCell threshold of 0.05, NES threshold of 3.0, `min_genes=20`, up to 50 targets and regulator thresholds of 5, 10 and 50. AUCell activity was calculated for 422 canonical regulons, of which 319 had at least 10 targets and 103 were retained as small-regulon review objects. CellRank association used the 5,000-cell fate subset. The canonical branch used pySCENIC 0.12.1, ctxcore 0.2.0, anndata 0.12.16, NumPy 1.26.4, pandas 2.2.2, SciPy 1.13.1, PyArrow 24.0.0, Matplotlib 3.10.8 and statsmodels 0.14.6. Because only GRNBoost2 seed 777 was run in the canonical branch, regulon activity is reported as association evidence and not as direct TF binding or causal activity. The upstream environment snapshot lists a different worker-count field; this administrative discrepancy is flagged for author reconciliation before the final Methods is frozen.

### CellOracle

CellOracle used a curated hg38 promoter prior network and the 15-TF panel defined by the upstream candidate-selection record. State-specific GRNs were fitted across five recorded states using CellOracle 0.20.0 in a Python 3.10.20 WSL2 environment. The recorded GRN parameters were `alpha_links=10`, `alpha_simulation=1`, 20 bagging iterations, 30 PCA dimensions, `k=30`, `p_threshold=0.001`, `threshold_number=10000` and eight jobs. Virtual perturbations were knockout simulations with three propagation steps, delta clipping, 100 neighbours, a sampled fraction of 0.3, `sigma_corr=0.05`, 20 grid steps, 50 grid neighbours, 50 top genes per state and seed 15071990. Negative-control analyses used five recorded control TFs. Restore or overexpression simulations were not available. AP-1 summaries were derived from member-wise perturbation results and were not treated as a combined suppression experiment. CellOracle output was therefore interpreted as directional computational prioritization.

### scTenifoldKnk

The primary Module 7 workflow retained the recorded 3,000-gene network background and 45 TF-level runs with `qc=false`, `nc_nNet=1`, `nc_nCells=100`, one core and seed 11. A separate additive reproducibility audit standardized the selected HNF4A, EGR1 and SOX4 runs to `nc_nNet=10`, `nc_nCells=500`, `nc_nComp=3`, `ma_nDim=2`, `qc=false`, `qc_minCells=3` and three seeds (15071990, 15071991 and 15071992); HNF4A and SOX4 were rerun and EGR1 was retained from a matching three-seed run. The primary Figure 2/3 frozen panel counts remain those specified by the corresponding frozen figure/source-data records. scTenifoldKnk manifold-alignment distances were treated as unsigned network displacement. They were not converted into signed activation, suppression, upregulation, downregulation or experimental knockout effects. The exact historical package version is not fully recoverable and is not replaced with a current installation.

## 10. Differential expression and pathway analysis

Cell-level differential testing was prespecified for exploratory contrasts using Wilcoxon tests implemented through Seurat `FindMarkers`, `min.pct=0.10`, absolute log2 fold change at least 0.25 and BH-FDR below 0.05. The frozen inventory does not contain a complete global traditional cell-level DEG table for all primary contrasts; this specification must not be presented as a completed universal DEG result unless the author supplies the corresponding result table.

For sample/patient-level inference, cells from the same dataset, patient and condition/state were aggregated into pseudobulk profiles. Cells or nuclei were never treated as biological replicates. Complete paired designs used `~ patient_id + condition`; unpaired designs used `~ condition` only where the source identity and estimability rules allowed it. The edgeR workflow was `DGEList` construction, `filterByExpr(y, design=design)`, TMM normalization with `calcNormFactors`, robust dispersion estimation, robust quasi-likelihood fitting with `glmQLFit`, and `glmQLFTest`. Multiple testing used the Benjamini-Hochberg procedure. The prespecified pseudobulk reporting threshold was FDR<0.05 with absolute log2 fold change at least 0.50. Historical eligible comparisons and GSE326201 retained all result rows and significant subsets; any incomplete or non-estimable comparison was kept as such.

The GSE189175 matched analysis used six raw integer patient-by-condition profiles from three complete patients, a fixed hepatocyte/epithelial marker compartment and the contrast Tumour versus Adjacent under `~ patient_id_original + condition`. The exact recorded runtime was R 4.3.1 with edgeR 3.42.4 and limma 3.56.2. GSE326201 used direct GEO patient/sample mapping, InferCNV recall and final eligibility thresholds of at least 50 reference cells and 20 malignant cells per patient; three complete pairs remained and the cohort was retained as Tier 1 exploratory. The GSE326201 InferCNV record used R 4.6.1, infercnv 1.28.0 and Matrix 1.7.5; its edgeR record used R 4.5.0, edgeR 4.8.2 and limma 3.66.0.

Pathway analysis used one-sided hypergeometric ORA for GO biological process, KEGG and Reactome when the corresponding input and matched background were available. Network-perturbation enrichment used the source-specific network background, including the recorded 3,000-gene background for the Module 7 network. Preranked GSEA or GSVA was used only where the corresponding frozen record identifies the ranking and database. BH correction was applied within prespecified analysis families. Network-displacement enrichment was kept separate from traditional DEG enrichment, and the absence of an FDR-supported pathway was retained rather than filled by nominally significant or untested pathways.

## 11. Relative network direction and mediation summaries

CellOracle and scTenifoldKnk outputs were integrated through programme-level effect summaries, gene-set overlap, rank concordance and forward/reverse perturbation comparisons. Module 9.2 recorded eight signature sets, 48 asymmetry rows and 180 concordance rows. Module 9.3 used 1,000 bootstrap replicates for partial mediation summaries where outcome and path variables were available. These analyses were used to grade computational direction and partial coupling. They did not fit genetic epistasis, did not establish a direct causal cascade and did not convert the historical competing-model label `Model 1: Linear cascade` into the biological interpretation; the frozen biological architecture remains partially ordered.

## 12. External expression recurrence and exploratory clinical association

Figure 7 v2 used programme definitions frozen from the single-cell/network analyses before bulk outcomes were inspected. Patient-level scores were obtained by averaging eligible samples within each patient, and cohort-specific tumour-normal expression recurrence was summarized with Hedges' g and 95% confidence intervals. A two-cohort random-effects synthesis was treated as exploratory. Matched random signatures preserved gene number, mean expression, variance and detection-rate strata and were evaluated using the same recurrence and TCGA-only secondary clinical models. The Figure 7 v2 script used seed 20260812 and 1,000 matched random signatures per cohort-axis combination. The primary programme score is an unsigned associated target-programme score because the prespecified signed-target coverage rule was not met for all axes.

TCGA-LIHC was eligible for secondary clinical association analyses when the frozen clinical fields were defined and complete-case criteria were met. ICGC-LIRI-JP expression-to-donor mapping was retained for expression recurrence. Its OS branch was analysed only under the separately audited exploratory policy: a continuous frozen programme score per one within-ICGC tumour standard deviation, 231 donor-linked records and 42 events, with no data-driven cut-off search and with PH/nonlinearity diagnostics retained. This branch is `ESTIMABLE_BUT_NOT_VALIDATED` and is restricted to Supplementary/Extended Data. ICGC Age, Gender, Stage, `fustat` and `futime` semantics and release-specific provenance remain incomplete; no current codebook or release is guessed. No claim of independent prognostic validation, clinical prediction or clinical utility is made.

## 13. Exploratory transcriptomic rescue analysis

Figure 8 used a frozen continuous landmark-space representation of 978 genes. Each evidence coordinate was robust-scaled to [-1,1], combined using the frozen weights 0.35, 0.15, 0.10, 0.10, 0.10, 0.15 and 0.05, and shrunk according to directional agreement; unsupported coordinates remained zero. DrugReflector V3.5 outputs and frozen checkpoints were compared with related L1000FWD/CLUE resources, curated exact identifier mappings, Figure 6 network compatibility and PRISM cancer-cell phenotype summaries. The mapping hierarchy was InChIKey, canonical standardized name and exact BRD identifier; fuzzy name matching was not used. Two thousand matched null signatures were used for specificity benchmarking, and 15 prespecified signature variants were retained.

Candidate tiers required the recorded combinations of rank stability, three-fold agreement, matched-null specificity, CMap-family corroboration, curated mechanism/network support, non-broad-cytotoxic phenotype and coverage. The analysis is currently `EXTENDED_DATA_ONLY`; no Tier A, B or C candidate passed the joint gate. The output is a hypothesis-generating prioritization and does not establish efficacy, safety, normal-cell selectivity, treatment recommendation or clinical actionability.

## 14. Software environment and reproducibility

The primary local analysis records identify Windows 10/11, Python 3.11.5 and R 4.5.0 for the major project workflows. The scVI/scANVI runs record scvi-tools 1.3.3 and PyTorch 2.4.1+cu121. CellRank records CellRank 2.0.7, Scanpy 1.11.5 and anndata 0.12.16. Canonical SCENIC records pySCENIC 0.12.1, ctxcore 0.2.0, NumPy 1.26.4, pandas 2.2.2, SciPy 1.13.1, PyArrow 24.0.0 and statsmodels 0.14.6. CellOracle records version 0.20.0 and a Python 3.10.20 WSL2 runtime. Figure 5–8 R records include ggplot2 4.0.3, ggsci 5.0.0, patchwork 1.3.2, data.table 1.18.4, mgcv 1.9.4, metafor 5.0.1, UCell 2.14.0, survival 3.8-6, jsonlite 2.0.0, MASS 7.3-65, ggrepel 0.9.8 and related packages as listed in the software table.

The project dependency scan identified 72 R package entries and 28 Python distribution/module entries from actual source-code calls, including `library`, `require`, `requireNamespace`, namespace calls and Bioconductor signals. The full scan table records source files, current loadability, current versions, historical versions and closure status. When the exact historical version cannot be recovered, the record is explicitly labelled `historical_exact_version_not_recoverable`; a current package version is not used as a historical substitute. Current environment drift between system Python and `.venv-scvi` is retained as provenance and does not imply that the current runtime reproduced a historical result.

All Methods claims are traceable to the upstream manifests, source-data records, scripts, run reports and frozen gate records listed in the Stage21 input audit. Stage21 generated no new biological analysis, did not modify Figures 1–8, did not alter Results or Discussion, and did not reopen Stage19.

## 15. Statistical reporting and analysis boundaries

All inferential p-values were adjusted with BH within the prespecified analysis family when the relevant test family was available. Effect sizes, confidence intervals, patient/sample units, missingness and evidence status were retained in the source tables. The analysis distinguishes cells, samples, patients/donors and datasets. Machine learning was not activated as a new task; any historical prediction artefact remains a legacy audit object and is not used to support a new clinical model.

The study is observational and computational. CNV-supported malignant-like calls, relative pseudotime, CellRank fate probabilities, regulon activity, virtual perturbations, external expression recurrence and transcriptomic rescue scores support association or hypothesis prioritization. They do not alone demonstrate direct TF binding, genetic dependency, a strict temporal cascade, patient-level causality, clinical prediction, treatment efficacy or safety.

## 16. Data and code availability

The exact submission-ready Data Availability and Code Availability wording, including unresolved repository and licence fields, is provided in `输出/数据可用性与代码可用性段落.md`. Before submission, the author must insert the verified persistent code/data repository identifier, final access dates, journal-specific acknowledgement text, original-article citations and institutionally approved secondary-use/ethics wording. These administrative items do not reopen Stage19 and do not change the scientific evidence boundaries.

"""


METHODS_TEXT = METHODS_TEXT.replace(
    "leaving unnormalized raw-count provenance unresolved.",
    "leaving unnormalized raw-count provenance unresolved; its current status is `BLOCKED_PROVENANCE_UNRESOLVED`.",
    1,
).replace(
    "The analysis is currently `EXTENDED_DATA_ONLY`; no Tier A, B or C candidate passed the joint gate.",
    "The current scope retains the frozen `EXTENDED_DATA_ONLY` gate; Tier A/B/C eligibility is recorded in the upstream Figure 8 audit and is not reassessed by this Methods draft.",
    1,
)


DATA_CODE_TEXT = r"""# Data Availability and Code Availability

## Data Availability

The public datasets analysed in this study are available from the following repositories: GEO accessions GSE149614, GSE151530, GSE174748, GSE185477, GSE202379, GSE212046, GSE189175 and GSE326201; the HCC atlas is available from figshare, version 1, DOI: 10.6084/m9.figshare.22332568.v1; TCGA-LIHC data are available through the Genomic Data Commons; and ICGC-LIRI-JP data are available through the ICGC project records. GSE189175 was used as a small matched patient-level supportive cohort. GSE326201 was used as an exploratory Tier 1 cohort. GSE282701 was evaluated as a design-level rescue candidate but was not used for formal validation because the provenance of unnormalized raw counts remained unresolved.

The source data underlying the figures, derived statistics and analysis-ready tables are retained in the project source-data and metadata directories and will be deposited in `[AUTHOR ACTION REQUIRED: insert verified repository and persistent identifier]` before submission. Readers should obtain the original public files from the repositories and accessions above and cite the corresponding data-generating articles. The final manuscript must include the exact access dates and journal-required dataset citations.

The current audit did not confirm an explicit dataset-specific licence for the GEO Series records; public repository access must not be interpreted as CC0 or unrestricted reuse. The HCC atlas figshare record is identified as CC BY 4.0 in the project audit and requires attribution. The GSE149614 record links raw-data access to EGA managed access, and GSE212046 documents an IRB-based non-upload condition for human raw data. Human-derived raw matrices and non-public raw reads are not redistributed by this project. ICGC-LIRI-JP OS analyses remain exploratory and `ESTIMABLE_BUT_NOT_VALIDATED`.

## Code Availability

The analysis scripts, configuration records, run reports, source-data traceability records and reproducibility manifests are maintained in the project `scripts/`, `metadata/`, `reports/` and stage-specific directories. The final public code release will be deposited at `[AUTHOR ACTION REQUIRED: insert verified code repository, release tag and DOI/Zenodo or other persistent identifier]`. The release should include the exact source-data manifest, software environment manifest, parameter table and instructions for reproducing the figure-level source data. External raw human files subject to source-repository or consent restrictions will not be redistributed; the code release will provide accession-based retrieval instructions and derived non-sensitive files where permitted.

## Ethics / secondary-use wording for author confirmation

Candidate wording: This study is a secondary analysis of de-identified public datasets. No new human participants were recruited or contacted and no new biospecimens were collected. Ethics approvals and consent procedures were those of the original data-generating studies. The authors will confirm that this secondary use is covered by applicable institutional policy before submission.

`MANUAL_CONFIRMATION_REQUIRED`: institutional IRB/ethics wording, dataset-specific licence or waiver language, final Zotero reference matching, repository acknowledgement, access dates and persistent code/data identifiers.

## 中文核对

- 公开数据库的访问权不等于对第三方权利的全面放弃；GEO、GSE189175 和 ICGC-LIRI-JP 的最终使用措辞仍需作者确认。
- GSE282701 只能写为未纳入正式验证的 provenance candidate；不能在 Methods、Results 或 Discussion 中写入其 CNV、pseudobulk 或三轴结果。
- GSE326201 保持 Tier 1 exploratory；ICGC OS 保持 `ESTIMABLE_BUT_NOT_VALIDATED`；Figure 8 保持 `EXTENDED_DATA_ONLY`。
- 代码仓库永久标识、最终访问日期和伦理/二次使用措辞需要在投稿前补齐。
"""


def build_input_audit() -> str:
    return r"""# Stage21 输入材料审计

## 进入依据

- 当前阶段：`21-SCI生信研究Methods撰写器`。
- 用户已明确授权正式进入 Stage21，并限定为基于冻结的真实分析流程、参数、软件版本、数据来源和审计记录撰写 Methods。
- Stage19：`STAGE19_CLOSED_WITH_LIMITATIONS`；Stage20 eligibility：`YES`；Stage20 已完成冻结证据 Discussion 草稿；Stage20.5 合规状态：`PASS_WITH_LIMITATION`。
- Stage21 不重跑生物学分析，不修改 Figure 1–8、Results 或 Discussion，不升级 GSE326201、GSE282701、ICGC OS 或 Figure 8 证据等级，不重新打开 Stage19。

## 已读取的冻结输入

1. Stage08 数据集纳排、下载、样本分组和公共数据 accession/version/licence 记录。
2. Stage10 分析流程图、模块-软件表和关键参数预设。
3. Stage11 counts/QC/标准化/scVI/批次诊断记录。
4. Stage12 差异分析、富集、阈值、对比设计和传统 DEG 状态记录。
5. Stage13 TF-gene 网络、多证据候选筛选和三轴候选缺口记录。
6. Stage14 机器学习禁用/条件性启用边界。
7. Stage15 外部验证结果、ICGC OS 状态和未来实验计划。
8. Stage16 图件契约、图注、source-data 索引和 Figure 7/8 证据边界。
9. Stage17 代码依赖扫描、版本环境、运行顺序和历史版本不可恢复标记。
10. Stage18 Results 主线与证据等级；Stage19 Results 审计和最终闭合判定。
11. Stage19B GSE189175 matched patient-level validation；GSE326201 Tier 1 exploratory；GSE282701 raw-count provenance rescue。
12. Stage20 Discussion 冻结输入和 Stage20.5 Public Data Compliance Audit。

## 冻结科学边界

- 三轴采用部分有序架构，不写严格 HNF4A/PPARA → AP-1/CEBPB/EGR1 → SOX4 因果级联。
- HNF4A/PPARA、AP-1/CEBPB/EGR1 和 SOX4 均为关联/候选调控层；计算网络不能替代直接因果或实验扰动。
- GSE326201 为 Tier 1 exploratory，正式 Tier 2+ 不可用。
- GSE282701 因未归一化 raw-count provenance unresolved 未用于正式验证。
- ICGC OS 为 `ESTIMABLE_BUT_NOT_VALIDATED`，只允许预先定义的 donor-level 连续变量探索，并限 Supplementary/Extended Data。
- Figure 8 为 `EXTENDED_DATA_ONLY`，不得写疗效、安全性、临床可操作性或治疗建议。
- 成人 HCC 细胞实验是未来计划，不能写成已完成实验。

## Methods 书写审计结论

- 可直接写入：已由 run report、脚本或冻结审计明确记录的输入、QC、模型、参数、软件和统计流程。
- 必须标记：Monocle3/Slingshot 历史精确版本、CopyKAT 历史精确版本、部分历史依赖包精确版本、SCENIC worker-count 记录冲突、TCGA/ICGC release/access date/codebook、最终 repository identifier、dataset-specific licence/waiver 和机构伦理措辞。
- 传统 cell-level DEG 的全局完成状态不可从现有清单确认；Methods 将其写成预设的探索性方案并明确 `SPECIFIED_NOT_COMPLETE_AS_GLOBAL_TABLE`。

## 输出策略

本阶段生成独立 Stage21 目录，所有新增 Markdown 首行带项目版权水印。参数表同时保存实际代码依赖和方法参数；输入清单保存来源文件的存在性、大小、修改时间及可行时的 SHA-256。Stage21 不覆盖既有上游文件。
"""


def build_quality_rows() -> list[dict[str, object]]:
    return [
        {"check_id": "21-01", "item": "Stage19/20/20.5 entry authorization", "status": "PASS", "evidence": "STAGE19_FINAL_GATE.json; stage20_run_record.json; PUBLIC_DATA_COMPLIANCE_GATE.json", "risk": "Stage21 is limited to frozen-evidence Methods."},
        {"check_id": "21-02", "item": "Required Stage08–17 input records", "status": "PASS", "evidence": "Stage21 input manifest; upstream stage audits", "risk": "Some historical records retain legacy wording."},
        {"check_id": "21-03", "item": "Primary data and accession roles", "status": "PASS_WITH_LIMITATION", "evidence": "Stage08 intake and public-data provenance", "risk": "Per-file download/release fields are incomplete for some external caches."},
        {"check_id": "21-04", "item": "Raw counts/QC/integration parameters", "status": "PASS", "evidence": "Stage11 audit; scVI/scANVI run reports; Stage10 parameter table", "risk": "Dataset-level mixing does not remove sample-level residual risk."},
        {"check_id": "21-05", "item": "Annotation/CNV/malignant-like definition", "status": "PASS_WITH_LIMITATION", "evidence": "CellTypist/SingleR reports; CopyKAT/proxy scripts and reports", "risk": "CopyKAT historical exact version unavailable; labels are computational."},
        {"check_id": "21-06", "item": "Trajectory/CellRank method traceability", "status": "PASS_WITH_LIMITATION", "evidence": "trajectory and CellRank reports; Figure5 report", "risk": "Monocle3/Slingshot historical exact versions unavailable; pseudotime is relative."},
        {"check_id": "21-07", "item": "SCENIC/cisTarget parameters and boundary", "status": "PASS_WITH_LIMITATION", "evidence": "canonical SCENIC reports and environment", "risk": "single GRNBoost2 seed; worker-count reconciliation required."},
        {"check_id": "21-08", "item": "CellOracle/scTenifoldKnk traceability", "status": "PASS_WITH_LIMITATION", "evidence": "CellOracle reports; scTenifoldKnk v2 audit", "risk": "Virtual perturbation is not direct causality; scTenifold historical package exact version incomplete."},
        {"check_id": "21-09", "item": "Differential/enrichment methods", "status": "PASS_WITH_LIMITATION", "evidence": "Stage12 plans; edgeR scripts and run logs", "risk": "Global traditional cell-level DEG table is not confirmed complete."},
        {"check_id": "21-10", "item": "GSE189175/GSE326201 patient-level methods", "status": "PASS_WITH_LIMITATION", "evidence": "Stage19B and GSE326201 final reports/audits", "risk": "Small exploratory cohorts; Tier 2+ remains unavailable."},
        {"check_id": "21-11", "item": "Figure7 external recurrence and ICGC OS boundary", "status": "PASS_WITH_LIMITATION", "evidence": "Figure7 v2 report; ICGC OS policy and final gate", "risk": "Release/codebook and clinical semantics remain incomplete; no prognostic-validation claim."},
        {"check_id": "21-12", "item": "Figure8 rescue-method boundary", "status": "PASS", "evidence": "Figure8 v2 report and integrated score definition", "risk": "EXTENDED_DATA_ONLY remains mandatory."},
        {"check_id": "21-13", "item": "Software dependency and version table", "status": "PARTIALLY_CLOSED", "evidence": "metadata/reproducibility/software_environment_manifest.json", "risk": "Historical exact versions unavailable for selected packages; current versions are not substitutes."},
        {"check_id": "21-14", "item": "Data/Code availability wording", "status": "MANUAL_CONFIRMATION_REQUIRED", "evidence": "Stage20.5 compliance audit and Data Availability draft", "risk": "Repository identifier, access dates, licence/waiver and institutional ethics wording pending."},
        {"check_id": "21-15", "item": "No biological rerun or upstream modification", "status": "PASS", "evidence": "Stage21 run record; output namespace isolated", "risk": "No upstream scientific state is changed by this stage."},
        {"check_id": "21-16", "item": "Methods reproducibility gate", "status": "BLOCKED_FOR_FINAL_MANUSCRIPT", "evidence": "Missing historical exact versions and author-controlled provenance fields", "risk": "Draft is usable for author review but not final submission text until action items close."},
    ]


def build_score_rows() -> list[dict[str, object]]:
    return [
        {"评价项": "证据完整性", "评分0-5": 4, "证据来源": "Stage08–20.5 input audit and Stage21 parameter table", "问题说明": "核心流程可追溯；若干历史版本和外部 release 字段仍待补充。", "是否触发硬阻断": "否", "修正建议": "补齐作者控制的 release/access-date/licence/codebook 字段。"},
        {"评价项": "结果-图表一致性", "评分0-5": 5, "证据来源": "Stage19 final gate; Stage20 frozen Discussion; Figure contract", "问题说明": "Methods preserves Figure7/8 and three-axis boundaries; no figure was changed.", "是否触发硬阻断": "否", "修正建议": "整合全文时继续使用同一证据术语。"},
        {"评价项": "方法可重复性", "评分0-5": 3, "证据来源": "Stage10–17 records; GSE189175/GSE326201 run manifests", "问题说明": "主要参数已冻结；部分历史 exact version、worker-count record 和外部 release 尚未闭合。", "是否触发硬阻断": "是（最终投稿文本）", "修正建议": "保留 AUTHOR ACTION REQUIRED，禁止猜测或用当前版本冒充历史版本。"},
        {"评价项": "逻辑连贯性", "评分0-5": 5, "证据来源": "Stage18 storyline; Stage20 Discussion", "问题说明": "从数据、QC、状态、网络到外部验证的顺序清晰。", "是否触发硬阻断": "否", "修正建议": "保持 discovery/validation/clinical exploratory 分层。"},
        {"评价项": "路线专属规范符合度", "评分0-5": 4, "证据来源": "Stage21 skill; nature-data audit", "问题说明": "Methods/Data/Code 可用性均有专门段落；伦理和许可仍需作者确认。", "是否触发硬阻断": "否", "修正建议": "投稿前按目标期刊更新文字。"},
        {"评价项": "目标期刊适配度", "评分0-5": 3, "证据来源": "用户高影响力综合期刊目标；期刊尚未最终确定", "问题说明": "科学方法草稿可作为 Nature-leaning 基础，但格式、长度和具体数据政策未适配。", "是否触发硬阻断": "否", "修正建议": "选刊后做期刊级格式审校。"},
        {"评价项": "夸大或编造风险", "评分0-5": 5, "证据来源": "Stage19 final gate; Figure7/8 policy; Stage20.5 audit", "问题说明": "保留 partial、heterogeneous、ESTIMABLE_BUT_NOT_VALIDATED 和 EXTENDED_DATA_ONLY。", "是否触发硬阻断": "否", "修正建议": "整合时不把候选或计算关联改写为因果。"},
        {"评价项": "是否允许进入下一写作环节", "评分0-5": 0, "证据来源": "Stage21 user confirmation gate", "问题说明": "本轮仅完成 Stage21 草稿；需作者确认 Methods 草稿及人工核查项后，才可考虑 Stage22。", "是否触发硬阻断": "是", "修正建议": "等待作者逐项确认；本轮不自动进入 Stage22。"},
    ]


def build_human_review() -> str:
    return r"""# Stage21 需要人工核查

以下事项必须由作者确认或补充后，Methods 才能作为最终投稿文本使用。它们不改变 Stage19 的关闭状态，也不要求在本轮重跑生物学分析。

1. `AUTHOR ACTION REQUIRED`：确认 GSE149614、GSE151530、GSE174748、GSE185477、GSE202379、GSE212046、GSE189175、GSE326201、TCGA-LIHC 和 ICGC-LIRI-JP 的最终 accession、版本/release、访问日期、下载记录和使用资格。
2. `MANUAL_CONFIRMATION_REQUIRED`：确认公共人类数据二次分析的机构/伦理负责人措辞；当前候选文本不能替代机构确认或 IRB/豁免文件。
3. `MANUAL_CONFIRMATION_REQUIRED`：完成 GSE326201、GSE282701 原始论文与最终 Zotero 条目的匹配；不能根据缺失 GEO citation 猜测文章。
4. `MANUAL_CONFIRMATION_REQUIRED`：确认 GEO、GSE189175、ICGC-LIRI-JP 在没有明确 dataset-specific licence 时的最终使用、引用和 acknowledgement 措辞。
5. `AUTHOR ACTION REQUIRED`：确认目标期刊要求的 repository acknowledgement、数据访问日期、代码仓库永久标识、release tag/DOI 和最终 Data Availability/Code Availability 文字。
6. `AUTHOR ACTION REQUIRED`：确认 CopyKAT 的历史实际包版本/运行环境是否有项目外记录；若无，保留 `historical_exact_version_not_recoverable`。
7. `AUTHOR ACTION REQUIRED`：确认 Monocle3、Slingshot 和 scTenifoldKnk 的历史实际包版本；若无，保留不可恢复标记，不能使用当前版本代替。
8. `AUTHOR ACTION REQUIRED`：核对 canonical SCENIC 运行记录中 GRNBoost2 worker 数（run report=8；环境快照另有 worker 字段）后再冻结最终 Methods 表述。
9. `MANUAL_CONFIRMATION_REQUIRED`：确认 cell-level FindMarkers/Wilcoxon 是否有完整全局 DEG 输出；当前草稿仅将其列为预设探索性方法，不宣称已完成全局结果。
10. `AUTHOR ACTION REQUIRED`：确认是否接受本 Stage21 Methods 草稿的英文投稿语体和其中所有 `AUTHOR ACTION REQUIRED` 标记；未确认前不得自动进入 Stage22。

## 科学边界复核

- GSE326201：`Tier 1 exploratory`；正式 Tier 2+ 不可用。
- GSE282701：`BLOCKED_PROVENANCE_UNRESOLVED`，未用于正式 validation、InferCNV、pseudobulk 或 edgeR 结论。
- ICGC OS：`ESTIMABLE_BUT_NOT_VALIDATED`，仅 donor-level 连续变量探索，Supplementary/Extended Data only。
- Figure 8：`EXTENDED_DATA_ONLY`，无疗效、安全性或临床可操作性表述。
- 成人 HCC 细胞实验：未来计划，当前无已完成实验结果。

下面是基于当前材料的推荐复制回答；你可以直接复制发送，也可以改动其中选项或补充内容。

```text
1 确认：接受 Stage21 Methods 草稿作为冻结证据基础；缺失字段保留 AUTHOR ACTION REQUIRED。
2 确认：伦理/二次使用措辞待机构或伦理负责人最终确认。
3 确认：GSE326201、GSE282701 原始论文与 Zotero matching 待我人工完成。
4 确认：GEO、GSE189175、ICGC-LIRI-JP 的 dataset-specific licence/使用措辞待我最终确认。
5 确认：repository acknowledgement、访问日期、代码仓库永久标识和 Data Availability 文字待补充。
6 确认：若无 CopyKAT 历史精确版本证据，保留 historical_exact_version_not_recoverable。
7 确认：若无 Monocle3、Slingshot、scTenifoldKnk 历史精确版本证据，保留 historical_exact_version_not_recoverable。
8 确认：我将核对 SCENIC worker 数冲突后再冻结最终 Methods。
9 补充：cell-level FindMarkers/Wilcoxon 完整全局 DEG 输出路径为【请填写；没有则写“未找到”】。
10 确认：Stage21 当前不自动进入 Stage22。
"""


def build_handoff() -> str:
    return r"""# Stage21 交接记录

- 当前阶段：`21-SCI生信研究Methods撰写器`。
- 当前状态：`METHODS_DRAFT_GENERATED_PENDING_USER_CONFIRMATION`。
- 进入授权：用户已明确要求正式进入 Stage21，并限定为冻结证据 Methods 撰写。
- Stage19：`STAGE19_CLOSED_WITH_LIMITATIONS`，不重新打开。
- Stage20：Discussion 已作为冻结输入；本阶段不改写 Discussion。
- Stage20.5：Public Data Compliance Audit 为 `PASS_WITH_LIMITATION`；伦理/二次使用为 `MANUAL_CONFIRMATION_REQUIRED`；Data Availability 为 `MANUSCRIPT_READY_WITH_PLACEHOLDERS`。

## 本阶段已生成

- `输出/Methods草稿.md`
- `输出/软件包与参数表.csv`
- `输出/R实际依赖包版本表.csv`
- `输出/Python实际依赖包版本表.csv`
- `输出/Methods参数追溯表.csv`
- `输出/数据可用性与代码可用性段落.md`
- `输入/输入材料审计.md`
- `质量核查/质量核查表.csv`
- `质量核查/初稿质量评分表.csv`
- `质量核查/质量核查报告.md`
- `需要人工核查.md`
- `过程记录/stage21_input_manifest.tsv`
- `过程记录/stage21_run_record.json`

## 仍未闭合

- 历史 Monocle3、Slingshot、CopyKAT、scTenifoldKnk 以及部分其他依赖的精确版本；不得猜测。
- canonical SCENIC worker-count 记录冲突，需要作者核对。
- TCGA/ICGC cache release/access date、ICGC codebook 语义和最终 dataset-specific licence/waiver wording。
- 伦理/二次使用措辞、最终 Zotero matching、目标期刊 acknowledgement、代码仓库永久标识和 Data Availability placeholders。
- cell-level traditional DEG 是否有完整全局输出；当前 Methods 不把它写成已完成结果。

## 禁止事项

- 不重跑任何生物学分析。
- 不修改 Figure 1–8、Results、Discussion、Stage19 gate 或既有 source-of-truth 文件。
- 不升级 GSE326201、GSE282701、ICGC OS 或 Figure 8 的证据等级。
- 不自动进入 `22-SCI生信研究Introduction撰写器`；必须等待作者逐项确认 Stage21 草稿和人工核查清单。

## 下一阶段

理论下一阶段为 `22-SCI生信研究Introduction撰写器`，当前状态为 `NOT_ENTERED`，须经作者确认后再启动。
"""


def build_qa_report(rows: list[dict[str, object]]) -> str:
    hard = [row for row in rows if row["status"] == "BLOCKED_FOR_FINAL_MANUSCRIPT"]
    blocked = [row for row in rows if row["status"] == "BLOCKED_FOR_FINAL_MANUSCRIPT"]
    return f"""# Stage21 Methods 质量核查报告

## 结论

- Stage21 草稿生成：`PASS_WITH_LIMITATIONS`。
- 最终投稿 Methods 门禁：`BLOCKED_FOR_FINAL_MANUSCRIPT`，原因是若干作者控制的 provenance、伦理/许可和历史精确版本字段仍未闭合。
- 本阶段未重跑生物学分析，未修改 Figure 1–8、Results、Discussion 或 Stage19 状态。
- 自动交接到 Stage22：`NOT_ENTERED`。

## 核查统计

- 质量核查项目：{len(rows)}。
- 需要作者确认或补充的质量项目：{sum(row['status'] in {'PASS_WITH_LIMITATION', 'MANUAL_CONFIRMATION_REQUIRED', 'PARTIALLY_CLOSED', 'BLOCKED_FOR_FINAL_MANUSCRIPT'} for row in rows)}。
- Methods 最终投稿硬阻断记录：{len(blocked)}。
- 质量表中触发最终投稿硬阻断的条目：{len(hard)}。

## 已通过的部分

- Stage19/20/20.5 入口授权与边界读取。
- 主数据、参考数据、raw-count/QC/scVI/scanVI、CellRank、网络、外部复现和 Figure 8 方法链有上游证据。
- GSE326201、GSE282701、ICGC OS 和 Figure 8 的证据等级在 Methods 中保持原状。
- 新文件使用独立 Stage21 目录，不覆盖上游结果。

## 必须保留的限制

1. `historical_exact_version_not_recoverable`：无法从现有证据恢复的历史版本不以当前版本替代。
2. CopyKAT、Monocle3、Slingshot 和 scTenifoldKnk 的历史精确版本需作者补充或保留不可恢复标记。
3. canonical SCENIC worker-count 记录需要核对。
4. TCGA/ICGC release/access date、ICGC codebook、dataset-specific licence/waiver、伦理措辞和 repository identifier 需人工补齐。
5. 传统 cell-level DEG 全局完成状态未确认，草稿仅按预设方案记录。

## 科学边界

- 三轴是部分有序候选架构，不是严格因果级联。
- 计算网络、CNV-supported malignant-like state、轨迹和外部 recurrence 不等同于直接因果。
- ICGC OS = `ESTIMABLE_BUT_NOT_VALIDATED`，Supplementary/Extended Data only。
- Figure 8 = `EXTENDED_DATA_ONLY`，不写疗效、安全性或临床可操作性。
"""


def main() -> None:
    for folder in (OUT, INPUT, LOG, QA):
        folder.mkdir(parents=True, exist_ok=True)

    input_rows = build_input_manifest()
    write_csv(LOG / "stage21_input_manifest.tsv", input_rows, list(input_rows[0].keys()), delimiter="\t")

    r_rows, py_rows, software_manifest = build_dependency_outputs()
    write_csv(OUT / "R实际依赖包版本表.csv", r_rows, list(r_rows[0].keys()))
    write_csv(OUT / "Python实际依赖包版本表.csv", py_rows, list(py_rows[0].keys()))

    parameter_rows = methods_parameter_rows()
    write_csv(OUT / "Methods参数追溯表.csv", parameter_rows, list(parameter_rows[0].keys()))

    combined_rows: list[dict[str, object]] = []
    for row in r_rows:
        combined_rows.append(
            {
                "record_type": "R_dependency",
                "module_or_figure": "code_scan",
                "analysis_role": row["package_class"],
                "parameter_or_package": row["package"],
                "recorded_value": f"current={row['current_version']}; status={row['current_status']}; calls={row['call_types']}",
                "source_path": row["source_files"],
                "status": row["closure_status"],
                "historical_version_status": row["historical_version_status"],
                "notes": f"historical_versions_seen={row['historical_versions_seen']}; evidence={row['historical_evidence_files']}",
            }
        )
    for row in py_rows:
        combined_rows.append(
            {
                "record_type": "Python_dependency",
                "module_or_figure": "code_scan",
                "analysis_role": row["roles"],
                "parameter_or_package": row["distribution"],
                "recorded_value": f"active={row['active_python_version']}; project_venv={row['project_venv_version']}; active_status={row['active_status']}",
                "source_path": row["source_files"],
                "status": row["closure_status"],
                "historical_version_status": row["historical_version_status"],
                "notes": f"module={row['module']}; historical_versions_seen={row['historical_versions_seen']}; evidence={row['historical_evidence_files']}",
            }
        )
    for row in parameter_rows:
        combined_rows.append(
            {
                "record_type": "method_parameter",
                "module_or_figure": row["module_or_figure"],
                "analysis_role": row["analysis_role"],
                "parameter_or_package": row["parameter_or_item"],
                "recorded_value": row["recorded_value"],
                "source_path": row["source_path"],
                "status": row["status"],
                "historical_version_status": row["provenance_status"],
                "notes": row["notes"],
            }
        )
    write_csv(OUT / "软件包与参数表.csv", combined_rows, [
        "record_type", "module_or_figure", "analysis_role", "parameter_or_package", "recorded_value", "source_path", "status", "historical_version_status", "notes"
    ])

    write_text(INPUT / "输入材料审计.md", build_input_audit(), watermark=True)
    write_text(OUT / "Methods草稿.md", METHODS_TEXT, watermark=True)
    write_text(OUT / "数据可用性与代码可用性段落.md", DATA_CODE_TEXT, watermark=True)

    quality_rows = build_quality_rows()
    write_csv(QA / "质量核查表.csv", quality_rows, list(quality_rows[0].keys()))
    score_rows = build_score_rows()
    write_csv(QA / "初稿质量评分表.csv", score_rows, list(score_rows[0].keys()))
    write_text(QA / "质量核查报告.md", build_qa_report(quality_rows), watermark=True)
    write_text(STAGE / "需要人工核查.md", build_human_review(), watermark=True)
    write_text(STAGE / "下一步交接记录.md", build_handoff(), watermark=True)

    missing_history = [
        "Monocle3", "Slingshot", "CopyKAT", "scTenifoldKnk", "AnnotationDbi", "BiocManager", "Matrix", "Seurat", "edgeR", "limma", "monocle3", "msigdbr", "slingshot", "survminer", "cellrank", "arboreto", "ctxcore", "pyscenic", "scipy", "matplotlib", "seaborn"
    ]
    run_record = {
        "created_at_local": now_local(),
        "stage": 21,
        "stage_name": "21-SCI生信研究Methods撰写器",
        "status": "METHODS_DRAFT_GENERATED_PENDING_USER_CONFIRMATION",
        "entry_authorized_by_user": True,
        "stage19_status": "STAGE19_CLOSED_WITH_LIMITATIONS",
        "stage20_status": "DISCUSSION_COMPLETED_FROZEN_INPUT",
        "stage20_5_status": "PASS_WITH_LIMITATION",
        "stage21_scope": "frozen_methods_writing_only",
        "biological_rerun": False,
        "figures_modified": False,
        "results_modified": False,
        "discussion_modified": False,
        "stage19_reopened": False,
        "input_manifest_rows": len(input_rows),
        "input_missing_rows": sum(row["exists"] == "FALSE" for row in input_rows),
        "r_dependency_count": len(r_rows),
        "python_dependency_count": len(py_rows),
        "method_parameter_count": len(parameter_rows),
        "historical_exact_version_not_recoverable_examples": missing_history,
        "scientific_boundaries": {
            "three_axis_architecture": "partial_order_only",
            "gse326201": "Tier 1 exploratory",
            "gse282701": "BLOCKED_PROVENANCE_UNRESOLVED; not used for formal validation",
            "icgc_os": "ESTIMABLE_BUT_NOT_VALIDATED; Supplementary/Extended Data only",
            "figure8": "EXTENDED_DATA_ONLY",
            "adult_hcc_experiment": "future_plan_no_completed_result",
        },
        "manual_confirmation_required": True,
        "next_stage": "22-SCI生信研究Introduction撰写器",
        "auto_handoff": False,
        "output_files": [
            rel(OUT / "Methods草稿.md"),
            rel(OUT / "软件包与参数表.csv"),
            rel(OUT / "R实际依赖包版本表.csv"),
            rel(OUT / "Python实际依赖包版本表.csv"),
            rel(OUT / "Methods参数追溯表.csv"),
            rel(OUT / "数据可用性与代码可用性段落.md"),
            rel(INPUT / "输入材料审计.md"),
            rel(QA / "质量核查表.csv"),
            rel(QA / "初稿质量评分表.csv"),
            rel(QA / "质量核查报告.md"),
            rel(STAGE / "需要人工核查.md"),
            rel(STAGE / "下一步交接记录.md"),
        ],
        "software_environment_manifest_source": rel(ROOT / "metadata/reproducibility/software_environment_manifest.json"),
    }
    write_text(LOG / "stage21_run_record.json", json.dumps(run_record, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
