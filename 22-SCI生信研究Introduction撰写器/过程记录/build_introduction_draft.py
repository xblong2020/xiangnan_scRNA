from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
STAGE = ROOT / "22-SCI生信研究Introduction撰写器"
WATERMARK = (
    "> 版权声明：本文件由杨师兄原创“研究型论文 Skill 系统”生成。  \n"
    "> 未经书面授权，禁止复制、传播、改编、转售、商用或用于第三方交付。  \n"
    "> 授权请联系杨师兄。\n"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def write_md(path: Path, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(WATERMARK + "\n" + body.rstrip() + "\n", encoding="utf-8")


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def input_materials() -> list[tuple[str, Path, str]]:
    return [
        (
            "research_question",
            ROOT / "04-SCI生信研究问题与疾病背景构建器/输出/生信研究问题框架.md",
            "Required input: disease background, research question and frozen boundaries",
        ),
        (
            "mechanism_map",
            ROOT / "04-SCI生信研究问题与疾病背景构建器/输出/疾病背景与机制线索表.csv",
            "Required input: axis-specific biological roles and limitations",
        ),
        (
            "literature_scope",
            ROOT / "06-SCI生信文献筛选与证据提取器/输出/高质量对标文献清单.md",
            "Required input: screened literature scope and access status",
        ),
        (
            "literature_evidence_summary",
            ROOT / "06-SCI生信文献筛选与证据提取器/输出/全文筛选证据总结.md",
            "Required input: literature evidence limitations and citation boundary",
        ),
        (
            "pipeline_comparator",
            ROOT / "06-SCI生信文献筛选与证据提取器/输出/分析流程对标证据表.csv",
            "Required input: methods and external-validation comparator evidence",
        ),
        (
            "research_gaps",
            ROOT / "09-SCI生信研究空白识别器/输出/生信研究空白矩阵.csv",
            "Required input: prioritized research gaps",
        ),
        (
            "innovation_candidates",
            ROOT / "09-SCI生信研究空白识别器/输出/创新点候选清单.md",
            "Required input: frozen innovation and claim boundaries",
        ),
        (
            "priority_questions",
            ROOT / "09-SCI生信研究空白识别器/输出/可分析问题优先级.md",
            "Required input: P0/P1/P2 priority and most executable question",
        ),
        (
            "discussion_context",
            ROOT / "20-SCI生信研究Discussion撰写器/输出/Discussion草稿.md",
            "Optional input: frozen Discussion interpretation and limitations",
        ),
        (
            "limitations",
            ROOT / "20-SCI生信研究Discussion撰写器/输出/局限性与未来验证.md",
            "Optional input: limitations and future validation boundary",
        ),
        (
            "discussion_handoff",
            ROOT / "20-SCI生信研究Discussion撰写器/下一步交接记录.md",
            "Optional input: Stage20 handoff record",
        ),
        (
            "methods_draft",
            ROOT / "21-SCI生信研究Methods撰写器/输出/Methods草稿.md",
            "Required input: frozen Methods wording and provenance limits",
        ),
        (
            "methods_handoff",
            ROOT / "21-SCI生信研究Methods撰写器/下一步交接记录.md",
            "Required input: Stage21 technical closure and external pending items",
        ),
        (
            "stage21_gate",
            ROOT / "21-SCI生信研究Methods撰写器/final_closure_audit_v1/06_stage21_gate/STAGE21_FINAL_GATE.json",
            "Required input: Stage21 final gate",
        ),
        (
            "author_confirmation",
            ROOT / "21-SCI生信研究Methods撰写器/final_closure_audit_v1/09_logs/STAGE21_AUTHOR_CONFIRMATION_RECORD.md",
            "Required input: author-confirmed conditional Stage22 authorization",
        ),
        (
            "public_data_record",
            ROOT / "metadata/public_data_accession_version_license.md",
            "Required input: public-data accession, version and access record; final wording remains pending",
        ),
    ]


def build_input_audit() -> tuple[list[dict[str, str]], str]:
    rows: list[dict[str, str]] = []
    for role, path, description in input_materials():
        exists = path.is_file()
        rows.append(
            {
                "role": role,
                "path": str(path),
                "relative_path": rel(path) if exists else "",
                "exists": "TRUE" if exists else "FALSE",
                "file_size_bytes": str(path.stat().st_size) if exists else "",
                "sha256": sha256(path) if exists else "",
                "description": description,
                "status": "AVAILABLE" if exists else "MISSING_OR_NOT_PROVIDED",
            }
        )
    missing = [row["role"] for row in rows if row["exists"] != "TRUE"]
    audit = [
        "# Stage22 输入材料审计",
        "",
        "- Stage22 status: `STAGE22_MANUSCRIPT_INTEGRATION_IN_PROGRESS`.",
        "- Stage21 status: `STAGE21_CLOSED_WITH_PENDING_EXTERNAL_CONFIRMATIONS`.",
        "- Stage19 status: `STAGE19_CLOSED_WITH_LIMITATIONS`; Stage19 was not reopened.",
        "- Scope: manuscript integration and Introduction drafting from frozen evidence only.",
        "- Biological rerun: `FALSE`.",
        "- Figure 1–8 modification: `FALSE`.",
        "- Results/Discussion modification: `FALSE`.",
        "- Stage20.5: author-confirmed `PASS_WITH_LIMITATION` for citation/licence/access/repository acknowledgement, `MANUAL_CONFIRMATION_REQUIRED` for ethics/secondary use, and `MANUSCRIPT_READY_WITH_PLACEHOLDERS` for Data Availability. No separate Stage20.5 directory was assumed.",
        "",
        "## Decision-question status",
        "",
        "Stage22-specific A/B/C/D answers were not newly collected in this run. The draft uses a provisional research-objective-led structure, high-impact cautious density, mechanism-exploration objective and partial-evidence integration profile inferred from the frozen project context. These are drafting defaults, not new author confirmations.",
        "",
        "## Input completeness",
        "",
        f"- Materials checked: `{len(rows)}`.",
        f"- Missing materials: `{len(missing)}` ({', '.join(missing) if missing else 'none'}).",
        "- Missing or uncertain external facts remain in `需要人工核查.md`; they are not filled from inference.",
        "",
        "## Evidence boundary carried forward",
        "",
        "- GSE326201: `Tier 1 exploratory`.",
        "- GSE282701: `BLOCKED_PROVENANCE_UNRESOLVED`; not formal validation.",
        "- ICGC OS: `ESTIMABLE_BUT_NOT_VALIDATED`; exploratory donor-level continuous analysis and Supplementary/Extended Data only.",
        "- Figure 8: `EXTENDED_DATA_ONLY`; hypothesis-generating transcriptomic rescue candidates only.",
        "- Adult HCC cell experiments: future plan; no completed experimental result.",
        "- Historical exact software versions not recoverable remain explicitly unrecoverable.",
        "",
        "## Audit table",
        "",
        "See `过程记录/stage22_input_manifest.tsv` for size, timestamp-independent SHA-256 and availability status.",
    ]
    return rows, "\n".join(audit)


INTRODUCTION = r'''# Introduction

## A state-based view of hepatocellular carcinoma development

Hepatocellular carcinoma (HCC) develops within a liver shaped by chronic injury, regenerative pressure, metabolic remodeling and heterogeneous malignant cell states. This context makes hepatocyte transformation difficult to describe as a single molecular event: the same tissue can contain mature hepatocytes, injury-associated or regenerative states, transitional populations and cells with malignant-like transcriptional features. A useful disease model therefore needs to resolve how lineage identity, stress adaptation and malignant-state plasticity are distributed across cellular states while preserving the etiological and patient-specific context in which those states arise. [REF-PENDING: HCC disease background and chronic injury; verify against the project’s Zotero records before submission]

Hepatocyte identity is a biologically meaningful reference axis in this setting. HNF4A and PPARA are associated with differentiated hepatocyte transcriptional and metabolic programmes, and attenuation of these programmes can accompany loss of mature cellular features during injury, dedifferentiation or neoplastic progression. The direction and functional meaning of such attenuation may differ between regenerative and malignant contexts, however. In particular, expression-level identity loss does not by itself establish that HNF4A or PPARA is the initiating driver of malignant transformation. The project’s literature screen contains HNF4A-linked HCC mechanism studies and hepatocyte-focused single-cell/bulk analyses, but the exact article-to-claim mapping remains pending manual Zotero verification. [REF-PENDING: `SQDMHNEF`; verify original article, DOI/PMID and permitted claim]

An injury-responsive stress layer provides a second source of complexity. AP-1 family factors, CEBPB and EGR1 can participate in inflammatory, regenerative and immediate-early transcriptional responses, and those responses may overlap with the plasticity of malignant cells. Consequently, activation of JUN, FOS, JUND, CEBPB or EGR1 should not be equated automatically with cancer-specific progression. The relevant question is whether a stress-transition programme can be separated from generic dissociation or inflammatory responses and positioned relative to hepatocyte identity loss and malignant-like states across datasets with different causes of liver disease. This distinction is particularly important for single-cell and single-nucleus studies, in which tissue handling and cell-state composition can influence the observed stress landscape. [REF-PENDING: project methods/context records and relevant Zotero entries; verify exact primary citations]

Single-cell and single-nucleus transcriptomics provide an opportunity to examine this problem at the level of cellular state rather than only at the level of bulk tissue averages. CNV inference, trajectory reconstruction and fate-probability models can add complementary structure to expression-defined states, while regulon and prior-network approaches can prioritize candidate transcriptional relationships. Yet these layers represent different kinds of evidence. Inferred CNV is not equivalent to orthogonally measured genomic alteration; pseudotime is not physical time; SCENIC/cisTarget, CellOracle and scTenifoldKnk provide network or perturbation hypotheses rather than direct binding or functional causality. An evidence-layered design is therefore needed to prevent a large collection of computational outputs from being misread as a linear mechanism. [REF-PENDING: `Y9E8YGYC`, `52VCXU3P`; verify bibliographic mapping and scope]

The project’s screened literature indicates that HCC studies increasingly combine single-cell atlases, trajectory analysis, transcription-factor networks and TCGA/ICGC or other external cohorts. This methodological convergence creates an opportunity, but it also leaves several unresolved issues. First, identity erosion, injury-associated stress transition and malignant-state plasticity are often analysed in separate frameworks rather than as partially overlapping axes. Second, cross-dataset differences in disease aetiology, platform, tissue handling, sample structure and clinical annotation can make an apparently stable programme context-dependent. Third, external expression recurrence and clinical prediction are frequently discussed together even though they answer different questions. Finally, computational perturbation can prioritize candidates, but it cannot replace perturbation, rescue or orthogonal chromatin and protein evidence. These observations motivate an integration strategy that treats reproducibility, heterogeneity and evidence level as part of the biological question. [REF-PENDING: `8YW6KXL4`, `BJ9R7EE6`, `GYQ8H5TX`; verify original papers and exact claims]

## Study objective

Here we develop a cross-dataset single-cell/single-nucleus framework to examine hepatocyte states spanning reference or injury-associated populations and a CNV-supported malignant-like state. The central question is whether three transcriptional-regulatory programmes can be resolved and compared within a common state space: (i) HNF4A/PPARA-associated hepatocyte identity loss, (ii) AP-1/CEBPB/EGR1-associated stress transition, and (iii) a later SOX4-associated malignant-state or plasticity programme. We use trajectory and CellRank analyses to estimate relative state positioning, canonical SCENIC/cisTarget to assess regulon-level support, CellOracle and scTenifoldKnk to generate directional perturbation hypotheses, and independent bulk or single-cell resources to evaluate axis-specific recurrence. The analyses are organized as complementary evidence layers, with sample- and patient-level limitations retained wherever they affect interpretation.

The intended contribution is a partially ordered, overlapping three-axis architecture rather than a proven linear cascade. The framework is designed to test whether identity loss, stress transition and SOX4-associated stabilization show reproducible but asymmetric relationships, and whether those relationships remain interpretable when CNV inference, fate probability, dataset effects and external-cohort heterogeneity are considered separately. External data are used for population-level corroboration and axis-specific recurrence; they are not used here to construct a new clinical prediction model. Computational network outputs are treated as candidate regulatory evidence, and any future adult HCC cell experiments are planned as additional validation rather than represented as completed results.

Finally, the framework can generate exploratory transcriptomic-rescue hypotheses for later testing. Such candidates remain subject to an independent mechanism gate, orthogonal evidence and experimental rescue. Figure 8 therefore remains `EXTENDED_DATA_ONLY` in the project’s evidence ledger, and no efficacy, safety or clinical-actionability claim follows from the computational ranking. Likewise, ICGC survival information remains `ESTIMABLE_BUT_NOT_VALIDATED` and restricted to exploratory Supplementary/Extended Data use until the outstanding provenance and clinical-codebook requirements are resolved.

'''


LOGIC_CHAIN = r'''# 引言逻辑链

本稿采用“疾病问题 → 三轴生物学张力 → 单细胞/网络方法机会 → 文献与证据缺口 → 研究目标 → 贡献与边界”的递进结构。Stage22 专属四组选项本轮未重新收集，以下结构属于基于冻结项目上下文的临时写作配置，待作者审阅。

| 顺序 | 引言功能 | 本项目对应内容 | 允许的证据强度 | 主要来源 |
|---|---|---|---|---|
| 1 | 疾病问题 | HCC发生于慢性损伤、再生、代谢重塑和恶性样异质性共存的环境 | 疾病背景；具体文献待匹配 | 04阶段问题框架；06阶段文献库 |
| 2 | 身份轴 | HNF4A/PPARA-associated hepatocyte identity loss | identity-associated state programme；不写启动因果 | 04机制线索表；06文献对标 |
| 3 | 应激轴 | AP-1/CEBPB/EGR1 stress-transition | injury/stress-transition candidate；不等同恶性驱动 | 04机制线索表；09空白矩阵 |
| 4 | 恶性样轴 | CNV-supported malignant-like state与later SOX4-associated plasticity | computationally supported malignant-like state；SOX4为candidate | 04机制线索表；09创新点清单 |
| 5 | 方法机会 | 单细胞/单核、CNV、trajectory、CellRank、SCENIC/cisTarget、CellOracle、scTenifoldKnk | complementary evidence layers；不把方法收敛写成因果 | 06流程对标；21 Methods |
| 6 | 核心缺口 | 三轴整合、跨队列异质性、患者/样本结构、功能验证和provenance | 部分有序问题可检验；严格级联未证实 | 09空白矩阵；20 Discussion |
| 7 | 研究目标 | 在统一状态空间中评估三轴相对位置、方向性候选和外部recurrence | population-level corroboration与hypothesis generation | 04问题框架；09可分析问题优先级 |
| 8 | 边界 | GSE326201 Tier 1 exploratory；GSE282701 blocked；ICGC OS未验证；Figure8 Extended Data only | 保留pending/blocked状态 | Stage19/21 closure；作者确认 |

## 需要避免的引言表述

- `proven linear causal cascade`、`master regulator`、`independent prognostic validation`：当前证据不支持。
- `therapeutic efficacy`、`safety`、`clinical actionability`：Figure 8仅为探索性转录救援候选。
- `completed adult HCC cell validation`：实验尚未完成。
- 将GSE326201或GSE282701写成同等等级的正式独立验证：分别受Tier 1 exploratory和raw-count provenance unresolved限制。
- 将ICGC OS写成正式外部临床预测验证：当前状态为`ESTIMABLE_BUT_NOT_VALIDATED`。

'''


def citation_rows() -> list[dict[str, str]]:
    return [
        {
            "citation_id": "REF-PENDING-01",
            "project_reference_key": "SQDMHNEF",
            "title_or_source_description": "Targeting FDFT1 Reduces Cholesterol and Bile Acid Production and Delays Hepatocellular Carcinoma Progression Through the HNF4A/ALDOB/AKT1 Axis.",
            "evidence_use": "HNF4A-linked hepatocyte identity/metabolic context; use only after article-level verification",
            "source_file": "06-SCI生信文献筛选与证据提取器/输出/分析流程对标证据表.csv",
            "verification_status": "PENDING_REFERENCE_CONFIRMATION",
            "allowed_claim": "Contextual association only; no direct causal claim for this project",
            "manual_action": "Match original article, DOI/PMID and Zotero item before manuscript citation",
        },
        {
            "citation_id": "REF-PENDING-02",
            "project_reference_key": "Y9E8YGYC",
            "title_or_source_description": "Single-Cell Characterization of Terminal States and State-Specific Transcriptional Regulatory Networks in Hepatocellular Carcinoma.",
            "evidence_use": "HCC single-cell terminal-state, trajectory and regulatory-network methodological context",
            "source_file": "06-SCI生信文献筛选与证据提取器/输出/分析流程对标证据表.csv",
            "verification_status": "PENDING_REFERENCE_CONFIRMATION",
            "allowed_claim": "Methodological precedent; does not validate the present three-axis model",
            "manual_action": "Match original article, DOI/PMID and Zotero item before manuscript citation",
        },
        {
            "citation_id": "REF-PENDING-03",
            "project_reference_key": "52VCXU3P",
            "title_or_source_description": "Large-scale single-cell analysis and in silico perturbation reveal dynamic evolution of HCC: from initiation to therapeutic targeting.",
            "evidence_use": "Single-cell, trajectory, SCENIC/cisTarget and in silico perturbation comparator",
            "source_file": "06-SCI生信文献筛选与证据提取器/输出/分析流程对标证据表.csv",
            "verification_status": "PENDING_REFERENCE_CONFIRMATION",
            "allowed_claim": "Methodological comparator; computational perturbation remains hypothesis-generating",
            "manual_action": "Match original article, DOI/PMID and Zotero item before manuscript citation",
        },
        {
            "citation_id": "REF-PENDING-04",
            "project_reference_key": "8YW6KXL4",
            "title_or_source_description": "Deep dissection of stemness-related hierarchies in hepatocellular carcinoma.",
            "evidence_use": "HCC state hierarchy and external-cohort methodological context",
            "source_file": "06-SCI生信文献筛选与证据提取器/输出/分析流程对标证据表.csv",
            "verification_status": "PENDING_REFERENCE_CONFIRMATION",
            "allowed_claim": "Contextual/methodological precedent only",
            "manual_action": "Match original article, DOI/PMID and Zotero item before manuscript citation",
        },
        {
            "citation_id": "REF-PENDING-05",
            "project_reference_key": "BJ9R7EE6",
            "title_or_source_description": "Single-cell multi-omics reveals DUSP9 as a key regulator of cancer stemness and a potential therapeutic target in hepatocellular carcinoma.",
            "evidence_use": "HCC single-cell multi-omics, SCENIC and external validation comparator",
            "source_file": "06-SCI生信文献筛选与证据提取器/输出/分析流程对标证据表.csv",
            "verification_status": "PENDING_REFERENCE_CONFIRMATION",
            "allowed_claim": "Comparator only; does not establish a therapeutic claim for this project",
            "manual_action": "Match original article, DOI/PMID and Zotero item before manuscript citation",
        },
        {
            "citation_id": "REF-PENDING-06",
            "project_reference_key": "GYQ8H5TX",
            "title_or_source_description": "Deciphering the Oncogenic Landscape of Hepatocytes Through Integrated Single-Nucleus and Bulk RNA-Seq of Hepatocellular Carcinoma.",
            "evidence_use": "Integrated single-nucleus/bulk HCC comparator",
            "source_file": "06-SCI生信文献筛选与证据提取器/输出/分析流程对标证据表.csv",
            "verification_status": "PENDING_REFERENCE_CONFIRMATION",
            "allowed_claim": "Comparator only; exact bibliographic and figure-level scope requires verification",
            "manual_action": "Match original article, DOI/PMID and Zotero item before manuscript citation",
        },
        {
            "citation_id": "PROJECT-EVIDENCE-04",
            "project_reference_key": "PROJECT-04",
            "title_or_source_description": "Frozen disease background and mechanism-line evidence table",
            "evidence_use": "Defines the three axes and their explicit claim limits",
            "source_file": "04-SCI生信研究问题与疾病背景构建器/输出/疾病背景与机制线索表.csv",
            "verification_status": "PROJECT_INTERNAL_SOURCE_VERIFIED",
            "allowed_claim": "Internal design evidence, not an external citation",
            "manual_action": "Keep as provenance; replace manuscript citation placeholders with verified primary literature",
        },
        {
            "citation_id": "PROJECT-EVIDENCE-09",
            "project_reference_key": "PROJECT-09",
            "title_or_source_description": "Frozen research-gap matrix and innovation candidates",
            "evidence_use": "Defines the partially ordered architecture and evidence boundaries",
            "source_file": "09-SCI生信研究空白识别器/输出/生信研究空白矩阵.csv",
            "verification_status": "PROJECT_INTERNAL_SOURCE_VERIFIED",
            "allowed_claim": "Internal study rationale, not an external citation",
            "manual_action": "Use for audit trail; do not cite as literature",
        },
        {
            "citation_id": "PUBLIC-DATA-01",
            "project_reference_key": "PUBLIC_DATA_RECORD",
            "title_or_source_description": "Public data accession, version, download and reuse record",
            "evidence_use": "Identifies public data provenance and unresolved licence/ethics wording",
            "source_file": "metadata/public_data_accession_version_license.md",
            "verification_status": "PASS_WITH_LIMITATION; MANUAL_CONFIRMATION_REQUIRED",
            "allowed_claim": "Data-source provenance only; no blanket licence or ethics waiver claim",
            "manual_action": "Complete dataset-specific citation, licence, repository acknowledgement and ethics wording before submission",
        },
    ]


def quality_rows() -> list[dict[str, str]]:
    return [
        {"check_id": "Q01", "check_item": "Required upstream materials located and hashed", "status": "PASS_WITH_LIMITATIONS", "evidence_source": "过程记录/stage22_input_manifest.tsv", "notes": "All required local inputs were available in this run; separate Stage20.5 directory was not assumed.", "hard_blocker": "NO", "action": "Retain the manifest and author confirmation record."},
        {"check_id": "Q02", "check_item": "Disease problem, three-axis rationale, gap and objective form a progression", "status": "PASS", "evidence_source": "输出/引言逻辑链.md; 04; 09", "notes": "Objective-led structure is coherent and explicitly provisional pending author review.", "hard_blocker": "NO", "action": "Author review for emphasis and journal-specific density."},
        {"check_id": "Q03", "check_item": "Introduction does not report new Results or numerical findings", "status": "PASS", "evidence_source": "输出/Introduction草稿.md", "notes": "The text states questions, design and boundaries without introducing Figure results.", "hard_blocker": "NO", "action": "Check against the frozen Results during full integration."},
        {"check_id": "Q04", "check_item": "External citations are traceable to project evidence", "status": "PARTIALLY_CLOSED", "evidence_source": "输出/引用证据表.csv; 06", "notes": "Project keys and titles are recorded, but DOI/PMID/Zotero matching remains pending.", "hard_blocker": "YES", "action": "Manual Zotero matching is required before final manuscript readiness."},
        {"check_id": "Q05", "check_item": "Data provenance and public-data compliance wording are final", "status": "PARTIALLY_CLOSED", "evidence_source": "metadata/public_data_accession_version_license.md; Stage21 author confirmation", "notes": "Citation/licence/access/repository acknowledgement and ethics wording remain author-controlled pending items.", "hard_blocker": "YES", "action": "Complete the external compliance checklist before submission."},
        {"check_id": "Q06", "check_item": "Causal, temporal and clinical claims remain within evidence level", "status": "PASS_WITH_LIMITATIONS", "evidence_source": "输出/Introduction草稿.md; Stage21 final gate", "notes": "Partial order, computational hypotheses, ICGC and Figure8 boundaries are retained.", "hard_blocker": "NO", "action": "Recheck during Results/Discussion/Methods integration."},
        {"check_id": "Q07", "check_item": "Figure8 evidence boundary retained", "status": "PASS", "evidence_source": "输出/Introduction草稿.md; 09 innovation candidates", "notes": "Figure8 remains EXTENDED_DATA_ONLY and hypothesis-generating.", "hard_blocker": "NO", "action": "Do not upgrade during later manuscript integration."},
        {"check_id": "Q08", "check_item": "GSE326201, GSE282701 and ICGC evidence levels retained", "status": "PASS", "evidence_source": "输入/输入材料审计.md; Stage21 final gate", "notes": "No cohort substitution or evidence-tier upgrade was made.", "hard_blocker": "NO", "action": "Keep cohort labels synchronized in later sections."},
        {"check_id": "Q09", "check_item": "Target-journal adaptation is complete", "status": "PENDING", "evidence_source": "User project record", "notes": "Target journal, word limits and exact reference style remain unconfirmed.", "hard_blocker": "NO", "action": "Handle during journal selection and full-text integration."},
        {"check_id": "Q10", "check_item": "Stage22 initial-draft gate permits automatic next-stage handoff", "status": "NO_AUTOMATIC_HANDOFF", "evidence_source": "Q04-Q05; 需要人工核查.md", "notes": "The draft is usable for review but external citation/compliance confirmations remain open.", "hard_blocker": "YES", "action": "Keep Stage22 in progress; do not enter Stage23 automatically."},
    ]


def score_rows() -> list[dict[str, str]]:
    return [
        {"评价项": "证据完整性", "评分0-5": "4", "证据来源": "04/06/09/20/21 input audit", "问题说明": "Rationale and boundaries are traceable; bibliographic matching remains pending.", "是否触发硬阻断": "是", "修正建议": "Complete DOI/PMID/Zotero mapping."},
        {"评价项": "结果-图表一致性", "评分0-5": "4", "证据来源": "Stage21 final gate; frozen Figure boundary record", "问题说明": "No Results or new figure claims were added; full manuscript cross-check remains.", "是否触发硬阻断": "否", "修正建议": "Cross-check during Stage22 full integration."},
        {"评价项": "方法可重复性", "评分0-5": "4", "证据来源": "Stage21 Methods draft and closure audit", "问题说明": "Technical closure accepted; historical exact versions remain unrecoverable by policy.", "是否触发硬阻断": "否", "修正建议": "Preserve explicit unrecoverable-version labels."},
        {"评价项": "逻辑连贯性", "评分0-5": "4", "证据来源": "输出/引言逻辑链.md", "问题说明": "Disease-to-gap-to-objective sequence is explicit.", "是否触发硬阻断": "否", "修正建议": "Author review of emphasis and length."},
        {"评价项": "路线专属规范符合度", "评分0-5": "4", "证据来源": "Stage22 skill; Stage21 boundary records", "问题说明": "Single-cell, network and external-validation roles are separated; final citation wording remains open.", "是否触发硬阻断": "是", "修正建议": "Complete external compliance and citation review."},
        {"评价项": "目标期刊适配度", "评分0-5": "2", "证据来源": "User project record", "问题说明": "Nature or equivalent is a preference; final journal and format are not fixed.", "是否触发硬阻断": "否", "修正建议": "Apply journal-specific style only after target selection."},
        {"评价项": "夸大或编造风险", "评分0-5": "5", "证据来源": "输出/Introduction草稿.md; boundary audit", "问题说明": "Causal, prognostic, efficacy, safety and completed-experiment claims are explicitly excluded.", "是否触发硬阻断": "否", "修正建议": "Maintain the same language gate in later sections."},
        {"评价项": "是否允许进入下一写作环节", "评分0-5": "0", "证据来源": "需要人工核查.md; Q04-Q05", "问题说明": "Stage22 remains in progress and must not auto-handoff to Stage23.", "是否触发硬阻断": "是", "修正建议": "Obtain itemized author confirmation after review."},
    ]


def main() -> None:
    dirs = [STAGE / name for name in ["输入", "输出", "过程记录", "质量核查"]]
    for directory in dirs:
        directory.mkdir(parents=True, exist_ok=True)

    rows, audit = build_input_audit()
    write_md(STAGE / "输入/输入材料审计.md", audit)
    write_tsv(
        STAGE / "过程记录/stage22_input_manifest.tsv",
        list(rows[0].keys()),
        rows,
    )

    write_md(STAGE / "输出/Introduction草稿.md", INTRODUCTION)
    write_md(STAGE / "输出/引言逻辑链.md", LOGIC_CHAIN)
    write_csv(
        STAGE / "输出/引用证据表.csv",
        [
            "citation_id",
            "project_reference_key",
            "title_or_source_description",
            "evidence_use",
            "source_file",
            "verification_status",
            "allowed_claim",
            "manual_action",
        ],
        citation_rows(),
    )
    write_csv(
        STAGE / "质量核查/质量核查表.csv",
        ["check_id", "check_item", "status", "evidence_source", "notes", "hard_blocker", "action"],
        quality_rows(),
    )
    write_csv(
        STAGE / "质量核查/初稿质量评分表.csv",
        ["评价项", "评分0-5", "证据来源", "问题说明", "是否触发硬阻断", "修正建议"],
        score_rows(),
    )

    human_review = r'''# Stage22 需要人工核查

Stage22 当前状态为 `STAGE22_MANUSCRIPT_INTEGRATION_IN_PROGRESS`。本稿是基于冻结材料的 Introduction 初稿，外部事实和最终投稿合规事项仍不能自动确认。

1. `PENDING_REFERENCE_CONFIRMATION`：逐项将 GSE326201、GSE282701、GSE189175、ICGC-LIRI-JP 及引言引用候选与原始论文、DOI/PMID和Zotero key匹配。
2. `MANUAL_CONFIRMATION_REQUIRED`：由机构伦理负责人最终审核公共人类数据二次使用/伦理措辞。
3. `MANUAL_CONFIRMATION_REQUIRED`：确认 GEO、GSE189175、GSE326201、GSE282701、ICGC-LIRI-JP 的最终 licence、repository acknowledgement 和 Data Availability wording；不得从公共可访问性推断统一许可。
4. `PENDING_REPOSITORY_ARCHIVAL`：补充代码仓库 remote、release tag、archival、DOI或其他永久标识；该项阻断 FINAL_MANUSCRIPT_READY，但不阻断 Stage22整合。
5. `PENDING`：确认本 Introduction 的疾病背景、三轴定义、研究空白和目标表述准确，尤其是是否接受“partially ordered, overlapping three-axis architecture”的措辞。
6. `PENDING`：确认临时写作配置（研究目标导向、高影响力谨慎密度、机制探索目标）和Introduction篇幅/语体，或补充修改意见。
7. `PENDING`：保留 `historical_exact_version_not_recoverable`；不得把当前软件版本替代为历史版本。
8. `PENDING`：确认 Stage22 不自动进入 Stage23，直到本轮稿件和人工核查项逐项处理。

## 固定科学边界

- Stage19 保持 `STAGE19_CLOSED_WITH_LIMITATIONS`，不重新打开。
- GSE326201 保持 `Tier 1 exploratory`；GSE282701 保持 `BLOCKED_PROVENANCE_UNRESOLVED`。
- ICGC OS 保持 `ESTIMABLE_BUT_NOT_VALIDATED`，仅 Supplementary/Extended Data 探索性使用。
- Figure 8 保持 `EXTENDED_DATA_ONLY`，不写疗效、安全性或临床可操作性。
- 成人 HCC 细胞实验是未来计划，当前没有完成的实验结果。
- 不新增生物学分析、不修改 Figure 1–8、不改变冻结假设、不更换验证队列、不伪造引用/许可/伦理/仓储信息。

下面是基于当前材料的推荐复制回答；你可以直接复制发送，也可以改动其中选项或补充内容。

```text
1 确认：我已逐项核对Introduction草稿中的疾病背景、三轴定义、研究空白和研究目标；需要修改的地方为【无/请填写】。
2 确认：接受“partially ordered, overlapping three-axis architecture”作为当前谨慎表述，不解释为严格级联。
3 补充：引言候选文献的 DOI/PMID/Zotero key 待我人工匹配；完成后提供映射记录。
4 补充：公共数据 licence、repository acknowledgement、Data Availability 和伦理/二次使用措辞待我最终审核。
5 确认：代码仓库永久标识仍为 PENDING_REPOSITORY_ARCHIVAL，不将当前状态写成最终可投稿。
6 确认：GSE326201、GSE282701、ICGC OS、Figure 8 和成人HCC实验边界保持不变。
7 确认：本轮不自动进入Stage23；Stage22继续保持 STAGE22_MANUSCRIPT_INTEGRATION_IN_PROGRESS。
```
'''
    write_md(STAGE / "需要人工核查.md", human_review)

    handoff = r'''# Stage22 交接记录

- 当前阶段：`22-SCI生信研究Introduction撰写器`。
- 当前状态：`STAGE22_MANUSCRIPT_INTEGRATION_IN_PROGRESS`。
- 进入依据：作者已确认 Stage21 `STAGE21_CLOSED_WITH_PENDING_EXTERNAL_CONFIRMATIONS`，并授权条件进入 Stage22。
- 当前范围：仅基于冻结的疾病背景、研究空白、项目内文献证据、Discussion/Methods和合规审计记录生成 Introduction 初稿。
- Stage19：`STAGE19_CLOSED_WITH_LIMITATIONS`，未重新打开。
- Stage21：技术闭合已接受；外部引用、伦理、许可、仓储永久标识和最终 Data Availability 仍为人工事项。
- 生物学分析重跑：`FALSE`。
- Figure 1–8修改：`FALSE`。
- Results/Discussion修改：`FALSE`。
- Stage23自动交接：`FALSE`。
- `FINAL_MANUSCRIPT_READY`: `FALSE`。
- `FINAL_SUBMISSION_READY`: `FALSE`。

## 本阶段产物

- `输入/输入材料审计.md`
- `输出/Introduction草稿.md`
- `输出/引言逻辑链.md`
- `输出/引用证据表.csv`
- `质量核查/质量核查表.csv`
- `质量核查/初稿质量评分表.csv`
- `过程记录/stage22_input_manifest.tsv`
- `需要人工核查.md`

## 仍然未闭合

- `PENDING_REFERENCE_CONFIRMATION`：原始论文与 DOI/PMID/Zotero key 匹配。
- `MANUAL_CONFIRMATION_REQUIRED`：伦理/公共数据二次使用措辞。
- `PASS_WITH_LIMITATION` 待最终人工关闭：公共数据 licence、repository acknowledgement和Data Availability wording。
- `PENDING_REPOSITORY_ARCHIVAL`：代码仓库永久标识。
- `historical_exact_version_not_recoverable`：继续保留，不能猜测。

## 下一阶段

按 Stage22 技能的正常路线，下一阶段为 `23-SCI生信标题摘要关键词撰写器`。本记录不构成自动交接；在作者逐项核查本阶段草稿前，不进入 Stage23。
'''
    write_md(STAGE / "下一步交接记录.md", handoff)

    run_record = {
        "created_at_local": datetime.now().astimezone().isoformat(timespec="seconds"),
        "stage": 22,
        "status": "STAGE22_MANUSCRIPT_INTEGRATION_IN_PROGRESS",
        "decision_questions_status": "NOT_NEWLY_COLLECTED_IN_THIS_RUN",
        "draft_profile": {
            "structure": "provisional_research_objective_led",
            "literature_density": "provisional_high_impact_cautious",
            "research_goal": "mechanism_exploration",
            "evidence_profile": "partial_evidence_with_explicit_pending_items",
            "profile_status": "INFERRED_FROM_FROZEN_PROJECT_CONTEXT_NOT_NEW_AUTHOR_CONFIRMATION",
        },
        "stage19_status": "STAGE19_CLOSED_WITH_LIMITATIONS",
        "stage21_status": "STAGE21_CLOSED_WITH_PENDING_EXTERNAL_CONFIRMATIONS",
        "stage20_5_status": "PASS_WITH_LIMITATION_AND_MANUAL_CONFIRMATIONS_PENDING",
        "stage22_scope": "MANUSCRIPT_INTEGRATION_ONLY",
        "biological_rerun": False,
        "new_cohort": False,
        "new_results": False,
        "figure1_8_modified": False,
        "results_modified": False,
        "discussion_modified": False,
        "stage19_reopened": False,
        "stage23_auto_handoff": False,
        "final_manuscript_ready": False,
        "final_submission_ready": False,
        "evidence_boundaries": {
            "gse326201": "Tier 1 exploratory",
            "gse282701": "BLOCKED_PROVENANCE_UNRESOLVED",
            "icgc_os": "ESTIMABLE_BUT_NOT_VALIDATED; Supplementary/Extended Data only",
            "figure8": "EXTENDED_DATA_ONLY",
            "adult_hcc_experiment": "future_plan_no_completed_result",
        },
        "pending_external_confirmations": [
            "PENDING_REFERENCE_CONFIRMATION",
            "MANUAL_CONFIRMATION_REQUIRED_ETHICS_SECONDARY_USE",
            "MANUAL_CONFIRMATION_REQUIRED_DATASET_LICENCE_ACKNOWLEDGEMENT",
            "PENDING_REPOSITORY_ARCHIVAL",
        ],
        "next_stage": "23-SCI生信标题摘要关键词撰写器",
        "user_confirmation_required_before_next_stage": True,
    }
    (STAGE / "过程记录/stage22_run_record.json").write_text(
        json.dumps(run_record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"stage": 22, "status": run_record["status"], "input_count": len(rows)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
