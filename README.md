# HCC hepatocyte state-transition single-cell analysis code

## Project scope

This repository contains the analysis and reproducibility code for a hepatocellular carcinoma (HCC) single-cell study of hepatocyte state change across a reference/normal-like state, chronic injury or cirrhosis-associated states, transition programmes, and a CNV-supported malignant-like state. The working biological architecture has three partially ordered, overlapping axes:

1. HNF4A/PPARA identity loss;
2. AP-1/CEBPB/EGR1 stress-transition;
3. later SOX4-associated malignant-state stabilization.

The project records are evidence-led. Computational association, trajectory positioning, network displacement, and external recurrence are retained with their documented limitations and are not converted into direct causal or clinical claims.

## Analysis overview

The source tree covers public-data preparation, count/QC harmonisation, scVI/scanVI integration, cell-type annotation, CNV/malignant-state processing, trajectory analysis, CellOracle/SCENIC and scTenifoldKnk network workflows, external validation, and figure-source-data generation. The top-level scripts/ directory is the release source tree. Stage-specific writing and audit records remain in the project workspace; the Stage22 archival closure package is preserved under the Introduction stage directory.

Raw human data, restricted patient-level files, large sequencing objects, derived figure files, and temporary outputs are not redistributed by this code release. Reproduction uses accession-based retrieval and the project manifests where the source repositories and usage terms permit it.

## Repository layout

- scripts/: Python, R, and supporting workflow entry points.
- 22-SCI生信研究Introduction撰写器/code_repository_archival_closure_v1/: Git provenance, release, remote, archival, manuscript wording, QA, and Stage22 gate records.
- data/, figures/, metadata/, reports/, tmp/: protected local research assets or generated outputs; excluded from Git tracking by policy.
- tests/: existing project logic tests; running tests does not substitute for rerunning the biological pipeline.
- REPRODUCIBILITY_NOTES.md: environment and reproduction boundary notes.
- CODE_AVAILABILITY.md: conservative repository-availability wording.
- CITATION.cff: citation metadata with unresolved external fields kept explicit.

## Prerequisites and major workflows

The project evidence records a Windows environment with Python 3.11.5, pytest 9.0.3, and R 4.5.0. Key packages recorded in the project include anndata 0.12.16, scanpy 1.11.5, scvi-tools 1.3.3, torch 2.4.1+cu121, numpy 1.26.4, pandas 2.2.2, statsmodels 0.14.6, data.table 1.18.4, survival 3.8.6, jsonlite 2.0.0, ggplot2 4.0.3, patchwork 1.3.2, and ragg 1.5.2. These are project-recorded environment facts; exact historical versions for every analysis step are not fully recoverable.

The principal workflow order is data intake and download records → preprocessing/QC and batch handling → differential and enrichment analysis → network and candidate analysis → conditional modelling → external validation planning → figure contracts/source data → reproducibility audit → results storyline → manuscript integration. Each expensive biological entry point requires an explicit input, parameter, seed, output, and audit record; this Stage22 closure does not invoke any of them.

## Frozen scientific boundaries

- GSE326201: Tier 1 exploratory cohort; it must not be described as formal Tier 2+ validation.
- GSE282701: BLOCKED_PROVENANCE_UNRESOLVED; it is not an eligible formal validation input.
- ICGC-LIRI-JP overall survival: ESTIMABLE_BUT_NOT_VALIDATED; supplementary/Extended Data use only, with no clinical utility or independent prognostic-validation claim.
- Figure 8: EXTENDED_DATA_ONLY under the objective evidence gate, regardless of the user preference for a main-text candidate; no efficacy, safety, treatment, or clinical-actionability claim.
- Adult HCC cell experiments: future plan only; no completed experimental result is represented here.
- Stage19 remains STAGE19_CLOSED_WITH_LIMITATIONS.
- historical_exact_version_not_recoverable remains explicit wherever the historical runtime cannot be independently reconstructed.

## Scope of this Stage22 closure

This closure performs Git, filesystem, provenance, release-preparation, and documentation work only. No biological analysis was rerun; no new cohort or biological result was created; Figure 1–8, Results, and Discussion were not modified; Stage19 was not reopened; and Stage23 is not entered automatically. The release tag and commit recorded in the closure package describe repository provenance, not a new scientific analysis.

The public GitHub repository is https://github.com/xblong2020/xiangnan_scRNA. The verified remote main publication commit is a7ee1edb8c139e9dca94a3e3702c750701e019e0, and the frozen annotated tag v1.0.0 points to 1c0303049ba629b9b986cb9a7e088384f15f5a87. The GitHub Release page still requires manual creation because the GitHub CLI is unavailable in this environment. Zenodo archival and its permanent identifier remain pending.

## Citation

The analysis code is publicly available at https://github.com/xblong2020/xiangnan_scRNA, release v1.0.0. A Zenodo permanent archival identifier remains pending. After Zenodo returns a real DOI, replace only the pending archive fields with the verified DOI and record URL. See CITATION.cff and the Stage22 closure package for the evidence trail.
