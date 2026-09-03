# Reproducibility Notes

## Release boundary

The release scope is source-first: the top-level scripts/ tree, selected manuscript-writing source scripts, root repository metadata, and the Stage22 archival closure package. The following remain outside Git tracking: raw or restricted human data, FASTQ/BAM and other large objects, data/, figures/, metadata/, reports/, tmp/, local virtual environments, caches, .Rhistory, generated plots, Zotero local databases/attachments, and credential-bearing files.

## Reproduction principle

Reproduction begins by checking out the verified immutable GitHub Release tag v1.0.0, confirming the commit and manifest hashes, obtaining public datasets from their authoritative accession/repository under the applicable terms, and running only the explicitly documented workflow entry point with its recorded parameters, seed, input manifest, output directory, and software environment. A release tag is a provenance anchor; it is not evidence that every historical runtime can be reconstructed.

## Recorded environment

Project records report Windows, Python 3.11.5, pytest 9.0.3, and R 4.5.0. They also record key Python packages anndata 0.12.16, scanpy 1.11.5, scvi-tools 1.3.3, torch 2.4.1+cu121, numpy 1.26.4, pandas 2.2.2, and statsmodels 0.14.6, plus R packages data.table 1.18.4, survival 3.8.6, jsonlite 2.0.0, ggplot2 4.0.3, patchwork 1.3.2, and ragg 1.5.2. Several exact historical package versions remain unavailable and are recorded as historical_exact_version_not_recoverable.

## Scientific evidence boundaries

GSE326201 is Tier 1 exploratory; GSE282701 remains BLOCKED_PROVENANCE_UNRESOLVED; ICGC OS remains ESTIMABLE_BUT_NOT_VALIDATED and supplementary/Extended Data only; Figure 8 remains EXTENDED_DATA_ONLY. These labels are part of the reproducibility contract.

## Stage22 non-rerun statement

This repository closure used Git, GitHub/Zenodo metadata inspection, filesystem inspection, and documentation generation. It did not rerun scVI/scanVI, CopyKAT, trajectory, SCENIC, CellOracle, scTenifoldKnk, external validation, or any other biological analysis. It did not modify Figure 1–8, Results, Discussion, or Stage19.
