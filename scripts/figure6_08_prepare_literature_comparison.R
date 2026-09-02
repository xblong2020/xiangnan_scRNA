#!/usr/bin/env Rscript

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figure6_common.R"))

paper <- list(
  citation = "Zhang XF et al. Journal of Hepatology. 2026;84(2):370-384.",
  title = "Targeting cell-state plasticity driven by FOXM1/CEBPB axis disrupts developmental heterogeneity and therapeutic resistance in hepatocellular carcinoma",
  doi = "10.1016/j.jhep.2025.09.022", pubmed = "41043722", url = "https://pubmed.ncbi.nlm.nih.gov/41043722/"
)
tbl <- data.table(
  dimension = c("Biological scope","State continuum","Principal regulators","Inferred function","Plasticity context",
    "Malignant-state specificity","Proliferation dependency","Perturbation evidence","External validation","Experimental validation"),
  published_foxm1_cebpb = c(
    "HCC developmental heterogeneity and therapeutic resistance", "Developmental state plasticity", "FOXM1 and CEBPB toggle",
    "Mutual suppression regulates cell-state plasticity", "Developmental heterogeneity and resistance",
    "HCC experimental models", "FOXM1 is proliferation/plasticity-associated", "Experimental genetic perturbation",
    "Clinical and preclinical HCC context", "In vitro, in vivo and preclinical validation"
  ),
  current_three_axis = c(
    "Computational HCC hepatocyte-state acquisition architecture", "Identity loss → stress transition → SOX4-associated state",
    "HNF4A/PPARA; AP-1/CEBPB/EGR1; SOX4", "Computationally inferred coupling among three frozen programmes",
    "State-specific extension across identity, injury and malignant-like programmes", "SOX4-associated malignant-state stabilization hypothesis",
    "Proliferation measured as a separate response; sparse frozen gene coverage", "CellOracle KO plus unsigned scTenifoldKnk magnitude",
    "Frozen module8 cohorts provide association-level context", "No new wet-lab validation in Figure 6"
  ),
  relationship = c("complementary","partially overlapping","partially overlapping","partially overlapping","complementary",
    "distinct analytical scope","unresolved","distinct evidence level","complementary","distinct evidence level"),
  evidence_source = c(rep(paste0("PubMed ", paper$pubmed, "; project Figure 2-6 frozen outputs"), 10))
)
figure6_fwrite(tbl, file.path(FIGURE6_METADATA_DIR, "figure6h_comparison_table.tsv"), compress = FALSE)
figure6_write_json(list(
  panel = "Figure 6H", paper = paper, relationship_terms = c("complementary","partially overlapping","distinct analytical scope","state-specific extension","computationally inferred architecture"),
  foxm1_project_perturbation_available = FALSE, cebpb_project_ko_available = TRUE,
  interpretation = "The models are complementary and partially overlapping, with distinct analytical and validation scope.",
  prohibited_claims = c("our model supersedes the published model", "FOXM1/CEBPB is incorrect", "SOX4 is the only master regulator")
), file.path(FIGURE6_METADATA_DIR, "figure6h_comparison_report.json"))

