root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
script <- file.path(root, "scripts", "figure8_v2_07_prism.R")
if (!file.exists(script)) stop("Expected RED failure: PRISM module is not implemented")
Sys.setenv(FIGURE8_V2_TEST_MODE = "1")
source(script, local = FALSE)

model <- data.table::data.table(
  CellLineName = c("HCC515", "HA1E", "Hep G2", "SNU-398", "KMCH-1", "OTHER"),
  OncotreeCode = c("LUAD", NA, "LIHB", "HCC", "HCCIHCH", "BRCA"),
  OncotreePrimaryDisease = c("Non-Small Cell Lung Cancer", "Non-Cancerous", "Hepatoblastoma", "Hepatocellular Carcinoma", "Hepatocellular Carcinoma plus Intrahepatic Cholangiocarcinoma", "Breast Cancer"),
  OncotreeLineage = c("Lung", "Kidney", "Liver", "Liver", "Liver", "Breast")
)
labels <- mapply(figure8_v2_context_label, model$CellLineName, model$OncotreeCode, model$OncotreePrimaryDisease, model$OncotreeLineage)
stopifnot(labels[[1]] == "lung_adenocarcinoma_non_liver")
stopifnot(labels[[2]] == "kidney_derived_non_liver")
stopifnot(labels[[3]] == "liver_derived_hepatoblastoma_like_caveat")
stopifnot(labels[[4]] == "adult_hepatocellular_carcinoma")
stopifnot(labels[[5]] == "mixed_hcc_intrahepatic_cholangiocarcinoma")
stopifnot(labels[[6]] == "other_cancer_non_liver")

stopifnot(figure8_v2_choose_secondary_screen(c("HTS002", "MTS010")) == "MTS010")
stopifnot(figure8_v2_choose_secondary_screen(c("HTS002")) == "HTS002")

secondary_toy <- figure8_v2_add_secondary_sensitivity(data.table::data.table(auc = c(0.2, 1.1)))
stopifnot(identical(secondary_toy$secondary_sensitivity, c(-0.2, -1.1)))
stopifnot(figure8_v2_combine_prism_classes("hcc_liver_enriched", NA_character_) == "hcc_liver_enriched_single_release")
stopifnot(figure8_v2_combine_prism_classes("hcc_liver_enriched", "broad_cytotoxicity") == "discordant_phenotype")
stopifnot(figure8_v2_combine_prism_classes("no_enriched_support", "broad_cytotoxicity") == "broad_cytotoxicity")

cat("figure8_v2 PRISM logic tests passed\n")
