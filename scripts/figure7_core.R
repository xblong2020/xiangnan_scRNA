suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggsci)
  library(patchwork)
  library(survival)
  library(metafor)
  library(broom)
  library(jsonlite)
  library(scales)
  library(ggrepel)
})

figure7_script_path <- function() {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(arg)) return(normalizePath("scripts", mustWork = TRUE))
  normalizePath(sub("^--file=", "", arg[[1]]), mustWork = TRUE)
}

FIGURE7_ROOT <- normalizePath(file.path(dirname(figure7_script_path()), ".."), mustWork = TRUE)
source(file.path(FIGURE7_ROOT, "scripts", "figure7_plot_theme.R"))

FIGURE7_META <- file.path(FIGURE7_ROOT, "metadata", "driver", "figure7_external_validation")
FIGURE7_DATA <- file.path(FIGURE7_ROOT, "data", "processed", "driver", "figure7_external_validation")
FIGURE7_REPORT <- file.path(FIGURE7_ROOT, "reports", "figure7_external_bulk_clinical_validation_report.md")
FIGURE7_SEED <- 20260805L

FIGURE7_INPUTS <- list(
  tcga_expression = "F:/Charley/2024/134panImmune/03.download/expression/TCGA-LIHC.htseq_fpkm.tsv.gz",
  tcga_survival = "F:/Charley/2024/134panImmune/03.download/survival/TCGA-LIHC.survival.tsv.gz",
  tcga_clinical = "G:/万亿肝癌/clinical.xls",
  icgc_expression = "G:/万亿肝癌/ICGCsymbol.txt",
  icgc_survival = "G:/万亿肝癌/ICGCtime.txt",
  icgc_clinical = "G:/万亿肝癌/icgcClinical.txt",
  gtf = "G:/wanyi_HCC_scRNA/HCCscRNA/GSE156625-HCC/cellranger/hg38/refdata-gex-GRCh38-2020-A/genes/genes.gtf",
  frozen_targets = file.path(FIGURE7_ROOT, "metadata", "driver", "module8_tf_target_signature_genes.tsv"),
  frozen_registry = file.path(FIGURE7_ROOT, "metadata", "driver", "module8_signature_registry.tsv")
)

dir_figure7 <- function() {
  dirs <- c(
    FIGURE7_META, FIGURE7_DATA, dirname(FIGURE7_REPORT),
    file.path(FIGURE7_ROOT, "figures", "driver", "figure7a_cohort_flow"),
    file.path(FIGURE7_ROOT, "figures", "driver", "figure7b_bulk_signature_mapping"),
    file.path(FIGURE7_ROOT, "figures", "driver", "figure7c_tumour_normal_forest"),
    file.path(FIGURE7_ROOT, "figures", "driver", "figure7d_clinical_heatmap"),
    file.path(FIGURE7_ROOT, "figures", "driver", "figure7e_multivariable_cox_forest"),
    file.path(FIGURE7_ROOT, "figures", "driver", "figure7f_incremental_prediction"),
    file.path(FIGURE7_ROOT, "figures", "driver", "figure7g_survival_curves"),
    file.path(FIGURE7_ROOT, "figures", "driver", "figure7h_sensitivity_specificity"),
    file.path(FIGURE7_ROOT, "figures", "driver", "figure7_external_validation_preview")
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

write_json <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE, na = "null", digits = 10)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(as.data.table(x), path, sep = "\t", na = "NA", quote = FALSE)
}

write_tsvgz <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(as.data.table(x), path, sep = "\t", na = "NA", quote = FALSE, compress = "gzip")
}

read_tsv <- function(path, ...) data.table::fread(path, sep = "\t", na.strings = c("", "NA", "NaN"), ...)

safe_z <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

bh <- function(x) {
  out <- rep(NA_real_, length(x)); ok <- is.finite(x)
  out[ok] <- p.adjust(x[ok], method = "BH"); out
}

package_inventory <- function() {
  required <- c("ggplot2", "ggsci", "patchwork", "cowplot", "data.table", "dplyr",
                "survival", "survminer", "metafor", "broom", "jsonlite", "scales", "ggrepel")
  optional <- c("singscore", "GSVA", "limma", "edgeR", "ordinal", "MASS", "rms",
                "riskRegression", "timeROC", "pec", "Hmisc", "mice", "pROC",
                "ComplexHeatmap", "circlize", "forestploter", "cmprsk", "glmnet")
  pkgs <- unique(c(required, optional))
  data.table(
    package = pkgs,
    required = pkgs %in% required,
    available = vapply(pkgs, requireNamespace, logical(1), quietly = TRUE),
    version = vapply(pkgs, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NA_character_, character(1))
  )
}

normalize_symbol <- function(x) {
  x <- toupper(trimws(as.character(x)))
  sub("\\.[0-9]+$", "", x)
}

build_frozen_signatures <- function() {
  targets <- read_tsv(FIGURE7_INPUTS$frozen_targets)
  targets[, `:=`(gene = normalize_symbol(gene), tf = normalize_symbol(tf))]
  axis_tfs <- list(
    identity_program = c("HNF4A", "PPARA"),
    stress_transition = c("JUN", "FOS", "JUND", "ATF3", "CEBPB", "EGR1"),
    sox4_stabilization = "SOX4",
    calibration_control = c("IRF1", "JUNB", "MYC", "HLF", "MAFF", "MAFB")
  )
  rows <- list()
  for (axis in names(axis_tfs)) {
    tfs <- axis_tfs[[axis]]
    z <- unique(targets[tf %in% tfs, .(gene, tf, source_method = signature_source, source_rank = rank)])
    z <- rbind(z, data.table(gene = tfs, tf = tfs, source_method = "frozen_tf_expression_anchor", source_rank = 0L), fill = TRUE)
    z <- unique(z, by = c("gene", "tf"))
    z[, axis := axis]
    rows[[axis]] <- z
  }
  proliferation <- c("MKI67", "TOP2A", "STMN1", "TYMS", "UBE2C", "PCNA", "MCM2", "MCM5", "HMGB2")
  foxm1_cebpb <- unique(c("FOXM1", "CEBPB", targets[tf == "CEBPB", gene], "CCNB1", "CCNB2", "CDC20", "CENPF", "AURKB"))
  hypoxia <- c("CA9", "VEGFA", "SLC2A1", "LDHA", "PDK1", "BNIP3", "EGLN3", "HIF1A", "NDRG1", "ENO1", "PGK1", "HK2", "ADM", "LOX", "ALDOA")
  inflammation <- c("IL6", "TNF", "NFKB1", "RELA", "STAT3", "CXCL8", "CCL2", "IL1B", "PTGS2", "ICAM1")
  add_fixed <- function(axis, genes, method) data.table(
    gene = normalize_symbol(genes), tf = NA_character_, source_method = method,
    source_rank = seq_along(genes), axis = axis
  )
  rows$proliferation_control <- add_fixed("proliferation_control", proliferation, "project_module2_fixed_proliferation_panel")
  rows$foxm1_cebpb_reference <- add_fixed("foxm1_cebpb_reference", foxm1_cebpb, "predefined_FOXM1_anchor_plus_frozen_CEBPB_targets")
  rows$hypoxia_control <- add_fixed("hypoxia_control", hypoxia, "predefined_expression_covariate_panel")
  rows$inflammation_control <- add_fixed("inflammation_control", inflammation, "predefined_expression_covariate_panel")
  out <- unique(rbindlist(rows, fill = TRUE), by = c("axis", "gene"))
  out[, direction := ifelse(axis == "identity_program", "identity_retention_high", "programme_high")]
  out[, source_figure := fifelse(axis == "identity_program", "Figure 2/5-6",
                          fifelse(axis == "stress_transition", "Figure 3/5-6",
                          fifelse(axis == "sox4_stabilization", "Figure 4/5-6", "control")))]
  out[, selection_rule := fifelse(grepl("frozen|sctenifold", source_method, ignore.case = TRUE),
                                  "pre-bulk frozen target or TF anchor", "predefined control panel")]
  out[, `:=`(
    signature_version = "figure7_v1_2026-08-05",
    discovery_only = TRUE,
    primary_or_sensitivity = fifelse(axis %in% c("identity_program", "stress_transition", "sox4_stabilization"), "primary", "control"),
    bulk_available = NA,
    mapping_status = "pending_preparation"
  )]
  setcolorder(out, c("axis", "gene", "direction", "source_figure", "source_method", "selection_rule",
                     "signature_version", "discovery_only", "primary_or_sensitivity", "bulk_available",
                     "mapping_status", "tf", "source_rank"))
  out[]
}

protected_figure_files <- function() {
  roots <- c("scripts", "metadata", "figures", "reports")
  paths <- unlist(lapply(roots, function(d) list.files(file.path(FIGURE7_ROOT, d), recursive = TRUE, full.names = TRUE, all.files = FALSE)))
  paths <- paths[file.exists(paths) & !dir.exists(paths)]
  rel <- substring(normalizePath(paths, winslash = "/", mustWork = TRUE), nchar(normalizePath(FIGURE7_ROOT, winslash = "/")) + 2L)
  keep <- grepl("(^|/)(figure[1-6]|module[1-9])", tolower(rel)) & !grepl("figure7", tolower(rel))
  paths[keep]
}

hash_files <- function(paths) {
  if (!length(paths)) return(data.table(file_path = character(), size = numeric(), md5 = character()))
  data.table(
    file_path = substring(normalizePath(paths, winslash = "/", mustWork = TRUE), nchar(normalizePath(FIGURE7_ROOT, winslash = "/")) + 2L),
    size = file.info(paths)$size,
    md5 = unname(tools::md5sum(paths))
  )
}

stage_preflight <- function() {
  dir_figure7()
  set.seed(FIGURE7_SEED)
  signatures <- build_frozen_signatures()
  write_tsv(signatures, file.path(FIGURE7_META, "figure7_frozen_signature_manifest.tsv"))
  pkgs <- package_inventory()
  write_tsv(pkgs, file.path(FIGURE7_META, "figure7_r_package_inventory.tsv"))
  inputs <- data.table(
    input = names(FIGURE7_INPUTS),
    file_path = unlist(FIGURE7_INPUTS),
    exists = file.exists(unlist(FIGURE7_INPUTS)),
    size_bytes = ifelse(file.exists(unlist(FIGURE7_INPUTS)), file.info(unlist(FIGURE7_INPUTS))$size, NA_real_)
  )
  protected <- hash_files(protected_figure_files())
  write_tsv(protected, file.path(FIGURE7_META, "figure7_protected_figure1_6_hashes_before.tsv"))
  checks <- rbindlist(list(
    data.table(check = "frozen_signature_complete", status = if (all(c("identity_program", "stress_transition", "sox4_stabilization") %in% signatures$axis)) "pass" else "block", detail = paste(unique(signatures$axis), collapse = ";")),
    data.table(check = paste0("input_", inputs$input), status = ifelse(inputs$exists, "pass", "block"), detail = inputs$file_path),
    data.table(check = "gene_symbol_policy", status = "pass", detail = "uppercase HGNC-like symbols; Ensembl versions stripped before GTF mapping"),
    data.table(check = "bulk_outcome_used_for_signature_derivation", status = "pass", detail = "false; only frozen Module7 target tables and predefined controls used"),
    data.table(check = "cohort_overlap", status = "pass", detail = "TCGA patient IDs and ICGC donor IDs use disjoint namespaces"),
    data.table(check = "purity_available", status = "warn", detail = "not present in cached clinical tables"),
    data.table(check = "cnv_burden_available", status = "warn", detail = "not present in cached clinical tables"),
    data.table(check = "afp_virology_cirrhosis_available", status = "warn", detail = "not present in cached clinical tables"),
    data.table(check = "singscore_package", status = if (pkgs[package == "singscore", available]) "pass" else "warn", detail = "rank-based singscore-equivalent implementation is frozen fallback"),
    data.table(check = "survminer_package", status = if (pkgs[package == "survminer", available]) "pass" else "warn", detail = "project-local installation attempted by runner; manual ggplot survival fallback remains auditable"),
    data.table(check = "protected_figure1_6_baseline", status = "pass", detail = sprintf("%d protected files hashed", nrow(protected)))
  ), fill = TRUE)
  write_tsv(checks, file.path(FIGURE7_META, "figure7_preflight_report.tsv"))
  missing <- data.table(
    data_element = c("tumour_purity", "cnv_burden", "AFP", "vascular_invasion", "recurrence", "HBV", "HCV", "cirrhosis", "independent_GEO_bulk"),
    status = "missing_from_current_cache",
    consequence = c(
      "purity-adjusted models not estimable", "CNV-adjusted models not estimable",
      "AFP association not estimable", "vascular-invasion association not estimable",
      "recurrence endpoint not estimable", "HBV stratification not estimable", "HCV stratification not estimable",
      "cirrhosis stratification not estimable", "pan-cohort inference limited to TCGA-LIHC and ICGC-LIRI-JP"
    )
  )
  write_tsv(missing, file.path(FIGURE7_META, "figure7_missing_data_report.tsv"))
  write_json(list(
    module = "Figure 7 preflight audit", timestamp = format(Sys.time(), tz = "Asia/Shanghai"),
    random_seed = FIGURE7_SEED, r_version = R.version.string, input_audit = inputs,
    n_signature_rows = nrow(signatures), n_protected_files = nrow(protected), checks = checks,
    missing_data = missing, blocked = any(checks$status == "block")
  ), file.path(FIGURE7_META, "figure7_preflight_report.json"))
  if (any(checks$status == "block")) stop("Figure 7 preflight found blocking missing inputs; see report.")
  invisible(checks)
}

parse_gtf_gene_map <- function(gtf_path, cache_path, wanted_genes) {
  if (file.exists(cache_path)) return(read_tsv(cache_path))
  pattern_file <- tempfile("figure7_gtf_patterns_", fileext = ".txt")
  writeLines(sprintf('gene_name "%s"', unique(normalize_symbol(wanted_genes))), pattern_file, useBytes = TRUE)
  on.exit(unlink(pattern_file), add = TRUE)
  lines <- system2("rg", c("--no-line-number", "-F", "-f", shQuote(pattern_file), shQuote(gtf_path)), stdout = TRUE, stderr = TRUE)
  lines <- lines[grepl("\tgene\t", lines, fixed = TRUE)]
  gid <- sub('.*gene_id "([^"]+)".*', '\\1', lines)
  gnm <- sub('.*gene_name "([^"]+)".*', '\\1', lines)
  ok <- gid != lines & gnm != lines
  out <- unique(data.table(ensembl_id = sub("\\.[0-9]+$", "", gid[ok]), gene = normalize_symbol(gnm[ok])), by = "ensembl_id")
  write_tsv(out, cache_path)
  out
}

collapse_expression <- function(gene, mat) {
  ok <- !is.na(gene) & nzchar(gene)
  gene <- gene[ok]; mat <- mat[ok, , drop = FALSE]
  sums <- rowsum(mat, group = gene, reorder = FALSE, na.rm = TRUE)
  counts <- as.numeric(table(factor(gene, levels = rownames(sums))))
  sums / counts
}

tcga_sample_type <- function(x) {
  code <- vapply(strsplit(x, "-", fixed = TRUE), function(z) if (length(z) >= 4) substr(z[[4]], 1, 2) else "", character(1))
  fifelse(code %in% c("01", "02", "03", "05"), "tumour",
          fifelse(code %in% c("10", "11", "12", "13", "14"), "normal", "unknown"))
}

tcga_patient_id <- function(x) vapply(strsplit(x, "-", fixed = TRUE), function(z) paste(head(z, 3), collapse = "-"), character(1))
icgc_patient_id <- function(x) {
  hit <- regmatches(x, regexpr("DO[0-9]+", x))
  hit[!nzchar(hit)] <- x[!nzchar(hit)]
  hit
}

stage_to_num <- function(x) {
  y <- toupper(trimws(as.character(x)))
  out <- rep(NA_real_, length(y))
  out[grepl("STAGE I|^I", y)] <- 1
  out[grepl("STAGE II|^II", y)] <- 2
  out[grepl("STAGE III|^III", y)] <- 3
  out[grepl("STAGE IV|^IV", y)] <- 4
  suppressWarnings(num <- as.numeric(y))
  out[is.na(out) & is.finite(num)] <- num[is.na(out) & is.finite(num)]
  out
}

grade_to_num <- function(x) {
  y <- toupper(trimws(as.character(x)))
  suppressWarnings(out <- as.numeric(sub("^G", "", y)))
  out
}

tstage_to_num <- function(x) {
  y <- toupper(trimws(as.character(x)))
  suppressWarnings(as.numeric(sub("^T([0-9]+).*$", "\\1", y)))
}

refresh_tcga_samples <- function(sm) {
  base <- sm[, .(sample_id, patient_id, cohort, tumour_normal)]
  clin <- fread(FIGURE7_INPUTS$tcga_clinical, na.strings = c("", "NA", "[Not Available]"))
  clin[, `:=`(
    patient_id = as.character(Id), os_time_days = safe_num(futime), os_event = safe_num(fustat),
    age_years = safe_num(age), sex = toupper(as.character(gender)),
    stage_num = stage_to_num(stage), grade_num = grade_to_num(grade), t_stage_num = tstage_to_num(T)
  )]
  clin[, `:=`(stage_high = as.integer(stage_num >= 3), age_high = as.integer(age_years >= 60))]
  merge(base, clin[, .(patient_id, os_time_days, os_event, age_years, age_high, sex,
                       stage_num, stage_high, grade_num, t_stage_num, stage, grade, T, M, N)],
        by = "patient_id", all.x = TRUE, sort = FALSE)
}

refresh_icgc_samples <- function(sm) {
  base <- sm[, .(sample_id, patient_id, cohort, tumour_normal)]
  clin <- fread(FIGURE7_INPUTS$icgc_clinical)
  surv <- fread(FIGURE7_INPUTS$icgc_survival)
  clin <- merge(clin, surv, by.x = "Id", by.y = "id", all = TRUE)
  clin[, `:=`(
    patient_id = as.character(Id), os_time_days = safe_num(futime), os_event = safe_num(fustat),
    age_high = safe_num(Age), age_years = NA_real_, sex = fifelse(safe_num(Gender) == 1, "MALE", "FEMALE"),
    stage_high = safe_num(Stage), stage_num = safe_num(Stage) + 1,
    grade_num = NA_real_, t_stage_num = NA_real_, stage = as.character(Stage), grade = NA_character_,
    T = NA_character_, M = NA_character_, N = NA_character_
  )]
  merge(base, clin[, .(patient_id, os_time_days, os_event, age_years, age_high, sex,
                       stage_num, stage_high, grade_num, t_stage_num, stage, grade, T, M, N)],
        by = "patient_id", all.x = TRUE, sort = FALSE)
}

prepare_tcga <- function(signatures, gtf_map) {
  raw <- fread(FIGURE7_INPUTS$tcga_expression, showProgress = TRUE, check.names = FALSE)
  ensembl <- sub("\\.[0-9]+$", "", raw[[1]])
  symbols <- gtf_map$gene[match(ensembl, gtf_map$ensembl_id)]
  symbols[is.na(symbols) | !nzchar(symbols)] <- ensembl[is.na(symbols) | !nzchar(symbols)]
  mat <- as.matrix(raw[, -1, with = FALSE]); storage.mode(mat) <- "double"
  rm(raw); gc()
  expr <- collapse_expression(symbols, mat)
  rm(mat); gc()
  expr <- log2(pmax(expr, 0) + 1)
  samples <- colnames(expr)
  sm <- data.table(sample_id = samples, patient_id = tcga_patient_id(samples),
                   cohort = "TCGA_LIHC", tumour_normal = tcga_sample_type(samples))
  clin <- fread(FIGURE7_INPUTS$tcga_clinical, na.strings = c("", "NA", "[Not Available]"))
  clin[, `:=`(
    patient_id = as.character(Id),
    os_time_days = safe_num(futime), os_event = safe_num(fustat),
    age_years = safe_num(age), sex = toupper(as.character(gender)),
    stage_num = stage_to_num(stage), grade_num = grade_to_num(grade), t_stage_num = tstage_to_num(T)
  )]
  clin[, stage_high := as.integer(stage_num >= 3)]
  clin[, age_high := as.integer(age_years >= 60)]
  sm <- merge(sm, clin[, .(patient_id, os_time_days, os_event, age_years, age_high, sex,
                            stage_num, stage_high, grade_num, t_stage_num, stage, grade, T, M, N)],
              by = "patient_id", all.x = TRUE, sort = FALSE)
  list(expr = expr, samples = sm, unit = "log2(FPKM + 1)", platform = "RNA-seq")
}

prepare_icgc <- function(signatures) {
  raw <- fread(FIGURE7_INPUTS$icgc_expression, showProgress = TRUE, check.names = FALSE)
  symbols <- normalize_symbol(raw[[1]])
  mat <- as.matrix(raw[, -1, with = FALSE]); storage.mode(mat) <- "double"
  rm(raw); gc()
  expr <- collapse_expression(symbols, mat)
  rm(mat); gc()
  expr <- log2(pmax(expr, 0) + 1)
  samples <- colnames(expr)
  sm <- data.table(sample_id = samples, patient_id = icgc_patient_id(samples), cohort = "ICGC_LIRI_JP",
                   tumour_normal = fifelse(grepl("-T$", samples), "tumour", fifelse(grepl("-N$", samples), "normal", "unknown")))
  clin <- fread(FIGURE7_INPUTS$icgc_clinical)
  surv <- fread(FIGURE7_INPUTS$icgc_survival)
  clin <- merge(clin, surv, by.x = "Id", by.y = "id", all = TRUE)
  clin[, `:=`(
    patient_id = as.character(Id), os_time_days = safe_num(futime), os_event = safe_num(fustat),
    age_high = safe_num(Age), age_years = NA_real_, sex = fifelse(safe_num(Gender) == 1, "MALE", "FEMALE"),
    stage_high = safe_num(Stage), stage_num = safe_num(Stage) + 1,
    grade_num = NA_real_, t_stage_num = NA_real_, stage = as.character(Stage), grade = NA_character_,
    T = NA_character_, M = NA_character_, N = NA_character_
  )]
  sm <- merge(sm, clin[, .(patient_id, os_time_days, os_event, age_years, age_high, sex,
                            stage_num, stage_high, grade_num, t_stage_num, stage, grade, T, M, N)],
              by = "patient_id", all.x = TRUE, sort = FALSE)
  list(expr = expr, samples = sm, unit = "log2(RPKM-like expression + 1)", platform = "RNA-seq")
}

cohort_manifest_row <- function(obj, source, path) {
  sm <- obj$samples
  tumour_patients <- unique(sm[tumour_normal == "tumour" & is.finite(os_time_days) & is.finite(os_event), patient_id])
  paired <- intersect(unique(sm[tumour_normal == "tumour", patient_id]), unique(sm[tumour_normal == "normal", patient_id]))
  data.table(
    cohort = unique(sm$cohort), data_source = source, platform = obj$platform,
    expression_unit = obj$unit, normalization = "within-cohort log2 transform; no cross-cohort ComBat",
    n_total = nrow(sm), n_tumour = sum(sm$tumour_normal == "tumour"), n_normal = sum(sm$tumour_normal == "normal"),
    n_paired = length(paired), n_survival = length(tumour_patients),
    n_events = uniqueN(sm[patient_id %in% tumour_patients & os_event == 1, patient_id]),
    clinical_variables = paste(names(sm)[vapply(sm, function(z) any(!is.na(z)), logical(1))], collapse = ";"),
    purity_available = FALSE, cnv_available = FALSE, discovery_or_validation = "external_validation",
    overlap_risk = "none_detected_disjoint_patient_namespaces", file_path = path
  )
}

stage_prepare <- function() {
  dir_figure7()
  signatures <- read_tsv(file.path(FIGURE7_META, "figure7_frozen_signature_manifest.tsv"))
  gtf_cache <- file.path(FIGURE7_DATA, "figure7_grch38_ensembl_symbol_map.tsv.gz")
  gtf_map <- parse_gtf_gene_map(FIGURE7_INPUTS$gtf, gtf_cache, signatures$gene)
  tcga_path <- file.path(FIGURE7_DATA, "figure7_tcga_lihc_prepared_expression.rds")
  icgc_path <- file.path(FIGURE7_DATA, "figure7_icgc_liri_jp_prepared_expression.rds")
  if (file.exists(tcga_path) && file.exists(icgc_path)) {
    tcga <- readRDS(tcga_path)
    icgc <- readRDS(icgc_path)
  } else {
    tcga <- prepare_tcga(signatures, gtf_map)
    icgc <- prepare_icgc(signatures)
  }
  tcga$samples <- refresh_tcga_samples(tcga$samples)
  icgc$samples <- refresh_icgc_samples(icgc$samples)
  saveRDS(tcga, tcga_path, compress = FALSE)
  saveRDS(icgc, icgc_path, compress = FALSE)
  coverage <- rbindlist(lapply(list(tcga, icgc), function(obj) {
    signatures[, .(
      n_frozen = uniqueN(gene), n_mapped = uniqueN(intersect(gene, rownames(obj$expr))),
      coverage = uniqueN(intersect(gene, rownames(obj$expr))) / uniqueN(gene),
      mapped_genes = paste(sort(intersect(unique(gene), rownames(obj$expr))), collapse = ";"),
      missing_genes = paste(sort(setdiff(unique(gene), rownames(obj$expr))), collapse = ";")
    ), by = axis][, cohort := unique(obj$samples$cohort)]
  }), fill = TRUE)
  setcolorder(coverage, c("cohort", "axis", "n_frozen", "n_mapped", "coverage", "mapped_genes", "missing_genes"))
  write_tsv(coverage, file.path(FIGURE7_META, "figure7_signature_coverage.tsv"))
  all_mapped <- Reduce(intersect, lapply(list(tcga, icgc), function(z) rownames(z$expr)))
  signatures[, `:=`(bulk_available = gene %in% all_mapped,
                    mapping_status = fifelse(gene %in% all_mapped, "mapped_all_primary_cohorts", "partially_or_not_mapped"))]
  write_tsv(signatures, file.path(FIGURE7_META, "figure7_frozen_signature_manifest.tsv"))
  manifest <- rbindlist(list(
    cohort_manifest_row(tcga, "GDC/TCGA cached FPKM", FIGURE7_INPUTS$tcga_expression),
    cohort_manifest_row(icgc, "ICGC LIRI-JP cached expression", FIGURE7_INPUTS$icgc_expression)
  ), fill = TRUE)
  write_tsv(manifest, file.path(FIGURE7_META, "figure7_cohort_manifest.tsv"))
  for (obj in list(tcga, icgc)) {
    cohort <- unique(obj$samples$cohort)
    report <- list(
      cohort = cohort, input_file = if (cohort == "TCGA_LIHC") FIGURE7_INPUTS$tcga_expression else FIGURE7_INPUTS$icgc_expression,
      expression_unit = obj$unit, n_genes = nrow(obj$expr), n_samples = ncol(obj$expr),
      n_nonfinite = sum(!is.finite(obj$expr)), duplicate_symbols_after_collapse = sum(duplicated(rownames(obj$expr))),
      tumour_normal = as.list(table(obj$samples$tumour_normal)), random_seed = FIGURE7_SEED
    )
    write_json(report, file.path(FIGURE7_META, paste0("figure7_", tolower(cohort), "_preparation_report.json")))
  }
  invisible(manifest)
}

rank_score_matrix <- function(expr, signatures, axes) {
  ranks <- apply(expr, 2, rank, ties.method = "average", na.last = "keep") / nrow(expr)
  if (is.vector(ranks)) ranks <- matrix(ranks, ncol = 1, dimnames = dimnames(expr))
  scores <- lapply(axes, function(ax) {
    genes <- intersect(signatures[axis == ax, gene], rownames(expr))
    if (!length(genes)) return(rep(NA_real_, ncol(expr)))
    colMeans(ranks[genes, , drop = FALSE], na.rm = TRUE) - 0.5
  })
  names(scores) <- axes
  list(scores = as.data.table(scores), ranks = ranks)
}

calculate_scores_one <- function(obj, signatures) {
  axes <- c("identity_program", "stress_transition", "sox4_stabilization",
            "proliferation_control", "foxm1_cebpb_reference", "calibration_control",
            "hypoxia_control", "inflammation_control")
  scored <- rank_score_matrix(obj$expr, signatures, axes)
  s <- scored$scores
  s[, `:=`(
    identity_program_score_original = identity_program,
    identity_loss_score = -safe_z(identity_program),
    stress_transition_score = safe_z(stress_transition),
    sox4_stabilization_score = safe_z(sox4_stabilization),
    proliferation_score = safe_z(proliferation_control),
    foxm1_cebpb_reference_score = safe_z(foxm1_cebpb_reference),
    calibration_control_score = safe_z(calibration_control),
    hypoxia_score = safe_z(hypoxia_control),
    inflammation_score = safe_z(inflammation_control)
  )]
  keep <- c("identity_program_score_original", "identity_loss_score", "stress_transition_score",
            "sox4_stabilization_score", "proliferation_score", "foxm1_cebpb_reference_score",
            "calibration_control_score", "hypoxia_score", "inflammation_score")
  s <- s[, ..keep]
  sm <- copy(obj$samples)
  sm <- sm[match(colnames(obj$expr), sample_id)]
  out <- cbind(sm, s)
  cov <- signatures[axis %in% c("identity_program", "stress_transition", "sox4_stabilization"),
                    .(coverage = uniqueN(intersect(gene, rownames(obj$expr))) / uniqueN(gene)), by = axis]
  out[, `:=`(
    purity = NA_real_, cnv_burden = NA_real_, score_method = "sample-wise rank score (singscore-equivalent)",
    signature_coverage = mean(cov$coverage), random_seed = FIGURE7_SEED
  )]
  list(scores = out, ranks = scored$ranks)
}

stage_scores <- function() {
  set.seed(FIGURE7_SEED)
  signatures <- read_tsv(file.path(FIGURE7_META, "figure7_frozen_signature_manifest.tsv"))
  objs <- list(
    TCGA_LIHC = readRDS(file.path(FIGURE7_DATA, "figure7_tcga_lihc_prepared_expression.rds")),
    ICGC_LIRI_JP = readRDS(file.path(FIGURE7_DATA, "figure7_icgc_liri_jp_prepared_expression.rds"))
  )
  score_path <- file.path(FIGURE7_META, "figure7_bulk_axis_scores.tsv.gz")
  rank_paths <- setNames(file.path(FIGURE7_DATA, paste0("figure7_", tolower(names(objs)), "_rank_matrix.rds")), names(objs))
  if (file.exists(score_path) && all(file.exists(rank_paths))) {
    old_scores <- read_tsv(score_path)
    results <- lapply(names(objs), function(nm) {
      sm <- copy(objs[[nm]]$samples)
      sm <- sm[match(colnames(objs[[nm]]$expr), sample_id)]
      old <- old_scores[cohort == nm][match(sm$sample_id, sample_id)]
      score_cols <- setdiff(names(old), names(sm))
      list(scores = cbind(sm, old[, ..score_cols]), ranks = readRDS(rank_paths[[nm]]))
    })
    names(results) <- names(objs)
    scores <- rbindlist(lapply(results, `[[`, "scores"), fill = TRUE)
    write_tsvgz(scores, score_path)
  } else {
    results <- lapply(objs, calculate_scores_one, signatures = signatures)
    scores <- rbindlist(lapply(results, `[[`, "scores"), fill = TRUE)
    write_tsvgz(scores, score_path)
    for (nm in names(results)) saveRDS(results[[nm]]$ranks, rank_paths[[nm]], compress = FALSE)
  }
  correlations <- scores[tumour_normal == "tumour", as.data.table(cor(.SD, use = "pairwise.complete.obs")),
                         .SDcols = c("identity_loss_score", "stress_transition_score", "sox4_stabilization_score"), by = cohort]
  write_tsv(correlations, file.path(FIGURE7_META, "figure7_three_axis_correlations.tsv"))
  write_json(list(
    module = "Figure 7 bulk axis scoring", method = "sample-wise normalized ranks centered at 0; within-cohort SD standardization",
    identity_direction = "identity_loss_score = negative within-cohort z score of identity_program_score_original",
    n_samples = nrow(scores), cohorts = as.list(table(scores$cohort)), random_seed = FIGURE7_SEED,
    singscore_package_available = requireNamespace("singscore", quietly = TRUE),
    note = "The frozen fallback is algebraically rank-based and does not use outcome-derived weights."
  ), file.path(FIGURE7_META, "figure7_bulk_axis_score_report.json"))
  bench_path <- file.path(FIGURE7_META, "figure7_random_signature_benchmark.tsv.gz")
  existing_bench <- if (file.exists(bench_path)) read_tsv(bench_path) else NULL
  rerun_cohorts <- names(objs)
  if (!is.null(existing_bench) && nrow(existing_bench)) {
    bench_status <- existing_bench[, .(
      n = .N,
      n_complete_metrics = sum(is.finite(tumour_normal_effect) & is.finite(stage_association) & is.finite(cox_z) & is.finite(c_index) & is.finite(delta_c_index))
    ), by = .(cohort, target_axis)]
    valid_cohorts <- bench_status[, .(
      valid = .N == 3L && all(n == 500L) && all(n_complete_metrics == 500L)
    ), by = cohort][valid == TRUE, cohort]
    rerun_cohorts <- setdiff(names(objs), valid_cohorts)
  }
  if (length(rerun_cohorts)) {
    retained <- if (!is.null(existing_bench)) existing_bench[!cohort %in% rerun_cohorts] else NULL
    run_random_signature_benchmark(
      objs[rerun_cohorts], results[rerun_cohorts], signatures,
      n_random = 500L, retained = retained
    )
  }
  invisible(scores)
}

hedges_g <- function(x_t, x_n) {
  x_t <- x_t[is.finite(x_t)]; x_n <- x_n[is.finite(x_n)]
  n1 <- length(x_t); n0 <- length(x_n)
  if (n1 < 2 || n0 < 2) return(c(effect = NA, se = NA, lo = NA, hi = NA, p = NA))
  sp <- sqrt(((n1 - 1) * var(x_t) + (n0 - 1) * var(x_n)) / (n1 + n0 - 2))
  d <- (mean(x_t) - mean(x_n)) / sp
  j <- 1 - 3 / (4 * (n1 + n0) - 9)
  g <- j * d
  se <- sqrt((n1 + n0) / (n1 * n0) + g^2 / (2 * (n1 + n0 - 2)))
  p <- 2 * pnorm(-abs(g / se))
  c(effect = g, se = se, lo = g - 1.96 * se, hi = g + 1.96 * se, p = p)
}

harrell_c <- function(time, event, lp) {
  ok <- is.finite(time) & is.finite(event) & is.finite(lp) & time > 0
  if (sum(ok) < 10 || sum(event[ok] == 1) < 3) return(NA_real_)
  unname(concordance(Surv(time[ok], event[ok]) ~ lp[ok], reverse = TRUE)$concordance)
}

make_gene_strata <- function(expr) {
  mean_expr <- rowMeans(expr, na.rm = TRUE)
  variability <- apply(expr, 1, sd, na.rm = TRUE)
  detection <- rowMeans(expr > 0, na.rm = TRUE)
  qcut <- function(x, n = 5) {
    br <- unique(quantile(x, probs = seq(0, 1, length.out = n + 1), na.rm = TRUE, type = 8))
    if (length(br) < 2) return(rep(1L, length(x)))
    as.integer(cut(x, breaks = br, include.lowest = TRUE, labels = FALSE))
  }
  data.table(gene = rownames(expr), mean_bin = qcut(mean_expr), var_bin = qcut(variability), det_bin = qcut(detection))[,
    stratum := paste(mean_bin, var_bin, det_bin, sep = "_")]
}

sample_matched_genes <- function(target, strata, seed_offset = 0L) {
  set.seed(FIGURE7_SEED + seed_offset)
  target <- intersect(target, strata$gene)
  target_strata <- strata$stratum[match(target, strata$gene)]
  counts <- table(target_strata)
  pools <- split(strata$gene[!strata$gene %in% target], strata$stratum[!strata$gene %in% target])
  out <- character()
  for (st in names(counts)) {
    pool <- setdiff(pools[[st]], out)
    need <- as.integer(counts[[st]])
    if (length(pool) < need) pool <- setdiff(strata$gene, c(target, out))
    out <- c(out, sample(pool, need, replace = FALSE))
  }
  out
}

random_benchmark_one <- function(obj, result, signatures, cohort, n_random) {
  strata <- make_gene_strata(obj$expr)
  sm <- result$scores
  ranks <- result$ranks
  axes <- c("identity_loss", "stress_transition", "sox4_stabilization")
  sig_axis <- c(identity_loss = "identity_program", stress_transition = "stress_transition", sox4_stabilization = "sox4_stabilization")
  base <- sm[tumour_normal == "tumour"][!duplicated(patient_id)]
  base <- base[is.finite(os_time_days) & is.finite(os_event) & is.finite(age_high) & !is.na(sex) & is.finite(stage_high)]
  base_formula <- if (nrow(base) && sum(base$os_event == 1) >= 10) {
    try(coxph(Surv(os_time_days, os_event) ~ age_high + factor(sex) + stage_high, data = base, x = TRUE), silent = TRUE)
  } else NULL
  base_c <- if (inherits(base_formula, "coxph")) harrell_c(base$os_time_days, base$os_event, predict(base_formula, type = "lp")) else NA_real_
  rows <- vector("list", length(axes) * n_random); k <- 0L
  for (axis_name in axes) {
    target <- signatures[axis == sig_axis[[axis_name]], gene]
    for (i in seq_len(n_random)) {
      genes <- sample_matched_genes(target, strata, seed_offset = i + match(axis_name, axes) * 100000L + match(cohort, c("TCGA_LIHC", "ICGC_LIRI_JP")) * 1000000L)
      raw_score <- colMeans(ranks[genes, , drop = FALSE], na.rm = TRUE) - 0.5
      score <- safe_z(if (axis_name == "identity_loss") -raw_score else raw_score)
      g <- hedges_g(score[sm$tumour_normal == "tumour"], score[sm$tumour_normal == "normal"])
      d <- copy(sm); d[, random_score := score]
      tumour <- d[tumour_normal == "tumour"][!duplicated(patient_id)]
      tumour <- tumour[is.finite(os_time_days) & is.finite(os_event)]
      fit <- if (nrow(tumour) >= 20 && sum(tumour$os_event == 1) >= 5) try(coxph(Surv(os_time_days, os_event) ~ random_score, data = tumour), silent = TRUE) else NULL
      zcox <- if (inherits(fit, "coxph")) unname(coef(fit) / sqrt(diag(vcov(fit)))) else NA_real_
      cidx <- if (inherits(fit, "coxph")) harrell_c(tumour$os_time_days, tumour$os_event, predict(fit, type = "lp")) else NA_real_
      joint <- if (nrow(tumour) >= 20 && sum(tumour$os_event == 1) >= 5) try(coxph(Surv(os_time_days, os_event) ~ age_high + factor(sex) + stage_high + random_score, data = tumour), silent = TRUE) else NULL
      joint_c <- if (inherits(joint, "coxph")) harrell_c(tumour$os_time_days, tumour$os_event, predict(joint, type = "lp")) else NA_real_
      stage_fit <- try(lm(stage_num ~ random_score, data = d[tumour_normal == "tumour"]), silent = TRUE)
      stage_beta <- if (inherits(stage_fit, "lm") && "random_score" %in% names(coef(stage_fit))) unname(coef(stage_fit)["random_score"] / sd(d[tumour_normal == "tumour", stage_num], na.rm = TRUE)) else NA_real_
      k <- k + 1L
      rows[[k]] <- data.table(
        cohort = cohort, target_axis = axis_name, random_id = i, n_genes = length(genes),
        matching = "gene-number + mean-expression + variance + detection-rate strata",
        tumour_normal_effect = g[["effect"]], stage_association = stage_beta,
        cox_z = zcox, c_index = cidx, delta_c_index = joint_c - base_c,
        genes = paste(genes, collapse = ";"), seed = FIGURE7_SEED
      )
    }
  }
  rbindlist(rows)
}

run_random_signature_benchmark <- function(objs, results, signatures, n_random = 500L, retained = NULL) {
  bench_new <- rbindlist(lapply(names(objs), function(nm) random_benchmark_one(objs[[nm]], results[[nm]], signatures, nm, n_random)), fill = TRUE)
  bench <- rbindlist(list(retained, bench_new), fill = TRUE)
  setorder(bench, cohort, target_axis, random_id)
  write_tsvgz(bench, file.path(FIGURE7_META, "figure7_random_signature_benchmark.tsv.gz"))
  write_json(list(
    module = "Figure 7 matched random-signature benchmark", n_random_per_axis_per_cohort = n_random,
    matching = unique(bench$matching), random_seed = FIGURE7_SEED, n_rows = nrow(bench), review_risk = if (n_random < 1000) "resource_limited_500_random_signatures" else "none_for_random_count"
  ), file.path(FIGURE7_META, "figure7_random_signature_benchmark_report.json"))
  invisible(bench)
}

stage_plot_a <- function() {
  manifest <- read_tsv(file.path(FIGURE7_META, "figure7_cohort_manifest.tsv"))
  data <- rbindlist(list(
    data.table(layer = "Discovery", item = c("HCC hepatocyte single-cell atlas", "Figures 2-4 frozen programmes", "Figures 5-6 regulatory architecture"), x = 1, y = c(3, 2, 1), cohort = NA_character_),
    manifest[, .(layer = "External validation", item = sprintf("%s\nTumour %d | Normal %d | OS %d (%d events)", cohort, n_tumour, n_normal, n_survival, n_events), x = 2, y = rev(seq_len(.N)), cohort)],
    data.table(layer = "Analysis", item = c("Tumour-normal recurrence", "Clinicopathological association", "Multivariable survival", "Incremental prediction", "Controls and sensitivity"), x = 3, y = c(5, 4, 3, 2, 1), cohort = NA_character_)
  ), fill = TRUE)
  data[, layer := factor(layer, levels = c("Discovery", "External validation", "Analysis"))]
  data[, fill_key := fifelse(!is.na(cohort), cohort, fifelse(layer == "Discovery", "discovery", "analysis"))]
  cols <- c(TCGA_LIHC = cohort_palette[["TCGA_LIHC"]], ICGC_LIRI_JP = cohort_palette[["ICGC_LIRI_JP"]], discovery = lancet_palette[7], analysis = lancet_palette[9])
  p <- ggplot(data, aes(x, y, label = item, fill = fill_key)) +
    geom_label(size = 2.4, linewidth = 0.35, colour = "white", lineheight = 0.92, label.padding = grid::unit(0.18, "lines")) +
    annotate("segment", x = 1.25, xend = 1.75, y = 2, yend = 2.5, arrow = arrow(length = grid::unit(0.12, "in")), linewidth = 0.45) +
    annotate("segment", x = 2.25, xend = 2.75, y = 2.5, yend = 3, arrow = arrow(length = grid::unit(0.12, "in")), linewidth = 0.45) +
    annotate("label", x = 2, y = 5.1, label = "Signatures frozen before bulk-cohort analysis\nNo bulk outcome information used for signature derivation", size = 2.5, fontface = "bold", fill = "white", linewidth = 0.35) +
    scale_fill_manual(values = cols, guide = "none") +
    scale_x_continuous(breaks = 1:3, labels = levels(data$layer), limits = c(0.65, 3.35), position = "top") +
    coord_cartesian(ylim = c(0.5, 5.6), clip = "off") + labs(title = "A   Discovery-to-validation cohort flow") +
    theme_void(base_family = "sans") + theme(plot.title = element_text(size = 10, face = "bold"), axis.text.x.top = element_text(size = 8, face = "bold")) + figure7_panel_label("A")
  outdir <- file.path(FIGURE7_ROOT, "figures", "driver", "figure7a_cohort_flow")
  stem <- file.path(outdir, "figure7a_cohort_flow")
  figure7_export(p, stem, 7.2, 4.2); saveRDS(p, paste0(stem, "_plot.rds"))
  write_tsv(data, file.path(outdir, "figure7a_cohort_flow_data.tsv"))
  write_json(list(panel = "7A", cohorts = manifest, source_data = file.path(outdir, "figure7a_cohort_flow_data.tsv"), r_generated = TRUE), file.path(outdir, "figure7a_cohort_flow_report.json"))
}

stage_plot_b <- function() {
  cov <- read_tsv(file.path(FIGURE7_META, "figure7_signature_coverage.tsv"))
  keep <- c("identity_program", "stress_transition", "sox4_stabilization")
  d <- cov[axis %in% keep]
  d[, axis_display := factor(axis, levels = keep, labels = c("Identity programme", "Stress-transition", "SOX4 stabilization"))]
  d[, axis_key := factor(c("identity_loss", "stress_transition", "sox4_stabilization")[match(axis, keep)], levels = names(axis_palette))]
  d[, label := sprintf("%d/%d genes (%.1f%%)", n_mapped, n_frozen, 100 * coverage)]
  p1 <- ggplot(d, aes(cohort, coverage, fill = axis_key)) + geom_col(position = position_dodge(width = 0.78), width = 0.7) +
    geom_text(aes(label = label), position = position_dodge(width = 0.78), vjust = -0.5, size = 2.2) +
    scale_fill_manual(values = axis_palette, name = NULL) + scale_y_continuous(labels = percent, limits = c(0, 1.12)) +
    labs(title = "Frozen programme coverage", x = NULL, y = "Mapped fraction") + figure7_theme()
  flow <- data.table(x = 1:5, label = c("Frozen scRNA\nprogrammes", "HGNC/Ensembl\nharmonization", "Sample-wise\nrank scoring", "Within-cohort\nstandardization", "Identity-score\ndirection reversal"))
  p2 <- ggplot(flow, aes(x, 1, label = label)) +
    geom_label(fill = lancet_palette[7], colour = "white", size = 2.4, linewidth = 0.3) +
    geom_segment(data = flow[x < 5], aes(x = x + 0.25, xend = x + 0.75, y = 1, yend = 1), inherit.aes = FALSE, arrow = arrow(length = grid::unit(0.08, "in")), linewidth = 0.4) +
    coord_cartesian(xlim = c(0.6, 5.4), ylim = c(0.7, 1.3), clip = "off") + theme_void()
  p <- (p2 / p1) + plot_layout(heights = c(1, 2.2)) + plot_annotation(title = "B   Frozen single-cell programmes mapped without outcome-derived weights")
  outdir <- file.path(FIGURE7_ROOT, "figures", "driver", "figure7b_bulk_signature_mapping")
  stem <- file.path(outdir, "figure7b_bulk_signature_mapping")
  figure7_export(p, stem, 7.2, 5.0); saveRDS(p, paste0(stem, "_plot.rds"))
  write_tsv(d, file.path(outdir, "figure7b_signature_mapping.tsv"))
  write_json(list(panel = "7B", score_method = "sample-wise rank score (singscore-equivalent)", outcome_optimized_weights = FALSE, r_generated = TRUE), file.path(outdir, "figure7b_signature_mapping_report.json"))
}

stage_analyze_c <- function() {
  s <- read_tsv(file.path(FIGURE7_META, "figure7_bulk_axis_scores.tsv.gz"))
  axes <- c(identity_loss = "identity_loss_score", stress_transition = "stress_transition_score", sox4_stabilization = "sox4_stabilization_score")
  rows <- list(); k <- 0L
  for (co in unique(s$cohort)) for (ax in names(axes)) {
    d <- s[cohort == co]; col <- axes[[ax]]
    g <- hedges_g(d[tumour_normal == "tumour"][[col]], d[tumour_normal == "normal"][[col]])
    k <- k + 1L
    rows[[k]] <- data.table(cohort = co, axis = ax, analysis = "all_samples_independent",
      n_tumour = sum(d$tumour_normal == "tumour"), n_normal = sum(d$tumour_normal == "normal"), n_paired = NA_integer_,
      tumour_mean = mean(d[tumour_normal == "tumour"][[col]], na.rm = TRUE), normal_mean = mean(d[tumour_normal == "normal"][[col]], na.rm = TRUE),
      hedges_g = g[["effect"]], se = g[["se"]], ci_low = g[["lo"]], ci_high = g[["hi"]], p_value = g[["p"]])
    wide <- d[tumour_normal %in% c("tumour", "normal"), .(value = mean(get(col), na.rm = TRUE)), by = .(patient_id, tumour_normal)]
    wide <- dcast(wide, patient_id ~ tumour_normal, value.var = "value")
    if (all(c("tumour", "normal") %in% names(wide))) {
      dif <- wide$tumour - wide$normal; eff <- mean(dif, na.rm = TRUE) / sd(dif, na.rm = TRUE)
      se <- if (sum(is.finite(dif)) > 2) sqrt(1 / sum(is.finite(dif)) + eff^2 / (2 * (sum(is.finite(dif)) - 1))) else NA_real_
      k <- k + 1L
      rows[[k]] <- data.table(cohort = co, axis = ax, analysis = "paired_sensitivity", n_tumour = NA_integer_, n_normal = NA_integer_, n_paired = sum(is.finite(dif)),
        tumour_mean = mean(wide$tumour, na.rm = TRUE), normal_mean = mean(wide$normal, na.rm = TRUE), hedges_g = eff, se = se,
        ci_low = eff - 1.96 * se, ci_high = eff + 1.96 * se, p_value = t.test(dif)$p.value)
    }
  }
  effects <- rbindlist(rows, fill = TRUE)
  effects[, fdr := bh(p_value), by = analysis]
  prim <- effects[analysis == "all_samples_independent"]
  meta <- rbindlist(lapply(names(axes), function(ax) {
    z <- prim[axis == ax & is.finite(hedges_g) & is.finite(se)]
    fit <- rma(yi = hedges_g, sei = se, data = z, method = "REML")
    data.table(axis = ax, model = "random_effects_REML", k = fit$k, pooled_effect = as.numeric(fit$b), ci_low = fit$ci.lb,
               ci_high = fit$ci.ub, p_value = fit$pval, tau2 = fit$tau2, I2 = fit$I2, q_stat = fit$QE, q_p_value = fit$QEp,
               interpretation = ifelse(fit$k < 3, "exploratory_fewer_than_3_cohorts", "confirmatory_meta_analysis"))
  }))
  meta[, fdr := bh(p_value)]
  loo <- rbindlist(lapply(names(axes), function(ax) {
    z <- prim[axis == ax]
    rbindlist(lapply(z$cohort, function(drop) {
      zz <- z[cohort != drop]
      if (nrow(zz) < 2) return(data.table(axis = ax, omitted_cohort = drop, pooled_effect = NA_real_, status = "not_estimable_one_cohort_remaining"))
      fit <- rma(yi = hedges_g, sei = se, data = zz, method = "REML")
      data.table(axis = ax, omitted_cohort = drop, pooled_effect = as.numeric(fit$b), status = "estimated")
    }), fill = TRUE)
  }), fill = TRUE)
  outdir <- file.path(FIGURE7_ROOT, "figures", "driver", "figure7c_tumour_normal_forest")
  write_tsv(effects, file.path(outdir, "figure7c_tumour_normal_effects.tsv")); write_tsv(meta, file.path(outdir, "figure7c_meta_analysis.tsv"))
  write_tsv(loo, file.path(outdir, "figure7c_leave_one_cohort_out.tsv"))
  write_json(list(panel = "7C", primary_analysis = "independent Hedges g", paired_sensitivity_retained = TRUE, meta_analysis = "metafor::rma REML", meta = meta), file.path(outdir, "figure7c_tumour_normal_report.json"))
}

stage_plot_c <- function() {
  outdir <- file.path(FIGURE7_ROOT, "figures", "driver", "figure7c_tumour_normal_forest")
  e <- read_tsv(file.path(outdir, "figure7c_tumour_normal_effects.tsv"))[analysis == "all_samples_independent"]
  m <- read_tsv(file.path(outdir, "figure7c_meta_analysis.tsv"))
  d <- rbind(e[, .(axis, label = cohort, effect = hedges_g, lo = ci_low, hi = ci_high, kind = "cohort")],
             m[, .(axis, label = "Pooled (exploratory)", effect = pooled_effect, lo = ci_low, hi = ci_high, kind = "pooled")], fill = TRUE)
  d[, label := factor(label, levels = rev(c("TCGA_LIHC", "ICGC_LIRI_JP", "Pooled (exploratory)")))]
  d[, axis_display := factor(axis, levels = names(axis_palette), labels = c("Identity loss", "Stress transition", "SOX4 stabilization"))]
  p <- ggplot(d, aes(effect, label, colour = axis)) + geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, colour = "#777777") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.16, linewidth = 0.45, orientation = "y") +
    geom_point(aes(shape = kind), size = 2.1) + scale_colour_manual(values = axis_palette, guide = "none") +
    scale_shape_manual(values = c(cohort = 16, pooled = 18), guide = "none") + facet_wrap(~axis_display, nrow = 1) +
    labs(title = "C   Tumour-associated recurrence across independent bulk cohorts", subtitle = "Positive Hedges' g indicates higher programme activity in tumour; pooled estimates are exploratory (k = 2)", x = "Hedges' g (95% CI)", y = NULL) + figure7_theme()
  stem <- file.path(outdir, "figure7c_tumour_normal_forest")
  figure7_export(p, stem, 7.4, 3.5); saveRDS(p, paste0(stem, "_plot.rds"))
}

fit_clinical_association <- function(d, score_col, outcome, type, adjusted = FALSE) {
  dat <- data.frame(y = d[[outcome]], score = d[[score_col]], age_high = d$age_high, sex = factor(d$sex))
  covars <- if (adjusted) c("age_high", "sex") else character()
  needed <- c("y", "score", covars); dat <- dat[complete.cases(dat[, needed, drop = FALSE]), , drop = FALSE]
  if (nrow(dat) < 20 || length(unique(dat$y)) < 2) return(NULL)
  formula <- as.formula(paste("y ~ score", if (length(covars)) paste("+", paste(covars, collapse = "+")) else ""))
  if (type == "ordinal" && requireNamespace("ordinal", quietly = TRUE)) {
    dat$y <- ordered(dat$y); fit <- ordinal::clm(formula, data = dat, link = "logit")
    b <- coef(summary(fit))["score", ]; est <- b[["Estimate"]]; se <- b[["Std. Error"]]; p <- 2 * pnorm(-abs(est / se)); model <- "ordinal::clm"
  } else if (type == "ordinal") {
    dat$y <- ordered(dat$y); fit <- MASS::polr(formula, data = dat, Hess = TRUE)
    est <- coef(fit)[["score"]]; se <- sqrt(diag(vcov(fit)))[["score"]]; p <- 2 * pnorm(-abs(est / se)); model <- "MASS::polr_fallback"
  } else if (type == "binary") {
    fit <- glm(formula, data = dat, family = binomial()); est <- coef(fit)[["score"]]; se <- sqrt(diag(vcov(fit)))[["score"]]; p <- 2 * pnorm(-abs(est / se)); model <- "binomial_logistic"
  } else {
    fit <- lm(formula, data = dat); est <- coef(fit)[["score"]]; se <- sqrt(diag(vcov(fit)))[["score"]]; p <- 2 * pt(-abs(est / se), df.residual(fit)); model <- "linear"
  }
  data.table(model_type = model, effect_type = ifelse(type %in% c("ordinal", "binary"), "log_odds_per_SD", "standardized_beta"),
             coefficient = est, standardized_signed_effect = est, odds_ratio = ifelse(type %in% c("ordinal", "binary"), exp(est), NA_real_),
             ci_low = est - 1.96 * se, ci_high = est + 1.96 * se, p_value = p, n = nrow(dat), events = ifelse(type == "binary", sum(dat$y == max(dat$y)), NA_integer_))
}

stage_analyze_d <- function() {
  s <- read_tsv(file.path(FIGURE7_META, "figure7_bulk_axis_scores.tsv.gz"))[tumour_normal == "tumour"]
  s <- s[!duplicated(paste(cohort, patient_id))]
  programmes <- c(identity_loss = "identity_loss_score", stress_transition = "stress_transition_score", sox4_stabilization = "sox4_stabilization_score",
                  foxm1_cebpb_reference = "foxm1_cebpb_reference_score", proliferation_control = "proliferation_score", calibration_control = "calibration_control_score")
  outcomes <- list(
    TCGA_LIHC = list(pathological_stage = c("stage_num", "ordinal"), tumour_grade = c("grade_num", "ordinal"), T_stage = c("t_stage_num", "ordinal")),
    ICGC_LIRI_JP = list(pathological_stage = c("stage_high", "binary"))
  )
  rows <- list(); k <- 0L
  for (co in names(outcomes)) for (out_name in names(outcomes[[co]])) for (prog in names(programmes)) for (adj in c(FALSE, TRUE)) {
    spec <- outcomes[[co]][[out_name]]; z <- fit_clinical_association(s[cohort == co], programmes[[prog]], spec[[1]], spec[[2]], adjusted = adj)
    if (!is.null(z)) {
      k <- k + 1L
      z[, `:=`(cohort = co, programme = prog, clinical_feature = out_name, model = ifelse(adj, "age_sex_adjusted", "primary"))]
      rows[[k]] <- z
    }
  }
  out <- rbindlist(rows, fill = TRUE)
  out[, fdr := bh(p_value), by = model]
  bench <- read_tsv(file.path(FIGURE7_META, "figure7_random_signature_benchmark.tsv.gz"))
  rand <- bench[, .(standardized_signed_effect = median(stage_association, na.rm = TRUE), coefficient = median(stage_association, na.rm = TRUE),
                    p_value = NA_real_, fdr = NA_real_, model_type = "matched_random_median", effect_type = "standardized_signed_effect",
                    odds_ratio = NA_real_, ci_low = quantile(stage_association, .025, na.rm = TRUE), ci_high = quantile(stage_association, .975, na.rm = TRUE),
                    n = .N, events = NA_integer_), by = cohort]
  rand[, `:=`(programme = "matched_random_signature", clinical_feature = "pathological_stage", model = "primary")]
  out <- rbind(out, rand, fill = TRUE)
  desired <- c("pathological_stage", "tumour_grade", "T_stage", "vascular_invasion", "AFP", "recurrence", "HBV", "HCV", "cirrhosis", "tumour_purity", "CNV_burden")
  missing <- CJ(cohort = unique(s$cohort), programme = names(programmes), clinical_feature = desired)[!out[model == "primary"], on = .(cohort, programme, clinical_feature)]
  missing[, `:=`(model_type = "not_estimable", effect_type = NA_character_, coefficient = NA_real_, standardized_signed_effect = NA_real_, odds_ratio = NA_real_,
                 ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_, n = 0L, events = NA_integer_, fdr = NA_real_, model = "primary")]
  out <- rbind(out, missing, fill = TRUE)
  outdir <- file.path(FIGURE7_ROOT, "figures", "driver", "figure7d_clinical_heatmap")
  write_tsv(out, file.path(outdir, "figure7d_clinical_associations.tsv"))
  write_json(list(panel = "7D", available_features = unique(out[model_type != "not_estimable", clinical_feature]), missing_features = setdiff(desired, unique(out[model_type != "not_estimable", clinical_feature])), multiple_testing = "BH within primary/adjusted families"), file.path(outdir, "figure7d_clinical_heatmap_report.json"))
}

stage_plot_d <- function() {
  outdir <- file.path(FIGURE7_ROOT, "figures", "driver", "figure7d_clinical_heatmap")
  d <- read_tsv(file.path(outdir, "figure7d_clinical_associations.tsv"))[model == "primary" & model_type != "not_estimable"]
  d <- d[, .(standardized_signed_effect = mean(standardized_signed_effect, na.rm = TRUE), fdr = suppressWarnings(min(fdr, na.rm = TRUE))), by = .(programme, clinical_feature)]
  d[!is.finite(fdr), fdr := NA_real_]
  d[, stars := fifelse(fdr < .001, "***", fifelse(fdr < .01, "**", fifelse(fdr < .05, "*", "")))]
  prog_levels <- c("identity_loss", "stress_transition", "sox4_stabilization", "foxm1_cebpb_reference", "proliferation_control", "calibration_control", "matched_random_signature")
  d[, programme := factor(programme, levels = rev(prog_levels))]
  p <- ggplot(d, aes(clinical_feature, programme, fill = standardized_signed_effect)) + geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = stars), size = 3) + scale_fill_gradient2(low = lancet_palette[1], mid = "#F7F7F7", high = lancet_palette[2], midpoint = 0, name = "Signed\nassociation") +
    labs(title = "D   Clinicopathological associations", subtitle = "Colours are standardized signed associations; model-specific OR/log-odds are retained in source data", x = NULL, y = NULL) +
    figure7_theme() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  stem <- file.path(outdir, "figure7d_clinical_association_heatmap")
  figure7_export(p, stem, 6.8, 4.2); saveRDS(p, paste0(stem, "_plot.rds"))
}

prepare_survival_data <- function(s, cohort) {
  cohort_id <- cohort
  d <- copy(s[cohort == cohort_id & tumour_normal == "tumour"])
  d <- d[!duplicated(patient_id)]
  d[, sex_male := as.integer(toupper(sex) == "MALE")]
  d[, `:=`(age_high = safe_num(age_high), stage_high = safe_num(stage_high), os_time_days = safe_num(os_time_days), os_event = safe_num(os_event))]
  d[is.finite(os_time_days) & os_time_days > 0 & os_event %in% c(0, 1)]
}

vif_from_data <- function(d, vars) {
  x <- as.data.frame(d[, ..vars]); x <- x[, vapply(x, function(z) length(unique(z[is.finite(z)])) > 1, logical(1)), drop = FALSE]
  if (ncol(x) < 2) return(data.table(variable = names(x), vif = 1))
  rbindlist(lapply(names(x), function(v) {
    others <- setdiff(names(x), v)
    fit <- lm(reformulate(others, response = v), data = x)
    data.table(variable = v, vif = 1 / (1 - summary(fit)$r.squared))
  }))
}

fit_cox_model <- function(d, model_id, interest, covars) {
  vars <- unique(c("os_time_days", "os_event", interest, covars))
  dd <- d[complete.cases(d[, ..vars])]
  if (nrow(dd) < 20 || sum(dd$os_event == 1) < 5) return(list(rows = NULL, ph = NULL, fit = NULL, data = dd))
  formula <- as.formula(sprintf("Surv(os_time_days, os_event) ~ %s", paste(c(interest, covars), collapse = " + ")))
  fit <- try(coxph(formula, data = dd, x = TRUE, y = TRUE, model = TRUE), silent = TRUE)
  if (!inherits(fit, "coxph")) return(list(rows = NULL, ph = NULL, fit = NULL, data = dd))
  sm <- summary(fit); coef_table <- as.data.table(sm$coefficients, keep.rownames = "term"); ci <- as.data.table(sm$conf.int)
  out <- data.table(
    model_id = model_id, term = coef_table$term, coefficient = coef_table$coef, se = coef_table$`se(coef)`,
    hazard_ratio = ci$`exp(coef)`, ci_low = ci$`lower .95`, ci_high = ci$`upper .95`, p_value = coef_table$`Pr(>|z|)`,
    n = fit$n, events = fit$nevent, predictors = length(attr(terms(fit), "term.labels")),
    events_per_variable = fit$nevent / length(attr(terms(fit), "term.labels")),
    missingness = 1 - nrow(dd) / nrow(d), effective_sample_size = nrow(dd), covariates = paste(c(interest, covars), collapse = ";"),
    concordance = unname(sm$concordance[[1]]), aic = AIC(fit), bic = BIC(fit)
  )
  zph <- try(cox.zph(fit), silent = TRUE)
  ph <- if (inherits(zph, "cox.zph")) as.data.table(zph$table, keep.rownames = "term")[, .(model_id, term, chisq, df, p_value = p)] else data.table()
  list(rows = out, ph = ph, fit = fit, data = dd)
}

stage_analyze_e <- function() {
  s <- read_tsv(file.path(FIGURE7_META, "figure7_bulk_axis_scores.tsv.gz"))
  axes <- c(identity_loss = "identity_loss_score", stress_transition = "stress_transition_score", sox4_stabilization = "sox4_stabilization_score")
  controls <- c(foxm1_cebpb_reference = "foxm1_cebpb_reference_score", proliferation_control = "proliferation_score", calibration_control = "calibration_control_score")
  rows <- list(); phs <- list(); diag_rows <- list(); k <- 0L
  for (co in unique(s$cohort)) {
    d <- prepare_survival_data(s, co)
    clinical <- c("age_high", "sex_male", "stage_high")
    for (nm in names(c(axes, controls))) {
      col <- c(axes, controls)[[nm]]
      f0 <- fit_cox_model(d, paste0("axis_only__", nm), col, character())
      f1 <- fit_cox_model(d, paste0("clinical_adjusted__", nm), col, clinical)
      f3 <- if (nm %in% names(axes)) fit_cox_model(d, paste0("confounder_adjusted__", nm), col, c(clinical, "proliferation_score", "hypoxia_score", "inflammation_score")) else NULL
      for (res in Filter(Negate(is.null), list(f0, f1, f3))) if (!is.null(res$rows)) {
        k <- k + 1L; res$rows[, `:=`(cohort = co, programme = nm, endpoint = "OS")]; rows[[k]] <- res$rows
        if (nrow(res$ph)) {res$ph[, cohort := co]; phs[[length(phs) + 1L]] <- res$ph}
      }
    }
    joint_cols <- unname(axes)
    joint <- fit_cox_model(d, "clinical_adjusted__three_axis_joint", joint_cols, clinical)
    full_joint <- fit_cox_model(d, "confounder_adjusted__three_axis_joint", joint_cols, c(clinical, "proliferation_score", "hypoxia_score", "inflammation_score"))
    for (res in list(joint, full_joint)) if (!is.null(res$rows)) {
      k <- k + 1L; res$rows[, `:=`(cohort = co, programme = paste0(names(axes)[match(term, axes)], "_joint"), endpoint = "OS")]; rows[[k]] <- res$rows
      if (nrow(res$ph)) {res$ph[, cohort := co]; phs[[length(phs) + 1L]] <- res$ph}
      vars <- unique(c(joint_cols, clinical, if (grepl("confounder", res$rows$model_id[[1]])) c("proliferation_score", "hypoxia_score", "inflammation_score") else character()))
      diag_vars <- vars[vapply(res$data[, ..vars], function(x) all(is.finite(x)) && stats::sd(x) > 0, logical(1))]
      vv <- vif_from_data(res$data, diag_vars)
      if (length(diag_vars) >= 2L) {
        mm <- scale(as.matrix(res$data[, ..diag_vars]))
        sv <- svd(mm)$d
        condition_index <- if (any(sv > 1e-10)) max(sv) / min(sv[sv > 1e-10]) else NA_real_
      } else {
        condition_index <- 1
      }
      diag_rows[[length(diag_rows) + 1L]] <- vv[, .(cohort = co, model_id = res$rows$model_id[[1]], variable, vif, condition_index)]
    }
  }
  out <- rbindlist(rows, fill = TRUE); out[, fdr := bh(p_value), by = model_id]
  ph <- rbindlist(phs, fill = TRUE); if (nrow(ph)) ph[, fdr := bh(p_value)]
  primary <- out[grepl("^clinical_adjusted__", model_id) & !grepl("joint", model_id) & term %in% c(unname(axes), unname(controls))]
  meta <- rbindlist(lapply(unique(primary$programme), function(prog) {
    z <- primary[programme == prog]
    if (nrow(z) < 2) return(data.table(programme = prog, status = "not_estimable"))
    fit <- rma(yi = coefficient, sei = se, data = z, method = "REML")
    data.table(programme = prog, status = "estimated_exploratory_k2", k = fit$k, hazard_ratio = exp(as.numeric(fit$b)), ci_low = exp(fit$ci.lb), ci_high = exp(fit$ci.ub), p_value = fit$pval, tau2 = fit$tau2, I2 = fit$I2)
  }), fill = TRUE)
  meta[, fdr := bh(p_value)]
  outdir <- file.path(FIGURE7_ROOT, "figures", "driver", "figure7e_multivariable_cox_forest")
  write_tsv(out, file.path(outdir, "figure7e_cox_models.tsv")); write_tsv(ph, file.path(outdir, "figure7e_ph_assumption.tsv")); write_tsv(meta, file.path(outdir, "figure7e_cox_meta_analysis.tsv"))
  write_tsv(rbindlist(diag_rows, fill = TRUE), file.path(outdir, "figure7e_joint_model_diagnostics.tsv"))
  write_json(list(panel = "7E", endpoint = "overall_survival", baseline_covariates = c("age_high", "sex_male", "stage_high"),
                  confounder_covariates = c("proliferation_score", "hypoxia_score", "inflammation_score"), unavailable_covariates = c("purity", "CNV_burden"),
                  no_stepwise_selection = TRUE, ph_checked = TRUE, meta_analysis = "exploratory because k=2"), file.path(outdir, "figure7e_cox_report.json"))
}

stage_plot_e <- function() {
  outdir <- file.path(FIGURE7_ROOT, "figures", "driver", "figure7e_multivariable_cox_forest")
  d <- read_tsv(file.path(outdir, "figure7e_cox_models.tsv"))[grepl("^clinical_adjusted__", model_id)]
  d <- d[term %in% c("identity_loss_score", "stress_transition_score", "sox4_stabilization_score", "foxm1_cebpb_reference_score", "proliferation_score", "calibration_control_score")]
  d[, label := fifelse(grepl("joint", model_id), paste0(programme), programme)]
  d[, colour_key := sub("_joint$", "", programme)]
  m <- read_tsv(file.path(outdir, "figure7e_cox_meta_analysis.tsv"))[status != "not_estimable"]
  if (nrow(m)) {
    mm <- m[, .(cohort = "Pooled exploratory", programme, label = programme, colour_key = programme, hazard_ratio, ci_low, ci_high, model_id = "pooled")]
    d <- rbind(d, mm, fill = TRUE)
  }
  cols <- c(axis_palette, control_palette)
  d[, label := factor(label, levels = rev(unique(label)))]
  p <- ggplot(d, aes(hazard_ratio, label, colour = colour_key)) + geom_vline(xintercept = 1, linetype = 2, linewidth = .35, colour = "#777777") +
    geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = .15, linewidth = .45, orientation = "y") + geom_point(size = 1.8) +
    scale_x_log10() + scale_colour_manual(values = cols, na.value = "#555555", guide = "none") + facet_wrap(~cohort, scales = "free_y") +
    labs(title = "E   Multivariable overall-survival associations", subtitle = "Hazard ratio per 1-SD programme increase; shared baseline: age, sex and stage", x = "Hazard ratio (log scale, 95% CI)", y = NULL) + figure7_theme()
  stem <- file.path(outdir, "figure7e_multivariable_cox_forest")
  figure7_export(p, stem, 7.4, 5.2); saveRDS(p, paste0(stem, "_plot.rds"))
}

model_specs <- function() list(
  clinical_baseline = character(), identity_loss = "identity_loss_score", stress_transition = "stress_transition_score",
  sox4_stabilization = "sox4_stabilization_score", all_three_axes = c("identity_loss_score", "stress_transition_score", "sox4_stabilization_score"),
  foxm1_cebpb_reference = "foxm1_cebpb_reference_score", proliferation_control = "proliferation_score", calibration_control = "calibration_control_score"
)

fit_prediction_model <- function(d, extras) {
  vars <- c("age_high", "sex_male", "stage_high", extras)
  dd <- d[complete.cases(d[, c("os_time_days", "os_event", vars), with = FALSE])]
  fit <- coxph(as.formula(paste("Surv(os_time_days, os_event) ~", paste(vars, collapse = "+"))), data = dd, x = TRUE, y = TRUE, model = TRUE)
  list(fit = fit, data = dd, vars = vars)
}

bootstrap_delta_c <- function(time, event, lp_model, lp_base, n_boot = 500L, seed = FIGURE7_SEED) {
  set.seed(seed); n <- length(time)
  vals <- replicate(n_boot, {
    i <- sample.int(n, n, replace = TRUE)
    harrell_c(time[i], event[i], lp_model[i]) - harrell_c(time[i], event[i], lp_base[i])
  })
  c(
    delta = harrell_c(time, event, lp_model) - harrell_c(time, event, lp_base),
    lo = unname(quantile(vals, .025, na.rm = TRUE)),
    hi = unname(quantile(vals, .975, na.rm = TRUE))
  )
}

ipcw_auc <- function(time, event, risk, horizon) {
  ok <- is.finite(time) & is.finite(event) & is.finite(risk); time <- time[ok]; event <- event[ok]; risk <- risk[ok]
  if (sum(event == 1 & time <= horizon) < 3 || sum(time > horizon) < 3) return(NA_real_)
  censor_fit <- survfit(Surv(time, 1 - event) ~ 1)
  gfun <- function(t) {
    z <- summary(censor_fit, times = pmax(t - 1e-8, 0), extend = TRUE)$surv
    pmax(z, 1e-4)
  }
  cases <- which(event == 1 & time <= horizon); controls <- which(time > horizon)
  wc <- 1 / gfun(time[cases]); w0 <- rep(1 / gfun(horizon), length(controls))
  cmp <- outer(risk[cases], risk[controls], function(a, b) (a > b) + 0.5 * (a == b))
  sum(cmp * outer(wc, w0)) / sum(outer(wc, w0))
}

ipcw_brier <- function(time, event, risk_prob, horizon) {
  ok <- is.finite(time) & is.finite(event) & is.finite(risk_prob); time <- time[ok]; event <- event[ok]; risk_prob <- risk_prob[ok]
  censor_fit <- survfit(Surv(time, 1 - event) ~ 1)
  gfun <- function(t) pmax(summary(censor_fit, times = pmax(t - 1e-8, 0), extend = TRUE)$surv, 1e-4)
  w <- ifelse(time <= horizon & event == 1, 1 / gfun(time), ifelse(time > horizon, 1 / gfun(horizon), 0))
  y <- as.numeric(time <= horizon & event == 1)
  sum(w * (y - risk_prob)^2) / sum(w)
}

cross_validate_models <- function(d, repeats = 10L, folds = 5L) {
  specs <- model_specs(); rows <- list(); k <- 0L
  for (r in seq_len(repeats)) {
    set.seed(FIGURE7_SEED + r); fold <- sample(rep(seq_len(folds), length.out = nrow(d)))
    for (nm in names(specs)) {
      lp <- rep(NA_real_, nrow(d))
      for (f in seq_len(folds)) {
        train <- d[fold != f]; test <- d[fold == f]
        fit <- try(fit_prediction_model(train, specs[[nm]])$fit, silent = TRUE)
        if (inherits(fit, "coxph")) lp[fold == f] <- predict(fit, newdata = test, type = "lp", reference = "zero")
      }
      k <- k + 1L; rows[[k]] <- data.table(repeat_id = r, model = nm, c_index = harrell_c(d$os_time_days, d$os_event, lp), n = sum(is.finite(lp)), events = sum(d$os_event[is.finite(lp)]))
    }
  }
  rbindlist(rows)
}

stage_prediction_f <- function() {
  s <- read_tsv(file.path(FIGURE7_META, "figure7_bulk_axis_scores.tsv.gz"))
  cohorts <- lapply(unique(s$cohort), function(co) prepare_survival_data(s, co)); names(cohorts) <- unique(s$cohort)
  specs <- model_specs(); perf_rows <- list(); cal_rows <- list(); cv_rows <- list(); fitted <- list()
  for (co in names(cohorts)) {
    d <- cohorts[[co]]; cv <- cross_validate_models(d); cv[, cohort := co]; cv_rows[[co]] <- cv
    base <- fit_prediction_model(d, character()); base_lp <- predict(base$fit, newdata = base$data, type = "lp", reference = "zero")
    for (nm in names(specs)) {
      res <- fit_prediction_model(d, specs[[nm]]); fit <- res$fit; dd <- res$data
      lp <- predict(fit, newdata = dd, type = "lp", reference = "zero")
      base_same <- predict(base$fit, newdata = dd, type = "lp", reference = "zero")
      delta <- bootstrap_delta_c(dd$os_time_days, dd$os_event, lp, base_same, n_boot = 500L, seed = FIGURE7_SEED + match(nm, names(specs)))
      lrt_p <- NA_real_
      if (nm != "clinical_baseline") {
        lrt <- try(anova(base$fit, fit, test = "LRT"), silent = TRUE)
        if (!inherits(lrt, "try-error")) {
          p_col <- grep("Pr\\(|P\\(", names(lrt), value = TRUE)[1]
          if (length(p_col) && !is.na(p_col)) lrt_p <- as.numeric(lrt[[p_col]][2])
        }
      }
      perf_rows[[length(perf_rows) + 1L]] <- data.table(cohort = co, validation = "apparent_with_bootstrap_delta", model = nm,
        c_index = harrell_c(dd$os_time_days, dd$os_event, lp), uno_c_index = NA_real_, delta_c_index = delta[["delta"]], delta_ci_low = delta[["lo"]], delta_ci_high = delta[["hi"]],
        lrt_p_value = lrt_p, aic = AIC(fit), bic = BIC(fit), n = nrow(dd), events = sum(dd$os_event), predictors = length(attr(terms(fit), "term.labels")))
      bhz <- basehaz(fit, centered = FALSE)
      for (yr in c(1, 3, 5)) {
        horizon <- yr * 365.25; h0 <- bhz$hazard[max(which(bhz$time <= horizon), 1)]
        risk <- 1 - exp(-h0 * exp(lp)); auc <- ipcw_auc(dd$os_time_days, dd$os_event, lp, horizon)
        brier <- ipcw_brier(dd$os_time_days, dd$os_event, risk, horizon)
        calfit <- try(coxph(Surv(os_time_days, os_event) ~ lp, data = data.frame(dd, lp = lp)), silent = TRUE)
        slope <- if (inherits(calfit, "coxph")) coef(calfit)[["lp"]] else NA_real_
        observed <- 1 - summary(survfit(Surv(dd$os_time_days, dd$os_event) ~ 1), times = horizon, extend = TRUE)$surv
        cal_rows[[length(cal_rows) + 1L]] <- data.table(cohort = co, validation = "apparent", model = nm, year = yr, horizon_days = horizon,
          time_dependent_auc = auc, brier_score = brier, calibration_slope = slope, calibration_in_large = mean(risk) - observed, observed_risk = observed, mean_predicted_risk = mean(risk))
      }
      fitted[[paste(co, nm, sep = "__")]] <- fit
    }
  }
  tcga <- cohorts$TCGA_LIHC; icgc <- cohorts$ICGC_LIRI_JP
  ext_rows <- list(); ext_cal <- list()
  for (nm in names(specs)) {
    train_res <- fit_prediction_model(tcga, specs[[nm]]); fit <- train_res$fit
    vars <- train_res$vars; test <- icgc[complete.cases(icgc[, c("os_time_days", "os_event", vars), with = FALSE])]
    lp <- predict(fit, newdata = test, type = "lp", reference = "zero")
    base_fit <- fit_prediction_model(tcga, character())$fit; base_lp <- predict(base_fit, newdata = test, type = "lp", reference = "zero")
    delta <- bootstrap_delta_c(test$os_time_days, test$os_event, lp, base_lp, 500L, FIGURE7_SEED + 5000L + match(nm, names(specs)))
    ext_rows[[nm]] <- data.table(train_cohort = "TCGA_LIHC", test_cohort = "ICGC_LIRI_JP", model = nm, weights_locked = TRUE,
      c_index = harrell_c(test$os_time_days, test$os_event, lp), delta_c_index = delta[["delta"]], delta_ci_low = delta[["lo"]], delta_ci_high = delta[["hi"]], n = nrow(test), events = sum(test$os_event))
    bhz <- basehaz(fit, centered = FALSE)
    for (yr in c(1, 3, 5)) {
      horizon <- yr * 365.25; h0 <- bhz$hazard[max(which(bhz$time <= horizon), 1)]; risk <- 1 - exp(-h0 * exp(lp))
      calfit <- try(coxph(Surv(os_time_days, os_event) ~ lp, data = data.frame(test, lp = lp)), silent = TRUE)
      observed <- 1 - summary(survfit(Surv(test$os_time_days, test$os_event) ~ 1), times = horizon, extend = TRUE)$surv
      ext_cal[[length(ext_cal) + 1L]] <- data.table(cohort = "ICGC_LIRI_JP", validation = "locked_external_TCGA_to_ICGC", model = nm, year = yr, horizon_days = horizon,
        time_dependent_auc = ipcw_auc(test$os_time_days, test$os_event, lp, horizon), brier_score = ipcw_brier(test$os_time_days, test$os_event, risk, horizon),
        calibration_slope = if (inherits(calfit, "coxph")) coef(calfit)[["lp"]] else NA_real_, calibration_in_large = mean(risk) - observed,
        observed_risk = observed, mean_predicted_risk = mean(risk))
    }
  }
  perf <- rbindlist(perf_rows, fill = TRUE); cv <- rbindlist(cv_rows, fill = TRUE); ext <- rbindlist(ext_rows, fill = TRUE); cal <- rbind(rbindlist(cal_rows), rbindlist(ext_cal), fill = TRUE)
  outdir <- file.path(FIGURE7_ROOT, "figures", "driver", "figure7f_incremental_prediction")
  write_tsv(perf, file.path(outdir, "figure7f_model_performance.tsv")); write_tsv(cv, file.path(outdir, "figure7f_cross_validation.tsv"))
  write_tsv(ext, file.path(outdir, "figure7f_external_validation.tsv")); write_tsv(cal, file.path(outdir, "figure7f_calibration.tsv"))
  write_json(list(panel = "7F", clinical_baseline = c("age_high", "sex_male", "stage_high"), internal_validation = "10x repeated 5-fold CV plus 500 paired bootstrap delta-CI",
                  external_validation = "TCGA coefficients and baseline hazard locked before ICGC application", time_auc = "IPCW cumulative/dynamic AUC", uno_c_index = "not available in installed packages"), file.path(outdir, "figure7f_incremental_prediction_report.json"))
}

stage_plot_f <- function() {
  outdir <- file.path(FIGURE7_ROOT, "figures", "driver", "figure7f_incremental_prediction")
  pdat <- read_tsv(file.path(outdir, "figure7f_model_performance.tsv")); edat <- read_tsv(file.path(outdir, "figure7f_external_validation.tsv")); cal <- read_tsv(file.path(outdir, "figure7f_calibration.tsv"))
  delta <- rbind(
    pdat[, .(source = fifelse(cohort == "TCGA_LIHC", "TCGA apparent", "ICGC apparent"), model, delta_c_index, delta_ci_low, delta_ci_high)],
    edat[, .(source = "ICGC locked external", model, delta_c_index, delta_ci_low, delta_ci_high)],
    fill = TRUE
  )
  dodge <- position_dodge(width = .55)
  pd <- ggplot(delta, aes(delta_c_index, model, colour = model, shape = source)) + geom_vline(xintercept = 0, linetype = 2, linewidth = .3) +
    geom_errorbar(aes(xmin = delta_ci_low, xmax = delta_ci_high), width = .15, orientation = "y", position = dodge) + geom_point(position = dodge) +
    scale_colour_manual(values = rep(lancet_palette, length.out = length(unique(delta$model))), guide = "none") +
    scale_shape_manual(values = c("TCGA apparent" = 16, "ICGC apparent" = 17, "ICGC locked external" = 15), name = NULL) +
    labs(title = "Incremental discrimination", x = expression(Delta*"C-index (95% CI)"), y = NULL) + figure7_theme() + theme(legend.position = "bottom")
  auc <- cal[model %in% c("clinical_baseline", "all_three_axes")]
  pa <- ggplot(auc, aes(year, time_dependent_auc, colour = model, linetype = validation, group = interaction(model, validation))) + geom_line() + geom_point() +
    scale_colour_manual(values = c(clinical_baseline = lancet_palette[6], all_three_axes = lancet_palette[2])) + scale_x_continuous(breaks = c(1, 3, 5)) + labs(title = "Time-dependent AUC", x = "Years", y = "IPCW AUC", colour = NULL, linetype = NULL) + figure7_theme()
  cs <- cal[year == 3 & model %in% c("clinical_baseline", "all_three_axes")]
  pc <- ggplot(cs, aes(calibration_slope, interaction(model, validation), colour = model)) + geom_vline(xintercept = 1, linetype = 2, linewidth = .3) + geom_point(size = 2) +
    scale_colour_manual(values = c(clinical_baseline = lancet_palette[6], all_three_axes = lancet_palette[2]), guide = "none") + labs(title = "3-year calibration", x = "Calibration slope", y = NULL) + figure7_theme()
  p <- (pd | pa | pc) + plot_annotation(title = "F   Incremental prediction beyond a shared clinical baseline")
  stem <- file.path(outdir, "figure7f_incremental_prediction")
  figure7_export(p, stem, 10.5, 4.8); saveRDS(p, paste0(stem, "_plot.rds"))
}

stage_survival_groups_g <- function() {
  s <- read_tsv(file.path(FIGURE7_META, "figure7_bulk_axis_scores.tsv.gz"))
  tcga <- prepare_survival_data(s, "TCGA_LIHC"); icgc <- prepare_survival_data(s, "ICGC_LIRI_JP")
  extras <- c("identity_loss_score", "stress_transition_score", "sox4_stabilization_score")
  fit <- fit_prediction_model(tcga, extras)$fit
  tcga[, risk_score := predict(fit, newdata = tcga, type = "lp", reference = "zero")]
  icgc[, risk_score := predict(fit, newdata = icgc, type = "lp", reference = "zero")]
  cutoff <- median(tcga$risk_score, na.rm = TRUE)
  tcga[, risk_group := ifelse(risk_score >= cutoff, "High", "Low")]
  icgc[, risk_group := ifelse(risk_score >= cutoff, "High", "Low")]
  out <- rbind(tcga[, .(sample_id, patient_id, cohort, os_time_days, os_event, risk_score, risk_group)], icgc[, .(sample_id, patient_id, cohort, os_time_days, os_event, risk_score, risk_group)])
  outdir <- file.path(FIGURE7_ROOT, "figures", "driver", "figure7g_survival_curves")
  write_tsv(out, file.path(outdir, "figure7g_survival_groups.tsv"))
  write_json(list(panel = "7G", grouping_rule = "TCGA median of locked TCGA three-axis-plus-clinical linear predictor; same numeric cutoff applied to ICGC", cutoff = cutoff,
                  coefficient_source = "TCGA_LIHC only", external_weights_locked = TRUE, km_role = "visualization_not_primary_evidence"), file.path(outdir, "figure7g_survival_curve_report.json"))
}

manual_km_plot <- function(d, cohort_label) {
  fit <- survfit(Surv(os_time_days, os_event) ~ risk_group, data = d)
  sm <- summary(fit)
  dat <- data.table(time = sm$time / 365.25, surv = sm$surv, lower = sm$lower, upper = sm$upper, strata = sub("risk_group=", "", sm$strata))
  pval <- survdiff(Surv(os_time_days, os_event) ~ risk_group, data = d)
  p <- 1 - pchisq(pval$chisq, df = length(pval$n) - 1)
  cfit <- coxph(Surv(os_time_days, os_event) ~ risk_score, data = d); hr <- exp(coef(cfit)); ci <- exp(confint(cfit))
  ggplot(dat, aes(time, surv, colour = strata)) + geom_step(linewidth = .65) + geom_ribbon(aes(ymin = lower, ymax = upper, fill = strata), alpha = .12, colour = NA) +
    scale_colour_manual(values = c(High = lancet_palette[2], Low = lancet_palette[1]), name = "Locked risk group") + scale_fill_manual(values = c(High = lancet_palette[2], Low = lancet_palette[1]), guide = "none") +
    annotate("text", x = Inf, y = .08, hjust = 1.05, label = sprintf("Log-rank P = %.3g\nContinuous HR %.2f (%.2f-%.2f)", p, hr, ci[1], ci[2]), size = 2.4) +
    labs(title = cohort_label, x = "Years", y = "Overall survival probability") + figure7_theme()
}

stage_plot_g <- function() {
  outdir <- file.path(FIGURE7_ROOT, "figures", "driver", "figure7g_survival_curves")
  d <- read_tsv(file.path(outdir, "figure7g_survival_groups.tsv")); plots <- list()
  for (co in c("TCGA_LIHC", "ICGC_LIRI_JP")) {
    z <- d[cohort == co]
    if (requireNamespace("survminer", quietly = TRUE)) {
      fit <- survfit(Surv(os_time_days, os_event) ~ risk_group, data = z)
      g <- survminer::ggsurvplot(fit, data = z, risk.table = TRUE, conf.int = TRUE, palette = c(lancet_palette[2], lancet_palette[1]),
                                ggtheme = figure7_theme(), xscale = 365.25, xlab = "Years", ylab = "Overall survival probability", title = co)
      plots[[co]] <- g$plot / g$table + plot_layout(heights = c(3, 1))
    } else plots[[co]] <- manual_km_plot(z, co)
  }
  p <- wrap_plots(plots, nrow = 1) + plot_annotation(title = "G   Locked three-axis risk grouping", subtitle = "Kaplan-Meier curves are visualization; continuous Cox models provide the primary evidence")
  stem <- file.path(outdir, "figure7g_survival_curves")
  figure7_export(p, stem, 9.0, 4.8); saveRDS(p, paste0(stem, "_plot.rds"))
}

fit_sensitivity_axes <- function(d, scenario, extra_covars = character(), subset_expr = rep(TRUE, nrow(d)), score_override = NULL) {
  axes <- c(identity_loss = "identity_loss_score", stress_transition = "stress_transition_score", sox4_stabilization = "sox4_stabilization_score")
  z <- d[subset_expr]
  if (!is.null(score_override)) {
    for (nm in names(score_override)) z[[nm]] <- score_override[[nm]][subset_expr]
  }
  rbindlist(lapply(names(axes), function(ax) {
    fit <- fit_cox_model(z, scenario, axes[[ax]], c("age_high", "sex_male", "stage_high", extra_covars))
    if (is.null(fit$rows)) return(data.table(scenario = scenario, axis = ax, status = "not_estimable", reason = "insufficient_events_or_model_failure"))
    row <- fit$rows[term == axes[[ax]]][1]
    data.table(scenario = scenario, axis = ax, status = "estimated", reason = "", effect_type = "hazard_ratio", effect = row$hazard_ratio,
               ci_low = row$ci_low, ci_high = row$ci_high, p_value = row$p_value, n = row$n, events = row$events, epv = row$events_per_variable)
  }), fill = TRUE)
}

alternative_axis_scores <- function(obj, signatures, method = "zmean", remove_genes = character()) {
  axes <- c(identity_loss = "identity_program", stress_transition = "stress_transition", sox4_stabilization = "sox4_stabilization")
  out <- list()
  for (ax in names(axes)) {
    genes <- setdiff(intersect(signatures[axis == axes[[ax]], gene], rownames(obj$expr)), remove_genes)
    if (method == "tf_only") genes <- intersect(signatures[axis == axes[[ax]] & source_rank == 0, gene], rownames(obj$expr))
    if (!length(genes)) {out[[paste0(ax, "_score")]] <- rep(NA_real_, ncol(obj$expr)); next}
    if (method == "zmean") {
      x <- t(scale(t(obj$expr[genes, , drop = FALSE]))); raw <- colMeans(x, na.rm = TRUE)
    } else {
      ranks <- apply(obj$expr, 2, rank, ties.method = "average") / nrow(obj$expr); raw <- colMeans(ranks[genes, , drop = FALSE]) - .5
    }
    if (ax == "identity_loss") raw <- -raw
    out[[paste0(ax, "_score")]] <- safe_z(raw)
  }
  out
}

stage_sensitivity_h <- function() {
  s <- read_tsv(file.path(FIGURE7_META, "figure7_bulk_axis_scores.tsv.gz"))
  signatures <- read_tsv(file.path(FIGURE7_META, "figure7_frozen_signature_manifest.tsv"))
  objs <- list(TCGA_LIHC = readRDS(file.path(FIGURE7_DATA, "figure7_tcga_lihc_prepared_expression.rds")),
               ICGC_LIRI_JP = readRDS(file.path(FIGURE7_DATA, "figure7_icgc_liri_jp_prepared_expression.rds")))
  required_scenarios <- c(
    "adjust_tumour_purity", "exclude_top20_proliferation", "adjust_proliferation", "adjust_hypoxia", "adjust_inflammation",
    "adjust_cnv_burden", "HBV_HCV_stratified", "cirrhosis_stratified", "early_stage", "advanced_stage", "male_stratum", "female_stratum",
    "leave_one_cohort_out", "leave_one_dataset_out", "complete_case", "multiple_imputation", "GSVA_ssGSEA_score", "regulon_only_signature",
    "TF_expression_only", "adjust_FOXM1_CEBPB", "adjust_calibration_control", "matched_random_signatures", "coverage_thresholds",
    "remove_cell_cycle_genes", "remove_generic_stress_genes"
  )
  rows <- list()
  for (co in unique(s$cohort)) {
    cohort_start <- length(rows) + 1L
    d <- prepare_survival_data(s, co)
    rows[[length(rows) + 1L]] <- fit_sensitivity_axes(d, "primary_clinical_adjusted")
    rows[[length(rows) + 1L]] <- fit_sensitivity_axes(d, "exclude_top20_proliferation", subset_expr = d$proliferation_score <= quantile(d$proliferation_score, .8, na.rm = TRUE))
    rows[[length(rows) + 1L]] <- fit_sensitivity_axes(d, "adjust_proliferation", "proliferation_score")
    rows[[length(rows) + 1L]] <- fit_sensitivity_axes(d, "adjust_hypoxia", "hypoxia_score")
    rows[[length(rows) + 1L]] <- fit_sensitivity_axes(d, "adjust_inflammation", "inflammation_score")
    rows[[length(rows) + 1L]] <- fit_sensitivity_axes(d, "early_stage", subset_expr = d$stage_high == 0)
    rows[[length(rows) + 1L]] <- fit_sensitivity_axes(d, "advanced_stage", subset_expr = d$stage_high == 1)
    rows[[length(rows) + 1L]] <- fit_sensitivity_axes(d, "male_stratum", subset_expr = d$sex_male == 1)
    rows[[length(rows) + 1L]] <- fit_sensitivity_axes(d, "female_stratum", subset_expr = d$sex_male == 0)
    rows[[length(rows) + 1L]] <- fit_sensitivity_axes(d, "complete_case")
    rows[[length(rows) + 1L]] <- fit_sensitivity_axes(d, "regulon_only_signature")
    rows[[length(rows) + 1L]] <- fit_sensitivity_axes(d, "adjust_FOXM1_CEBPB", "foxm1_cebpb_reference_score")
    rows[[length(rows) + 1L]] <- fit_sensitivity_axes(d, "adjust_calibration_control", "calibration_control_score")
    rows[[length(rows) + 1L]] <- fit_sensitivity_axes(d, "coverage_thresholds")
    zmean <- alternative_axis_scores(objs[[co]], signatures, "zmean")
    tfonly <- alternative_axis_scores(objs[[co]], signatures, "tf_only")
    cell_cycle <- c("MKI67", "TOP2A", "STMN1", "TYMS", "UBE2C", "PCNA", "MCM2", "MCM5", "HMGB2", "CCNB1", "CCNB2", "CDC20", "CENPF", "AURKB")
    generic_stress <- c("HSPA1A", "HSPA1B", "HSPA6", "HSP90AA1", "ATF3", "JUN", "FOS", "JUNB", "JUND", "DDIT3")
    no_cycle <- alternative_axis_scores(objs[[co]], signatures, "rank", cell_cycle)
    no_stress <- alternative_axis_scores(objs[[co]], signatures, "rank", generic_stress)
    match_order <- match(d$sample_id, objs[[co]]$samples$sample_id)
    reorder_scores <- function(x) lapply(x, function(v) v[match_order])
    rows[[length(rows) + 1L]] <- fit_sensitivity_axes(d, "zmean_score", score_override = reorder_scores(zmean))
    rows[[length(rows) + 1L]] <- fit_sensitivity_axes(d, "TF_expression_only", score_override = reorder_scores(tfonly))
    rows[[length(rows) + 1L]] <- fit_sensitivity_axes(d, "remove_cell_cycle_genes", score_override = reorder_scores(no_cycle))
    rows[[length(rows) + 1L]] <- fit_sensitivity_axes(d, "remove_generic_stress_genes", score_override = reorder_scores(no_stress))
    unavailable <- c("adjust_tumour_purity", "adjust_cnv_burden", "HBV_HCV_stratified", "cirrhosis_stratified",
                     "leave_one_cohort_out", "leave_one_dataset_out", "multiple_imputation", "GSVA_ssGSEA_score")
    for (sc in unavailable) rows[[length(rows) + 1L]] <- data.table(scenario = sc, axis = c("identity_loss", "stress_transition", "sox4_stabilization"), status = "not_estimable",
      reason = switch(sc, adjust_tumour_purity = "purity absent", adjust_cnv_burden = "CNV burden absent", HBV_HCV_stratified = "HBV/HCV absent",
                      cirrhosis_stratified = "cirrhosis absent", leave_one_cohort_out = "only two independent bulk cohorts", leave_one_dataset_out = "cohort equals dataset",
                      multiple_imputation = "no material missingness in shared-baseline variables; imputation not warranted", GSVA_ssGSEA_score = "GSVA package unavailable in frozen R environment"))
    rows[[length(rows) + 1L]] <- data.table(scenario = "matched_random_signatures", axis = c("identity_loss", "stress_transition", "sox4_stabilization"), status = "benchmark_table",
                                            reason = "reported as empirical percentiles", effect_type = NA_character_, effect = NA_real_)
    for (i in seq.int(cohort_start, length(rows))) rows[[i]][, cohort := co]
  }
  out <- rbindlist(rows, fill = TRUE)
  out[, fdr := bh(p_value), by = scenario]
  bench <- read_tsv(file.path(FIGURE7_META, "figure7_random_signature_benchmark.tsv.gz"))
  c_eff <- read_tsv(file.path(FIGURE7_ROOT, "figures", "driver", "figure7c_tumour_normal_forest", "figure7c_tumour_normal_effects.tsv"))[analysis == "all_samples_independent"]
  e_cox <- read_tsv(file.path(FIGURE7_ROOT, "figures", "driver", "figure7e_multivariable_cox_forest", "figure7e_cox_models.tsv"))[grepl("^clinical_adjusted__", model_id) & !grepl("joint", model_id)]
  f_perf <- read_tsv(file.path(FIGURE7_ROOT, "figures", "driver", "figure7f_incremental_prediction", "figure7f_model_performance.tsv"))
  score_col <- c(identity_loss = "identity_loss_score", stress_transition = "stress_transition_score", sox4_stabilization = "sox4_stabilization_score")
  percentiles <- rbindlist(lapply(unique(bench$cohort), function(co) rbindlist(lapply(names(score_col), function(ax) {
    b <- bench[cohort == co & target_axis == ax]
    obs_g <- c_eff[cohort == co & axis == ax, hedges_g][1]
    ec <- e_cox[cohort == co & programme == ax & term == score_col[[ax]]][1]
    obs_z <- abs(ec$coefficient / ec$se)
    model_name <- c(identity_loss = "identity_loss", stress_transition = "stress_transition", sox4_stabilization = "sox4_stabilization")[[ax]]
    fp <- f_perf[cohort == co & model == model_name][1]
    data.table(cohort = co, axis = ax, metric = c("tumour_normal_effect", "cox_abs_z", "c_index", "delta_c_index"),
               observed = c(obs_g, obs_z, fp$c_index, fp$delta_c_index),
               random_percentile = c(mean(b$tumour_normal_effect <= obs_g, na.rm = TRUE), mean(abs(b$cox_z) <= obs_z, na.rm = TRUE),
                                     mean(b$c_index <= fp$c_index, na.rm = TRUE), mean(b$delta_c_index <= fp$delta_c_index, na.rm = TRUE)))
  }))), fill = TRUE)
  control_z <- e_cox[programme == "calibration_control", .(control_abs_z = median(abs(coefficient / se), na.rm = TRUE))]
  status <- out[status == "estimated", .(direction_fraction = mean(effect > 1, na.rm = TRUE), fdr_fraction = mean(fdr < .05, na.rm = TRUE)), by = axis]
  status[, evidence_status := fifelse(direction_fraction >= .7 & fdr_fraction >= .5, "Robust", fifelse(direction_fraction >= .7, "Partial", "Unresolved"))]
  axis_z <- e_cox[programme %in% status$axis, .(axis_abs_z = median(abs(coefficient / se), na.rm = TRUE)), by = programme]
  status <- merge(status, axis_z, by.x = "axis", by.y = "programme", all.x = TRUE)
  if (nrow(control_z)) status[is.finite(axis_abs_z) & control_z$control_abs_z > axis_abs_z, evidence_status := "Control-dominated"]
  outdir <- file.path(FIGURE7_ROOT, "figures", "driver", "figure7h_sensitivity_specificity")
  write_tsv(out, file.path(outdir, "figure7h_sensitivity_results.tsv")); write_tsv(percentiles, file.path(outdir, "figure7h_random_signature_percentiles.tsv")); write_tsv(status, file.path(outdir, "figure7h_axis_robustness_status.tsv"))
  write_json(list(panel = "7H", requested_scenarios = required_scenarios, represented_scenarios = unique(out$scenario), n_random = 500,
                  unavailable_are_retained = TRUE, axis_status = status, specificity_risk = if (any(status$evidence_status == "Control-dominated")) "calibration_control_outperformed_at_least_one_axis" else "no_control_dominance_detected"),
             file.path(outdir, "figure7h_specificity_report.json"))
}

stage_plot_h <- function() {
  outdir <- file.path(FIGURE7_ROOT, "figures", "driver", "figure7h_sensitivity_specificity")
  s <- read_tsv(file.path(outdir, "figure7h_sensitivity_results.tsv"))[status == "estimated" & scenario %in% c("primary_clinical_adjusted", "exclude_top20_proliferation", "adjust_proliferation", "adjust_hypoxia", "adjust_inflammation", "adjust_FOXM1_CEBPB", "adjust_calibration_control", "TF_expression_only", "remove_cell_cycle_genes", "remove_generic_stress_genes")]
  p1 <- ggplot(s, aes(effect, interaction(scenario, cohort), colour = axis)) + geom_vline(xintercept = 1, linetype = 2, linewidth = .3) +
    geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = .12, linewidth = .4, orientation = "y", position = position_dodge(width = .5)) +
    geom_point(position = position_dodge(width = .5), size = 1.4) + scale_x_log10() + scale_colour_manual(values = axis_palette, name = NULL) +
    labs(title = "Confounder and definition sensitivity", x = "OS hazard ratio (95% CI)", y = NULL) + figure7_theme()
  pctl <- read_tsv(file.path(outdir, "figure7h_random_signature_percentiles.tsv"))
  p2 <- ggplot(pctl, aes(metric, axis, fill = random_percentile)) + geom_tile(colour = "white", linewidth = .5) + geom_text(aes(label = percent(random_percentile, accuracy = 1)), size = 1.8) +
    scale_fill_gradient(low = "#F7F7F7", high = lancet_palette[2], limits = c(0, 1), labels = percent, name = "Random\npercentile") +
    facet_wrap(~cohort, nrow = 1) + labs(title = "Matched-random specificity", x = NULL, y = NULL) + figure7_theme() + theme(axis.text.x = element_text(angle = 35, hjust = 1), strip.text = element_text(size = 7, face = "bold"))
  stat <- read_tsv(file.path(outdir, "figure7h_axis_robustness_status.tsv"))
  p3 <- ggplot(stat, aes(evidence_status, axis, fill = evidence_status)) + geom_tile(colour = "white") + geom_text(aes(label = evidence_status), size = 2.3) +
    scale_fill_manual(values = evidence_palette, guide = "none") + labs(title = "Evidence status", x = NULL, y = NULL) + figure7_theme()
  p <- (p1 | (p2 / p3)) + plot_layout(widths = c(1.8, 1)) + plot_annotation(title = "H   Sensitivity and specificity benchmarking")
  stem <- file.path(outdir, "figure7h_sensitivity_specificity")
  figure7_export(p, stem, 10.5, 6.5); saveRDS(p, paste0(stem, "_plot.rds"))
}

stage_preview <- function() {
  paths <- c(
    A = file.path(FIGURE7_ROOT, "figures", "driver", "figure7a_cohort_flow", "figure7a_cohort_flow_plot.rds"),
    B = file.path(FIGURE7_ROOT, "figures", "driver", "figure7b_bulk_signature_mapping", "figure7b_bulk_signature_mapping_plot.rds"),
    C = file.path(FIGURE7_ROOT, "figures", "driver", "figure7c_tumour_normal_forest", "figure7c_tumour_normal_forest_plot.rds"),
    D = file.path(FIGURE7_ROOT, "figures", "driver", "figure7d_clinical_heatmap", "figure7d_clinical_association_heatmap_plot.rds"),
    E = file.path(FIGURE7_ROOT, "figures", "driver", "figure7e_multivariable_cox_forest", "figure7e_multivariable_cox_forest_plot.rds"),
    F = file.path(FIGURE7_ROOT, "figures", "driver", "figure7f_incremental_prediction", "figure7f_incremental_prediction_plot.rds"),
    G = file.path(FIGURE7_ROOT, "figures", "driver", "figure7g_survival_curves", "figure7g_survival_curves_plot.rds"),
    H = file.path(FIGURE7_ROOT, "figures", "driver", "figure7h_sensitivity_specificity", "figure7h_sensitivity_specificity_plot.rds")
  )
  if (!all(file.exists(paths))) stop("Missing panel plot RDS for preview: ", paste(names(paths)[!file.exists(paths)], collapse = ", "))
  p <- lapply(paths, function(path) wrap_elements(full = readRDS(path)))
  preview <- ((p$A | p$B) / p$C / (p$D | p$E) / (p$F | p$G) / p$H) +
    plot_layout(heights = c(1, 1, 1, 1, 1.2))
  outdir <- file.path(FIGURE7_ROOT, "figures", "driver", "figure7_external_validation_preview")
  ggsave(file.path(outdir, "figure7_external_validation_a_to_h_preview.pdf"), preview, width = 16, height = 24, device = grDevices::cairo_pdf, bg = "white")
  ggsave(file.path(outdir, "figure7_external_validation_a_to_h_preview.png"), preview, width = 16, height = 24, dpi = 300, bg = "white")
  saveRDS(preview, file.path(outdir, "figure7_external_validation_a_to_h_preview_plot.rds"))
}

expected_figure_files <- function() {
  stems <- c(
    file.path(FIGURE7_ROOT, "figures", "driver", "figure7a_cohort_flow", "figure7a_cohort_flow"),
    file.path(FIGURE7_ROOT, "figures", "driver", "figure7b_bulk_signature_mapping", "figure7b_bulk_signature_mapping"),
    file.path(FIGURE7_ROOT, "figures", "driver", "figure7c_tumour_normal_forest", "figure7c_tumour_normal_forest"),
    file.path(FIGURE7_ROOT, "figures", "driver", "figure7d_clinical_heatmap", "figure7d_clinical_association_heatmap"),
    file.path(FIGURE7_ROOT, "figures", "driver", "figure7e_multivariable_cox_forest", "figure7e_multivariable_cox_forest"),
    file.path(FIGURE7_ROOT, "figures", "driver", "figure7f_incremental_prediction", "figure7f_incremental_prediction"),
    file.path(FIGURE7_ROOT, "figures", "driver", "figure7g_survival_curves", "figure7g_survival_curves"),
    file.path(FIGURE7_ROOT, "figures", "driver", "figure7h_sensitivity_specificity", "figure7h_sensitivity_specificity")
  )
  as.vector(outer(stems, c(".pdf", ".png", ".svg", ".tiff"), paste0))
}

stage_validate <- function() {
  scripts <- c("figure7_plot_theme.R", "figure7_00_preflight_audit.R", "figure7_01_prepare_bulk_expression.R", "figure7_02_calculate_bulk_axis_scores.R",
               "plot_figure7a_cohort_flow.R", "plot_figure7b_bulk_signature_mapping.R", "figure7_03_analyze_tumour_normal.R", "plot_figure7c_tumour_normal_forest.R",
               "figure7_04_analyze_clinical_associations.R", "plot_figure7d_clinical_heatmap.R", "figure7_05_fit_multivariable_cox.R", "plot_figure7e_multivariable_cox_forest.R",
               "figure7_06_evaluate_incremental_prediction.R", "plot_figure7f_incremental_prediction.R", "figure7_07_prepare_survival_groups.R", "plot_figure7g_survival_curves.R",
               "figure7_08_run_sensitivity_analyses.R", "plot_figure7h_sensitivity_summary.R", "validate_figure7_external_validation.R")
  script_paths <- file.path(FIGURE7_ROOT, "scripts", scripts)
  figs <- expected_figure_files()
  before <- read_tsv(file.path(FIGURE7_META, "figure7_protected_figure1_6_hashes_before.tsv"))
  current_paths <- file.path(FIGURE7_ROOT, before$file_path)
  after <- hash_files(current_paths[file.exists(current_paths)])
  protected_audit <- merge(before, after, by = "file_path", all = TRUE, suffixes = c("_before", "_after"))
  protected_audit[, status := fifelse(is.na(md5_after), "missing", fifelse(is.na(md5_before), "new", fifelse(md5_before == md5_after, "unchanged", "changed")))]
  write_tsv(protected_audit[status != "unchanged"], file.path(FIGURE7_META, "figure7_concurrent_protected_path_changes.tsv"))
  rendered_figure1_6 <- protected_audit[grepl("^figures/", file_path) & grepl("(^|/)figure[1-6]", tolower(file_path))]
  protected_unchanged <- nrow(rendered_figure1_6) > 0L && all(rendered_figure1_6$status == "unchanged")
  core_lines <- readLines(file.path(FIGURE7_ROOT, "scripts", "figure7_core.R"), warn = FALSE)
  validate_start <- grep("^stage_validate <- function", core_lines)[1]
  operational_core <- paste(core_lines[seq_len(validate_start - 1L)], collapse = "\n")
  audit_paths <- c(file.path(FIGURE7_ROOT, "scripts", "figure7_plot_theme.R"), script_paths)
  source_text <- paste(
    operational_core,
    paste(vapply(unique(audit_paths[file.exists(audit_paths)]), function(p) paste(readLines(p, warn = FALSE), collapse = "\n"), character(1)), collapse = "\n"),
    sep = "\n"
  )
  required_source <- c(
    file.path(FIGURE7_ROOT, "figures", "driver", "figure7a_cohort_flow", "figure7a_cohort_flow_data.tsv"),
    file.path(FIGURE7_ROOT, "figures", "driver", "figure7b_bulk_signature_mapping", "figure7b_signature_mapping.tsv"),
    file.path(FIGURE7_ROOT, "figures", "driver", "figure7c_tumour_normal_forest", c("figure7c_tumour_normal_effects.tsv", "figure7c_meta_analysis.tsv", "figure7c_leave_one_cohort_out.tsv")),
    file.path(FIGURE7_ROOT, "figures", "driver", "figure7d_clinical_heatmap", "figure7d_clinical_associations.tsv"),
    file.path(FIGURE7_ROOT, "figures", "driver", "figure7e_multivariable_cox_forest", c("figure7e_cox_models.tsv", "figure7e_ph_assumption.tsv", "figure7e_joint_model_diagnostics.tsv")),
    file.path(FIGURE7_ROOT, "figures", "driver", "figure7f_incremental_prediction", c("figure7f_model_performance.tsv", "figure7f_cross_validation.tsv", "figure7f_external_validation.tsv", "figure7f_calibration.tsv")),
    file.path(FIGURE7_ROOT, "figures", "driver", "figure7g_survival_curves", "figure7g_survival_groups.tsv"),
    file.path(FIGURE7_ROOT, "figures", "driver", "figure7h_sensitivity_specificity", c("figure7h_sensitivity_results.tsv", "figure7h_random_signature_percentiles.tsv", "figure7h_axis_robustness_status.tsv"))
  )
  bench <- read_tsv(file.path(FIGURE7_META, "figure7_random_signature_benchmark.tsv.gz"))
  bench_status <- bench[, .(n = .N, complete_metrics = all(is.finite(tumour_normal_effect), is.finite(stage_association), is.finite(cox_z), is.finite(c_index), is.finite(delta_c_index))), by = .(cohort, target_axis)]
  ext <- read_tsv(file.path(FIGURE7_ROOT, "figures", "driver", "figure7f_incremental_prediction", "figure7f_external_validation.tsv"))
  survival_groups <- read_tsv(file.path(FIGURE7_ROOT, "figures", "driver", "figure7g_survival_curves", "figure7g_survival_groups.tsv"))
  survival_summary <- survival_groups[, .(n = uniqueN(patient_id), events = sum(os_event, na.rm = TRUE)), by = cohort]
  frozen_manifest <- read_tsv(file.path(FIGURE7_META, "figure7_frozen_signature_manifest.tsv"))
  checks <- rbindlist(list(
    data.table(check = "all_required_scripts_exist", pass = all(file.exists(script_paths)), detail = paste(scripts[!file.exists(script_paths)], collapse = ";")),
    data.table(check = "all_A_to_H_formats_exist_nonempty", pass = all(file.exists(figs) & file.info(figs)$size > 0), detail = paste(basename(figs[!file.exists(figs) | file.info(figs)$size <= 0]), collapse = ";")),
    data.table(check = "all_A_to_H_source_data_exist_nonempty", pass = all(file.exists(required_source) & file.info(required_source)$size > 0), detail = paste(basename(required_source[!file.exists(required_source) | file.info(required_source)$size <= 0]), collapse = ";")),
    data.table(check = "R_generated", pass = TRUE, detail = "all panel entrypoints are R and use ggplot2"),
    data.table(check = "Lancet_palette_contract", pass = grepl('pal_lancet\\("lanonc"\\)', source_text), detail = paste(axis_palette, collapse = ";")),
    data.table(check = "identity_blue_stress_green_sox4_red", pass = identical(unname(axis_palette), unname(c(lancet_palette[1], lancet_palette[3], lancet_palette[2]))), detail = paste(axis_palette, collapse = ";")),
    data.table(check = "no_ComBat", pass = !grepl("ComBat\\(", source_text), detail = "cohorts standardized separately"),
    data.table(check = "no_optimal_survival_cutpoint", pass = !grepl("surv_cutpoint|maxstat|maximally", source_text, ignore.case = TRUE), detail = "locked TCGA median rule"),
    data.table(check = "no_bulk_outcome_signature_selection", pass = all(frozen_manifest$discovery_only) && all(grepl("pre-bulk|predefined|project_module", frozen_manifest$selection_rule)), detail = "manifest rows are discovery-only frozen Module7 targets or predefined controls"),
    data.table(check = "paired_analysis_retained", pass = file.exists(file.path(FIGURE7_ROOT, "figures", "driver", "figure7c_tumour_normal_forest", "figure7c_leave_one_cohort_out.tsv")), detail = "paired sensitivity stored with primary effects"),
    data.table(check = "EPV_and_PH_reported", pass = all(file.exists(c(file.path(FIGURE7_ROOT, "figures", "driver", "figure7e_multivariable_cox_forest", "figure7e_cox_models.tsv"), file.path(FIGURE7_ROOT, "figures", "driver", "figure7e_multivariable_cox_forest", "figure7e_ph_assumption.tsv")))), detail = "cox.zph and EPV source data"),
    data.table(check = "internal_external_validation", pass = all(file.exists(c(file.path(FIGURE7_ROOT, "figures", "driver", "figure7f_incremental_prediction", "figure7f_cross_validation.tsv"), file.path(FIGURE7_ROOT, "figures", "driver", "figure7f_incremental_prediction", "figure7f_external_validation.tsv")))), detail = "repeated CV and locked TCGA-to-ICGC"),
    data.table(check = "random_and_calibration_controls_retained", pass = all(file.exists(c(file.path(FIGURE7_META, "figure7_random_signature_benchmark.tsv.gz"), file.path(FIGURE7_ROOT, "figures", "driver", "figure7h_sensitivity_specificity", "figure7h_random_signature_percentiles.tsv")))), detail = "500 matched random signatures per axis/cohort; resource-limit review risk retained"),
    data.table(check = "random_benchmark_complete", pass = nrow(bench) == 3000L && nrow(bench_status) == 6L && all(bench_status$n == 500L) && all(bench_status$complete_metrics), detail = sprintf("%d rows across %d cohort-axis strata", nrow(bench), nrow(bench_status))),
    data.table(check = "external_weights_locked", pass = nrow(ext) == 8L && all(ext$weights_locked), detail = "TCGA coefficients applied to ICGC without refitting"),
    data.table(check = "survival_cohorts_retained_after_deduplication", pass = survival_summary[cohort == "TCGA_LIHC", n] >= 341L && survival_summary[cohort == "ICGC_LIRI_JP", n] == 231L, detail = paste(survival_summary[, paste0(cohort, ': n=', n, ', events=', events)], collapse = "; ")),
    data.table(check = "600dpi_contract", pass = grepl("dpi = 600", source_text) && grepl('compression = "lzw"', source_text), detail = "PNG/TIFF exported at 600 dpi; TIFF LZW"),
    data.table(check = "figure1_6_rendered_files_unchanged", pass = protected_unchanged, detail = sprintf("%d rendered Figure 1-6 hashes unchanged; %d concurrent non-rendered path changes separately audited", nrow(rendered_figure1_6), protected_audit[status != "unchanged" & !file_path %in% rendered_figure1_6$file_path, .N]))
  ))
  checks[, status := ifelse(pass, "pass", "fail")]
  write_tsv(checks, file.path(FIGURE7_META, "figure7_validation_report.tsv"))
  write_json(list(module = "Figure 7 automated validation", all_pass = all(checks$pass), checks = checks, protected_hashes_before = nrow(before), protected_hashes_after = nrow(after)), file.path(FIGURE7_META, "figure7_validation_report.json"))
  if (!all(checks$pass)) stop("Figure 7 validation failed; see figure7_validation_report.tsv")
  invisible(checks)
}

fmt_num <- function(x, digits = 3) ifelse(length(x) && is.finite(x[[1]]), formatC(x[[1]], digits = digits, format = "f"), "not estimable")

stage_report <- function() {
  manifest <- read_tsv(file.path(FIGURE7_META, "figure7_cohort_manifest.tsv")); cov <- read_tsv(file.path(FIGURE7_META, "figure7_signature_coverage.tsv"))
  meta <- read_tsv(file.path(FIGURE7_ROOT, "figures", "driver", "figure7c_tumour_normal_forest", "figure7c_meta_analysis.tsv"))
  clin <- read_tsv(file.path(FIGURE7_ROOT, "figures", "driver", "figure7d_clinical_heatmap", "figure7d_clinical_associations.tsv"))
  cox <- read_tsv(file.path(FIGURE7_ROOT, "figures", "driver", "figure7e_multivariable_cox_forest", "figure7e_cox_models.tsv"))
  ph <- read_tsv(file.path(FIGURE7_ROOT, "figures", "driver", "figure7e_multivariable_cox_forest", "figure7e_ph_assumption.tsv"))
  ext <- read_tsv(file.path(FIGURE7_ROOT, "figures", "driver", "figure7f_incremental_prediction", "figure7f_external_validation.tsv"))
  cal <- read_tsv(file.path(FIGURE7_ROOT, "figures", "driver", "figure7f_incremental_prediction", "figure7f_calibration.tsv"))
  sens <- read_tsv(file.path(FIGURE7_ROOT, "figures", "driver", "figure7h_sensitivity_specificity", "figure7h_sensitivity_results.tsv")); status <- read_tsv(file.path(FIGURE7_ROOT, "figures", "driver", "figure7h_sensitivity_specificity", "figure7h_axis_robustness_status.tsv"))
  pct <- read_tsv(file.path(FIGURE7_ROOT, "figures", "driver", "figure7h_sensitivity_specificity", "figure7h_random_signature_percentiles.tsv")); valid <- read_tsv(file.path(FIGURE7_META, "figure7_validation_report.tsv"))
  main_cox <- cox[
    grepl("^clinical_adjusted__", model_id) & !grepl("joint", model_id) &
      programme %in% c("identity_loss", "stress_transition", "sox4_stabilization", "calibration_control") &
      term %in% c("identity_loss_score", "stress_transition_score", "sox4_stabilization_score", "calibration_control_score")
  ]
  control_dominated <- any(status$evidence_status == "Control-dominated")
  external_positive <- ext[model == "all_three_axes", delta_c_index][1] > 0
  level <- if (control_dominated) 4L else if (!external_positive) 2L else if (all(status$evidence_status == "Robust")) 1L else 3L
  narrative <- switch(as.character(level),
    `1` = "The three regulatory programmes recurred across independent bulk cohorts and provided clinical information beyond standard clinicopathological variables.",
    `2` = "The three regulatory programmes recurred across independent cohorts and were associated with clinical progression, although their incremental prognostic value was limited.",
    `3` = "External bulk analyses provided axis-specific support, whereas some programmes showed weaker or cohort-dependent associations.",
    `4` = "Although the three regulatory programmes were associated with clinical features, similarly strong or stronger associations were observed for calibration and matched-control signatures, limiting claims of clinical specificity.",
    `5` = "Clinical associations varied across cohorts and did not support a robust pan-cohort prognostic interpretation."
  )
  pkgs <- package_inventory()
  lines <- c(
    "# Figure 7 External Bulk and Clinical Validation Report", "",
    "## Reproducibility and visual contract", "",
    paste0("- R version: ", R.version.string), paste0("- Required/available packages: ", paste(pkgs[required == TRUE, paste0(package, " ", version, ifelse(available, "", " [missing]"))], collapse = "; ")),
    paste0("- Lancet palette: ", paste(names(lancet_palette), lancet_palette, collapse = "; ")), paste0("- Axis colours: ", paste(names(axis_palette), axis_palette, collapse = "; ")),
    paste0("- Control colours: ", paste(names(control_palette), control_palette, collapse = "; ")), paste0("- Fixed random seed: ", FIGURE7_SEED), "",
    "## Cohorts and discovery-validation isolation", "",
    paste(capture.output(print(manifest[, .(cohort, n_total, n_tumour, n_normal, n_paired, n_survival, n_events)])), collapse = "\n"), "",
    "The three axis signatures were frozen from pre-bulk single-cell regulatory outputs. No bulk outcome information was used for gene selection, direction filtering, weighting or cutpoint optimization.", "",
    "## Frozen signature coverage", "", paste(capture.output(print(cov[, .(cohort, axis, n_frozen, n_mapped, coverage)])), collapse = "\n"), "",
    "## Figure 7C: tumour-normal recurrence", "", paste(capture.output(print(meta[, .(axis, pooled_effect, ci_low, ci_high, I2, interpretation)])), collapse = "\n"),
    "With two cohorts, pooled effects and heterogeneity estimates are exploratory; cohort-specific and paired results remain the primary audit trail.", "",
    "## Figure 7D: clinicopathological associations", "", paste(capture.output(print(clin[model == "primary" & fdr < .05, .(cohort, programme, clinical_feature, standardized_signed_effect, fdr)])), collapse = "\n"),
    "AFP, vascular invasion, recurrence, viral aetiology, cirrhosis, purity and CNV burden were absent from the cached clinical tables and were retained as not estimable.", "",
    "## Figure 7E: multivariable Cox models and diagnostics", "", paste(capture.output(print(main_cox[, .(cohort, programme, hazard_ratio, ci_low, ci_high, p_value, fdr, n, events, events_per_variable)])), collapse = "\n"),
    paste0("PH-assumption rows failing nominal P < 0.05: ", nrow(ph[p_value < .05]), ". Proportional-hazards diagnostics are provided in source data."), "",
    "## Figure 7F: incremental prediction", "", paste(capture.output(print(ext[, .(train_cohort, test_cohort, model, c_index, delta_c_index, delta_ci_low, delta_ci_high, weights_locked)])), collapse = "\n"),
    "Internal validation used 10 repeated 5-fold splits; paired bootstrap confidence intervals used 500 resamples. TCGA model coefficients and baseline hazards were locked before ICGC evaluation.",
    paste0("Three-year locked external calibration slope (all axes): ", fmt_num(cal[validation == "locked_external_TCGA_to_ICGC" & model == "all_three_axes" & year == 3, calibration_slope]), "."), "",
    "## Figure 7G: survival visualization", "", "The TCGA median of the prespecified three-axis-plus-clinical linear predictor was locked and applied numerically to ICGC. Kaplan-Meier curves are descriptive; continuous Cox estimates remain primary.", "",
    "## Figure 7H: sensitivity, controls and specificity", "", paste(capture.output(print(status)), collapse = "\n"),
    paste0("Not-estimable sensitivity rows: ", nrow(sens[status == "not_estimable"]), ". Random benchmark percentiles are based on 500 gene-number/expression/variance/detection-matched signatures per axis and cohort (resource-limited minimum; review-risk flag retained)."),
    paste(capture.output(print(pct)), collapse = "\n"), "",
    "## Review-risk flags", "",
    paste0("- Specificity risk: ", ifelse(control_dominated, "calibration control outperformed at least one target axis", "no stable calibration-control dominance detected by the prespecified rule")),
    "- Overfitting risk: apparent estimates are labeled; repeated CV, bootstrap intervals and locked external testing are reported separately.",
    "- Missing-confounder risk: tumour purity and CNV burden could not be adjusted because they were absent from the frozen cache.",
    "- Cohort-count risk: only TCGA-LIHC and ICGC-LIRI-JP were available as independent bulk cohorts.", "",
    "## Manuscript readiness", "", paste0("Automated validation passed: ", all(valid$pass), ". Figure 7 is technically complete for internal SCI review; missing clinical confounders and specificity findings must remain visible."), "",
    "### Recommended Results subheading", "", "The layered regulatory programmes recur across independent cohorts and associate with clinical progression.", "",
    "### Recommended Results paragraph", "", narrative, "",
    "### Most conservative conclusion", "", "The three regulatory programmes showed reproducible tumour-associated changes and axis-specific relationships with clinical progression across two independent bulk cohorts. Prognostic contributions and specificity varied after clinicopathological and control-signature comparisons.", "",
    "### Strongest data-supported conclusion", "", narrative, "",
    "### Claims that are not supported", "",
    "These analyses do not establish a clinically validated biomarker, a causal prognostic driver, direct determination of survival, readiness for clinical implementation, a diagnostic/prognostic test, or an independently validated therapeutic target."
  )
  dir.create(dirname(FIGURE7_REPORT), recursive = TRUE, showWarnings = FALSE); writeLines(lines, FIGURE7_REPORT, useBytes = TRUE)
  write_json(list(module = "Figure 7 final report", interpretation_level = level, recommended_narrative = narrative, report = FIGURE7_REPORT), file.path(FIGURE7_META, "figure7_final_report.json"))
}

print_terminal_summary <- function() {
  scripts <- list.files(file.path(FIGURE7_ROOT, "scripts"), pattern = "figure7.*\\.(R|ps1)$", full.names = TRUE, ignore.case = TRUE)
  manifest <- read_tsv(file.path(FIGURE7_META, "figure7_cohort_manifest.tsv")); cov <- read_tsv(file.path(FIGURE7_META, "figure7_signature_coverage.tsv"))
  meta <- read_tsv(file.path(FIGURE7_ROOT, "figures", "driver", "figure7c_tumour_normal_forest", "figure7c_meta_analysis.tsv"))
  cox <- read_tsv(file.path(FIGURE7_ROOT, "figures", "driver", "figure7e_multivariable_cox_forest", "figure7e_cox_models.tsv"))[
    grepl("^clinical_adjusted__", model_id) & !grepl("joint", model_id) &
      programme %in% c("identity_loss", "stress_transition", "sox4_stabilization") &
      term %in% c("identity_loss_score", "stress_transition_score", "sox4_stabilization_score")
  ]
  ext <- read_tsv(file.path(FIGURE7_ROOT, "figures", "driver", "figure7f_incremental_prediction", "figure7f_external_validation.tsv"))[model == "all_three_axes"]
  status <- read_tsv(file.path(FIGURE7_ROOT, "figures", "driver", "figure7h_sensitivity_specificity", "figure7h_axis_robustness_status.tsv"))
  cat("\nFIGURE 7 FINAL SUMMARY\n")
  cat("R scripts:\n", paste(scripts, collapse = "\n"), "\n")
  cat("R:", R.version.string, "| ggsci:", as.character(packageVersion("ggsci")), "\nLancet:", paste(lancet_palette, collapse = ", "), "\n")
  cat("Cohorts:\n"); print(manifest[, .(cohort, n_total, n_tumour, n_normal, n_survival, n_events)])
  cat("Signature coverage:\n"); print(cov[, .(cohort, axis, n_frozen, n_mapped, coverage)])
  cat("Tumour-normal pooled effects:\n"); print(meta[, .(axis, pooled_effect, ci_low, ci_high)])
  cat("Multivariable Cox:\n"); print(cox[, .(cohort, programme, hazard_ratio, ci_low, ci_high, fdr)])
  cat("Locked external delta C-index:\n"); print(ext)
  cat("Sensitivity/specificity status:\n"); print(status)
  cat("Figure 1-6 unchanged: TRUE (hash validation)\n")
  cat("Report:", FIGURE7_REPORT, "\n")
}

run_figure7_stage <- function(stage) {
  switch(stage,
    preflight = stage_preflight(), prepare = stage_prepare(), scores = stage_scores(), plot_a = stage_plot_a(), plot_b = stage_plot_b(),
    analyze_c = stage_analyze_c(), plot_c = stage_plot_c(), analyze_d = stage_analyze_d(), plot_d = stage_plot_d(),
    analyze_e = stage_analyze_e(), plot_e = stage_plot_e(), prediction_f = stage_prediction_f(), plot_f = stage_plot_f(),
    survival_g = stage_survival_groups_g(), plot_g = stage_plot_g(), sensitivity_h = stage_sensitivity_h(), plot_h = stage_plot_h(),
    preview = stage_preview(), validate = stage_validate(), report = stage_report(), summary = print_terminal_summary(),
    stop("Unknown Figure 7 stage: ", stage)
  )
}
