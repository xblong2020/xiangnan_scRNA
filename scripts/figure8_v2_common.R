#!/usr/bin/env Rscript

if (.Platform$OS.type == "windows") {
  suppressWarnings(try(Sys.setlocale("LC_CTYPE", "Chinese"), silent = TRUE))
  suppressWarnings(try(Sys.setlocale("LC_COLLATE", "Chinese"), silent = TRUE))
}

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x[[1]])) y else x

figure8_v2_project_root <- function() {
  explicit <- Sys.getenv("FIGURE8_PROJECT_ROOT", unset = "")
  if (nzchar(explicit)) return(normalizePath(explicit, winslash = "/", mustWork = TRUE))
  cwd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists("AGENTS.md") && dir.exists("scripts")) return(".")
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(file_arg)) stop("Cannot determine project root; run from the project root or set FIGURE8_PROJECT_ROOT")
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), winslash = "/", mustWork = TRUE)
}

FIGURE8_V2_ROOT <- figure8_v2_project_root()
FIGURE8_V2_METADATA <- file.path(FIGURE8_V2_ROOT, "metadata/driver/figure8_transcriptomic_reversal_v2_mainfigure")
FIGURE8_V2_DATA <- file.path(FIGURE8_V2_ROOT, "data/processed/driver/figure8_transcriptomic_reversal_v2_mainfigure")
FIGURE8_V2_FIGURES <- file.path(FIGURE8_V2_ROOT, "figures/driver/figure8_transcriptomic_reversal_v2_mainfigure")
FIGURE8_V2_REPORTS <- file.path(FIGURE8_V2_ROOT, "reports/figure8_transcriptomic_reversal_v2_mainfigure")
FIGURE8_V2_SEED <- 20260805L

figure8_v2_existing_r_library <- function() {
  candidates <- c(
    file.path(Sys.getenv("USERPROFILE"), ".codex/r-libs/figure8-v2"),
    "data/processed/driver/figure8_transcriptomic_reversal/r_library",
    "data/processed/driver/figure6_directional_network/r_library"
  )
  hit <- candidates[dir.exists(candidates)]
  if (!length(hit)) stop("No existing frozen project R library is available")
  hit[[1]]
}

figure8_v2_init_dirs <- function() {
  dirs <- c(FIGURE8_V2_METADATA, FIGURE8_V2_DATA, FIGURE8_V2_FIGURES, FIGURE8_V2_REPORTS)
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

figure8_v2_read_tsv <- function(path, ...) {
  if (!file.exists(path)) stop("Required v2 input is missing: ", path)
  fread(path, ...)
}

figure8_v2_write_tsv <- function(x, filename, directory = FIGURE8_V2_METADATA, compress = grepl("[.]gz$", filename)) {
  figure8_v2_init_dirs()
  path <- file.path(directory, filename)
  fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA", compress = if (compress) "gzip" else "none")
  path
}

figure8_v2_write_json <- function(x, filename, directory = FIGURE8_V2_METADATA) {
  figure8_v2_init_dirs()
  path <- file.path(directory, filename)
  write_json(x, path, pretty = TRUE, auto_unbox = TRUE, na = "null", digits = NA)
  path
}

figure8_v2_robust_unit <- function(x, clip = 3) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  finite <- is.finite(x)
  if (!any(finite)) return(out)
  center <- median(x[finite])
  spread <- mad(x[finite], center = center, constant = 1)
  if (!is.finite(spread) || spread <= .Machine$double.eps) {
    spread <- max(abs(x[finite] - center))
  }
  if (!is.finite(spread) || spread <= .Machine$double.eps) {
    out[finite] <- 0
    return(out)
  }
  z <- (x[finite] - center) / spread
  out[finite] <- pmax(-clip, pmin(clip, z)) / clip
  out
}

figure8_v2_rank_percentile <- function(rank, n_compounds) {
  pmax(0, pmin(1, 1 - (as.numeric(rank) - 1) / max(1, n_compounds - 1)))
}

figure8_v2_agreement_shrunk <- function(components, weights) {
  components <- as.numeric(components)
  weights <- as.numeric(weights)
  keep <- is.finite(components) & is.finite(weights) & weights > 0
  if (!any(keep)) return(list(weighted_mean = NA_real_, agreement = NA_real_, score = NA_real_, n_components = 0L))
  components <- pmax(-1, pmin(1, components[keep]))
  weights <- weights[keep]
  signed_sum <- sum(weights * components)
  weighted_mean <- signed_sum / sum(weights)
  absolute_mass <- sum(weights * abs(components))
  agreement <- if (absolute_mass <= .Machine$double.eps) 0 else abs(signed_sum) / absolute_mass
  list(
    weighted_mean = weighted_mean,
    agreement = agreement,
    score = weighted_mean * agreement,
    n_components = length(components)
  )
}

figure8_v2_safe_name <- function(x) {
  x <- iconv(as.character(x), from = "", to = "ASCII//TRANSLIT", sub = "")
  tolower(gsub("[^a-zA-Z0-9]+", "", trimws(x)))
}

figure8_v2_entity_key <- function(inchi_key, canonical_name, brd_id) {
  inchi <- trimws(as.character(inchi_key))
  name <- figure8_v2_safe_name(canonical_name)
  brd <- trimws(as.character(brd_id))
  invalid_inchi <- tolower(inchi) %in% c("-666", "na", "nan", "none", "unknown", "restricted", "not available", "unavailable")
  has_inchi <- !is.na(inchi) & nzchar(inchi) & !invalid_inchi
  has_name <- !is.na(name) & nzchar(name)
  ifelse(has_inchi, paste0("INCHI:", inchi), ifelse(has_name, paste0("NAME:", name), paste0("BRD:", brd)))
}

figure8_v2_empirical_p <- function(observed, null) {
  null <- as.numeric(null)
  null <- null[is.finite(null)]
  if (!is.finite(observed) || !length(null)) return(list(p_upper = NA_real_, p_lower = NA_real_, p_two_sided = NA_real_, n_null = length(null)))
  p_upper <- (1 + sum(null >= observed)) / (1 + length(null))
  p_lower <- (1 + sum(null <= observed)) / (1 + length(null))
  list(p_upper = p_upper, p_lower = p_lower, p_two_sided = min(1, 2 * min(p_upper, p_lower)), n_null = length(null))
}

figure8_v2_specificity_label <- function(p) {
  if (!is.finite(p)) return("unavailable")
  if (p < 0.05) "strong" else if (p <= 0.10) "suggestive" else "not_specific"
}

figure8_v2_directional_specificity_label <- function(p_upper, p_lower, p_two_sided) {
  if (!all(is.finite(c(p_upper, p_lower, p_two_sided)))) return("unavailable")
  if (p_two_sided < 0.05 && p_upper < 0.05) return("strong")
  if (p_two_sided <= 0.10 && p_upper <= 0.10) return("suggestive")
  if (p_two_sided < 0.05 && p_lower < 0.05) return("significantly_worse")
  "not_specific"
}

figure8_v2_evidence_summary <- function(values, weights) {
  values <- as.numeric(values)
  weights <- as.numeric(weights)
  if (length(values) != length(weights)) stop("values and weights must have equal length")
  valid_weight <- is.finite(weights) & weights > 0
  total_weight <- sum(weights[valid_weight])
  observed <- valid_weight & is.finite(values)
  observed_weight <- sum(weights[observed])
  conservative <- if (total_weight > 0) sum(values[observed] * weights[observed]) / total_weight else NA_real_
  coverage_aware <- if (observed_weight > 0) sum(values[observed] * weights[observed]) / observed_weight else NA_real_
  list(
    conservative_score = conservative,
    coverage_aware_score = coverage_aware,
    coverage_confidence = if (total_weight > 0) observed_weight / total_weight else NA_real_,
    missing_dimensions = which(!is.finite(values) & valid_weight)
  )
}

figure8_v2_prism_class <- function(liver_percentile, pan_percentile, liver_minus_pan, n_liver_lines, enrichment_delta = 0.10) {
  if (!is.finite(n_liver_lines) || n_liver_lines < 3) return("insufficient_liver_lines")
  if (!all(is.finite(c(liver_percentile, pan_percentile, liver_minus_pan)))) return("unavailable")
  if (liver_percentile >= 0.80 && pan_percentile >= 0.80 && liver_minus_pan < enrichment_delta) return("broad_cytotoxicity")
  if (liver_percentile >= 0.80 && liver_minus_pan >= enrichment_delta) return("hcc_liver_enriched")
  if (pan_percentile >= 0.80) return("pan_cancer_activity")
  "no_enriched_support"
}

figure8_v2_prism_score <- function(liver_percentile, pan_percentile, liver_minus_pan) {
  if (!all(is.finite(c(liver_percentile, pan_percentile, liver_minus_pan)))) return(NA_real_)
  enrichment <- pmax(0, pmin(1, (liver_minus_pan + 0.25) / 0.50))
  pmax(0, pmin(1, 0.70 * liver_percentile + 0.30 * enrichment - 0.25 * pan_percentile))
}

figure8_v2_assign_tier <- function(
  rank_stability, fold_agreement, specificity_p, cmap_support, cmap_profiled,
  network_moa_score, nuisance_penalty, prism_class, coverage, strong_opposition
) {
  failed <- character()
  if (isTRUE(strong_opposition)) failed <- c(failed, "strong_opposition")
  if (is.finite(nuisance_penalty) && nuisance_penalty >= 0.75) failed <- c(failed, "nuisance_penalty")
  if (!is.na(prism_class) && prism_class == "broad_cytotoxicity") failed <- c(failed, "broad_cytotoxicity")
  if (length(failed)) return(list(tier = "discordant", failed_gates = paste(failed, collapse = ",")))

  stable <- is.finite(rank_stability) && rank_stability >= 0.75
  folds <- is.finite(fold_agreement) && fold_agreement >= 0.75
  specific <- is.finite(specificity_p) && specificity_p <= 0.10
  cmap <- isTRUE(cmap_support)
  mechanism <- is.finite(network_moa_score) && network_moa_score >= 0.50
  phenotype <- !is.na(prism_class) && prism_class == "hcc_liver_enriched"

  if (stable && folds && specific && cmap && mechanism && phenotype && is.finite(coverage) && coverage >= 0.75) {
    return(list(tier = "tier_A", failed_gates = ""))
  }
  if (stable && folds && specific && cmap && mechanism && is.finite(coverage) && coverage >= 0.60) {
    return(list(tier = "tier_B", failed_gates = if (phenotype) "" else "orthogonal_phenotype"))
  }
  if (stable && folds && is.finite(coverage) && coverage >= 0.50) {
    return(list(tier = "tier_C", failed_gates = paste(c(if (!specific) "specificity", if (!cmap && isTRUE(cmap_profiled)) "cmap_corroboration", if (!mechanism) "mechanism"), collapse = ",")))
  }
  list(tier = "unresolved", failed_gates = "insufficient_internal_support_or_coverage")
}
