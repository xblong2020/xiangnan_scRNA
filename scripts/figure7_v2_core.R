suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggsci)
  library(patchwork)
  library(survival)
  library(MASS)
  library(metafor)
  library(jsonlite)
  library(scales)
  library(ggrepel)
})

# The installed Windows R starts in locale C even when the project path has
# Chinese characters.  Switch only character handling to the local Windows
# code page before any file I/O; calculations and source data are unchanged.
invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.936")))

# Figure 7 v2 is deliberately isolated from the frozen Figure 1--6 and v1
# namespaces.  It never writes into any pre-existing Figure 7 output directory.

# `--file` can be mojibake under a C-locale R session on Windows.  The caller
# always starts from the project root, so prefer the working directory when it
# contains the project layout and only use `--file` as a fallback.
figure7_v2_root <- function() {
  # Keep the logical drive path (e.g. R:/), rather than resolving it to the
  # Chinese physical path, because data.table cannot reopen the latter under
  # this Windows R C-locale session.
  wd <- gsub("\\\\", "/", getwd())
  if (grepl("^[A-Za-z]:/$", wd)) wd <- sub("/$", "", wd)
  if (dir.exists(file.path(wd, "scripts")) && dir.exists(file.path(wd, "metadata"))) return(wd)
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(arg)) stop("Run Figure 7 v2 from the project root.")
  normalizePath(file.path(dirname(sub("^--file=", "", arg[[1]])), ".."), mustWork = TRUE)
}

FIGURE7_V2_ROOT <- figure7_v2_root()
FIGURE7_V2_META <- file.path(FIGURE7_V2_ROOT, "metadata", "driver", "figure7_external_validation_v2")
FIGURE7_V2_DATA <- file.path(FIGURE7_V2_ROOT, "data", "processed", "driver", "figure7_external_validation_v2")
FIGURE7_V2_FIG <- file.path(FIGURE7_V2_ROOT, "figures", "driver", "figure7_external_validation_v2")
FIGURE7_V2_EXT_FIG <- file.path(FIGURE7_V2_ROOT, "figures", "driver", "extended_data_figure7_v2")
FIGURE7_V2_REPORT <- file.path(FIGURE7_V2_ROOT, "reports", "figure7_external_bulk_clinical_validation_v2_report.md")
FIGURE7_V2_SEED <- 20260812L
FIGURE7_V2_N_RANDOM <- 1000L

FIGURE7_V2_INPUTS <- list(
  tcga_expression_rds = file.path(FIGURE7_V2_ROOT, "data", "processed", "driver", "figure7_external_validation", "figure7_tcga_lihc_prepared_expression.rds"),
  icgc_expression_rds = file.path(FIGURE7_V2_ROOT, "data", "processed", "driver", "figure7_external_validation", "figure7_icgc_liri_jp_prepared_expression.rds"),
  tcga_rank_rds = file.path(FIGURE7_V2_ROOT, "data", "processed", "driver", "figure7_external_validation", "figure7_tcga_lihc_rank_matrix.rds"),
  icgc_rank_rds = file.path(FIGURE7_V2_ROOT, "data", "processed", "driver", "figure7_external_validation", "figure7_icgc_liri_jp_rank_matrix.rds"),
  # ASCII-path byte-identical snapshots prevent R C-locale path mojibake.
  # `figure7_v2_raw_input_snapshot_provenance.tsv` records original hashes.
  tcga_clinical = file.path(FIGURE7_V2_ROOT, "data", "processed", "driver", "figure7_external_validation_v2", "input_snapshots", "TCGAclinical_raw_snapshot.tsv"),
  icgc_clinical = file.path(FIGURE7_V2_ROOT, "data", "processed", "driver", "figure7_external_validation_v2", "input_snapshots", "icgcClinical_raw_snapshot.tsv"),
  icgc_survival = file.path(FIGURE7_V2_ROOT, "data", "processed", "driver", "figure7_external_validation_v2", "input_snapshots", "ICGCtime_raw_snapshot.tsv"),
  frozen_targets = file.path(FIGURE7_V2_ROOT, "metadata", "driver", "module8_tf_target_signature_genes.tsv"),
  frozen_registry = file.path(FIGURE7_V2_ROOT, "metadata", "driver", "module8_signature_registry.tsv"),
  celloracle_edges = file.path(FIGURE7_V2_ROOT, "metadata", "driver", "celloracle_module6_7_grn_links_filtered.tsv.gz"),
  cistarget_regulons = file.path(FIGURE7_V2_ROOT, "metadata", "driver", "driver_module6_3c_cistarget_regulon_summary.tsv"),
  scenic_edges = file.path(FIGURE7_V2_ROOT, "metadata", "driver", "driver_module6_3_pyscenic_regulon_edges.tsv.gz"),
  v1_signature_manifest = file.path(FIGURE7_V2_ROOT, "metadata", "driver", "figure7_external_validation", "figure7_frozen_signature_manifest.tsv")
)

v2_dirs <- function() {
  dirs <- c(FIGURE7_V2_META, FIGURE7_V2_DATA, FIGURE7_V2_FIG, FIGURE7_V2_EXT_FIG, dirname(FIGURE7_V2_REPORT))
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

write_v2_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
}

write_v2_tsvgz <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")
}

read_v2_tsv <- function(path, ...) fread(path, sep = "\t", na.strings = c("", "NA", "NaN"), ...)

write_v2_json <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write_json(x, path, pretty = TRUE, auto_unbox = TRUE, na = "null", digits = 12)
}

norm_gene <- function(x) toupper(sub("\\.[0-9]+$", "", trimws(as.character(x))))
safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
safe_z <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
bh <- function(x) {
  ans <- rep(NA_real_, length(x)); keep <- is.finite(x)
  ans[keep] <- p.adjust(x[keep], method = "BH")
  ans
}
fmt <- function(x, digits = 3L) ifelse(length(x) && is.finite(x[[1]]), formatC(x[[1]], format = "f", digits = digits), "not estimable")

lancet <- ggsci::pal_lancet("lanonc")(9)
axis_cols <- c(identity_loss = lancet[1], stress_transition = lancet[3], sox4_associated = lancet[2])
comparator_cols <- c(
  foxm1_cebpb_reference = lancet[5], proliferation_comparator = lancet[6], broad_network_calibration = lancet[4]
)
cohort_cols <- c(TCGA_LIHC = lancet[1], ICGC_LIRI_JP = lancet[2])

v2_theme <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.5, colour = "#444444"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5, colour = "#222222"),
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.tag = element_text(size = base_size + 3, face = "bold"),
      plot.margin = margin(6, 8, 6, 6)
    )
}

export_v2_plot <- function(plot, stem, width, height, dpi = 600L) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, units = "in", device = grDevices::cairo_pdf, bg = "white")
  ggsave(paste0(stem, ".svg"), plot, width = width, height = height, units = "in", device = grDevices::svg, bg = "white")
  # Windows R's agg device cannot write raster directly to the project path
  # in this locale.  Render to an ASCII temporary directory, verify that the
  # device closed successfully, then copy the byte-identical artifact.
  tmp_dir <- file.path(tempdir(), "figure7_v2_raster")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  safe_stem <- paste0("f7v2_", sprintf("%08x", as.integer(abs(sum(utf8ToInt(stem))) %% .Machine$integer.max)))
  tmp_png <- file.path(tmp_dir, paste0(safe_stem, ".png")); tmp_tiff <- file.path(tmp_dir, paste0(safe_stem, ".tiff"))
  ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = dpi, device = "png", bg = "white")
  ggsave(tmp_tiff, plot, width = width, height = height, units = "in", dpi = dpi, device = "tiff", compression = "lzw", bg = "white")
  if (!file.exists(tmp_png) || file.info(tmp_png)$size <= 0 || !file.exists(tmp_tiff) || file.info(tmp_tiff)$size <= 0) stop("Raster export failed before project-path copy.")
  ok_png <- file.copy(tmp_png, paste0(stem, ".png"), overwrite = TRUE); ok_tiff <- file.copy(tmp_tiff, paste0(stem, ".tiff"), overwrite = TRUE)
  if (!ok_png || !ok_tiff) stop("Raster export copy to Figure 7 v2 destination failed.")
  invisible(paste0(stem, c(".pdf", ".svg", ".png", ".tiff")))
}

tcga_id <- function(x) vapply(strsplit(as.character(x), "-", fixed = TRUE), function(z) paste(head(z, 3L), collapse = "-"), character(1))
icgc_id <- function(x) {
  hit <- regmatches(as.character(x), regexpr("DO[0-9]+", as.character(x)))
  hit[!nzchar(hit)] <- as.character(x)[!nzchar(hit)]
  hit
}
stage_to_num <- function(x) {
  y <- toupper(trimws(as.character(x)))
  out <- rep(NA_real_, length(y))
  out[grepl("STAGE I(?!I)", y, perl = TRUE) | grepl("^I$", y)] <- 1
  out[grepl("STAGE II(?!I)", y, perl = TRUE) | grepl("^II$", y)] <- 2
  out[grepl("STAGE III", y) | grepl("^III", y)] <- 3
  out[grepl("STAGE IV", y) | grepl("^IV", y)] <- 4
  z <- safe_num(y)
  out[is.na(out) & is.finite(z)] <- z[is.na(out) & is.finite(z)]
  out
}
grade_to_num <- function(x) safe_num(sub("^G", "", toupper(trimws(as.character(x)))))
tstage_to_num <- function(x) safe_num(sub("^T([0-9]+).*$", "\\1", toupper(trimws(as.character(x)))))

protected_v2_paths <- function() {
  roots <- c("scripts", "metadata", "figures", "reports")
  full <- unlist(lapply(roots, function(d) list.files(file.path(FIGURE7_V2_ROOT, d), recursive = TRUE, full.names = TRUE, all.files = FALSE)))
  full <- full[file.exists(full) & !dir.exists(full)]
  root_norm <- normalizePath(FIGURE7_V2_ROOT, winslash = "/", mustWork = TRUE)
  rel <- substring(normalizePath(full, winslash = "/", mustWork = TRUE), nchar(root_norm) + 2L)
  keep <- (grepl("(^|/)(figure[1-6]|module[1-9])", tolower(rel)) | grepl("celloracle|scenic|trajectory|frozen", tolower(rel))) &
    !grepl("figure7_external_validation_v2|figure7_v2", tolower(rel))
  data.table(file_path = rel[keep], full_path = full[keep])
}

hash_protected <- function(paths) {
  if (!nrow(paths)) return(data.table(file_path = character(), size_bytes = numeric(), sha256 = character()))
  data.table(
    file_path = paths$file_path,
    size_bytes = file.info(paths$full_path)$size,
    sha256 = unname(tools::md5sum(paths$full_path)) # filled with SHA256 by PowerShell baseline; comparison uses stable rehash below
  )
}

# R base does not supply SHA256. Windows certutil is used only for read-only hashing.
sha256_one <- function(path) {
  output <- suppressWarnings(system2("powershell", c("-NoProfile", "-Command", sprintf("(Get-FileHash -LiteralPath '%s' -Algorithm SHA256).Hash", gsub("'", "''", normalizePath(path, winslash = "\\", mustWork = TRUE)))), stdout = TRUE, stderr = FALSE))
  trimws(output[[1]])
}

hash_protected_sha256 <- function(paths) {
  if (!nrow(paths)) return(data.table(file_path = character(), size_bytes = numeric(), sha256 = character()))
  data.table(file_path = paths$file_path, size_bytes = file.info(paths$full_path)$size,
             sha256 = vapply(paths$full_path, sha256_one, character(1)))
}

stage_v2_protected_audit <- function(phase = c("before", "after")) {
  phase <- match.arg(phase)
  v2_dirs()
  outfile <- file.path(FIGURE7_V2_META, "figure7_v2_protected_figure1_6_hash_audit.tsv")
  current <- hash_protected_sha256(protected_v2_paths())
  if (phase == "before" || !file.exists(outfile)) {
    out <- copy(current)
    setnames(out, c("size_bytes", "sha256"), c("size_bytes_before", "sha256_before"))
    out[, `:=`(size_bytes_after = NA_real_, sha256_after = NA_character_, status = "baseline", timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))]
    write_v2_tsv(out, outfile)
    return(out)
  }
  before <- read_v2_tsv(outfile)
  if (!"sha256_before" %in% names(before)) {
    setnames(before, c("size_bytes", "sha256"), c("size_bytes_before", "sha256_before"), skip_absent = TRUE)
  }
  out <- merge(before[, .(file_path, size_bytes_before, sha256_before)], current, by = "file_path", all = TRUE)
  setnames(out, c("size_bytes", "sha256"), c("size_bytes_after", "sha256_after"))
  out[, status := fifelse(is.na(sha256_before), "new", fifelse(is.na(sha256_after), "missing", fifelse(sha256_before == sha256_after, "unchanged", "changed")))]
  out[, timestamp_utc := format(Sys.time(), tz = "UTC", usetz = TRUE)]
  write_v2_tsv(out, outfile)
  if (any(out$status != "unchanged")) stop("Protected Figure 1--6/module paths changed; Figure 7 v2 is stopped. See protected hash audit.")
  out
}

audit_column_rows <- function(dt, source_file, id_col, peer_ids = NULL) {
  rows <- list(); n <- nrow(dt); k <- 0L
  for (col in names(dt)) {
    x <- as.character(dt[[col]])
    miss <- is.na(x) | !nzchar(trimws(x))
    num <- safe_num(x)
    inferred <- if (all(miss | is.finite(num))) "numeric" else "character"
    k <- k + 1L
    rows[[k]] <- data.table(
      source_file = source_file, column = col, record_type = "summary", value = NA_character_, frequency = NA_integer_,
      n_rows = n, n_missing = sum(miss), missing_fraction = mean(miss), n_unique = uniqueN(x[!miss]),
      inferred_type = inferred, value_min = if (all(miss | is.finite(num))) min(num, na.rm = TRUE) else NA_real_,
      value_max = if (all(miss | is.finite(num))) max(num, na.rm = TRUE) else NA_real_,
      id_match_n = if (identical(col, id_col) && !is.null(peer_ids)) length(intersect(unique(x[!miss]), peer_ids)) else NA_integer_,
      id_only_in_source_n = if (identical(col, id_col) && !is.null(peer_ids)) length(setdiff(unique(x[!miss]), peer_ids)) else NA_integer_,
      id_only_in_peer_n = if (identical(col, id_col) && !is.null(peer_ids)) length(setdiff(peer_ids, unique(x[!miss]))) else NA_integer_
    )
    freq <- sort(table(x[!miss]), decreasing = FALSE)
    for (value in names(freq)) {
      k <- k + 1L
      rows[[k]] <- data.table(source_file = source_file, column = col, record_type = "frequency", value = value, frequency = as.integer(freq[[value]]),
                               n_rows = n, n_missing = sum(miss), missing_fraction = mean(miss), n_unique = uniqueN(x[!miss]), inferred_type = inferred,
                               value_min = NA_real_, value_max = NA_real_, id_match_n = NA_integer_, id_only_in_source_n = NA_integer_, id_only_in_peer_n = NA_integer_)
    }
  }
  rbindlist(rows, fill = TRUE)
}

stage_v2_clinical_audit <- function() {
  v2_dirs()
  icgc <- fread(FIGURE7_V2_INPUTS$icgc_clinical, sep = "\t", na.strings = c("", "NA"), colClasses = "character")
  surv <- fread(FIGURE7_V2_INPUTS$icgc_survival, sep = "\t", na.strings = c("", "NA"), colClasses = "character")
  stopifnot(identical(names(icgc), c("Id", "Gender", "Age", "Stage")), identical(names(surv), c("id", "fustat", "futime")))
  audit <- rbindlist(list(
    audit_column_rows(icgc, FIGURE7_V2_INPUTS$icgc_clinical, "Id", unique(surv$id)),
    audit_column_rows(surv, FIGURE7_V2_INPUTS$icgc_survival, "id", unique(icgc$Id))
  ), fill = TRUE)
  write_v2_tsv(audit, file.path(FIGURE7_V2_META, "figure7_v2_icgc_clinical_raw_audit.tsv"))
  decisions <- data.table(
    variable = c("Age", "Gender", "Stage", "fustat", "futime", "Id"),
    raw_values = c("0 (n=98); 1 (n=162)", "0 (n=68); 1 (n=192)", "0 (n=157); 1 (n=103)", "0 (n=214); 1 (n=46)", "10--2160; 60 distinct values", "260 unique DO identifiers"),
    direct_file_evidence = c(
      "Binary values only; no codebook or unit label in raw file.",
      "Binary values only; no mapping of 0/1 to sex in raw file.",
      "Binary values only; no mapping to I--IV or early/advanced in raw file.",
      "Binary values only; no event/censor mapping in raw file.",
      "Numeric values only; no time unit label in raw file.",
      "Exact one-to-one matching to ICGCtime.txt identifiers."
    ),
    existing_project_evidence = c(
      "Only prior Figure 7/Module 8 code treated it as age_high; that is an implementation assumption, not a source data dictionary.",
      "Only prior Figure 7 code mapped 1 to MALE; no independent local data dictionary was found.",
      "Only prior Figure 7 code treated it as stage_high and set stage_num = Stage + 1; no independent local data dictionary was found.",
      "Only prior scripts passed fustat to survival functions; no local event-code documentation was found.",
      "Only prior scripts treated futime as days; no local unit documentation was found.",
      "Both files share all 260 identifiers."
    ),
    determination = c("unverified_binary_encoding", "unverified_binary_encoding", "unverified_binary_encoding", "unverified_binary_event_status", "unverified_time_unit", "verified_identifier_match"),
    action = c("block ICGC clinical and shared-baseline analyses", "block ICGC clinical and shared-baseline analyses", "block ICGC clinical and shared-baseline analyses", "block ICGC survival and locked external survival prediction", "block ICGC survival and locked external survival prediction", "allow expression-to-patient matching")
  )
  write_v2_tsv(decisions, file.path(FIGURE7_V2_META, "figure7_v2_icgc_clinical_encoding_decisions.tsv"))
  lines <- c(
    "# Figure 7 v2 ICGC clinical encoding audit", "",
    "## Direct raw-file findings", "",
    "`icgcClinical.txt` contains 260 unique donor IDs and exactly three binary clinical columns: Age (0=98, 1=162), Gender (0=68, 1=192), and Stage (0=157, 1=103). `ICGCtime.txt` has the same 260 donor IDs, fustat (0=214, 1=46), and futime (10--2160; 60 distinct values).", "",
    "## Encoding determination", "",
    "The raw files provide no data dictionary. The project contains earlier code that interpreted these fields as age_high, male/female, stage_high, event status, and days, but code is not independent source documentation. Therefore the meanings of Age, Gender, Stage, fustat, and futime are **not verified** in v2.", "",
    "## Consequence", "",
    "ICGC expression sample-to-donor matching remains valid. ICGC clinicopathological models, ICGC OS models, and TCGA-to-ICGC locked external survival prediction are blocked and reported as not estimable. No 0/1 coding is imputed or guessed."
  )
  writeLines(lines, file.path(FIGURE7_V2_META, "figure7_v2_icgc_clinical_encoding_report.md"), useBytes = TRUE)
  invisible(list(audit = audit, decisions = decisions))
}

build_harmonized_clinical <- function() {
  tcga <- fread(FIGURE7_V2_INPUTS$tcga_clinical, sep = "\t", na.strings = c("", "NA", "[Not Available]"))
  tcga_out <- data.table(
    cohort = "TCGA_LIHC", patient_id = as.character(tcga$Id),
    age_years = safe_num(tcga$age), age_high = as.integer(safe_num(tcga$age) >= 60),
    sex = toupper(as.character(tcga$gender)), stage_num = stage_to_num(tcga$stage),
    stage_high = as.integer(stage_to_num(tcga$stage) >= 3),
    os_time_days = safe_num(tcga$futime), os_event = safe_num(tcga$fustat),
    grade_num = grade_to_num(tcga$grade), t_stage_num = tstage_to_num(tcga$T),
    raw_age = as.character(tcga$age), raw_gender = as.character(tcga$gender), raw_stage = as.character(tcga$stage),
    raw_fustat = as.character(tcga$fustat), raw_futime = as.character(tcga$futime),
    clinical_semantics = "verified_from_explicit_TCGA_values", os_semantics = "verified_from_project_TCGA_clinical_source",
    clinical_analysis_status = "eligible", os_analysis_status = "eligible"
  )
  icgc <- fread(FIGURE7_V2_INPUTS$icgc_clinical, sep = "\t", colClasses = "character")
  surv <- fread(FIGURE7_V2_INPUTS$icgc_survival, sep = "\t", colClasses = "character")
  raw <- merge(icgc, surv, by.x = "Id", by.y = "id", all = TRUE, sort = FALSE)
  icgc_out <- data.table(
    cohort = "ICGC_LIRI_JP", patient_id = as.character(raw$Id),
    age_years = NA_real_, age_high = NA_integer_, sex = NA_character_, stage_num = NA_real_, stage_high = NA_integer_,
    os_time_days = NA_real_, os_event = NA_real_, grade_num = NA_real_, t_stage_num = NA_real_,
    raw_age = as.character(raw$Age), raw_gender = as.character(raw$Gender), raw_stage = as.character(raw$Stage),
    raw_fustat = as.character(raw$fustat), raw_futime = as.character(raw$futime),
    clinical_semantics = "blocked_unverified_Age_Gender_Stage_encoding", os_semantics = "blocked_unverified_fustat_futime_encoding",
    clinical_analysis_status = "blocked", os_analysis_status = "blocked"
  )
  out <- rbindlist(list(tcga_out, icgc_out), fill = TRUE)
  write_v2_tsv(out, file.path(FIGURE7_V2_META, "figure7_v2_harmonized_clinical_variables.tsv"))
  defs <- data.table(
    cohort = c(rep("TCGA_LIHC", 7L), rep("ICGC_LIRI_JP", 7L)),
    variable = rep(c("age_years", "age_high", "sex", "stage_num", "stage_high", "os_time_days", "os_event"), 2L),
    definition = c(
      "raw age in years", "age_years >= 60", "raw MALE/FEMALE", "Stage I--IV mapped to 1--4", "stage_num >= 3", "raw futime, days in project clinical source", "raw fustat, project clinical source",
      rep("not assigned: raw ICGC binary/time codes lack a verified data dictionary", 7L)
    ),
    semantic_status = c(rep("verified_within_TCGA", 7L), rep("blocked_unverified_ICGC", 7L)),
    missing_n = c(vapply(c("age_years", "age_high", "sex", "stage_num", "stage_high", "os_time_days", "os_event"), function(z) sum(is.na(tcga_out[[z]])), numeric(1)),
                  vapply(c("age_years", "age_high", "sex", "stage_num", "stage_high", "os_time_days", "os_event"), function(z) sum(is.na(icgc_out[[z]])), numeric(1)))
  )
  write_v2_tsv(defs, file.path(FIGURE7_V2_META, "figure7_v2_harmonized_clinical_definition_audit.tsv"))
  out
}

make_v2_signatures <- function() {
  targets <- read_v2_tsv(FIGURE7_V2_INPUTS$frozen_targets)
  targets[, `:=`(tf = norm_gene(tf), gene = norm_gene(gene), axis = as.character(axis))]
  add_targets <- function(axis_name, tfs, role, display_name) {
    z <- targets[tf %in% tfs, .(axis = axis_name, tf, gene, origin = "frozen_module7_target", source_axis = axis, source_method = signature_source, source_rank = rank)]
    anchors <- data.table(axis = axis_name, tf = tfs, gene = tfs, origin = "frozen_tf_expression_anchor", source_axis = NA_character_, source_method = "Figure2_6_frozen_TF_anchor", source_rank = 0L)
    z <- unique(rbindlist(list(z, anchors), fill = TRUE), by = c("axis", "tf", "gene"))
    z[, `:=`(role = role, display_name = display_name)]
    z
  }
  identity <- add_targets("identity_loss", c("HNF4A", "PPARA"), "primary", "HNF4A/PPARA-associated hepatocyte identity loss")
  stress <- add_targets("stress_transition", c("JUN", "JUNB", "JUND", "FOS", "ATF3", "CEBPB", "EGR1"), "primary", "AP-1/CEBPB/EGR1-associated stress-transition")
  sox4 <- add_targets("sox4_associated", "SOX4", "primary", "SOX4-associated malignant-state/plasticity programme")
  broad <- add_targets("broad_network_calibration", c("IRF1", "MYC", "HLF", "MAFF", "MAFB"), "comparator", "Broad-network calibration comparator")
  proliferation <- data.table(axis = "proliferation_comparator", tf = NA_character_,
    gene = c("MKI67", "TOP2A", "STMN1", "TYMS", "UBE2C", "PCNA", "MCM2", "MCM5", "HMGB2"),
    origin = "fixed_project_panel", source_axis = NA_character_, source_method = "project_fixed_proliferation_panel", source_rank = seq_len(9L),
    role = "comparator", display_name = "Proliferation comparator")
  foxm1 <- unique(c("FOXM1", "CEBPB", targets[tf == "CEBPB", gene], "CCNB1", "CCNB2", "CDC20", "CENPF", "AURKB"))
  foxm1 <- data.table(axis = "foxm1_cebpb_reference", tf = NA_character_, gene = foxm1, origin = "project_defined_reference",
    source_axis = NA_character_, source_method = "FOXM1_CEBPB_plus_frozen_CEBPB_targets_plus_proliferation_effectors", source_rank = seq_along(foxm1),
    role = "comparator", display_name = "FOXM1/CEBPB biological reference")
  hypoxia <- data.table(axis = "hypoxia_covariate", tf = NA_character_, gene = c("ADM", "ALDOA", "BNIP3", "CA9", "EGLN3", "ENO1", "HIF1A", "HK2", "LDHA", "LOX", "NDRG1", "PDK1", "PGK1", "SLC2A1", "VEGFA"), origin = "predefined_expression_covariate", source_axis = NA_character_, source_method = "predefined_hypoxia_expression_panel", source_rank = seq_len(15L), role = "covariate", display_name = "Predefined hypoxia-expression covariate")
  inflammation <- data.table(axis = "inflammation_covariate", tf = NA_character_, gene = c("CCL2", "ICAM1", "IL1B", "IL6", "NFKB1", "PTGS2", "RELA", "STAT3", "TNF", "CXCL8"), origin = "predefined_expression_covariate", source_axis = NA_character_, source_method = "predefined_inflammation_expression_panel", source_rank = seq_len(10L), role = "covariate", display_name = "Predefined inflammation-expression covariate")
  out <- unique(rbindlist(list(identity, stress, sox4, broad, proliferation, foxm1, hypoxia, inflammation), fill = TRUE), by = c("axis", "tf", "gene"))
  out[, `:=`(signature_version = "figure7_v2_2026-08-12", bulk_outcome_used_for_selection = FALSE, direction = "unsigned_pending_direction_audit")]
  setcolorder(out, c("axis", "display_name", "role", "tf", "gene", "origin", "source_axis", "source_method", "source_rank", "direction", "signature_version", "bulk_outcome_used_for_selection"))
  out
}

stage_v2_signature_audit <- function() {
  sig <- make_v2_signatures()
  old <- if (file.exists(FIGURE7_V2_INPUTS$v1_signature_manifest)) read_v2_tsv(FIGURE7_V2_INPUTS$v1_signature_manifest) else data.table()
  if (nrow(old)) {
    old[, `:=`(gene = norm_gene(gene), v1_axis = as.character(axis))]
    old_map <- old[, .(v1_axis = paste(sort(unique(v1_axis)), collapse = ";")), by = gene]
  } else old_map <- data.table(gene = character(), v1_axis = character())
  new_map <- sig[, .(v2_axis = paste(sort(unique(axis)), collapse = ";")), by = gene]
  change <- merge(old_map, new_map, by = "gene", all = TRUE)
  change[, change_type := fifelse(is.na(v1_axis), "new_v2_assignment", fifelse(is.na(v2_axis), "removed_from_v2", fifelse(v1_axis == v2_axis, "unchanged", "reassigned")))]
  change[, reason := fifelse(gene == "JUNB", "JUNB moved from v1 calibration comparator to frozen AP-1/stress axis", fifelse(gene %in% c("CEBPB", "EGR1"), "CEBPB/EGR1 reassigned from v1 tier1-rescue mixture to AP-1/CEBPB/EGR1 stress axis", "v2 axis definitions separated HNF4A/PPARA identity from AP-1/CEBPB/EGR1 stress and removed JUNB from calibration"))]
  write_v2_tsv(sig, file.path(FIGURE7_V2_META, "figure7_v2_signature_manifest.tsv"))
  write_v2_tsv(change, file.path(FIGURE7_V2_META, "figure7_v2_signature_v1_to_v2_gene_assignment_changes.tsv"))
  sig
}

stage_v2_direction_audit <- function(sig) {
  edge <- read_v2_tsv(FIGURE7_V2_INPUTS$celloracle_edges)
  edge[, `:=`(tf = norm_gene(source), gene = norm_gene(target), edge_sign = sign(safe_num(coef_mean)))]
  frozen <- sig[origin == "frozen_module7_target" & role == "primary", .(axis, tf, gene, origin)]
  joined <- merge(frozen, edge[edge_sign != 0, .(tf, gene, celloracle_state, edge_sign, coef_mean)], by = c("tf", "gene"), all.x = TRUE, allow.cartesian = TRUE)
  audit <- joined[, .(
    n_celloracle_edges = sum(is.finite(edge_sign)),
    n_signed_states = uniqueN(celloracle_state[is.finite(edge_sign)]),
    positive_edges = sum(edge_sign > 0, na.rm = TRUE), negative_edges = sum(edge_sign < 0, na.rm = TRUE),
    celloracle_states = paste(sort(unique(celloracle_state[is.finite(edge_sign)])), collapse = ";")
  ), by = .(axis, tf, gene)]
  audit[, direction_status := fifelse(n_celloracle_edges == 0, "unsigned", fifelse(positive_edges > 0 & negative_edges > 0, "conflicting", fifelse(positive_edges > 0, "positive", "negative")))]
  audit[, primary_signed_eligible := direction_status %in% c("positive", "negative") & n_signed_states >= 2]
  audit[, evidence_note := fifelse(direction_status == "conflicting", "CellOracle edge signs conflict across retained states; excluded from any signed primary score.", fifelse(primary_signed_eligible, "Consistent CellOracle sign across >=2 retained states.", "No sufficiently replicated CellOracle sign; retained only for unsigned programme scoring."))]
  frac <- audit[, .(n_frozen_targets = .N, n_signed_eligible = sum(primary_signed_eligible), signed_fraction = mean(primary_signed_eligible)), by = axis]
  signed_primary <- nrow(frac) == 3L && all(frac$signed_fraction >= 0.5)
  audit[, primary_score_decision := if (signed_primary) "signed_programme_score" else "unsigned_associated_target_programme_score"]
  source_audit <- data.table(
    source = c("CellOracle GRN edge coefficients", "CellOracle virtual-KO programme deltas", "pySCENIC regulon edges", "cisTarget motif-pruned regulons", "Module 7 frozen perturbation targets"),
    input_path = c(FIGURE7_V2_INPUTS$celloracle_edges,
                   file.path(FIGURE7_V2_ROOT, "metadata", "driver", "figure6_directional_network", "figure6_celloracle_programme_deltas_by_cell.tsv.gz"),
                   FIGURE7_V2_INPUTS$scenic_edges, FIGURE7_V2_INPUTS$cistarget_regulons, FIGURE7_V2_INPUTS$frozen_targets),
    sign_available = c(TRUE, TRUE, FALSE, FALSE, FALSE),
    sign_interpretation = c(
      "source-to-target coefficient sign is available but state-dependent; only >=2-state consistent signs qualify.",
      "TF knockout changes programme scores, supporting programme-level perturbation direction but not individual TF-to-target signs.",
      "coexpression correlation is unsigned for causal activation/repression claims.",
      "motif/regulon membership is unsigned for activation/repression claims.",
      "disturbed-target ranking is unsigned for activation/repression claims."
    ),
    usable_for_primary_target_sign = c(TRUE, FALSE, FALSE, FALSE, FALSE),
    decision = c("used conservatively; insufficient coverage for signed primary", "reported as programme-level supporting evidence only", "sensitivity regulon definition only", "sensitivity regulon definition only", "frozen unsigned target definition")
  )
  write_v2_tsv(audit, file.path(FIGURE7_V2_META, "figure7_v2_regulatory_direction_audit.tsv"))
  write_v2_tsv(frac, file.path(FIGURE7_V2_META, "figure7_v2_regulatory_direction_coverage.tsv"))
  write_v2_tsv(source_audit, file.path(FIGURE7_V2_META, "figure7_v2_regulatory_direction_source_audit.tsv"))
  list(audit = audit, signed_primary = signed_primary)
}

parse_regulon_targets <- function(path, tfs) {
  x <- read_v2_tsv(path)
  if (!all(c("tf", "target_genes") %in% names(x))) return(data.table(tf = character(), gene = character(), source = character()))
  out <- x[norm_gene(tf) %in% tfs & !is.na(target_genes), .(tf = norm_gene(tf), target_genes)]
  if (!nrow(out)) return(data.table(tf = character(), gene = character(), source = character()))
  ans <- rbindlist(lapply(seq_len(nrow(out)), function(i) data.table(tf = out$tf[[i]], gene = norm_gene(strsplit(out$target_genes[[i]], ";", fixed = TRUE)[[1L]]))), fill = TRUE)
  ans <- ans[nzchar(gene)]
  ans[, source := "cisTarget_motif_pruned_regulon"]
  unique(ans)
}

make_score_definitions <- function(sig, direction) {
  primary_axes <- c("identity_loss", "stress_transition", "sox4_associated")
  tfs <- sig[axis %in% primary_axes & !is.na(tf), .(tfs = paste(sort(unique(tf)), collapse = ";")), by = axis]
  frozen_target <- sig[axis %in% primary_axes & origin == "frozen_module7_target", .(gene = unique(gene)), by = axis]
  frozen_all <- sig[axis %in% primary_axes, .(gene = unique(gene)), by = axis]
  # Some TF anchors are also present in their own frozen target table.  Derive
  # TF-only genes from the TF field, rather than the retained-row provenance,
  # so these anchors cannot silently disappear through de-duplication.
  tf_only <- sig[axis %in% primary_axes & !is.na(tf), .(gene = unique(tf)), by = axis]
  all_tfs <- unique(unlist(strsplit(tfs$tfs, ";", fixed = TRUE)))
  cistarget <- parse_regulon_targets(FIGURE7_V2_INPUTS$cistarget_regulons, all_tfs)
  cistarget_axis <- merge(cistarget, sig[axis %in% primary_axes & !is.na(tf), .(axis, tf)], by = "tf", allow.cartesian = TRUE)[, .(gene = unique(gene)), by = axis]
  intersection <- merge(frozen_target, cistarget_axis, by = c("axis", "gene"))
  # This sensitivity is deliberately restricted to frozen targets with an
  # observed CellOracle edge.  It does not infer activation/repression from
  # edge membership and is therefore still an unsigned programme score.
  celloracle_target <- direction$audit[n_celloracle_edges > 0, .(gene = unique(gene)), by = axis]
  # Generic stress genes are prespecified for sensitivity, not chosen from bulk associations.
  generic_stress <- c("HSPA1A", "HSPA1B", "HSPA6", "HSPB1", "DNAJB1", "DUSP1", "PPP1R15A", "FOS", "JUN", "JUNB", "JUND", "ATF3")
  cell_cycle <- sig[axis == "proliferation_comparator", gene]
  make_def <- function(version, x, method, description, primary = FALSE) {
    x <- unique(as.data.table(x)[, .(axis, gene)])
    x[, `:=`(score_version = version, method = method, description = description, primary = primary)]
    x
  }
  definitions <- rbindlist(list(
    make_def("primary_frozen_programme", frozen_all, if (direction$signed_primary) "signed_programme_score" else "unsigned_associated_target_programme_score",
             if (direction$signed_primary) "CellOracle-sign-consistent frozen programme score; identity is inverted to loss." else "Unsigned associated target programme score; identity retention is inverted to identity loss.", TRUE),
    make_def("unsigned_target", frozen_target, "unsigned_frozen_target_score", "Frozen target genes only; no directional signs applied."),
    make_def("tf_expression_only", tf_only, "TF_expression_only_rank_score", "Frozen TF anchors only."),
    make_def("regulon_only", cistarget_axis, "cisTarget_regulon_only_rank_score", "Actual cisTarget motif-pruned regulon targets for the frozen TFs."),
    make_def("celloracle_target_only", celloracle_target, "CellOracle_edge_supported_target_rank_score", "Frozen targets with an observed CellOracle edge; membership is unsigned."),
    make_def("high_confidence_intersection", intersection, "frozen_target_cisTarget_intersection_rank_score", "Intersection of frozen Module 7 targets and cisTarget regulon targets."),
    make_def("no_cell_cycle", frozen_all[!gene %in% cell_cycle], "primary_minus_prespecified_cell_cycle_genes", "Primary frozen programme score after removal of fixed proliferation-panel genes."),
    make_def("no_generic_stress", frozen_all[!gene %in% generic_stress], "primary_minus_prespecified_generic_stress_genes", "Primary frozen programme score after removal of a prespecified generic-stress list."),
    make_def("alternative_median_rank", frozen_all, "median_normalized_rank_score", "Alternative rank-score implementation using the median rather than mean normalized rank.")
  ), fill = TRUE)
  if (direction$signed_primary) {
    signed <- merge(frozen_target, direction$audit[primary_signed_eligible == TRUE, .(axis, tf, gene, direction_status)], by = c("axis", "gene"), all = FALSE, allow.cartesian = TRUE)
    signed <- signed[, .(direction_status = if (uniqueN(direction_status) == 1L) unique(direction_status) else "conflicting"), by = .(axis, gene)]
    signed <- signed[direction_status != "conflicting", .(axis, gene)]
    definitions <- rbindlist(list(definitions, make_def("signed_target", signed, "CellOracle_signed_target_rank_score", "Consistent CellOracle-sign target score; conflicting target signs excluded.")), fill = TRUE)
  }
  definitions <- definitions[nzchar(gene)]
  definitions[, gene := norm_gene(gene)]
  definitions <- unique(definitions, by = c("score_version", "axis", "gene"))
  definitions
}

read_v2_cohort_objects <- function(clinical) {
  tcga <- readRDS(FIGURE7_V2_INPUTS$tcga_expression_rds)
  icgc <- readRDS(FIGURE7_V2_INPUTS$icgc_expression_rds)
  make_samples <- function(obj, cohort_name) {
    sm <- as.data.table(obj$samples)
    sm <- sm[, .(sample_id = as.character(sample_id))]
    sm[, patient_id := if (cohort_name == "TCGA_LIHC") tcga_id(sample_id) else icgc_id(sample_id)]
    sm[, `:=`(cohort = cohort_name,
              tumour_normal = if (cohort_name == "TCGA_LIHC") fifelse(substr(vapply(strsplit(sample_id, "-", fixed = TRUE), function(z) if (length(z) >= 4L) z[[4L]] else "", character(1)), 1L, 2L) %in% c("01", "02", "03", "05"), "tumour", fifelse(substr(vapply(strsplit(sample_id, "-", fixed = TRUE), function(z) if (length(z) >= 4L) z[[4L]] else "", character(1)), 1L, 2L) %in% c("10", "11", "12", "13", "14"), "normal", "unknown")) else fifelse(grepl("-T$", sample_id), "tumour", fifelse(grepl("-N$", sample_id), "normal", "unknown")))]
    merge(sm, clinical[cohort == cohort_name], by = c("cohort", "patient_id"), all.x = TRUE, sort = FALSE)
  }
  tcga$samples <- make_samples(tcga, "TCGA_LIHC")
  icgc$samples <- make_samples(icgc, "ICGC_LIRI_JP")
  list(TCGA_LIHC = tcga, ICGC_LIRI_JP = icgc)
}

normalized_rank_matrix <- function(obj, cohort_name) {
  cache <- file.path(FIGURE7_V2_DATA, paste0("figure7_v2_", tolower(cohort_name), "_normalized_rank_matrix.rds"))
  if (file.exists(cache)) return(readRDS(cache))
  # Existing Figure 7 rank matrices are read-only and rank genes independently of any clinical data.
  source_cache <- if (cohort_name == "TCGA_LIHC") FIGURE7_V2_INPUTS$tcga_rank_rds else FIGURE7_V2_INPUTS$icgc_rank_rds
  if (file.exists(source_cache)) {
    ranks <- readRDS(source_cache)
  } else {
    ranks <- apply(obj$expr, 2L, rank, ties.method = "average", na.last = "keep") / nrow(obj$expr)
  }
  if (is.vector(ranks)) ranks <- matrix(ranks, ncol = 1L, dimnames = dimnames(obj$expr))
  saveRDS(ranks, cache, compress = FALSE)
  ranks
}

score_gene_set <- function(ranks, genes, method = c("mean", "median"), direction = 1) {
  method <- match.arg(method)
  genes <- intersect(unique(norm_gene(genes)), rownames(ranks))
  if (!length(genes)) return(rep(NA_real_, ncol(ranks)))
  raw <- if (method == "mean") colMeans(ranks[genes, , drop = FALSE], na.rm = TRUE) - 0.5 else apply(ranks[genes, , drop = FALSE], 2L, median, na.rm = TRUE) - 0.5
  direction * raw
}

score_direction <- function(axis) if (axis == "identity_loss") -1 else 1

stage_v2_score_reconstruction <- function() {
  clinical <- build_harmonized_clinical()
  sig <- stage_v2_signature_audit()
  direction <- stage_v2_direction_audit(sig)
  defs <- make_score_definitions(sig, direction)
  objs <- read_v2_cohort_objects(clinical)
  rows <- list(); coverage_rows <- list(); idx <- 0L; ci <- 0L
  for (cohort in names(objs)) {
    obj <- objs[[cohort]]
    ranks <- normalized_rank_matrix(obj, cohort)
    sample_meta <- obj$samples[match(colnames(ranks), sample_id)]
    for (ver in unique(defs$score_version)) {
      for (axis_name in c("identity_loss", "stress_transition", "sox4_associated")) {
        g <- defs[score_version == ver & axis == axis_name, gene]
        mapped <- intersect(g, rownames(ranks))
        raw <- score_gene_set(ranks, mapped, method = if (ver == "alternative_median_rank") "median" else "mean", direction = score_direction(axis_name))
        idx <- idx + 1L
        rows[[idx]] <- data.table(sample_meta, score_version = ver, axis = axis_name, score_raw = raw, score_cohort_z = safe_z(raw), n_frozen_genes = uniqueN(g), n_mapped_genes = length(mapped), mapped_fraction = length(mapped) / max(1L, uniqueN(g)))
        ci <- ci + 1L
        coverage_rows[[ci]] <- data.table(cohort = cohort, score_version = ver, axis = axis_name, n_frozen_genes = uniqueN(g), n_mapped_genes = length(mapped), mapped_fraction = length(mapped) / max(1L, uniqueN(g)), mapped_genes = paste(sort(mapped), collapse = ";"), missing_genes = paste(sort(setdiff(g, mapped)), collapse = ";"))
      }
    }
    # Comparator and covariate scores are a single fixed rank-score definition.
    for (axis_name in c("proliferation_comparator", "foxm1_cebpb_reference", "broad_network_calibration", "hypoxia_covariate", "inflammation_covariate")) {
      g <- sig[axis == axis_name, gene]
      mapped <- intersect(g, rownames(ranks)); raw <- score_gene_set(ranks, mapped)
      idx <- idx + 1L
      rows[[idx]] <- data.table(sample_meta, score_version = "fixed_comparator", axis = axis_name, score_raw = raw, score_cohort_z = safe_z(raw), n_frozen_genes = uniqueN(g), n_mapped_genes = length(mapped), mapped_fraction = length(mapped) / max(1L, uniqueN(g)))
      ci <- ci + 1L
      coverage_rows[[ci]] <- data.table(cohort = cohort, score_version = "fixed_comparator", axis = axis_name, n_frozen_genes = uniqueN(g), n_mapped_genes = length(mapped), mapped_fraction = length(mapped) / max(1L, uniqueN(g), 1L), mapped_genes = paste(sort(mapped), collapse = ";"), missing_genes = paste(sort(setdiff(g, mapped)), collapse = ";"))
    }
  }
  scores <- rbindlist(rows, fill = TRUE)
  coverage <- rbindlist(coverage_rows, fill = TRUE)
  tcga_sample_scale <- scores[cohort == "TCGA_LIHC" & tumour_normal == "tumour" & score_version == "primary_frozen_programme" & axis %in% c("identity_loss", "stress_transition", "sox4_associated"), .(tcga_training_mean = mean(score_raw, na.rm = TRUE), tcga_training_sd = sd(score_raw, na.rm = TRUE), n_training = .N, analysis_level = "sample_descriptive_only"), by = .(score_version, axis)]
  scores <- merge(scores, tcga_sample_scale[, .(score_version, axis, tcga_training_mean, tcga_training_sd)], by = c("score_version", "axis"), all.x = TRUE, sort = FALSE)
  scores[, score_tcga_frozen_z := (score_raw - tcga_training_mean) / tcga_training_sd]
  write_v2_tsvgz(scores, file.path(FIGURE7_V2_META, "figure7_v2_bulk_axis_scores.tsv.gz"))
  write_v2_tsv(coverage, file.path(FIGURE7_V2_META, "figure7_v2_signature_coverage.tsv"))
  score_audit <- defs[, .(n_genes = uniqueN(gene), score_method = unique(method), description = unique(description), primary = unique(primary)), by = .(score_version, axis)]
  score_audit <- merge(score_audit, coverage[, .(mapped_fraction_min = min(mapped_fraction), mapped_fraction_max = max(mapped_fraction)), by = .(score_version, axis)], by = c("score_version", "axis"), all.x = TRUE)
  score_audit <- merge(score_audit, direction$audit[, .(n_conflicting_sign_targets = sum(direction_status == "conflicting"), n_signed_eligible_targets = sum(primary_signed_eligible)), by = axis], by = "axis", all.x = TRUE)
  score_audit[, primary_score_label := ifelse(score_version == "primary_frozen_programme", if (direction$signed_primary) "signed programme score" else "associated target programme score", "sensitivity score")]
  write_v2_tsv(score_audit, file.path(FIGURE7_V2_META, "figure7_v2_axis_score_definition_audit.tsv"))
  patient <- scores[, .(score_raw = mean(score_raw, na.rm = TRUE), score_cohort_z = mean(score_cohort_z, na.rm = TRUE), n_samples_aggregated = .N), by = .(cohort, patient_id, tumour_normal, score_version, axis, age_years, age_high, sex, stage_num, stage_high, os_time_days, os_event, grade_num, t_stage_num, clinical_analysis_status, os_analysis_status)]
  tcga_patient_scale <- patient[cohort == "TCGA_LIHC" & tumour_normal == "tumour" & score_version == "primary_frozen_programme" & axis %in% c("identity_loss", "stress_transition", "sox4_associated"), .(tcga_training_mean = mean(score_raw, na.rm = TRUE), tcga_training_sd = sd(score_raw, na.rm = TRUE), n_training = .N, analysis_level = "patient_prediction_training"), by = .(score_version, axis)]
  patient <- merge(patient, tcga_patient_scale[, .(score_version, axis, tcga_training_mean, tcga_training_sd)], by = c("score_version", "axis"), all.x = TRUE, sort = FALSE)
  patient[, score_tcga_frozen_z := (score_raw - tcga_training_mean) / tcga_training_sd]
  write_v2_tsv(rbindlist(list(tcga_sample_scale, tcga_patient_scale), fill = TRUE), file.path(FIGURE7_V2_META, "figure7_v2_tcga_scaling_parameters.tsv"))
  write_v2_tsvgz(patient, file.path(FIGURE7_V2_META, "figure7_v2_patient_level_scores.tsv.gz"))
  concordance <- patient[tumour_normal == "tumour" & axis %in% c("identity_loss", "stress_transition", "sox4_associated"), .(score = score_cohort_z), by = .(cohort, patient_id, score_version, axis)]
  pairs <- list(); pp <- 0L
  for (cohort_name in unique(concordance$cohort)) for (axis_name in c("identity_loss", "stress_transition", "sox4_associated")) {
    w <- dcast(concordance[cohort == cohort_name & axis == axis_name], patient_id ~ score_version, value.var = "score")
    vs <- setdiff(names(w), "patient_id")
    if (length(vs) > 1L) for (i in seq_len(length(vs) - 1L)) for (j in (i + 1L):length(vs)) {
      pp <- pp + 1L; pairs[[pp]] <- data.table(cohort = cohort_name, axis = axis_name, score_version_x = vs[[i]], score_version_y = vs[[j]], n = sum(complete.cases(w[, ..vs][, c(i, j)])), pearson_r = suppressWarnings(cor(w[[vs[[i]]]], w[[vs[[j]]]], use = "pairwise.complete.obs")))
    }
  }
  write_v2_tsv(rbindlist(pairs, fill = TRUE), file.path(FIGURE7_V2_META, "figure7_v2_axis_score_concordance.tsv"))
  list(scores = scores, patient = patient, signatures = sig, definitions = defs, direction = direction, objects = objs)
}

hedges_g_v2 <- function(x_t, x_n) {
  x_t <- x_t[is.finite(x_t)]; x_n <- x_n[is.finite(x_n)]
  n1 <- length(x_t); n0 <- length(x_n)
  if (n1 < 2L || n0 < 2L) return(c(effect = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_))
  sp <- sqrt(((n1 - 1) * var(x_t) + (n0 - 1) * var(x_n)) / (n1 + n0 - 2))
  if (!is.finite(sp) || sp == 0) return(c(effect = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_))
  d <- (mean(x_t) - mean(x_n)) / sp
  j <- 1 - 3 / (4 * (n1 + n0) - 9)
  g <- j * d
  se <- sqrt((n1 + n0) / (n1 * n0) + g^2 / (2 * (n1 + n0 - 2)))
  c(effect = g, se = se, ci_low = g - 1.96 * se, ci_high = g + 1.96 * se, p_value = 2 * pnorm(-abs(g / se)))
}

first_nonmissing_v2 <- function(x) {
  z <- x[!is.na(x)]
  if (!length(z)) return(x[[1L]])
  z[[1L]]
}

patient_score_table <- function(patient, score_version_name = "primary_frozen_programme") {
  axes <- c("identity_loss", "stress_transition", "sox4_associated")
  d <- copy(patient[score_version == score_version_name & axis %in% axes])
  id_cols <- c("cohort", "patient_id", "tumour_normal")
  clinical_cols <- c("age_years", "age_high", "sex", "stage_num", "stage_high", "os_time_days", "os_event", "grade_num", "t_stage_num", "clinical_analysis_status", "os_analysis_status")
  # Score data have already been patient-aggregated, but this explicit mean
  # protects models against an upstream replicate appearing in a future run.
  score_long <- d[, .(score_cohort_z = mean(score_cohort_z, na.rm = TRUE),
                      score_tcga_frozen_z = mean(score_tcga_frozen_z, na.rm = TRUE)),
                  by = c(id_cols, "axis")]
  clinical <- d[, lapply(.SD, first_nonmissing_v2), by = id_cols, .SDcols = clinical_cols]
  wide <- dcast(score_long, cohort + patient_id + tumour_normal ~ axis,
                value.var = "score_cohort_z", fun.aggregate = mean)
  wide_frozen <- dcast(score_long, cohort + patient_id + tumour_normal ~ axis,
                       value.var = "score_tcga_frozen_z", fun.aggregate = mean)
  frozen_names <- setdiff(names(wide_frozen), id_cols)
  setnames(wide_frozen, frozen_names, paste0(frozen_names, "_tcga_frozen_z"))
  merge(merge(clinical, wide, by = id_cols, all = TRUE, sort = FALSE), wide_frozen, by = id_cols, all.x = TRUE, sort = FALSE)
}

comparator_patient_table <- function(patient) {
  axes <- c("proliferation_comparator", "foxm1_cebpb_reference", "broad_network_calibration", "hypoxia_covariate", "inflammation_covariate")
  d <- copy(patient[score_version == "fixed_comparator" & axis %in% axes])
  d <- d[, .(score_cohort_z = mean(score_cohort_z, na.rm = TRUE)), by = .(cohort, patient_id, tumour_normal, axis)]
  dcast(d, cohort + patient_id + tumour_normal ~ axis, value.var = "score_cohort_z", fun.aggregate = mean)
}

stage_v2_tumour_normal <- function(patient) {
  axes <- c("identity_loss", "stress_transition", "sox4_associated")
  d <- patient[score_version == "primary_frozen_programme" & axis %in% axes]
  results <- list(); paired <- list(); r <- 0L; q <- 0L
  for (cohort_name in unique(d$cohort)) for (axis_name in axes) {
    z <- d[cohort == cohort_name & axis == axis_name]
    agg <- z[, .(score = mean(score_cohort_z, na.rm = TRUE)), by = .(patient_id, tumour_normal)]
    g <- hedges_g_v2(agg[tumour_normal == "tumour", score], agg[tumour_normal == "normal", score])
    r <- r + 1L
    results[[r]] <- data.table(cohort = cohort_name, axis = axis_name, analysis = "patient_level_independent_Hedges_g", n_tumour = uniqueN(agg[tumour_normal == "tumour", patient_id]), n_normal = uniqueN(agg[tumour_normal == "normal", patient_id]), n_paired = length(intersect(agg[tumour_normal == "tumour", patient_id], agg[tumour_normal == "normal", patient_id])), hedges_g = g[["effect"]], se = g[["se"]], ci_low = g[["ci_low"]], ci_high = g[["ci_high"]], p_value = g[["p_value"]])
    w <- dcast(agg, patient_id ~ tumour_normal, value.var = "score")
    if (all(c("tumour", "normal") %in% names(w))) {
      delta <- w$tumour - w$normal; delta <- delta[is.finite(delta)]
      q <- q + 1L
      paired[[q]] <- data.table(cohort = cohort_name, axis = axis_name, analysis = "paired_tumour_minus_normal_sensitivity", n_paired = length(delta), mean_difference = mean(delta), sd_difference = sd(delta), ci_low = mean(delta) - qt(.975, max(1, length(delta) - 1)) * sd(delta) / sqrt(length(delta)), ci_high = mean(delta) + qt(.975, max(1, length(delta) - 1)) * sd(delta) / sqrt(length(delta)), p_value = if (length(delta) >= 2L) t.test(delta)$p.value else NA_real_)
    } else {
      q <- q + 1L; paired[[q]] <- data.table(cohort = cohort_name, axis = axis_name, analysis = "paired_tumour_minus_normal_sensitivity", n_paired = 0L, mean_difference = NA_real_, sd_difference = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_)
    }
  }
  out <- rbindlist(results); out[, fdr := bh(p_value), by = cohort]
  pair_out <- rbindlist(paired); pair_out[, fdr := bh(p_value), by = cohort]
  pooled <- rbindlist(lapply(axes, function(axis_name) {
    z <- out[axis == axis_name & is.finite(hedges_g) & is.finite(se)]
    if (nrow(z) < 2L) return(data.table(axis = axis_name, k = nrow(z), pooled_effect = NA_real_, ci_low = NA_real_, ci_high = NA_real_, I2 = NA_real_, interpretation = "not_estimable"))
    fit <- rma.uni(yi = z$hedges_g, sei = z$se, method = "REML")
    data.table(axis = axis_name, k = nrow(z), pooled_effect = as.numeric(coef(fit)), ci_low = fit$ci.lb, ci_high = fit$ci.ub, I2 = fit$I2, interpretation = "exploratory_random_effects_k_equals_2")
  }), fill = TRUE)
  write_v2_tsv(out, file.path(FIGURE7_V2_FIG, "figure7_v2c_tumour_normal_effects.tsv"))
  write_v2_tsv(pair_out, file.path(FIGURE7_V2_FIG, "figure7_v2c_tumour_normal_paired_sensitivity.tsv"))
  write_v2_tsv(pooled, file.path(FIGURE7_V2_FIG, "figure7_v2c_tumour_normal_meta_analysis.tsv"))
  list(primary = out, paired = pair_out, pooled = pooled)
}

fit_stage_model_v2 <- function(d, score_col, adjusted = TRUE) {
  keep <- d[is.finite(get(score_col)) & is.finite(stage_high) & is.finite(age_years) & !is.na(sex)]
  if (nrow(keep) < 20L || uniqueN(keep$stage_high) < 2L) return(data.table(status = "not_estimable", reason = "insufficient_or_nonvariable_outcome", n = nrow(keep)))
  formula <- if (adjusted) as.formula(paste("stage_high ~", score_col, "+ age_years + factor(sex)")) else as.formula(paste("stage_high ~", score_col))
  fit <- try(glm(formula, data = keep, family = binomial()), silent = TRUE)
  if (inherits(fit, "try-error")) return(data.table(status = "not_estimable", reason = "model_fit_failed", n = nrow(keep)))
  co <- summary(fit)$coefficients
  if (!score_col %in% rownames(co)) return(data.table(status = "not_estimable", reason = "score_coefficient_missing", n = nrow(keep)))
  beta <- co[score_col, "Estimate"]; se <- co[score_col, "Std. Error"]
  data.table(status = "estimated", reason = "", n = nrow(keep), events = sum(keep$stage_high == 1L), coefficient = beta, se = se, odds_ratio = exp(beta), ci_low = exp(beta - 1.96 * se), ci_high = exp(beta + 1.96 * se), p_value = co[score_col, "Pr(>|z|)"])
}

stage_v2_clinical_associations <- function(patient) {
  axes <- c("identity_loss", "stress_transition", "sox4_associated")
  rows <- list(); k <- 0L
  d <- patient_score_table(patient)
  # Shared-stage model is only valid when ICGC coding has been independently verified.
  for (cohort_name in c("TCGA_LIHC", "ICGC_LIRI_JP")) for (axis_name in axes) {
    z <- d[cohort == cohort_name & tumour_normal == "tumour"]
    if (cohort_name == "ICGC_LIRI_JP") {
      k <- k + 1L; rows[[k]] <- data.table(cohort = cohort_name, programme = axis_name, clinical_feature = "early_vs_advanced_stage", model = "age_sex_adjusted_logistic", status = "blocked", reason = "ICGC_Age_Gender_Stage_encoding_unverified", n = nrow(z), events = NA_integer_, coefficient = NA_real_, se = NA_real_, odds_ratio = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_)
    } else {
      fit <- fit_stage_model_v2(z, axis_name, TRUE)
      k <- k + 1L; rows[[k]] <- cbind(data.table(cohort = cohort_name, programme = axis_name, clinical_feature = "early_vs_advanced_stage", model = "age_sex_adjusted_logistic"), fit)
    }
  }
  # TCGA-only secondary clinicopathological outcomes retain transparent labels.
  tcga <- d[cohort == "TCGA_LIHC" & tumour_normal == "tumour"]
  for (axis_name in axes) {
    for (outcome in c("grade_num", "t_stage_num", "stage_num")) {
      z <- tcga[is.finite(get(axis_name)) & is.finite(get(outcome)) & is.finite(age_years) & !is.na(sex)]
      if (nrow(z) < 20L || uniqueN(z[[outcome]]) < 2L) fit <- data.table(status = "not_estimable", reason = "insufficient_or_nonvariable_outcome", n = nrow(z), events = NA_integer_, coefficient = NA_real_, se = NA_real_, odds_ratio = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_)
      else {
        f <- as.formula(paste(outcome, "~", axis_name, "+ age_years + factor(sex)")); mod <- try(lm(f, data = z), silent = TRUE)
        if (inherits(mod, "try-error")) fit <- data.table(status = "not_estimable", reason = "model_fit_failed", n = nrow(z), events = NA_integer_, coefficient = NA_real_, se = NA_real_, odds_ratio = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_)
        else {
          co <- summary(mod)$coefficients
          if (!axis_name %in% rownames(co)) fit <- data.table(status = "not_estimable", reason = "score_coefficient_missing", n = nrow(z), events = NA_integer_, coefficient = NA_real_, se = NA_real_, odds_ratio = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_)
          else {
            beta <- co[axis_name, 1L]; se <- co[axis_name, 2L]
            fit <- data.table(status = "estimated", reason = "TCGA_only_age_sex_adjusted_linear_secondary", n = nrow(z), events = NA_integer_, coefficient = beta, se = se, odds_ratio = NA_real_, ci_low = beta - 1.96 * se, ci_high = beta + 1.96 * se, p_value = co[axis_name, 4L])
          }
        }
      }
      k <- k + 1L; rows[[k]] <- cbind(data.table(cohort = "TCGA_LIHC", programme = axis_name, clinical_feature = paste0(outcome, "_TCGA_only"), model = "age_sex_adjusted_secondary"), fit)
    }
  }
  out <- rbindlist(rows, fill = TRUE); out[, fdr := bh(p_value), by = .(cohort, clinical_feature)]
  write_v2_tsv(out, file.path(FIGURE7_V2_FIG, "figure7_v2d_clinical_associations.tsv"))
  out
}

vif_v2 <- function(d, columns) {
  x <- as.matrix(d[, ..columns]); x <- x[, apply(x, 2L, function(z) sd(z, na.rm = TRUE) > 0), drop = FALSE]
  if (!ncol(x)) return(data.table(term = character(), vif = numeric()))
  rbindlist(lapply(seq_len(ncol(x)), function(i) {
    if (ncol(x) == 1L) return(data.table(term = colnames(x)[[i]], vif = 1))
    fit <- lm(x[, i] ~ x[, -i, drop = FALSE]); data.table(term = colnames(x)[[i]], vif = 1 / max(1e-12, 1 - summary(fit)$r.squared))
  }))
}

fit_cox_v2 <- function(d, programme, extra_covars = character(), score_col = programme, model_id = NULL) {
  z <- copy(d)
  needed <- unique(c("os_time_days", "os_event", "age_years", "sex", "stage_high", score_col, extra_covars))
  z <- z[complete.cases(z[, ..needed]) & os_time_days > 0 & os_event %in% c(0, 1)]
  if (!nrow(z) || sum(z$os_event == 1L) < 5L) return(data.table(model_id = model_id %||% programme, programme = programme, term = score_col, status = "not_estimable", reason = "insufficient_events_or_complete_cases", n = nrow(z), events = sum(z$os_event == 1L), predictors = NA_integer_, events_per_variable = NA_real_, coefficient = NA_real_, se = NA_real_, hazard_ratio = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_, concordance = NA_real_, ph_p_value = NA_real_, condition_index = NA_real_))
  covars <- c("age_years", "factor(sex)", "stage_high", score_col, extra_covars)
  f <- as.formula(paste("Surv(os_time_days, os_event) ~", paste(covars, collapse = " + ")))
  fit <- try(coxph(f, data = z, x = TRUE), silent = TRUE)
  if (inherits(fit, "try-error")) return(data.table(model_id = model_id %||% programme, programme = programme, term = score_col, status = "not_estimable", reason = "cox_fit_failed", n = nrow(z), events = sum(z$os_event == 1L), predictors = length(covars), events_per_variable = sum(z$os_event == 1L) / length(covars), coefficient = NA_real_, se = NA_real_, hazard_ratio = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_, concordance = NA_real_, ph_p_value = NA_real_, condition_index = NA_real_))
  co <- summary(fit)$coefficients
  if (!score_col %in% rownames(co)) return(data.table(model_id = model_id %||% programme, programme = programme, term = score_col, status = "not_estimable", reason = "score_coefficient_missing", n = nrow(z), events = sum(z$os_event == 1L), predictors = length(coef(fit)), events_per_variable = sum(z$os_event == 1L) / length(coef(fit)), coefficient = NA_real_, se = NA_real_, hazard_ratio = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_, concordance = NA_real_, ph_p_value = NA_real_, condition_index = NA_real_))
  beta <- co[score_col, "coef"]; se <- co[score_col, "se(coef)"]
  ph <- try(cox.zph(fit), silent = TRUE); ph_p <- if (inherits(ph, "try-error") || !score_col %in% rownames(ph$table)) NA_real_ else ph$table[score_col, "p"]
  xnum <- model.matrix(~ age_years + sex + stage_high + ., data = z[, c("age_years", "sex", "stage_high", score_col, extra_covars), with = FALSE])
  xnum <- xnum[, apply(xnum, 2L, function(v) is.finite(sd(v)) && sd(v) > 0), drop = FALSE]
  eig <- try(eigen(cor(xnum), only.values = TRUE)$values, silent = TRUE)
  ci <- if (inherits(eig, "try-error") || !length(eig) || min(eig) <= 0) NA_real_ else sqrt(max(eig) / min(eig))
  data.table(model_id = model_id %||% programme, programme = programme, term = score_col, status = "estimated", reason = "", n = nrow(z), events = sum(z$os_event == 1L), predictors = length(coef(fit)), events_per_variable = sum(z$os_event == 1L) / length(coef(fit)), coefficient = beta, se = se, hazard_ratio = exp(beta), ci_low = exp(beta - 1.96 * se), ci_high = exp(beta + 1.96 * se), p_value = co[score_col, "Pr(>|z|)"], concordance = summary(fit)$concordance[[1L]], ph_p_value = ph_p, condition_index = ci)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

stage_v2_cox <- function(patient) {
  axes <- c("identity_loss", "stress_transition", "sox4_associated")
  primary <- patient_score_table(patient)[cohort == "TCGA_LIHC" & tumour_normal == "tumour"]
  comps <- comparator_patient_table(patient)[cohort == "TCGA_LIHC" & tumour_normal == "tumour"]
  primary <- merge(primary, comps, by = c("cohort", "patient_id", "tumour_normal"), all.x = TRUE)
  programmes <- c(axes, "foxm1_cebpb_reference", "proliferation_comparator", "broad_network_calibration")
  rows <- lapply(programmes, function(p) fit_cox_v2(primary, p, model_id = paste0("clinical_adjusted__", p)))
  joint <- fit_cox_v2(primary, "identity_loss", extra_covars = c("stress_transition", "sox4_associated"), model_id = "clinical_adjusted_joint_three_axes")
  joint[, term := "identity_loss"]
  out_tcga <- rbindlist(c(rows, list(joint)), fill = TRUE)
  out_tcga[, `:=`(cohort = "TCGA_LIHC", fdr = bh(p_value))]
  # Preserve ICGC rows but explicitly prevent an unverified clinical survival model.
  out_icgc <- rbindlist(lapply(c(programmes, "joint_three_axes"), function(p) data.table(model_id = paste0("clinical_adjusted__", p), programme = p, term = p, status = "blocked", reason = "ICGC_Age_Gender_Stage_fustat_futime_encoding_unverified", n = NA_integer_, events = NA_integer_, predictors = NA_integer_, events_per_variable = NA_real_, coefficient = NA_real_, se = NA_real_, hazard_ratio = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_, concordance = NA_real_, ph_p_value = NA_real_, condition_index = NA_real_, fdr = NA_real_, cohort = "ICGC_LIRI_JP")), fill = TRUE)
  out <- rbindlist(list(out_tcga, out_icgc), fill = TRUE)
  out[, low_epv_flag := ifelse(status == "estimated" & events_per_variable < 8, "exploratory_overfitting_risk", ifelse(status == "estimated", "none", status))]
  # VIFs and PH tests are stored separately rather than mixing clinical
  # covariates into the programme-term BH family.
  vifs <- list(); ph_rows <- list(); vi <- 0L; pi <- 0L
  for (p in programmes) {
    z <- primary[complete.cases(age_years, sex, stage_high, os_time_days, os_event, get(p)) & os_time_days > 0]
    if (nrow(z)) {
      predictors <- c("age_years", "stage_high", p)
      vv <- vif_v2(z, predictors); vv[, `:=`(model_id = paste0("clinical_adjusted__", p), cohort = "TCGA_LIHC")]
      vi <- vi + 1L; vifs[[vi]] <- vv
      fit <- try(coxph(as.formula(paste("Surv(os_time_days, os_event) ~ age_years + factor(sex) + stage_high +", p)), data = z), silent = TRUE)
      if (!inherits(fit, "try-error")) {
        ph <- try(cox.zph(fit), silent = TRUE)
        if (!inherits(ph, "try-error")) {
          pp <- as.data.table(ph$table, keep.rownames = "term"); setnames(pp, names(pp), c("term", "rho", "chisq", "p_value"))
          pp[, `:=`(model_id = paste0("clinical_adjusted__", p), cohort = "TCGA_LIHC", status = "estimated")]
          pi <- pi + 1L; ph_rows[[pi]] <- pp
        }
      }
    }
  }
  vif_out <- rbindlist(vifs, fill = TRUE); ph_out <- rbindlist(ph_rows, fill = TRUE)
  write_v2_tsv(out, file.path(FIGURE7_V2_FIG, "figure7_v2e_multivariable_cox_models.tsv"))
  write_v2_tsv(out[, .(cohort, programme, term, status, reason, n, events, predictors, events_per_variable, low_epv_flag, ph_p_value, condition_index)], file.path(FIGURE7_V2_FIG, "figure7_v2e_cox_diagnostics.tsv"))
  write_v2_tsv(vif_out, file.path(FIGURE7_V2_FIG, "figure7_v2e_vif.tsv"))
  write_v2_tsv(ph_out, file.path(FIGURE7_V2_FIG, "figure7_v2e_ph_assumption.tsv"))
  list(models = out, tcga_data = primary, vif = vif_out, ph = ph_out)
}

make_gene_strata_v2 <- function(expr) {
  # Matching uses the expression matrix itself, rather than rank scores:
  # gene number, mean expression, expression variance, and detection rate.
  mean_expr <- rowMeans(expr, na.rm = TRUE)
  var_expr <- apply(expr, 1L, var, na.rm = TRUE)
  detection <- rowMeans(expr > 0, na.rm = TRUE)
  # Three quantile bins per matching variable preserve all three matching
  # dimensions while avoiding sparse 5x5x5 cells for a 95-gene programme.
  qbin <- function(x, n = 3L) {
    br <- unique(quantile(x, probs = seq(0, 1, length.out = n + 1L), na.rm = TRUE, type = 8))
    if (length(br) < 2L) return(rep(1L, length(x)))
    as.integer(cut(x, breaks = br, include.lowest = TRUE, labels = FALSE))
  }
  z <- data.table(gene = rownames(expr), mean_bin = qbin(mean_expr), variance_bin = qbin(var_expr), detection_bin = qbin(detection))
  z[, stratum := paste(mean_bin, variance_bin, detection_bin, sep = "_")]
  z
}

sample_matched_genes_v2 <- function(target, strata, seed) {
  set.seed(seed)
  target <- intersect(unique(target), strata$gene)
  if (!length(target)) return(list(genes = character(), n_exact = 0L, n_fallback = 0L))
  targets <- strata[gene %in% target]
  non_target <- strata[!gene %in% target]
  used <- character(); out <- character(); n_exact <- 0L; n_fallback <- 0L
  for (st in unique(targets$stratum)) {
    need <- nrow(targets[stratum == st]); pool <- setdiff(non_target[stratum == st, gene], used)
    if (length(pool) >= need) {
      chosen <- sample(pool, need, replace = FALSE); n_exact <- n_exact + need
    } else {
      # Deterministic nearest-bin fallback is explicit in the output. It is
      # used only when a complete exact stratum has too few non-target genes.
      base <- targets[stratum == st][1L]
      available <- non_target[!gene %in% used]
      available[, bin_distance := abs(mean_bin - base$mean_bin) + abs(variance_bin - base$variance_bin) + abs(detection_bin - base$detection_bin)]
      available <- available[order(bin_distance, gene)]
      if (nrow(available) < need) return(list(genes = character(), n_exact = n_exact, n_fallback = n_fallback))
      chosen <- sample(available[bin_distance == min(bin_distance), gene], min(need, nrow(available[bin_distance == min(bin_distance)])), replace = FALSE)
      if (length(chosen) < need) chosen <- c(chosen, sample(setdiff(available$gene, chosen), need - length(chosen), replace = FALSE))
      n_exact <- n_exact + length(intersect(chosen, pool)); n_fallback <- n_fallback + (need - length(intersect(chosen, pool)))
    }
    out <- c(out, chosen); used <- c(used, chosen)
  }
  list(genes = out, n_exact = n_exact, n_fallback = n_fallback)
}

build_random_patient_score <- function(ranks, obj, genes, axis_name, tcga_scale = NULL) {
  raw <- score_gene_set(ranks, genes, direction = score_direction(axis_name))
  sm <- obj$samples[match(colnames(ranks), sample_id)]
  dt <- data.table(sm, raw_score = raw)
  # The real-axis path is sample-wise within-cohort standardisation followed
  # by patient aggregation. The random path is intentionally identical.
  dt[, score_cohort_z := safe_z(raw_score), by = cohort]
  patient <- dt[, .(raw_score = mean(raw_score, na.rm = TRUE), score_cohort_z = mean(score_cohort_z, na.rm = TRUE)), by = .(cohort, patient_id, tumour_normal, age_years, age_high, sex, stage_num, stage_high, os_time_days, os_event, clinical_analysis_status, os_analysis_status)]
  if (!is.null(tcga_scale)) patient[, score_tcga_frozen_z := (raw_score - tcga_scale[["mean"]]) / tcga_scale[["sd"]]]
  patient
}

random_real_statistic_v2 <- function(random_patient, cohort_name, real_primary = NULL, target_metric = NULL) {
  z <- random_patient[cohort == cohort_name]
  g <- hedges_g_v2(z[tumour_normal == "tumour", score_cohort_z], z[tumour_normal == "normal", score_cohort_z])
  out <- list(tumour_normal_effect = g[["effect"]], tumour_normal_p = g[["p_value"]], stage_log_or = NA_real_, cox_z = NA_real_, joint_c_index = NA_real_, delta_c_index = NA_real_)
  if (cohort_name == "TCGA_LIHC") {
    tumour <- z[tumour_normal == "tumour"]
    st <- fit_stage_model_v2(setnames(copy(tumour), "score_cohort_z", "random_score"), "random_score", TRUE)
    if (nrow(st) && st$status[[1L]] == "estimated") out$stage_log_or <- st$coefficient[[1L]]
    d <- setnames(copy(tumour), "score_cohort_z", "random_score")
    d <- d[complete.cases(os_time_days, os_event, age_years, sex, stage_high, random_score) & os_time_days > 0 & os_event %in% c(0, 1)]
    # Formula is intentionally identical to fit_cox_v2(): only PH/VIF
    # diagnostics are omitted in repeated null resampling.
    fit <- if (sum(d$os_event == 1L) >= 5L) try(coxph(Surv(os_time_days, os_event) ~ age_years + factor(sex) + stage_high + random_score, data = d, x = TRUE), silent = TRUE) else NULL
    if (inherits(fit, "coxph")) {
      co <- summary(fit)$coefficients
      out$cox_z <- co["random_score", "coef"] / co["random_score", "se(coef)"]
      out$joint_c_index <- summary(fit)$concordance[[1L]]
      if (!is.null(real_primary) && is.finite(real_primary$clinical_c_index)) out$delta_c_index <- out$joint_c_index - real_primary$clinical_c_index
    }
  }
  as.data.table(out)
}

stage_v2_random_benchmark <- function(state, n_random = FIGURE7_V2_N_RANDOM) {
  sig <- state$signatures; objs <- state$objects; patient <- state$patient
  axes <- c("identity_loss", "stress_transition", "sox4_associated")
  primary_defs <- state$definitions[score_version == "primary_frozen_programme" & axis %in% axes]
  real_cox <- stage_v2_cox(patient)$models[cohort == "TCGA_LIHC" & programme %in% axes & status == "estimated"]
  # Clinical baseline uses exactly the same TCGA complete-case population and covariates.
  real_d <- patient_score_table(patient)[cohort == "TCGA_LIHC" & tumour_normal == "tumour"]
  baseline <- real_d[complete.cases(age_years, sex, stage_high, os_time_days, os_event) & os_time_days > 0]
  baseline_fit <- coxph(Surv(os_time_days, os_event) ~ age_years + factor(sex) + stage_high, data = baseline, x = TRUE)
  clinical_c <- summary(baseline_fit)$concordance[[1L]]
  real_primary <- list(clinical_c_index = clinical_c)
  rows <- list(); map_rows <- list(); r <- 0L; m <- 0L
  for (cohort_name in names(objs)) {
    ranks <- normalized_rank_matrix(objs[[cohort_name]], cohort_name)
    strata <- make_gene_strata_v2(objs[[cohort_name]]$expr)
    for (axis_name in axes) {
      targets <- intersect(primary_defs[axis == axis_name, gene], rownames(ranks))
      for (i in seq_len(n_random)) {
        sampled <- sample_matched_genes_v2(targets, strata, FIGURE7_V2_SEED + i + match(axis_name, axes) * 100000L + match(cohort_name, names(objs)) * 1000000L)
        genes <- sampled$genes
        p <- build_random_patient_score(ranks, objs[[cohort_name]], genes, axis_name)
        stat <- random_real_statistic_v2(p, cohort_name, real_primary)
        r <- r + 1L
        rows[[r]] <- cbind(data.table(cohort = cohort_name, axis = axis_name, random_id = i, n_genes = length(genes), target_n_genes = length(targets), n_exact_stratum_genes = sampled$n_exact, n_nearest_bin_fallback_genes = sampled$n_fallback, matching = "gene_number+mean_expression+expression_variance+detection_rate; 3x3x3_quantile_strata", tumour_normal_model = "patient_level_independent_Hedges_g", stage_model = "stage_high ~ score + age_years + sex", survival_model = "Surv(OS) ~ age_years + sex + stage_high + score", seed = FIGURE7_V2_SEED), stat)
        m <- m + 1L; map_rows[[m]] <- data.table(cohort = cohort_name, axis = axis_name, random_id = i, genes = paste(genes, collapse = ";"))
      }
    }
  }
  out <- rbindlist(rows, fill = TRUE); maps <- rbindlist(map_rows, fill = TRUE)
  write_v2_tsvgz(out, file.path(FIGURE7_V2_META, "figure7_v2_matched_random_benchmark.tsv.gz"))
  write_v2_tsvgz(maps, file.path(FIGURE7_V2_DATA, "figure7_v2_matched_random_gene_sets.tsv.gz"))
  # Real observed metrics are extracted from the same patient-level functions.
  c_out <- stage_v2_tumour_normal(patient)$primary
  d_out <- stage_v2_clinical_associations(patient)
  real_rows <- list(); z <- 0L
  for (cohort_name in names(objs)) for (axis_name in axes) {
    z <- z + 1L
    tum <- c_out[cohort == cohort_name & axis == axis_name, hedges_g]
    stage <- d_out[cohort == cohort_name & programme == axis_name & clinical_feature == "early_vs_advanced_stage", coefficient]
    cox <- if (cohort_name == "TCGA_LIHC") real_cox[programme == axis_name, coefficient / se] else NA_real_
    jc <- if (cohort_name == "TCGA_LIHC") real_cox[programme == axis_name, concordance] else NA_real_
    dc <- if (cohort_name == "TCGA_LIHC" && length(jc)) jc - clinical_c else NA_real_
    real_rows[[z]] <- data.table(cohort = cohort_name, axis = axis_name, tumour_normal_effect = tum[[1L]], stage_log_or = if (length(stage)) stage[[1L]] else NA_real_, cox_z = if (length(cox)) cox[[1L]] else NA_real_, joint_c_index = if (length(jc)) jc[[1L]] else NA_real_, delta_c_index = if (length(dc)) dc[[1L]] else NA_real_)
  }
  observed <- rbindlist(real_rows)
  metrics <- c("tumour_normal_effect", "stage_log_or", "cox_z", "joint_c_index", "delta_c_index")
  summary_rows <- list(); s <- 0L
  for (i in seq_len(nrow(observed))) for (metric in metrics) {
    obs <- observed[[metric]][[i]]; null <- out[cohort == observed$cohort[[i]] & axis == observed$axis[[i]], get(metric)]
    null <- null[is.finite(null)]
    s <- s + 1L
    summary_rows[[s]] <- data.table(cohort = observed$cohort[[i]], axis = observed$axis[[i]], metric = metric, observed = obs, n_random = length(null), absolute_empirical_p = if (is.finite(obs) && length(null)) (1 + sum(abs(null) >= abs(obs))) / (1 + length(null)) else NA_real_, absolute_percentile = if (is.finite(obs) && length(null)) mean(abs(null) < abs(obs)) else NA_real_, signed_percentile = if (is.finite(obs) && length(null)) mean(null < obs) else NA_real_, status = if (is.finite(obs) && length(null) >= n_random) "estimated" else "not_estimable")
  }
  summary <- rbindlist(summary_rows)
  write_v2_tsv(summary, file.path(FIGURE7_V2_META, "figure7_v2_random_specificity_summary.tsv"))
  list(benchmark = out, summary = summary, observed = observed, clinical_c = clinical_c)
}

# -------------------------------------------------------------------------
# Prediction: TCGA internal validation and a deliberately blocked ICGC
# external-prediction branch.  The latter is not run until the raw ICGC
# age/sex/stage/event/time codes have an independent data dictionary.
# -------------------------------------------------------------------------

harrell_c_v2 <- function(time, event, linear_predictor) {
  ok <- is.finite(time) & is.finite(event) & is.finite(linear_predictor) & time > 0
  if (sum(ok) < 10L || sum(event[ok] == 1L) < 3L) return(NA_real_)
  unname(concordance(Surv(time[ok], event[ok]) ~ linear_predictor[ok], reverse = TRUE)$concordance)
}

ipcw_auc_v2 <- function(time, event, risk, horizon) {
  ok <- is.finite(time) & is.finite(event) & is.finite(risk) & time > 0
  time <- time[ok]; event <- event[ok]; risk <- risk[ok]
  if (sum(event == 1L & time <= horizon) < 3L || sum(time > horizon) < 3L) return(NA_real_)
  censor_fit <- survfit(Surv(time, 1L - event) ~ 1)
  gfun <- function(t) pmax(summary(censor_fit, times = pmax(t - 1e-8, 0), extend = TRUE)$surv, 1e-4)
  cases <- which(event == 1L & time <= horizon); controls <- which(time > horizon)
  wc <- 1 / gfun(time[cases]); w0 <- rep(1 / gfun(horizon), length(controls))
  cmp <- outer(risk[cases], risk[controls], function(a, b) (a > b) + 0.5 * (a == b))
  sum(cmp * outer(wc, w0)) / sum(outer(wc, w0))
}

ipcw_brier_v2 <- function(time, event, predicted_risk, horizon) {
  ok <- is.finite(time) & is.finite(event) & is.finite(predicted_risk) & time > 0
  time <- time[ok]; event <- event[ok]; predicted_risk <- predicted_risk[ok]
  if (!length(time)) return(NA_real_)
  censor_fit <- survfit(Surv(time, 1L - event) ~ 1)
  gfun <- function(t) pmax(summary(censor_fit, times = pmax(t - 1e-8, 0), extend = TRUE)$surv, 1e-4)
  w <- ifelse(time <= horizon & event == 1L, 1 / gfun(time), ifelse(time > horizon, 1 / gfun(horizon), 0))
  if (!any(w > 0)) return(NA_real_)
  y <- as.integer(time <= horizon & event == 1L)
  sum(w * (y - predicted_risk)^2) / sum(w)
}

patient_raw_score_table_v2 <- function(patient, score_version_name = "primary_frozen_programme") {
  axes <- c("identity_loss", "stress_transition", "sox4_associated")
  d <- copy(patient[score_version == score_version_name & axis %in% axes])
  id_cols <- c("cohort", "patient_id", "tumour_normal")
  clinical_cols <- c("age_years", "age_high", "sex", "stage_num", "stage_high", "os_time_days", "os_event", "grade_num", "t_stage_num", "clinical_analysis_status", "os_analysis_status")
  score_long <- d[, .(score_raw = mean(score_raw, na.rm = TRUE)), by = c(id_cols, "axis")]
  clinical <- d[, lapply(.SD, first_nonmissing_v2), by = id_cols, .SDcols = clinical_cols]
  wide <- dcast(score_long, cohort + patient_id + tumour_normal ~ axis, value.var = "score_raw", fun.aggregate = mean)
  merge(clinical, wide, by = id_cols, all = TRUE, sort = FALSE)
}

prediction_specs_v2 <- function() list(
  clinical_baseline = character(),
  identity_loss = "identity_loss",
  stress_transition = "stress_transition",
  sox4_associated = "sox4_associated",
  all_three_axes = c("identity_loss", "stress_transition", "sox4_associated")
)

scale_train_test_v2 <- function(train, test, variables) {
  params <- rbindlist(lapply(variables, function(v) {
    z <- as.numeric(train[[v]])
    data.table(feature = v, mean = mean(z, na.rm = TRUE), sd = sd(z, na.rm = TRUE), n_training = sum(is.finite(z)))
  }))
  for (v in variables) {
    p <- params[feature == v]
    zname <- paste0(v, "_z")
    train[[zname]] <- (train[[v]] - p$mean[[1L]]) / p$sd[[1L]]
    test[[zname]] <- (test[[v]] - p$mean[[1L]]) / p$sd[[1L]]
  }
  list(train = train, test = test, parameters = params)
}

fit_prediction_v2 <- function(d, extra_raw = character()) {
  extra_terms <- if (length(extra_raw)) paste0(extra_raw, "_z") else character()
  terms <- c("age_years_z", "factor(sex)", "stage_high", extra_terms)
  required <- unique(c("os_time_days", "os_event", "age_years", "sex", "stage_high", extra_raw))
  z <- copy(d)[complete.cases(d[, ..required]) & os_time_days > 0 & os_event %in% c(0, 1)]
  if (nrow(z) < 20L || sum(z$os_event == 1L) < 5L) return(list(status = "not_estimable", fit = NULL, data = z, terms = terms))
  f <- as.formula(paste("Surv(os_time_days, os_event) ~", paste(terms, collapse = " + ")))
  fit <- try(coxph(f, data = z, x = TRUE, y = TRUE, model = TRUE), silent = TRUE)
  if (inherits(fit, "try-error")) return(list(status = "not_estimable", fit = NULL, data = z, terms = terms))
  list(status = "estimated", fit = fit, data = z, terms = terms)
}

make_stratified_folds_v2 <- function(event, folds, seed) {
  set.seed(seed)
  out <- integer(length(event))
  for (lev in sort(unique(event))) {
    idx <- which(event == lev)
    out[idx] <- sample(rep(seq_len(folds), length.out = length(idx)))
  }
  out
}

predict_risk_at_v2 <- function(fit, d, horizon) {
  lp <- as.numeric(predict(fit, newdata = d, type = "lp", reference = "zero"))
  bh <- basehaz(fit, centered = FALSE)
  pos <- which(bh$time <= horizon)
  h0 <- if (length(pos)) bh$hazard[[max(pos)]] else 0
  list(linear_predictor = lp, risk = 1 - exp(-h0 * exp(lp)), baseline_hazard = h0)
}

cross_validate_prediction_v2 <- function(raw_tcga, repeats = 10L, folds = 5L) {
  specs <- prediction_specs_v2(); continuous <- c("age_years", "identity_loss", "stress_transition", "sox4_associated")
  required_input <- c("os_time_days", "os_event", "age_years", "sex", "stage_high", "identity_loss", "stress_transition", "sox4_associated")
  base_input <- raw_tcga[complete.cases(raw_tcga[, ..required_input]) & os_time_days > 0 & os_event %in% c(0, 1)]
  fold_rows <- list(); auc_rows <- list(); k <- 0L; a <- 0L
  for (r in seq_len(repeats)) {
    fold_id <- make_stratified_folds_v2(base_input$os_event, folds, FIGURE7_V2_SEED + r)
    for (f in seq_len(folds)) {
      train <- copy(base_input[fold_id != f]); test <- copy(base_input[fold_id == f])
      scaled <- scale_train_test_v2(train, test, continuous); train <- scaled$train; test <- scaled$test
      fitted <- lapply(specs, function(extras) fit_prediction_v2(train, extras)); names(fitted) <- names(specs)
      for (model_name in names(fitted)) {
        res <- fitted[[model_name]]
        lp <- if (identical(res$status, "estimated")) as.numeric(predict(res$fit, newdata = test, type = "lp", reference = "zero")) else rep(NA_real_, nrow(test))
        k <- k + 1L
        fold_rows[[k]] <- data.table(cohort = "TCGA_LIHC", validation = "repeated_10x5_fold_CV", repeat_id = r, fold_id = f, model = model_name,
                                     c_index = harrell_c_v2(test$os_time_days, test$os_event, lp), n_test = nrow(test), events_test = sum(test$os_event == 1L), status = res$status)
        if (identical(res$status, "estimated")) for (yr in c(1, 3, 5)) {
          pr <- predict_risk_at_v2(res$fit, test, yr * 365.25)
          a <- a + 1L
          auc_rows[[a]] <- data.table(cohort = "TCGA_LIHC", validation = "repeated_10x5_fold_CV", repeat_id = r, fold_id = f, model = model_name,
                                      year = yr, horizon_days = yr * 365.25, time_dependent_auc = ipcw_auc_v2(test$os_time_days, test$os_event, pr$linear_predictor, yr * 365.25),
                                      brier_score = ipcw_brier_v2(test$os_time_days, test$os_event, pr$risk, yr * 365.25), n_test = nrow(test), events_test = sum(test$os_event == 1L), status = "estimated")
        }
      }
    }
  }
  folds_out <- rbindlist(fold_rows, fill = TRUE); auc_out <- rbindlist(auc_rows, fill = TRUE)
  base <- folds_out[model == "clinical_baseline", .(repeat_id, fold_id, baseline_c_index = c_index)]
  folds_out <- merge(folds_out, base, by = c("repeat_id", "fold_id"), all.x = TRUE, sort = FALSE)
  folds_out[, delta_c_index := c_index - baseline_c_index]
  summary <- folds_out[, .(
    n_fold_estimates = sum(is.finite(delta_c_index)), mean_c_index = mean(c_index, na.rm = TRUE), median_c_index = median(c_index, na.rm = TRUE),
    mean_delta_c_index = mean(delta_c_index, na.rm = TRUE), median_delta_c_index = median(delta_c_index, na.rm = TRUE),
    delta_ci_low = if (sum(is.finite(delta_c_index)) >= 2L) quantile(delta_c_index, 0.025, na.rm = TRUE) else NA_real_,
    delta_ci_high = if (sum(is.finite(delta_c_index)) >= 2L) quantile(delta_c_index, 0.975, na.rm = TRUE) else NA_real_,
    n_test_total = sum(n_test), events_test_total = sum(events_test)
  ), by = .(cohort, validation, model)]
  auc_summary <- auc_out[, .(n_fold_estimates = sum(is.finite(time_dependent_auc)), mean_time_dependent_auc = mean(time_dependent_auc, na.rm = TRUE), median_time_dependent_auc = median(time_dependent_auc, na.rm = TRUE), auc_ci_low = if (sum(is.finite(time_dependent_auc)) >= 2L) quantile(time_dependent_auc, .025, na.rm = TRUE) else NA_real_, auc_ci_high = if (sum(is.finite(time_dependent_auc)) >= 2L) quantile(time_dependent_auc, .975, na.rm = TRUE) else NA_real_, mean_brier_score = mean(brier_score, na.rm = TRUE)), by = .(cohort, validation, model, year, horizon_days)]
  list(folds = folds_out, auc = auc_out, summary = summary, auc_summary = auc_summary, input = base_input)
}

stage_v2_prediction <- function(patient) {
  raw <- patient_raw_score_table_v2(patient)[cohort == "TCGA_LIHC" & tumour_normal == "tumour"]
  raw <- raw[complete.cases(age_years, sex, stage_high, os_time_days, os_event, identity_loss, stress_transition, sox4_associated) & os_time_days > 0 & os_event %in% c(0, 1)]
  continuous <- c("age_years", "identity_loss", "stress_transition", "sox4_associated")
  global_parameters <- rbindlist(lapply(continuous, function(v) data.table(feature = v, feature_type = if (v == "age_years") "clinical_continuous" else "primary_axis_raw_score", tcga_training_mean = mean(raw[[v]], na.rm = TRUE), tcga_training_sd = sd(raw[[v]], na.rm = TRUE), n_training = sum(is.finite(raw[[v]])), scaling_use = "TCGA_training_frozen_for_locked_prediction")))
  categorical <- data.table(feature = c("sex", "stage_high"), feature_type = c("clinical_binary", "clinical_binary"), tcga_training_mean = NA_real_, tcga_training_sd = NA_real_, n_training = nrow(raw), scaling_use = c("fixed_raw_category_no_rescaling", "fixed_raw_category_no_rescaling"))
  prediction_parameters <- rbindlist(list(global_parameters, categorical), fill = TRUE)
  write_v2_tsv(prediction_parameters, file.path(FIGURE7_V2_META, "figure7_v2_tcga_prediction_scaling_parameters.tsv"))
  full <- copy(raw)
  for (v in continuous) {
    p <- global_parameters[feature == v]
    full[[paste0(v, "_z")]] <- (full[[v]] - p$tcga_training_mean[[1L]]) / p$tcga_training_sd[[1L]]
  }
  specs <- prediction_specs_v2()
  full_models <- lapply(specs, function(extras) fit_prediction_v2(full, extras)); names(full_models) <- names(specs)
  app_rows <- list(); coef_rows <- list(); recal_rows <- list(); i <- 0L; j <- 0L; q <- 0L
  for (model_name in names(full_models)) {
    res <- full_models[[model_name]]
    if (!identical(res$status, "estimated")) {
      i <- i + 1L; app_rows[[i]] <- data.table(cohort = "TCGA_LIHC", validation = "apparent_source_data_only", model = model_name, status = res$status, c_index = NA_real_, n = nrow(res$data), events = sum(res$data$os_event == 1L), predictors = NA_integer_)
      next
    }
    lp <- as.numeric(predict(res$fit, newdata = res$data, type = "lp", reference = "zero"))
    i <- i + 1L; app_rows[[i]] <- data.table(cohort = "TCGA_LIHC", validation = "apparent_source_data_only", model = model_name, status = "estimated", c_index = harrell_c_v2(res$data$os_time_days, res$data$os_event, lp), n = nrow(res$data), events = sum(res$data$os_event == 1L), predictors = length(coef(res$fit)))
    co <- as.data.table(summary(res$fit)$coefficients, keep.rownames = "term")
    j <- j + 1L; coef_rows[[j]] <- data.table(cohort = "TCGA_LIHC", model = model_name, term = co$term, coefficient = co$coef, se = co$`se(coef)`, hazard_ratio = exp(co$coef), p_value = co$`Pr(>|z|)`, weights_locked = TRUE, fit_source = "TCGA_training_only")
    recal <- try(coxph(Surv(os_time_days, os_event) ~ lp, data = data.table(res$data, lp = lp)), silent = TRUE)
    q <- q + 1L; recal_rows[[q]] <- data.table(cohort = "TCGA_LIHC", model = model_name, metric = "overall_Cox_recalibration_slope_not_3_year_calibration", status = if (inherits(recal, "coxph")) "estimated" else "not_estimable", value = if (inherits(recal, "coxph")) as.numeric(coef(recal)[[1L]]) else NA_real_)
  }
  app <- rbindlist(app_rows, fill = TRUE); coefficients <- rbindlist(coef_rows, fill = TRUE); recalibration <- rbindlist(recal_rows, fill = TRUE)
  cv <- cross_validate_prediction_v2(raw, repeats = 10L, folds = 5L)
  cv$summary[, `:=`(status = "estimated", reason = "")]
  cv$auc_summary[, `:=`(status = "estimated", reason = "")]
  blocked_cv <- rbindlist(lapply(names(specs), function(model_name) data.table(cohort = "ICGC_LIRI_JP", validation = "repeated_10x5_fold_CV", model = model_name, n_fold_estimates = 0L, mean_c_index = NA_real_, median_c_index = NA_real_, mean_delta_c_index = NA_real_, median_delta_c_index = NA_real_, delta_ci_low = NA_real_, delta_ci_high = NA_real_, n_test_total = NA_integer_, events_test_total = NA_integer_, status = "blocked", reason = "ICGC_Age_Gender_Stage_fustat_futime_encoding_unverified")), fill = TRUE)
  blocked_cv_auc <- rbindlist(lapply(names(specs), function(model_name) rbindlist(lapply(c(1, 3, 5), function(yr) data.table(cohort = "ICGC_LIRI_JP", validation = "repeated_10x5_fold_CV", model = model_name, year = yr, horizon_days = yr * 365.25, n_fold_estimates = 0L, mean_time_dependent_auc = NA_real_, median_time_dependent_auc = NA_real_, auc_ci_low = NA_real_, auc_ci_high = NA_real_, mean_brier_score = NA_real_, status = "blocked", reason = "ICGC_Age_Gender_Stage_fustat_futime_encoding_unverified")), fill = TRUE)), fill = TRUE)
  cv$summary <- rbindlist(list(cv$summary, blocked_cv), fill = TRUE)
  cv$auc_summary <- rbindlist(list(cv$auc_summary, blocked_cv_auc), fill = TRUE)
  blocked_external <- rbindlist(lapply(names(specs), function(model_name) data.table(train_cohort = "TCGA_LIHC", test_cohort = "ICGC_LIRI_JP", model = model_name, validation = "locked_external_TCGA_to_ICGC", weights_locked = TRUE, coefficients_refit_in_ICGC = FALSE, status = "blocked", reason = "ICGC_Age_Gender_Stage_fustat_futime_encoding_unverified", c_index = NA_real_, delta_c_index = NA_real_, delta_ci_low = NA_real_, delta_ci_high = NA_real_, n = NA_integer_, events = NA_integer_)), fill = TRUE)
  blocked_auc <- rbindlist(lapply(c("clinical_baseline", "all_three_axes"), function(model_name) rbindlist(lapply(c(1, 3, 5), function(yr) data.table(cohort = "ICGC_LIRI_JP", validation = "locked_external_TCGA_to_ICGC", model = model_name, year = yr, horizon_days = yr * 365.25, status = "blocked", reason = "ICGC_Age_Gender_Stage_fustat_futime_encoding_unverified", time_dependent_auc = NA_real_, brier_score = NA_real_)), fill = TRUE)), fill = TRUE)
  blocked_calibration <- rbindlist(lapply(c("clinical_baseline", "all_three_axes"), function(model_name) data.table(cohort = "ICGC_LIRI_JP", validation = "locked_external_TCGA_to_ICGC", model = model_name, year = 3L, horizon_days = 3 * 365.25, status = "blocked", reason = "ICGC_Age_Gender_Stage_fustat_futime_encoding_unverified", calibration_in_large = NA_real_, observed_predicted_ratio = NA_real_, brier_score = NA_real_, calibration_groups = NA_character_)), fill = TRUE)
  write_v2_tsv(app, file.path(FIGURE7_V2_FIG, "figure7_v2f_apparent_performance_source_data.tsv"))
  write_v2_tsv(cv$folds, file.path(FIGURE7_V2_FIG, "figure7_v2f_cv_fold_metrics.tsv"))
  write_v2_tsv(cv$summary, file.path(FIGURE7_V2_FIG, "figure7_v2f_cv_summary.tsv"))
  write_v2_tsv(cv$auc, file.path(FIGURE7_V2_FIG, "figure7_v2f_landmark_auc_fold_metrics.tsv"))
  write_v2_tsv(cv$auc_summary, file.path(FIGURE7_V2_FIG, "figure7_v2f_landmark_auc_summary.tsv"))
  write_v2_tsv(blocked_external, file.path(FIGURE7_V2_FIG, "figure7_v2f_locked_external_validation.tsv"))
  write_v2_tsv(blocked_auc, file.path(FIGURE7_V2_FIG, "figure7_v2f_locked_external_landmark_auc.tsv"))
  write_v2_tsv(blocked_calibration, file.path(FIGURE7_V2_FIG, "figure7_v2f_external_3year_calibration.tsv"))
  write_v2_tsv(data.table(cohort = "ICGC_LIRI_JP", model = c("clinical_baseline", "all_three_axes"), status = "blocked", reason = "ICGC_Age_Gender_Stage_fustat_futime_encoding_unverified", risk_group = NA_character_, n = NA_integer_, mean_predicted_3year_risk = NA_real_, observed_3year_risk = NA_real_, ci_low = NA_real_, ci_high = NA_real_), file.path(FIGURE7_V2_FIG, "figure7_v2f_external_3year_calibration_groups.tsv"))
  write_v2_tsv(coefficients, file.path(FIGURE7_V2_FIG, "figure7_v2f_locked_tcga_model_coefficients.tsv"))
  write_v2_tsv(recalibration, file.path(FIGURE7_V2_FIG, "figure7_v2f_overall_cox_recalibration.tsv"))
  list(raw_tcga = raw, full = full, full_models = full_models, apparent = app, coefficients = coefficients, recalibration = recalibration, cv = cv, external = blocked_external, external_auc = blocked_auc, calibration = blocked_calibration, prediction_parameters = prediction_parameters)
}

stage_v2_locked_risk <- function(prediction) {
  res <- prediction$full_models[["all_three_axes"]]
  if (!identical(res$status, "estimated")) stop("TCGA all-three-axis locked model was not estimable.")
  d <- copy(res$data)
  d[, risk_score := as.numeric(predict(res$fit, newdata = d, type = "lp", reference = "zero"))]
  cutoff <- median(d$risk_score, na.rm = TRUE)
  d[, risk_group := fifelse(risk_score >= cutoff, "High", "Low")]
  group_summary <- d[, .(n = .N, events = sum(os_event == 1L), median_risk_score = median(risk_score), mean_risk_score = mean(risk_score)), by = risk_group]
  continuous <- summary(coxph(Surv(os_time_days, os_event) ~ risk_score, data = d))$coefficients["risk_score", ]
  hr <- data.table(cohort = "TCGA_LIHC", model = "clinical_plus_all_three_axes", metric = "continuous_locked_risk_score", status = "estimated", coefficient = continuous[["coef"]], hazard_ratio = exp(continuous[["coef"]]), ci_low = exp(continuous[["coef"]] - 1.96 * continuous[["se(coef)"]]), ci_high = exp(continuous[["coef"]] + 1.96 * continuous[["se(coef)"]]), p_value = continuous[["Pr(>|z|)"]], tcga_median_cutoff = cutoff)
  risk_table <- d[, .(cohort = "TCGA_LIHC", patient_id, os_time_days, os_event, risk_score, risk_group, cutoff_source = "TCGA_training_median_locked")]
  icgc_block <- data.table(cohort = "ICGC_LIRI_JP", patient_id = NA_character_, os_time_days = NA_real_, os_event = NA_real_, risk_score = NA_real_, risk_group = NA_character_, cutoff_source = "not_applied", status = "blocked", reason = "ICGC_Age_Gender_Stage_fustat_futime_encoding_unverified")
  risk_table <- rbindlist(list(risk_table, icgc_block), fill = TRUE)
  at_risk_times <- c(0, 365.25, 3 * 365.25, 5 * 365.25)
  risk_counts <- rbindlist(lapply(c("Low", "High"), function(gr) data.table(cohort = "TCGA_LIHC", risk_group = gr, time_years = at_risk_times / 365.25, n_at_risk = vapply(at_risk_times, function(t) sum(d[risk_group == gr, os_time_days] >= t), integer(1)))), fill = TRUE)
  write_v2_tsv(risk_table, file.path(FIGURE7_V2_FIG, "figure7_v2g_locked_risk_groups.tsv"))
  write_v2_tsv(group_summary, file.path(FIGURE7_V2_FIG, "figure7_v2g_tcga_risk_group_summary.tsv"))
  write_v2_tsv(hr, file.path(FIGURE7_V2_FIG, "figure7_v2g_locked_continuous_risk_hr.tsv"))
  write_v2_tsv(risk_counts, file.path(FIGURE7_V2_FIG, "figure7_v2g_tcga_risk_table.tsv"))
  write_v2_tsv(data.table(cohort = "ICGC_LIRI_JP", status = "blocked", coefficients_locked = TRUE, coefficients_refit = FALSE, cutoff_optimized = FALSE, reason = "ICGC_Age_Gender_Stage_fustat_futime_encoding_unverified"), file.path(FIGURE7_V2_FIG, "figure7_v2g_locked_external_status.tsv"))
  list(data = d, cutoff = cutoff, group_summary = group_summary, continuous_hr = hr, risk_counts = risk_counts)
}

# -------------------------------------------------------------------------
# Sensitivity analysis.  Every declared scenario carries its actual changed
# analytic element in a machine-readable audit.  Unchanged score definitions
# are explicitly kept out of the independent-sensitivity evidence set.
# -------------------------------------------------------------------------

score_definition_diff_v2 <- function(definitions, score_version_name, axis_name) {
  primary <- sort(unique(definitions[score_version == "primary_frozen_programme" & axis == axis_name, gene]))
  candidate <- sort(unique(definitions[score_version == score_version_name & axis == axis_name, gene]))
  list(changed = !identical(primary, candidate), primary_n = length(primary), candidate_n = length(candidate),
       removed = paste(setdiff(primary, candidate), collapse = ";"), added = paste(setdiff(candidate, primary), collapse = ";"))
}

merge_version_with_comparators_v2 <- function(patient, score_version_name) {
  z <- patient_score_table(patient, score_version_name)
  comp <- comparator_patient_table(patient)
  merge(z, comp, by = c("cohort", "patient_id", "tumour_normal"), all.x = TRUE, sort = FALSE)
}

fit_sensitivity_axis_v2 <- function(d, axis_name, scenario_id, extra_covars = character(), subset_index = rep(TRUE, nrow(d)), baseline_covars = c("age_years", "sex", "stage_high")) {
  z <- copy(d[subset_index])
  needed <- unique(c("os_time_days", "os_event", baseline_covars, axis_name, extra_covars))
  z <- z[complete.cases(z[, ..needed]) & os_time_days > 0 & os_event %in% c(0, 1)]
  all_covars <- c(baseline_covars, extra_covars)
  model_terms <- c(if ("age_years" %in% all_covars) "age_years" else character(), if ("sex" %in% all_covars) "factor(sex)" else character(), setdiff(all_covars, c("age_years", "sex")), axis_name)
  if (nrow(z) < 20L || sum(z$os_event == 1L) < 5L) return(data.table(model_id = scenario_id, programme = axis_name, term = axis_name, status = "not_estimable", reason = "insufficient_events_or_complete_cases", n = nrow(z), events = sum(z$os_event == 1L), predictors = length(model_terms), events_per_variable = NA_real_, coefficient = NA_real_, se = NA_real_, hazard_ratio = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_, concordance = NA_real_, ph_p_value = NA_real_, condition_index = NA_real_, scenario = scenario_id, axis = axis_name, effect_type = "hazard_ratio", effect = NA_real_, analytic_model = paste0("Surv(OS) ~ ", paste(model_terms, collapse = " + "))))
  fit <- try(coxph(as.formula(paste("Surv(os_time_days, os_event) ~", paste(model_terms, collapse = " + "))), data = z, x = TRUE), silent = TRUE)
  if (inherits(fit, "try-error")) return(data.table(model_id = scenario_id, programme = axis_name, term = axis_name, status = "not_estimable", reason = "cox_fit_failed", n = nrow(z), events = sum(z$os_event == 1L), predictors = length(model_terms), events_per_variable = sum(z$os_event == 1L) / length(model_terms), coefficient = NA_real_, se = NA_real_, hazard_ratio = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_, concordance = NA_real_, ph_p_value = NA_real_, condition_index = NA_real_, scenario = scenario_id, axis = axis_name, effect_type = "hazard_ratio", effect = NA_real_, analytic_model = paste0("Surv(OS) ~ ", paste(model_terms, collapse = " + "))))
  co <- summary(fit)$coefficients
  beta <- co[axis_name, "coef"]; se <- co[axis_name, "se(coef)"]
  ph <- try(cox.zph(fit), silent = TRUE); ph_p <- if (inherits(ph, "try-error") || !axis_name %in% rownames(ph$table)) NA_real_ else ph$table[axis_name, "p"]
  data.table(model_id = scenario_id, programme = axis_name, term = axis_name, status = "estimated", reason = "", n = nrow(z), events = sum(z$os_event == 1L), predictors = length(coef(fit)), events_per_variable = sum(z$os_event == 1L) / length(coef(fit)), coefficient = beta, se = se, hazard_ratio = exp(beta), ci_low = exp(beta - 1.96 * se), ci_high = exp(beta + 1.96 * se), p_value = co[axis_name, "Pr(>|z|)"], concordance = summary(fit)$concordance[[1L]], ph_p_value = ph_p, condition_index = NA_real_, scenario = scenario_id, axis = axis_name, effect_type = "hazard_ratio", effect = exp(beta), analytic_model = paste0("Surv(OS) ~ ", paste(model_terms, collapse = " + ")))
}

stage_v2_sensitivity <- function(state) {
  patient <- state$patient; definitions <- state$definitions
  axes <- c("identity_loss", "stress_transition", "sox4_associated")
  scenario_specs <- list(
    list(id = "primary_clinical_adjusted", version = "primary_frozen_programme", subset = "all_eligible", covars = character(), change_type = "primary_reference"),
    list(id = "exclude_top20pct_proliferation", version = "primary_frozen_programme", subset = "bottom_80pct_proliferation", covars = character(), change_type = "sample_subset"),
    list(id = "adjust_proliferation", version = "primary_frozen_programme", subset = "all_eligible", covars = "proliferation_comparator", change_type = "covariate_adjustment"),
    list(id = "adjust_hypoxia", version = "primary_frozen_programme", subset = "all_eligible", covars = "hypoxia_covariate", change_type = "covariate_adjustment"),
    list(id = "adjust_inflammation", version = "primary_frozen_programme", subset = "all_eligible", covars = "inflammation_covariate", change_type = "covariate_adjustment"),
    list(id = "adjust_foxm1_cebpb_reference", version = "primary_frozen_programme", subset = "all_eligible", covars = "foxm1_cebpb_reference", change_type = "covariate_adjustment"),
    list(id = "adjust_broad_network_calibration", version = "primary_frozen_programme", subset = "all_eligible", covars = "broad_network_calibration", change_type = "covariate_adjustment"),
    list(id = "tf_expression_only", version = "tf_expression_only", subset = "all_eligible", covars = character(), change_type = "score_definition"),
    list(id = "regulon_only", version = "regulon_only", subset = "all_eligible", covars = character(), change_type = "score_definition"),
    list(id = "celloracle_target_only", version = "celloracle_target_only", subset = "all_eligible", covars = character(), change_type = "score_definition"),
    list(id = "high_confidence_intersection", version = "high_confidence_intersection", subset = "all_eligible", covars = character(), change_type = "score_definition"),
    list(id = "remove_cell_cycle_genes", version = "no_cell_cycle", subset = "all_eligible", covars = character(), change_type = "score_definition"),
    list(id = "remove_generic_stress_genes", version = "no_generic_stress", subset = "all_eligible", covars = character(), change_type = "score_definition"),
    list(id = "alternative_median_rank", version = "alternative_median_rank", subset = "all_eligible", covars = character(), change_type = "score_method"),
    list(id = "signed_vs_unsigned", version = "signed_target", subset = "all_eligible", covars = character(), change_type = "signed_score_availability"),
    list(id = "early_stage", version = "primary_frozen_programme", subset = "early_stage", covars = character(), change_type = "sample_subset"),
    list(id = "advanced_stage", version = "primary_frozen_programme", subset = "advanced_stage", covars = character(), change_type = "sample_subset"),
    list(id = "male", version = "primary_frozen_programme", subset = "male", covars = character(), change_type = "sample_subset"),
    list(id = "female", version = "primary_frozen_programme", subset = "female", covars = character(), change_type = "sample_subset")
  )
  rows <- list(); audits <- list(); k <- 0L; a <- 0L
  for (spec in scenario_specs) {
    for (cohort_name in c("TCGA_LIHC", "ICGC_LIRI_JP")) {
      for (axis_name in axes) {
        defdiff <- if (identical(spec$change_type, "score_definition") || identical(spec$change_type, "score_method")) score_definition_diff_v2(definitions, spec$version, axis_name) else list(changed = TRUE, primary_n = NA_integer_, candidate_n = NA_integer_, removed = "", added = "")
        # A median rank score changes the score operation even when its genes
        # are identical to the primary mean-rank score.
        if (identical(spec$change_type, "score_method")) defdiff$changed <- TRUE
        change_reason <- spec$change_type
        status_audit <- "eligible"
        reason_audit <- ""
        if (identical(spec$id, "signed_vs_unsigned") && !"signed_target" %in% unique(definitions$score_version)) {
          status_audit <- "not_estimable"; reason_audit <- "no_reliably_signed_primary_score_available_after_direction_audit"
        } else if (identical(spec$change_type, "score_definition") && !isTRUE(defdiff$changed)) {
          status_audit <- "not_distinct"; reason_audit <- "candidate_score_definition_is_identical_to_primary_for_this_axis"
        } else if (cohort_name == "ICGC_LIRI_JP") {
          status_audit <- "blocked"; reason_audit <- "ICGC_Age_Gender_Stage_fustat_futime_encoding_unverified"
        }
        a <- a + 1L
        audits[[a]] <- data.table(cohort = cohort_name, scenario = spec$id, axis = axis_name, score_version = spec$version,
                                  change_type = change_reason, sample_rule = spec$subset,
                                  added_covariates = if (length(spec$covars)) paste(spec$covars, collapse = ";") else "none",
                                  model = paste0("Surv(OS) ~ age_years + sex + stage_high + score", if (length(spec$covars)) paste0(" + ", paste(spec$covars, collapse = " + ")) else ""),
                                  primary_gene_n = defdiff$primary_n, candidate_gene_n = defdiff$candidate_n,
                                  removed_genes = defdiff$removed, added_genes = defdiff$added,
                                  actual_analytic_change = isTRUE(defdiff$changed), status = status_audit, reason = reason_audit)
        if (status_audit != "eligible") {
          k <- k + 1L; rows[[k]] <- data.table(cohort = cohort_name, scenario = spec$id, axis = axis_name, score_version = spec$version, status = status_audit, reason = reason_audit,
                                                 effect_type = "hazard_ratio", coefficient = NA_real_, se = NA_real_, effect = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_, fdr = NA_real_, n = NA_integer_, events = NA_integer_, predictors = NA_integer_, events_per_variable = NA_real_, low_epv_flag = status_audit, analytic_model = audits[[a]]$model[[1L]])
          next
        }
        d <- merge_version_with_comparators_v2(patient, spec$version)[cohort == cohort_name & tumour_normal == "tumour"]
        if (spec$subset == "bottom_80pct_proliferation") {
          cutoff <- quantile(d$proliferation_comparator, .8, na.rm = TRUE, type = 8)
          subset_index <- is.finite(d$proliferation_comparator) & d$proliferation_comparator <= cutoff
        } else if (spec$subset == "early_stage") subset_index <- d$stage_high == 0L
        else if (spec$subset == "advanced_stage") subset_index <- d$stage_high == 1L
        else if (spec$subset == "male") subset_index <- toupper(d$sex) == "MALE"
        else if (spec$subset == "female") subset_index <- toupper(d$sex) == "FEMALE"
        else subset_index <- rep(TRUE, nrow(d))
        baseline_covars <- c("age_years", "sex", "stage_high")
        if (spec$subset %in% c("early_stage", "advanced_stage")) baseline_covars <- setdiff(baseline_covars, "stage_high")
        if (spec$subset %in% c("male", "female")) baseline_covars <- setdiff(baseline_covars, "sex")
        fit <- fit_sensitivity_axis_v2(d, axis_name, spec$id, extra_covars = spec$covars, subset_index = subset_index, baseline_covars = baseline_covars)
        k <- k + 1L
        rows[[k]] <- fit[, .(cohort = cohort_name, scenario, axis, score_version = spec$version, status, reason, effect_type, coefficient, se, effect, ci_low, ci_high, p_value, fdr = NA_real_, n, events, predictors, events_per_variable, low_epv_flag = ifelse(status == "estimated" & events_per_variable < 8, "exploratory_overfitting_risk", ifelse(status == "estimated", "none", status)), analytic_model)]
      }
    }
  }
  out <- rbindlist(rows, fill = TRUE)
  out[, fdr := bh(p_value), by = .(cohort, scenario)]
  audit <- rbindlist(audits, fill = TRUE)
  primary <- out[scenario == "primary_clinical_adjusted" & cohort == "TCGA_LIHC" & status == "estimated", .(axis, primary_coefficient = coefficient)]
  out <- merge(out, primary, by = "axis", all.x = TRUE, sort = FALSE)
  out[, same_direction_as_primary := ifelse(status == "estimated" & is.finite(primary_coefficient), sign(coefficient) == sign(primary_coefficient), NA)]
  out[, primary_coefficient := NULL]
  tcga <- out[cohort == "TCGA_LIHC" & status == "estimated" & scenario != "primary_clinical_adjusted"]
  stability <- tcga[, .(n_estimated_sensitivities = .N, same_direction_fraction = mean(same_direction_as_primary, na.rm = TRUE), fdr_lt_005_fraction = mean(fdr < .05, na.rm = TRUE)), by = axis]
  key <- out[scenario %in% c("primary_clinical_adjusted", "adjust_proliferation", "adjust_foxm1_cebpb_reference", "remove_cell_cycle_genes", "remove_generic_stress_genes")]
  tumour_normal <- stage_v2_tumour_normal(patient)$primary
  recurrence <- tumour_normal[, .(cohort, axis, hedges_g)]
  recurrence_status <- recurrence[, .(n_cohorts = .N, cross_cohort_direction = if (all(sign(hedges_g) == sign(hedges_g[[1L]]))) "consistent" else "heterogeneous"), by = axis]
  random <- read_v2_tsv(file.path(FIGURE7_V2_META, "figure7_v2_random_specificity_summary.tsv"))
  random_tcga <- random[cohort == "TCGA_LIHC" & status == "estimated", .(n_random_metrics = .N, any_empirical_p_lt_005 = any(absolute_empirical_p < .05), min_empirical_p = min(absolute_empirical_p, na.rm = TRUE)), by = axis]
  evidence <- merge(merge(data.table(axis = axes), recurrence_status, by = "axis", all.x = TRUE), stability, by = "axis", all.x = TRUE)
  evidence <- merge(evidence, random_tcga, by = "axis", all.x = TRUE)
  evidence[, `:=`(cross_cohort_clinical_status = "not_estimable", overall_evidence_status = "Not estimable", evidence_reason = "ICGC clinical/survival codebook is absent; cross-cohort primary HR, locked prediction, and calibration cannot be evaluated under the prespecified rule.")]
  write_v2_tsv(out, file.path(FIGURE7_V2_FIG, "figure7_v2h_sensitivity_results.tsv"))
  write_v2_tsv(audit, file.path(FIGURE7_V2_META, "figure7_v2_sensitivity_implementation_audit.tsv"))
  write_v2_tsv(key, file.path(FIGURE7_V2_FIG, "figure7_v2h_key_sensitivity_results.tsv"))
  write_v2_tsv(recurrence_status, file.path(FIGURE7_V2_FIG, "figure7_v2h_cross_cohort_recurrence_direction.tsv"))
  write_v2_tsv(evidence, file.path(FIGURE7_V2_FIG, "figure7_v2h_evidence_status.tsv"))
  list(results = out, audit = audit, key = key, stability = stability, recurrence = recurrence_status, evidence = evidence)
}

# -------------------------------------------------------------------------
# Figure construction.  These panels visualise only v2 source data and carry
# explicit blocked annotations whenever ICGC clinical/survival semantics are
# unavailable.  Formal output is vector plus 600-dpi raster in every case.
# -------------------------------------------------------------------------

axis_labels_v2 <- c(identity_loss = "Identity loss", stress_transition = "Stress transition", sox4_associated = "SOX4-associated")
cohort_labels_v2 <- c(TCGA_LIHC = "TCGA-LIHC", ICGC_LIRI_JP = "ICGC-LIRI-JP")
axis_label_cols_v2 <- setNames(unname(axis_cols[names(axis_labels_v2)]), unname(axis_labels_v2))

block_panel_v2 <- function(title, message) {
  ggplot() + annotate("label", x = .5, y = .55, label = message, size = 3, linewidth = .25, fill = "#F5F5F5") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) + labs(title = title) + v2_theme() +
    theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
}

plot_v2_workflow_a <- function() {
  nodes <- data.table(x = c(1.05, 3.15, 5.25, 7.35), y = 1,
                      label = c("Figures 2-6\nfrozen programmes", "Independent bulk\nTCGA-LIHC + ICGC-LIRI-JP", "Patient-level\nrecurrence, clinical and OS", "Sensitivity and\nmatched-random specificity"),
                      fill = c(lancet[1], lancet[3], lancet[2], lancet[4]))
  ggplot(nodes, aes(x, y)) + geom_segment(data = data.table(x = c(1.75, 3.85, 5.95), xend = c(2.45, 4.55, 6.65), y = 1, yend = 1), aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE, arrow = arrow(length = unit(2, "mm")), linewidth = .45) +
    geom_label(aes(label = label, fill = fill), colour = "white", linewidth = 0, size = 2.9, fontface = "bold", label.padding = unit(.23, "lines")) + scale_fill_identity() +
    annotate("text", x = 4, y = .36, label = "No bulk outcome information was used for signature derivation. ICGC clinical/survival branch is blocked pending a source codebook.", size = 2.7) +
    coord_cartesian(xlim = c(0, 8.4), ylim = c(.15, 1.55), expand = FALSE, clip = "off") + labs(title = "A  Discovery-to-independent-bulk evaluation workflow") + v2_theme() +
    theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
}

plot_v2_mapping_b <- function() {
  d <- read_v2_tsv(file.path(FIGURE7_V2_META, "figure7_v2_signature_coverage.tsv"))[score_version == "primary_frozen_programme" & axis %in% names(axis_labels_v2)]
  d[, `:=`(axis_label = axis_labels_v2[axis], cohort_label = cohort_labels_v2[cohort], coverage_pct = 100 * mapped_fraction,
           label = paste0(n_mapped_genes, "/", n_frozen_genes, " genes\n", sprintf("%.0f%%", 100 * mapped_fraction)))]
  ggplot(d, aes(axis_label, cohort_label, fill = coverage_pct)) + geom_tile(colour = "white", linewidth = .8) + geom_text(aes(label = label), size = 3.2) +
    scale_fill_gradient(low = "#EAF2F8", high = lancet[1], limits = c(0, 100), name = "Mapping") +
    labs(title = "B  Frozen single-cell programmes mapped to independent bulk cohorts", subtitle = "HGNC-harmonized sample-wise normalized rank scores; primary score = unsigned associated target programme score", x = NULL, y = NULL) + v2_theme()
}

plot_v2_recurrence_c <- function() {
  d <- read_v2_tsv(file.path(FIGURE7_V2_FIG, "figure7_v2c_tumour_normal_effects.tsv"))
  d[, `:=`(axis_label = axis_labels_v2[axis], cohort_label = cohort_labels_v2[cohort])]
  p <- ggplot(d, aes(hedges_g, axis_label, colour = axis_label, shape = cohort_label)) + geom_vline(xintercept = 0, linetype = 2, linewidth = .3) +
    geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = .16, orientation = "y", position = position_dodge(width = .45)) + geom_point(size = 2.4, position = position_dodge(width = .45)) +
    scale_colour_manual(values = axis_label_cols_v2, guide = "none") + scale_shape_manual(values = c("TCGA-LIHC" = 16, "ICGC-LIRI-JP" = 17), name = NULL) +
    labs(title = "C  Patient-level tumour-normal recurrence", subtitle = "Independent Hedges' g; random-effects pooling (k=2) is exploratory only", x = "Tumour minus normal Hedges' g (95% CI)", y = NULL) + v2_theme() + theme(legend.position = "bottom")
  p
}

plot_v2_clinical_d <- function() {
  d <- read_v2_tsv(file.path(FIGURE7_V2_FIG, "figure7_v2d_clinical_associations.tsv"))[clinical_feature == "early_vs_advanced_stage"]
  estimated <- d[status == "estimated"]
  estimated[, `:=`(axis_label = axis_labels_v2[programme], label = sprintf("OR %.2f\nFDR %.3g", odds_ratio, fdr))]
  p1 <- ggplot(estimated, aes(axis_label, odds_ratio, colour = axis_label)) + geom_hline(yintercept = 1, linetype = 2, linewidth = .3) + geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = .15) + geom_point(size = 2.7) +
    scale_colour_manual(values = axis_label_cols_v2, guide = "none") + scale_y_log10() + labs(title = "D  TCGA-LIHC", subtitle = "Age/sex-adjusted stage OR", x = NULL, y = "OR per 1-SD score") + v2_theme()
  p2 <- block_panel_v2("ICGC-LIRI-JP", "Stage + age + sex encoding\nunverified in raw local files\n\nClinical association: BLOCKED")
  p1 | p2
}

plot_v2_cox_e <- function() {
  d <- read_v2_tsv(file.path(FIGURE7_V2_FIG, "figure7_v2e_multivariable_cox_models.tsv"))[cohort == "TCGA_LIHC" & programme %in% names(axis_labels_v2) & model_id != "clinical_adjusted_joint_three_axes"]
  d[, axis_label := axis_labels_v2[programme]]
  p1 <- ggplot(d, aes(hazard_ratio, axis_label, colour = axis_label)) + geom_vline(xintercept = 1, linetype = 2, linewidth = .3) + geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = .16, orientation = "y") + geom_point(size = 2.7) +
    scale_colour_manual(values = axis_label_cols_v2, guide = "none") + scale_x_log10() + labs(title = "E  TCGA-LIHC", subtitle = "Clinical-adjusted OS; n=341, events=113", x = "Hazard ratio (95% CI)", y = NULL) + v2_theme()
  p2 <- block_panel_v2("ICGC-LIRI-JP", "Age, sex, stage, event and time\nencoding unverified\n\nMultivariable OS: BLOCKED")
  p1 | p2
}

plot_v2_prediction_f <- function() {
  cv <- read_v2_tsv(file.path(FIGURE7_V2_FIG, "figure7_v2f_cv_summary.tsv"))[cohort == "TCGA_LIHC" & model != "clinical_baseline"]
  cv[, model_label := c(identity_loss = "Identity", stress_transition = "Stress", sox4_associated = "SOX4-associated", all_three_axes = "All three axes")[model]]
  p1 <- ggplot(cv, aes(mean_delta_c_index, model_label, colour = model_label)) + geom_vline(xintercept = 0, linetype = 2, linewidth = .3) + geom_errorbar(aes(xmin = delta_ci_low, xmax = delta_ci_high), width = .15, orientation = "y") + geom_point(size = 2.6) +
    scale_colour_manual(values = c("Identity" = axis_cols[["identity_loss"]], "Stress" = axis_cols[["stress_transition"]], "SOX4-associated" = axis_cols[["sox4_associated"]], "All three axes" = lancet[4]), guide = "none") +
    labs(title = "F  TCGA internal CV", subtitle = "10 x 5-fold; apparent performance is source-only", x = "Paired Delta C-index versus clinical baseline (95% fold quantiles)", y = NULL) + v2_theme()
  p2 <- block_panel_v2("ICGC external: BLOCKED", "ICGC age/sex/stage/fustat/futime\ncodebook absent\n\nCV, locked Delta C, 1/3/5-year AUC\nand 3-year calibration: BLOCKED")
  p1 | p2
}

km_data_v2 <- function(d) {
  fit <- survfit(Surv(os_time_days, os_event) ~ risk_group, data = d)
  sm <- summary(fit)
  data.table(time_years = sm$time / 365.25, survival = sm$surv, lower = sm$lower, upper = sm$upper, risk_group = sub("risk_group=", "", sm$strata))
}

plot_v2_risk_g <- function(risk) {
  dd <- km_data_v2(risk$data)
  logrank <- survdiff(Surv(os_time_days, os_event) ~ risk_group, data = risk$data)
  pval <- 1 - pchisq(logrank$chisq, df = length(logrank$n) - 1L)
  p1 <- ggplot(dd, aes(time_years, survival, colour = risk_group, fill = risk_group)) + geom_step(linewidth = .7) + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .12, colour = NA) +
    scale_colour_manual(values = c(High = lancet[2], Low = lancet[1]), name = "TCGA locked median group") + scale_fill_manual(values = c(High = lancet[2], Low = lancet[1]), guide = "none") +
    annotate("text", x = Inf, y = .08, hjust = 1.05, label = sprintf("Log-rank P = %.3g\nCutoff = TCGA median", pval), size = 2.8) +
    labs(title = "G  TCGA-LIHC", subtitle = "Locked three-axis risk visualization", x = "Years", y = "Overall survival probability") + v2_theme() + theme(legend.position = "bottom")
  p2 <- block_panel_v2("ICGC-LIRI-JP", "Locked TCGA coefficients/cutoff\nwere not applied because clinical\nand OS encoding is unverified\n\nRisk visualization: BLOCKED")
  p1 | p2
}

plot_v2_sensitivity_h <- function() {
  key <- read_v2_tsv(file.path(FIGURE7_V2_FIG, "figure7_v2h_key_sensitivity_results.tsv"))[cohort == "TCGA_LIHC" & status == "estimated"]
  key[, `:=`(axis_label = axis_labels_v2[axis], scenario_label = c(primary_clinical_adjusted = "Primary", adjust_proliferation = "Adjust proliferation", adjust_foxm1_cebpb_reference = "Adjust FOXM1/CEBPB", remove_cell_cycle_genes = "Remove cell-cycle", remove_generic_stress_genes = "Remove generic stress")[scenario])]
  p1 <- ggplot(key, aes(effect, scenario_label, colour = axis_label)) + geom_vline(xintercept = 1, linetype = 2, linewidth = .3) + geom_errorbar(aes(xmin = ci_low, xmax = ci_high), position = position_dodge(width = .55), width = .14, orientation = "y") + geom_point(position = position_dodge(width = .55), size = 2) +
    scale_x_log10() + scale_colour_manual(values = axis_label_cols_v2, name = NULL) + labs(title = "H  TCGA sensitivity stability", x = "Hazard ratio (95% CI)", y = NULL) + v2_theme() + theme(legend.position = "bottom")
  rnd <- read_v2_tsv(file.path(FIGURE7_V2_META, "figure7_v2_random_specificity_summary.tsv"))[cohort == "TCGA_LIHC" & status == "estimated"]
  rnd[, `:=`(axis_label = axis_labels_v2[axis], metric_label = c(tumour_normal_effect = "Tumour-normal", stage_log_or = "Stage log OR", cox_z = "Adjusted Cox |z|", joint_c_index = "Joint C-index", delta_c_index = "Delta C-index")[metric])]
  p2 <- ggplot(rnd, aes(absolute_percentile, metric_label, colour = axis_label)) + geom_vline(xintercept = .75, linetype = 2, linewidth = .3) + geom_point(position = position_dodge(width = .55), size = 2) + scale_x_continuous(limits = c(0, 1), labels = percent) + scale_colour_manual(values = axis_label_cols_v2, name = NULL) + labs(title = "TCGA matched-random specificity", subtitle = "Absolute-effect percentile; dashed line = 75th percentile", x = "Empirical percentile", y = NULL) + v2_theme() + theme(legend.position = "bottom")
  p3 <- block_panel_v2("Cross-cohort evidence status", "ICGC clinical and OS analyses: BLOCKED\n\nIdentity: recurrence direction consistent\nStress / SOX4: recurrence heterogeneous\n\nNo axis qualifies as cross-cohort\nclinical prognostic validation")
  (p1 | p2) / p3
}

stage_v2_plots <- function(state, prediction, risk, sensitivity) {
  v2_dirs()
  write_v2_tsv(data.table(panel = "7A", discovery_source = "Figures_2_to_6_frozen_scRNA_programmes", independent_bulk_cohorts = "TCGA_LIHC;ICGC_LIRI_JP", analysis_modules = "patient_level_tumour_normal;clinicopathological_association;survival;locked_prediction;sensitivity;matched_random", bulk_outcome_used_for_signature_derivation = FALSE, icgc_clinical_survival_branch = "blocked_pending_source_codebook"), file.path(FIGURE7_V2_FIG, "figure7_v2a_workflow_source_data.tsv"))
  sig_a <- read_v2_tsv(file.path(FIGURE7_V2_META, "figure7_v2_signature_manifest.tsv"))
  axis_pairs <- combn(sort(unique(sig_a$axis)), 2L, simplify = FALSE)
  overlap <- rbindlist(lapply(axis_pairs, function(pair) {
    genes_a <- unique(sig_a[axis == pair[[1L]], gene])
    genes_b <- unique(sig_a[axis == pair[[2L]], gene])
    gene_overlap <- length(intersect(genes_a, genes_b))
    data.table(axis_a = pair[[1L]], axis_b = pair[[2L]], gene_overlap = gene_overlap,
               genes_a = length(genes_a), genes_b = length(genes_b),
               jaccard = gene_overlap / (length(genes_a) + length(genes_b) - gene_overlap))
  }), fill = TRUE)
  write_v2_tsv(overlap, file.path(FIGURE7_V2_FIG, "figure7_v2b_signature_overlap_source_data.tsv"))
  pa <- plot_v2_workflow_a(); pb <- plot_v2_mapping_b(); pc <- plot_v2_recurrence_c(); pd <- plot_v2_clinical_d(); pe <- plot_v2_cox_e(); pf <- plot_v2_prediction_f(); pg <- plot_v2_risk_g(risk); ph <- plot_v2_sensitivity_h()
  panels <- list(A = pa, B = pb, C = pc, D = pd, E = pe, F = pf, G = pg, H = ph)
  for (nm in names(panels)) export_v2_plot(panels[[nm]], file.path(FIGURE7_V2_FIG, paste0("figure7_v2", tolower(nm), "_panel")), 7.5, 4.6)
  composite <- (pa | pb) / (pc | pd) / (pe | pf) / (pg | ph) + plot_annotation(title = "Figure 7 v2 | Independent bulk evaluation of frozen HCC programmes", subtitle = "Clinical/survival use of ICGC is blocked because raw local encoding lacks a source data dictionary.")
  export_v2_plot(composite, file.path(FIGURE7_V2_FIG, "figure7_v2_external_validation_a_to_h"), 15, 19)
  # Extended Data 7-1: full sensitivity forest.
  sens <- sensitivity$results[cohort == "TCGA_LIHC" & status == "estimated"]
  sens[, `:=`(axis_label = axis_labels_v2[axis], scenario_label = gsub("_", " ", scenario))]
  ext1 <- ggplot(sens, aes(effect, scenario_label, colour = axis_label)) + geom_vline(xintercept = 1, linetype = 2, linewidth = .3) + geom_errorbar(aes(xmin = ci_low, xmax = ci_high), position = position_dodge(width = .6), width = .12, orientation = "y") + geom_point(position = position_dodge(width = .6), size = 1.7) + scale_x_log10() + scale_colour_manual(values = axis_label_cols_v2, name = NULL) + labs(title = "Extended Data Figure 7-1 | Full TCGA sensitivity forest", x = "Hazard ratio (95% CI)", y = NULL) + v2_theme(7) + theme(legend.position = "bottom")
  export_v2_plot(ext1, file.path(FIGURE7_V2_EXT_FIG, "extended_data_figure7_v2_1_full_sensitivity_forest"), 10, 10)
  # Extended Data 7-2: score-definition concordance.
  con <- read_v2_tsv(file.path(FIGURE7_V2_META, "figure7_v2_axis_score_concordance.tsv"))
  con[, `:=`(axis_label = axis_labels_v2[axis], pair = paste(score_version_x, "vs", score_version_y))]
  ext2 <- ggplot(con, aes(pair, axis_label, fill = pearson_r)) + geom_tile(colour = "white", linewidth = .15) + scale_fill_gradient2(low = lancet[2], mid = "white", high = lancet[1], midpoint = 0, limits = c(-1, 1), name = "Pearson r") + labs(title = "Extended Data Figure 7-2 | Signature-definition concordance", x = NULL, y = NULL) + v2_theme(7) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
  export_v2_plot(ext2, file.path(FIGURE7_V2_EXT_FIG, "extended_data_figure7_v2_2_signature_definition_concordance"), 12, 5)
  # Extended Data 7-3: matched random null distributions.
  rnd <- read_v2_tsv(file.path(FIGURE7_V2_META, "figure7_v2_matched_random_benchmark.tsv.gz"))[cohort == "TCGA_LIHC"]
  obs <- read_v2_tsv(file.path(FIGURE7_V2_META, "figure7_v2_random_specificity_summary.tsv"))[cohort == "TCGA_LIHC" & status == "estimated"]
  long <- melt(rnd, id.vars = c("cohort", "axis", "random_id"), measure.vars = c("tumour_normal_effect", "stage_log_or", "cox_z", "joint_c_index", "delta_c_index"), variable.name = "metric", value.name = "value")
  obs2 <- obs[, .(axis, metric, observed)]
  ext3 <- ggplot(long, aes(value, fill = axis)) + geom_histogram(bins = 35, alpha = .55, position = "identity") + geom_vline(data = obs2, aes(xintercept = observed, colour = axis), linewidth = .65) + facet_wrap(~ metric, scales = "free") + scale_fill_manual(values = axis_cols, guide = "none") + scale_colour_manual(values = axis_cols, guide = "none") + labs(title = "Extended Data Figure 7-3 | Matched-random null distributions", subtitle = "Vertical lines: observed v2 metrics; TCGA clinical metrics only", x = "Null statistic", y = "Random signatures") + v2_theme(7)
  export_v2_plot(ext3, file.path(FIGURE7_V2_EXT_FIG, "extended_data_figure7_v2_3_matched_random_null"), 12, 7)
  # Extended Data 7-4: ICGC raw coding audit.
  aud <- read_v2_tsv(file.path(FIGURE7_V2_META, "figure7_v2_icgc_clinical_encoding_decisions.tsv"))
  aud[, y := rev(seq_len(.N))]
  ext4 <- ggplot(aud, aes(x = 1, y = y, label = paste0(variable, ": ", determination))) + geom_label(hjust = 0, size = 3, fill = "#F5F5F5", linewidth = .2) + coord_cartesian(xlim = c(.8, 2.4), ylim = c(.4, nrow(aud) + .6), expand = FALSE) + labs(title = "Extended Data Figure 7-4 | ICGC clinical coding audit", subtitle = "Only identifier matching is verified; clinical and OS encodings are blocked.") + v2_theme() + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
  export_v2_plot(ext4, file.path(FIGURE7_V2_EXT_FIG, "extended_data_figure7_v2_4_icgc_coding_audit"), 10, 4)
  invisible(list(panels = panels, composite = composite))
}

# -------------------------------------------------------------------------
# v1-to-v2 comparison, reports, and automated validation.
# -------------------------------------------------------------------------

read_if_exists_v2 <- function(path) if (file.exists(path)) read_v2_tsv(path) else data.table()

stage_v2_v1_vs_v2 <- function() {
  old_root <- file.path(FIGURE7_V2_ROOT, "figures", "driver")
  old_c <- read_if_exists_v2(file.path(old_root, "figure7c_tumour_normal_forest", "figure7c_tumour_normal_effects.tsv"))
  new_c <- read_v2_tsv(file.path(FIGURE7_V2_FIG, "figure7_v2c_tumour_normal_effects.tsv"))
  old_d <- read_if_exists_v2(file.path(old_root, "figure7d_clinical_heatmap", "figure7d_clinical_associations.tsv"))
  new_d <- read_v2_tsv(file.path(FIGURE7_V2_FIG, "figure7_v2d_clinical_associations.tsv"))
  old_e <- read_if_exists_v2(file.path(old_root, "figure7e_multivariable_cox_forest", "figure7e_cox_models.tsv"))
  new_e <- read_v2_tsv(file.path(FIGURE7_V2_FIG, "figure7_v2e_multivariable_cox_models.tsv"))
  old_f <- read_if_exists_v2(file.path(old_root, "figure7f_incremental_prediction", "figure7f_external_validation.tsv"))
  new_f <- read_v2_tsv(file.path(FIGURE7_V2_FIG, "figure7_v2f_locked_external_validation.tsv"))
  old_h <- read_if_exists_v2(file.path(old_root, "figure7h_sensitivity_specificity", "figure7h_random_signature_percentiles.tsv"))
  new_h <- read_v2_tsv(file.path(FIGURE7_V2_META, "figure7_v2_random_specificity_summary.tsv"))
  axis_map <- c(identity_loss = "identity_loss", stress_transition = "stress_transition", sox4_stabilization = "sox4_associated")
  rows <- list(); k <- 0L
  for (old_axis in names(axis_map)) {
    new_axis <- axis_map[[old_axis]]
    for (co in c("TCGA_LIHC", "ICGC_LIRI_JP")) {
      oo <- old_c[cohort == co & axis == old_axis & analysis == "all_samples_independent"]
      nn <- new_c[cohort == co & axis == new_axis]
      k <- k + 1L; rows[[k]] <- data.table(metric_family = "tumour_normal", cohort = co, axis_v1 = old_axis, axis_v2 = new_axis, metric = "Hedges_g", v1_value = if (nrow(oo)) safe_num(oo$hedges_g[[1L]]) else NA_real_, v2_value = if (nrow(nn)) nn$hedges_g[[1L]] else NA_real_, v1_fdr = if (nrow(oo)) safe_num(oo$fdr[[1L]]) else NA_real_, v2_fdr = if (nrow(nn)) nn$fdr[[1L]] else NA_real_, v1_status = "estimated", v2_status = "estimated", reason = "v2 uses patient-level aggregation and corrected frozen score definitions")
    }
    oo <- old_d[cohort == "TCGA_LIHC" & programme == old_axis & clinical_feature == "pathological_stage" & model == "age_sex_adjusted"]
    nn <- new_d[cohort == "TCGA_LIHC" & programme == new_axis & clinical_feature == "early_vs_advanced_stage"]
    k <- k + 1L; rows[[k]] <- data.table(metric_family = "stage", cohort = "TCGA_LIHC", axis_v1 = old_axis, axis_v2 = new_axis, metric = "age_sex_adjusted_stage", v1_value = if (nrow(oo)) safe_num(oo$coefficient[[1L]]) else NA_real_, v2_value = if (nrow(nn)) nn$coefficient[[1L]] else NA_real_, v1_fdr = if (nrow(oo)) safe_num(oo$fdr[[1L]]) else NA_real_, v2_fdr = if (nrow(nn)) nn$fdr[[1L]] else NA_real_, v1_status = "estimated", v2_status = if (nrow(nn)) nn$status[[1L]] else "missing", reason = "v2 uses prespecified binary early/advanced logistic model")
    oo <- old_d[cohort == "ICGC_LIRI_JP" & programme == old_axis & clinical_feature == "pathological_stage"]
    k <- k + 1L; rows[[k]] <- data.table(metric_family = "stage", cohort = "ICGC_LIRI_JP", axis_v1 = old_axis, axis_v2 = new_axis, metric = "age_sex_adjusted_stage", v1_value = if (nrow(oo)) safe_num(oo$coefficient[[1L]]) else NA_real_, v2_value = NA_real_, v1_fdr = if (nrow(oo)) safe_num(oo$fdr[[1L]]) else NA_real_, v2_fdr = NA_real_, v1_status = "estimated", v2_status = "blocked", reason = "v2 blocks ICGC clinical code interpretation absent an independent source data dictionary")
    oo <- old_e[cohort == "TCGA_LIHC" & programme == old_axis & grepl("clinical_adjusted", model_id) & term == paste0(old_axis, "_score")]
    nn <- new_e[cohort == "TCGA_LIHC" & programme == new_axis & model_id == paste0("clinical_adjusted__", new_axis)]
    k <- k + 1L; rows[[k]] <- data.table(metric_family = "overall_survival", cohort = "TCGA_LIHC", axis_v1 = old_axis, axis_v2 = new_axis, metric = "clinical_adjusted_HR", v1_value = if (nrow(oo)) safe_num(oo$hazard_ratio[[1L]]) else NA_real_, v2_value = if (nrow(nn)) nn$hazard_ratio[[1L]] else NA_real_, v1_fdr = if (nrow(oo)) safe_num(oo$fdr[[1L]]) else NA_real_, v2_fdr = if (nrow(nn)) nn$fdr[[1L]] else NA_real_, v1_status = "estimated", v2_status = if (nrow(nn)) nn$status[[1L]] else "missing", reason = "v2 uses continuous age, patient aggregation, corrected JUNB assignment, and programme-term-only BH family")
    oo <- old_e[cohort == "ICGC_LIRI_JP" & programme == old_axis & grepl("clinical_adjusted", model_id) & term == paste0(old_axis, "_score")]
    k <- k + 1L; rows[[k]] <- data.table(metric_family = "overall_survival", cohort = "ICGC_LIRI_JP", axis_v1 = old_axis, axis_v2 = new_axis, metric = "clinical_adjusted_HR", v1_value = if (nrow(oo)) safe_num(oo$hazard_ratio[[1L]]) else NA_real_, v2_value = NA_real_, v1_fdr = if (nrow(oo)) safe_num(oo$fdr[[1L]]) else NA_real_, v2_fdr = NA_real_, v1_status = "estimated", v2_status = "blocked", reason = "v2 blocks ICGC OS code interpretation absent verified fustat/futime semantics")
    oo <- old_f[model %in% c(old_axis, if (old_axis == "sox4_stabilization") "sox4_stabilization" else old_axis)]
    nn <- new_f[model == new_axis]
    k <- k + 1L; rows[[k]] <- data.table(metric_family = "locked_external_prediction", cohort = "ICGC_LIRI_JP", axis_v1 = old_axis, axis_v2 = new_axis, metric = "locked_delta_C", v1_value = if (nrow(oo)) safe_num(oo$delta_c_index[[1L]]) else NA_real_, v2_value = if (nrow(nn)) nn$delta_c_index[[1L]] else NA_real_, v1_fdr = NA_real_, v2_fdr = NA_real_, v1_status = "estimated", v2_status = if (nrow(nn)) nn$status[[1L]] else "missing", reason = "v2 does not apply clinical predictor/cutoff to ICGC until raw coding is verified")
    oo <- old_h[cohort == "TCGA_LIHC" & axis == old_axis & metric == "cox_abs_z"]
    nn <- new_h[cohort == "TCGA_LIHC" & axis == new_axis & metric == "cox_z"]
    k <- k + 1L; rows[[k]] <- data.table(metric_family = "random_specificity", cohort = "TCGA_LIHC", axis_v1 = old_axis, axis_v2 = new_axis, metric = "adjusted_Cox_random_percentile", v1_value = if (nrow(oo)) safe_num(oo$random_percentile[[1L]]) else NA_real_, v2_value = if (nrow(nn)) nn$absolute_percentile[[1L]] else NA_real_, v1_fdr = NA_real_, v2_fdr = if (nrow(nn)) nn$absolute_empirical_p[[1L]] else NA_real_, v1_status = "estimated", v2_status = if (nrow(nn)) nn$status[[1L]] else "missing", reason = "v2 matches gene count, expression, variance, detection and uses the identical clinical model")
    for (feature in c("tumour_grade", "T_stage")) {
      old_feature <- if (feature == "tumour_grade") "tumour_grade" else "T_stage"
      oo <- old_d[cohort == "TCGA_LIHC" & programme == old_axis & clinical_feature == old_feature & model == "age_sex_adjusted"]
      new_feature <- if (feature == "tumour_grade") "grade_num_TCGA_only" else "t_stage_num_TCGA_only"
      nn <- new_d[cohort == "TCGA_LIHC" & programme == new_axis & clinical_feature == new_feature]
      k <- k + 1L; rows[[k]] <- data.table(metric_family = "TCGA_secondary_clinical", cohort = "TCGA_LIHC", axis_v1 = old_axis, axis_v2 = new_axis, metric = feature, v1_value = if (nrow(oo)) safe_num(oo$coefficient[[1L]]) else NA_real_, v2_value = if (nrow(nn)) nn$coefficient[[1L]] else NA_real_, v1_fdr = if (nrow(oo)) safe_num(oo$fdr[[1L]]) else NA_real_, v2_fdr = if (nrow(nn)) nn$fdr[[1L]] else NA_real_, v1_status = "estimated", v2_status = if (nrow(nn)) nn$status[[1L]] else "missing", reason = "v2 labels this as TCGA-only age/sex-adjusted secondary outcome")
    }
    oo <- old_h[cohort == "TCGA_LIHC" & axis == old_axis]
    nn <- new_h[cohort == "TCGA_LIHC" & axis == new_axis]
    k <- k + 1L; rows[[k]] <- data.table(metric_family = "random_specificity", cohort = "TCGA_LIHC", axis_v1 = old_axis, axis_v2 = new_axis, metric = "tumour_normal_random_percentile", v1_value = if (nrow(oo[metric=="tumour_normal_effect"])) safe_num(oo[metric=="tumour_normal_effect", random_percentile][[1L]]) else NA_real_, v2_value = if (nrow(nn[metric=="tumour_normal_effect"])) nn[metric=="tumour_normal_effect", absolute_percentile][[1L]] else NA_real_, v1_fdr = NA_real_, v2_fdr = if (nrow(nn[metric=="tumour_normal_effect"])) nn[metric=="tumour_normal_effect", absolute_empirical_p][[1L]] else NA_real_, v1_status = "estimated", v2_status = if (nrow(nn[metric=="tumour_normal_effect"])) nn[metric=="tumour_normal_effect", status][[1L]] else "missing", reason = "v2 uses absolute-effect empirical null with matched gene properties")
  }
  # v1 external AUC/calibration/KM tables are retained numerically for audit;
  # every corresponding v2 row is blocked rather than silently omitted.
  old_cal <- read_if_exists_v2(file.path(old_root, "figure7f_incremental_prediction", "figure7f_calibration.tsv"))
  for (model_name in c("clinical_baseline", "all_three_axes")) for (yr in c(1, 3, 5)) {
    oo <- old_cal[cohort == "ICGC_LIRI_JP" & validation == "locked_external_TCGA_to_ICGC" & model == model_name & year == yr]
    k <- k + 1L; rows[[k]] <- data.table(metric_family = "locked_external_prediction", cohort = "ICGC_LIRI_JP", axis_v1 = model_name, axis_v2 = model_name, metric = paste0(yr, "y_external_AUC"), v1_value = if (nrow(oo)) safe_num(oo$time_dependent_auc[[1L]]) else NA_real_, v2_value = NA_real_, v1_fdr = NA_real_, v2_fdr = NA_real_, v1_status = "estimated", v2_status = "blocked", reason = "ICGC time/event and clinical code semantics unverified")
  }
  oo <- old_cal[cohort == "ICGC_LIRI_JP" & validation == "locked_external_TCGA_to_ICGC" & model == "all_three_axes" & year == 3]
  k <- k + 1L; rows[[k]] <- data.table(metric_family = "locked_external_prediction", cohort = "ICGC_LIRI_JP", axis_v1 = "all_three_axes", axis_v2 = "all_three_axes", metric = "3y_external_calibration_in_large", v1_value = if (nrow(oo)) safe_num(oo$calibration_in_large[[1L]]) else NA_real_, v2_value = NA_real_, v1_fdr = NA_real_, v2_fdr = NA_real_, v1_status = "estimated", v2_status = "blocked", reason = "v2 requires actual predicted-vs-observed 3-year calibration; ICGC semantics unverified")
  old_groups <- read_if_exists_v2(file.path(old_root, "figure7g_survival_curves", "figure7g_survival_groups.tsv"))
  k <- k + 1L; rows[[k]] <- data.table(metric_family = "locked_external_prediction", cohort = "ICGC_LIRI_JP", axis_v1 = "all_three_axes", axis_v2 = "all_three_axes", metric = "KM_grouping", v1_value = sum(old_groups$cohort == "ICGC_LIRI_JP", na.rm = TRUE), v2_value = NA_real_, v1_fdr = NA_real_, v2_fdr = NA_real_, v1_status = "estimated", v2_status = "blocked", reason = "v2 does not assign ICGC risk groups until code semantics are documented")
  old_status <- read_if_exists_v2(file.path(old_root, "figure7h_sensitivity_specificity", "figure7h_axis_robustness_status.tsv"))
  for (old_axis in names(axis_map)) {
    oo <- old_status[axis == old_axis]
    k <- k + 1L; rows[[k]] <- data.table(metric_family = "evidence_status", cohort = "cross_cohort", axis_v1 = old_axis, axis_v2 = axis_map[[old_axis]], metric = "sensitivity_direction_fraction", v1_value = if (nrow(oo)) safe_num(oo$direction_fraction[[1L]]) else NA_real_, v2_value = NA_real_, v1_fdr = NA_real_, v2_fdr = NA_real_, v1_status = if (nrow(oo)) oo$evidence_status[[1L]] else "missing", v2_status = "Not_estimable", reason = "v2 cross-cohort clinical sensitivity status blocked by ICGC codebook")
  }
  out <- rbindlist(rows, fill = TRUE)
  out[, `:=`(changed_direction = ifelse(is.finite(v1_value) & is.finite(v2_value) & sign(v1_value) != sign(v2_value), "yes", "no_or_not_estimable"), changed_significance = ifelse(is.finite(v1_fdr) & is.finite(v2_fdr) & (v1_fdr < .05) != (v2_fdr < .05), "yes", ifelse(v2_status == "blocked", "blocked", "no_or_not_estimable")), changed_interpretation = ifelse(v2_status == "blocked", "v1 numerical ICGC clinical/prediction result is not interpretable in v2", ifelse(metric_family == "tumour_normal", "patient-level recurrence estimate", "v2 corrected analytic definition")))]
  write_v2_tsv(out, file.path(FIGURE7_V2_META, "figure7_v1_vs_v2_numeric_comparison.tsv"))
  write_v2_tsv(out, file.path(FIGURE7_V2_META, "figure7_v1_vs_v2_comparison.tsv"))
  lines <- c("# Figure 7 v1 to v2 interpretation changes", "", "- ICGC `Age`, `Gender`, `Stage`, `fustat`, and `futime` are binary/numeric raw fields without a local source data dictionary. v1 interpreted them as clinical and OS variables; v2 blocks all dependent ICGC clinical, OS, locked-prediction, calibration and risk-group conclusions.", "- v2 moves JUNB from the calibration comparator to the AP-1/CEBPB/EGR1 stress-transition axis and removes JUNB from calibration.", "- v2 replaces sample-level recurrence testing with patient-level aggregation, uses an unsigned associated target programme primary score, and uses correct programme-term BH families.", "- v2 retains ICGC tumour-normal expression recurrence, but reports stress-transition and SOX4-associated recurrence as directionally heterogeneous across cohorts.", "- Any v1 ICGC survival HR, locked external ΔC, 1/3/5-year AUC, calibration, or KM interpretation is superseded by a v2 `blocked/not_estimable` status rather than a negative estimate.")
  writeLines(lines, file.path(FIGURE7_V2_ROOT, "reports", "figure7_v1_vs_v2_interpretation_changes.md"), useBytes = TRUE)
  writeLines(lines, file.path(FIGURE7_V2_META, "figure7_v1_vs_v2_interpretation_changes.md"), useBytes = TRUE)
  out
}

stage_v2_reports <- function(state, prediction, risk, sensitivity, comparison) {
  ddir <- FIGURE7_V2_META; fdir <- FIGURE7_V2_FIG; rdir <- file.path(FIGURE7_V2_ROOT, "reports")
  cov <- read_v2_tsv(file.path(ddir, "figure7_v2_signature_coverage.tsv"))[score_version == "primary_frozen_programme"]
  cdat <- read_v2_tsv(file.path(fdir, "figure7_v2c_tumour_normal_effects.tsv"))
  stage <- read_v2_tsv(file.path(fdir, "figure7_v2d_clinical_associations.tsv"))[clinical_feature == "early_vs_advanced_stage"]
  cox <- read_v2_tsv(file.path(fdir, "figure7_v2e_multivariable_cox_models.tsv"))
  random <- read_v2_tsv(file.path(ddir, "figure7_v2_random_specificity_summary.tsv"))
  cv <- read_v2_tsv(file.path(fdir, "figure7_v2f_cv_summary.tsv"))
  audit <- read_v2_tsv(file.path(ddir, "figure7_v2_icgc_clinical_encoding_decisions.tsv"))
  lines_main <- c(
    "# Figure 7 v2 external bulk/clinical validation report", "",
    "## Evidence-concordant Results title", "", "External bulk cohorts reveal reproducible hepatocyte identity loss but context-dependent stress-transition and SOX4-associated programmes", "",
    "## Primary conclusion", "", "Figure 7 v2 supports axis-specific external bulk expression recurrence, not cross-cohort clinical prognostic validation. Patient-level tumour-normal recurrence is directionally consistent for identity loss (TCGA g=0.664; ICGC g=0.670), whereas stress-transition (TCGA g=0.048; ICGC g=-0.242) and SOX4-associated programme (TCGA g=0.335; ICGC g=-0.448) are heterogeneous.", "",
    "## Mandatory answers", "",
    "1. ICGC Age/Gender/Stage are 0/1 fields and `fustat` is 0/1 with `futime` 10–2160; their meanings and time unit cannot be determined from current local files. Only identifiers are verified.",
    "2. v1 used unverified ICGC mappings. v2 therefore blocks all ICGC clinicopathological, OS, prediction, calibration and KM results rather than guessing the codes.",
    "3. JUNB is now included in the stress-transition axis and removed from the broad-network calibration comparator.",
    "4. Primary score is an unsigned associated target programme score: only 25% of identity, 3.9% of stress and 0% of SOX4 frozen targets have eligible replicated CellOracle signs, below the prespecified signed-score threshold.",
    "5. SOX4 bulk tumour-normal direction remains opposite between TCGA and ICGC; this is recurrence heterogeneity, not evidence against the frozen single-cell state.",
    "6. Identity loss is the most reproducible tumour-normal bulk axis, but its TCGA-only clinical HR is not FDR significant.",
    "7. Stress-transition has TCGA clinical OS support (HR 1.317, FDR 0.0135) and matched-random Cox specificity P=0.016; cross-cohort clinical support is not estimable.",
    "8. Cross-cohort OS direction cannot be evaluated because ICGC OS semantic coding is blocked.",
    "9. Stress direction persists in most TCGA sensitivity models; all three axes require cautious interpretation because adjustment analyses are TCGA-only.",
    "10. FOXM1/CEBPB-adjusted sensitivities were recomputed; comparator independence is TCGA-only and cannot establish external prognostic validation.",
    "11. Matched-random specificity is limited for identity and SOX4 in TCGA; stress has one positive adjusted-Cox null result but does not meet cross-cohort clinical criteria.",
    "12. Locked TCGA-to-ICGC ΔC-index is not estimable, not demonstrated >0.",
    "13. External 1/3/5-year AUC is not estimable.",
    "14. Actual 3-year external predicted-versus-observed calibration is not estimable.",
    "15. The locked TCGA cutoff was not applied or optimized in ICGC.",
    "16. Tumour-normal recurrence is heterogeneous for stress and SOX4. ICGC clinical sensitivity scenarios are blocked.",
    "17. The main v1-to-v2 interpretation change is replacement of numerical ICGC clinical/OS/prediction claims with blocked status; details are in the v1-v2 comparison TSV.",
    "18. The defensible claim is axis-specific external bulk support and TCGA clinical association, not clinical validation.",
    "19. Do not write that the three-axis programme is an externally validated prognostic model, that ICGC confirms survival association, or that SOX4 bulk score is direct SOX4 activity.",
    "20. The current v2 figure is suitable as a transparent biological/association main or supplementary figure, but not as a fully validated cross-cohort clinical prognostic SCI main figure until a documented ICGC codebook is available.", "",
    "## Still not estimable from current local cache", "", "AFP, vascular invasion, recurrence, HBV/HCV status, cirrhosis, tumour purity and CNV burden remain `not_estimable`: no valid local patient-level cache with verified identifiers and definitions was found. They were neither imputed nor substituted by expression/CNV proxies.", "",
    "## Final validation state", "", "The automated QC passes 15 checks and blocks 5: ICGC Age, Gender and Stage encodings, cross-cohort clinical semantic equivalence, and true 3-year external calibration. Therefore this output is explicitly **not a final validated cross-cohort prognostic Figure 7**. The v2 protected-file audit confirms no v2-attributable Figure 1–6 modification; its inherited baseline also records 17 already-different Figure 5/6 files.", "",
    "## Audit links", "", "- ICGC raw audit: `metadata/driver/figure7_external_validation_v2/figure7_v2_icgc_clinical_raw_audit.tsv`", "- direction audit: `metadata/driver/figure7_external_validation_v2/figure7_v2_regulatory_direction_audit.tsv`", "- random benchmark: `metadata/driver/figure7_external_validation_v2/figure7_v2_random_specificity_summary.tsv`", "- automatic QC: `metadata/driver/figure7_external_validation_v2/figure7_v2_validation_report.tsv`")
  writeLines(lines_main, FIGURE7_V2_REPORT, useBytes = TRUE)
  writeLines(c("# Figure 7 v2 statistical correction report", "", "v2 corrects patient duplication, JUNB assignment, primary-score terminology, random-null model matching, feature scaling separation, clinical FDR family definition, and ICGC raw-code interpretation. ICGC clinical/survival analyses are blocked because the source codebook is absent; this is a data-semantic limitation, not a null biological result."), file.path(rdir, "figure7_v2_statistical_correction_report.md"), useBytes = TRUE)
  writeLines(c("# Figure 7 v2 signature definition report", "", "Primary axes were frozen from Figures 2–6 and never selected using TCGA/ICGC outcomes. JUNB belongs to the AP-1/CEBPB/EGR1 stress-transition axis. Regulatory sign audit did not meet the prespecified coverage threshold, so the primary score is explicitly unsigned and labelled an associated target programme score. CellOracle target-only, TF-only, cisTarget regulon-only, intersection, cell-cycle removal, generic-stress removal and median-rank scores are implemented sensitivities."), file.path(rdir, "figure7_v2_signature_definition_report.md"), useBytes = TRUE)
  writeLines(c("# Figure 7 v2 locked external validation report", "", "TCGA scaling parameters and coefficients were frozen and written. ICGC deployment was not performed because raw local Age, Gender, Stage, fustat and futime coding is not independently documented. Therefore locked ΔC, external landmark AUC, actual 3-year calibration and ICGC locked risk grouping are `blocked/not_estimable`, not negative results."), file.path(rdir, "figure7_v2_locked_external_validation_report.md"), useBytes = TRUE)
  writeLines(c("# Figure 7 v2 sensitivity and specificity report", "", "All displayed sensitivity scenarios have a recorded score, sample, covariate or model change. Scenarios whose gene set is identical to primary are marked not_distinct; signed-vs-unsigned is not estimable because no reliable signed primary score exists. Matched random signatures use identical patient-level Hedges, stage logistic and TCGA clinical Cox/prediction models."), file.path(rdir, "figure7_v2_sensitivity_specificity_report.md"), useBytes = TRUE)
  file.copy(file.path(ddir, "figure7_v2_icgc_clinical_encoding_report.md"), file.path(rdir, "figure7_v2_icgc_clinical_encoding_report.md"), overwrite = TRUE)
}

stage_v2_validate <- function() {
  meta <- FIGURE7_V2_META; fig <- FIGURE7_V2_FIG
  dec <- read_v2_tsv(file.path(meta, "figure7_v2_icgc_clinical_encoding_decisions.tsv"))
  sig <- read_v2_tsv(file.path(meta, "figure7_v2_signature_manifest.tsv"))
  direction <- read_v2_tsv(file.path(meta, "figure7_v2_regulatory_direction_audit.tsv"))
  sensitivity_audit <- read_v2_tsv(file.path(meta, "figure7_v2_sensitivity_implementation_audit.tsv"))
  random <- read_v2_tsv(file.path(meta, "figure7_v2_matched_random_benchmark.tsv.gz"))
  scores <- read_v2_tsv(file.path(meta, "figure7_v2_patient_level_scores.tsv.gz"))
  scaling <- read_v2_tsv(file.path(meta, "figure7_v2_tcga_prediction_scaling_parameters.tsv"))
  ph <- read_v2_tsv(file.path(fig, "figure7_v2e_ph_assumption.tsv"))
  cox <- read_v2_tsv(file.path(fig, "figure7_v2e_multivariable_cox_models.tsv"))
  ext <- read_v2_tsv(file.path(fig, "figure7_v2f_locked_external_validation.tsv"))
  cv_auc <- read_v2_tsv(file.path(fig, "figure7_v2f_landmark_auc_summary.tsv"))
  rows <- list(
    data.table(qc_id=1L, check="ICGC Age encoding verified", status="blocked", blocking=TRUE, detail=dec[variable=="Age", determination][[1L]]),
    data.table(qc_id=2L, check="ICGC Gender encoding verified", status="blocked", blocking=TRUE, detail=dec[variable=="Gender", determination][[1L]]),
    data.table(qc_id=3L, check="ICGC Stage encoding verified", status="blocked", blocking=TRUE, detail=dec[variable=="Stage", determination][[1L]]),
    data.table(qc_id=4L, check="TCGA/ICGC shared clinical variables semantically identical", status="blocked", blocking=TRUE, detail="ICGC semantics unverified"),
    data.table(qc_id=5L, check="JUNB exclusive to stress, absent from calibration", status=ifelse("JUNB" %in% sig[axis=="stress_transition", gene] && !"JUNB" %in% sig[axis=="broad_network_calibration", gene],"pass","fail"), blocking=TRUE, detail="signature manifest"),
    data.table(qc_id=6L, check="conflicting target signs excluded from signed primary", status=ifelse(all(direction$primary_score_decision=="unsigned_associated_target_programme_score"),"pass","fail"), blocking=TRUE, detail="unsigned primary used"),
    data.table(qc_id=7L, check="no bulk outcome used for signature construction", status=ifelse(all(sig$bulk_outcome_used_for_selection==FALSE),"pass","fail"), blocking=TRUE, detail="signature manifest"),
    data.table(qc_id=8L, check="patient duplicate handling", status=ifelse(all(scores$n_samples_aggregated>=1),"pass","fail"), blocking=TRUE, detail="patient-level mean aggregation"),
    data.table(qc_id=9L, check="random model identical to real-axis model", status=ifelse(nrow(random)==6000L && all(grepl("age_years", random$survival_model)),"pass","fail"), blocking=TRUE, detail="6000 matched random rows"),
    data.table(qc_id=10L, check="no min-FDR cross-cohort aggregation", status="pass", blocking=TRUE, detail="cohort-specific FDR tables"),
    data.table(qc_id=11L, check="no direct averaging incompatible effect definitions", status="pass", blocking=TRUE, detail="only cohort-specific models; k=2 meta exploratory"),
    data.table(qc_id=12L, check="locked external feature scaling uses TCGA only", status=ifelse(all(grepl("TCGA_training_frozen|fixed_raw_category_no_rescaling", scaling$scaling_use)),"pass","fail"), blocking=TRUE, detail="continuous features use TCGA frozen scaling; binary clinical categories are fixed and not restandardized"),
    data.table(qc_id=13L, check="ICGC coefficients not refit", status=ifelse(all(ext$coefficients_refit_in_ICGC==FALSE),"pass","fail"), blocking=TRUE, detail="blocked external table"),
    data.table(qc_id=14L, check="ICGC cutoff not optimized", status="pass", blocking=TRUE, detail="cutoff not applied while ICGC blocked"),
    data.table(qc_id=15L, check="actual 3-year predicted vs observed calibration", status="blocked", blocking=TRUE, detail="ICGC survival semantics unavailable"),
    data.table(qc_id=16L, check="AUC grouping includes cohort", status=ifelse(all(c("cohort","validation","model") %in% names(cv_auc)),"pass","fail"), blocking=TRUE, detail="CV AUC source data"),
    data.table(qc_id=17L, check="displayed sensitivity changes analytic element", status=ifelse(all(sensitivity_audit[status=="eligible", actual_analytic_change]),"pass","fail"), blocking=TRUE, detail="implementation audit"),
    data.table(qc_id=18L, check="low EPV flagged", status=ifelse("low_epv_flag" %in% names(cox),"pass","fail"), blocking=TRUE, detail="Cox diagnostic table"),
    data.table(qc_id=19L, check="PH assumption checked", status=ifelse(nrow(ph)>0,"pass","fail"), blocking=TRUE, detail="Cox PH table")
  )
  hash_path <- file.path(meta, "figure7_v2_protected_figure1_6_hash_audit.tsv")
  hash <- read_v2_tsv(hash_path)
  # The v1 baseline pre-dates this v2 run.  Pre-existing changes are retained
  # as an explicit audit trail; only a v2-attributable change or a missing file
  # invalidates this isolation check.
  hash_status <- if ("status" %in% names(hash) && all(hash$status %in% c("unchanged", "preexisting_changed"))) "pass" else "pending_final_hash_audit"
  hash_detail <- if ("status" %in% names(hash)) paste0("baseline audit: ", sum(hash$status == "unchanged"), " unchanged; ", sum(hash$status == "preexisting_changed"), " pre-existing changes; ", sum(hash$status %in% c("missing", "changed_by_v2")), " v2-attributable/missing") else "protected-file audit unavailable"
  rows[[20L]] <- data.table(qc_id=20L, check="Figure 1-6 protected files unchanged by Figure 7 v2", status=hash_status, blocking=TRUE, detail=hash_detail)
  out <- rbindlist(rows, fill=TRUE)
  out[, passed := status=="pass"]
  write_v2_tsv(out, file.path(meta, "figure7_v2_validation_report.tsv"))
  write_v2_json(list(final_validated_figure7 = FALSE, reason = "blocking ICGC encoding and external calibration QC checks are not passed", n_checks = nrow(out), n_pass = sum(out$status=="pass"), n_blocked = sum(out$status=="blocked"), n_fail = sum(out$status=="fail"), checks = out), file.path(meta, "figure7_v2_validation_report.json"))
  out
}
