#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
suppressPackageStartupMessages(library(data.table))

paths <- figure5_six_panel_paths(root, create = FALSE)
checks <- data.table(id = integer(), check = character(), passed = logical(), detail = character())
add_check <- function(id, check, passed, detail = "") {
  checks <<- rbind(checks, data.table(id = as.integer(id), check = check, passed = isTRUE(passed), detail = as.character(detail)))
}

report_path <- function(name) file.path(paths$metadata, name)
read_report <- function(name) jsonlite::read_json(report_path(name), simplifyVector = TRUE)
file_exists <- function(x) all(file.exists(x))
all_formats <- function(stem) all(file.exists(paste0(stem, c(".pdf", ".png", ".svg", ".tiff"))))

plot_objects <- file.path(paths$processed, "plot_objects", paste0("figure5_panel_", LETTERS[1:6], ".rds"))
add_check(1, "six main panel objects A-F exist", file_exists(plot_objects), paste(plot_objects[!file.exists(plot_objects)], collapse = ";"))
add_check(2, "six-panel main namespace exists", dir.exists(paths$figures) && dir.exists(paths$preview), paths$figures)

a_text <- paste(readLines(file.path(root, "scripts", "plot_figure5a_corrected_workflow.R"), warn = FALSE), collapse = "\n")
a_report <- read_report("figure5_six_panel_A_report.json")
a_nodes_text <- gsub("[[:space:]]+", " ", paste(a_report$workflow_nodes, collapse = " "))
a_formal <- as.character(a_report$workflow_nodes_formal)
a_expected_formal <- c(
  "Frozen hepatocyte atlas and metadata",
  "Frozen three-axis programmes: HNF4A/PPARA identity loss; AP-1/CEBPB/EGR1 stress transition; SOX4 malignant-state stabilization",
  "Oriented consensus pseudotime: 0 = normal/reference; 1 = CNV-supported malignant/malignant-like",
  "CNV-supported malignant endpoint and CellRank malignant-fate probability",
  "Coverage-qualified patient/sample-token pseudobulk GAM",
  "Corrected temporal landmarks: onset; t10; t50; maximum slope",
  "Dataset-stratified bootstrap",
  "Tie-aware precedence probability",
  "Conservative overlapping regulatory-activity model"
)
a_required <- c("Frozen hepatocyte atlas", "Frozen three-axis programmes", "Oriented consensus pseudotime",
                "CNV-supported malignant endpoint", "CellRank malignant-fate probability",
                "Coverage-qualified patient/sample-token", "pseudobulk GAM", "Corrected temporal landmarks",
                "Dataset-stratified bootstrap", "Tie-aware precedence probability",
                "Conservative overlapping", "regulatory-activity model")
add_check(3, "5A contains all nine formal workflow stages", identical(a_formal, a_expected_formal), paste(setdiff(a_expected_formal, a_formal), collapse = ";"))
a_forbidden <- c("DPT (not available)", "CellRank pseudotime", "Patient pseudobulk", "Leave-one-dataset-out")
add_check(4, "5A omits forbidden independent-method/stage wording", !any(vapply(a_forbidden, grepl, logical(1), x = a_text, fixed = TRUE)),
          paste(a_forbidden[vapply(a_forbidden, grepl, logical(1), x = a_text, fixed = TRUE)], collapse = ";"))
add_check(5, "5A descriptive state-anchor note is exact", identical(a_report$state_anchor_note, "State labels provide descriptive anchors and do not define discrete temporal transitions."), a_report$state_anchor_note)
add_check(6, "5A bootstrap fallback note is exact", identical(a_report$bootstrap_coverage_note, "Because eligible patient-token data did not span the complete pseudotemporal range, the primary bootstrap used the prespecified sample-token coverage fallback."), a_report$bootstrap_coverage_note)

panel_reports <- c(A = "figure5_six_panel_A_report.json", B = "figure5_six_panel_B_report.json", C = "figure5_six_panel_C_report.json",
                   D = "figure5_six_panel_D_report.json", E = "figure5_six_panel_E_report.json", F = "figure5_six_panel_F_report.json")
add_check(7, "panel A-F reports exist", file_exists(report_path(panel_reports)), paste(panel_reports[!file.exists(report_path(panel_reports))], collapse = ";"))
expected_titles <- c(
  A = "Corrected temporal-positioning analysis workflow",
  B = "Coverage-corrected regulatory programmes along unified pseudotime",
  C = "Temporal organization of frozen TF, regulon and target-gene programmes",
  D = "Corrected bootstrap temporal landmarks",
  E = "Tie-aware temporal precedence across complementary landmarks",
  F = "Conservative overlapping regulatory-activity model"
)
observed_titles <- vapply(panel_reports, function(name) {
  if (!file.exists(report_path(name))) return(NA_character_)
  as.character(read_report(name)$title)
}, character(1))
add_check(8, "panel A-F titles match final specification", all(observed_titles == expected_titles, na.rm = TRUE), paste(names(observed_titles)[observed_titles != expected_titles], collapse = ";"))

d_plot <- fread(file.path(paths$metadata, "figure5_six_panel_D_temporal_landmarks.tsv"))
add_check(9, "5D retains peak and core landmarks", all(c("onset_time", "t10", "t50", "maximum_slope_time", "peak_time") %chin% d_plot$landmark), paste(unique(d_plot$landmark), collapse = ";"))

e_data <- fread(file.path(paths$metadata, "figure5_six_panel_E_precedence_probabilities.tsv"))
e_levels <- c("onset", "t10", "t50", "maximum_slope")
add_check(10, "5E contains corrected onset, t10, t50 and maximum slope", all(e_levels %chin% e_data$landmark), paste(setdiff(e_levels, unique(e_data$landmark)), collapse = ";"))
add_check(11, "5E excludes peak", !any(e_data$landmark %chin% c("peak", "peak_time")), paste(unique(e_data$landmark), collapse = ";"))
add_check(12, "5E labels precedence as Pr(earlier)", all(grepl("Pr\\(earlier\\)", e_data$label, fixed = FALSE) | e_data$label == "Not available"), "label audit")
e_text <- paste(readLines(file.path(root, "scripts", "plot_figure5e_six_panel_precedence.R"), warn = FALSE), collapse = "\n")
add_check(13, "5E does not use P= notation", !grepl("P=", e_text, fixed = TRUE), "static script audit")

f_data <- fread(file.path(paths$metadata, "figure5_six_panel_F_activity_band_summary.tsv"))
f_report <- read_report("figure5_six_panel_F_report.json")
f_text <- paste(readLines(file.path(root, "scripts", "plot_figure5f_six_panel_activity_model.R"), warn = FALSE), collapse = "\n")
add_check(14, "5F excludes onset_time from formal boundary", identical(f_report$uses_onset_time_for_formal_band, FALSE), as.character(f_report$uses_onset_time_for_formal_band))
add_check(15, "5F stable starts use bootstrap t10", all(f_data$boundary_method[f_data$t10_stable] == "bootstrap_t10"), paste(unique(f_data$boundary_method), collapse = ";"))
add_check(16, "5F preserves corrected onset in source data", all(is.finite(f_data$onset)), paste(f_data$onset, collapse = ";"))
add_check(17, "5F has exact activity-band caption", identical(f_report$caption, "Activity bands represent relative programme prominence derived from smoothed scores and bootstrap temporal landmarks; they do not indicate discrete activation or termination events."), f_report$caption)
add_check(18, "5F uses polygons/alpha and no hard rectangles", grepl("geom_polygon", f_text, fixed = TRUE) && grepl("scale_alpha_identity", f_text, fixed = TRUE) && !grepl("geom_rect", f_text, fixed = TRUE), "continuous alpha-gradient band audit")
add_check(19, "5F uses visual right-edge fade without programme end", all(is.finite(f_data$fade_start) & is.finite(f_data$fade_end)) && all(grepl("fade", f_data$right_edge_rule, fixed = TRUE)), paste(unique(f_data$right_edge_rule), collapse = ";"))
add_check(20, "5F reports unresolved boundaries when needed", is.character(f_data$boundary_status) && all(f_data$boundary_status %chin% c("resolved", "boundary unresolved")), paste(unique(f_data$boundary_status), collapse = ";"))

x_data <- fread(file.path(paths$extended_metadata, "extended_data_figureX_temporal_method_concordance.tsv"))
old_x <- fread(file.path(figure5_onset_fix_paths(root, create = FALSE)$metadata, "figure5f_method_concordance.tsv"))
old_x <- old_x[!is.na(method) & !is.na(comparison)]
common_x <- intersect(names(old_x), names(x_data))
setorderv(old_x, common_x); setorderv(x_data, common_x)
same_x <- nrow(old_x) == nrow(x_data) && isTRUE(all.equal(old_x[, ..common_x], x_data[, ..common_x], check.attributes = FALSE))
add_check(21, "Extended Data Figure X retains complete old method rows", same_x, sprintf("old=%d new=%d", nrow(old_x), nrow(x_data)))
add_check(22, "Extended Data Figure X retains Not available rows", any(x_data$status == "Not available") && any(x_data$n_valid == 0), paste(unique(x_data$status), collapse = ";"))
add_check(23, "Extended Data Figure X title is exact", identical(jsonlite::read_json(file.path(paths$extended_metadata, "extended_data_figureX_report.json"), simplifyVector = TRUE)$title, "Availability and sensitivity of temporal-order estimates across trajectory methods and resampling schemes"), "title audit")

y_data <- fread(file.path(paths$extended_metadata, "extended_data_figureY_patient_token_coverage.tsv"))
y_report <- jsonlite::read_json(file.path(paths$extended_metadata, "extended_data_figureY_report.json"), simplifyVector = TRUE)
y_required <- c("patient_sample_token", "pseudotime_minimum", "pseudotime_maximum", "number_of_bins", "coverage_status", "exclusion_reason")
add_check(24, "Extended Data Figure Y has required coverage columns", all(y_required %chin% names(y_data)), paste(setdiff(y_required, names(y_data)), collapse = ";"))
add_check(25, "Extended Data Figure Y retains exact coverage limitation statement", identical(y_report$body, "No individual patient token covered the full pseudotemporal range required for independent within-patient ordering."), y_report$body)
add_check(26, "Extended Data Figure Y has no qualified patient token", identical(as.integer(y_report$n_patient_tokens_coverage_qualified), 0L), as.character(y_report$n_patient_tokens_coverage_qualified))

montage_text <- paste(readLines(file.path(root, "scripts", "plot_figure5_six_panel_montage.R"), warn = FALSE), collapse = "\n")
add_check(27, "montage uses requested four-row layout", all(vapply(c("plots$A /", "(plots$B | plots$C)", "(plots$D | plots$E)", "plots$F +"), grepl, logical(1), x = montage_text, fixed = TRUE)), "row layout audit")
layout_map <- fread(file.path(paths$metadata, "figure5_six_panel_layout_old_vs_new.tsv"))
add_check(28, "old-versus-new mapping contains all old panels", nrow(layout_map) == 8L && all(paste0("5", LETTERS[1:8]) %chin% layout_map$old_panel), paste(layout_map$old_panel, collapse = ";"))

main_stems <- c(
  file.path(paths$figures, "panel_A_corrected_workflow", "figure5_six_panel_A_corrected_workflow"),
  file.path(paths$figures, "panel_B_programmes", "figure5_six_panel_B_programmes"),
  file.path(paths$figures, "panel_C_programme_heatmap", "figure5_six_panel_C_programme_heatmap"),
  file.path(paths$figures, "panel_D_landmarks", "figure5_six_panel_D_landmarks"),
  file.path(paths$figures, "panel_E_precedence", "figure5_six_panel_E_precedence"),
  file.path(paths$figures, "panel_F_activity_model", "figure5_six_panel_F_activity_model"),
  file.path(paths$preview, "figure5_six_panel_main")
)
extended_stems <- c(file.path(paths$extended_figures, "figureX_temporal_method_concordance", "extended_data_figureX_temporal_method_concordance"),
                    file.path(paths$extended_figures, "figureY_patient_token_coverage", "extended_data_figureY_patient_token_coverage"))
add_check(29, "main and Extended Data exports have PDF/SVG/600-dpi PNG/TIFF stems", all(vapply(c(main_stems, extended_stems), all_formats, logical(1))), paste(c(main_stems, extended_stems)[!vapply(c(main_stems, extended_stems), all_formats, logical(1))], collapse = ";"))
theme_text <- paste(readLines(file.path(root, "scripts", "figure5_plot_theme.R"), warn = FALSE), collapse = "\n")
add_check(30, "export contract requests 600 dpi PNG/TIFF", grepl("dpi = 600", theme_text, fixed = TRUE) && grepl("res = dpi", theme_text, fixed = TRUE), "figure5_plot_theme.R")
dpi_audit_path <- file.path(paths$metadata, "figure5_six_panel_export_dpi_audit.tsv")
dpi_audit <- if (file.exists(dpi_audit_path)) fread(dpi_audit_path) else data.table()
add_check(31, "export DPI audit records 600 dpi raster and vector companions", nrow(dpi_audit) == length(c(main_stems, extended_stems)) &&
            all(dpi_audit$png_dpi_x >= 599 & dpi_audit$png_dpi_y >= 599 & dpi_audit$tiff_dpi_x >= 600 & dpi_audit$tiff_dpi_y >= 600 &
                dpi_audit$pdf_exists & dpi_audit$svg_exists), sprintf("rows=%d", nrow(dpi_audit)))

protected <- fread(file.path(paths$metadata, "figure5_six_panel_protected_current_manifest.tsv"))
add_check(32, "old corrected result artefacts remain byte-identical", nrow(protected) > 0L && all(protected$exists) && all(protected$unchanged), sprintf("%d/%d unchanged", sum(protected$unchanged, na.rm = TRUE), nrow(protected)))

figure5_write_tsv(checks, file.path(paths$metadata, "figure5_six_panel_validation_report.tsv"))
figure5_write_json(list(
  validation = "Figure 5 six-panel refactor",
  checks = split(checks, seq_len(nrow(checks))),
  passed = all(checks$passed),
  n_checks = nrow(checks),
  n_passed = sum(checks$passed)
), file.path(paths$metadata, "figure5_six_panel_validation_report.json"))
if (!all(checks$passed)) {
  print(checks[passed == FALSE])
  stop(sprintf("Figure 5 six-panel validation failed: %d/%d checks passed", sum(checks$passed), nrow(checks)), call. = FALSE)
}
cat(sprintf("Figure 5 six-panel validation passed: %d/%d checks\n", sum(checks$passed), nrow(checks)))
