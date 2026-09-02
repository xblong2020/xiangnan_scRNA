#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
paths <- figure5_onset_fix_paths(root)
old_paths <- figure5_paths(root, create = FALSE)

checks <- list()
add_check <- function(name, passed, details) {
  checks[[length(checks) + 1L]] <<- data.table(
    check = name,
    status = if (isTRUE(passed)) "PASS" else "FAIL",
    details = as.character(details)
  )
}

required_metadata <- c(
  "figure5_onset_fix_original_file_md5.tsv",
  "figure5b_gam_model_summary.tsv", "figure5b_gam_predictions.tsv.gz", "figure5b_gam_report.json",
  "figure5c_heatmap_matrix.tsv.gz", "figure5c_temporal_landmarks.tsv",
  "figure5_temporal_landmarks.tsv", "figure5_temporal_landmark_diagnostics.tsv",
  "figure5_bootstrap_temporal_landmarks.tsv.gz", "figure5_bootstrap_temporal_landmark_diagnostics.tsv.gz",
  "figure5_bootstrap_temporal_landmarks_report.json", "figure5e_precedence_probabilities.tsv",
  "figure5f_method_concordance.tsv", "figure5g_patient_temporal_differences.tsv",
  "figure5g_patient_meta_analysis.tsv", "figure5h_phase_boundaries.tsv", "figure5h_activity_band_summary.tsv",
  "figure5_landmarks_old_vs_new.tsv", "figure5_precedence_old_vs_new.tsv",
  "figure5_original_onset_sensitivity.tsv", "figure5_temporal_landmark_qc_failures.tsv",
  "figure5_export_dpi_audit.tsv"
)
missing_metadata <- required_metadata[!file.exists(file.path(paths$metadata, required_metadata))]
add_check("required corrected metadata outputs", !length(missing_metadata),
          if (length(missing_metadata)) paste(missing_metadata, collapse = ", ") else paste(length(required_metadata), "files present"))

boot <- as.data.table(readRDS(file.path(paths$processed, "figure5_bootstrap_temporal_landmarks.rds")))
diagnostics <- fread(file.path(paths$metadata, "figure5_bootstrap_temporal_landmark_diagnostics.tsv.gz"))
add_check("formal bootstrap iteration count", uniqueN(boot$iteration) == 1000L, paste("iterations:", uniqueN(boot$iteration)))
add_check("formal bootstrap landmark row count", nrow(boot) == 1000L * 3L * 8L, paste("rows:", nrow(boot)))
add_check("formal bootstrap diagnostic row count", nrow(diagnostics) == 1000L * 3L, paste("rows:", nrow(diagnostics)))
seed_map <- unique(boot[, .(iteration, iteration_seed, base_seed)])
add_check("actual iteration seeds recorded", nrow(seed_map) == 1000L && all(seed_map$iteration_seed == seed_map$base_seed + seed_map$iteration),
          "iteration_seed equals base_seed + iteration")
add_check("bootstrap model is REML", all(boot$gam_method == "REML"), paste(unique(boot$gam_method), collapse = ", "))
add_check("bootstrap unit is not cell", all(boot$bootstrap_unit %chin% c("patient", "sample")), paste(unique(boot$bootstrap_unit), collapse = ", "))
add_check("coverage failures do not extrapolate", {
  failure_keys <- diagnostics[coverage_ok == FALSE, .(iteration, axis)]
  !nrow(failure_keys) || !nrow(merge(boot[is.finite(time)], failure_keys, by = c("iteration", "axis")))
}, paste("coverage failure fraction:", mean(!diagnostics$coverage_ok)))

main <- fread(file.path(paths$metadata, "figure5_temporal_landmarks.tsv"))
main_wide <- dcast(main, axis ~ landmark, value.var = "time")
add_check("main onset excludes baseline", all(main_wide[is.finite(onset_time), onset_time > 0.10]),
          paste(main_wide$axis, main_wide$onset_time, sep = "=", collapse = "; "))
add_check("main t10 not later than t50", all(main_wide[is.finite(t10) & is.finite(t50), t10 <= t50]),
          paste(main_wide$axis, main_wide$t10, main_wide$t50, sep = ":", collapse = "; "))
add_check("right-boundary maximum slope unresolved", !is.finite(main_wide[axis == "sox4_stabilization", maximum_slope_time]) &&
            diagnostics[axis == "sox4_stabilization", all(maximum_slope_boundary_hit)],
          "SOX4 maximum slope is NA because all formal fits are boundary-limited")

onset_fraction <- boot[landmark == "onset_time", .(finite_fraction = mean(is.finite(time))), by = axis]
add_check("onset stability classified from finite fraction", all(onset_fraction$finite_fraction >= 0.80),
          paste(onset_fraction$axis, onset_fraction$finite_fraction, sep = "=", collapse = "; "))

precedence <- fread(file.path(paths$metadata, "figure5e_precedence_probabilities.tsv"))
add_check("tie-aware precedence counts reconcile", all(precedence$n_earlier + precedence$n_tied + precedence$n_later == precedence$n_valid),
          "earlier + tied + later equals n_valid")
expected_probability <- with(precedence, (n_earlier + 0.5 * n_tied) / n_valid)
probability_ok <- !is.finite(expected_probability) | abs(expected_probability - precedence$probability) < 1e-12
add_check("ties contribute exactly 0.5", all(probability_ok), "probability=(n_earlier+0.5*n_tied)/n_valid")
add_check("precedence exports delta quantiles and Monte Carlo SE",
          all(c("delta_q025", "delta_q975", "Monte_Carlo_SE", "mc_interval_label") %in% names(precedence)) &&
            all(precedence$mc_interval_label == "Monte Carlo interval only"),
          "bootstrap delta interval is separated from Monte Carlo probability error")

methods <- fread(file.path(paths$metadata, "figure5f_method_concordance.tsv"))
add_check("method concordance uses valid-unit counts", all(c("n_units", "n_valid", "evidence_weight") %in% names(methods)),
          "n_units, n_valid and evidence_weight present")
add_check("main and patient-pseudobulk marked non-independent",
          all(grepl("non-independent", methods[method == "patient-pseudobulk", independence_status], fixed = TRUE)),
          "patient-pseudobulk duplicate is explicitly non-independent")
add_check("DPT and CellRank remain unavailable",
          all(methods[method %chin% c("DPT", "CellRank pseudotime"), status == "Not available"]),
          "no unavailable pseudotime method was fabricated")

patient_meta <- fread(file.path(paths$metadata, "figure5g_patient_meta_analysis.tsv"))
add_check("patient-token analysis respects coverage", all(patient_meta$status == "Not available"),
          "no token met both >=5 bins and full coverage; Figure 5G reports Not available")

h_summary <- fread(file.path(paths$metadata, "figure5h_activity_band_summary.tsv"))
h_report <- jsonlite::read_json(file.path(paths$metadata, "figure5h_overlapping_phase_model_report.json"), simplifyVector = TRUE)
caption <- "Activity bands represent relative programme prominence derived from smoothed scores and bootstrap temporal landmarks; they do not indicate discrete activation or termination events."
add_check("Figure 5H exact caption", identical(h_report$caption, caption), h_report$caption)
add_check("Figure 5H does not use old onset", identical(h_report$uses_old_onset_for_formal_band, FALSE),
          "old onset is retained only in sensitivity/source data")
add_check("Figure 5H has no forced fade", all(!h_summary$decline_stable) && all(h_summary$right_edge_rule == "continued prominence to observed pseudotime end"),
          "all three decline landmarks are unstable and all bands continue to the endpoint")
add_check("Figure 5H has no discrete rectangles", !grepl("geom_rect", paste(readLines(file.path(root, "scripts", "plot_figure5h_overlapping_phase_model.R"), warn = FALSE), collapse = "\n"), fixed = TRUE),
          "continuous alpha-gradient polygons used")

panel_dirs <- file.path(paths$figures, c(
  "figure5c_temporal_heatmap", "figure5d_temporal_landmarks", "figure5e_precedence_matrix",
  "figure5f_method_concordance", "figure5g_patient_temporal_forest", "figure5h_overlapping_phase_model"
))
panel_stems <- c(
  "figure5c_temporal_heatmap", "figure5d_temporal_landmarks", "figure5e_precedence_matrix",
  "figure5f_method_concordance", "figure5g_patient_temporal_forest", "figure5h_overlapping_phase_model"
)
panel_files <- unlist(Map(function(directory, stem) file.path(directory, paste0(stem, c(".pdf", ".png", ".svg", ".tiff"))), panel_dirs, panel_stems))
add_check("Figure 5C-H publication formats", all(file.exists(panel_files)) && all(file.info(panel_files)$size > 0),
          paste(length(panel_files), "non-empty panel files"))

dpi_audit <- fread(file.path(paths$metadata, "figure5_export_dpi_audit.tsv"))
add_check("publication PNG resolution", nrow(dpi_audit) >= 6L && all(dpi_audit$dpi_x >= 599) && all(dpi_audit$dpi_y >= 599),
          paste("audited PNG files:", nrow(dpi_audit), "; minimum DPI:", min(dpi_audit$dpi_x, dpi_audit$dpi_y)))

preview_files <- file.path(paths$preview, paste0("figure5_temporal_positioning_onset_fix_a_to_h_preview", c(".pdf", ".png")))
add_check("corrected montage", all(file.exists(preview_files)) && all(file.info(preview_files)$size > 0),
          paste(basename(preview_files), collapse = ", "))

correction_report <- file.path(paths$reports, "figure5_onset_correction_report.md")
report_text <- if (file.exists(correction_report)) paste(readLines(correction_report, warn = FALSE), collapse = "\n") else ""
section_matches <- gregexpr("(?m)^## ([1-9]|1[0-9]|20)\\.", report_text, perl = TRUE)[[1L]]
section_count <- if (identical(section_matches, -1L)) 0L else length(section_matches)
add_check("20-section correction report", file.exists(correction_report) && section_count == 20L,
          paste("numbered sections:", section_count))

old_new_landmarks <- fread(file.path(paths$metadata, "figure5_landmarks_old_vs_new.tsv"))
old_new_precedence <- fread(file.path(paths$metadata, "figure5_precedence_old_vs_new.tsv"))
add_check("old-vs-new landmark comparison complete", nrow(old_new_landmarks) == 24L && all(c("old_estimate", "new_estimate", "interpretation_changed") %in% names(old_new_landmarks)),
          paste("rows:", nrow(old_new_landmarks)))
add_check("old-vs-new precedence comparison complete", nrow(old_new_precedence) == 9L,
          paste("rows:", nrow(old_new_precedence)))

qc <- fread(file.path(paths$metadata, "figure5_temporal_landmark_qc.tsv"))
add_check("formal 15-rule QC", nrow(qc) == 15L && all(qc$status == "PASS"),
          paste("passes:", sum(qc$status == "PASS"), "of", nrow(qc)))
add_check("Figure 1-4 frozen baseline integrity", qc[check_number == 15L, status] == "PASS",
          qc[check_number == 15L, details])

validation <- rbindlist(checks)
figure5_write_tsv(validation, file.path(paths$metadata, "figure5_validation_report.tsv"))
figure5_write_json(list(
  n_checks = nrow(validation),
  n_pass = sum(validation$status == "PASS"),
  n_fail = sum(validation$status == "FAIL"),
  checks = split(validation, seq_len(nrow(validation)))
), file.path(paths$metadata, "figure5_validation_report.json"))
print(validation)
if (any(validation$status == "FAIL")) stop("Corrected Figure 5 validation failed")
