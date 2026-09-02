#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
source_paths <- figure5_paths(root, create = FALSE)
paths <- figure5_onset_fix_paths(root)
cli <- figure5_cli(list(bootstrap_suffix = "development_100"))
suffix <- if (identical(cli$bootstrap_suffix, "formal") || !nzchar(cli$bootstrap_suffix)) "" else paste0("_", cli$bootstrap_suffix)

landmarks <- fread(file.path(paths$metadata, "figure5_temporal_landmarks.tsv"))
main_diagnostics <- fread(file.path(paths$metadata, "figure5_temporal_landmark_diagnostics.tsv"))
boot <- as.data.table(readRDS(file.path(paths$processed, paste0("figure5_bootstrap_temporal_landmarks", suffix, ".rds"))))
boot_diagnostics <- fread(file.path(paths$metadata, paste0("figure5_bootstrap_temporal_landmark_diagnostics", suffix, ".tsv.gz")))

qc <- list()
add_check <- function(number, rule, passed, details, status_if_not_passed = "FAIL") {
  qc[[length(qc) + 1L]] <<- data.table(
    check_number = as.integer(number), rule = rule,
    status = if (isTRUE(passed)) "PASS" else status_if_not_passed,
    details = as.character(details)
  )
}

finite_onset <- landmarks[landmark == "onset_time" & is.finite(time), time]
add_check(1, "onset strictly exceeds baseline_end", all(finite_onset > figure5_temporal_parameters$baseline_end),
          paste("minimum finite onset:", if (length(finite_onset)) min(finite_onset) else "none"))

main_onset <- merge(
  landmarks[landmark == "onset_time", .(axis, onset_time = time)],
  main_diagnostics[, .(axis, onset_found)], by = "axis", all = TRUE
)
add_check(2, "no sustained activation returns NA onset", all(main_onset[onset_found == FALSE, !is.finite(onset_time)]),
          paste("axes without sustained onset:", paste(main_onset[onset_found == FALSE, axis], collapse = ", ")))

rounded_onset <- round(main_onset[is.finite(onset_time), onset_time], 6)
add_check(3, "axes do not receive one common default onset", length(rounded_onset) < 2L || uniqueN(rounded_onset) > 1L,
          paste("main onsets:", paste(main_onset$axis, main_onset$onset_time, sep = "=", collapse = "; ")))

core_text <- paste(readLines(file.path(root, "scripts", "figure5_temporal_core.R"), warn = FALSE), collapse = "\n")
add_check(4, "t10/t50 use baseline-to-peak persistent crossings",
          grepl("baseline_value + 0.10 * diagnostics$total_rise", core_text, fixed = TRUE) &&
            grepl("baseline_value + 0.50 * diagnostics$total_rise", core_text, fixed = TRUE) &&
            grepl("crossing_min_run", core_text, fixed = TRUE),
          "source audit of fixed baseline-to-peak targets and persistent crossing rule")

finite_slope <- landmarks[landmark == "maximum_slope_time" & is.finite(time), time]
add_check(5, "maximum slope lies after baseline window", all(finite_slope > figure5_temporal_parameters$baseline_end),
          paste("finite maximum slopes:", paste(finite_slope, collapse = ", ")))
add_check(6, "maximum slope is not assigned at the right search boundary",
          all(finite_slope < figure5_temporal_parameters$slope_search_end) &&
            all(main_diagnostics[maximum_slope_boundary_hit == TRUE, grepl("maximum_slope_at_search_boundary", failure_reason)]),
          paste("boundary-limited axes set unresolved:", paste(main_diagnostics[maximum_slope_boundary_hit == TRUE, axis], collapse = ", ")))

wide_main <- dcast(landmarks, axis ~ landmark, value.var = "time")
add_check(7, "t10 is not later than t50", all(wide_main[is.finite(t10) & is.finite(t50), t10 <= t50]),
          paste(wide_main$axis, wide_main$t10, wide_main$t50, sep = ":", collapse = "; "))
relationship_ok <- wide_main[, all((!is.finite(t10) | !is.finite(peak_time) | t10 <= peak_time) &
                                     (!is.finite(t50) | !is.finite(peak_time) | t50 <= peak_time) &
                                     (!is.finite(onset_time) | !is.finite(peak_time) | onset_time <= peak_time))]
add_check(8, "onset/t10/t50/peak relationships are audited", relationship_ok,
          "onset may follow t10 when the fixed 0.25 baseline-SD threshold exceeds the 10% rise threshold; all finite starts/crossings must precede peak")

tie_test <- tie_aware_precedence_probability(c(0.3, 0.4), c(0.3, 0.4))
add_check(9, "tied precedence contributes 0.5", identical(tie_test$n_tied, 2L) && tie_test$probability == 0.5,
          paste("test probability:", tie_test$probability))

add_check(10, "bootstrap uses REML", all(boot$gam_method == "REML"),
          paste("methods:", paste(unique(boot$gam_method), collapse = ", ")))
main_script <- paste(readLines(file.path(root, "scripts", "figure5_04_fit_three_axis_gam.R"), warn = FALSE), collapse = "\n")
bootstrap_script <- paste(readLines(file.path(root, "scripts", "figure5_07_bootstrap_temporal_landmarks.R"), warn = FALSE), collapse = "\n")
add_check(11, "main and bootstrap reuse one GAM/landmark structure",
          grepl("landmarks_from_table", main_script, fixed = TRUE) &&
            grepl("landmarks_from_table", bootstrap_script, fixed = TRUE) &&
            !grepl("GCV.Cp", bootstrap_script, fixed = TRUE),
          "both paths call landmarks_from_table/fit_landmark_gam; bootstrap contains no GCV.Cp")

coverage_failures <- boot_diagnostics[coverage_ok == FALSE, .(iteration, axis)]
if (nrow(coverage_failures)) {
  failed_times <- merge(boot, coverage_failures, by = c("iteration", "axis"))[is.finite(time)]
  coverage_no_extrapolation <- nrow(failed_times) == 0L
} else {
  coverage_no_extrapolation <- TRUE
}
add_check(12, "insufficient coverage returns NA without extrapolation", coverage_no_extrapolation,
          paste("coverage-failed axis-iterations:", nrow(coverage_failures)))

profile_grid <- seq(0, 1, length.out = 201)
profile <- build_figure5h_activity_profile(profile_grid, plogis((profile_grid - 0.4) / 0.06),
                                            start = 0.2, t50 = 0.4, maximum_slope = 0.4,
                                            decline_onset = NA_real_)
add_check(13, "Figure 5H has no forced fade without stable decline", all(profile$fade == 1),
          "synthetic monotonic activity retains fade factor 1 through pseudotime end")

plot_files <- file.path(root, "scripts", c(
  "plot_figure5c_temporal_heatmap.R", "plot_figure5d_temporal_landmarks.R",
  "plot_figure5e_precedence_matrix.R", "plot_figure5f_method_concordance.R",
  "plot_figure5g_patient_forest.R", "plot_figure5h_overlapping_phase_model.R"
))
plot_text <- vapply(plot_files, function(path) paste(readLines(path, warn = FALSE), collapse = "\n"), character(1))
theme_text <- paste(readLines(file.path(root, "scripts", "figure5_plot_theme.R"), warn = FALSE), collapse = "\n")
add_check(14, "all revised plots use R/ggplot2 and ggsci Lancet palette",
          all(grepl("figure5_plot_theme.R", plot_text, fixed = TRUE)) &&
            grepl("ggsci::pal_lancet(\"lanonc\")", theme_text, fixed = TRUE),
          "six plot scripts source the shared ggplot2/ggsci Lancet theme")

baseline_path <- file.path(source_paths$metadata, "figure5_figure1_4_baseline.tsv")
baseline <- fread(baseline_path)
existing <- baseline[file.exists(file_path)]
current_md5 <- unname(tools::md5sum(existing$file_path))
changed <- existing[tolower(md5) != tolower(current_md5), relative_path]
missing <- baseline[!file.exists(file_path), relative_path]
add_check(15, "Figure 1-4 and other frozen old results remain unchanged",
          !length(changed) && !length(missing),
          sprintf("checked %d baseline files; changed=%d; missing=%d", nrow(baseline), length(changed), length(missing)))

qc_table <- rbindlist(qc)
failures <- qc_table[status != "PASS"]
figure5_write_tsv(qc_table, file.path(paths$metadata, paste0("figure5_temporal_landmark_qc", suffix, ".tsv")))
figure5_write_tsv(failures, file.path(paths$metadata, "figure5_temporal_landmark_qc_failures.tsv"))
figure5_write_json(list(
  bootstrap_run = if (nzchar(suffix)) substring(suffix, 2L) else "formal",
  n_checks = nrow(qc_table),
  n_pass = sum(qc_table$status == "PASS"),
  n_fail = sum(qc_table$status == "FAIL"),
  n_exception = sum(qc_table$status == "EXCEPTION"),
  checks = split(qc_table, seq_len(nrow(qc_table)))
), file.path(paths$metadata, paste0("figure5_temporal_landmark_qc", suffix, ".json")))

print(qc_table)
if (any(qc_table$status == "FAIL")) stop("Figure 5 onset-fix QC has fatal failures")
