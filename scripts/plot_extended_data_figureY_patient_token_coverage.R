#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
source(file.path(root, "scripts", "figure5_plot_theme.R"))
suppressPackageStartupMessages(library(data.table))

paths <- figure5_six_panel_paths(root)
source_paths <- figure5_paths(root, create = FALSE)
pseudobulk <- as.data.table(readRDS(file.path(source_paths$processed, "figure5_patient_pseudotime_pseudobulk.rds")))
pseudobulk <- pseudobulk[method == "main/consensus pseudotime"]

make_audit <- function(dt, unit_name, token_col) {
  if (!nrow(dt)) return(data.table())
  dt[, token_tmp := as.character(get(token_col))]
  dt <- dt[!is.na(token_tmp) & nzchar(token_tmp)]
  if (!nrow(dt)) return(data.table())
  out <- dt[, .(
    dataset_id = paste(sort(unique(dataset_id)), collapse = ";"),
    pseudotime_minimum = if (any(is.finite(pseudotime))) min(pseudotime, na.rm = TRUE) else NA_real_,
    pseudotime_maximum = if (any(is.finite(pseudotime))) max(pseudotime, na.rm = TRUE) else NA_real_,
    number_of_bins = uniqueN(pseudotime_bin[is.finite(pseudotime_bin)]),
    n_rows = .N,
    n_cells = sum(n_cells[is.finite(n_cells)], na.rm = TRUE)
  ), by = .(patient_sample_token = token_tmp)]
  out[, aggregation_unit := unit_name]
  out[, coverage_status := fcase(
    is.finite(pseudotime_minimum) & is.finite(pseudotime_maximum) &
      pseudotime_minimum <= 0.10 & pseudotime_maximum >= 0.90 & number_of_bins >= 5,
    "coverage-qualified",
    number_of_bins < 5 & (pseudotime_minimum > 0.10 | pseudotime_maximum < 0.90),
    "insufficient-bins;incomplete-range",
    number_of_bins < 5, "insufficient-bins",
    !is.finite(pseudotime_minimum) | !is.finite(pseudotime_maximum), "no-finite-pseudotime",
    default = "incomplete-range"
  )]
  out[, exclusion_reason := fcase(
    coverage_status == "coverage-qualified", "",
    grepl("insufficient-bins", coverage_status, fixed = TRUE) & grepl("incomplete-range", coverage_status, fixed = TRUE),
    "fewer than 5 pseudotime bins; incomplete pseudotime range (required <=0.10 to >=0.90)",
    grepl("insufficient-bins", coverage_status, fixed = TRUE), "fewer than 5 pseudotime bins",
    coverage_status == "no-finite-pseudotime", "no finite pseudotime values",
    default = "incomplete pseudotime range (required <=0.10 to >=0.90)"
  )]
  out[, `:=`(coverage_min_rule = 0.10, coverage_max_rule = 0.90)]
  out[]
}

patient_dt <- pseudobulk[aggregation_unit == "patient"]
sample_dt <- pseudobulk[aggregation_unit == "sample"]
if ("sample_id" %in% names(sample_dt)) sample_dt[, token_for_audit := sample_id] else sample_dt[, token_for_audit := patient_id]
audit <- rbindlist(list(
  make_audit(patient_dt, "patient", "patient_id"),
  make_audit(sample_dt, "sample", "token_for_audit")
), fill = TRUE)
setorder(audit, aggregation_unit, dataset_id, patient_sample_token)
audit[, token_label := paste0(aggregation_unit, "::", patient_sample_token)]
audit[, token_label := factor(token_label, levels = rev(unique(token_label)))]

coverage_note <- "No individual patient token covered the full pseudotemporal range required for independent within-patient ordering."
p <- ggplot(audit, aes(y = token_label, colour = coverage_status)) +
  geom_vline(xintercept = c(0.10, 0.90), linewidth = 0.45, linetype = "dashed", colour = lancet_palette[8]) +
  geom_segment(data = audit[is.finite(pseudotime_minimum) & is.finite(pseudotime_maximum)],
               aes(x = pseudotime_minimum, xend = pseudotime_maximum, y = token_label, yend = token_label), linewidth = 0.75) +
  geom_point(data = audit[is.finite(pseudotime_minimum)], aes(x = pseudotime_minimum), size = 1.8) +
  geom_point(data = audit[is.finite(pseudotime_maximum)], aes(x = pseudotime_maximum), size = 1.8, shape = 17) +
  facet_grid(aggregation_unit ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c("coverage-qualified" = lancet_palette[3],
                                  "incomplete-range" = lancet_palette[6],
                                  "insufficient-bins" = lancet_palette[5],
                                  "insufficient-bins;incomplete-range" = lancet_palette[8],
                                  "no-finite-pseudotime" = "#FFFFFF"), drop = FALSE) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.10, 0.5, 0.90, 1),
                     labels = c("0", "0.10", "0.50", "0.90", "1")) +
  labs(title = "Extended Data Figure Y  Patient/sample-token pseudotime coverage audit",
       subtitle = coverage_note,
       x = "Unified pseudotime coverage", y = NULL, colour = "Coverage status") +
  theme_figure5(6.8) +
  theme(axis.text.y = element_text(size = 5.7), legend.position = "bottom",
        strip.text = element_text(size = 7.5, face = "bold"), plot.title = element_text(size = 10, face = "bold"))

out_dir <- file.path(paths$extended_figures, "figureY_patient_token_coverage")
outputs <- export_figure5_plot(p, file.path(out_dir, "extended_data_figureY_patient_token_coverage"),
                               9.0, max(5.8, nrow(audit) * 0.18))
figure5_write_tsv(audit[, .(aggregation_unit, dataset_id, patient_sample_token, pseudotime_minimum,
                            pseudotime_maximum, number_of_bins, coverage_status, exclusion_reason,
                            n_rows, n_cells, coverage_min_rule, coverage_max_rule)],
                  file.path(paths$extended_metadata, "extended_data_figureY_patient_token_coverage.tsv"))
figure5_write_json(list(
  figure = "Extended Data Figure Y",
  title = "Patient/sample-token pseudotime coverage audit",
  required_columns = c("patient_sample_token", "pseudotime_minimum", "pseudotime_maximum", "number_of_bins", "coverage_status", "exclusion_reason"),
  body = coverage_note,
  qualification_rule = "pseudotime minimum <=0.10, maximum >=0.90 and >=5 bins",
  n_tokens = nrow(audit),
  n_patient_tokens_coverage_qualified = audit[aggregation_unit == "patient" & coverage_status == "coverage-qualified", .N],
  outputs = as.list(outputs)
), file.path(paths$extended_metadata, "extended_data_figureY_report.json"))
dir.create(paths$extended_processed, recursive = TRUE, showWarnings = FALSE)
saveRDS(p, file.path(paths$extended_processed, "extended_data_figureY_patient_token_coverage.rds"))

