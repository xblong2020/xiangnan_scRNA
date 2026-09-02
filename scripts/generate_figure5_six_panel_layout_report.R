#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure5_temporal_core.R"))
suppressPackageStartupMessages(library(data.table))

paths <- figure5_six_panel_paths(root)

mapping <- data.table(
  old_panel = c("5A", "5B", "5C", "5D", "5E", "5F", "5G", "5H"),
  new_panel = c("5A", "5B", "5C", "5D", "5E", "Extended Data Figure X", "Extended Data Figure Y", "5F"),
  status = c("redesigned", "retained/refactored", "retained/refactored", "retained/refactored", "retained/refactored", "moved", "moved", "renumbered"),
  output = c(
    "panel_A_corrected_workflow/figure5_six_panel_A_corrected_workflow",
    "panel_B_programmes/figure5_six_panel_B_programmes",
    "panel_C_programme_heatmap/figure5_six_panel_C_programme_heatmap",
    "panel_D_landmarks/figure5_six_panel_D_landmarks",
    "panel_E_precedence/figure5_six_panel_E_precedence",
    "extended_data/figureX_temporal_method_concordance/extended_data_figureX_temporal_method_concordance",
    "extended_data/figureY_patient_token_coverage/extended_data_figureY_patient_token_coverage",
    "panel_F_activity_model/figure5_six_panel_F_activity_model"
  ),
  notes = c(
    "formal left-to-right frozen-programme workflow; no independent pseudotime methods or pre-drawn phase bands",
    "coverage-corrected primary programme GAM",
    "frozen TF, regulon and target-gene programme organization",
    "peak retained here only",
    "corrected onset, t10, t50 and maximum slope; Pr(earlier)",
    "all NA, small-n and Not available method/resampling rows retained",
    "patient/sample-token coverage audit with exclusion reasons",
    "t10-priority continuous activity bands; onset retained only in source/sensitivity data"
  )
)
figure5_write_tsv(mapping, file.path(paths$metadata, "figure5_six_panel_layout_old_vs_new.tsv"))

baseline_path <- file.path(paths$metadata, "figure5_six_panel_protected_original_manifest.tsv")
if (!file.exists(baseline_path)) stop("Protected original manifest is missing", call. = FALSE)
baseline <- fread(baseline_path)
current <- rbindlist(lapply(seq_len(nrow(baseline)), function(i) {
  path <- baseline$file_path[[i]]
  exists <- file.exists(path)
  data.table(file_path = path, exists = exists,
             size = if (exists) as.numeric(file.info(path)$size) else NA_real_,
             md5 = if (exists) unname(tools::md5sum(path)) else NA_character_,
             unchanged = exists && identical(unname(tools::md5sum(path)), baseline$md5[[i]]))
}), fill = TRUE)
figure5_write_tsv(current, file.path(paths$metadata, "figure5_six_panel_protected_current_manifest.tsv"))

caption <- "Activity bands represent relative programme prominence derived from smoothed scores and bootstrap temporal landmarks; they do not indicate discrete activation or termination events."
fallback_note <- "Because eligible patient-token data did not span the complete pseudotemporal range, the primary bootstrap used the prespecified sample-token coverage fallback."
state_note <- "State labels provide descriptive anchors and do not define discrete temporal transitions."
coverage_note <- "No individual patient token covered the full pseudotemporal range required for independent within-patient ordering."

md <- c(
  "# Figure 5 six-panel refactor and old-versus-new layout report",
  "",
  "## Final main-figure title",
  "",
  "**Corrected temporal positioning reveals overlapping regulatory programmes with later SOX4 prominence.**",
  "",
  "## Main Figure 5 panel map",
  "",
  "| Panel | Final content | Source-data namespace |",
  "|---|---|---|",
  "| 5A | Corrected temporal-positioning analysis workflow | `metadata/driver/figure5_temporal_positioning_six_panel/` |",
  "| 5B | Coverage-corrected regulatory programmes along unified pseudotime | six-panel namespace |",
  "| 5C | Temporal organization of frozen TF, regulon and target-gene programmes | six-panel namespace |",
  "| 5D | Corrected bootstrap temporal landmarks | six-panel namespace |",
  "| 5E | Tie-aware temporal precedence across complementary landmarks | six-panel namespace |",
  "| 5F | Conservative overlapping regulatory-activity model (former corrected 5H) | six-panel namespace |",
  "",
  "The layout is: row 1 = 5A full width; row 2 = 5B + 5C; row 3 = 5D + 5E; row 4 = 5F full width.",
  "",
  "## Figure 5A workflow contract",
  "",
  "The redesigned workflow is left-to-right: frozen hepatocyte atlas and metadata; frozen HNF4A/PPARA identity-loss, AP-1/CEBPB/EGR1 stress-transition and SOX4 malignant-state-stabilization programmes; oriented consensus pseudotime (0 = normal/reference, 1 = CNV-supported malignant/malignant-like); CNV-supported malignant endpoint and CellRank malignant-fate probability; coverage-qualified patient/sample-token pseudobulk GAM; corrected onset/t10/t50/maximum-slope landmarks; dataset-stratified bootstrap; tie-aware precedence probability; conservative overlapping regulatory-activity model.",
  "",
  paste0("State-label note: ", state_note),
  "",
  "Sensitivity and coverage audits—Extended Data: Monocle3, Slingshot, CytoTRACE2, LODO, LOSO, CNV-strict, no-proliferation and no-generic-stress are Extended Data audits rather than independent primary evidence in 5A.",
  "",
  paste0("Coverage note: ", fallback_note),
  "",
  "## Landmark and precedence policy",
  "",
  "Figure 5D retains corrected onset, t10, t50, maximum slope and peak. Figure 5E reports corrected onset, t10, t50 and maximum slope using `Pr(earlier)`, with tie-aware bootstrap probabilities; peak is excluded because no prespecified interior-stable peak audit passed.",
  "",
  "Formal activity-band starts in 5F use stable bootstrap t10 first. A fallback is accepted only when the first persistent five-grid-point run has the prespecified biological derivative direction, derivative 95% CI excluding zero and at least 0.25 SD change from the early baseline. Otherwise the label is `boundary unresolved`. Corrected onset remains in source-data and sensitivity outputs and is excluded from the formal band boundary.",
  "",
  paste0("Figure 5F caption: ", caption),
  "",
  "## Extended Data relocation",
  "",
  "- Extended Data Figure X: **Availability and sensitivity of temporal-order estimates across trajectory methods and resampling schemes**. All NA, small-n and Not available rows are retained.",
  paste0("- Extended Data Figure Y: patient/sample-token coverage audit with pseudotime minimum, pseudotime maximum, number of bins, coverage status and exclusion reason. ", coverage_note),
  "",
  "## Results and source-data references",
  "",
  "Results text should cite the six-panel main figure as Figure 5A–F. Cross-method availability/sensitivity is cited as Extended Data Figure X; patient/sample-token coverage is cited as Extended Data Figure Y. The former main-panel 5F/5G claims are not used as independent primary evidence.",
  "",
  "The verbatim Results/legend replacement text is also recorded in `reports/figure5_six_panel_results_reference_update.md`.",
  "",
  "All new source-data files use the `figure5_six_panel_*` or `extended_data_figure*` namespace. The previous corrected onset-fix namespace is read-only and remains available for sensitivity/source-data traceability.",
  "",
  "## Protected-result audit",
  "",
  sprintf("The protected manifest contains %d pre-existing corrected artefacts; %d remain byte-identical after six-panel generation.", nrow(baseline), sum(current$unchanged, na.rm = TRUE)),
  "",
  "See `metadata/driver/figure5_temporal_positioning_six_panel/figure5_six_panel_layout_old_vs_new.tsv` and the validation report for machine-readable checks.",
  ""
)
report_path <- file.path(paths$reports, "figure5_six_panel_layout_report.md")
writeLines(md, report_path, useBytes = TRUE)
figure5_write_json(list(
  title = "Corrected temporal positioning reveals overlapping regulatory programmes with later SOX4 prominence.",
  mapping = split(mapping, seq_len(nrow(mapping))),
  caption = caption,
  fallback_note = fallback_note,
  state_note = state_note,
  extended_data = list(figureX = "Availability and sensitivity of temporal-order estimates across trajectory methods and resampling schemes",
                       figureY = coverage_note),
  protected_artifact_count = nrow(baseline),
  protected_artifact_unchanged = sum(current$unchanged, na.rm = TRUE)
), file.path(paths$metadata, "figure5_six_panel_layout_report.json"))
cat(sprintf("Wrote %s\n", report_path))
