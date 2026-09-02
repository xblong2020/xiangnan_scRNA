#!/usr/bin/env Rscript

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figure6_common.R"))

delta_path <- file.path(FIGURE6_METADATA_DIR, "figure6_celloracle_programme_deltas_by_cell.tsv.gz")
shift_path <- file.path(FIGURE6_PROJECT_ROOT, "metadata", "driver", "celloracle_module6_8_cell_shift_summary.tsv.gz")
if (!file.exists(delta_path)) stop("Missing raw CellOracle programme export: ", delta_path)

delta <- figure6_fread(delta_path)
shift <- figure6_fread(shift_path)[, .(tf, cell_id, malignant_axis_projection)]
delta <- merge(delta, shift, by = c("tf", "cell_id"), all.x = TRUE)
delta <- delta[as.character(celloracle_main_strict) %in% c("TRUE", "T", "1", "true")]
delta[, dataset := as.character(dataset)]
delta[, sample_id_fig6 := fifelse(!is.na(cnv_sample) & nzchar(cnv_sample), as.character(cnv_sample),
  fifelse(!is.na(study_sample) & nzchar(study_sample), as.character(study_sample), as.character(sample_id)))]

programme_cols <- setdiff(FIGURE6_CORE_OUTPUTS, "malignant_fate_change")
baseline_cols <- paste0("baseline_", programme_cols)
missing_cols <- setdiff(c(programme_cols, baseline_cols), names(delta))
if (length(missing_cols)) stop("Programme export lacks: ", paste(missing_cols, collapse = ", "))

long_parts <- lapply(programme_cols, function(outcome) {
  baseline <- paste0("baseline_", outcome)
  tmp <- delta[, .(
    n_cells = .N,
    raw_mean_change = mean(get(outcome), na.rm = TRUE),
    baseline_sd = stats::sd(get(baseline), na.rm = TRUE),
    baseline_mean = mean(get(baseline), na.rm = TRUE),
    baseline_stress = mean(driver_main_strict__module_Stressed_Injured, na.rm = TRUE),
    baseline_proliferation = mean(driver_main_strict__module_Proliferation, na.rm = TRUE),
    baseline_cnv_proxy = mean(cnv_proxy_z, na.rm = TRUE)
  ), by = .(tf, dataset, sample_id = sample_id_fig6)]
  fallback_sd <- stats::sd(delta[[baseline]], na.rm = TRUE)
  tmp[!is.finite(baseline_sd) | baseline_sd <= 0, baseline_sd := fallback_sd]
  tmp[, standardized_effect := raw_mean_change / baseline_sd]
  tmp[, output := outcome]
  tmp
})

fate <- delta[, .(
  n_cells = .N,
  raw_mean_change = mean(malignant_axis_projection, na.rm = TRUE),
  projection_sd = stats::sd(malignant_axis_projection, na.rm = TRUE),
  baseline_stress = mean(driver_main_strict__module_Stressed_Injured, na.rm = TRUE),
  baseline_proliferation = mean(driver_main_strict__module_Proliferation, na.rm = TRUE),
  baseline_cnv_proxy = mean(cnv_proxy_z, na.rm = TRUE)
), by = .(tf, dataset, sample_id = sample_id_fig6)]
fallback_fate_sd <- stats::sd(delta$malignant_axis_projection, na.rm = TRUE)
fate[!is.finite(projection_sd) | projection_sd <= 0, projection_sd := fallback_fate_sd]
fate[, `:=`(
  standardized_effect = raw_mean_change / projection_sd,
  baseline_sd = projection_sd,
  baseline_mean = NA_real_,
  output = "malignant_fate_change"
)]

sample_effects <- rbindlist(c(long_parts, list(fate)), fill = TRUE)
sample_effects[, `:=`(
  perturbation = paste(tf, "KO"),
  axis = figure6_axis_for_tf(tf),
  response_availability = "available",
  effect_interpretation = fifelse(output == "cnv_malignant_signature_change",
    "CNV-associated malignant expression signature change", "signed predicted change")
)]

ap1 <- sample_effects[tf %in% FIGURE6_AP1_MEMBERS,
  .(
    n_cells = max(n_cells), raw_mean_change = stats::median(raw_mean_change, na.rm = TRUE),
    baseline_sd = stats::median(baseline_sd, na.rm = TRUE), baseline_mean = stats::median(baseline_mean, na.rm = TRUE),
    standardized_effect = stats::median(standardized_effect, na.rm = TRUE),
    baseline_stress = stats::median(baseline_stress, na.rm = TRUE),
    baseline_proliferation = stats::median(baseline_proliferation, na.rm = TRUE),
    baseline_cnv_proxy = stats::median(baseline_cnv_proxy, na.rm = TRUE)
  ), by = .(dataset, sample_id, output)]
ap1[, `:=`(
  tf = "AP1_AGGREGATE", perturbation = "AP-1 member-KO aggregate", axis = "stress_axis",
  response_availability = "available", effect_interpretation = "median across JUN/JUNB/JUND/FOS/ATF3 KO sample effects"
)]
sample_effects <- rbindlist(list(sample_effects, ap1), fill = TRUE)
figure6_fwrite(sample_effects, file.path(FIGURE6_METADATA_DIR, "figure6_sample_level_effects.tsv.gz"))

rows <- unique(sample_effects[, .(tf, perturbation, axis, output)])
effects <- rows[, {
  z <- sample_effects[tf == .BY$tf & output == .BY$output]
  b <- figure6_stratified_bootstrap(z, "standardized_effect", n_boot = 1000L, seed = 20260805L + .GRP)
  dataset_est <- z[, .(x = mean(standardized_effect, na.rm = TRUE)), by = dataset]$x
  list(
    effect_estimate = b$estimate, ci_low = b$ci_low, ci_high = b$ci_high, pvalue = b$pvalue,
    sign = sign(b$estimate), n_cells = sum(z$n_cells), n_samples = uniqueN(z$sample_id),
    n_datasets = uniqueN(z$dataset), stability = mean(sign(dataset_est) == sign(b$estimate), na.rm = TRUE)
  )
}, by = .(tf, perturbation, axis, output)]
effects[, fdr := p.adjust(pvalue, method = "BH")]
effects[, `:=`(method = "CellOracle KO; R sample-level dataset-stratified bootstrap", state_subset = "global strict-main")]

not_available <- CJ(
  perturbation = c("HNF4A restore/OE", "PPARA restore/OE", "SOX4 OE"),
  output = FIGURE6_CORE_OUTPUTS
)
not_available[, `:=`(
  tf = sub(" .*", "", perturbation),
  axis = c(HNF4A = "identity_axis", PPARA = "identity_axis", SOX4 = "sox4_axis")[sub(" .*", "", perturbation)],
  effect_estimate = NA_real_, ci_low = NA_real_, ci_high = NA_real_, pvalue = NA_real_, fdr = NA_real_, sign = NA_real_,
  n_cells = 0L, n_samples = 0L, n_datasets = 0L, stability = NA_real_, method = "not available",
  state_subset = "not available", availability = "Not available"
)]
effects[, availability := "Available"]
effects <- rbindlist(list(effects, not_available), fill = TRUE, use.names = TRUE)
setcolorder(effects, c("perturbation", "tf", "axis", "output", setdiff(names(effects), c("perturbation", "tf", "axis", "output"))))
figure6_fwrite(effects, file.path(FIGURE6_METADATA_DIR, "figure6_perturbation_response_effects.tsv.gz"))

report <- list(
  analysis = "Figure 6 perturbation-response matrix",
  n_cell_rows = nrow(delta), n_sample_effect_rows = nrow(sample_effects), n_effect_rows = nrow(effects),
  bootstrap = list(unit = "sample", n = 1000, stratification = "dataset", seed = 20260805),
  programme_standardization = "sample mean predicted delta expression divided by within-sample baseline programme SD; global SD fallback",
  malignant_fate_definition = "sample mean CellOracle KO vector projection onto frozen malignant axis divided by within-sample projection SD",
  ap1_definition = "median sample effect across JUN, JUNB, JUND, FOS and ATF3 individual KO simulations",
  unavailable = c("HNF4A restore/OE", "PPARA restore/OE", "SOX4 OE"),
  guardrails = c("No cell-level p-values", "CNV output is an expression signature", "No KO was relabelled as restoration")
)
figure6_write_json(report, file.path(FIGURE6_METADATA_DIR, "figure6_perturbation_response_report.json"))
