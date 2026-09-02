#!/usr/bin/env Rscript

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figure6_common.R"))

s <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6_sample_level_effects.tsv.gz"))
defs <- data.table(
  comparison = c("Axis A → Axis B", "Axis A → Axis C", "Axis B → Axis C"),
  source_axis = c("identity_axis", "identity_axis", "stress_axis"),
  target_axis = c("stress_axis", "sox4_axis", "sox4_axis"),
  forward_tfs = c("HNF4A;PPARA", "HNF4A;PPARA", "EGR1;CEBPB;AP1_AGGREGATE"),
  reverse_tfs = c("EGR1;CEBPB;AP1_AGGREGATE", "SOX4", "SOX4"),
  forward_output = c("stress_transition_change", "sox4_programme_change", "sox4_programme_change"),
  reverse_output = c("identity_program_change", "identity_program_change", "stress_transition_change")
)
sample_pairs <- list(); results <- list()
for (i in seq_len(nrow(defs))) {
  d <- defs[i]
  ft <- strsplit(d$forward_tfs, ";", fixed = TRUE)[[1]]
  rt <- strsplit(d$reverse_tfs, ";", fixed = TRUE)[[1]]
  f <- s[tf %in% ft & output == d$forward_output,
    .(forward_signed_effect = mean(standardized_effect, na.rm = TRUE)), by = .(dataset, sample_id)]
  r <- s[tf %in% rt & output == d$reverse_output,
    .(reverse_signed_effect = mean(standardized_effect, na.rm = TRUE)), by = .(dataset, sample_id)]
  z <- merge(f, r, by = c("dataset", "sample_id"))
  z[, `:=`(
    comparison = d$comparison, source_axis = d$source_axis, target_axis = d$target_axis,
    forward_absolute_effect = abs(forward_signed_effect), reverse_absolute_effect = abs(reverse_signed_effect)
  )]
  z[, directional_asymmetry_score := forward_absolute_effect - reverse_absolute_effect]
  sample_pairs[[i]] <- z
  b <- figure6_stratified_bootstrap(z, "directional_asymmetry_score", 1000L, 20260825L + i)
  f_est <- figure6_dataset_balanced_mean(z$forward_signed_effect, z$dataset)
  r_est <- figure6_dataset_balanced_mean(z$reverse_signed_effect, z$dataset)
  dataset_asym <- z[, .(x = mean(directional_asymmetry_score, na.rm = TRUE)), by = dataset]$x
  results[[i]] <- data.table(
    comparison = d$comparison, source_axis = d$source_axis, target_axis = d$target_axis,
    forward_signed_effect = f_est, reverse_signed_effect = r_est,
    forward_absolute_effect = abs(f_est), reverse_absolute_effect = abs(r_est),
    directional_asymmetry_score = b$estimate, ci_low = b$ci_low, ci_high = b$ci_high, pvalue = b$pvalue,
    n_samples = uniqueN(z$sample_id), n_datasets = uniqueN(z$dataset),
    dataset_consistency = mean(sign(dataset_asym) == sign(b$estimate), na.rm = TRUE)
  )
}
sample_pairs <- rbindlist(sample_pairs)
out <- rbindlist(results)
out[, fdr := p.adjust(pvalue, "BH")]
out[, classification := fifelse(ci_low > 0, "Forward-dominant", fifelse(ci_high < 0, "Reverse-dominant", "Symmetric/unresolved"))]
figure6_fwrite(sample_pairs, file.path(FIGURE6_METADATA_DIR, "figure6c_directional_asymmetry_sample_data.tsv.gz"))
figure6_fwrite(out, file.path(FIGURE6_METADATA_DIR, "figure6c_directional_asymmetry.tsv"), compress = FALSE)
figure6_write_json(list(
  panel = "Figure 6C", definition = "abs(forward signed effect) - abs(reverse signed effect)",
  inference = "1000 sample-level bootstrap replicates stratified by dataset", comparisons = nrow(out),
  guardrail = "Directional asymmetry score is descriptive computational evidence, not a causal index or epistasis score."
), file.path(FIGURE6_METADATA_DIR, "figure6c_directional_asymmetry_report.json"))

