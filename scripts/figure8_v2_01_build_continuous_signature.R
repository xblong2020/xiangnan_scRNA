#!/usr/bin/env Rscript

source(file.path("scripts", "figure8_v2_common.R"))

v1_r_library <- file.path(FIGURE8_V2_ROOT, "data/processed/driver/figure8_transcriptomic_reversal/r_library")
if (dir.exists(v1_r_library)) .libPaths(c(v1_r_library, .libPaths()))

suppressPackageStartupMessages({
  library(Matrix)
  library(hdf5r)
})

figure8_v2_component_weights <- function() {
  c(
    state_component = 0.35,
    trajectory_component = 0.15,
    axis_A_component = 0.10,
    axis_B_component = 0.10,
    axis_C_component = 0.10,
    malignant_state_component = 0.15,
    perturbation_component = 0.05
  )
}

figure8_v2_model_genes <- function(root = FIGURE8_V2_ROOT) {
  if (file.exists("AGENTS.md") && dir.exists("scripts")) root <- "."
  path <- file.path(root, "metadata/driver/figure8_transcriptomic_reversal/figure8_drugreflector_model_genes.tsv")
  model <- fread(path)
  if (!identical(names(model), c("model_gene_order", "gene"))) stop("Unexpected DrugReflector model-gene schema")
  model[, model_gene_order := as.integer(model_gene_order)]
  model[, gene := as.character(gene)]
  if (nrow(model) != 978L || uniqueN(model$gene) != 978L || !identical(model$model_gene_order, seq_len(978L))) {
    stop("Frozen DrugReflector input order is not the required 978 unique landmarks")
  }
  model
}

figure8_v2_combine_components <- function(x, weights = figure8_v2_component_weights()) {
  x <- as.data.table(copy(x))
  missing <- setdiff(names(weights), names(x))
  if (length(missing)) stop("Missing rescue-score components: ", paste(missing, collapse = ", "))
  component_names <- names(weights)
  component_matrix <- as.matrix(x[, ..component_names])
  component_matrix[!is.finite(component_matrix)] <- 0
  component_matrix[component_matrix > 1] <- 1
  component_matrix[component_matrix < -1] <- -1
  signed_sum <- as.numeric(component_matrix %*% weights)
  absolute_mass <- as.numeric(abs(component_matrix) %*% weights)
  agreement <- fifelse(absolute_mass > .Machine$double.eps, abs(signed_sum) / absolute_mass, 0)
  positive_n <- rowSums(component_matrix > 1e-12)
  negative_n <- rowSums(component_matrix < -1e-12)
  x[, `:=`(
    weighted_rescue_mean = signed_sum,
    directional_agreement = agreement,
    final_rescue_vscore = pmax(-1, pmin(1, signed_sum * agreement)),
    conflict_flag = positive_n > 0 & negative_n > 0,
    n_nonzero_components = positive_n + negative_n
  )]
  x[, final_rescue_direction := fifelse(final_rescue_vscore > 0, "up", fifelse(final_rescue_vscore < 0, "down", "zero"))]
  x
}

figure8_v2_axis_balance <- function(x) {
  x <- as.data.table(x)
  cols <- c(axis_A = "axis_A_component", axis_B = "axis_B_component", axis_C = "axis_C_component")
  mass <- vapply(cols, function(nm) sum(abs(x[[nm]]), na.rm = TRUE), numeric(1))
  total <- sum(mass)
  frac <- if (total > 0) mass / total else rep(NA_real_, length(mass))
  entropy <- if (all(is.finite(frac)) && sum(frac > 0)) -sum(frac[frac > 0] * log(frac[frac > 0])) else NA_real_
  effective_axes <- if (is.finite(entropy)) exp(entropy) else NA_real_
  severe <- any(frac > 0.70, na.rm = TRUE) || (is.finite(effective_axes) && effective_axes < 1.8)
  data.table(
    axis = names(cols), component_column = unname(cols), absolute_mass = unname(mass),
    absolute_mass_fraction = unname(frac), effective_axis_number = effective_axes,
    severe_axis_domination = severe,
    domination_rule = "any axis >70% of A/B/C absolute mass or effective-axis number <1.8"
  )
}

figure8_v2_continuous_long <- function(score) {
  score <- as.data.table(score)
  score[final_rescue_direction != "zero", .(
    signature_id = "landmark_continuous_three_axis_rescue_vscore",
    gene,
    desired_direction = final_rescue_direction,
    v_score = final_rescue_vscore,
    axis = fifelse(
      abs(axis_A_component) >= pmax(abs(axis_B_component), abs(axis_C_component)),
      "axis_A_identity",
      fifelse(abs(axis_B_component) >= abs(axis_C_component), "axis_B_stress", "axis_C_sox4")
    ),
    is_landmark = TRUE
  )]
}

figure8_v2_h5_string <- function(h5, path) {
  as.character(h5[[path]][])
}

figure8_v2_h5_obs <- function(h5, name) {
  obj <- h5[[paste0("obs/", name)]]
  encoding <- tryCatch(as.character(obj$attr_open("encoding-type")$read()), error = function(e) "array")
  if (identical(encoding, "categorical")) {
    categories <- as.character(obj[["categories"]][])
    codes <- as.integer(obj[["codes"]][])
    out <- rep(NA_character_, length(codes))
    keep <- codes >= 0L
    out[keep] <- categories[codes[keep] + 1L]
    return(out)
  }
  as.vector(obj[])
}

figure8_v2_read_h5ad_landmarks <- function(path, model_genes) {
  h5 <- H5File$new(path, mode = "r")
  on.exit(h5$close_all(), add = TRUE)
  var_names <- figure8_v2_h5_string(h5, "var/_index")
  selected_position <- match(model_genes, var_names)
  expression_covered <- !is.na(selected_position)

  x_group <- h5[["X"]]
  shape <- as.integer(x_group$attr_open("shape")$read())
  data <- as.numeric(x_group[["data"]][])
  indices <- as.integer(x_group[["indices"]][])
  indptr <- as.integer(x_group[["indptr"]][])
  selected_zero <- selected_position[expression_covered] - 1L
  mapped_col <- match(indices, selected_zero)
  keep <- !is.na(mapped_col)
  row_index <- rep.int(seq_len(shape[[1]]), diff(indptr))
  matrix <- sparseMatrix(
    i = row_index[keep], j = mapped_col[keep], x = data[keep],
    dims = c(shape[[1]], sum(expression_covered)), giveCsparse = TRUE
  )
  colnames(matrix) <- model_genes[expression_covered]
  list(
    matrix = matrix,
    expression_covered = expression_covered,
    obs = data.table(
      cell_id = figure8_v2_h5_string(h5, "obs/_index"),
      trajectory_role = figure8_v2_h5_obs(h5, "trajectory_role"),
      trajectory_root_end_role = figure8_v2_h5_obs(h5, "trajectory_root_end_role"),
      sample_id = figure8_v2_h5_obs(h5, "sample_id"),
      pseudotime = as.numeric(figure8_v2_h5_obs(h5, "driver_main_strict__pseudotime_mean")),
      main_strict_eligible = as.logical(figure8_v2_h5_obs(h5, "driver_main_strict__eligible")),
      malignant_fate_probability = as.numeric(figure8_v2_h5_obs(h5, "cellrank_fate_prob_cnv_supported_malignant"))
    )
  )
}

figure8_v2_sample_means <- function(matrix, samples, minimum_cells = 20L) {
  samples <- as.character(samples)
  sample_levels <- sort(unique(samples[!is.na(samples) & nzchar(samples)]))
  rows <- lapply(sample_levels, function(sample) {
    idx <- which(samples == sample)
    if (length(idx) < minimum_cells) return(NULL)
    means <- as.list(Matrix::colMeans(matrix[idx, , drop = FALSE]))
    as.data.table(c(list(sample_id = sample, n_cells = length(idx)), means))
  })
  rbindlist(rows, fill = TRUE)
}

figure8_v2_correlation_component <- function(matrix, values, desired_sign = -1) {
  keep <- is.finite(values)
  if (sum(keep) < 20L) return(data.table(gene = colnames(matrix), rho = NA_real_, p = NA_real_, q = NA_real_, component = NA_real_))
  dense <- as.matrix(matrix[keep, , drop = FALSE])
  rho <- suppressWarnings(as.numeric(cor(dense, values[keep], method = "spearman", use = "pairwise.complete.obs")))
  n <- sum(keep)
  t_stat <- rho * sqrt(pmax(0, n - 2) / pmax(1e-12, 1 - rho^2))
  p <- 2 * pt(-abs(t_stat), df = n - 2)
  p[!is.finite(p)] <- 1
  q <- p.adjust(p, method = "BH")
  component <- desired_sign * rho * (1 - pmin(1, q))
  component[!is.finite(component)] <- 0
  data.table(gene = colnames(matrix), rho = rho, p = p, q = q, component = component)
}

figure8_v2_axis_evidence <- function(model_genes, root = FIGURE8_V2_ROOT) {
  if (file.exists("AGENTS.md") && dir.exists("scripts")) root <- "."
  result <- data.table(gene = model_genes, axis_A_component = 0, axis_B_component = 0, axis_C_component = 0, malignant_state_component_prior = 0)
  signature <- fread(file.path(root, "metadata/driver/module9_4_drug_reversal_signature.tsv"))
  signature[, normalized_weight := abs(as.numeric(final_weight)) / max(abs(as.numeric(final_weight)), na.rm = TRUE), by = component]
  sig_axis <- signature[, .(
    axis_A_component = max(fifelse(component %in% c("hnf4a_ppara_rescue", "mature_hepatocyte", "tier1_rescue") & desired_direction == "up", normalized_weight, 0), na.rm = TRUE),
    axis_B_component = -max(fifelse(component == "ap1_stress_proliferation" & desired_direction == "down", normalized_weight, 0), na.rm = TRUE),
    axis_C_component = -max(fifelse(component == "sox4_state_specific" & desired_direction == "down", normalized_weight, 0), na.rm = TRUE),
    malignant_state_component_prior = -max(fifelse(component == "c_malignant_like_fate" & desired_direction == "down", normalized_weight, 0), na.rm = TRUE)
  ), by = gene]

  target <- fread(file.path(root, "metadata/driver/module8_tf_target_signature_genes.tsv"))
  target[, target_weight := 1 / log2(pmax(2, as.numeric(rank) + 1))]
  target_axis <- target[axis %in% c("tier1_rescue", "ap1_stress_proliferation", "sox4_state_specific"), .(
    axis_A_component = max(fifelse(axis == "tier1_rescue", target_weight, 0), na.rm = TRUE),
    axis_B_component = -max(fifelse(axis == "ap1_stress_proliferation", target_weight, 0), na.rm = TRUE),
    axis_C_component = -max(fifelse(axis == "sox4_state_specific", target_weight, 0), na.rm = TRUE)
  ), by = gene]

  pathway <- fread(file.path(root, "metadata/driver/module8_pathway_signature_genes.tsv"))
  pathway[, pathway_weight := 1 / log2(pmax(2, as.numeric(term_rank) + 1))]
  pathway_axis <- pathway[axis %in% c("tier1_rescue", "ap1_stress_proliferation", "sox4_state_specific"), .(
    axis_A_component = max(fifelse(axis == "tier1_rescue", pathway_weight, 0), na.rm = TRUE),
    axis_B_component = -max(fifelse(axis == "ap1_stress_proliferation", pathway_weight, 0), na.rm = TRUE),
    axis_C_component = -max(fifelse(axis == "sox4_state_specific", pathway_weight, 0), na.rm = TRUE)
  ), by = gene]

  merge_axis <- function(base, add) {
    z <- merge(base, add, by = "gene", all.x = TRUE, suffixes = c("", ".add"))
    for (nm in c("axis_A_component", "axis_B_component", "axis_C_component")) {
      add_nm <- paste0(nm, ".add")
      if (add_nm %in% names(z)) {
        positive <- nm == "axis_A_component"
        z[, (nm) := if (positive) pmax(get(nm), fifelse(is.na(get(add_nm)), 0, get(add_nm))) else pmin(get(nm), fifelse(is.na(get(add_nm)), 0, get(add_nm)))]
        z[, (add_nm) := NULL]
      }
    }
    z
  }
  result <- merge_axis(result, sig_axis)
  result <- merge_axis(result, target_axis)
  result <- merge_axis(result, pathway_axis)
  result[, `:=`(
    axis_A_component = pmax(0, pmin(1, axis_A_component)),
    axis_B_component = pmax(-1, pmin(0, axis_B_component)),
    axis_C_component = pmax(-1, pmin(0, axis_C_component)),
    malignant_state_component_prior = pmax(-1, pmin(0, malignant_state_component_prior))
  )]
  result
}

figure8_v2_perturbation_component <- function(model_genes, root = FIGURE8_V2_ROOT) {
  if (file.exists("AGENTS.md") && dir.exists("scripts")) root <- "."
  path <- file.path(root, "metadata/driver/celloracle_module6_8_top_gene_delta_by_state.tsv.gz")
  x <- fread(path)
  x <- x[celloracle_state == "malignant_or_malignant_like" & tf %in% c("HNF4A", "PPARA", "HLF", "JUN", "JUNB", "JUND", "FOS", "CEBPB", "EGR1", "ATF3", "SOX4")]
  x[, scaled_delta := figure8_v2_robust_unit(mean_delta_x), by = tf]
  x[, desired_delta := fifelse(tf %in% c("HNF4A", "PPARA", "HLF"), -scaled_delta, scaled_delta)]
  agg <- x[, .(perturbation_component = mean(desired_delta, na.rm = TRUE), perturbation_tfs = paste(sort(unique(tf)), collapse = ";")), by = gene]
  out <- merge(data.table(gene = model_genes), agg, by = "gene", all.x = TRUE)
  out[!is.finite(perturbation_component), perturbation_component := 0]
  out[is.na(perturbation_tfs), perturbation_tfs := ""]
  out[, perturbation_component := pmax(-1, pmin(1, perturbation_component))]
  out
}

figure8_v2_build_signature <- function(root = FIGURE8_V2_ROOT) {
  if (file.exists("AGENTS.md") && dir.exists("scripts")) root <- "."
  set.seed(FIGURE8_V2_SEED)
  figure8_v2_init_dirs()
  model <- figure8_v2_model_genes(root)
  h5ad_path <- file.path(root, "data/processed/driver/driver_union_full_expression.module6_3b.h5ad")
  h5 <- figure8_v2_read_h5ad_landmarks(h5ad_path, model$gene)
  matrix <- h5$matrix
  obs <- h5$obs

  normal_idx <- which(obs$trajectory_role == "normal_reference")
  malignant_idx <- which(obs$trajectory_root_end_role == "end_malignant_cnv" | obs$trajectory_role == "malignant_cnv_supported")
  normal_samples <- figure8_v2_sample_means(matrix[normal_idx, , drop = FALSE], obs$sample_id[normal_idx])
  malignant_samples <- figure8_v2_sample_means(matrix[malignant_idx, , drop = FALSE], obs$sample_id[malignant_idx])
  gene_cols <- colnames(matrix)
  normal_mean <- if (nrow(normal_samples)) colMeans(as.matrix(normal_samples[, ..gene_cols])) else rep(NA_real_, length(gene_cols))
  malignant_mean <- if (nrow(malignant_samples)) colMeans(as.matrix(malignant_samples[, ..gene_cols])) else rep(NA_real_, length(gene_cols))
  state_effect <- normal_mean - malignant_mean
  state_component <- figure8_v2_robust_unit(state_effect)

  eligible <- obs$main_strict_eligible & is.finite(obs$pseudotime)
  trajectory <- figure8_v2_correlation_component(matrix[eligible, , drop = FALSE], obs$pseudotime[eligible], desired_sign = -1)
  fate <- figure8_v2_correlation_component(matrix, obs$malignant_fate_probability, desired_sign = -1)

  expr_mean <- Matrix::colMeans(matrix)
  expr_sq_mean <- Matrix::colMeans(matrix^2)
  expr_variance <- pmax(0, expr_sq_mean - expr_mean^2)
  detection_rate <- Matrix::colMeans(matrix > 0)

  expression <- data.table(
    gene = gene_cols, expression_covered = TRUE, mean_expression = expr_mean,
    expression_variance = expr_variance, detection_rate = detection_rate,
    state_effect_normal_minus_malignant = state_effect, state_component = state_component
  )
  expression <- merge(expression, trajectory[, .(gene, trajectory_rho = rho, trajectory_q = q, trajectory_component = component)], by = "gene", all.x = TRUE)
  expression <- merge(expression, fate[, .(gene, malignant_fate_rho = rho, malignant_fate_q = q, malignant_state_component_data = component)], by = "gene", all.x = TRUE)

  score <- merge(model, expression, by = "gene", all.x = TRUE, sort = FALSE)
  setorder(score, model_gene_order)
  score[, expression_covered := !is.na(expression_covered) & expression_covered]
  for (nm in c("mean_expression", "expression_variance", "detection_rate")) score[!is.finite(get(nm)), (nm) := 0]

  axis <- figure8_v2_axis_evidence(model$gene, root)
  perturb <- figure8_v2_perturbation_component(model$gene, root)
  score <- merge(score, axis, by = "gene", all.x = TRUE, sort = FALSE)
  score <- merge(score, perturb, by = "gene", all.x = TRUE, sort = FALSE)
  setorder(score, model_gene_order)
  for (nm in c("state_component", "trajectory_component", "axis_A_component", "axis_B_component", "axis_C_component", "malignant_state_component_data", "malignant_state_component_prior", "perturbation_component")) {
    score[!is.finite(get(nm)), (nm) := 0]
  }
  score[, malignant_state_component := pmin(malignant_state_component_data, malignant_state_component_prior)]
  score <- figure8_v2_combine_components(score)

  v1 <- fread(file.path(root, "metadata/driver/module9_4_drug_reversal_signature.tsv"))
  qc <- v1[, .(qc_flag = any(housekeeping_or_qc_flag), v1_signature_member = TRUE), by = gene]
  score <- merge(score, qc, by = "gene", all.x = TRUE, sort = FALSE)
  setorder(score, model_gene_order)
  score[is.na(qc_flag), qc_flag := grepl("^MT-|^RPL|^RPS|^HIST|MALAT1$|NEAT1$", gene)]
  score[is.na(v1_signature_member), v1_signature_member := FALSE]
  score[, regulatory_sign := sign(axis_A_component + axis_B_component + axis_C_component)]
  component_cols <- names(figure8_v2_component_weights())
  score[, evidence_sources := apply(.SD, 1, function(z) paste(sub("_component$", "", component_cols[abs(as.numeric(z)) > 1e-12]), collapse = ";")), .SDcols = component_cols]
  score[, landmark_status := "DrugReflector_frozen_model_landmark"]

  required <- c(
    "gene", "landmark_status", "axis_A_component", "axis_B_component", "axis_C_component",
    "malignant_state_component", "perturbation_component", "regulatory_sign", "final_rescue_direction",
    "final_rescue_vscore", "evidence_sources", "conflict_flag", "qc_flag"
  )
  setcolorder(score, c("model_gene_order", required, setdiff(names(score), c("model_gene_order", required))))
  figure8_v2_write_tsv(score, "figure8_v2_gene_level_rescue_vscore.tsv")

  coverage <- score[, .(
    model_gene_order, gene, expression_covered, usable_rescue_score = abs(final_rescue_vscore) > 0,
    score_sign = final_rescue_direction, final_rescue_vscore, mean_expression, expression_variance,
    detection_rate, axis_A_component, axis_B_component, axis_C_component,
    malignant_state_component, perturbation_component, conflict_flag, qc_flag
  )]
  figure8_v2_write_tsv(coverage, "figure8_v2_978_landmark_coverage.tsv")
  balance <- figure8_v2_axis_balance(score)
  figure8_v2_write_tsv(balance, "figure8_v2_landmark_axis_balance.tsv")

  v1_manifest <- fread(file.path(root, "metadata/driver/figure8_transcriptomic_reversal/figure8_signature_variant_manifest.tsv"))
  continuous_row <- data.table(
    signature_id = "landmark_continuous_three_axis_rescue_vscore",
    signature_label = "Continuous landmark-space three-axis rescue v-score",
    up_gene_count = score[final_rescue_vscore > 0, .N],
    down_gene_count = score[final_rescue_vscore < 0, .N],
    landmark_up_count = score[final_rescue_vscore > 0, .N],
    landmark_down_count = score[final_rescue_vscore < 0, .N],
    component_weights = paste(names(figure8_v2_component_weights()), figure8_v2_component_weights(), sep = "=", collapse = ";"),
    conflict_rule = "agreement-shrunk continuous score; unsupported coordinates zero",
    normalization = "component robust scaling to [-1,1]; frozen weighted agreement shrinkage",
    model_input_dimension = 978L,
    primary_or_sensitivity = "v2_primary"
  )
  manifest <- rbindlist(list(v1_manifest, continuous_row), fill = TRUE, use.names = TRUE)
  figure8_v2_write_tsv(manifest, "figure8_v2_signature_variant_manifest.tsv")

  v1_long <- fread(file.path(root, "metadata/driver/figure8_transcriptomic_reversal/figure8_signature_variants_long.tsv.gz"))
  continuous_long <- figure8_v2_continuous_long(score)
  all_long <- rbindlist(list(v1_long, continuous_long), fill = TRUE, use.names = TRUE)
  figure8_v2_write_tsv(all_long, "figure8_v2_signature_variants_long.tsv.gz", compress = TRUE)
  model_long <- all_long[gene %in% model$gene, .(v_score = sum(as.numeric(v_score), na.rm = TRUE)), by = .(signature_id, gene)]
  wide <- dcast(model_long, signature_id ~ gene, value.var = "v_score", fill = 0)
  missing_model <- setdiff(model$gene, names(wide))
  for (gene in missing_model) wide[, (gene) := 0]
  setcolorder(wide, c("signature_id", model$gene))
  figure8_v2_write_tsv(wide, "figure8_v2_signature_variants_wide.tsv.gz", directory = FIGURE8_V2_DATA, compress = TRUE)

  report_lines <- c(
    "# Figure 8 v2 landmark input QC", "",
    paste0("- Frozen landmark universe: ", nrow(score)),
    paste0("- Expression-covered landmarks: ", sum(score$expression_covered)),
    paste0("- Non-zero rescue scores: ", sum(abs(score$final_rescue_vscore) > 0)),
    paste0("- Positive/negative/zero: ", sum(score$final_rescue_vscore > 0), "/", sum(score$final_rescue_vscore < 0), "/", sum(score$final_rescue_vscore == 0)),
    paste0("- v1 sparse reference: 47/300 (15.7%)"),
    paste0("- Effective axis number: ", format(balance$effective_axis_number[[1]], digits = 4)),
    paste0("- Severe axis domination: ", balance$severe_axis_domination[[1]]), "",
    "The 25% v1 review-risk threshold is not used as a v2 target. Unsupported model coordinates remain zero and no landmark was added manually."
  )
  writeLines(report_lines, file.path(FIGURE8_V2_METADATA, "figure8_v2_landmark_input_qc_report.md"), useBytes = TRUE)
  figure8_v2_write_json(list(
    module = "figure8_v2_continuous_signature", seed = FIGURE8_V2_SEED,
    source_h5ad = "data/processed/driver/driver_union_full_expression.module6_3b.h5ad",
    normal_cells = length(normal_idx), malignant_cells = length(malignant_idx),
    normal_samples = nrow(normal_samples), malignant_samples = nrow(malignant_samples),
    model_landmarks = nrow(score), expression_covered = sum(score$expression_covered),
    usable_nonzero = sum(abs(score$final_rescue_vscore) > 0), axis_balance = balance,
    component_weights = as.list(figure8_v2_component_weights()),
    compound_informed_tuning = FALSE, outcome_informed_selection = FALSE
  ), "figure8_v2_continuous_signature_report.json")
  invisible(score)
}

if (sys.nframe() == 0L && Sys.getenv("FIGURE8_V2_TEST_MODE") != "1") {
  result <- figure8_v2_build_signature()
  cat("FIGURE8_V2_SIGNATURE landmarks=", nrow(result), " nonzero=", sum(abs(result$final_rescue_vscore) > 0), "\n", sep = "")
}
