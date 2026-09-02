#!/usr/bin/env Rscript

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figure6_common.R"))
suppressPackageStartupMessages({library(lavaan); library(rsample)})

raw <- figure6_fread(file.path(FIGURE6_PROJECT_ROOT, "metadata", "driver", "module9_3_scrna_pseudobulk_axis_scores.tsv.gz"))
dat <- raw[run_id == "main_strict", .(
  A = mean(A_hnf4a_ppara_loss, na.rm = TRUE), B = mean(B_transition_activation, na.rm = TRUE),
  C = mean(C_sox4_axis, na.rm = TRUE), fate = mean(mean_C_malignant_like_fate, na.rm = TRUE),
  n_cells = max(n_cells, na.rm = TRUE)
), by = .(dataset, cnv_sample)]
dat[, c("A", "B", "C", "fate") := lapply(.SD, figure6_safe_scale), .SDcols = c("A", "B", "C", "fate")]
dat <- dat[complete.cases(dat[, .(A, B, C, fate)])]
figure6_fwrite(dat, file.path(FIGURE6_METADATA_DIR, "figure6f_model_input_sample_pseudobulk.tsv"), compress = FALSE)

models <- list(
  "Model 1: Linear cascade" = "B ~ A\nC ~ B\nfate ~ C",
  "Model 2: Parallel convergence" = "fate ~ A + B + C\nA ~~ B\nA ~~ C\nB ~~ C",
  "Model 3: Partially ordered" = "B ~ A\nC ~ A + B\nfate ~ A + B + C"
)
fits <- lapply(models, function(syntax) lavaan::sem(syntax, data = dat, meanstructure = TRUE, fixed.x = FALSE, estimator = "ML"))
metric_names <- c("chisq", "df", "pvalue", "aic", "bic", "cfi", "tli", "rmsea", "srmr")
fit_summary <- rbindlist(lapply(seq_along(fits), function(i) {
  f <- fits[[i]]; m <- suppressWarnings(lavaan::fitMeasures(f, metric_names))
  as.data.table(as.list(m))[, `:=`(model = names(fits)[i], converged = lavaan::inspect(f, "converged"), n_samples = nrow(dat))]
}), fill = TRUE)
params <- rbindlist(lapply(seq_along(fits), function(i) {
  p <- lavaan::parameterEstimates(fits[[i]], standardized = TRUE, ci = TRUE)
  as.data.table(p)[op == "~", .(model = names(fits)[i], lhs, op, rhs, estimate = est, se, z, pvalue = pvalue,
    ci_low = ci.lower, ci_high = ci.upper, standardized_estimate = std.all)]
}))
params[, fdr := p.adjust(pvalue, "BH"), by = model]

predict_fate <- function(model_name, train, test) {
  form <- if (grepl("Linear", model_name)) fate ~ C else fate ~ A + B + C
  fit <- stats::lm(form, data = train)
  as.numeric(stats::predict(fit, newdata = test))
}
set.seed(20260805)
folds <- rsample::vfold_cv(as.data.frame(dat), v = min(5L, nrow(dat)), repeats = 10L, strata = dataset)
cv <- rbindlist(lapply(seq_len(nrow(folds)), function(i) {
  train <- as.data.table(rsample::analysis(folds$splits[[i]])); test <- as.data.table(rsample::assessment(folds$splits[[i]]))
  rbindlist(lapply(names(models), function(m) {
    pred <- predict_fate(m, train, test)
    data.table(model = m, resample_id = paste(folds$id[i], folds$id2[i], sep = "_"),
      rmse = sqrt(mean((test$fate - pred)^2, na.rm = TRUE)), mae = mean(abs(test$fate - pred), na.rm = TRUE), n_test = nrow(test))
  }))
}))
figure6_fwrite(cv, file.path(FIGURE6_METADATA_DIR, "figure6f_cross_validation.tsv"), compress = FALSE)
cv_summary <- cv[, .(cv_rmse = mean(rmse), cv_rmse_sd = sd(rmse), cv_mae = mean(mae)), by = model]

edge_defs <- list(
  "Model 1: Linear cascade" = list(c("B", "A"), c("C", "B"), c("fate", "C")),
  "Model 2: Parallel convergence" = list(c("fate", "A"), c("fate", "B"), c("fate", "C")),
  "Model 3: Partially ordered" = list(c("B", "A"), c("C", "A"), c("C", "B"), c("fate", "A"), c("fate", "B"), c("fate", "C"))
)
set.seed(20260815)
boot_rows <- list(); k <- 1L
for (m in names(edge_defs)) for (edge in edge_defs[[m]]) {
  lhs <- edge[1]; rhs <- edge[2]
  formula_txt <- if (m == "Model 3: Partially ordered" && lhs %in% c("C", "fate")) {
    if (lhs == "C") "C ~ A + B" else "fate ~ A + B + C"
  } else if (m == "Model 2: Parallel convergence") "fate ~ A + B + C" else paste(lhs, "~", rhs)
  original <- coef(lm(as.formula(formula_txt), data = dat))[rhs]
  boot_coef <- replicate(1000L, {
    ix <- unlist(lapply(split(seq_len(nrow(dat)), dat$dataset), function(z) sample(z, length(z), replace = TRUE)), use.names = FALSE)
    unname(coef(lm(as.formula(formula_txt), data = dat[ix]))[rhs])
  })
  lodo <- vapply(unique(dat$dataset), function(ds) {
    tr <- dat[dataset != ds]
    if (nrow(tr) < 6) return(NA_real_)
    unname(coef(lm(as.formula(formula_txt), data = tr))[rhs])
  }, numeric(1))
  boot_rows[[k]] <- data.table(model = m, lhs, rhs, edge = paste(rhs, "→", lhs), original_estimate = original,
    bootstrap_sign_stability = mean(sign(boot_coef) == sign(original), na.rm = TRUE),
    bootstrap_ci_low = quantile(boot_coef, .025, na.rm = TRUE), bootstrap_ci_high = quantile(boot_coef, .975, na.rm = TRUE),
    lodo_sign_stability = mean(sign(lodo) == sign(original), na.rm = TRUE), n_bootstrap = 1000L, n_lodo = sum(is.finite(lodo)))
  k <- k + 1L
}
stability <- rbindlist(boot_rows)
figure6_fwrite(stability, file.path(FIGURE6_METADATA_DIR, "figure6f_bootstrap_edge_stability.tsv"), compress = FALSE)

fit_summary <- merge(fit_summary, cv_summary, by = "model", all.x = TRUE)
st <- stability[, .(mean_edge_stability = mean(bootstrap_sign_stability), mean_lodo_stability = mean(lodo_sign_stability)), by = model]
fit_summary <- merge(fit_summary, st, by = "model", all.x = TRUE)
fit_summary[, `:=`(rank_bic = frank(bic, ties.method = "min"), rank_aic = frank(aic, ties.method = "min"),
  rank_cv = frank(cv_rmse, ties.method = "min"), rank_stability = frank(-mean_edge_stability, ties.method = "min"))]
fit_summary[, selection_rank_sum := rank_bic + rank_aic + rank_cv + rank_stability]
fit_summary[, selected := selection_rank_sum == min(selection_rank_sum)]
figure6_fwrite(fit_summary, file.path(FIGURE6_METADATA_DIR, "figure6f_model_fit_summary.tsv"), compress = FALSE)
figure6_fwrite(params, file.path(FIGURE6_METADATA_DIR, "figure6f_model_parameters.tsv"), compress = FALSE)
selected <- fit_summary[selected == TRUE, model]
figure6_write_json(list(
  panel = "Figure 6F", analysis_unit = "sample-level pseudobulk; three trajectory methods averaged within dataset/sample",
  n_samples = nrow(dat), n_datasets = uniqueN(dat$dataset), models = as.list(models),
  selection_rule = "minimum sum of ranks across AIC, BIC, repeated 10x5-fold fate-prediction RMSE and mean bootstrap edge sign stability",
  selected_model = as.list(selected), bootstrap = 1000, cv = "10 repeats of 5-fold CV, dataset-stratified",
  identifiability = "Model 3 is an acyclic partially ordered SEM; no simultaneous reciprocal causal paths were fitted.",
  limitation = "Models 2 and 3 share the same fate regression and can tie in fate-only CV; saturated models are not favored by CFI/RMSEA alone."
), file.path(FIGURE6_METADATA_DIR, "figure6f_model_report.json"))
