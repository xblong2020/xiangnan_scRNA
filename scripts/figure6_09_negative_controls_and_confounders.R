#!/usr/bin/env Rscript

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figure6_common.R"))
suppressPackageStartupMessages(library(broom))

effects <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6_perturbation_response_effects.tsv.gz"))[availability == "Available"]
sample <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6_sample_level_effects.tsv.gz"))
candidates <- c("HNF4A", "PPARA", "EGR1", "CEBPB", "AP1_AGGREGATE", "SOX4")
controls <- FIGURE6_AXIS_TFS$control

neg <- effects[tf %in% candidates, {
  null <- effects[tf %in% controls & output == .BY$output, abs(effect_estimate)]
  obs <- abs(effect_estimate)
  list(
    effect_estimate = effect_estimate, absolute_effect = obs,
    control_median_absolute_effect = median(null, na.rm = TRUE),
    control_max_absolute_effect = max(null, na.rm = TRUE),
    empirical_control_pvalue = (1 + sum(null >= obs, na.rm = TRUE)) / (1 + sum(is.finite(null))),
    n_control_tfs = sum(is.finite(null)), candidate_exceeds_control_median = obs > median(null, na.rm = TRUE),
    candidate_exceeds_all_controls = obs > max(null, na.rm = TRUE)
  )
}, by = .(tf, perturbation, axis, output)]
neg[, empirical_fdr := p.adjust(empirical_control_pvalue, "BH")]
neg[, specificity_status := fifelse(n_control_tfs < 10, "specificity_risk_limited_control_pool",
  fifelse(empirical_fdr < .05, "specific", "not_specific"))]
figure6_fwrite(neg, file.path(FIGURE6_METADATA_DIR, "figure6_negative_control_results.tsv.gz"))

adjust_one <- function(z) {
  needed <- c("standardized_effect", "baseline_stress", "baseline_proliferation", "dataset")
  z <- z[complete.cases(z[, ..needed])]
  if (nrow(z) < 10 || uniqueN(z$dataset) < 2) return(data.table(
    adjusted_effect = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_, pvalue = NA_real_, n_samples = nrow(z),
    n_datasets = uniqueN(z$dataset), cnv_adjusted_effect = NA_real_, cnv_adjusted_pvalue = NA_real_, n_samples_cnv = 0L,
    status = "not_estimable"))
  z[, c("baseline_stress", "baseline_proliferation") := lapply(.SD, figure6_safe_scale),
    .SDcols = c("baseline_stress", "baseline_proliferation")]
  fit <- tryCatch(lm(standardized_effect ~ baseline_stress + baseline_proliferation + factor(dataset), data = z), error = function(e) NULL)
  if (is.null(fit)) return(data.table(adjusted_effect=NA_real_,se=NA_real_,ci_low=NA_real_,ci_high=NA_real_,pvalue=NA_real_,n_samples=nrow(z),n_datasets=uniqueN(z$dataset),status="fit_failed"))
  # Equal-weight dataset marginal mean at centered covariates.
  newdata <- data.frame(
    baseline_stress = 0, baseline_proliferation = 0,
    dataset = factor(sort(unique(z$dataset)), levels = levels(factor(z$dataset)))
  )
  mm <- model.matrix(delete.response(terms(fit)), newdata)
  contrast <- colMeans(mm)
  beta <- coef(fit); vc <- vcov(fit)
  estimate <- sum(contrast * beta)
  se <- sqrt(as.numeric(t(contrast) %*% vc %*% contrast))
  z_value <- estimate / se
  z_cnv <- z[is.finite(baseline_cnv_proxy)]
  cnv_estimate <- cnv_pvalue <- NA_real_
  if (nrow(z_cnv) >= 10 && uniqueN(z_cnv$dataset) >= 2) {
    z_cnv[, baseline_cnv_proxy := figure6_safe_scale(baseline_cnv_proxy)]
    fit_cnv <- tryCatch(lm(standardized_effect ~ baseline_stress + baseline_proliferation + baseline_cnv_proxy + factor(dataset), data=z_cnv), error=function(e) NULL)
    if (!is.null(fit_cnv)) {
      nd_cnv <- data.frame(baseline_stress=0, baseline_proliferation=0, baseline_cnv_proxy=0,
        dataset=factor(sort(unique(z_cnv$dataset)), levels=levels(factor(z_cnv$dataset))))
      mm_cnv <- model.matrix(delete.response(terms(fit_cnv)), nd_cnv); contrast_cnv <- colMeans(mm_cnv)
      cnv_estimate <- sum(contrast_cnv * coef(fit_cnv))
      cnv_se <- sqrt(as.numeric(t(contrast_cnv) %*% vcov(fit_cnv) %*% contrast_cnv))
      cnv_pvalue <- 2*pnorm(-abs(cnv_estimate/cnv_se))
    }
  }
  data.table(adjusted_effect=estimate, se=se, ci_low=estimate-1.96*se, ci_high=estimate+1.96*se,
    pvalue=2*pnorm(-abs(z_value)), cnv_adjusted_effect=cnv_estimate, cnv_adjusted_pvalue=cnv_pvalue,
    n_samples=nrow(z), n_samples_cnv=nrow(z_cnv), n_datasets=uniqueN(z$dataset), status="estimable")
}
conf <- sample[tf %in% candidates, adjust_one(.SD), by = .(tf, perturbation, axis, output)]
conf[, fdr := p.adjust(pvalue, "BH")]
conf[, `:=`(
  generic_stress_adjustment = "baseline frozen Stressed_Injured module score",
  proliferation_adjustment = "baseline frozen Proliferation module score",
  cnv_adjustment = "extended sensitivity model using baseline cnv_proxy_z where >=10 samples",
  dataset_adjustment = "fixed effect",
  hypoxia_adjustment = "not testable: no frozen hypoxia score",
  s_g2m_adjustment = "not testable: separate S/G2M scores unavailable"
)]
figure6_fwrite(conf, file.path(FIGURE6_METADATA_DIR, "figure6_confounder_adjustment_summary.tsv"), compress = FALSE)
figure6_write_json(list(
  analysis = "Figure 6 negative-control and confounder sensitivity",
  negative_controls = as.list(controls), n_controls = length(controls), requested_minimum = 10,
  specificity_risk = "Only five eligible frozen TF controls were available; empirical p-values cannot reach 0.05.",
  adjusted_covariates = c("frozen Stressed_Injured module", "frozen Proliferation module", "dataset fixed effect", "extended cnv_proxy_z sensitivity"),
  unavailable_covariates = c("hypoxia score", "separate S score", "separate G2M score"),
  interpretation = "Sensitivity analyses assess robustness of sample-level predicted effects and do not establish causal specificity."
), file.path(FIGURE6_METADATA_DIR, "figure6_sensitivity_report.json"))
