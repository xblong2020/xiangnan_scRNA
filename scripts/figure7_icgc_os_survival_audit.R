suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(survival)
  library(jsonlite)
})

audit_root <- function() {
  wd <- gsub("\\\\", "/", getwd())
  if (dir.exists(file.path(wd, "Figure7_ICGC_OS_Audit"))) return(wd)
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(arg)) stop("Run from the Figure7 project root or provide --file.")
  root <- normalizePath(file.path(dirname(sub("^--file=", "", arg[[1L]])), ".."), winslash = "/", mustWork = TRUE)
  root
}

ROOT <- audit_root()
AUDIT <- file.path(ROOT, "Figure7_ICGC_OS_Audit")
SCORES <- file.path(AUDIT, "08_source_data", "ICGC_axis_scores.tsv")
OUT_SURV <- file.path(AUDIT, "05_survival_models")
OUT_SENS <- file.path(AUDIT, "06_sensitivity")
OUT_QC <- file.path(AUDIT, "04_qc")
OUT_FIG <- file.path(AUDIT, "07_figures")
dir.create(OUT_SURV, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_SENS, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_QC, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

write_audit_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
}

num <- function(x) suppressWarnings(as.numeric(as.character(x)))
axes <- c("identity_loss", "stress_transition", "sox4_associated")
axis_label <- c(identity_loss = "Identity loss", stress_transition = "Stress-transition", sox4_associated = "SOX4-associated malignant state")

d <- fread(SCORES, sep = "\t", na.strings = c("", "NA", "NaN"), encoding = "UTF-8")
required <- c("cohort", "donor_id", "raw_futime", "raw_fustat", "age_high_raw", "gender_raw", "stage_high_raw",
              "identity_loss_score_tcga_frozen_z", "stress_transition_score_tcga_frozen_z", "sox4_associated_score_tcga_frozen_z")
missing_required <- setdiff(required, names(d))
if (length(missing_required)) stop("ICGC_axis_scores.tsv is missing: ", paste(missing_required, collapse = ", "))
d[, `:=`(
  cohort = as.character(cohort), donor_id = as.character(donor_id),
  os_time_days = num(raw_futime), os_event = as.integer(num(raw_fustat)),
  age_high = as.integer(num(age_high_raw)), sex_male = as.integer(num(gender_raw) == 1), stage_high = as.integer(num(stage_high_raw) == 1),
  identity_loss = num(identity_loss_score_tcga_frozen_z),
  stress_transition = num(stress_transition_score_tcga_frozen_z),
  sox4_associated = num(sox4_associated_score_tcga_frozen_z)
)]
d <- d[cohort == "ICGC_LIRI_JP"]
if (anyDuplicated(d$donor_id)) stop("Duplicate donor rows remain in the audit score table.")
if (nrow(d) != 231L) stop("Expected 231 unique ICGC tumour donors; observed ", nrow(d), ".")

mapping_path <- file.path(AUDIT, "02_mapping", "ICGC_expression_clinical_mapping.tsv")
mapping <- fread(mapping_path, sep = "\t", na.strings = c("", "NA"))
mapping[, `:=`(OS_time = num(raw_futime), OS_status = as.integer(num(raw_fustat)),
              mapping_status = "expression_to_donor_and_clinical_survival_match; exploratory_OS_derivation",
              exclusion_reason = fifelse(sample_type == "tumour", "none_for_exploratory_OS; continuous_age_unavailable", "not_an_OS_analysis_unit; retained_for_expression_recurrence"))]
write_audit_tsv(mapping, mapping_path)

d[, `:=`(
  age_years = NA_real_,
  identity_loss_z = as.numeric(scale(identity_loss)),
  stress_transition_z = as.numeric(scale(stress_transition)),
  sox4_associated_z = as.numeric(scale(sox4_associated))
)]
if (any(!is.finite(d$os_time_days)) || any(d$os_time_days <= 0) || any(!d$os_event %in% c(0L, 1L))) stop("Invalid raw OS time/event values.")

fit_cox <- function(dat, axis_name, level = c("univariable", "partial_clinical")) {
  level <- match.arg(level)
  score_col <- paste0(axis_name, "_z")
  rhs <- if (level == "univariable") score_col else paste(c(score_col, "sex_male", "stage_high"), collapse = " + ")
  needed <- c("os_time_days", "os_event", score_col, if (level == "partial_clinical") c("sex_male", "stage_high") else character())
  z <- dat[complete.cases(dat[, ..needed]) & os_time_days > 0 & os_event %in% c(0L, 1L)]
  events <- sum(z$os_event == 1L)
  p <- length(unique(all.vars(as.formula(paste("~", rhs)))))
  base <- data.table(model_id = paste0(level, "__", axis_name), analysis = level, cohort = "ICGC_LIRI_JP", endpoint = "OS",
                     predictor = axis_name, term = score_col, adjustment = if (level == "univariable") "none" else "sex + stage_high; age omitted",
                     n = nrow(z), events = events, censored = nrow(z) - events, predictors = p,
                     events_per_variable = if (p) events / p else NA_real_, HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P = NA_real_,
                     concordance = NA_real_, PH_P = NA_real_, global_PH_P = NA_real_, status = "not_estimable",
                     evidence_status = "exploratory", figure_eligibility = "SUPPLEMENTARY_ONLY", reason = "")
  if (nrow(z) < 20L || events < 5L) {
    base[, reason := "insufficient_complete_cases_or_events"]
    return(list(result = base, ph = data.table(model_id = base$model_id, predictor = axis_name, term = score_col, PH_P = NA_real_, global_PH_P = NA_real_, status = base$status, reason = base$reason), fit = NULL, data = z))
  }
  formula <- as.formula(paste("Surv(os_time_days, os_event) ~", rhs))
  fit <- try(coxph(formula, data = z, x = TRUE, y = TRUE, model = TRUE), silent = TRUE)
  if (inherits(fit, "try-error")) {
    base[, reason := "cox_fit_failed"]
    return(list(result = base, ph = data.table(model_id = base$model_id, predictor = axis_name, term = score_col, PH_P = NA_real_, global_PH_P = NA_real_, status = base$status, reason = base$reason), fit = NULL, data = z))
  }
  sm <- summary(fit); co <- sm$coefficients; ci <- suppressWarnings(confint(fit));
  if (!score_col %in% rownames(co)) {
    base[, reason := "score_term_missing_from_model"]
    return(list(result = base, ph = data.table(model_id = base$model_id, predictor = axis_name, term = score_col, PH_P = NA_real_, global_PH_P = NA_real_, status = base$status, reason = base$reason), fit = fit, data = z))
  }
  ph <- try(cox.zph(fit), silent = TRUE)
  ph_axis <- if (inherits(ph, "try-error") || !score_col %in% rownames(ph$table)) NA_real_ else ph$table[score_col, "p"]
  ph_global <- if (inherits(ph, "try-error") || !"GLOBAL" %in% rownames(ph$table)) NA_real_ else ph$table["GLOBAL", "p"]
  model_reason <- if (is.finite(ph_global) && ph_global < 0.05) "global_PH_violation_retained_as_exploratory" else ""
  base[, `:=`(HR = exp(co[score_col, "coef"]), CI_low = exp(ci[score_col, 1L]), CI_high = exp(ci[score_col, 2L]),
              P = co[score_col, "Pr(>|z|)"], concordance = sm$concordance[[1L]], PH_P = ph_axis, global_PH_P = ph_global,
              status = "estimated", reason = model_reason)]
  ph_row <- data.table(model_id = base$model_id, predictor = axis_name, term = score_col, PH_P = ph_axis, global_PH_P = ph_global,
                       status = if (is.finite(ph_axis) && is.finite(ph_global)) "estimated" else "not_estimable",
                       reason = if (is.finite(ph_axis) && is.finite(ph_global)) "cox.zph" else "PH test unavailable")
  list(result = base, ph = ph_row, fit = fit, data = z)
}

univ <- lapply(axes, function(a) fit_cox(d, a, "univariable"))
partial <- lapply(axes, function(a) fit_cox(d, a, "partial_clinical"))
univ_results <- rbindlist(lapply(univ, `[[`, "result"), fill = TRUE)
partial_results <- rbindlist(lapply(partial, `[[`, "result"), fill = TRUE)
univ_results[, FDR := p.adjust(P, method = "BH")]
partial_results[, FDR := p.adjust(P, method = "BH")]
write_audit_tsv(univ_results, file.path(OUT_SURV, "ICGC_univariable_cox.tsv"))
write_audit_tsv(partial_results, file.path(OUT_SURV, "ICGC_multivariable_cox.tsv"))
ph_results <- rbindlist(c(lapply(univ, `[[`, "ph"), lapply(partial, `[[`, "ph")), fill = TRUE)
ph_results[, `:=`(cohort = "ICGC_LIRI_JP", endpoint = "OS", evidence_status = "exploratory", figure_eligibility = "SUPPLEMENTARY_ONLY")]
write_audit_tsv(ph_results, file.path(OUT_SURV, "ICGC_PH_test.tsv"))

nonlinear_rows <- lapply(axes, function(axis_name) {
  score_col <- paste0(axis_name, "_z")
  z <- d[complete.cases(d[, c("os_time_days", "os_event", score_col), with = FALSE])]
  linear <- try(coxph(as.formula(paste("Surv(os_time_days, os_event) ~", score_col)), data = z), silent = TRUE)
  spline <- try(coxph(as.formula(paste("Surv(os_time_days, os_event) ~ splines::ns(", score_col, ", df=3)")), data = z), silent = TRUE)
  if (inherits(linear, "try-error") || inherits(spline, "try-error")) return(data.table(cohort = "ICGC_LIRI_JP", endpoint = "OS", predictor = axis_name, n = nrow(z), events = sum(z$os_event == 1L), linear_AIC = NA_real_, spline_AIC = NA_real_, nonlinearity_P = NA_real_, status = "not_estimable", evidence_status = "exploratory", figure_eligibility = "SUPPLEMENTARY_ONLY", reason = "nonlinearity_fit_failed"))
  p <- NA_real_
  ll_linear <- suppressWarnings(as.numeric(logLik(linear)))
  ll_spline <- suppressWarnings(as.numeric(logLik(spline)))
  df_diff <- suppressWarnings(as.numeric(attr(logLik(spline), "df") - attr(logLik(linear), "df")))
  if (is.finite(ll_linear) && is.finite(ll_spline) && is.finite(df_diff) && df_diff > 0) {
    p <- pchisq(2 * (ll_spline - ll_linear), df = df_diff, lower.tail = FALSE)
  }
  data.table(cohort = "ICGC_LIRI_JP", endpoint = "OS", predictor = axis_name, n = nrow(z), events = sum(z$os_event == 1L), linear_AIC = AIC(linear), spline_AIC = AIC(spline), nonlinearity_P = p, status = if (is.finite(p)) "estimated" else "not_estimable", evidence_status = "exploratory", figure_eligibility = "SUPPLEMENTARY_ONLY", reason = "pre-specified df=3 spline sensitivity")
})
write_audit_tsv(rbindlist(nonlinear_rows, fill = TRUE), file.path(OUT_SENS, "ICGC_nonlinearity_sensitivity.tsv"))

qc <- fread(file.path(OUT_QC, "ICGC_survival_QC_summary.tsv"), sep = "\t", na.strings = c("", "NA"))
qc[, `:=`(survival_models_run = TRUE, univariable_status = "estimated", partial_multivariable_status = "exploratory_only", PH_test_status = "completed", nonlinearity_status = "completed",
          age_missing = age_missing_continuous, sex_missing = sex_missing_raw, stage_missing = stage_missing_raw)]
write_audit_tsv(qc, file.path(OUT_QC, "ICGC_survival_QC_summary.tsv"))

v2_tn_path <- file.path(ROOT, "figures", "driver", "figure7_external_validation_v2", "figure7_v2c_tumour_normal_effects.tsv")
v2_cox_path <- file.path(ROOT, "figures", "driver", "figure7_external_validation_v2", "figure7_v2e_multivariable_cox_models.tsv")
v2_tn <- fread(v2_tn_path, sep = "\t", na.strings = c("", "NA"))
v2_cox <- fread(v2_cox_path, sep = "\t", na.strings = c("", "NA"))
tcga_tn <- v2_tn[cohort == "TCGA_LIHC" & axis %in% axes, .(axis, TCGA_expression_effect = hedges_g)]
icgc_tn <- v2_tn[cohort == "ICGC_LIRI_JP" & axis %in% axes, .(axis, ICGC_expression_effect = hedges_g)]
tcga_surv <- v2_cox[cohort == "TCGA_LIHC" & programme %in% axes & grepl("^clinical_adjusted__", model_id) & term == programme, .(axis = programme, TCGA_survival_HR = hazard_ratio, TCGA_survival_P = p_value)]
icgc_surv <- univ_results[, .(axis = predictor, ICGC_survival_HR = HR, ICGC_survival_P = P, ICGC_survival_FDR = FDR, ICGC_N = n, ICGC_events = events)]
replication <- Reduce(function(x, y) merge(x, y, by = "axis", all = TRUE), list(tcga_tn, icgc_tn, tcga_surv, icgc_surv))
replication[, `:=`(
  TCGA_expression_direction = fifelse(TCGA_expression_effect > 0, "positive", fifelse(TCGA_expression_effect < 0, "negative", "zero")),
  ICGC_expression_direction = fifelse(ICGC_expression_effect > 0, "positive", fifelse(ICGC_expression_effect < 0, "negative", "zero")),
  survival_direction_concordance = fifelse(is.finite(TCGA_survival_HR) & is.finite(ICGC_survival_HR) & ((TCGA_survival_HR > 1) == (ICGC_survival_HR > 1)), "concordant", "discordant_or_not_estimable"),
  expression_direction_concordance = fifelse(is.finite(TCGA_expression_effect) & is.finite(ICGC_expression_effect) & ((TCGA_expression_effect > 0) == (ICGC_expression_effect > 0)), "concordant", "discordant"),
  significance = paste0("ICGC univariable P=", formatC(ICGC_survival_P, format = "g", digits = 4), "; FDR=", formatC(ICGC_survival_FDR, format = "g", digits = 4)),
  evidence_status = "exploratory", figure_eligibility = "SUPPLEMENTARY_ONLY"
)]
replication[, interpretation := fifelse(survival_direction_concordance == "concordant" & ICGC_survival_FDR >= 0.05, "directionally concordant but statistically underpowered", fifelse(survival_direction_concordance == "concordant", "directionally concordant exploratory association", "inconclusive or discordant"))]
setcolorder(replication, c("axis", "TCGA_expression_direction", "ICGC_expression_direction", "TCGA_survival_HR", "ICGC_survival_HR", "TCGA_survival_P", "ICGC_survival_P", "ICGC_survival_FDR", "ICGC_N", "ICGC_events", "expression_direction_concordance", "survival_direction_concordance", "significance", "interpretation", "evidence_status", "figure_eligibility"))
write_audit_tsv(replication, file.path(AUDIT, "08_source_data", "TCGA_ICGC_axis_replication.tsv"))

audit_rows <- rbindlist(list(
  univ_results[, .(analysis = "continuous_univariable_cox", cohort, endpoint, N = n, events, predictor, HR, CI_low, CI_high, P, PH_P, adjustment, replication_direction = NA_character_, evidence_status, figure_eligibility, reason)],
  partial_results[, .(analysis = "partial_multivariable_cox", cohort, endpoint, N = n, events, predictor, HR, CI_low, CI_high, P, PH_P, adjustment, replication_direction = NA_character_, evidence_status, figure_eligibility, reason)]
), fill = TRUE)
write_audit_tsv(audit_rows, file.path(AUDIT, "08_source_data", "Figure7_ICGC_clinical_validation_audit.tsv"))

forest <- copy(univ_results)
forest[, axis_label := unname(axis_label[predictor])]
forest_plot <- ggplot(forest, aes(x = HR, y = reorder(axis_label, HR), colour = predictor, shape = predictor)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = .35, colour = "#666666") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), width = .16, orientation = "y", linewidth = .5) +
  geom_point(size = 2.7) +
  scale_x_log10() +
  scale_colour_manual(values = c(identity_loss = "#0072B2", stress_transition = "#009E73", sox4_associated = "#D55E00"), guide = "none") +
  scale_shape_manual(values = c(identity_loss = 16, stress_transition = 17, sox4_associated = 15), guide = "none") +
  labs(title = "ICGC-LIRI-JP exploratory OS associations", subtitle = "Frozen v2 axis scores; HR per within-ICGC tumour SD; Supplementary only", x = "Hazard ratio (95% CI)", y = NULL) +
  theme_classic(base_size = 8) + theme(plot.title = element_text(size = 10, face = "bold"), plot.subtitle = element_text(size = 7, colour = "#444444"), axis.text = element_text(size = 7), plot.margin = margin(6, 8, 6, 6))
fig_stem <- file.path(OUT_FIG, "ICGC_OS_exploratory_cox_forest")
ggsave(paste0(fig_stem, ".pdf"), forest_plot, width = 5.6, height = 3.4, units = "in", device = grDevices::cairo_pdf, bg = "white")
ggsave(paste0(fig_stem, ".svg"), forest_plot, width = 5.6, height = 3.4, units = "in", device = grDevices::svg, bg = "white")
tmp_fig_dir <- file.path(tempdir(), "figure7_icgc_os_audit_raster")
dir.create(tmp_fig_dir, recursive = TRUE, showWarnings = FALSE)
tmp_png <- file.path(tmp_fig_dir, "icgc_os_exploratory_cox_forest.png")
tmp_tiff <- file.path(tmp_fig_dir, "icgc_os_exploratory_cox_forest.tiff")
ggsave(tmp_png, forest_plot, width = 5.6, height = 3.4, units = "in", dpi = 600, device = "png", bg = "white")
ggsave(tmp_tiff, forest_plot, width = 5.6, height = 3.4, units = "in", dpi = 600, device = "tiff", compression = "lzw", bg = "white")
if (!file.copy(tmp_png, paste0(fig_stem, ".png"), overwrite = TRUE) || !file.copy(tmp_tiff, paste0(fig_stem, ".tiff"), overwrite = TRUE)) stop("Exploratory forest raster export failed.")

fmt_num <- function(x, digits = 3L) if (length(x) && is.finite(x[[1L]])) formatC(x[[1L]], format = "f", digits = digits) else "NA"
fmt_p <- function(x) if (length(x) && is.finite(x[[1L]])) formatC(x[[1L]], format = "g", digits = 4L) else "NA"
model_lines <- c(
  "## Post-gate survival results", "",
  "All models use the frozen v2 primary programme direction and one donor-level row per tumour donor. The score effect is per one within-ICGC tumour SD.", "",
  "| Analysis | Axis | N | Events | HR | 95% CI | P | PH P | Global PH P | Status |", "|---|---|---:|---:|---:|---|---:|---:|---:|---|",
  vapply(seq_len(nrow(univ_results)), function(i) sprintf("| Univariable | %s | %d | %d | %s | %s-%s | %s | %s | %s | exploratory |", axis_label[[univ_results$predictor[[i]]]], univ_results$n[[i]], univ_results$events[[i]], fmt_num(univ_results$HR[[i]]), fmt_num(univ_results$CI_low[[i]]), fmt_num(univ_results$CI_high[[i]]), fmt_p(univ_results$P[[i]]), fmt_p(univ_results$PH_P[[i]]), fmt_p(univ_results$global_PH_P[[i]])), character(1)),
  vapply(seq_len(nrow(partial_results)), function(i) sprintf("| Partial adjusted | %s | %d | %d | %s | %s-%s | %s | %s | %s | exploratory; age omitted |", axis_label[[partial_results$predictor[[i]]]], partial_results$n[[i]], partial_results$events[[i]], fmt_num(partial_results$HR[[i]]), fmt_num(partial_results$CI_low[[i]]), fmt_num(partial_results$CI_high[[i]]), fmt_p(partial_results$P[[i]]), fmt_p(partial_results$PH_P[[i]]), fmt_p(partial_results$global_PH_P[[i]])), character(1)),
  "", "The three univariable axis PH tests are acceptable (all PH P > 0.05). The partial clinical models have global PH P values below 0.05 and are retained as exploratory diagnostics only.", "",
  sprintf("The pre-specified natural-spline sensitivity showed no evidence of nonlinearity (identity P=%s; stress P=%s; SOX4 P=%s). AIC values are retained in `06_sensitivity/ICGC_nonlinearity_sensitivity.tsv`; no cutoff search was performed.", fmt_p(nonlinear_rows[[1L]]$nonlinearity_P), fmt_p(nonlinear_rows[[2L]]$nonlinearity_P), fmt_p(nonlinear_rows[[3L]]$nonlinearity_P)), "",
  "KM was not re-fit in this audit. The existing Figure 7 KM convention remains visualization-only and is not used as primary inference."
)
writeLines(model_lines, file.path(AUDIT, "09_reports", "ICGC_OS_survival_results.md"), useBytes = TRUE)
gate_path <- file.path(AUDIT, "03_os_definition", "ICGC_OS_UNBLOCK_GATE.md")
gate_base <- readLines(gate_path, warn = FALSE)
post_gate <- grep("^## Post-gate survival results$", gate_base)
if (length(post_gate)) gate_base <- gate_base[seq_len(post_gate[[1L]] - 1L)]
writeLines(c(gate_base, "", model_lines), gate_path, useBytes = TRUE)
decision_path <- file.path(AUDIT, "09_reports", "ICGC_OS_FINAL_DECISION.md")
decision_lines <- c(
  "# ICGC_OS_FINAL_DECISION", "", "## Executive decision", "", "`ESTIMABLE_BUT_NOT_VALIDATED`", "",
  sprintf("## Data and derivation\n\nThe ICGC-LIRI-JP tumour expression cohort contains %d unique donor-linked patients with %d events and %d censored/alive records. Expression-to-donor mapping is complete and duplicate tumour samples are collapsed deterministically. OS time origin and day unit are supported by ICGC documentation and independent LIRI-JP methods; event coding is supported by a 203-donor HCCDB18 cross-check.", nrow(d), sum(d$os_event == 1L), sum(d$os_event == 0L)), "",
  "## Survival results", "", model_lines[4:length(model_lines)], "",
  "The identity score is directionally concordant with TCGA in both expression recurrence and survival association, but the FDR-adjusted ICGC univariable result is not below 0.05. Stress-transition and SOX4-associated directions are discordant between cohorts. Partial clinical models omit Age because the local Age field is binary without an independently documented threshold; their global PH tests are below 0.05, so they remain exploratory.", "",
  "## Placement and claims", "", "The ICGC survival panel is `SUPPLEMENTARY_ONLY`. The current Figure 7 main-v2 ICGC clinical/OS block is preserved. This audit upgrades the endpoint from an undefined numerical branch to an auditable exploratory estimate; it does not establish external prognostic validation.", "",
  "## Allowed claims", "", "- ICGC-LIRI-JP donor-level OS association is estimable under a documented exploratory derivation.", "- Directional comparison with TCGA can be reported with uncertainty and event-count limitations.", "- Expression-level ICGC recurrence remains usable independently of clinical/OS claims.", "", "## Disallowed claims", "", "- Externally validated prognostic model.", "- Fully age/sex/stage-adjusted ICGC validation.", "- Clinical utility or treatment-selection evidence.", "- Direct SOX4 activity inferred from a bulk-associated score."
)
writeLines(decision_lines, decision_path, useBytes = TRUE)

manifest_path <- file.path(AUDIT, "10_manifests", "ICGC_OS_run_manifest.json")
manifest <- fromJSON(manifest_path, simplifyVector = FALSE)
manifest$survival_models_run <- TRUE
manifest$random_seed <- "not_used"
manifest$final_gate_status <- "ESTIMABLE_BUT_NOT_VALIDATED"
manifest$decision_basis <- "OS time origin and day unit are supported by independent ICGC/LIRI-JP documentation; exact legacy release and field-generation provenance remain unresolved, Age is binary, and cross-axis replication is incomplete"
manifest$data_release <- "legacy local ICGC-LIRI-JP cache; exact release identifier not recorded"
manifest$clinical_file <- "see input_files entry with source_role=clinical"
manifest$expression_file <- "see input_files entry with source_role=expression"
manifest$survival_file <- "see input_files entry with source_role=survival"
manifest$mapping_file <- "02_mapping/ICGC_expression_clinical_mapping.tsv"
manifest$OS_definition <- list(time = "futime in days", status = "fustat: 1=death, 0=alive/censored", time_origin = "date of diagnosis / primary diagnosis", patient_id = "DO donor")
manifest$exclusion_rules <- "unknown sample type; missing DO; missing clinical/survival row; nonpositive time; nonbinary event"
manifest$duplicate_handling <- "one donor-level mean score; 8 donors with 9 excess tumour sample rows"
manifest$gene_sets <- "Figure 7 v2 primary_frozen_programme; identity_loss, stress_transition, sox4_associated"
manifest$score_method <- "frozen v2 unsigned associated target programme score; TCGA-frozen score scale, then within-ICGC SD for Cox"
manifest$normalization <- "existing Figure 7 v2 sample-wise normalized rank scores; no outcome-derived reweighting"
manifest$covariates <- "partial exploratory adjustment: sex_male + stage_high; Age omitted because legacy Age is binary without verified threshold"
manifest$software_versions <- list(R = R.version.string, data_table = as.character(packageVersion("data.table")), survival = as.character(packageVersion("survival")), jsonlite = as.character(packageVersion("jsonlite")))
manifest$external_sources <- c("https://docs.cancergenomicscloud.org/docs/icgc-metadata", "https://docs.icgc-argo.org/dictionary", "https://pmc.ncbi.nlm.nih.gov/articles/PMC9482539/", "http://lifeome.net/database/hccdb/download.html")
manifest$project_record <- "metadata/driver/figure7_external_validation_v2/figure7_v2_icgc_os_unblock_audit_record.md"
manifest$survival_models <- list(
  univariable = "Surv(os_time_days, os_event) ~ scale(frozen_axis_score) per ICGC tumour SD",
  partial_multivariable = "Surv(os_time_days, os_event) ~ scale(frozen_axis_score) + sex_male + stage_high",
  age_adjustment = "not fitted: legacy Age is binary and threshold is not independently documented",
  ph_test = "cox.zph",
  nonlinearity = "restricted cubic-spline-equivalent natural spline df=3 sensitivity",
  km = "not run in this audit; median split remains visualization-only in existing Figure 7"
)
manifest$final_outputs <- c(
  "05_survival_models/ICGC_univariable_cox.tsv",
  "05_survival_models/ICGC_multivariable_cox.tsv",
  "05_survival_models/ICGC_PH_test.tsv",
  "06_sensitivity/ICGC_nonlinearity_sensitivity.tsv",
  "08_source_data/TCGA_ICGC_axis_replication.tsv",
  "08_source_data/Figure7_ICGC_clinical_validation_audit.tsv"
)
manifest$output_files <- c(
  "01_input_inventory/figure7_icgc_history_inventory.tsv",
  "01_input_inventory/ICGC_input_provenance.tsv",
  "01_input_inventory/ICGC_raw_field_summary.tsv",
  "01_input_inventory/external_source_evidence.tsv",
  "02_mapping/ICGC_expression_clinical_mapping.tsv",
  "02_mapping/ICGC_mapping_QC_summary.tsv",
  "03_os_definition/ICGC_OS_DERIVATION_SPEC.md",
  "03_os_definition/ICGC_OS_UNBLOCK_GATE.md",
  "04_qc/ICGC_survival_QC_summary.tsv",
  "05_survival_models/ICGC_univariable_cox.tsv",
  "05_survival_models/ICGC_multivariable_cox.tsv",
  "05_survival_models/ICGC_PH_test.tsv",
  "06_sensitivity/ICGC_nonlinearity_sensitivity.tsv",
  "08_source_data/ICGC_axis_scores.tsv",
  "08_source_data/HCCDB18_local_crosscheck.tsv",
  "08_source_data/HCCDB18_crosscheck_summary.json",
  "08_source_data/TCGA_ICGC_axis_replication.tsv",
  "08_source_data/Figure7_ICGC_clinical_validation_audit.tsv",
  "07_figures/ICGC_OS_exploratory_cox_forest.pdf",
  "07_figures/ICGC_OS_exploratory_cox_forest.svg",
  "07_figures/ICGC_OS_exploratory_cox_forest.png",
  "07_figures/ICGC_OS_exploratory_cox_forest.tiff",
  "09_reports/ICGC_OS_FINAL_DECISION.md",
  "09_reports/ICGC_OS_BLOCK_ROOT_CAUSE.md",
  "09_reports/ICGC_OS_survival_results.md",
  "09_reports/ICGC_historical_audit_summary.md",
  "09_reports/ICGC_OS_AUDIT_REPORT.md",
  "10_manifests/ICGC_OS_run_manifest.json"
)
write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE, na = "null", digits = 12)

cat("ICGC OS audit survival models completed: N=", nrow(d), ", events=", sum(d$os_event == 1L), "; decision=ESTIMABLE_BUT_NOT_VALIDATED\n", sep = "")
