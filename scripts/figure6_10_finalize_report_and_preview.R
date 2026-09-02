#!/usr/bin/env Rscript

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figure6_common.R"))
suppressPackageStartupMessages(library(patchwork))

plot_paths <- file.path(FIGURE6_METADATA_DIR, paste0("figure6", letters[1:8], "_plot.rds"))
if (!all(file.exists(plot_paths))) stop("Missing panel RDS: ", paste(basename(plot_paths[!file.exists(plot_paths)]), collapse = ", "))
plots <- lapply(plot_paths, readRDS)
preview_plots <- lapply(plots, function(p) wrap_elements(full = p & theme(legend.position = "none", plot.caption = element_blank())))
preview <- wrap_plots(preview_plots, design = "AAAA\nBBCC\nDDEE\nFFGG\nHHHH") +
  plot_annotation(title = "Figure 6 | Directional perturbation of the three regulatory axes",
    tag_levels = "A", theme = theme(plot.title = element_text(size = 13, hjust = .5, family = "sans"),
      plot.tag = element_text(size = 12, face = "bold", family = "sans")))
preview_dir <- file.path(FIGURE6_PROJECT_ROOT, "figures", "driver", "figure6_directional_network_preview")
dir.create(preview_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(preview_dir, "figure6_directional_network_a_to_h_preview.pdf"), preview, width = 16, height = 21, device = cairo_pdf, bg = "white", limitsize = FALSE)
ggsave(file.path(preview_dir, "figure6_directional_network_a_to_h_preview.png"), preview, width = 16, height = 21, dpi = 300, bg = "white", limitsize = FALSE)

effects <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6_perturbation_response_effects.tsv.gz"))
asym <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6c_directional_asymmetry.tsv"))
conc <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6d_cross_method_concordance.tsv"))
fit <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6f_model_fit_summary.tsv"))
edges <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6g_edge_evidence.tsv"))
neg <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6_negative_control_results.tsv.gz"))
conf <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6_confounder_adjustment_summary.tsv"))
path_jac <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6e_pathway_similarity.tsv"))
pkg <- figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6_r_package_versions.tsv"))
selected <- paste(fit[selected == TRUE, model], collapse = "; ")
protected_audit_path <- file.path(FIGURE6_METADATA_DIR, "figure6_protected_asset_change_audit.tsv")
protected_changes <- if (file.exists(protected_audit_path)) figure6_fread(protected_audit_path)[unchanged == FALSE] else data.table()
main_effects <- effects[availability == "Available" & tf %in% c("HNF4A","PPARA","EGR1","CEBPB","AP1_AGGREGATE","SOX4")][order(fdr)][1:min(.N, 12)]

lines <- c(
  "# Figure 6 directional perturbation network report", "",
  "## Reproducible environment", "",
  paste0("- R: ", R.version.string),
  paste0("- Packages: ", paste(paste0(pkg$package, " ", pkg$version), collapse = "; ")),
  paste0("- ggsci Lancet palette: ", paste(lancet_palette, collapse = ", ")),
  paste0("- Axis colours: A ", axis_palette["identity_axis"], "; B ", axis_palette["stress_axis"], "; C ", axis_palette["sox4_axis"]), "",
  "## Frozen interventions and response definitions", "",
  paste0("Available KO perturbations: ", paste(sort(unique(effects[availability == "Available", perturbation])), collapse = ", "), "."),
  "Unavailable interventions: HNF4A restore/OE, PPARA restore/OE and SOX4 OE. No KO was inverted or relabelled as restoration.",
  "AP-1 is the sample-level median of JUN, JUNB, JUND, FOS and ATF3 individual KO effects; it is not a combined suppression simulation.",
  "Programme responses are sample means of CellOracle predicted delta expression over frozen Figure 5 genes, standardized by the matching baseline programme SD.",
  "Malignant fate is the standardized KO-vector projection onto the frozen malignant progression axis. CNV output denotes an expression signature rather than CNV burden.", "",
  "## Figure 6A–E", "",
  paste0("Figure 6A contains ", nrow(effects[availability == "Available"]), " available perturbation-output estimates and ", nrow(effects[availability != "Available"]), " explicitly unavailable cells."),
  paste0("Leading estimates by sample-bootstrap FDR: ", paste(sprintf("%s/%s %.2f (FDR %.3g)", main_effects$perturbation, main_effects$output, main_effects$effect_estimate, main_effects$fdr), collapse = "; "), "."),
  "Figure 6B reuses the frozen audited HNF4A, EGR1 and SOX4 knockout UMAP vector fields with common coordinate and arrow scaling.",
  paste0("Figure 6C classifications: ", paste(paste(asym$comparison, asym$classification, sprintf("D=%.2f [%.2f, %.2f]", asym$directional_asymmetry_score, asym$ci_low, asym$ci_high)), collapse = "; "), "."),
  paste0("Figure 6D median top-50 Jaccard = ", sprintf("%.3f", median(conc$top50_jaccard, na.rm=TRUE)), "; median magnitude-rank correlation = ", sprintf("%.3f", median(conc$spearman_rank_correlation, na.rm=TRUE)), "."),
  "scTenifoldKnk manifold distances are unsigned; Figure 6D uses ranks, magnitudes and gene overlap. Method disagreement is retained.",
  paste0("Figure 6E uses the shared module7.2 main-strict 3,000-gene background. Pathway similarity was ", if (any(path_jac$status == "estimable")) "estimable for perturbation pairs with FDR-significant pathway sets" else "not estimable because no FDR-significant pathway set was available", "."), "",
  "## Figure 6F–H", "",
  paste0("Competing SEMs use ", fit$n_samples[1], " sample-level pseudobulks from ", uniqueN(figure6_fread(file.path(FIGURE6_METADATA_DIR, "figure6f_model_input_sample_pseudobulk.tsv"))$dataset), " datasets. Selected model by the preregistered four-metric rank-sum rule: ", selected, "."),
  "Model assessment combines AIC, BIC, repeated 10×5-fold fate-prediction RMSE, 1,000-bootstrap edge sign stability and leave-one-dataset-out stability. Saturated-model fit indices are not treated as decisive.",
  paste0("Figure 6G grades: ", paste(sprintf("%s→%s %s (%.2f)", edges$source_code, edges$target_code, edges$evidence_grade, edges$evidence_score), collapse = "; "), "."),
  "Edge direction represents computational support and does not establish direct causality.",
  "Figure 6H positions this computational three-axis architecture as complementary and partially overlapping with the experimentally validated FOXM1/CEBPB plasticity axis (Zhang et al., J Hepatol 2026; PMID 41043722). FOXM1 was not perturbed in this project.", "",
  "## Controls, sensitivity and review-risk flags", "",
  paste0("Negative-control pool: ", paste(FIGURE6_AXIS_TFS$control, collapse = ", "), " (n=5). This is below the requested 10–50 controls and creates a specificity risk; minimum empirical p is 1/6."),
  paste0("Candidate effects exceeding all available controls: ", sum(neg$candidate_exceeds_all_controls, na.rm=TRUE), "/", nrow(neg), "."),
  paste0("Confounder-adjusted estimable contrasts: ", sum(conf$status == "estimable"), "/", nrow(conf), ". Dataset, frozen stress, proliferation and CNV-proxy scores were included where estimable."),
  "Hypoxia and separate S/G2M covariates were unavailable and remain untested. Restore/OE evidence is absent. Frozen HCC-malignant programme coverage is sparse. Several temporal-order tests were unsupported.",
  paste0("Protected-asset audit: ", nrow(protected_changes), " Figure 2–5 files differed from the initial baseline. These external/concurrent changes were preserved; Figure 6 outputs remained within dedicated paths. See figure6_protected_asset_change_audit.tsv."),
  "Virtual perturbation methods identified overlapping but incompletely concordant network responses, precluding definitive inference of inter-axis directionality.", "",
  "## Manuscript-ready language", "",
  "Recommended Results subtitle: “Directional perturbation supports a partially ordered regulatory architecture of malignant-state acquisition.”", "",
  "Recommended legend: Figure 6. Directional virtual perturbation of frozen hepatocyte-identity, stress-transition and SOX4-associated programmes. CellOracle KO responses were summarized at sample level with dataset-stratified bootstrap inference; scTenifoldKnk was used only for unsigned magnitude/rank concordance. Competing identifiable SEMs and prespecified evidence weights yielded conservative edge grades. Missing restore/OE interventions and unavailable covariates are shown explicitly.", "",
  "Recommended Results paragraph: CellOracle perturbations produced programme-specific responses across the three frozen axes, while forward–reverse comparisons showed heterogeneous directional asymmetry. Cross-model gene-response concordance was incomplete, and sample-level competing models differed under the preregistered multi-metric selection rule. Integrated evidence therefore supports a conservative, partially ordered computational architecture rather than a strict causal cascade.", "",
  "Most conservative conclusion: Cross-model perturbation analyses supported asymmetric coupling among hepatocyte identity, stress-transition and SOX4-associated malignant-state programmes, although several inter-axis directions remained incompletely resolved.", "",
  "Strongest data-supported conclusion: Directional virtual perturbation supports programme coupling and a partially ordered candidate architecture at sample-pseudobulk resolution.", "",
  "Claims excluded: direct causality; genetic epistasis; actual CNV burden change; true restoration/OE; combined AP-1 suppression; signed scTenifoldKnk effect; superiority over the FOXM1/CEBPB model; SOX4 as the only master regulator.", "",
  "## Main-figure readiness", "",
  "Status: conditionally suitable for an SCI main figure. Visual and computational reproducibility requirements are met; causal wording and specificity claims must remain conservative because restore/OE, a large control pool, hypoxia/S-G2M adjustment and wet-lab validation are unavailable."
)
dir.create(dirname(FIGURE6_REPORT_PATH), recursive = TRUE, showWarnings = FALSE)
writeLines(lines, FIGURE6_REPORT_PATH, useBytes = TRUE)
figure6_write_json(list(
  preview_pdf = figure6_norm_path(file.path(preview_dir, "figure6_directional_network_a_to_h_preview.pdf")),
  preview_png = figure6_norm_path(file.path(preview_dir, "figure6_directional_network_a_to_h_preview.png")),
  report = figure6_norm_path(FIGURE6_REPORT_PATH), selected_model = selected,
  sci_main_figure_status = "conditional", causal_claim_allowed = FALSE
), file.path(FIGURE6_METADATA_DIR, "figure6_finalization_report.json"))
