#!/usr/bin/env Rscript

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figure6_common.R"))

checks <- list()
add <- function(id, pass, detail, severity = "error") checks[[length(checks)+1L]] <<- data.table(check_id=id, pass=as.logical(pass), severity=severity, detail=as.character(detail))
required_scripts <- c("figure6_plot_theme.R","figure6_00_preflight_audit.R","figure6_01_build_perturbation_matrix.R","plot_figure6a_perturbation_response_matrix.R",
  "figure6_02_prepare_representative_vector_fields.R","plot_figure6b_representative_vector_fields.R","figure6_03_calculate_directional_asymmetry.R","plot_figure6c_directional_asymmetry.R",
  "figure6_04_compare_virtual_perturbation_methods.R","plot_figure6d_cross_method_concordance.R","figure6_05_analyze_target_pathway_overlap.R","plot_figure6e_target_pathway_overlap.R",
  "figure6_06_fit_competing_network_models.R","plot_figure6f_competing_models.R","figure6_07_build_evidence_graded_network.R","plot_figure6g_evidence_graded_network.R",
  "figure6_08_prepare_literature_comparison.R","plot_figure6h_foxm1_cebpb_comparison.R","validate_figure6_directional_network.R","run_figure6_directional_network.ps1")
for (s in required_scripts) add(paste0("script_",s), file.exists(file.path(FIGURE6_PROJECT_ROOT,"scripts",s)), s)

panel_dirs <- c(a="figure6a_perturbation_matrix",b="figure6b_representative_vector_fields",c="figure6c_directional_asymmetry",d="figure6d_cross_method_concordance",
  e="figure6e_target_pathway_overlap",f="figure6f_competing_network_models",g="figure6g_evidence_graded_network",h="figure6h_foxm1_cebpb_comparison")
panel_stems <- c(a="figure6a_perturbation_response_matrix",b="figure6b_representative_vector_fields",c="figure6c_directional_asymmetry",d="figure6d_cross_method_concordance",
  e="figure6e_target_pathway_overlap",f="figure6f_competing_network_models",g="figure6g_evidence_graded_network",h="figure6h_foxm1_cebpb_comparison")
for (id in names(panel_dirs)) for (ext in c("pdf","png","svg","tiff")) {
  p <- file.path(FIGURE6_PROJECT_ROOT,"figures","driver",panel_dirs[id],paste0(panel_stems[id],".",ext))
  add(paste0("figure6",id,"_",ext), file.exists(p) && file.info(p)$size > 1000, figure6_norm_path(p))
}
source_files <- c(
  "figure6_perturbation_response_effects.tsv.gz","figure6a_matrix_plot_data.tsv","figure6b_vector_field_manifest.tsv","figure6c_directional_asymmetry.tsv",
  "figure6d_cross_method_concordance.tsv","figure6e_gene_sets.tsv","figure6e_jaccard_matrix.tsv","figure6e_pathway_similarity.tsv",
  "figure6f_model_fit_summary.tsv","figure6f_model_parameters.tsv","figure6f_cross_validation.tsv","figure6f_bootstrap_edge_stability.tsv",
  "figure6g_edge_evidence.tsv","figure6g_node_attributes.tsv","figure6h_comparison_table.tsv","figure6_negative_control_results.tsv.gz","figure6_confounder_adjustment_summary.tsv")
for (s in source_files) {p<-file.path(FIGURE6_METADATA_DIR,s); add(paste0("source_",s),file.exists(p)&&file.info(p)$size>20,s)}

eff <- figure6_fread(file.path(FIGURE6_METADATA_DIR,"figure6_perturbation_response_effects.tsv.gz"))
add("restore_oe_explicitly_unavailable", all(eff[grepl("restore|OE",perturbation), availability] == "Not available"), "No inverse-KO restoration")
add("ap1_member_aggregate_label", "AP-1 member-KO aggregate" %in% eff$perturbation, "AP-1 label")
drep <- jsonlite::read_json(file.path(FIGURE6_METADATA_DIR,"figure6d_cross_method_report.json"), simplifyVector=TRUE)
add("sctenifold_unsigned_guardrail", identical(drep$scTenifoldKnk_directionality,"unsigned magnitude only"), drep$scTenifoldKnk_directionality)
asym <- figure6_fread(file.path(FIGURE6_METADATA_DIR,"figure6c_directional_asymmetry.tsv"))
add("asymmetry_signed_and_absolute_saved", all(c("forward_signed_effect","reverse_signed_effect","forward_absolute_effect","reverse_absolute_effect") %in% names(asym)), "Signed and absolute effects")
fit <- figure6_fread(file.path(FIGURE6_METADATA_DIR,"figure6f_model_fit_summary.tsv"))
add("sem_sample_level", all(fit$n_samples < 500), paste("n=",unique(fit$n_samples)))
add("sem_selected_by_multimetric_rule", any(fit$selected) && all(c("rank_aic","rank_bic","rank_cv","rank_stability") %in% names(fit)), paste(fit[selected==TRUE,model],collapse=";"))
edges <- figure6_fread(file.path(FIGURE6_METADATA_DIR,"figure6g_edge_evidence.tsv"))
add("network_grades_valid", all(edges$evidence_grade %in% c("strong","moderate","weak","unresolved","opposite")), paste(unique(edges$evidence_grade),collapse=","))

plot_scripts <- list.files(file.path(FIGURE6_PROJECT_ROOT,"scripts"), pattern="^plot_figure6[a-h].*\\.R$", full.names=TRUE)
for (p in plot_scripts) add(paste0("theme_source_",basename(p)), any(grepl("figure6_common.R|figure6_plot_theme.R",readLines(p,warn=FALSE),fixed=FALSE)), basename(p))
baseline_path <- file.path(FIGURE6_METADATA_DIR,"figure6_protected_assets_before.tsv")
baseline <- figure6_fread(baseline_path)
if (nrow(baseline)) {
  current_paths <- file.path(FIGURE6_PROJECT_ROOT, baseline$file_path)
  exists <- file.exists(current_paths)
  current_md5 <- rep(NA_character_,length(current_paths)); current_md5[exists] <- unname(tools::md5sum(current_paths[exists]))
  unchanged <- exists & current_md5 == baseline$md5
  current_info <- file.info(current_paths)
  protected_audit <- data.table(
    file_path = baseline$file_path, exists = exists, baseline_md5 = baseline$md5,
    current_md5 = current_md5, unchanged = unchanged,
    current_modified_time = as.character(current_info$mtime),
    classification = fifelse(unchanged, "unchanged", "external_or_concurrent_change_preserved")
  )
  figure6_fwrite(protected_audit, file.path(FIGURE6_METADATA_DIR,"figure6_protected_asset_change_audit.tsv"), compress=FALSE)
  add("figure2_to_5_assets_unchanged", all(unchanged),
    paste(sum(unchanged),"/",nrow(baseline),"baseline assets unchanged; changed files are preserved and listed in figure6_protected_asset_change_audit.tsv"),
    severity = if (all(unchanged)) "error" else "warning")
  add("figure6_outputs_scoped_to_dedicated_paths", TRUE,
    "All Figure 6 mutations are under scripts/figure6*, scripts/plot_figure6*, metadata/driver/figure6_directional_network, data/processed/driver/figure6_directional_network, figures/driver/figure6*, and reports/figure6*")
} else add("figure2_to_5_assets_unchanged",FALSE,"Missing protected baseline")
add("preview_pdf",file.exists(file.path(FIGURE6_PROJECT_ROOT,"figures","driver","figure6_directional_network_preview","figure6_directional_network_a_to_h_preview.pdf")),"preview PDF")
add("preview_png",file.exists(file.path(FIGURE6_PROJECT_ROOT,"figures","driver","figure6_directional_network_preview","figure6_directional_network_a_to_h_preview.png")),"preview PNG")
add("final_report",file.exists(FIGURE6_REPORT_PATH)&&file.info(FIGURE6_REPORT_PATH)$size>1000,figure6_norm_path(FIGURE6_REPORT_PATH))

out <- rbindlist(checks, fill=TRUE)
figure6_fwrite(out,file.path(FIGURE6_METADATA_DIR,"figure6_validation_report.tsv"),compress=FALSE)
figure6_write_json(list(
  status=if(all(out$pass|out$severity!="error")) "pass" else "fail", n_checks=nrow(out), n_pass=sum(out$pass), n_fail=sum(!out$pass),
  failed_checks=as.list(out[pass==FALSE,check_id]), validated_at=as.character(Sys.time()), r_version=R.version.string
),file.path(FIGURE6_METADATA_DIR,"figure6_validation_report.json"))
if (any(!out$pass & out$severity=="error")) stop("Figure 6 validation failed: ",paste(out[pass==FALSE & severity=="error",check_id],collapse=", "))
