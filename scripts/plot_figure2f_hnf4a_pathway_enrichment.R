#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(ggplot2); library(ggsci); library(scales)})
source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])),
                 "figure2_hnf4a_common.R"))
root <- figure2_project_root()
input_path <- figure2_get_arg("--input", file.path(root,
  "metadata/driver/figure2e_hnf4a_sctenifoldknk/figure2e_hnf4a_normal_reference_perturbation_genes.tsv"))
background_path <- figure2_get_arg("--background", file.path(root,
  "data/processed/driver/figure2e_hnf4a_sctenifoldknk/normal_reference/figure2e_hnf4a_normal_reference_genes.tsv"))
gmt_dir <- figure2_get_arg("--gmt-dir", file.path(root, "metadata/driver/sctenifoldknk_module7_4_genesets"))
out_dir <- figure2_get_arg("--out-dir", file.path(root, "metadata/driver/figure2f_hnf4a"))
figure_dir <- figure2_get_arg("--figure-dir", file.path(root, "figures/driver/figure2f_hnf4a"))
target_tf <- figure2_get_arg("--target-tf", "HNF4A")
fdr <- as.numeric(figure2_get_arg("--fdr-cutoff", "0.05"))
top_n <- as.integer(figure2_get_arg("--top-n", "10"))
min_size <- as.integer(figure2_get_arg("--min-size", "5"))
max_size <- as.integer(figure2_get_arg("--max-size", "500"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

perturb <- read.delim(input_path, stringsAsFactors = FALSE, check.names = FALSE)
perturb$p.adj <- as.numeric(perturb$p.adj)
selected <- perturb[perturb$tf == target_tf & perturb$gene != target_tf &
  is.finite(perturb$p.adj) & perturb$p.adj < fdr, , drop = FALSE]
background <- unique(toupper(trimws(read.delim(background_path, stringsAsFactors = FALSE)[[1]])))
input_genes <- intersect(unique(toupper(trimws(selected$gene))), background)
gmt_files <- c(KEGG = "KEGG_2021_Human.gmt", Reactome = "Reactome_2022.gmt",
               GO_BP = "GO_Biological_Process_2023.gmt")

read_gmt <- function(path, database) {
  rows <- lapply(readLines(path, warn = FALSE), function(line) {
    f <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(f) < 3) return(NULL)
    genes <- intersect(unique(toupper(f[3:length(f)])), background)
    if (length(genes) < min_size || length(genes) > max_size) return(NULL)
    overlap <- intersect(input_genes, genes)
    data.frame(database = database, term = f[1], overlap_count = length(overlap),
      term_size = length(genes), input_gene_count = length(input_genes),
      background_gene_count = length(background),
      pvalue = if (length(input_genes)) phyper(length(overlap) - 1L, length(genes),
        length(background) - length(genes), length(input_genes), lower.tail = FALSE) else 1,
      overlap_genes = paste(sort(overlap), collapse = ";"), stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}
all_results <- do.call(rbind, lapply(names(gmt_files), function(db)
  read_gmt(file.path(gmt_dir, gmt_files[[db]]), db)))
if (is.null(all_results) || !nrow(all_results)) {
  all_results <- data.frame(database = character(), term = character(), overlap_count = integer(),
    term_size = integer(), input_gene_count = integer(), background_gene_count = integer(),
    pvalue = numeric(), overlap_genes = character())
}
all_results$p.adjust <- p.adjust(all_results$pvalue, method = "BH")
all_results$minus_log10_pvalue <- -log10(pmax(all_results$pvalue, .Machine$double.xmin))
all_results$significant <- all_results$p.adjust < fdr
all_results$tf <- target_tf; all_results$subset <- "normal_reference"
all_results <- all_results[order(all_results$p.adjust, all_results$pvalue, all_results$term), ]
sig <- all_results[all_results$significant, , drop = FALSE]
plot_dat <- head(sig, top_n)
if (nrow(plot_dat)) {
  plot_dat$term_label <- vapply(plot_dat$term, function(x) paste(strwrap(x, 38), collapse = "\n"), character(1))
  plot_dat$term_plot <- factor(plot_dat$term_label, levels = rev(plot_dat$term_label))
}
all_path <- file.path(out_dir, "figure2f_hnf4a_enrichment_all.tsv")
plot_path <- file.path(out_dir, "figure2f_hnf4a_plot_data.tsv")
write.table(all_results, all_path, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(plot_dat[, setdiff(names(plot_dat), c("term_label", "term_plot")), drop = FALSE],
            plot_path, sep = "\t", quote = FALSE, row.names = FALSE)

figure_generated <- FALSE
if (nrow(plot_dat)) {
  lancet <- ggsci::pal_lancet("lanonc")(9)
  cols <- c(lancet[5], lancet[1], lancet[4])
  p <- ggplot(plot_dat, aes(minus_log10_pvalue, term_plot)) +
    geom_segment(aes(x = 0, xend = minus_log10_pvalue, yend = term_plot),
                 linewidth = .65, colour = "#B8C2CC", lineend = "round") +
    geom_point(aes(size = overlap_count, colour = minus_log10_pvalue), alpha = .96) +
    scale_colour_gradientn(colours = cols, name = expression(-log[10](italic(P)))) +
    scale_size_continuous(range = c(3, 7), breaks = pretty_breaks(3), name = "Overlap count") +
    labs(x = expression(-log[10](italic(P)*"-value")), y = NULL,
         title = "HNF4A pathway enrichment", tag = "Figure 2F") +
    coord_cartesian(clip = "off") + figure2_theme() +
    theme(axis.text.y = element_text(size = 7.4, lineheight = .92))
  figure2_save(p, figure_dir, "figure2f_hnf4a_pathway_enrichment", 6.4, 4.8, tiff = TRUE)
  figure_generated <- TRUE
}
recommendation <- if (!nrow(sig))
  "Move pathway enrichment to Extended Data or replace Figure 2F with an HNF4A target-module score panel." else
  "Retain the FDR-significant ORA panel."
report <- list(
  module = "Figure 2F", target_tf = target_tf,
  method = "One-sided hypergeometric ORA with pooled BH correction across local GMT pathways",
  input = figure2_norm_path(input_path), background = figure2_norm_path(background_path),
  subset = "normal_reference", target_gene_excluded = TRUE, fdr_cutoff = fdr,
  n_significant_perturbed_genes = length(input_genes), n_background_genes = length(background),
  n_pathways_tested = nrow(all_results), n_significant_pathways = nrow(sig),
  display_rule = "Only BH FDR-significant pathways are eligible for the main plot",
  n_plotted = nrow(plot_dat), figure_generated = figure_generated,
  enrichment_interpretation = "ORA of non-directional network displacement genes",
  recommendation = recommendation,
  outputs = list(all_results = figure2_norm_path(all_path), plot_data = figure2_norm_path(plot_path)),
  caveat = "Manifold distance has no activation/suppression direction; nominal pathways are not displayed as a positive main result."
)
figure2_write_json(report, file.path(out_dir, "figure2f_hnf4a_report.json"))
