#!/usr/bin/env Rscript

## Figure 3F: strict ORA of FDR-significant EGR1-perturbed genes.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsci)
  library(scales)
  library(jsonlite)
})

file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- file_arg[grepl("^--file=", file_arg)]
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else ""
PROJECT_ROOT <- normalizePath(
  if (nzchar(script_path)) file.path(dirname(script_path), "..") else getwd(),
  mustWork = FALSE
)
source(file.path(PROJECT_ROOT, "scripts", "figure3_egr1_common.R"))

data_dir <- normalizePath(
  figure3_get_arg("--data-dir", file.path(PROJECT_ROOT, "metadata/driver/figure3f_egr1")),
  mustWork = FALSE
)
figure_dir <- normalizePath(
  figure3_get_arg("--figure-dir", file.path(PROJECT_ROOT, "figures/driver/figure3f_egr1")),
  mustWork = FALSE
)
subset <- figure3_get_arg("--subset", "stressed_regenerative")
input_path <- normalizePath(
  figure3_get_arg(
    "--input",
    file.path(
      PROJECT_ROOT,
      "metadata/driver/figure3e_egr1",
      paste0("figure3e_egr1_", subset, "_consensus_perturbation_genes.tsv")
    )
  ),
  mustWork = TRUE
)
background_path <- normalizePath(
  figure3_get_arg(
    "--background",
    file.path(
      PROJECT_ROOT,
      "data/processed/driver/figure3e_egr1_sctenifoldknk",
      subset,
      paste0("figure3e_egr1_", subset, "_genes.tsv")
    )
  ),
  mustWork = TRUE
)
gmt_dir <- normalizePath(
  figure3_get_arg(
    "--gmt-dir",
    file.path(PROJECT_ROOT, "metadata/driver/sctenifoldknk_module7_4_genesets")
  ),
  mustWork = TRUE
)
fdr_cutoff <- as.numeric(figure3_get_arg("--fdr-cutoff", "0.05"))
top_n <- as.integer(figure3_get_arg("--top-n", "10"))
min_size <- as.integer(figure3_get_arg("--min-size", "5"))
max_size <- as.integer(figure3_get_arg("--max-size", "500"))
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

perturb <- read.delim(input_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("tf", "gene", "p.adj", "distance", "subset")
missing <- setdiff(required, names(perturb))
if (length(missing)) stop("Missing Figure 3F input columns: ", paste(missing, collapse = ", "))
perturb$p.adj <- as.numeric(perturb$p.adj)
selected <- perturb[
  perturb$tf == "EGR1" & perturb$gene != "EGR1" &
    is.finite(perturb$p.adj) & perturb$p.adj < fdr_cutoff,
  ,
  drop = FALSE
]
input_genes_raw <- unique(toupper(trimws(selected$gene)))
background <- read.delim(background_path, stringsAsFactors = FALSE, check.names = FALSE)[[1]]
background <- unique(toupper(trimws(background)))
input_genes <- intersect(input_genes_raw, background)

gmt_files <- c(
  KEGG = "KEGG_2021_Human.gmt",
  Reactome = "Reactome_2022.gmt",
  GO_BP = "GO_Biological_Process_2023.gmt"
)

read_gmt <- function(path, database) {
  lines <- readLines(path, warn = FALSE)
  rows <- lapply(lines, function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 3L) return(NULL)
    genes <- unique(toupper(fields[3:length(fields)]))
    genes <- intersect(genes[nzchar(genes)], background)
    if (length(genes) < min_size || length(genes) > max_size) return(NULL)
    overlap <- intersect(input_genes, genes)
    pvalue <- if (length(input_genes)) {
      phyper(
        length(overlap) - 1L,
        length(genes),
        length(background) - length(genes),
        length(input_genes),
        lower.tail = FALSE
      )
    } else {
      1
    }
    data.frame(
      database = database,
      term = fields[1],
      overlap_count = length(overlap),
      term_size = length(genes),
      input_gene_count = length(input_genes),
      background_gene_count = length(background),
      pvalue = pvalue,
      overlap_genes = paste(sort(overlap), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

all_results <- do.call(rbind, lapply(names(gmt_files), function(database) {
  read_gmt(file.path(gmt_dir, gmt_files[[database]]), database)
}))
if (!nrow(all_results)) {
  all_results <- data.frame(
    database = character(),
    term = character(),
    overlap_count = integer(),
    term_size = integer(),
    input_gene_count = integer(),
    background_gene_count = integer(),
    pvalue = numeric(),
    overlap_genes = character(),
    stringsAsFactors = FALSE
  )
}
if (nrow(all_results)) {
  # A single BH family across all three databases is conservative and auditable.
  all_results$p.adjust <- p.adjust(all_results$pvalue, method = "BH")
  all_results$gene_ratio <- ifelse(
    all_results$input_gene_count > 0,
    all_results$overlap_count / all_results$input_gene_count,
    0
  )
  all_results$minus_log10_pvalue <- -log10(pmax(all_results$pvalue, .Machine$double.xmin))
  all_results$significant <- all_results$p.adjust < fdr_cutoff
  all_results$subset <- subset
  all_results$tf <- "EGR1"
  all_results <- all_results[order(all_results$p.adjust, all_results$pvalue, all_results$term), , drop = FALSE]
} else {
  all_results$p.adjust <- numeric()
  all_results$gene_ratio <- numeric()
  all_results$minus_log10_pvalue <- numeric()
  all_results$significant <- logical()
  all_results$subset <- character()
  all_results$tf <- character()
}

significant_pathways <- all_results[
  all_results$significant %in% TRUE & all_results$overlap_count > 0,
  ,
  drop = FALSE
]
plot_data <- head(significant_pathways, top_n)
if (nrow(plot_data)) {
  plot_data <- plot_data[
    order(plot_data$minus_log10_pvalue, decreasing = TRUE),
    ,
    drop = FALSE
  ]
  plot_data$rank <- seq_len(nrow(plot_data))
  plot_data$term_label <- vapply(
    paste0(plot_data$term, " [", plot_data$database, "]"),
    function(value) paste(strwrap(value, width = 42), collapse = "\n"),
    character(1)
  )
  plot_data$term_plot <- factor(plot_data$term_label, levels = rev(plot_data$term_label))
}

all_path <- file.path(data_dir, "figure3f_egr1_enrichment_all.tsv")
plot_path <- file.path(data_dir, "figure3f_egr1_plot_data.tsv")
write.table(all_results, all_path, sep = "\t", quote = FALSE, row.names = FALSE)
if (nrow(plot_data)) {
  plot_export <- plot_data[, setdiff(names(plot_data), c("term_label", "term_plot")), drop = FALSE]
  write.table(plot_export, plot_path, sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  write.table(significant_pathways, plot_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

formal_plot_generated <- nrow(plot_data) > 0
stem <- file.path(figure_dir, "figure3f_egr1_pathway_enrichment")
if (formal_plot_generated) {
  lancet <- ggsci::pal_lancet("lanonc")(9)
  gradient_colours <- c(lancet[5], lancet[1], lancet[4])
  x_max <- max(plot_data$minus_log10_pvalue)
  p <- ggplot(plot_data, aes(x = minus_log10_pvalue, y = term_plot)) +
    geom_segment(
      aes(x = 0, xend = minus_log10_pvalue, yend = term_plot),
      linewidth = 0.65,
      colour = "#B8C2CC",
      lineend = "round"
    ) +
    geom_point(aes(size = overlap_count, colour = minus_log10_pvalue), alpha = 0.96) +
    scale_colour_gradientn(
      colours = gradient_colours,
      name = expression(-log[10](italic(P)))
    ) +
    scale_size_continuous(
      range = c(3.0, 7.0),
      breaks = pretty_breaks(n = 3),
      name = "Overlap count"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.10))) +
    labs(
      x = expression(-log[10](italic(P)*"-value")),
      y = NULL,
      title = "Pathway enrichment",
      tag = "Figure 3F"
    ) +
    coord_cartesian(xlim = c(0, x_max * 1.05), clip = "off") +
    figure3_theme() +
    theme(
      axis.text.y = element_text(size = 7.4, colour = "black", lineheight = 0.92),
      legend.position = "right",
      plot.margin = margin(7, 9, 7, 7)
    )
  figure3_save(p, figure_dir, "figure3f_egr1_pathway_enrichment", 6.4, 4.8, tiff = TRUE)
}

review_risks <- list()
recommendation <- NULL
if (!formal_plot_generated) {
  review_risks <- list(list(
    flag = "no_fdr_significant_pathways",
    severity = "main_panel_blocking",
    detail = "No pathway passed the global BH p.adjust < 0.05 rule; the formal Figure 3F plot was suppressed."
  ))
  recommendation <- "Move pathway enrichment to Extended Data or replace Figure 3F with a stress-transition module score panel."
}
report <- list(
  module = "Figure 3F",
  method = "One-sided hypergeometric ORA using local GMT gene sets and global BH correction",
  target_tf = "EGR1",
  target_gene_excluded = TRUE,
  subset = subset,
  input = figure3_norm_path(input_path),
  background = figure3_norm_path(background_path),
  background_matches_figure3e_subset = TRUE,
  databases = as.list(names(gmt_files)),
  fdr_cutoff = fdr_cutoff,
  n_significant_perturbed_genes = length(input_genes),
  n_background_genes = length(background),
  n_pathways_tested = nrow(all_results),
  n_significant_pathways = nrow(significant_pathways),
  formal_plot_generated = formal_plot_generated,
  nominal_results_used_as_formal = FALSE,
  display_rule = "Up to 10 globally BH FDR-significant pathways ordered by adjusted P value, then nominal P value",
  n_plotted = nrow(plot_data),
  plotted_pathways = as.list(as.character(plot_data$term)),
  outputs = c(
    list(
      all_results = figure3_norm_path(all_path),
      plot_data = figure3_norm_path(plot_path)
    ),
    if (formal_plot_generated) list(
      pdf = figure3_norm_path(paste0(stem, ".pdf")),
      png = figure3_norm_path(paste0(stem, ".png")),
      svg = figure3_norm_path(paste0(stem, ".svg")),
      tiff = figure3_norm_path(paste0(stem, ".tiff"))
    ) else list()
  ),
  review_risk_flags = review_risks,
  recommendation_if_no_stable_pathways = recommendation,
  caveat = "ORA reports enrichment among unsigned network-displacement genes. It does not label pathways as activated, suppressed, upregulated, or downregulated."
)
figure3_write_json(report, file.path(data_dir, "figure3f_egr1_report.json"))
message(
  "Figure 3F report written; FDR-significant pathways=", nrow(significant_pathways),
  ", formal plot generated=", formal_plot_generated
)

