#!/usr/bin/env Rscript

## Figure 2F: pathway enrichment of malignant-like genes perturbed by SOX4 KO.
## Uses local GMT files and a one-sided hypergeometric test for reproducibility.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsci)
  library(scales)
  library(jsonlite)
})

file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- file_arg[grepl("^--file=", file_arg)]
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else ""
PROJECT_ROOT <- normalizePath(if (nzchar(script_path)) file.path(dirname(script_path), "..") else getwd(), mustWork = FALSE)
if (!dir.exists(file.path(PROJECT_ROOT, "scripts"))) PROJECT_ROOT <- normalizePath(getwd(), mustWork = TRUE)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  hit <- which(args == flag)
  if (!length(hit) || hit[1] == length(args)) return(default)
  args[hit[1] + 1L]
}

input_path <- normalizePath(get_arg("--input", file.path(PROJECT_ROOT,
  "metadata/driver/sctenifoldknk_module7_2_malignant_like_perturbation_genes.tsv")), mustWork = TRUE)
background_path <- normalizePath(get_arg("--background", file.path(PROJECT_ROOT,
  "data/processed/driver/sctenifoldknk_module7_1/malignant_like/sctenifoldknk_genes.tsv")), mustWork = TRUE)
gmt_dir <- normalizePath(get_arg("--gmt-dir", file.path(PROJECT_ROOT,
  "metadata/driver/sctenifoldknk_module7_4_genesets")), mustWork = TRUE)
out_dir <- normalizePath(get_arg("--out-dir", file.path(PROJECT_ROOT,
  "metadata/driver/figure2f_sox4")), mustWork = FALSE)
figure_dir <- normalizePath(get_arg("--figure-dir", file.path(PROJECT_ROOT,
  "figures/driver/figure2f_sox4")), mustWork = FALSE)
target_tf <- get_arg("--tf", "SOX4")
fdr_cutoff <- as.numeric(get_arg("--fdr-cutoff", "0.05"))
top_n <- as.integer(get_arg("--top-n", "10"))
min_size <- as.integer(get_arg("--min-size", "5"))
max_size <- as.integer(get_arg("--max-size", "500"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

perturb <- read.delim(input_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("tf", "gene", "p.adj", "subset")
missing <- setdiff(required, names(perturb))
if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))
perturb$p.adj <- as.numeric(perturb$p.adj)
selected <- perturb[
  perturb$tf == target_tf & perturb$gene != target_tf &
    is.finite(perturb$p.adj) & perturb$p.adj < fdr_cutoff,
  , drop = FALSE
]
input_genes <- unique(toupper(trimws(selected$gene)))
background <- read.delim(background_path, stringsAsFactors = FALSE, check.names = FALSE)[[1]]
background <- unique(toupper(trimws(background)))
input_genes <- intersect(input_genes, background)
if (length(input_genes) < min_size) stop("Too few significant genes in the network background: ", length(input_genes))

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
    if (!length(overlap)) return(NULL)
    data.frame(
      database = database,
      term = fields[1],
      overlap_count = length(overlap),
      term_size = length(genes),
      input_gene_count = length(input_genes),
      background_gene_count = length(background),
      pvalue = phyper(length(overlap) - 1L, length(genes),
        length(background) - length(genes), length(input_genes), lower.tail = FALSE),
      overlap_genes = paste(sort(overlap), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  out$p.adjust <- p.adjust(out$pvalue, method = "BH")
  out$gene_ratio <- out$overlap_count / out$input_gene_count
  out$minus_log10_pvalue <- -log10(pmax(out$pvalue, .Machine$double.xmin))
  out$significant <- out$p.adjust < fdr_cutoff
  out[order(out$pvalue, out$term), , drop = FALSE]
}

all_results <- do.call(rbind, lapply(names(gmt_files), function(db) {
  read_gmt(file.path(gmt_dir, gmt_files[[db]]), db)
}))
if (!nrow(all_results)) stop("No enriched pathways had at least one overlapping gene")
all_results$subset <- "malignant_like"
all_results$tf <- target_tf
all_results <- all_results[order(all_results$database, all_results$pvalue, all_results$term), ]

# Prefer significant KEGG pathways; supplement with significant Reactome/GO-BP.
priority <- c(KEGG = 1L, Reactome = 2L, GO_BP = 3L)
sig <- all_results[all_results$significant, , drop = FALSE]
candidates <- if (nrow(sig)) sig else all_results
candidates$db_priority <- unname(priority[candidates$database])
candidates <- candidates[order(candidates$db_priority, candidates$pvalue, candidates$term), , drop = FALSE]
plot_dat <- head(candidates, top_n)
plot_dat <- plot_dat[order(plot_dat$minus_log10_pvalue, decreasing = TRUE), , drop = FALSE]
plot_dat$rank <- seq_len(nrow(plot_dat))
plot_dat$term_label <- vapply(plot_dat$term, function(x) paste(strwrap(x, width = 38), collapse = "\n"), character(1))
plot_dat$term_plot <- factor(plot_dat$term_label, levels = rev(plot_dat$term_label))

all_path <- file.path(out_dir, "figure2f_sox4_enrichment_all.tsv")
plot_path <- file.path(out_dir, "figure2f_sox4_plot_data.tsv")
write.table(all_results, all_path, sep = "\t", quote = FALSE, row.names = FALSE)
plot_export <- plot_dat[, setdiff(names(plot_dat), c("term_label", "term_plot")), drop = FALSE]
write.table(plot_export, plot_path, sep = "\t", quote = FALSE, row.names = FALSE)

lancet <- ggsci::pal_lancet("lanonc")(9)
# Purple -> Lancet blue -> cyan; deliberately distinct from Figure 2E's multihue bars.
gradient_colours <- c(lancet[5], lancet[1], lancet[4])
x_max <- max(plot_dat$minus_log10_pvalue)

p <- ggplot(plot_dat, aes(x = minus_log10_pvalue, y = term_plot)) +
  geom_segment(aes(x = 0, xend = minus_log10_pvalue, yend = term_plot),
    linewidth = 0.65, colour = "#B8C2CC", lineend = "round") +
  geom_point(aes(size = overlap_count, colour = minus_log10_pvalue), alpha = 0.96) +
  scale_colour_gradientn(colours = gradient_colours, name = expression(-log[10](italic(P)))) +
  scale_size_continuous(range = c(3.0, 7.0), breaks = pretty_breaks(n = 3), name = "Gene count") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(x = expression(-log[10](italic(P)*"-value")), y = NULL, title = "Pathway enrichment") +
  coord_cartesian(xlim = c(0, x_max * 1.05), clip = "off") +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(size = 10, face = "plain", hjust = 0.5),
    axis.title.x = element_text(size = 8.5),
    axis.text.x = element_text(size = 7.5, colour = "black"),
    axis.text.y = element_text(size = 7.4, colour = "black", lineheight = 0.92),
    axis.line = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    legend.position = "right",
    legend.title = element_text(size = 7.5),
    legend.text = element_text(size = 7),
    legend.key.height = grid::unit(0.38, "cm"),
    plot.margin = margin(7, 9, 7, 7)
  )

stem <- file.path(figure_dir, "figure2f_sox4_pathway_enrichment")
ggsave(paste0(stem, ".pdf"), p, width = 6.4, height = 4.8, units = "in", device = cairo_pdf)
ggsave(paste0(stem, ".png"), p, width = 6.4, height = 4.8, units = "in", dpi = 600)
ggsave(paste0(stem, ".svg"), p, width = 6.4, height = 4.8, units = "in", device = grDevices::svg)
ggsave(paste0(stem, ".tiff"), p, width = 6.4, height = 4.8, units = "in", dpi = 600, compression = "lzw")

report <- list(
  module = "Figure 2F",
  method = "One-sided hypergeometric ORA using local GMT gene sets",
  plotting_language = "R",
  r_version = R.version.string,
  input = enc2utf8(gsub("\\\\", "/", input_path)),
  background = enc2utf8(gsub("\\\\", "/", background_path)),
  target_tf = target_tf,
  subset = "malignant_like",
  target_gene_excluded = TRUE,
  fdr_cutoff = fdr_cutoff,
  n_significant_perturbed_genes = length(input_genes),
  n_background_genes = length(background),
  n_pathways_tested = nrow(all_results),
  n_significant_pathways = sum(all_results$significant),
  display_rule = if (nrow(sig)) "Database priority KEGG > Reactome > GO_BP among BH FDR-significant pathways" else "No BH-significant pathways; database-prioritized nominal results shown",
  n_plotted = nrow(plot_dat),
  plotted_pathways = as.list(as.character(plot_dat$term)),
  palette = list(name = "ggsci Lancet lanonc purple-blue-cyan continuous gradient", colours = as.list(unname(gradient_colours))),
  outputs = list(
    all_results = enc2utf8(gsub("\\\\", "/", all_path)),
    plot_data = enc2utf8(gsub("\\\\", "/", plot_path)),
    pdf = enc2utf8(gsub("\\\\", "/", paste0(stem, ".pdf"))),
    png = enc2utf8(gsub("\\\\", "/", paste0(stem, ".png"))),
    svg = enc2utf8(gsub("\\\\", "/", paste0(stem, ".svg"))),
    tiff = enc2utf8(gsub("\\\\", "/", paste0(stem, ".tiff")))
  ),
  review_risk = "The source malignant-like scTenifoldKnk run used nc_nNet=1 and nc_nCells=100; publication-level claims require the planned higher-replicate rerun."
)
write_json(report, file.path(out_dir, "figure2f_sox4_report.json"), pretty = TRUE, auto_unbox = TRUE)
message("Figure 2F outputs written to: ", figure_dir)
