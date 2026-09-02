#!/usr/bin/env Rscript

## Figure 2E: Top scTenifoldKnk perturbed genes after SOX4 knockout.
## Main panel uses the malignant-like subset, excludes SOX4 itself, and ranks
## significant genes by manifold-alignment distance.

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
  if (length(hit) == 0L || hit[1] == length(args)) return(default)
  args[hit[1] + 1L]
}

input_path <- normalizePath(
  get_arg("--input", file.path(PROJECT_ROOT, "metadata/driver/sctenifoldknk_module7_2_malignant_like_perturbation_genes.tsv")),
  mustWork = TRUE
)
out_dir <- normalizePath(get_arg("--out-dir", file.path(PROJECT_ROOT, "metadata/driver/figure2e_sox4")), mustWork = FALSE)
figure_dir <- normalizePath(get_arg("--figure-dir", file.path(PROJECT_ROOT, "figures/driver/figure2e_sox4")), mustWork = FALSE)
target_tf <- get_arg("--tf", "SOX4")
top_n <- as.integer(get_arg("--top-n", "20"))
fdr_cutoff <- as.numeric(get_arg("--fdr-cutoff", "0.05"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.delim(input_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("tf", "gene", "distance", "p.adj", "p.value", "Z", "FC", "subset")
missing <- setdiff(required, names(dat))
if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))
dat$distance <- as.numeric(dat$distance)
dat$p.adj <- as.numeric(dat$p.adj)
dat$p.value <- as.numeric(dat$p.value)
dat$Z <- as.numeric(dat$Z)
dat$FC <- as.numeric(dat$FC)

tf_dat <- dat[dat$tf == target_tf & dat$gene != target_tf & is.finite(dat$distance), , drop = FALSE]
significant <- tf_dat[is.finite(tf_dat$p.adj) & tf_dat$p.adj < fdr_cutoff, , drop = FALSE]
if (nrow(significant) < top_n) {
  warning("Fewer than top_n genes pass FDR; filling remaining ranks by distance")
  ranked <- tf_dat[order(tf_dat$distance, decreasing = TRUE), , drop = FALSE]
} else {
  ranked <- significant[order(significant$distance, decreasing = TRUE), , drop = FALSE]
}
top <- head(ranked, top_n)
top$rank <- seq_len(nrow(top))
top$minus_log10_fdr <- -log10(pmax(top$p.adj, .Machine$double.xmin))
top$gene_plot <- factor(top$gene, levels = rev(top$gene))
top$significance <- ifelse(top$p.adj < 0.001, "***", ifelse(top$p.adj < 0.01, "**", "*"))

source_path <- file.path(out_dir, "figure2e_sox4_top20_perturbed_genes.tsv")
write.table(top[, c("rank", "tf", "gene", "distance", "Z", "FC", "p.value", "p.adj", "minus_log10_fdr", "subset")],
            source_path, sep = "\t", quote = FALSE, row.names = FALSE)

lancet <- ggsci::pal_lancet("lanonc")(9)
gradient_colours <- c(lancet[1], lancet[4], lancet[3], lancet[6], lancet[2])
x_max <- max(top$distance, na.rm = TRUE)

p <- ggplot(top, aes(x = distance, y = gene_plot)) +
  geom_col(aes(fill = distance), width = 0.68, alpha = 0.90, colour = NA) +
  geom_point(aes(fill = distance), shape = 21, size = 2.25, stroke = 0.45, colour = "white") +
  geom_text(
    aes(label = significance),
    x = top$distance + x_max * 0.025,
    hjust = 0,
    size = 2.8,
    colour = "grey20"
  ) +
  scale_x_continuous(
    labels = scales::label_scientific(digits = 2),
    expand = expansion(mult = c(0, 0.16))
  ) +
  scale_fill_gradientn(
    colours = gradient_colours,
    values = scales::rescale(c(min(top$distance),
                               quantile(top$distance, 0.25),
                               median(top$distance),
                               quantile(top$distance, 0.75),
                               max(top$distance))),
    labels = scales::label_scientific(digits = 2),
    name = "Distance"
  ) +
  labs(
    x = "Manifold alignment distance",
    y = NULL,
    title = paste0("Top 20 Perturbed Genes (No ", target_tf, ")")
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(size = 10, face = "plain", hjust = 0.5),
    axis.title.x = element_text(size = 8.5),
    axis.text.x = element_text(size = 7.5, colour = "black"),
    axis.text.y = element_text(size = 7.7, colour = "black", face = "plain"),
    axis.line = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.height = grid::unit(0.38, "cm"),
    plot.margin = margin(6, 15, 6, 6)
  )

stem <- file.path(figure_dir, "figure2e_sox4_top20_perturbed_genes")
ggsave(paste0(stem, ".pdf"), p, width = 5.4, height = 5.2, units = "in", device = cairo_pdf)
ggsave(paste0(stem, ".png"), p, width = 5.4, height = 5.2, units = "in", dpi = 600)
ggsave(paste0(stem, ".svg"), p, width = 5.4, height = 5.2, units = "in", device = grDevices::svg)
ggsave(paste0(stem, ".tiff"), p, width = 5.4, height = 5.2, units = "in", dpi = 600, compression = "lzw")

report <- list(
  module = "Figure 2E",
  method = "scTenifoldKnk SOX4 virtual knockout perturbed-gene ranking",
  plotting_language = "R",
  r_version = R.version.string,
  input = input_path,
  subset = unique(top$subset),
  target_tf = target_tf,
  target_gene_excluded = TRUE,
  selection = paste0("p.adj < ", fdr_cutoff, ", then descending manifold-alignment distance"),
  n_tested_excluding_target = nrow(tf_dat),
  n_significant_excluding_target = nrow(significant),
  n_plotted = nrow(top),
  top_genes = as.list(top$gene),
  palette = list(name = "ggsci Lancet lanonc continuous gradient", colours = as.list(unname(gradient_colours))),
  source_table = source_path,
  outputs = list(pdf = paste0(stem, ".pdf"), png = paste0(stem, ".png"), svg = paste0(stem, ".svg"), tiff = paste0(stem, ".tiff")),
  review_risk = "The existing malignant-like scTenifoldKnk run used nc_nNet=1 and nc_nCells=100. This figure is an auditable project-level reproduction; publication claims should be confirmed with nc_nNet=10, nc_nCells=500 and multiple seeds."
)
write_json(report, file.path(out_dir, "figure2e_sox4_report.json"), pretty = TRUE, auto_unbox = TRUE)
message("Figure 2E outputs written to: ", figure_dir)
