#!/usr/bin/env Rscript
root <- normalizePath(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts", "figure8_plot_theme.R"))
source(file.path(root, "scripts", "figure8_analysis_functions.R"))
run_figure8_prepare_signature()

