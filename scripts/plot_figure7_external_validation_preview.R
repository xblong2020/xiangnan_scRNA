script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "figure7_plot_theme.R")); source(file.path(script_dir, "figure7_core.R")); run_figure7_stage("preview")
