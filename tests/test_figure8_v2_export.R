root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path("scripts", "figure8_v2_theme.R"), local = FALSE)

plot <- ggplot2::ggplot(data.frame(x = 1:2, y = 1:2), ggplot2::aes(x, y)) + ggplot2::geom_line() + figure8_v2_theme()
paths <- figure8_v2_save_plot(plot, "figure8_v2_export_smoke_test", width_mm = 40, height_mm = 30, dpi = 600)
stopifnot(length(paths) == 4L)
stopifnot(all(file.exists(paths)))
stopifnot(all(file.info(paths)$size > 0))
cat("figure8_v2 export smoke test passed\n")
