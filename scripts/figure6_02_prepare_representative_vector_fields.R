#!/usr/bin/env Rscript

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figure6_common.R"))

spec <- data.table(
  axis = c("identity_axis", "stress_axis", "sox4_axis"),
  axis_label = c("Axis A", "Axis B", "Axis C"),
  perturbation = c("HNF4A knockout", "EGR1 knockout", "SOX4 knockout"),
  stem = c("figure2c_hnf4a", "figure3c_egr1", "figure2c_sox4"),
  folder = c("figure2c_hnf4a", "figure3c_egr1", "figure2c_sox4"),
  source_figure = c("Figure 2", "Figure 3", "Figure 4 (legacy Figure 2 filename)")
)
all_cells <- list(); all_grid <- list()
for (i in seq_len(nrow(spec))) {
  s <- spec[i]
  cells_path <- file.path(FIGURE6_PROJECT_ROOT, "metadata", "driver", s$folder, paste0(s$stem, "_matched_cells.tsv.gz"))
  grid_path <- file.path(FIGURE6_PROJECT_ROOT, "metadata", "driver", s$folder, paste0(s$stem, "_grid_umap.tsv.gz"))
  cells <- figure6_fread(cells_path)[, .(cell_id, celloracle_state, umap_1, umap_2)]
  grid <- figure6_fread(grid_path)
  cells[, `:=`(axis = s$axis, axis_label = s$axis_label, perturbation = s$perturbation)]
  grid[, `:=`(axis = s$axis, axis_label = s$axis_label, perturbation = s$perturbation)]
  all_cells[[i]] <- cells; all_grid[[i]] <- grid
  spec[i, `:=`(cells_path = figure6_norm_path(cells_path), grid_path = figure6_norm_path(grid_path),
    n_cells = nrow(cells), n_grid = nrow(grid), n_arrows = sum(grid$show %in% TRUE, na.rm = TRUE))]
}
cells <- rbindlist(all_cells, fill = TRUE)
grid <- rbindlist(all_grid, fill = TRUE)
shared_arrow_length <- stats::median(grid$arrow_length[is.finite(grid$arrow_length)], na.rm = TRUE)
grid[, `:=`(
  plot_xend = grid_x + flow_norm_x * shared_arrow_length,
  plot_yend = grid_y + flow_norm_y * shared_arrow_length,
  shared_arrow_length = shared_arrow_length
)]
figure6_fwrite(cells, file.path(FIGURE6_METADATA_DIR, "figure6b_vector_field_cells.tsv.gz"))
figure6_fwrite(grid, file.path(FIGURE6_METADATA_DIR, "figure6b_vector_field_grid.tsv.gz"))
figure6_fwrite(spec, file.path(FIGURE6_METADATA_DIR, "figure6b_vector_field_manifest.tsv"), compress = FALSE)
figure6_write_json(list(
  panel = "Figure 6B", coordinate = "CellOracle UMAP", representative_fields = as.list(spec$perturbation),
  shared_arrow_length = shared_arrow_length, restore_available = FALSE,
  guardrail = "HNF4A KO remains labelled and interpreted as knockout; no field was selected anew."
), file.path(FIGURE6_METADATA_DIR, "figure6b_vector_field_report.json"))

