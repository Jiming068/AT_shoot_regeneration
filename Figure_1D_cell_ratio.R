# Portable path configuration (override with environment variables).
project_dir <- Sys.getenv("SHOOT_REGEN_PROJECT_DIR", unset = ".")
output_root <- Sys.getenv("SHOOT_REGEN_OUTPUT_DIR", unset = file.path(project_dir, "02_results"))
external_data_dir <- Sys.getenv("SHOOT_REGEN_EXTERNAL_DATA_DIR", unset = file.path(project_dir, "external_data"))
work_dir <- Sys.getenv("SHOOT_REGEN_WORK_DIR", unset = ".")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

# Cell-ratio changes over time for WT clusters 0, 2 and 3.
# Proportions are calculated against all cells present at each time point,
# matching the denominator convention of cell_ratio_change_barplot.R.

input_rds <- file.path(project_dir, "01_obj/01_WT_snRNA/02_snRNA_subset_res0.2_SIM4-SIM16.rds")
reference_script <- file.path(project_dir, "03_script/cell_ratio_change_barplot.R")
default_outdir <- file.path(project_dir, "02_results/68_WT_C0_C2_C3_cell_ratio_by_time")
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1L) args[[1L]] else default_outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

target_clusters <- c("0", "2", "3")
time_levels <- c("SIM4", "SIM6", "SIM8", "SIM10", "SIM14", "SIM16")

# Exact colors retained from cell_ratio_change_barplot.R.
cluster_colors <- c(
  "0" = "#FFB74D",
  "2" = "#81C784",
  "3" = "#4C78A8"
)

save_plot <- function(plot, stem, width = 4.35, height = 3.45) {
  ggsave(
    file.path(outdir, paste0(stem, ".pdf")), plot,
    width = width, height = height, device = cairo_pdf,
    bg = "white", limitsize = FALSE
  )
  svg_device <- if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite
  } else {
    grDevices::svg
  }
  ggsave(
    file.path(outdir, paste0(stem, ".svg")), plot,
    width = width, height = height, device = svg_device,
    bg = "white", limitsize = FALSE
  )
  ggsave(
    file.path(outdir, paste0(stem, ".png")), plot,
    width = width, height = height, dpi = 450,
    bg = "white", limitsize = FALSE
  )
}

obj <- readRDS(input_rds)
required_metadata <- c("seurat_clusters", "date")
missing_metadata <- setdiff(required_metadata, colnames(obj@meta.data))
if (length(missing_metadata)) {
  stop("Missing required metadata: ", paste(missing_metadata, collapse = ", "))
}

cell_meta <- data.frame(
  cell = colnames(obj),
  cluster = as.character(obj$seurat_clusters),
  date_original = as.character(obj$date),
  stringsAsFactors = FALSE
) %>%
  mutate(
    date = date_original,
    date = factor(date, levels = time_levels)
  )

if (anyNA(cell_meta$date)) {
  bad_dates <- unique(cell_meta$date_original[is.na(cell_meta$date)])
  stop("Unexpected date labels: ", paste(bad_dates, collapse = ", "))
}
missing_clusters <- setdiff(target_clusters, unique(cell_meta$cluster))
if (length(missing_clusters)) {
  stop("Target clusters missing from input: ", paste(missing_clusters, collapse = ", "))
}

total_by_date <- cell_meta %>%
  count(date, name = "total_cells_all_clusters", .drop = FALSE)

ratio_long <- cell_meta %>%
  filter(cluster %in% target_clusters) %>%
  mutate(cluster = factor(cluster, levels = target_clusters)) %>%
  count(date, cluster, name = "n_cells", .drop = FALSE) %>%
  complete(
    date = factor(time_levels, levels = time_levels),
    cluster = factor(target_clusters, levels = target_clusters),
    fill = list(n_cells = 0L)
  ) %>%
  left_join(total_by_date, by = "date") %>%
  mutate(
    proportion_of_all_cells = n_cells / total_cells_all_clusters,
    percent_of_all_cells = 100 * proportion_of_all_cells
  ) %>%
  arrange(date, cluster)

selected_totals <- ratio_long %>%
  group_by(date) %>%
  summarise(
    selected_cells = sum(n_cells),
    total_cells_all_clusters = first(total_cells_all_clusters),
    selected_fraction_of_all_cells = selected_cells / total_cells_all_clusters,
    selected_percent_of_all_cells = 100 * selected_fraction_of_all_cells,
    .groups = "drop"
  )

ratio_within_selected <- ratio_long %>%
  left_join(selected_totals %>% select(date, selected_cells), by = "date") %>%
  mutate(
    proportion_within_C0_C2_C3 = ifelse(selected_cells > 0, n_cells / selected_cells, 0),
    percent_within_C0_C2_C3 = 100 * proportion_within_C0_C2_C3
  )

write.csv(
  ratio_within_selected,
  file.path(outdir, "02_C0_C2_C3_cell_counts_and_proportions_by_date_long.csv"),
  row.names = FALSE
)
write.csv(
  selected_totals,
  file.path(outdir, "03_C0_C2_C3_combined_fraction_by_date.csv"),
  row.names = FALSE
)

count_wide <- ratio_long %>%
  select(cluster, date, n_cells) %>%
  pivot_wider(names_from = date, values_from = n_cells)
proportion_wide <- ratio_long %>%
  select(cluster, date, percent_of_all_cells) %>%
  pivot_wider(names_from = date, values_from = percent_of_all_cells)
write.csv(count_wide,
          file.path(outdir, "04_C0_C2_C3_cell_counts_by_date_wide.csv"),
          row.names = FALSE)
write.csv(proportion_wide,
          file.path(outdir, "05_C0_C2_C3_percent_of_all_cells_by_date_wide.csv"),
          row.names = FALSE)

base_theme <- theme_classic(base_size = 10.5, base_family = "Arial") +
  theme(
    axis.title.y = element_text(face = "bold", size = 10.5),
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black", size = 9.5),
    axis.line = element_line(color = "black", linewidth = 0.45),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    legend.title = element_text(face = "bold", size = 9.5),
    legend.text = element_text(size = 9.5),
    legend.key.height = grid::unit(0.38, "cm"),
    legend.key.width = grid::unit(0.38, "cm"),
    plot.margin = margin(5, 5, 4, 5)
  )

# Main publication figure: the denominator is all cells at the corresponding
# time point. reverse=TRUE ensures the visual bottom-to-top order is 0, 2, 3.
p_main <- ggplot(
  ratio_long,
  aes(x = date, y = proportion_of_all_cells, fill = cluster)
) +
  geom_col(
    position = position_stack(reverse = TRUE),
    width = 0.70, color = "white", linewidth = 0.35
  ) +
  scale_fill_manual(
    values = cluster_colors, breaks = target_clusters,
    drop = FALSE, name = "Cluster"
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 0.80), breaks = seq(0, 0.8, 0.2),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(y = "Proportion of all cells") +
  base_theme

save_plot(p_main, "06_C0_C2_C3_cell_ratio_by_time_stacked_barplot")

# Companion labeled version for direct numerical reading.
p_labeled <- p_main +
  geom_text(
    aes(label = ifelse(percent_of_all_cells >= 2,
                       sprintf("%.1f%%", percent_of_all_cells), "")),
    position = position_stack(vjust = 0.5, reverse = TRUE),
    family = "Arial", fontface = "bold", size = 3.0,
    color = "#202020", show.legend = FALSE
  )
save_plot(p_labeled, "07_C0_C2_C3_cell_ratio_by_time_stacked_barplot_labeled")

# A normalized companion plot shows redistribution only within the selected
# C0/C2/C3 compartment; this is not used as the primary abundance figure.
p_selected <- ggplot(
  ratio_within_selected,
  aes(x = date, y = proportion_within_C0_C2_C3, fill = cluster)
) +
  geom_col(
    position = position_stack(reverse = TRUE),
    width = 0.70, color = "white", linewidth = 0.35
  ) +
  scale_fill_manual(
    values = cluster_colors, breaks = target_clusters,
    drop = FALSE, name = "Cluster"
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1), breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(y = "Composition within clusters 0, 2 and 3") +
  base_theme
save_plot(p_selected, "08_C0_C2_C3_within_selected_composition_by_time")

writeLines(c(
  paste("Input RDS:", input_rds),
  paste("Reference script:", reference_script),
  paste("Cells in input object:", ncol(obj)),
  paste("Cluster field: seurat_clusters"),
  paste("Target cluster order:", paste(target_clusters, collapse = ", ")),
  paste("Time order:", paste(time_levels, collapse = ", ")),
  "Input labels ending in '-root' were displayed without that suffix.",
  "Primary proportion denominator: all cells in the input RDS at each time point.",
  "The normalized-within-selected figure is provided only as a companion view.",
  "Clusters 0 and 2 retain the reference-script colors; cluster 3 uses #4C78A8 to improve separation from cluster 2.",
  "PDF and SVG are vector outputs; PNG previews are 450 dpi."
), file.path(outdir, "09_analysis_notes.txt"))

writeLines(capture.output(sessionInfo()), file.path(outdir, "10_sessionInfo.txt"))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg)) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
  script_destination <- file.path(outdir, "01_run_WT_C0_C2_C3_cell_ratio_by_time.R")
  if (normalizePath(script_path, mustWork = FALSE) !=
      normalizePath(script_destination, mustWork = FALSE)) {
    file.copy(script_path, script_destination, overwrite = TRUE)
  }
}

cat("Completed:", outdir, "\n")
print(ratio_long)
print(selected_totals)
