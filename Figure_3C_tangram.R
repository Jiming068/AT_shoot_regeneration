#!/usr/bin/env Rscript

# Reproduce the publication-ready Tangram row-Z-score heatmap.
#
# Input:
#   09_tangram_cluster_vs_time_celltype_raw_score.csv
#     Rows are snRNA-seq clusters and columns are spatial time/cell-type states.
#
# Calculation:
#   For each snRNA-seq cluster, the mean Tangram score across spatial states is
#   standardized using the population standard deviation (ddof = 0), matching
#   the original Python workflow exactly.
#
# Outputs:
#   20_tangram_row_zscore_heatmap_publication.pdf  (editable vector artwork)
#   20_tangram_row_zscore_heatmap_publication.png  (600 dpi raster preview)
#   20_tangram_row_zscore_heatmap_publication_source_data.csv
#   20_tangram_row_zscore_heatmap_publication_sessionInfo.txt

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

parse_option <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(prefix, "", hit[[1]], fixed = TRUE)
}

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()

input_file <- parse_option(
  "input",
  file.path(script_dir, "09_tangram_cluster_vs_time_celltype_raw_score.csv")
)
output_dir <- parse_option("outdir", script_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cluster_order <- c("0", "2", "4", "1", "8")
state_order <- c("SIM6_10", "SIM6_sub-1", "SIM6_sub-3", "SIM8_sub-1", "SIM10_sub-1")
color_limit <- 2

raw <- read.csv(input_file, row.names = 1, check.names = FALSE)
raw <- as.matrix(raw[cluster_order, state_order, drop = FALSE])
storage.mode(raw) <- "double"

population_zscore <- function(x) {
  center <- mean(x, na.rm = TRUE)
  sigma <- sqrt(mean((x - center)^2, na.rm = TRUE))
  if (!is.finite(sigma) || sigma == 0) return(rep(0, length(x)))
  (x - center) / sigma
}

zscore <- t(apply(raw, 1, population_zscore))
dimnames(zscore) <- list(cluster_order, state_order)

source_file <- file.path(
  output_dir,
  "20_tangram_row_zscore_heatmap_publication_source_data.csv"
)
write.csv(zscore, source_file, quote = FALSE)

plot_data <- as.data.frame(zscore) |>
  mutate(snRNA_cluster = rownames(zscore)) |>
  pivot_longer(
    cols = all_of(state_order),
    names_to = "spatial_state",
    values_to = "row_zscore"
  ) |>
  mutate(
    snRNA_cluster = factor(snRNA_cluster, levels = rev(cluster_order)),
    spatial_state = factor(spatial_state, levels = state_order),
    label_color = if_else(abs(row_zscore) >= 0.95, "white", "#222222")
  )

heatmap_plot <- ggplot(plot_data, aes(spatial_state, snRNA_cluster, fill = row_zscore)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(
    aes(label = sprintf("%.2f", row_zscore), color = label_color),
    family = "Arial",
    size = 3.25,
    show.legend = FALSE
  ) +
  scale_color_identity() +
  scale_fill_gradient2(
    low = "#4A90C2",
    mid = "#F7F7F7",
    high = "#C63D3B",
    midpoint = 0,
    limits = c(-color_limit, color_limit),
    breaks = c(-2, -1, 0, 1, 2),
    oob = squish,
    name = "Row Z-score"
  ) +
  coord_fixed(clip = "off") +
  labs(x = "Spatial state", y = "snRNA-seq cluster") +
  theme_classic(base_size = 9, base_family = "Arial") +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(face = "bold", margin = margin(t = 5)),
    axis.title.y = element_text(face = "bold", margin = margin(r = 5)),
    axis.text.x = element_text(angle = 42, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.key.height = grid::unit(15, "mm"),
    plot.margin = margin(4, 6, 4, 4)
  )

pdf_file <- file.path(output_dir, "20_tangram_row_zscore_heatmap_publication.pdf")
png_file <- file.path(output_dir, "20_tangram_row_zscore_heatmap_publication.png")

ggsave(
  pdf_file,
  plot = heatmap_plot,
  device = cairo_pdf,
  width = 4.65,
  height = 3.55,
  units = "in",
  bg = "white"
)
png(
  filename = png_file,
  width = 4.65,
  height = 3.55,
  units = "in",
  res = 600,
  type = "cairo",
  family = "Arial",
  bg = "white"
)
print(heatmap_plot)
invisible(dev.off())

session_file <- file.path(
  output_dir,
  "20_tangram_row_zscore_heatmap_publication_sessionInfo.txt"
)
capture.output(sessionInfo(), file = session_file)

message("Saved publication figure: ", pdf_file)
message("Saved 600-dpi preview: ", png_file)
message("Saved source matrix: ", source_file)
