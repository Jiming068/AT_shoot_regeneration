#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1]], fixed = TRUE)
}

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()

metadata_csv <- get_arg("metadata", file.path(script_dir, "00_spatial_metadata.csv.gz"))
classification_csv <- get_arg(
  "classification",
  file.path(script_dir, "SIM6_sub1_sub3_10_niche_classification_reproduced.csv")
)
output_dir <- get_arg("outdir", script_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

class_order <- c(
  "sub-1 core", "sub-1 bdry", "sub-3 shell",
  "sub-3 disp.", "10 shell", "10 disp."
)
class_colors <- c(
  "sub-1 core" = "#1B9E77",
  "sub-1 bdry" = "#66C2A5",
  "sub-3 shell" = "#D95F02",
  "sub-3 disp." = "#FDAE61",
  "10 shell" = "#7570B3",
  "10 disp." = "#BCBDDC"
)
background_color <- "#E4E4E4"
tissue_order <- c("SIM6_1-1", "SIM6_1-2", "SIM6_1-3", "SIM6_1-4", "SIM6_2-1", "SIM6_2-2", "SIM6_3-1", "SIM6_3-2")
padding_ratio <- 0.15

plot_data <- read.csv(
  metadata_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = c(cell_id = "character", time = "character", cluster = "character", tissue = "character")
)
required_meta <- c("cell_id", "time", "tissue", "x", "y")
missing_meta <- setdiff(required_meta, colnames(plot_data))
if (length(missing_meta)) stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "))
if (anyDuplicated(plot_data$cell_id)) stop("Metadata table contains duplicated cell IDs.")

plot_data <- plot_data |>
  filter(time == "SIM6", tissue %in% tissue_order)

classification <- read.csv(classification_csv, stringsAsFactors = FALSE, check.names = FALSE)
if (anyDuplicated(classification$cell_id)) stop("Classification table contains duplicated cell IDs.")
plot_data <- left_join(plot_data, classification[, c("cell_id", "niche_class")], by = "cell_id")

ranges <- plot_data |>
  group_by(tissue) |>
  summarise(x_range = max(x) - min(x), y_range = max(y) - min(y), .groups = "drop")
max_x_range <- max(ranges$x_range)
max_y_range <- max(ranges$y_range)
x_step <- max_x_range * (1 + padding_ratio)
y_step <- max_y_range * (1 + padding_ratio)

layout <- data.frame(
  tissue = tissue_order,
  column = rep(0:3, 2),
  row = rep(0:1, each = 4),
  stringsAsFactors = FALSE
)
plot_data <- plot_data |>
  group_by(tissue) |>
  mutate(local_x = x - min(x), local_y = y - min(y)) |>
  ungroup() |>
  left_join(layout, by = "tissue") |>
  mutate(
    plot_x = local_x + column * x_step,
    plot_y = local_y + row * y_step,
    niche_class = factor(niche_class, levels = class_order)
  )

background <- filter(plot_data, is.na(niche_class))
highlight <- filter(plot_data, !is.na(niche_class))

spatial_plot <- ggplot() +
  geom_point(
    data = background,
    aes(plot_x, plot_y),
    color = background_color,
    size = 0.38,
    shape = 16,
    stroke = 0
  ) +
  geom_point(
    data = highlight,
    aes(plot_x, plot_y, color = niche_class),
    size = 0.58,
    shape = 16,
    stroke = 0
  ) +
  scale_color_manual(values = class_colors, limits = class_order, drop = FALSE) +
  scale_y_reverse() +
  coord_fixed(expand = FALSE, clip = "off") +
  guides(color = guide_legend(
    title = "Core-shell class",
    override.aes = list(size = 2.2),
    title.position = "top"
  )) +
  theme_void(base_family = "Arial", base_size = 9) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8.5),
    legend.key.height = grid::unit(4.2, "mm"),
    legend.margin = margin(0, 0, 0, 5),
    plot.margin = margin(3, 3, 3, 3)
  )

pdf_file <- file.path(output_dir, "SIM6_niche_reclassification_spatial_publication.pdf")
png_file <- file.path(output_dir, "SIM6_niche_reclassification_spatial_publication.png")
source_file <- file.path(output_dir, "SIM6_niche_reclassification_spatial_plot_source_data.csv.gz")

ggsave(
  pdf_file,
  plot = spatial_plot,
  device = cairo_pdf,
  width = 13,
  height = 6.8,
  units = "in",
  bg = "white"
)
png(
  filename = png_file,
  width = 13,
  height = 6.8,
  units = "in",
  res = 600,
  type = "cairo",
  family = "Arial",
  bg = "white"
)
print(spatial_plot)
invisible(dev.off())

source_data <- plot_data[, c("cell_id", "tissue", "plot_x", "plot_y", "niche_class")]
con <- gzfile(source_file, open = "wt")
write.csv(source_data, con, row.names = FALSE, quote = TRUE)
close(con)

capture.output(sessionInfo(), file = file.path(output_dir, "SIM6_niche_reclassification_spatial_sessionInfo.txt"))
message("Saved publication PDF: ", pdf_file)
message("Saved 600-dpi PNG: ", png_file)
message("Saved plot source data: ", source_file)
