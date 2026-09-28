#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1]], fixed = TRUE)
}

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()

input_file <- get_arg(
  "input",
  file.path(script_dir, "cluster_comparison_heatmap_data_reproduced.csv")
)
output_dir <- get_arg("outdir", script_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cluster_order <- as.character(0:11)
wt_dates <- c("D4-root", "D6-root", "D8-root", "D10-root")
mutant_dates <- c("myb_SIM4", "myb_SIM6", "myb_SIM8", "myb_SIM10")
source_order <- c("WT observed", "myb31 predicted (in-reference)")

comparison <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("cluster", "display_time", "source", "proportion", "label")
if (length(setdiff(required, names(comparison)))) stop("Heatmap source columns are incomplete.")

comparison$cluster <- factor(as.character(comparison$cluster), levels = cluster_order)
comparison$display_time <- factor(
  as.character(comparison$display_time),
  levels = rev(c(wt_dates, mutant_dates))
)
comparison$source <- factor(as.character(comparison$source), levels = source_order)
comparison$label[is.na(comparison$label)] <- ""

max_proportion <- max(comparison$proportion, na.rm = TRUE)

heatmap_plot <- ggplot(comparison, aes(cluster, display_time, fill = proportion)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = label), family = "Arial", size = 3.15, color = "black") +
  scale_x_discrete(limits = cluster_order, drop = FALSE) +
  scale_fill_gradientn(
    colors = c("#D9D9D9", "#FFF7EC", "#FDD49E", "#FC8D59", "#B30000"),
    values = c(0, 0.001, 0.25, 0.60, 1),
    limits = c(0, max_proportion),
    breaks = c(0, 0.2, 0.4, 0.6),
    labels = percent,
    name = "Within-time\nproportion",
    guide = guide_colourbar(nbin = 256, display = "rectangles")
  ) +
  facet_grid(source ~ ., scales = "free_y", space = "free_y") +
  labs(x = "WT cluster", y = NULL) +
  theme_minimal(base_size = 13, base_family = "Arial") +
  theme(
    axis.title.x = element_text(size = 14, margin = margin(t = 4)),
    axis.text.x = element_text(size = 12, face = "bold", color = "#4D4D4D"),
    axis.text.y = element_text(size = 11.5, color = "#4D4D4D"),
    panel.grid = element_blank(),
    strip.text.y = element_text(size = 11.5, angle = 0, face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.margin = margin(4, 5, 3, 4)
  )

pdf_file <- file.path(
  output_dir,
  "05_WT_mutant_cluster_comparison_heatmap_publication_with_values.pdf"
)
png_file <- file.path(
  output_dir,
  "05_WT_mutant_cluster_comparison_heatmap_publication_with_values.png"
)

cairo_pdf(pdf_file, width = 9, height = 5.5, family = "Arial", onefile = FALSE)
print(heatmap_plot)
invisible(dev.off())

png(
  filename = png_file,
  width = 9,
  height = 5.5,
  units = "in",
  res = 600,
  type = "cairo",
  family = "Arial",
  bg = "white"
)
print(heatmap_plot)
invisible(dev.off())

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "05_WT_mutant_cluster_comparison_heatmap_sessionInfo.txt")
)
message("Saved editable PDF: ", pdf_file)
message("Saved 600-dpi PNG: ", png_file)
