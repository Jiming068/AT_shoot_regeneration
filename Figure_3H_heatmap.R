#!/usr/bin/env Rscript

# Portable path configuration (override with environment variables).
project_dir <- Sys.getenv("SHOOT_REGEN_PROJECT_DIR", unset = ".")
output_root <- Sys.getenv("SHOOT_REGEN_OUTPUT_DIR", unset = file.path(project_dir, "02_results"))
external_data_dir <- Sys.getenv("SHOOT_REGEN_EXTERNAL_DATA_DIR", unset = file.path(project_dir, "external_data"))
work_dir <- Sys.getenv("SHOOT_REGEN_WORK_DIR", unset = ".")

suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
set.seed(128)

input_rds <- file.path(project_dir, "01_obj/02_WT_stereo/Stereo-seq_SIM6-SIM10_final.rds")
gene_file <- file.path(project_dir, "03_script/fig2c_genes_selected_7_17.txt")
outdir <- file.path(output_root, "148_spatial_fig2c_7_17_expression_six_time_niche_classes")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
groups <- c("SIM6_sub-1 core", "SIM6_sub-1 bdry", "SIM6_sub-3 shell", "SIM6_10 shell", "SIM8_sub-1", "SIM10_sub-1")

genes <- read.delim(gene_file, header = FALSE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)
genes <- genes[nzchar(trimws(genes[[1]])), 1:3]
genes$AGI <- toupper(trimws(genes$AGI))
genes$symbol <- trimws(genes$symbol)
genes$program <- trimws(genes$program)
genes$program <- sub("Shoot founder/acquistion", "Shoot founder", genes$program, fixed = TRUE)
genes$program <- sub("Shoot founder/acquisition", "Shoot founder", genes$program, fixed = TRUE)
genes$program <- sub("Root/regeneration competence", "Root competence", genes$program, fixed = TRUE)

obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"
obj <- subset(obj, subset = time_niche_class %in% groups)
obj$time_niche_class <- factor(as.character(obj$time_niche_class), levels = groups)
Idents(obj) <- "time_niche_class"
obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)

feature_lookup <- setNames(rownames(obj), toupper(rownames(obj)))
present_agi <- genes$AGI[genes$AGI %in% names(feature_lookup)]
missing <- setdiff(genes$AGI, present_agi)
genes <- genes[match(present_agi, genes$AGI), , drop = FALSE]
features <- unname(feature_lookup[present_agi])

avg <- AverageExpression(
  obj, assays = "RNA", features = features, group.by = "time_niche_class",
  slot = "data", return.seurat = FALSE, verbose = FALSE
)$RNA[features, groups, drop = FALSE]
rownames(avg) <- present_agi
expr <- GetAssayData(obj, assay = "RNA", slot = "data")[features, , drop = FALSE]
pct <- sapply(groups, function(group) {
  cells <- colnames(obj)[obj$time_niche_class == group]
  Matrix::rowMeans(expr[, cells, drop = FALSE] > 0) * 100
})
rownames(pct) <- present_agi
colnames(pct) <- groups

z_raw <- t(scale(t(avg)))
z_raw[!is.finite(z_raw)] <- 0
z_plot <- pmax(pmin(z_raw, 1), -1)
minmax <- t(apply(avg, 1, function(x) {
  if (diff(range(x)) == 0) rep(0, length(x)) else (x - min(x)) / diff(range(x))
}))
colnames(minmax) <- groups

write.csv(data.frame(AGI = present_agi, symbol = genes$symbol, program = genes$program, avg, check.names = FALSE), file.path(outdir, "02_average_log_expression.csv"), row.names = FALSE)
write.csv(data.frame(AGI = present_agi, symbol = genes$symbol, program = genes$program, pct, check.names = FALSE), file.path(outdir, "03_percent_expressing.csv"), row.names = FALSE)
write.csv(data.frame(AGI = present_agi, symbol = genes$symbol, program = genes$program, z_raw, check.names = FALSE), file.path(outdir, "04_gene_Zscore_unclipped.csv"), row.names = FALSE)
write.csv(data.frame(AGI = present_agi, symbol = genes$symbol, program = genes$program, minmax, check.names = FALSE), file.path(outdir, "05_gene_minmax.csv"), row.names = FALSE)
writeLines(if (length(missing)) missing else "None", file.path(outdir, "06_missing_genes.txt"))
write.csv(as.data.frame(table(obj$time_niche_class)), file.path(outdir, "07_cell_counts_by_time_niche_class.csv"), row.names = FALSE)

to_long <- function(mat, value_name = "value") {
  d <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  d$symbol <- factor(genes$symbol[match(d$AGI, genes$AGI)], levels = rev(genes$symbol))
  d$program <- factor(genes$program[match(d$AGI, genes$AGI)], levels = unique(genes$program))
  d$time_niche_class <- factor(d$time_niche_class, levels = groups)
  d
}

theme_pub <- theme_classic(base_family = "Arial", base_size = 10) +
  theme(
    axis.title = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, size = 9.5, colour = "black"),
    axis.text.y = element_text(size = 9.5, colour = "black", face = "italic"),
    strip.text.y = element_text(size = 8.5, face = "bold"),
    strip.background = element_rect(fill = "#F2F2F2", colour = NA),
    panel.spacing.y = grid::unit(0.9, "mm"),
    legend.title = element_text(size = 8.5, face = "bold"),
    legend.text = element_text(size = 8),
    plot.margin = margin(3, 4, 3, 3)
  )

heatmap_plot <- function(mat, limits, colours, legend_title) {
  ggplot(to_long(mat), aes(time_niche_class, symbol, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    facet_grid(program ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_gradientn(colours = colours, limits = limits, oob = scales::squish, name = legend_title) +
    theme_pub
}

p_z <- heatmap_plot(z_plot, c(-1, 1), c("#A9C9E2", "#F7F7F7", "#E58B79"), "Expression\nZ-score")
p_mm <- heatmap_plot(minmax, c(0, 1), c("#F0F0F0", "#F7D7CC", "#E58B79"), "Relative\nexpression")
ggsave(file.path(outdir, "08_fig2c_7_17_spatial_six_classes_Zscore_heatmap.pdf"), p_z, width = 5.2, height = 7.7, device = cairo_pdf)
ggsave(file.path(outdir, "08_fig2c_7_17_spatial_six_classes_Zscore_heatmap.png"), p_z, width = 5.2, height = 7.7, dpi = 600, bg = "white")
ggsave(file.path(outdir, "09_fig2c_7_17_spatial_six_classes_minmax_heatmap.pdf"), p_mm, width = 5.2, height = 7.7, device = cairo_pdf)
ggsave(file.path(outdir, "09_fig2c_7_17_spatial_six_classes_minmax_heatmap.png"), p_mm, width = 5.2, height = 7.7, dpi = 600, bg = "white")

writeLines(c(
  "Spatial expression heatmaps for the current fig2c_genes_selected_7_17 gene set.",
  paste("Group order:", paste(groups, collapse = " -> ")),
  paste("Detected genes:", length(present_agi), "of", length(present_agi) + length(missing)),
  "RNA log-normalized expression was averaged within each time_niche_class.",
  "Z-score visualization was clipped to [-1, 1]; the complete unclipped matrix is retained in file 04.",
  "A second heatmap shows gene-wise 0-1 min-max relative expression."
), file.path(outdir, "10_method_information.txt"))
capture.output(sessionInfo(), file = file.path(outdir, "11_sessionInfo.txt"))
cat("Completed:", outdir, "\n")
