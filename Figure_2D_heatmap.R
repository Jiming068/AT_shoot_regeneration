# Portable path configuration (override with environment variables).
project_dir <- Sys.getenv("SHOOT_REGEN_PROJECT_DIR", unset = ".")
output_root <- Sys.getenv("SHOOT_REGEN_OUTPUT_DIR", unset = file.path(project_dir, "02_results"))
external_data_dir <- Sys.getenv("SHOOT_REGEN_EXTERNAL_DATA_DIR", unset = file.path(project_dir, "external_data"))
work_dir <- Sys.getenv("SHOOT_REGEN_WORK_DIR", unset = ".")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

input_rds <- file.path(project_dir, "01_obj/01_WT_snRNA/03_snRNA_C0-2-3-subclustered_res0.4_SIM4-SIM16.rds")
gene_file <- file.path(project_dir, "03_script/fig2c_genes_selected_7_11.txt")
outdir <- Sys.getenv(
  "EXPRESSION_OUTDIR",
  file.path(output_root, "145_WT_res04_C02418_updated_fig2c_7_11_expression")
)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
clusters <- c("0", "2", "4", "1", "8")

genes <- read.delim(
  gene_file, header = FALSE, sep = "\t", quote = "", comment.char = "",
  stringsAsFactors = FALSE
)
genes <- genes[nzchar(trimws(genes[[1]])), 1:3]
genes$AGI <- toupper(trimws(genes$AGI))
genes$symbol <- trimws(genes$symbol)
genes$program <- trimws(genes$program)
genes$program <- sub("Root/regeneration competence", "Root competence", genes$program, fixed = TRUE)
genes$program <- sub("oot/regeneration competence", "Root competence", genes$program, fixed = TRUE)
genes$program <- sub("Shoot founder/acquistion", "Shoot founder", genes$program, fixed = TRUE)
genes$program <- sub("Shoot founder/acquisition", "Shoot founder", genes$program, fixed = TRUE)

obj <- readRDS(input_rds)
obj <- subset(obj, subset = subcluster %in% clusters)
obj$subcluster <- factor(as.character(obj$subcluster), levels = clusters)
Idents(obj) <- "subcluster"
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, verbose = FALSE)

present <- genes$AGI[genes$AGI %in% rownames(obj)]
missing <- setdiff(genes$AGI, present)
genes <- genes[match(present, genes$AGI), , drop = FALSE]

avg <- AverageExpression(
  obj, assays = "RNA", features = present, group.by = "subcluster",
  slot = "data", return.seurat = FALSE, verbose = FALSE
)$RNA[present, clusters, drop = FALSE]
dat <- GetAssayData(obj, assay = "RNA", slot = "data")[present, , drop = FALSE]
pct <- sapply(clusters, function(cl) {
  cells <- colnames(obj)[obj$subcluster == cl]
  Matrix::rowMeans(dat[, cells, drop = FALSE] > 0) * 100
})
rownames(pct) <- present
colnames(pct) <- clusters

z <- t(scale(t(avg)))
z[!is.finite(z)] <- 0
z <- pmax(pmin(z, 1), -1)
mm <- t(apply(avg, 1, function(x) {
  if (diff(range(x)) == 0) rep(0, length(x)) else (x - min(x)) / diff(range(x))
}))
colnames(mm) <- clusters

write.csv(data.frame(AGI = present, symbol = genes$symbol, program = genes$program, avg, check.names = FALSE), file.path(outdir, "02_average_log_expression_C02418.csv"), row.names = FALSE)
write.csv(data.frame(AGI = present, symbol = genes$symbol, program = genes$program, pct, check.names = FALSE), file.path(outdir, "03_percent_expressing_C02418.csv"), row.names = FALSE)
write.csv(data.frame(AGI = present, symbol = genes$symbol, program = genes$program, z, check.names = FALSE), file.path(outdir, "04_gene_Zscore_C02418.csv"), row.names = FALSE)
write.csv(data.frame(AGI = present, symbol = genes$symbol, program = genes$program, mm, check.names = FALSE), file.path(outdir, "05_gene_minmax_C02418.csv"), row.names = FALSE)
writeLines(if (length(missing)) missing else "None", file.path(outdir, "06_missing_genes.txt"))

to_long <- function(mat, value = "value") {
  d <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  d$symbol <- factor(genes$symbol[match(d$AGI, genes$AGI)], levels = rev(genes$symbol))
  d$program <- factor(genes$program[match(d$AGI, genes$AGI)], levels = unique(genes$program))
  d$subcluster <- factor(d$subcluster, levels = clusters)
  d
}

theme_pub <- theme_classic(base_family = "Arial", base_size = 11) +
  theme(
    axis.title = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black", face = "italic"),
    strip.text.y = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "#F2F2F2", colour = NA),
    panel.spacing.y = grid::unit(1.0, "mm"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 9),
    plot.margin = margin(3, 4, 3, 3)
  )

save_plot <- function(p, stem, width = 4.15, height = 4.65) {
  ggsave(file.path(outdir, paste0(stem, ".pdf")), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(outdir, paste0(stem, ".png")), p, width = width, height = height, dpi = 600, bg = "white")
}

heatmap_plot <- function(mat, limits, colours, legend_title) {
  ggplot(to_long(mat), aes(subcluster, symbol, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.55) +
    facet_grid(program ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_gradientn(colours = colours, limits = limits, oob = scales::squish, name = legend_title) +
    theme_pub
}

p_z <- heatmap_plot(z, c(-1, 1), c("#A9C9E2", "#F7F7F7", "#E58B79"), "Expression\nZ-score")
p_mm <- heatmap_plot(mm, c(0, 1), c("#F0F0F0", "#F7D7CC", "#E58B79"), "Relative\nexpression")
save_plot(p_z, "07_updated_fig2c_genes_C02418_Zscore_heatmap")
save_plot(p_mm, "08_updated_fig2c_genes_C02418_minmax_heatmap")

avg_long <- to_long(avg, "average_expression")
pct_long <- to_long(pct, "percent_expressing")[, c("AGI", "subcluster", "percent_expressing")]
dot <- merge(avg_long, pct_long, by = c("AGI", "subcluster"), sort = FALSE)
dot$symbol <- factor(genes$symbol[match(dot$AGI, genes$AGI)], levels = rev(genes$symbol))
dot$program <- factor(genes$program[match(dot$AGI, genes$AGI)], levels = unique(genes$program))
dot$subcluster <- factor(dot$subcluster, levels = clusters)
p_dot <- ggplot(dot, aes(subcluster, symbol, size = percent_expressing, colour = average_expression)) +
  geom_point() +
  facet_grid(program ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_gradientn(colours = c("#D9E6F2", "#F4B6A6", "#C43C39"), name = "Average\nexpression") +
  scale_size(range = c(0.3, 5.0), limits = c(0, 100), name = "Expressing\ncells (%)") +
  theme_pub
save_plot(p_dot, "09_updated_fig2c_genes_C02418_dotplot", width = 4.75, height = 4.7)

writeLines(c(
  "Updated expression analysis using the current fig2c_genes_selected_7_11.txt.",
  paste("Subcluster order:", paste(clusters, collapse = " -> ")),
  "RNA log-normalized expression was averaged per subcluster.",
  "Gene-wise Z scores (capped at -2 to 2) and 0-1 min-max relative expression were calculated across the five subclusters.",
  sprintf("Detected genes: %d/%d.", length(present), length(present) + length(missing)),
  "Figures were exported as compact, publication-ready PDF and 600-dpi PNG files."
), file.path(outdir, "10_method_information.txt"))
capture.output(sessionInfo(), file = file.path(outdir, "11_sessionInfo.txt"))
cat("Completed:", outdir, "\n")
