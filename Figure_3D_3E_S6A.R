#!/usr/bin/env Rscript

# Portable path configuration (override with environment variables).
project_dir <- Sys.getenv("SHOOT_REGEN_PROJECT_DIR", unset = ".")
output_root <- Sys.getenv("SHOOT_REGEN_OUTPUT_DIR", unset = file.path(project_dir, "02_results"))
external_data_dir <- Sys.getenv("SHOOT_REGEN_EXTERNAL_DATA_DIR", unset = file.path(project_dir, "external_data"))
work_dir <- Sys.getenv("SHOOT_REGEN_WORK_DIR", unset = ".")

suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(scales) })
set.seed(128)

outdir <- file.path(output_root, "147_D6_spatial_cluster2_top200_module_score")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
marker_file <- file.path(project_dir, "02_results/146_Tangram_WT_C02418_to_spatial_5timecelltype/03_snRNA_top200_markers_per_cluster.csv")
spatial_file <- file.path(project_dir, "01_obj/02_WT_stereo/Stereo-seq_SIM6-SIM10_final.rds")

markers <- read.csv(marker_file, stringsAsFactors = FALSE)
genes <- unique(toupper(markers$gene[as.character(markers$cluster) == "2"]))
if (length(genes) != 200) warning("Cluster 2 marker list contains ", length(genes), " unique genes rather than 200.")

sp <- readRDS(spatial_file)
DefaultAssay(sp) <- "RNA"
sp <- subset(sp, subset = time == "SIM6")
sp <- NormalizeData(sp, assay = "RNA", verbose = FALSE)
genes_present <- intersect(genes, toupper(rownames(sp)))
feature_lookup <- setNames(rownames(sp), toupper(rownames(sp)))
genes_present_features <- unname(feature_lookup[genes_present])
write.csv(
  data.frame(gene = genes_present_features),
  file.path(outdir, "02_cluster2_top200_module_genes_present_in_spatial.csv"),
  row.names = FALSE
)
if (length(genes_present_features) < 20) stop("Fewer than 20 cluster 2 module genes are present in the spatial RNA assay.")

sp <- AddModuleScore(
  sp, features = list(genes_present_features), assay = "RNA",
  name = "C2_top200_module", ctrl = 100, nbin = 24, seed = 128
)
sp$C2_top200_module_score <- sp@meta.data[["C2_top200_module1"]]
md <- sp@meta.data
md$cell_id <- rownames(md)

# Publication map: scores below the D6 median are shown as neutral grey;
# the upper tail is capped at the 99th percentile.
lo <- unname(quantile(md$C2_top200_module_score, 0.50, na.rm = TRUE))
hi <- unname(quantile(md$C2_top200_module_score, 0.99, na.rm = TRUE))
md$score_plot <- pmin(pmax(md$C2_top200_module_score, lo), hi)
md$tissue <- factor(as.character(md$tissue), levels = unique(as.character(md$tissue)))
write.csv(
  md[, c("cell_id", "tissue", "time", "time_celltype", "cluster", "x", "y", "C2_top200_module_score", "score_plot")],
  gzfile(file.path(outdir, "03_D6_spatial_cluster2_top200_module_scores.csv.gz")), row.names = FALSE
)

theme_map <- theme_void(base_family = "Arial") +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text = element_text(size = 10, face = "bold", colour = "black"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    plot.margin = margin(3, 3, 3, 3)
  )

map <- ggplot(md, aes(x = x, y = y, colour = score_plot)) +
  geom_point(size = 0.38, alpha = 0.95) +
  scale_y_reverse() +
  scale_colour_gradientn(
    colours = c("#D9D9D9", "#F9E1D8", "#F5B8A2", "#E6836F", "#C94B4A"),
    limits = c(lo, hi), oob = squish, name = "Cluster 2 top200\nmodule score"
  ) +
  facet_wrap(~tissue, ncol = 4, scales = "free") + theme_map
ggsave(file.path(outdir, "04_D6_spatial_cluster2_top200_module_score_in_situ.pdf"), map, width = 8.0, height = 4.2, device = cairo_pdf)
ggsave(file.path(outdir, "04_D6_spatial_cluster2_top200_module_score_in_situ.png"), map, width = 8.0, height = 4.2, dpi = 600, bg = "white")

selected <- c("SIM6_10", "SIM6_sub-1", "SIM6_sub-3")
vdat <- md[md$time_celltype %in% selected, , drop = FALSE]
vdat$time_celltype <- factor(vdat$time_celltype, levels = selected)
summary_table <- do.call(rbind, lapply(selected, function(group) {
  values <- vdat$C2_top200_module_score[vdat$time_celltype == group]
  data.frame(time_celltype = group, n = length(values), mean = mean(values), median = median(values))
}))
write.csv(summary_table, file.path(outdir, "05_selected_D6_time_celltype_cluster2_top200_module_score_summary.csv"), row.names = FALSE)

# Sample-aware comparison: mean score per tissue and spatial class, followed by paired Wilcoxon tests.
pb <- aggregate(C2_top200_module_score ~ tissue + time_celltype, data = vdat, FUN = mean)
wide <- reshape(pb, idvar = "tissue", timevar = "time_celltype", direction = "wide")
score_names <- paste0("C2_top200_module_score.", selected)
for (nm in setdiff(score_names, colnames(wide))) wide[[nm]] <- NA_real_
complete <- wide[complete.cases(wide[, score_names, drop = FALSE]), , drop = FALSE]
comparisons <- data.frame(group1 = selected[c(1, 1, 2)], group2 = selected[c(2, 3, 3)], stringsAsFactors = FALSE)
stats <- do.call(rbind, lapply(seq_len(nrow(comparisons)), function(i) {
  a <- complete[[paste0("C2_top200_module_score.", comparisons$group1[i])]]
  b <- complete[[paste0("C2_top200_module_score.", comparisons$group2[i])]]
  p <- if (length(a) >= 2) suppressWarnings(wilcox.test(a, b, paired = TRUE, exact = FALSE)$p.value) else NA_real_
  data.frame(
    group1 = comparisons$group1[i], group2 = comparisons$group2[i],
    n_paired_samples = length(a), median_difference = if (length(a)) median(b - a) else NA_real_,
    p_value = p
  )
}))
stats$FDR_BH <- p.adjust(stats$p_value, method = "BH")
stats$label <- ifelse(
  is.na(stats$FDR_BH), "NA",
  ifelse(stats$FDR_BH < 0.001, "FDR < 0.001", paste0("FDR = ", formatC(stats$FDR_BH, format = "f", digits = 3)))
)
write.csv(pb, file.path(outdir, "06_D6_selected_groups_sample_mean_module_score.csv"), row.names = FALSE)
write.csv(stats, file.path(outdir, "07_D6_selected_groups_sample_aware_paired_Wilcoxon.csv"), row.names = FALSE)

yrange <- diff(range(vdat$C2_top200_module_score, na.rm = TRUE))
step <- ifelse(yrange > 0, yrange * 0.11, 0.05)
top <- max(vdat$C2_top200_module_score, na.rm = TRUE)
stats$y <- top + step * seq_len(nrow(stats))
stats$x1 <- match(stats$group1, selected)
stats$x2 <- match(stats$group2, selected)

violin <- ggplot(vdat, aes(x = time_celltype, y = C2_top200_module_score, fill = time_celltype)) +
  geom_violin(scale = "width", trim = TRUE, colour = NA, alpha = 0.90) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", colour = "#4A4A4A", linewidth = 0.32) +
  scale_fill_manual(values = c("SIM6_10" = "#D76E66", "SIM6_sub-1" = "#E69473", "SIM6_sub-3" = "#F0B67F"), guide = "none") +
  geom_segment(data = stats, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.35) +
  geom_segment(data = stats, aes(x = x1, xend = x1, y = y - step * 0.12, yend = y), inherit.aes = FALSE, linewidth = 0.35) +
  geom_segment(data = stats, aes(x = x2, xend = x2, y = y - step * 0.12, yend = y), inherit.aes = FALSE, linewidth = 0.35) +
  geom_text(data = stats, aes(x = (x1 + x2) / 2, y = y + step * 0.07, label = label), inherit.aes = FALSE, size = 3.1, fontface = "bold") +
  coord_cartesian(ylim = c(min(vdat$C2_top200_module_score, na.rm = TRUE), top + step * 3.7), clip = "off") +
  labs(x = "SIM6 spatial cell type", y = "Cluster 2 top200 module score") +
  theme_classic(base_family = "Arial", base_size = 11) +
  theme(axis.text.x = element_text(size = 10), axis.title = element_text(face = "bold"), plot.margin = margin(7, 5, 4, 4))
ggsave(file.path(outdir, "08_D6_time_celltype_cluster2_top200_module_score_violin.pdf"), violin, width = 3.9, height = 4.1, device = cairo_pdf)
ggsave(file.path(outdir, "08_D6_time_celltype_cluster2_top200_module_score_violin.png"), violin, width = 3.9, height = 4.1, dpi = 600, bg = "white")

writeLines(c(
  "Cluster 2 top200 marker module scored in all SIM6 spatial cells.",
  paste("Requested marker genes:", length(genes)),
  paste("Genes present in spatial RNA assay:", length(genes_present_features)),
  paste("SIM6 spatial cells:", nrow(md)),
  paste("Display lower limit: SIM6 score 50th percentile =", signif(lo, 5)),
  paste("Display upper limit: SIM6 score 99th percentile =", signif(hi, 5)),
  "Module score: Seurat AddModuleScore, ctrl=100, nbin=24, seed=128.",
  "Cells below the display median are rendered as neutral grey.",
  "Violin-plot significance is based on paired tissue-level mean scores with BH correction."
), file.path(outdir, "00_module_score_information.txt"))
capture.output(sessionInfo(), file = file.path(outdir, "09_sessionInfo.txt"))
cat("Completed:", outdir, "\n")
