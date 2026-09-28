# Portable path configuration (override with environment variables).
project_dir <- Sys.getenv("SHOOT_REGEN_PROJECT_DIR", unset = ".")
output_root <- Sys.getenv("SHOOT_REGEN_OUTPUT_DIR", unset = file.path(project_dir, "02_results"))
external_data_dir <- Sys.getenv("SHOOT_REGEN_EXTERNAL_DATA_DIR", unset = file.path(project_dir, "external_data"))
work_dir <- Sys.getenv("SHOOT_REGEN_WORK_DIR", unset = ".")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(pheatmap)
})

input_rds <- file.path(project_dir, "01_obj/02_WT_stereo/Stereo-seq_SIM6-SIM10_final.rds")
gene_xlsx <- file.path(project_dir, "02_results/01_cell_cycle_related/Arabidopsis_cell_cycle_genes_Kotaro_Torii.xlsx")
# This long-format table was generated from the Excel workbook in the reference
# analysis. SHA-256 validation confirmed that the workbook is byte-identical to
# the user-specified source workbook.
gene_csv <- file.path(project_dir, "02_results/44_WT_subclusters_0_2_7_4_8_1_6_cell_cycle_Kotaro_Torii/00_Kotaro_Torii_cell_cycle_gene_sets_long.csv")
default_outdir <- file.path(project_dir, "02_results/58_WT_stereo_SIM6_four_niches_cell_cycle_Kotaro_Torii")
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1L) args[[1L]] else default_outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
file.copy(gene_xlsx, file.path(outdir, basename(gene_xlsx)), overwrite = TRUE)

source_groups <- c(
  "SIM6_sub-1 core", "SIM6_sub-1 bdry", "SIM6_sub-3 shell", "SIM6_10 shell"
)
display_groups <- c("sub-1 core", "sub-1 bdry", "sub-3 shell", "10 shell")
group_map <- setNames(display_groups, source_groups)
phase_order <- c("G1", "S", "G2/M")
phase_cols <- c("G1" = "#E9B949", "S" = "#DE6E56", "G2/M" = "#6489B9")
group_cols <- c(
  "sub-1 core" = "#B65A4A", "sub-1 bdry" = "#E49A73",
  "sub-3 shell" = "#5F86A6", "10 shell" = "#93B8A0"
)

theme_pub <- function(base_size = 10.5) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.title = element_text(color = "#222222", face = "bold"),
      axis.text = element_text(color = "#222222"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(color = "#222222"),
      plot.title = element_text(face = "bold", size = rel(1.05), hjust = 0),
      plot.margin = margin(5, 7, 5, 5)
    )
}

save_plot <- function(plot, stem, width, height) {
  ggsave(
    file.path(outdir, paste0(stem, ".pdf")), plot,
    width = width, height = height, device = cairo_pdf
  )
  ggsave(
    file.path(outdir, paste0(stem, ".png")), plot,
    width = width, height = height, dpi = 450
  )
}

long <- read.csv(gene_csv, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("Gene", "Phase") %in% colnames(long)))
long$Gene <- toupper(trimws(long$Gene))
long$Phase <- trimws(long$Phase)
long <- long[
  long$Phase %in% phase_order & grepl("^AT[1-5CM]G[0-9]{5}$", long$Gene),
  , drop = FALSE
]
gene_sets <- setNames(
  lapply(phase_order, function(ph) unique(long$Gene[long$Phase == ph])),
  phase_order
)
dup_within <- long[duplicated(long[, c("Gene", "Phase")]), , drop = FALSE]
cross_phase <- aggregate(Phase ~ Gene, long, function(x) paste(unique(x), collapse = ";"))
cross_phase$n_phases <- lengths(strsplit(cross_phase$Phase, ";", fixed = TRUE))
cross_phase <- cross_phase[cross_phase$n_phases > 1, , drop = FALSE]
write.csv(long, file.path(outdir, "00_Kotaro_Torii_cell_cycle_gene_sets_long.csv"), row.names = FALSE)
write.csv(cross_phase, file.path(outdir, "01_cross_phase_gene_conflicts.csv"), row.names = FALSE)

obj <- readRDS(input_rds)
stopifnot("time_niche_class" %in% colnames(obj@meta.data))
missing_groups <- setdiff(source_groups, unique(as.character(obj$time_niche_class)))
if (length(missing_groups)) {
  stop("Missing requested time_niche_class labels: ", paste(missing_groups, collapse = ", "))
}
obj$time_niche_class <- as.character(obj$time_niche_class)
obj <- subset(obj, subset = time_niche_class %in% source_groups)
obj$niche_group <- factor(unname(group_map[obj$time_niche_class]), levels = display_groups)
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

present <- lapply(gene_sets, intersect, y = rownames(obj))
marker_qc <- data.frame(
  phase = phase_order,
  supplied = lengths(gene_sets),
  detected = lengths(present),
  detection_rate = lengths(present) / lengths(gene_sets),
  detected_genes = vapply(present, paste, collapse = ";", FUN.VALUE = character(1))
)
write.csv(marker_qc, file.path(outdir, "02_gene_set_detection_QC.csv"), row.names = FALSE)
stopifnot(all(lengths(present) >= 10), nrow(cross_phase) == 0, nrow(dup_within) == 0)

set.seed(20260825)
safe <- c("G1" = "G1", "S" = "S", "G2/M" = "G2M")
for (ph in phase_order) {
  nm <- paste0("KT_CC_", safe[[ph]])
  obj <- AddModuleScore(
    obj, features = list(present[[ph]]), assay = "RNA",
    name = nm, seed = 20260825
  )
  colnames(obj@meta.data)[colnames(obj@meta.data) == paste0(nm, "1")] <- nm
}
score_cols <- paste0("KT_CC_", unname(safe[phase_order]))
raw_scores <- as.matrix(obj@meta.data[, score_cols, drop = FALSE])
z_scores <- scale(raw_scores)
colnames(z_scores) <- phase_order
obj@meta.data[, paste0("KT_z_", unname(safe[phase_order]))] <- z_scores
ord <- t(apply(z_scores, 1, function(x) order(x, decreasing = TRUE)))
obj$KT_dominant_phase <- factor(phase_order[ord[, 1]], levels = phase_order)
obj$KT_phase_margin <-
  z_scores[cbind(seq_len(nrow(z_scores)), ord[, 1])] -
  z_scores[cbind(seq_len(nrow(z_scores)), ord[, 2])]
obj$KT_phase_confidence <- factor(
  ifelse(obj$KT_phase_margin >= 0.25, "Higher-confidence", "Boundary/low-margin"),
  levels = c("Higher-confidence", "Boundary/low-margin")
)
obj$KT_cycling_score <- rowMeans(z_scores[, c("S", "G2/M"), drop = FALSE]) - z_scores[, "G1"]

md <- obj@meta.data
md$cell <- rownames(md)
cell_cols <- c(
  "cell", "time_niche_class", "niche_group",
  intersect(c("time", "tissue", "x", "y"), colnames(md)),
  "KT_dominant_phase", "KT_phase_margin", "KT_phase_confidence",
  "KT_cycling_score", score_cols,
  paste0("KT_z_", unname(safe[phase_order]))
)
write.csv(md[, cell_cols], file.path(outdir, "03_cell_level_Kotaro_Torii_cycle_scores.csv"), row.names = FALSE)

phase_counts <- as.data.frame(table(niche_group = md$niche_group, phase = md$KT_dominant_phase))
phase_counts$proportion <- ave(
  phase_counts$Freq, phase_counts$niche_group,
  FUN = function(x) x / sum(x)
)
write.csv(phase_counts, file.path(outdir, "04_phase_counts_and_proportions.csv"), row.names = FALSE)

conf_counts <- as.data.frame(table(niche_group = md$niche_group, confidence = md$KT_phase_confidence))
conf_counts$proportion <- ave(
  conf_counts$Freq, conf_counts$niche_group,
  FUN = function(x) x / sum(x)
)
write.csv(conf_counts, file.path(outdir, "05_phase_assignment_confidence.csv"), row.names = FALSE)

mean_scores <- aggregate(md[, score_cols], list(niche_group = md$niche_group), mean)
mean_scores$niche_group <- factor(mean_scores$niche_group, levels = display_groups)
mean_scores <- mean_scores[order(mean_scores$niche_group), ]
write.csv(mean_scores, file.path(outdir, "06_mean_module_scores_by_niche.csv"), row.names = FALSE)

summary_rows <- do.call(rbind, lapply(display_groups, function(gr) {
  d <- md[as.character(md$niche_group) == gr, , drop = FALSE]
  pp <- prop.table(table(factor(d$KT_dominant_phase, levels = phase_order)))
  data.frame(
    niche_group = gr,
    n_cells = nrow(d),
    G1_pct = 100 * pp["G1"],
    S_pct = 100 * pp["S"],
    G2M_pct = 100 * pp["G2/M"],
    cycling_pct = 100 * sum(pp[c("S", "G2/M")]),
    mean_cycling_score = mean(d$KT_cycling_score),
    higher_confidence_pct = 100 * mean(d$KT_phase_confidence == "Higher-confidence"),
    median_phase_margin = median(d$KT_phase_margin),
    dominant_phase = names(which.max(pp))
  )
}))
rownames(summary_rows) <- NULL
write.csv(summary_rows, file.path(outdir, "07_niche_cycle_summary.csv"), row.names = FALSE)

tissue_counts <- as.data.frame(table(niche_group = md$niche_group, tissue = md$tissue))
tissue_counts$proportion <- ave(
  tissue_counts$Freq, tissue_counts$niche_group,
  FUN = function(x) x / sum(x)
)
write.csv(tissue_counts, file.path(outdir, "08_tissue_composition_by_niche.csv"), row.names = FALSE)
phase_tissue <- as.data.frame(table(
  niche_group = md$niche_group, tissue = md$tissue, phase = md$KT_dominant_phase
))
phase_tissue$proportion <- ave(
  phase_tissue$Freq,
  interaction(phase_tissue$niche_group, phase_tissue$tissue),
  FUN = function(x) if (sum(x) == 0) 0 else x / sum(x)
)
write.csv(phase_tissue, file.path(outdir, "09_phase_by_niche_and_tissue.csv"), row.names = FALSE)

p_phase <- ggplot(phase_counts, aes(niche_group, proportion, fill = phase)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = phase_cols, drop = FALSE) +
  scale_y_continuous(
    labels = function(x) paste0(round(100 * x), "%"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = NULL, y = "Fraction of spatial cells", fill = "Dominant module") +
  theme_pub(10.5) +
  theme(axis.text.x = element_text(face = "bold", angle = 25, hjust = 1))
save_plot(p_phase, "10_Kotaro_Torii_phase_composition", 4.45, 3.35)

mat <- as.matrix(mean_scores[, score_cols])
colnames(mat) <- phase_order
mat_z <- scale(mat)
for (ext in c("pdf", "png")) {
  path <- file.path(outdir, paste0("11_Kotaro_Torii_module_heatmap.", ext))
  if (ext == "pdf") {
    cairo_pdf(path, width = 3.55, height = 3.35, family = "Arial")
  } else {
    png(path, width = 1598, height = 1508, res = 450)
  }
  pheatmap(
    mat_z, cluster_rows = FALSE, cluster_cols = FALSE,
    color = colorRampPalette(c("#BFD7EA", "#F7F7F7", "#EFA7A0"))(101),
    breaks = seq(-2, 2, length.out = 102), border_color = "white",
    fontsize = 10.5, fontsize_row = 10.5, fontsize_col = 10.5,
    angle_col = 0
  )
  dev.off()
}

expr <- GetAssayData(obj, assay = "RNA", slot = "data")
rep_genes <- unlist(lapply(phase_order, function(ph) {
  gs <- present[[ph]]
  pct <- Matrix::rowMeans(expr[gs, , drop = FALSE] > 0)
  names(sort(pct, decreasing = TRUE))[seq_len(min(6, length(gs)))]
}), use.names = FALSE)
rep_table <- long[match(rep_genes, long$Gene), c("Gene", "Phase")]
write.csv(rep_table, file.path(outdir, "12_representative_genes_for_dotplot.csv"), row.names = FALSE)

p_dot <- DotPlot(obj, features = rep_genes, group.by = "niche_group", assay = "RNA") +
  scale_color_gradient2(low = "#BFD7EA", mid = "#F7F7F7", high = "#D95F59", midpoint = 0) +
  labs(
    x = "Representative cell-cycle gene (AGI)", y = NULL,
    color = "Average\nexpression", size = "% expressed"
  ) +
  theme_pub(9.2) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7.8),
    axis.text.y = element_text(face = "bold")
  )
save_plot(p_dot, "13_Kotaro_Torii_representative_gene_dotplot", 7.7, 3.55)

p_conf <- ggplot(conf_counts, aes(niche_group, proportion, fill = confidence)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = c(
    "Higher-confidence" = "#507C8C", "Boundary/low-margin" = "#D8D8D8"
  )) +
  scale_y_continuous(
    labels = function(x) paste0(round(100 * x), "%"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = NULL, y = "Fraction of spatial cells", fill = NULL) +
  theme_pub(10.5) +
  theme(axis.text.x = element_text(face = "bold", angle = 25, hjust = 1))
save_plot(p_conf, "14_phase_assignment_confidence", 4.6, 3.3)

# Use tissue-centered local coordinates so each spatial replicate can be faceted.
md$x_local <- ave(md$x, md$tissue, FUN = function(v) v - min(v))
md$y_local <- ave(md$y, md$tissue, FUN = function(v) v - min(v))
p_spatial <- ggplot(md, aes(x_local, y_local, color = KT_dominant_phase)) +
  geom_point(size = 0.55, alpha = 0.9) +
  facet_wrap(~tissue, nrow = 2) +
  scale_color_manual(values = phase_cols, drop = FALSE) +
  coord_equal() +
  labs(x = NULL, y = NULL, color = "Dominant module") +
  theme_void(base_size = 9.5, base_family = "Arial") +
  theme(
    strip.text = element_text(face = "bold", color = "#222222"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.spacing = grid::unit(0.7, "lines"),
    plot.margin = margin(4, 5, 4, 5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA)
  )
save_plot(p_spatial, "15_Kotaro_Torii_phase_spatial_map", 5.0, 4.5)

saveRDS(obj, file.path(outdir, "16_SIM6_four_niches_Kotaro_Torii_cycle_scores.rds"), compress = TRUE)
writeLines(capture.output(sessionInfo()), file.path(outdir, "17_sessionInfo.txt"))
writeLines(c(
  "Interpretation guardrails:",
  "1. Dominant module denotes relative transcriptomic resemblance, not a direct measurement of cell-cycle duration or division rate.",
  "2. G2 and M are combined because the Kotaro Torii workbook defines a single G2/M marker set.",
  "3. G1 is not relabeled as G0 because the source gene set does not define a quiescent program.",
  "4. Spatial cells are not independent biological replicates; tissue-stratified tables are supplied for diagnostic review.",
  "5. AddModuleScore values depend on detected expression and control-gene matching and should be interpreted comparatively among these four selected SIM6 niches."
), file.path(outdir, "18_interpretation_guardrails.txt"))

# Non-parametric descriptive comparisons across spatial groups. P values are
# cell-level and are therefore exploratory rather than replicate-aware inference.
kw_rows <- do.call(rbind, lapply(c(score_cols, "KT_cycling_score"), function(v) {
  kt <- kruskal.test(md[[v]] ~ md$niche_group)
  data.frame(feature = v, statistic = unname(kt$statistic), df = unname(kt$parameter), p_value = kt$p.value)
}))
kw_rows$p_adj_BH <- p.adjust(kw_rows$p_value, method = "BH")
write.csv(kw_rows, file.path(outdir, "19_global_score_comparisons_exploratory.csv"), row.names = FALSE)

pair_rows <- do.call(rbind, lapply(c(score_cols, "KT_cycling_score"), function(v) {
  pw <- pairwise.wilcox.test(md[[v]], md$niche_group, p.adjust.method = "BH", exact = FALSE)
  z <- as.data.frame(as.table(pw$p.value), stringsAsFactors = FALSE)
  z <- z[!is.na(z$Freq), , drop = FALSE]
  data.frame(feature = v, group_1 = z$Var1, group_2 = z$Var2, p_adj_BH = z$Freq)
}))
write.csv(pair_rows, file.path(outdir, "20_pairwise_score_comparisons_exploratory.csv"), row.names = FALSE)

score_long <- data.frame(
  niche_group = rep(md$niche_group, times = length(phase_order)),
  module = factor(rep(phase_order, each = nrow(md)), levels = phase_order),
  z_score = as.vector(z_scores[, phase_order])
)
p_violin <- ggplot(score_long, aes(niche_group, z_score, fill = niche_group)) +
  geom_violin(scale = "width", trim = TRUE, color = NA, alpha = 0.78) +
  geom_boxplot(width = 0.13, outlier.shape = NA, fill = "white", color = "#333333", linewidth = 0.3) +
  facet_wrap(~module, nrow = 1) +
  scale_fill_manual(values = group_cols, guide = "none") +
  labs(x = NULL, y = "Standardized module score") +
  theme_pub(9.5) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
    strip.text = element_text(face = "bold")
  )
save_plot(p_violin, "21_Kotaro_Torii_module_score_distributions", 7.4, 3.25)

methods_results <- c(
  "METHODS",
  "Spatial cells assigned to SIM6_sub-1 core, SIM6_sub-1 bdry, SIM6_sub-3 shell, or SIM6_10 shell in the time_niche_class metadata field were retained. RNA counts were log-normalized with a scale factor of 10,000. G1, S, and G2/M gene sets were taken from Arabidopsis_cell_cycle_genes_Kotaro_Torii.xlsx. Detected genes were scored independently with Seurat AddModuleScore using matched control genes and a fixed random seed. The three module scores were standardized across the selected spatial cells, and the phase with the largest standardized score was reported as the dominant transcriptomic module. Assignment margin was defined as the difference between the largest and second-largest standardized scores; margins >=0.25 were designated higher confidence. The cycling index was calculated as mean(zS, zG2/M) - zG1.",
  "",
  "INTERPRETATION",
  "The analysis compares relative cell-cycle transcriptional programs among the four selected SIM6 spatial niches. Dominant-module labels are not equivalent to direct measurements of cell-cycle duration, mitotic rate, or quiescence. Because individual spatial cells are nested within tissues, cell-level non-parametric P values are supplied only as exploratory summaries; biological inference should emphasize effect patterns, phase composition, and consistency across tissues."
)
writeLines(methods_results, file.path(outdir, "22_English_Methods_and_interpretation.txt"))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg)) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
  file.copy(
    script_path,
    file.path(outdir, "23_run_stereo_D6_four_niches_cell_cycle_kotaro_torii.R"),
    overwrite = TRUE
  )
}

cat("Completed:", outdir, "\n")
print(marker_qc[, 1:4])
print(summary_rows)
