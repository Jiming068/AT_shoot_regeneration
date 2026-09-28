# Portable path configuration (override with environment variables).
project_dir <- Sys.getenv("SHOOT_REGEN_PROJECT_DIR", unset = ".")
output_root <- Sys.getenv("SHOOT_REGEN_OUTPUT_DIR", unset = file.path(project_dir, "02_results"))
external_data_dir <- Sys.getenv("SHOOT_REGEN_EXTERNAL_DATA_DIR", unset = file.path(project_dir, "external_data"))
work_dir <- Sys.getenv("SHOOT_REGEN_WORK_DIR", unset = ".")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(RANN)
})

set.seed(128)

# -------------------------------------------------------------------------
# Input and output paths
# -------------------------------------------------------------------------
wt_file <- file.path(project_dir, "01_obj/01_WT_snRNA/seurat_subset_C0-2-3-reclustered_dim20_res0.4.rds")
mixed_file <- file.path(project_dir, "01_obj/03_WT_myb31_snRNA/01_all_subset_sct_dim30_res0.3.rds")
outdir <- file.path(project_dir, "02_results/06_WT_C023_D4-D10_subcluster_myb31_mapping")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

stopifnot(packageVersion("Seurat") == "4.4.0")

wt_dates <- c("D4-root", "D6-root", "D8-root", "D10-root")
mutant_dates <- c("myb_SIM4", "myb_SIM6", "myb_SIM8", "myb_SIM10")
cluster_order <- as.character(0:8)
cluster_colors <- setNames(scales::hue_pal()(9), cluster_order)
dims_use <- 1:30

# -------------------------------------------------------------------------
# Utility functions
# -------------------------------------------------------------------------
get_prediction_margin <- function(query, assay_name) {
  score_mat <- GetAssayData(query, assay = assay_name, slot = "data")
  apply(score_mat, 2, function(z) {
    z <- sort(z, decreasing = TRUE)
    if (length(z) < 2) return(NA_real_)
    z[1] - z[2]
  })
}

mean_knn_distance <- function(reference_embeddings, query_embeddings, k = 10) {
  RANN::nn2(
    data = reference_embeddings,
    query = query_embeddings,
    k = k
  )$nn.dists |>
    rowMeans()
}

internal_reference_knn <- function(reference_embeddings, k = 10) {
  # k + 1 is used because the first neighbor of each reference cell is itself.
  d <- RANN::nn2(
    data = reference_embeddings,
    query = reference_embeddings,
    k = k + 1
  )$nn.dists[, -1, drop = FALSE]
  rowMeans(d)
}

# -------------------------------------------------------------------------
# Construct the WT C0/2/3 D4-D10 reference
# -------------------------------------------------------------------------
message("Reading and subsetting the WT reference...")
wt_source <- readRDS(wt_file)
stopifnot(all(c("date", "subcluster") %in% colnames(wt_source@meta.data)))
DefaultAssay(wt_source) <- "RNA"

wt_cells <- rownames(wt_source@meta.data)[as.character(wt_source$date) %in% wt_dates]
wt_source <- subset(wt_source, cells = wt_cells)

wt_counts <- GetAssayData(wt_source, assay = "RNA", slot = "counts")
wt_meta <- wt_source@meta.data
saveRDS(
  wt_counts,
  file.path(outdir, "WT_D4-D10_RNA_counts_sparse_matrix.rds"),
  compress = TRUE
)
write.csv(wt_meta, file.path(outdir, "WT_D4-D10_metadata.csv"))

# Rebuild the reference from RNA counts so that mapping does not depend on SCT.
reference <- CreateSeuratObject(
  counts = wt_counts,
  meta.data = wt_meta,
  assay = "RNA",
  project = "WT_D4_D10_reference"
)
reference$date <- factor(as.character(reference$date), levels = wt_dates)
reference$subcluster <- factor(
  as.character(reference$subcluster),
  levels = cluster_order
)
Idents(reference) <- "subcluster"

reference <- NormalizeData(
  reference,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)
reference <- FindVariableFeatures(
  reference,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)
mapping_features <- VariableFeatures(reference)
reference <- ScaleData(reference, features = mapping_features, verbose = FALSE)
reference <- RunPCA(
  reference,
  features = mapping_features,
  npcs = max(dims_use),
  verbose = FALSE
)
reference <- RunUMAP(
  reference,
  reduction = "pca",
  dims = dims_use,
  return.model = TRUE,
  reduction.name = "wt.umap",
  reduction.key = "WTUMAP_",
  verbose = FALSE
)

# -------------------------------------------------------------------------
# Extract mutant RNA counts and construct a cluster-independent query
# -------------------------------------------------------------------------
message("Extracting mutant RNA counts...")
mixed <- readRDS(mixed_file)
stopifnot("batch" %in% colnames(mixed@meta.data))
mutant_cells <- rownames(mixed@meta.data)[
  tolower(as.character(mixed$batch)) == "mutant"
]

mutant_counts <- GetAssayData(
  mixed,
  assay = "RNA",
  slot = "counts"
)[, mutant_cells, drop = FALSE]
mutant_meta <- mixed@meta.data[mutant_cells, , drop = FALSE]

saveRDS(
  mutant_counts,
  file.path(outdir, "mutant_RNA_counts_sparse_matrix.rds"),
  compress = TRUE
)
write.csv(mutant_meta, file.path(outdir, "mutant_extracted_metadata.csv"))

query <- CreateSeuratObject(
  counts = mutant_counts,
  meta.data = mutant_meta,
  assay = "RNA",
  project = "myb31_query"
)
query <- NormalizeData(
  query,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

shared_genes <- intersect(rownames(reference), rownames(query))
mapping_features <- intersect(mapping_features, shared_genes)
if (length(mapping_features) < 500) {
  stop("Too few shared mapping features: ", length(mapping_features))
}

# -------------------------------------------------------------------------
# Anchor finding, label transfer and UMAP projection
# -------------------------------------------------------------------------
message("Finding transfer anchors...")
anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  normalization.method = "LogNormalize",
  reference.assay = "RNA",
  query.assay = "RNA",
  features = mapping_features,
  reference.reduction = "pca",
  reduction = "pcaproject",
  dims = dims_use,
  k.anchor = 5,
  verbose = TRUE
)

message("Transferring WT subcluster and time labels...")
query <- MapQuery(
  anchorset = anchors,
  reference = reference,
  query = query,
  refdata = list(
    subcluster = "subcluster",
    date = "date"
  ),
  reference.reduction = "pca",
  reduction.model = "wt.umap",
  verbose = TRUE
)

query$WT_pred_cluster <- factor(
  as.character(query$predicted.subcluster),
  levels = cluster_order
)
query$WT_score_cluster <- query$predicted.subcluster.score
query$WT_pred_date <- factor(
  as.character(query$predicted.date),
  levels = wt_dates
)
query$WT_score_date <- query$predicted.date.score

# -------------------------------------------------------------------------
# Reference-membership assessment
# -------------------------------------------------------------------------
message("Calculating prediction margins and reference-manifold distances...")
cluster_score_assay <- "prediction.score.subcluster"
query$WT_cluster_score_margin <- get_prediction_margin(
  query,
  cluster_score_assay
)[colnames(query)]

reference_pca <- Embeddings(reference, "pca")[, dims_use, drop = FALSE]
query_pca <- Embeddings(query, "ref.pca")[, dims_use, drop = FALSE]

wt_internal_distance <- internal_reference_knn(reference_pca, k = 10)
wt_q95 <- unname(quantile(wt_internal_distance, 0.95, na.rm = TRUE))
wt_q99 <- unname(quantile(wt_internal_distance, 0.99, na.rm = TRUE))
query$WT_mean_knn_distance <- mean_knn_distance(
  reference_embeddings = reference_pca,
  query_embeddings = query_pca,
  k = 10
)

query$WT_reference_status <- "uncertain"
query$WT_reference_status[
  query$WT_score_cluster >= 0.5 &
    query$WT_cluster_score_margin >= 0.1 &
    query$WT_mean_knn_distance <= wt_q95
] <- "in_reference"
query$WT_reference_status[
  query$WT_score_cluster < 0.4 |
    query$WT_mean_knn_distance > wt_q99
] <- "out_of_reference"
query$WT_reference_status <- factor(
  query$WT_reference_status,
  levels = c("in_reference", "out_of_reference", "uncertain")
)

saveRDS(
  reference,
  file.path(outdir, "WT_D4-D10_reference_processed.rds"),
  compress = TRUE
)
saveRDS(
  query,
  file.path(outdir, "mutant_mapped_to_WT_D4-D10.rds"),
  compress = TRUE
)
write.csv(query@meta.data, file.path(outdir, "mutant_mapping_metadata.csv"))

# -------------------------------------------------------------------------
# Summary tables
# -------------------------------------------------------------------------
mapping_summary <- as.data.frame(table(
  mutant_time = factor(query$date, levels = mutant_dates),
  predicted_cluster = factor(query$WT_pred_cluster, levels = cluster_order),
  predicted_WT_date = factor(query$WT_pred_date, levels = wt_dates),
  reference_status = query$WT_reference_status
))
mapping_summary <- mapping_summary[mapping_summary$Freq > 0, , drop = FALSE]
write.csv(
  mapping_summary,
  file.path(outdir, "mapping_summary.csv"),
  row.names = FALSE
)

status_summary <- as.data.frame(table(query$WT_reference_status))
status_summary$proportion <- status_summary$n_cells / sum(status_summary$n_cells)
write.csv(
  status_summary,
  file.path(outdir, "reference_status_summary.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Figure 01: WT reference and mapped mutant UMAP
# -------------------------------------------------------------------------
p_ref <- DimPlot(
  reference,
  reduction = "wt.umap",
  group.by = "subcluster",
  cols = cluster_colors,
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  scale_color_manual(
    values = cluster_colors,
    limits = cluster_order,
    drop = FALSE,
    name = "WT subcluster"
  ) +
  ggtitle("WT C0/2/3 D4-D10 reference")

p_query <- DimPlot(
  query,
  reduction = "ref.umap",
  group.by = "WT_pred_cluster",
  cols = cluster_colors,
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  scale_color_manual(
    values = cluster_colors,
    limits = cluster_order,
    drop = FALSE,
    name = "WT subcluster"
  ) +
  ggtitle("myb31 mapped to WT D4-D10")

p01 <- p_ref + p_query + plot_layout(guides = "collect")
ggsave(
  file.path(outdir, "01_WT_mutant_cluster_UMAP.pdf"),
  p01, width = 14, height = 6.5, device = cairo_pdf
)
ggsave(
  file.path(outdir, "01_WT_mutant_cluster_UMAP.png"),
  p01, width = 14, height = 6.5, dpi = 350, bg = "white"
)

# -------------------------------------------------------------------------
# Figure 02: mapped clusters split by reference status
# -------------------------------------------------------------------------
p02 <- DimPlot(
  query,
  reduction = "ref.umap",
  group.by = "WT_pred_cluster",
  split.by = "WT_reference_status",
  cols = cluster_colors,
  raster = TRUE,
  ncol = 3
) +
  scale_color_manual(
    values = cluster_colors,
    limits = cluster_order,
    drop = FALSE,
    name = "WT subcluster"
  )
ggsave(
  file.path(outdir, "02_cluster_by_reference_status.pdf"),
  p02, width = 15, height = 5.5, device = cairo_pdf
)
ggsave(
  file.path(outdir, "02_cluster_by_reference_status.png"),
  p02, width = 15, height = 5.5, dpi = 350, bg = "white"
)

# -------------------------------------------------------------------------
# Figures 03 and 04: predicted WT time and reference status
# -------------------------------------------------------------------------
p03 <- DimPlot(
  query,
  reduction = "ref.umap",
  group.by = "WT_pred_date",
  raster = TRUE
) +
  ggtitle("Predicted WT-equivalent time")
ggsave(
  file.path(outdir, "03_predicted_WT_date.pdf"),
  p03, width = 7, height = 6, device = cairo_pdf
)
ggsave(
  file.path(outdir, "03_predicted_WT_date.png"),
  p03, width = 7, height = 6, dpi = 350, bg = "white"
)

p04 <- DimPlot(
  query,
  reduction = "ref.umap",
  group.by = "WT_reference_status",
  cols = c(
    in_reference = "#D95F59",
    out_of_reference = "#91B9D2",
    uncertain = "#F3C677"
  ),
  raster = TRUE
) +
  ggtitle("WT C0/2/3 D4-D10 reference membership")
ggsave(
  file.path(outdir, "04_reference_status.pdf"),
  p04, width = 7, height = 6, device = cairo_pdf
)
ggsave(
  file.path(outdir, "04_reference_status.png"),
  p04, width = 7, height = 6, dpi = 350, bg = "white"
)

# -------------------------------------------------------------------------
# Figure 05: within-time WT versus mutant cluster composition
# -------------------------------------------------------------------------
wt_comp <- as.data.frame(table(
  dataset_time = factor(reference$date, levels = wt_dates),
  cluster = factor(reference$subcluster, levels = cluster_order)
))
wt_comp$source <- "WT observed"

query_in <- query@meta.data[
  query$WT_reference_status == "in_reference",
  ,
  drop = FALSE
]
mut_comp <- as.data.frame(table(
  dataset_time = factor(query_in$date, levels = mutant_dates),
  cluster = factor(query_in$WT_pred_cluster, levels = cluster_order)
))
mut_comp$source <- "myb31 predicted (in-reference)"

comp <- rbind(wt_comp, mut_comp)
comp$proportion <- ave(
  comp$Freq,
  comp$dataset_time,
  FUN = function(z) if (sum(z) > 0) z / sum(z) else rep(0, length(z))
)
comp$display_time <- as.character(comp$dataset_time)
comp$label <- ifelse(comp$proportion >= 0.01,
                     sprintf("%.1f", 100 * comp$proportion), "")
comp$cluster <- factor(comp$cluster, levels = cluster_order)
comp$source <- factor(
  comp$source,
  levels = c("WT observed", "myb31 predicted (in-reference)")
)
display_order <- c(wt_dates, mutant_dates)
comp$display_time <- factor(comp$display_time, levels = rev(display_order))

write.csv(
  comp,
  file.path(outdir, "cluster_comparison_heatmap_data.csv"),
  row.names = FALSE
)

p05 <- ggplot(comp, aes(cluster, display_time, fill = proportion)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = label), size = 3) +
  scale_x_discrete(limits = cluster_order, drop = FALSE) +
  scale_fill_gradientn(
    colors = c("#FFF7EC", "#FDD49E", "#FC8D59", "#B30000"),
    labels = scales::percent,
    name = "Within-time\nproportion"
  ) +
  facet_grid(source ~ ., scales = "free_y", space = "free_y") +
  labs(x = "WT subcluster", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(face = "bold"),
    panel.grid = element_blank(),
    strip.text.y = element_text(angle = 0, face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA)
  )
ggsave(
  file.path(outdir, "05_WT_mutant_cluster_comparison_heatmap.pdf"),
  p05, width = 10.5, height = 7, device = cairo_pdf
)
ggsave(
  file.path(outdir, "05_WT_mutant_cluster_comparison_heatmap.png"),
  p05, width = 10.5, height = 7, dpi = 350, bg = "white"
)

# -------------------------------------------------------------------------
# Figure 06: median mapping confidence by mutant time and predicted cluster
# -------------------------------------------------------------------------
confidence <- aggregate(
  query$WT_score_cluster,
  by = list(
    date = factor(query$date, levels = mutant_dates),
    WT_pred_cluster = factor(query$WT_pred_cluster, levels = cluster_order)
  ),
  FUN = function(z) c(median = median(z), n_cells = length(z))
)
confidence <- data.frame(
  date = confidence$date,
  WT_pred_cluster = confidence$WT_pred_cluster,
  WT_score_cluster = confidence$x[, "median"],
  n_cells = confidence$x[, "n_cells"]
)
write.csv(
  confidence,
  file.path(outdir, "cluster_confidence_heatmap_data.csv"),
  row.names = FALSE
)

confidence$date <- factor(confidence$date, levels = rev(mutant_dates))
confidence$WT_pred_cluster <- factor(
  confidence$WT_pred_cluster,
  levels = cluster_order
)

p06 <- ggplot(
  confidence,
  aes(WT_pred_cluster, date, fill = WT_score_cluster)
) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f", WT_score_cluster)), size = 3) +
  scale_x_discrete(limits = cluster_order, drop = FALSE) +
  scale_fill_gradientn(
    colors = c("#FFF7EC", "#FDD49E", "#FC8D59", "#B30000"),
    limits = c(0, 1),
    na.value = "grey92",
    name = "Median transfer\nscore"
  ) +
  labs(x = "Predicted WT subcluster", y = "Mutant time") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(face = "bold"),
    panel.grid = element_blank()
  )
ggsave(
  file.path(outdir, "06_cluster_mapping_confidence_heatmap.pdf"),
  p06, width = 10.5, height = 5, device = cairo_pdf
)
ggsave(
  file.path(outdir, "06_cluster_mapping_confidence_heatmap.png"),
  p06, width = 10.5, height = 5, dpi = 350, bg = "white"
)

# -------------------------------------------------------------------------
# Run information
# -------------------------------------------------------------------------
writeLines(c(
  paste("R version:", R.version.string),
  paste("Seurat version:", packageVersion("Seurat")),
  paste("Cleaned WT input:", wt_file),
  paste("Mixed input:", mixed_file),
  paste("WT dates retained:", paste(wt_dates, collapse = ", ")),
  paste("WT reference cells:", ncol(reference)),
  paste("Mutant cells:", ncol(query)),
  paste("Shared genes:", length(shared_genes)),
  paste("Mapping features:", length(mapping_features)),
  paste("WT internal 10-NN q95:", wt_q95),
  paste("WT internal 10-NN q99:", wt_q99),
  "Primary transferred label: subcluster.",
  "subcluster display order: 0,1,2,3,4,5,6,7,8.",
  "in_reference: cluster score>=0.5, margin>=0.1, distance<=WT q95.",
  "out_of_reference: cluster score<0.4 or distance>WT q99."
), file.path(outdir, "run_information.txt"))

writeLines(
  capture.output(sessionInfo()),
  file.path(outdir, "sessionInfo.txt")
)

message("Complete: ", outdir)

# Apply the standardized figure palette, cluster order, and publication theme.
# The companion script is saved in the same output directory.
source(file.path(outdir, "restyle_C023_subcluster_mapping_figures.R"))
