#!/usr/bin/env Rscript

# Portable path configuration (override with environment variables).
project_dir <- Sys.getenv("SHOOT_REGEN_PROJECT_DIR", unset = ".")
output_root <- Sys.getenv("SHOOT_REGEN_OUTPUT_DIR", unset = file.path(project_dir, "02_results"))
external_data_dir <- Sys.getenv("SHOOT_REGEN_EXTERNAL_DATA_DIR", unset = file.path(project_dir, "external_data"))
work_dir <- Sys.getenv("SHOOT_REGEN_WORK_DIR", unset = ".")

# Complete reproducible workflow for the SIM6 eight-class spatial-reference
# mapping and the time-niche composition heatmap in result directory 14.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

set.seed(128)
stopifnot(packageVersion("Seurat") == "4.4.0")

# -----------------------------------------------------------------------------
# Inputs and output directory
# -----------------------------------------------------------------------------
reference_file <- file.path(project_dir, "01_obj/02_WT_stereo/Stereo-seq_SIM6-SIM10_final.rds")
mixed_file <- paste0(
  file.path(project_dir, "01_obj/"),
  "03_WT_myb31_snRNA/01_all_subset_sct_dim30_res0.3.rds"
)
outdir <- paste0(
  file.path(project_dir, "02_results/"),
  "14_D6_8class_mutant_RNA_mapping_no_clusters"
)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ref_levels <- c(
  "SIM6_sub-1 core",
  "SIM6_sub-1 bdry",
  "SIM6_sub-3 shell",
  "SIM6_10 shell",
  "SIM6_3",
  "SIM6_4",
  "SIM6_5",
  "SIM6_8"
)
mutant_dates <- c("myb_SIM4", "myb_SIM6", "myb_SIM8", "myb_SIM10")
dims_use <- 1:30

class_colors <- c(
  "SIM6_sub-1 core" = "#D95F59",
  "SIM6_sub-1 bdry" = "#F3C677",
  "SIM6_sub-3 shell" = "#91B9D2",
  "SIM6_10 shell" = "#8DB58C",
  "SIM6_3" = "#E6A0C4",
  "SIM6_4" = "#B8A6D9",
  "SIM6_5" = "#E69F63",
  "SIM6_8" = "#79B7B2"
)
status_levels <- c("in_reference", "uncertain", "out_of_reference")
status_colors <- c(
  "in_reference" = "#D95F59",
  "uncertain" = "#F3C677",
  "out_of_reference" = "#91B9D2"
)
warm_gradient <- c("#D9D9D9", "#FFF7EC", "#FDD49E", "#FC8D59", "#B30000")

base_theme <- theme_minimal(base_size = 12, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    panel.grid = element_blank(),
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

save_pair <- function(plot, stem, width, height) {
  ggsave(
    file.path(outdir, paste0(stem, ".pdf")), plot,
    width = width, height = height, device = cairo_pdf
  )
  ggsave(
    file.path(outdir, paste0(stem, ".png")), plot,
    width = width, height = height, dpi = 400, bg = "white"
  )
}

# -----------------------------------------------------------------------------
# Construct the eight-class WT D6 spatial reference from RNA counts
# -----------------------------------------------------------------------------
message("Reading the WT spatial reference...")
spatial_source <- readRDS(reference_file)
stopifnot("time_niche_class" %in% colnames(spatial_source@meta.data))
DefaultAssay(spatial_source) <- "RNA"

ref_cells <- rownames(spatial_source@meta.data)[
  as.character(spatial_source$time_niche_class) %in% ref_levels
]
if (!length(ref_cells)) {
  stop("No requested D6 time_niche_class cells were found in the reference.")
}
spatial_source <- subset(spatial_source, cells = ref_cells)

ref_counts <- GetAssayData(spatial_source, assay = "RNA", slot = "counts")
ref_meta <- spatial_source@meta.data

reference <- CreateSeuratObject(
  counts = ref_counts,
  meta.data = ref_meta,
  assay = "RNA",
  project = "WT_D6_8class_spatial_reference"
)
reference$time_niche_class <- factor(
  as.character(reference$time_niche_class),
  levels = ref_levels
)
Idents(reference) <- "time_niche_class"

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

reference_counts <- as.data.frame(table(
  factor(reference$time_niche_class, levels = ref_levels)
))
write.csv(
  reference_counts,
  file.path(outdir, "WT_reference_class_counts.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Rebuild the mutant query from RNA counts without cluster information
# -----------------------------------------------------------------------------
message("Extracting mutant RNA counts without using mutant clusters...")
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

# Explicitly remove all pre-existing cluster and graph-resolution metadata.
drop_meta <- grep(
  "cluster|snn_res",
  colnames(mutant_meta),
  ignore.case = TRUE,
  value = TRUE
)
mutant_meta <- mutant_meta[, setdiff(colnames(mutant_meta), drop_meta), drop = FALSE]

saveRDS(
  mutant_counts,
  file.path(outdir, "mutant_RNA_counts_sparse_matrix.rds"),
  compress = TRUE
)

query <- CreateSeuratObject(
  counts = mutant_counts,
  meta.data = mutant_meta,
  assay = "RNA",
  project = "myb31_query_no_clusters"
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

# -----------------------------------------------------------------------------
# Anchor mapping and time-niche label transfer
# -----------------------------------------------------------------------------
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

message("Transferring the WT D6 time-niche labels...")
query <- MapQuery(
  anchorset = anchors,
  reference = reference,
  query = query,
  refdata = list(time_niche_class = "time_niche_class"),
  reference.reduction = "pca",
  reduction.model = "wt.umap",
  verbose = TRUE
)

query$WT_pred_class <- factor(
  as.character(query$predicted.time_niche_class),
  levels = ref_levels
)
query$WT_score_class <- query$predicted.time_niche_class.score

# Reproduce the score-only support categories recorded for this analysis.
query$reference_status <- "uncertain"
query$reference_status[query$WT_score_class >= 0.50] <- "in_reference"
query$reference_status[query$WT_score_class < 0.30] <- "out_of_reference"
query$reference_status <- factor(query$reference_status, levels = status_levels)

saveRDS(
  query,
  file.path(outdir, "mutant_RNA_mapped_D6_8class_no_clusters.rds"),
  compress = TRUE
)
write.csv(
  query@meta.data,
  file.path(outdir, "mutant_extracted_and_mapping_metadata_no_clusters.csv")
)

# -----------------------------------------------------------------------------
# Mapping and reference-support summaries
# -----------------------------------------------------------------------------
mapping_summary <- as.data.frame(table(
  mutant_time = factor(as.character(query$date), levels = mutant_dates),
  predicted_time_niche_class = factor(query$WT_pred_class, levels = ref_levels),
  reference_status = query$reference_status
))
mapping_summary <- mapping_summary[mapping_summary$Freq > 0, , drop = FALSE]
write.csv(
  mapping_summary,
  file.path(outdir, "mapping_summary.csv"),
  row.names = FALSE
)

status_summary <- as.data.frame(table(query$reference_status))
status_summary$proportion <- status_summary$n_cells / sum(status_summary$n_cells)
write.csv(
  status_summary,
  file.path(outdir, "reference_status_summary.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Figure 01: WT spatial reference and mapped mutant UMAP
# -----------------------------------------------------------------------------
p_ref <- DimPlot(
  reference,
  reduction = "wt.umap",
  group.by = "time_niche_class",
  cols = class_colors,
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  scale_colour_manual(values = class_colors, limits = ref_levels, drop = FALSE) +
  ggtitle("WT D6 eight-class spatial reference") +
  base_theme

p_query <- DimPlot(
  query,
  reduction = "ref.umap",
  group.by = "WT_pred_class",
  cols = class_colors,
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  scale_colour_manual(values = class_colors, limits = ref_levels, drop = FALSE) +
  ggtitle("myb31 RNA mapped to WT D6 spatial reference") +
  base_theme

save_pair(
  p_ref + p_query + plot_layout(guides = "collect") &
    theme(legend.position = "right"),
  "01_WT_D6_8class_mutant_RNA_UMAP",
  14,
  6.2
)

# -----------------------------------------------------------------------------
# Figure 02: reference-support category in the mapped mutant query
# -----------------------------------------------------------------------------
p02 <- DimPlot(
  query,
  reduction = "ref.umap",
  group.by = "reference_status",
  cols = status_colors,
  raster = TRUE
) +
  scale_colour_manual(values = status_colors, limits = status_levels, drop = FALSE) +
  ggtitle("WT D6 spatial reference support") +
  base_theme
save_pair(p02, "02_reference_status", 7.2, 6)

# -----------------------------------------------------------------------------
# Figure 03: WT spatial versus mutant time-niche composition
# -----------------------------------------------------------------------------
wt_comp <- as.data.frame(table(
  dataset_time = rep("WT_D6_spatial", ncol(reference)),
  time_niche_class = factor(reference$time_niche_class, levels = ref_levels)
))
wt_comp$source <- "WT observed"

query_in <- query@meta.data[
  query$reference_status == "in_reference",
  ,
  drop = FALSE
]
mut_comp <- as.data.frame(table(
  dataset_time = factor(as.character(query_in$date), levels = mutant_dates),
  time_niche_class = factor(query_in$WT_pred_class, levels = ref_levels)
))
mut_comp$source <- "myb31 predicted (in-reference)"

composition <- rbind(wt_comp, mut_comp)
composition$proportion <- ave(
  composition$Freq,
  composition$dataset_time,
  FUN = function(v) if (sum(v) > 0) v / sum(v) else rep(0, length(v))
)
composition$label <- ifelse(
  composition$proportion >= 0.01,
  sprintf("%.1f", 100 * composition$proportion),
  ""
)
composition$time_niche_class <- factor(
  composition$time_niche_class,
  levels = ref_levels
)
composition$dataset_time <- factor(
  as.character(composition$dataset_time),
  levels = c("WT_D6_spatial", mutant_dates)
)
write.csv(
  composition,
  file.path(outdir, "time_niche_comparison_heatmap_data.csv"),
  row.names = FALSE
)

p03 <- ggplot(
  composition,
  aes(time_niche_class, dataset_time, fill = proportion)
) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = label), family = "Arial", size = 3.2) +
  scale_x_discrete(limits = ref_levels, drop = FALSE) +
  scale_y_discrete(
    limits = c("WT_D6_spatial", mutant_dates),
    drop = FALSE
  ) +
  scale_fill_gradientn(
    colors = warm_gradient,
    values = c(0, 0.001, 0.25, 0.60, 1),
    limits = c(0, max(composition$proportion, na.rm = TRUE)),
    labels = percent,
    name = "Within-time\nproportion"
  ) +
  labs(x = "Reference time-niche class", y = NULL) +
  base_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
save_pair(p03, "03_WT_mutant_time_niche_comparison_heatmap", 9, 5.5)

# -----------------------------------------------------------------------------
# Figure 04: median transfer score by mutant time and predicted class
# -----------------------------------------------------------------------------
confidence <- aggregate(
  query$WT_score_class,
  by = list(
    date = factor(as.character(query$date), levels = mutant_dates),
    time_niche_class = factor(query$WT_pred_class, levels = ref_levels)
  ),
  FUN = function(z) c(median = median(z), n_cells = length(z))
)
confidence <- data.frame(
  date = confidence$date,
  time_niche_class = confidence$time_niche_class,
  prediction_score = confidence$x[, "median"],
  n_cells = confidence$x[, "n_cells"]
)
write.csv(
  confidence,
  file.path(outdir, "time_niche_confidence_heatmap_data.csv"),
  row.names = FALSE
)

p04 <- ggplot(
  confidence,
  aes(time_niche_class, date, fill = prediction_score)
) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f", prediction_score)), size = 3) +
  scale_x_discrete(limits = ref_levels, drop = FALSE) +
  scale_fill_gradientn(
    colors = c("#FFF7EC", "#FDD49E", "#FC8D59", "#B30000"),
    limits = c(0, 1),
    na.value = "grey92",
    name = "Median transfer\nscore"
  ) +
  labs(x = "Predicted WT D6 time-niche class", y = "Mutant time") +
  base_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
save_pair(p04, "04_time_niche_mapping_confidence_heatmap", 9, 5.4)

# Reproduce the final editable vector heatmaps with and without values.
vector_script <- file.path(
  outdir,
  "reproduce_D6_8class_time_niche_heatmap_vector.R"
)
if (file.exists(vector_script)) {
  source(vector_script, local = new.env(parent = globalenv()))
}

# -----------------------------------------------------------------------------
# Reproducibility record
# -----------------------------------------------------------------------------
n_anchors <- if ("anchors" %in% slotNames(anchors)) {
  nrow(anchors@anchors)
} else {
  NA_integer_
}
writeLines(c(
  paste("R version:", R.version.string),
  paste("Seurat version:", packageVersion("Seurat")),
  paste("Reference input:", reference_file),
  paste("Mutant input:", mixed_file),
  paste("WT reference cells:", ncol(reference)),
  paste("Reference labels:", length(ref_levels)),
  paste("Mutant RNA cells:", ncol(query)),
  "Mutant object contains RNA data and mapping results only; cluster metadata were removed.",
  paste("Shared genes:", length(shared_genes)),
  paste("Mapping features:", length(mapping_features)),
  paste("Anchors retained:", n_anchors),
  "Prediction score column: predicted.time_niche_class.score",
  paste(
    "Median prediction score:",
    signif(median(query$WT_score_class, na.rm = TRUE), 4)
  ),
  "Status thresholds: in_reference >= 0.5; uncertain 0.3-0.5; out_of_reference < 0.3.",
  "Interpretation: cross-platform expression similarity, not strict cell identity."
), file.path(outdir, "run_information.txt"))

writeLines(
  capture.output(sessionInfo()),
  file.path(outdir, "sessionInfo.txt")
)

message("Complete: ", outdir)
