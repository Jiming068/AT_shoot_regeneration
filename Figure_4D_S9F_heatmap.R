#!/usr/bin/env Rscript

# Reproduce WT–myb31 expression summaries and the two split publication heatmaps
# in a single run. The original two-step scripts are retained for provenance.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
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

input_rds <- get_arg("input-rds", NULL)
gene_file <- get_arg("gene-list", file.path(script_dir, "fig2c_genes_selected_7_13.tsv"))
output_dir <- get_arg("outdir", script_dir)

if (is.null(input_rds)) {
  stop(
    paste0(
      "Provide --input-rds=/path/to/integrated_seurat_object.rds\n",
      "Optional: --gene-list=/path/to/gene_list.tsv --outdir=/path/to/output"
    )
  )
}
if (!file.exists(input_rds)) stop("Input RDS does not exist: ", input_rds)
if (!file.exists(gene_file)) stop("Gene-list file does not exist: ", gene_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

days <- c("D4", "D6", "D8", "D10")
conditions <- c(paste0("WT_", days), paste0("myb31_", days))
condition_labels <- c(
  WT_D4 = "WT\nD4", WT_D6 = "WT\nD6", WT_D8 = "WT\nD8", WT_D10 = "WT\nD10",
  myb31_D4 = "myb31\nD4", myb31_D6 = "myb31\nD6",
  myb31_D8 = "myb31\nD8", myb31_D10 = "myb31\nD10"
)

# -----------------------------------------------------------------------------
# 1. Read and validate the selected gene list
# -----------------------------------------------------------------------------
genes <- read.delim(
  gene_file,
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  fill = TRUE
)
if (ncol(genes) < 3L) stop("Gene-list file must contain gene, symbol, and category columns.")
genes <- genes[, 1:3, drop = FALSE]
genes[] <- lapply(genes, trimws)
genes <- genes[nzchar(genes$gene), , drop = FALSE]
genes$gene <- toupper(genes$gene)
genes <- genes[!duplicated(genes$gene), , drop = FALSE]
genes$gene_order <- seq_len(nrow(genes))
genes$category <- sub(
  "^Rooot/regeneration competence$",
  "Root/regeneration competence",
  genes$category
)

# -----------------------------------------------------------------------------
# 2. Compute genotype-by-day mean expression, detection rate, and gene Z-scores
# -----------------------------------------------------------------------------
message("Reading Seurat object: ", input_rds)
obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"
expr <- if (packageVersion("SeuratObject") >= "5.0.0") {
  LayerData(obj, assay = "RNA", layer = "data")
} else {
  GetAssayData(obj, assay = "RNA", slot = "data")
}
meta <- obj@meta.data

if (!all(rownames(meta) == colnames(expr))) {
  stop("Metadata and expression-matrix cell order differ.")
}
missing_meta <- setdiff(c("batch", "date"), colnames(meta))
if (length(missing_meta)) {
  stop("Metadata is missing required columns: ", paste(missing_meta, collapse = ", "))
}

meta$genotype <- ifelse(tolower(as.character(meta$batch)) == "mutant", "myb31", "WT")
meta$day <- sub(".*SIM", "D", as.character(meta$date))
meta$condition <- paste(meta$genotype, meta$day, sep = "_")

present <- genes$gene %in% rownames(expr)
gene_detection <- transform(genes, detected_in_RNA_data = present)
write.csv(
  gene_detection,
  file.path(output_dir, "01_gene_detection_in_RNA_data.csv"),
  row.names = FALSE,
  quote = TRUE
)
if (!all(present)) {
  warning("Genes absent from RNA/data: ", paste(genes$gene[!present], collapse = ", "))
}

genes_use <- genes[present, , drop = FALSE]
if (!nrow(genes_use)) stop("None of the requested genes was detected in RNA/data.")
expr_use <- expr[genes_use$gene, , drop = FALSE]
mean_matrix <- matrix(
  NA_real_,
  nrow = nrow(genes_use),
  ncol = length(conditions),
  dimnames = list(genes_use$gene, conditions)
)
pct_matrix <- mean_matrix
cell_counts <- setNames(integer(length(conditions)), conditions)

for (condition in conditions) {
  cells <- rownames(meta)[meta$condition == condition]
  cell_counts[[condition]] <- length(cells)
  if (length(cells)) {
    mean_matrix[, condition] <- Matrix::rowMeans(expr_use[, cells, drop = FALSE])
    pct_matrix[, condition] <- 100 * Matrix::rowMeans(expr_use[, cells, drop = FALSE] > 0)
  }
}
if (any(cell_counts == 0L)) {
  stop(
    "The following genotype-by-day conditions contain no cells: ",
    paste(names(cell_counts)[cell_counts == 0L], collapse = ", ")
  )
}

# R scale() uses the sample standard deviation (denominator n - 1).
z_matrix <- t(scale(t(mean_matrix), center = TRUE, scale = TRUE))
z_matrix[!is.finite(z_matrix)] <- 0

mean_table <- data.frame(genes_use, mean_matrix, check.names = FALSE)
pct_table <- data.frame(genes_use, pct_matrix, check.names = FALSE)
zscore_table <- data.frame(genes_use, z_matrix, check.names = FALSE)

write.csv(
  mean_table,
  file.path(output_dir, "02_condition_mean_RNA_expression.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  pct_table,
  file.path(output_dir, "03_condition_pct_expressed.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  zscore_table,
  file.path(output_dir, "04_gene_wise_Zscore_matrix.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  data.frame(condition = conditions, n_cells = unname(cell_counts)),
  file.path(output_dir, "05_condition_cell_counts.csv"),
  row.names = FALSE,
  quote = TRUE
)

writeLines(c(
  paste0("Input RDS: ", normalizePath(input_rds)),
  paste0("Input RDS MD5: ", unname(tools::md5sum(input_rds))),
  paste0("Gene list: ", normalizePath(gene_file)),
  "Assay/layer: RNA/data (log-normalized expression).",
  paste0("Cells in expression matrix: ", ncol(expr)),
  paste0("Genes requested: ", nrow(genes)),
  paste0("Genes detected: ", nrow(genes_use)),
  "Aggregation: arithmetic mean of cell-level RNA/data values within each genotype-by-day condition.",
  "Scaling: gene-wise Z-score across eight conditions using R scale() and sample SD (n - 1).",
  paste0("Condition order: ", paste(conditions, collapse = ", "))
), file.path(output_dir, "06_expression_analysis_information.txt"))

# Release the large Seurat object before plotting.
rm(obj, expr, expr_use)
invisible(gc())

# -----------------------------------------------------------------------------
# 3. Split the genes and generate publication-ready PDF and PNG heatmaps
# -----------------------------------------------------------------------------
make_heatmap <- function(data, stem, width, height) {
  matrix_data <- as.matrix(data[, conditions, drop = FALSE])
  rownames(matrix_data) <- data$gene
  heat <- as.data.frame(as.table(matrix_data), stringsAsFactors = FALSE)
  heat <- merge(
    heat,
    data[, c("gene", "symbol", "category", "gene_order")],
    by = "gene",
    sort = FALSE
  )
  heat$symbol <- factor(heat$symbol, levels = rev(data$symbol))
  heat$category <- factor(heat$category, levels = unique(data$category))
  heat$condition <- factor(
    heat$condition,
    levels = conditions,
    labels = unname(condition_labels[conditions])
  )

  plot <- ggplot(heat, aes(condition, symbol, fill = z_score)) +
    geom_tile(color = "white", linewidth = 0.28) +
    geom_vline(xintercept = 4.5, color = "white", linewidth = 1.35) +
    facet_grid(category ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_gradient2(
      low = "#A8C7DF",
      mid = "#F7F7F7",
      high = "#E98C85",
      midpoint = 0,
      limits = c(-2.5, 2.5),
      oob = squish,
      name = "Gene-wise\nZ-score"
    ) +
    theme_minimal(base_size = 8.8, base_family = "Arial") +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(color = "#222222", face = "bold", size = 8.2),
      axis.text.y = element_text(color = "#222222", face = "italic", size = 7.8),
      strip.placement = "outside",
      strip.background.y = element_rect(fill = "#F0F2F4", color = NA),
      strip.text.y.left = element_text(
        angle = 0, face = "bold", color = "#333333", size = 7.2
      ),
      legend.title = element_text(face = "bold", size = 8.4),
      legend.text = element_text(size = 7.8),
      panel.spacing.y = grid::unit(0.65, "mm"),
      plot.margin = margin(4, 5, 4, 4),
      plot.background = element_rect(fill = "white", color = NA)
    )

  pdf_file <- file.path(output_dir, paste0(stem, ".pdf"))
  png_file <- file.path(output_dir, paste0(stem, ".png"))

  cairo_pdf(pdf_file, width = width, height = height, family = "Arial", onefile = FALSE)
  print(plot)
  invisible(dev.off())

  png(
    filename = png_file,
    width = width,
    height = height,
    units = "in",
    res = 600,
    type = "cairo",
    family = "Arial",
    bg = "white"
  )
  print(plot)
  invisible(dev.off())

  write.csv(
    data,
    file.path(output_dir, paste0(stem, "_source_matrix.csv")),
    row.names = FALSE,
    quote = TRUE
  )
}

hormone <- zscore_table[zscore_table$category %in% c("Auxin", "Cytokinin"), , drop = FALSE]
remaining <- zscore_table[!zscore_table$category %in% c("Auxin", "Cytokinin"), , drop = FALSE]
if (!nrow(hormone) || !nrow(remaining)) stop("Both split panels must contain genes.")

make_heatmap(
  hormone,
  "02_auxin_cytokinin_WT_myb31_expression_Zscore_heatmap_publication",
  width = 6.55,
  height = 4.25
)
make_heatmap(
  remaining,
  "03_remaining_genes_WT_myb31_expression_Zscore_heatmap_publication",
  width = 6.55,
  height = 6.10
)

writeLines(c(
  "Input Z-score matrix: generated in memory from the integrated Seurat object.",
  "Hormone panel definition: categories Auxin and Cytokinin.",
  "Remaining panel definition: all other categories.",
  paste0("Hormone genes: ", nrow(hormone)),
  paste0("Remaining genes: ", nrow(remaining)),
  "Column order: WT D4-D10 followed by myb31 D4-D10.",
  "Color limits: -2.5 to 2.5; values outside this range are squished.",
  "A white vertical divider separates WT and myb31 conditions."
), file.path(output_dir, "08_split_heatmap_information.txt"))

capture.output(sessionInfo(), file = file.path(output_dir, "07_expression_analysis_sessionInfo.txt"))
capture.output(sessionInfo(), file = file.path(output_dir, "09_split_heatmap_sessionInfo.txt"))

message("Completed expression summaries and both publication heatmaps in: ", output_dir)
