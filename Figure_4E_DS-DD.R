#!/usr/bin/env Rscript

# Reproduce Figure 2c selected-gene DS/DD heatmap from the genome-wide,
# sample-aware, state-resolved DS/DD result table.
#
# Required input:
#   --dsdd-results=/path/to/03_genomewide_state_resolved_DS_DD.csv.gz
# Optional inputs:
#   --gene-list=/path/to/fig2c_genes_selected_7_20.tsv
#   --outdir=/path/to/output_directory

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

dsdd_file <- get_arg("dsdd-results", NULL)
gene_file <- get_arg(
  "gene-list",
  file.path(script_dir, "fig2c_genes_selected_7_20.tsv")
)
output_dir <- get_arg("outdir", script_dir)

if (is.null(dsdd_file)) {
  stop(
    paste0(
      "Provide the genome-wide DS/DD table with:\n",
      "  --dsdd-results=/path/to/03_genomewide_state_resolved_DS_DD.csv.gz"
    )
  )
}
if (!file.exists(dsdd_file)) stop("DS/DD result file does not exist: ", dsdd_file)
if (!file.exists(gene_file)) stop("Gene-list file does not exist: ", gene_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

days <- c("D4", "D6", "D8", "D10")

# -----------------------------------------------------------------------------
# 1. Read the ordered gene set
# -----------------------------------------------------------------------------
gene_info <- read.delim(
  gene_file,
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "",
  fill = TRUE,
  stringsAsFactors = FALSE
)
if (ncol(gene_info) < 3L) {
  stop("Gene-list file must contain three columns: AGI, symbol, and category.")
}
gene_info <- gene_info[, 1:3, drop = FALSE]
gene_info[] <- lapply(gene_info, trimws)
gene_info <- gene_info[nzchar(gene_info$gene), , drop = FALSE]
gene_info$gene <- toupper(gene_info$gene)
gene_info <- gene_info[!duplicated(gene_info$gene), , drop = FALSE]
gene_info$gene_order <- seq_len(nrow(gene_info))

if (!nrow(gene_info)) stop("The gene-list file contains no valid genes.")
if (any(!nzchar(gene_info$symbol)) || any(!nzchar(gene_info$category))) {
  stop("Every selected gene must have a symbol and functional category.")
}
category_order <- unique(gene_info$category)

write.csv(
  gene_info,
  file.path(output_dir, "02_fig2c_7_20_gene_mapping.csv"),
  row.names = FALSE,
  quote = TRUE
)

# -----------------------------------------------------------------------------
# 2. Read and validate the genome-wide sample-aware, state-resolved DS/DD table
# -----------------------------------------------------------------------------
message("Reading genome-wide state-resolved DS/DD results...")
dsdd <- read.csv(
  if (grepl("\\.gz$", dsdd_file, ignore.case = TRUE)) gzfile(dsdd_file) else dsdd_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "day", "state", "gene", "strong_support", "signed_support",
  "DS_q_context", "DD_q_context"
)
missing_columns <- setdiff(required_columns, colnames(dsdd))
if (length(missing_columns)) {
  stop("DS/DD table is missing columns: ", paste(missing_columns, collapse = ", "))
}

dsdd$gene <- toupper(trimws(dsdd$gene))
dsdd$day <- trimws(dsdd$day)
if (!is.logical(dsdd$strong_support)) {
  dsdd$strong_support <- toupper(trimws(as.character(dsdd$strong_support))) == "TRUE"
}
dsdd$signed_support <- as.numeric(dsdd$signed_support)

selected <- dsdd[
  dsdd$gene %in% gene_info$gene & dsdd$day %in% days,
  ,
  drop = FALSE
]
if (!nrow(selected)) stop("No selected genes and dates were found in the DS/DD table.")

missing_genes <- setdiff(gene_info$gene, unique(selected$gene))
if (length(missing_genes)) {
  stop("Selected genes absent from DS/DD results: ", paste(missing_genes, collapse = ", "))
}

# A strongly supported context retains its signed DS/DD score; all unsupported
# contexts contribute zero. Thus each date-level mean summarizes the full set of
# matched query-state contexts and is not restricted to significant contexts.
selected$signed_zero <- ifelse(
  selected$strong_support & is.finite(selected$signed_support),
  selected$signed_support,
  0
)

# -----------------------------------------------------------------------------
# 3. Calculate the mean signed DS/DD score across matched states
# -----------------------------------------------------------------------------
mean_score <- aggregate(
  signed_zero ~ gene + day,
  selected,
  FUN = function(value) mean(value, na.rm = TRUE)
)
names(mean_score)[3] <- "mean_signed_score_across_states"

supported_count <- aggregate(
  strong_support ~ gene + day,
  selected,
  FUN = function(value) sum(value, na.rm = TRUE)
)
names(supported_count)[3] <- "n_strong_supported_states"

state_count <- aggregate(
  state ~ gene + day,
  selected,
  FUN = function(value) length(unique(value))
)
names(state_count)[3] <- "n_matched_states"

result <- Reduce(
  function(x, y) merge(x, y, by = c("gene", "day"), all = TRUE),
  list(mean_score, supported_count, state_count)
)

expected_grid <- expand.grid(
  gene = gene_info$gene,
  day = days,
  stringsAsFactors = FALSE
)
result <- merge(expected_grid, result, by = c("gene", "day"), all.x = TRUE)
result <- merge(result, gene_info, by = "gene", all.x = TRUE, sort = FALSE)
result <- result[order(result$gene_order, match(result$day, days)), , drop = FALSE]

if (anyNA(result$mean_signed_score_across_states)) {
  absent <- result[is.na(result$mean_signed_score_across_states), c("gene", "day")]
  stop(
    "Missing matched-state records for: ",
    paste(paste(absent$gene, absent$day, sep = "/"), collapse = ", ")
  )
}

write.csv(
  result,
  file.path(output_dir, "04_fig2c_7_20_date_mean_across_states_DS_DD.csv"),
  row.names = FALSE,
  quote = TRUE
)

# -----------------------------------------------------------------------------
# 4. Generate the publication-ready heatmap
# -----------------------------------------------------------------------------
plot_data <- result
plot_data$day <- factor(plot_data$day, levels = days)
plot_data$symbol <- factor(plot_data$symbol, levels = rev(gene_info$symbol))
plot_data$category <- factor(plot_data$category, levels = category_order)

palette <- c("#3B6FB6", "#A9BFDD", "#F3F3F3", "#F3A18F", "#D73027")

heatmap_plot <- ggplot(
  plot_data,
  aes(day, symbol, fill = mean_signed_score_across_states)
) +
  geom_tile(colour = "white", linewidth = 0.5) +
  facet_grid(
    category ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_fill_gradientn(
    colours = palette,
    limits = c(-2.5, 2.5),
    oob = squish,
    na.value = "#E7E7E7",
    name = "Mean signed\nDS/DD score",
    guide = guide_colourbar(nbin = 128, display = "rectangles")
  ) +
  labs(x = "Regeneration date") +
  theme_minimal(base_size = 11, base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 5)),
    axis.text.x = element_text(size = 11.5, face = "bold", colour = "black"),
    axis.text.y = element_text(size = 9.2, face = "italic", colour = "black"),
    panel.spacing.y = grid::unit(0.45, "lines"),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 8.8, hjust = 1),
    strip.background = element_blank(),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.margin = margin(4, 5, 3, 4),
    plot.background = element_rect(fill = "white", colour = NA)
  )

pdf_file <- file.path(
  output_dir,
  "06_fig2c_7_20_DS_DD_by_date_mean_states_heatmap.pdf"
)
png_file <- file.path(
  output_dir,
  "06_fig2c_7_20_DS_DD_by_date_mean_states_heatmap.png"
)

ggsave(
  pdf_file,
  heatmap_plot,
  width = 6.5,
  height = 5.2,
  units = "in",
  device = cairo_pdf
)
ggsave(
  png_file,
  heatmap_plot,
  width = 6.5,
  height = 5.2,
  units = "in",
  dpi = 450,
  bg = "white"
)

# -----------------------------------------------------------------------------
# 5. Record provenance and the software environment
# -----------------------------------------------------------------------------
writeLines(c(
  paste0("DS/DD input: ", normalizePath(dsdd_file)),
  paste0("DS/DD input MD5: ", unname(tools::md5sum(dsdd_file))),
  paste0("Gene list: ", normalizePath(gene_file)),
  paste0("Gene-list MD5: ", unname(tools::md5sum(gene_file))),
  paste0("Selected genes: ", nrow(gene_info)),
  paste0("Dates: ", paste(days, collapse = ", ")),
  paste0("Functional-category order: ", paste(category_order, collapse = "; ")),
  paste0("DS/DD rows used: ", nrow(selected)),
  paste0("Matched states per gene/date: ", paste(sort(unique(result$n_matched_states)), collapse = ", ")),
  paste0("Strong-supported states across gene/date combinations: ",
         paste(sort(unique(result$n_strong_supported_states)), collapse = ", ")),
  "Summary statistic: mean signed_support across all matched query-state contexts; unsupported contexts are set to zero.",
  "Direction: positive/red indicates mutant-high evidence; negative/blue indicates WT-high or reduced-in-mutant evidence.",
  "Color limits: -2.5 to 2.5; values beyond the limits are squished.",
  "Output size: 6.5 x 5.2 inches; PNG resolution: 450 dpi."
), file.path(output_dir, "07_fig2c_7_20_mean_states_analysis_information.txt"))

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "08_fig2c_7_20_mean_states_sessionInfo.txt")
)

message("Completed reproducible mean-across-states DS/DD heatmap: ", png_file)
