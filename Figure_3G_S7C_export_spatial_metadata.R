#!/usr/bin/env Rscript

# Portable path configuration (override with environment variables).
project_dir <- Sys.getenv("SHOOT_REGEN_PROJECT_DIR", unset = ".")
output_root <- Sys.getenv("SHOOT_REGEN_OUTPUT_DIR", unset = file.path(project_dir, "02_results"))
external_data_dir <- Sys.getenv("SHOOT_REGEN_EXTERNAL_DATA_DIR", unset = file.path(project_dir, "external_data"))
work_dir <- Sys.getenv("SHOOT_REGEN_WORK_DIR", unset = ".")

suppressPackageStartupMessages(library(Seurat))

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1]], fixed = TRUE)
}

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()

input_rds <- get_arg(
  "input-rds",
  file.path(project_dir, "01_obj/02_WT_stereo/Stereo-seq_SIM6-SIM10_final.rds")
)
output_csv <- get_arg("output", file.path(script_dir, "00_spatial_metadata.csv.gz"))

obj <- readRDS(input_rds)
required_meta <- c("time", "cluster", "tissue")
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta)) stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "))
if (!"spatial" %in% names(obj@reductions)) stop("The Seurat object has no 'spatial' reduction.")

xy <- Embeddings(obj, "spatial")
if (ncol(xy) < 2L) stop("The spatial reduction must contain at least two coordinates.")

metadata <- data.frame(
  cell_id = Cells(obj),
  time = as.character(obj$time),
  cluster = as.character(obj$cluster),
  tissue = as.character(obj$tissue),
  x = as.numeric(xy[, 1]),
  y = as.numeric(xy[, 2]),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
con <- gzfile(output_csv, open = "wt")
write.csv(metadata, con, row.names = FALSE, quote = TRUE)
close(con)

info_file <- file.path(dirname(output_csv), "00_input_export_information.txt")
writeLines(c(
  paste0("Input RDS: ", normalizePath(input_rds)),
  paste0("Input MD5: ", unname(tools::md5sum(input_rds))),
  paste0("Exported cells: ", nrow(metadata)),
  paste0("SIM6 cells: ", sum(metadata$time == "SIM6")),
  paste0("Output: ", normalizePath(output_csv, mustWork = FALSE)),
  "Coordinates: first two dimensions of the Seurat 'spatial' reduction"
), info_file)

message("Saved metadata: ", output_csv)
