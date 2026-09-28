#!/usr/bin/env Rscript

# Build the exact source matrix for the WT-versus-myb31 cluster-composition heatmap.
# The mutant denominator includes only query cells classified as in-reference.

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1]], fixed = TRUE)
}

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()

wt_file <- get_arg("wt-metadata", file.path(script_dir, "WT_cluster_date_metadata.csv.gz"))
mutant_file <- get_arg("mutant-metadata", file.path(script_dir, "mutant_mapping_for_heatmap.csv.gz"))
output_file <- get_arg(
  "output",
  file.path(script_dir, "cluster_comparison_heatmap_data_reproduced.csv")
)

wt_dates <- c("D4-root", "D6-root", "D8-root", "D10-root")
mutant_dates <- c("myb_SIM4", "myb_SIM6", "myb_SIM8", "myb_SIM10")
cluster_order <- as.character(0:11)

wt <- read.csv(
  wt_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = c(cell_id = "character", date = "character", seurat_clusters = "character")
)
mutant <- read.csv(
  mutant_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = c(
    cell_id = "character", date = "character", WT_pred_cluster = "character",
    WT_reference_status = "character"
  )
)

required_wt <- c("cell_id", "date", "seurat_clusters")
required_mutant <- c("cell_id", "date", "WT_pred_cluster", "WT_reference_status")
if (length(setdiff(required_wt, names(wt)))) stop("WT metadata columns are incomplete.")
if (length(setdiff(required_mutant, names(mutant)))) stop("Mutant metadata columns are incomplete.")
if (anyDuplicated(wt$cell_id) || anyDuplicated(mutant$cell_id)) stop("Cell IDs must be unique.")

complete_composition <- function(time, cluster, time_levels, cluster_levels, source_label) {
  observed <- aggregate(
    rep.int(1L, length(time)),
    by = list(dataset_time = as.character(time), cluster = as.character(cluster)),
    FUN = sum
  )
  names(observed)[3] <- "Freq"
  grid <- expand.grid(
    dataset_time = time_levels,
    cluster = cluster_levels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  result <- merge(grid, observed, by = c("dataset_time", "cluster"), all.x = TRUE, sort = FALSE)
  result$Freq[is.na(result$Freq)] <- 0L
  result$dataset_time <- factor(result$dataset_time, levels = time_levels)
  result$cluster <- factor(result$cluster, levels = cluster_levels)
  result <- result[order(result$cluster, result$dataset_time), , drop = FALSE]
  result$source <- source_label
  result$proportion <- ave(
    result$Freq,
    result$dataset_time,
    FUN = function(z) if (sum(z) > 0) z / sum(z) else rep(0, length(z))
  )
  result$display_time <- as.character(result$dataset_time)
  result$label <- ifelse(result$proportion >= 0.01, sprintf("%.1f", 100 * result$proportion), "")
  result
}

wt_keep <- wt$date %in% wt_dates & wt$seurat_clusters %in% cluster_order
wt_comp <- complete_composition(
  wt$date[wt_keep], wt$seurat_clusters[wt_keep], wt_dates, cluster_order, "WT observed"
)

mutant_keep <-
  mutant$date %in% mutant_dates &
  mutant$WT_pred_cluster %in% cluster_order &
  mutant$WT_reference_status == "in_reference"
mutant_comp <- complete_composition(
  mutant$date[mutant_keep], mutant$WT_pred_cluster[mutant_keep],
  mutant_dates, cluster_order, "myb31 predicted (in-reference)"
)

comparison <- rbind(wt_comp, mutant_comp)
comparison <- comparison[, c(
  "dataset_time", "cluster", "Freq", "source", "proportion", "display_time", "label"
)]
write.csv(comparison, output_file, row.names = FALSE, quote = TRUE)

sample_sizes <- rbind(
  aggregate(Freq ~ dataset_time + source, wt_comp, sum),
  aggregate(Freq ~ dataset_time + source, mutant_comp, sum)
)
write.csv(
  sample_sizes,
  file.path(dirname(output_file), "cluster_comparison_denominators_by_time.csv"),
  row.names = FALSE,
  quote = TRUE
)

checks <- aggregate(proportion ~ dataset_time + source, comparison, sum)
if (!all(abs(checks$proportion - 1) < 1e-12)) stop("Within-time proportions do not sum to one.")
if (nrow(comparison) != 96L) stop("Expected 96 rows (8 time points x 12 clusters).")

message("Saved heatmap source data: ", output_file)
message("WT cells included: ", sum(wt_keep))
message("Mutant in-reference cells included: ", sum(mutant_keep))
