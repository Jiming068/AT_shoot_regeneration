# Portable path configuration (override with environment variables).
project_dir <- Sys.getenv("SHOOT_REGEN_PROJECT_DIR", unset = ".")
output_root <- Sys.getenv("SHOOT_REGEN_OUTPUT_DIR", unset = file.path(project_dir, "02_results"))
external_data_dir <- Sys.getenv("SHOOT_REGEN_EXTERNAL_DATA_DIR", unset = file.path(project_dir, "external_data"))
work_dir <- Sys.getenv("SHOOT_REGEN_WORK_DIR", unset = ".")

suppressPackageStartupMessages(library(ggplot2))

input_dir <- file.path(project_dir, "02_results/143_WT_C02418_literature_S_G2M_gene_expression")
outdir <- Sys.getenv("DNA_S_G2M_OUTDIR", file.path(output_root, "143_DNAreplication_S_G2Mentry_heatmap"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

clusters <- c("0", "2", "4", "1", "8")
modules <- c("DNA_replication_S", "G2_M_entry")
module_labels <- c("DNA_replication_S" = "DNA replication/S", "G2_M_entry" = "G2/M entry")
excluded_symbols <- c(
  "CYCA3;1", "MCM7", "CDC6/CDC6A", "CDKB1;2",
  "MCM2", "CDT1B", "PCNA2"
)

dat <- read.csv(file.path(input_dir, "03_literature_gene_expression_and_Zscore.csv"),
                stringsAsFactors = FALSE, check.names = FALSE)
dat <- dat[dat$Module %in% modules, , drop = FALSE]
dat <- dat[!dat$Gene_symbol %in% excluded_symbols, , drop = FALSE]
dat$Module <- factor(dat$Module, levels = modules)
dat <- dat[order(dat$Module, match(dat$peak_state, clusters), -dat$dynamic_range), , drop = FALSE]
dat$display_gene <- make.unique(dat$display_gene)

z <- as.matrix(dat[, paste0("Z_", clusters), drop = FALSE])
rownames(z) <- dat$display_gene
colnames(z) <- clusters
heat <- as.data.frame(as.table(z), stringsAsFactors = FALSE)
heat$subcluster <- factor(heat$subcluster, levels = clusters)
heat$module <- factor(dat$Module[match(as.character(heat$gene), dat$display_gene)], levels = modules)
heat$gene <- factor(heat$gene, levels = rev(dat$display_gene))

write.csv(dat, file.path(outdir, "01_DNAreplication_S_G2Mentry_refined_gene_set.csv"), row.names = FALSE)
write.csv(heat, file.path(outdir, "02_DNAreplication_S_G2Mentry_refined_heatmap_data.csv"), row.names = FALSE)

p <- ggplot(heat, aes(x = subcluster, y = gene, fill = Zscore)) +
  geom_tile(color = "white", linewidth = 0.40) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y", switch = "y",
             labeller = as_labeller(module_labels)) +
  scale_fill_gradient2(low = "#89A9CE", mid = "#F5F5F5", high = "#DD776E",
                       limits = c(-2.5, 2.5), name = "Gene-wise\nZ-score") +
  labs(x = "Trajectory state", y = NULL) +
  theme_classic(base_family = "Arial", base_size = 9.5) +
  theme(
    axis.text.x = element_text(face = "bold", color = "#222222", size = 9),
    axis.text.y = element_text(face = "italic", color = "#222222", size = 8.2),
    axis.ticks = element_blank(),
    axis.title.x = element_text(face = "bold", size = 9.5),
    strip.text.y = element_text(face = "bold", size = 8.6),
    strip.background = element_rect(fill = "#F0F0F0", color = NA),
    panel.spacing.y = grid::unit(1.0, "mm"),
    legend.title = element_text(size = 8.5),
    legend.text = element_text(size = 8),
    plot.margin = margin(4, 5, 4, 4)
  )

pdf_file <- file.path(outdir, "03_DNAreplication_S_G2Mentry_refined_genes_expression_heatmap.pdf")
png_file <- file.path(outdir, "03_DNAreplication_S_G2Mentry_refined_genes_expression_heatmap.png")
ggsave(pdf_file, p, width = 4.8, height = 5.2, device = cairo_pdf, bg = "white")
ggsave(png_file, p, width = 4.8, height = 5.2, dpi = 500, bg = "white")

writeLines(c(
  "The heatmap retains only literature-curated DNA replication/S and G2/M entry genes.",
  paste0("Excluded genes: ", paste(excluded_symbols, collapse = ", "), "."),
  "Expression was averaged across WT subclusters in the trajectory order 0 -> 2 -> 4 -> 1 -> 8.",
  "Rows show gene-wise Z-scores calculated from these cluster-average RNA expression values; values are capped at [-2.5, 2.5]."
), file.path(outdir, "04_method_information.txt"))
