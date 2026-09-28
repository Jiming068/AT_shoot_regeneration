# Portable path configuration (override with environment variables).
project_dir <- Sys.getenv("SHOOT_REGEN_PROJECT_DIR", unset = ".")
output_root <- Sys.getenv("SHOOT_REGEN_OUTPUT_DIR", unset = file.path(project_dir, "02_results"))
external_data_dir <- Sys.getenv("SHOOT_REGEN_EXTERNAL_DATA_DIR", unset = file.path(project_dir, "external_data"))
work_dir <- Sys.getenv("SHOOT_REGEN_WORK_DIR", unset = ".")

dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)

library(Seurat)
library(patchwork)
library(pheatmap)
library(tidyverse)

set.seed(128)  #设置随机数种子，使结果可重复

# options(future.globals.maxSize = 5 * 1024^3)

snRNA2 <- readRDS(file.path(project_dir, "01_obj/01_WT_snRNA/01_snRNA_res1.2_SIM4-SIM16.rds"))
Idents(snRNA2) <- 'merged_cluster'

# subset 芽再生相关细胞类群

selected <- c('C1', 'C5', 'C8', 'C12', 'C18', 'C21')
subset <- subset(x = snRNA2, idents = selected, invert = TRUE)

DefaultAssay(subset) <- "RNA"

subset <- SCTransform(subset, vars.to.regress =c("percent.mt", 'percent.cp'), verbose = FALSE)

subset <- RunPCA(subset, verbose = FALSE)

ep <- ElbowPlot(subset, ndims = 50)
pdf("elbowplot.pdf", width = 4, height = 3)
print(ep)
dev.off()

subset <- FindNeighbors(subset, dims = 1:19, verbose = FALSE)
subset <- FindClusters(subset, verbose = FALSE, resolution = 0.2)

subset <- RunUMAP(subset, reduction = "pca", dims = 1:19)

subset$date <- factor(subset$date, levels =c('SIM4', 'SIM6', 'SIM8', 'SIM10', 'SIM14', 'SIM16'))

subset$orig.ident <- factor(subset$orig.ident, levels = c('SIM4_1', 'SIM4_2', 'SIM4_3', 'SIM6_1', 'SIM6_2', 'SIM6_3',
                 'SIM8_1', 'SIM8_2', 'SIM8_3', 'SIM10_1', 'SIM10_2', 'SIM10_3', 'SIM14_1', 'SIM14_2', 'SIM14_3',
                 'SIM16_1', 'SIM16_2', 'SIM16_3'))

plot1 <- DimPlot(subset, reduction = "umap", pt.size = 0.01, label=T, ncol = 3)+
    theme_bw()+
    theme(panel.grid =element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
pdf("01_recluster_sct_UMAP_dim19_res0.2_1.pdf", width = 4, height = 3.5)
print(plot1)
dev.off()

plot2 <- DimPlot(subset, reduction = "umap", pt.size = 0.01, split.by = 'date', label=F, ncol = 3)+
    theme_bw()+
    theme(panel.grid =element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
pdf("02_recluster_UMAP_dim19_res0.2_splitdate_1.pdf", width = 6, height = 4)
print(plot2)
dev.off()

# 删除列
subset@meta.data <- subset@meta.data[, !colnames(subset@meta.data) %in% c('SCT_snn_res.0.25','SCT_snn_res.1.2')]

# 验证删除结果

saveRDS(subset, file = "snRNA_subset_sct_dim19_res0.2.rds")

# 重新配色

subset <- readRDS(file.path(work_dir, "snRNA_subset_sct_dim19_res0.2.rds"))

Idents(object = subset) <- "seurat_clusters"

data_markers <- FindAllMarkers(subset, only.pos = TRUE, min.pct = 0, logfc.threshold = 0.1)
write.table(data_markers, file='recluster_DEGs_dim19_res0.2_20251225.txt', col.names = TRUE, sep = '\t')
write.csv(data_markers, file='recluster_DEGs_dim19_res0.2_20251225.csv', row.names = FALSE)

# DEGs for selected clusters

cluster0vs1_markers <- FindMarkers(subset,  ident.1 = 0, ident.2 = 1, min.pct = 0, logfc.threshold = 0.1, only.pos = TRUE)
write.table(cluster0vs1_markers, file='recluster_DEGs_dim20_res0.2_0vs1.txt', col.names = TRUE, sep = '\t')
write.csv(cluster0vs1_markers, file='recluster_DEGs_dim20_res0.2_0vs1.csv', row.names = TRUE)

cluster1vs0_markers <- FindMarkers(subset,  ident.1 = 1, ident.2 = 0, min.pct = 0, logfc.threshold = 0.1, only.pos = TRUE)
write.table(cluster1vs0_markers, file='recluster_DEGs_dim20_res0.2_1vs0.txt', col.names = TRUE, sep = '\t')
write.csv(cluster1vs0_markers, file='recluster_DEGs_dim20_res0.2_1vs0.csv', row.names = TRUE)

subset <- readRDS('snRNA_subset_sct_dim19_res0.2.rds')

selected <- c('1', '7', '9', '11')
subset_1 <- subset(x = subset, idents = selected, invert = TRUE)

plot1 <- DimPlot(subset_1, reduction = "umap", pt.size = 0.01, label=T, ncol = 3)+
    theme_bw()+
    theme(panel.grid =element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
pdf("01_recluster_subset_UMAP_dim19_res0.2_1.pdf", width = 4, height = 3.5)
print(plot1)
dev.off()

plot2 <- DimPlot(subset_1, reduction = "umap", pt.size = 0.01, split.by = 'date', label=F, ncol = 3)+
    theme_bw()+
    theme(panel.grid =element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
pdf("02_recluster_subset_UMAP_dim19_res0.2_splitdate_1.pdf", width = 6, height = 4)
print(plot2)
dev.off()
