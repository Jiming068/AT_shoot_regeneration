# Portable path configuration (override with environment variables).
project_dir <- Sys.getenv("SHOOT_REGEN_PROJECT_DIR", unset = ".")
output_root <- Sys.getenv("SHOOT_REGEN_OUTPUT_DIR", unset = file.path(project_dir, "02_results"))
external_data_dir <- Sys.getenv("SHOOT_REGEN_EXTERNAL_DATA_DIR", unset = file.path(project_dir, "external_data"))
work_dir <- Sys.getenv("SHOOT_REGEN_WORK_DIR", unset = ".")

# 2022-09-05
# correlation analysis

dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)

library(tidyverse)
library(Seurat)
library(patchwork)
set.seed(128)  #设置随机数种子，使结果可重复

library(pheatmap)

snRNA2 <- readRDS(file.path(external_data_dir, "01_sct_dim30_res1.2_d4-16_rm12_20230207.rds"))

# merged_meta <- snRNA2@meta.data
#write.csv(merged_meta, 'merged_meta_data_20220904.csv', row.names = F)

table(snRNA2$orig.ident)

Idents(object = snRNA2) <- "seurat_clusters"
Idents(object = snRNA2)

av <- AverageExpression(snRNA2, group.by = 'seurat_clusters', assays = 'RNA', slot = 'counts')

av = av[[1]]

# 选出标准差最大的1000个基因
cg <- names(tail(sort(apply(av, 1, sd)), 1000))

# 查看这1000个基因在各个细胞中的表达矩阵
view(av[cg,])

# 查看细胞群的相关性矩阵
view(cor(av[cg,], method = 'spearman'))

# pheatmap 绘制热图
pdf("Fig1_pheatmap_correlation_clusters_RNA_counts_spearman.pdf", width = 12, height = 12)
pheatmap::pheatmap(cor(av[cg,], method = 'spearman'), cluster_cols = TRUE, cluster_rows = TRUE)
dev.off()

pdf("Fig1_pheatmap_correlation_clusters_RNA_counts_pearson.pdf", width = 12, height = 12)
pheatmap::pheatmap(cor(av[cg,], method = 'pearson'), cluster_cols = TRUE, cluster_rows = TRUE, fontsize = 15)
dev.off()

pdf("Fig1_pheatmap_correlation_clusters_RNA_counts_pearson.pdf", width = 12, height = 12)
pheatmap::pheatmap(cor(av[cg,], method = 'pearson'), cluster_cols = TRUE, cluster_rows = TRUE, fontsize = 15)
dev.off()

av <- AverageExpression(snRNA2, group.by = 'seurat_clusters', assays = 'SCT', slot = 'counts')

av = av[[1]]

# 选出标准差最大的1000个基因
cg <- names(tail(sort(apply(av, 1, sd)), 1000))

# 查看这1000个基因在各个细胞中的表达矩阵
view(av[cg,])

# 查看细胞群的相关性矩阵
view(cor(av[cg,], method = 'spearman'))

# pheatmap 绘制热图
pdf("pheatmap_correlation_clusters_SCT_counts_spearman.pdf", width = 12, height = 12)
pheatmap::pheatmap(cor(av[cg,], method = 'spearman'), cluster_cols = TRUE, cluster_rows = TRUE)
dev.off()

pdf("pheatmap_correlation_clusters_SCT_counts_pearson.pdf", width = 12, height = 12)
pheatmap::pheatmap(cor(av[cg,], method = 'pearson'), cluster_cols = TRUE, cluster_rows = TRUE, fontsize = 15)
dev.off()
