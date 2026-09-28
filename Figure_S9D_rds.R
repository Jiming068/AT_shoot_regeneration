#!/usr/bin/env Rscript

# Converted from Figure_S9D_rds.ipynb.
# Notebook cell order is preserved; Markdown cells are retained as comments.

# ---- Code cell 1 ----
# Portable path configuration (override with environment variables).
project_dir <- Sys.getenv("SHOOT_REGEN_PROJECT_DIR", unset = ".")
output_root <- Sys.getenv("SHOOT_REGEN_OUTPUT_DIR", unset = file.path(project_dir, "02_results"))
external_data_dir <- Sys.getenv("SHOOT_REGEN_EXTERNAL_DATA_DIR", unset = file.path(project_dir, "external_data"))
work_dir <- Sys.getenv("SHOOT_REGEN_WORK_DIR", unset = ".")

# ---- Code cell 2 ----
library(Seurat)
library(patchwork)
set.seed(128)  #设置随机数种子，使结果可重复
library(reshape2)
#library(SeuratData)
library(tidyverse)
library(dplyr)
library(SeuratDisk)
library(harmony)
library(matrixStats)
packageVersion("matrixStats")

# ---- Code cell 3 ----
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)

# ---- Code cell 4 ----
snRNA2 <- readRDS("/data/work/05_cluster/01_all_sct_dim30_res0.5.rds")

# ---- Code cell 5 ----
subset <- subset(x=snRNA2, idents = c("0", '1', "2", "3", "4", '5', '7'), invert = TRUE)

# ---- Code cell 6 ----
subset <- SCTransform(subset, vars.to.regress = "percent.mt", verbose = FALSE)
subset <- RunPCA(subset, npcs = 50, verbose = FALSE)

ep <- ElbowPlot(subset, ndims = 50)
pdf("elbowplot.pdf", width = 4, height = 3)
print(ep)
dev.off()

# ---- Code cell 7 ----
saveRDS(object = subset, file = "01_WT_myb31_recluster_sct.rds")

# ---- Code cell 8 ----
subset <- readRDS("01_WT_myb31_recluster_sct.rds")

# ---- Code cell 9 ----
subset <- RunPCA(subset, npcs = 50, verbose = FALSE)

# ---- Code cell 10 ----
# Harmony整合（注意指定assay）
subset <- RunHarmony(subset,
                           group.by.vars = "batch",  # 批次变量名
                           theta = 2,                # 调整批次校正强度（默认2）
                           lambda = 2,
                           assay.use = "SCT")

# ---- Code cell 11 ----
#subset <- FindNeighbors(subset, dims = 1:30, reduction = "harmony", verbose = FALSE)
subset <- FindClusters(subset, verbose = FALSE, resolution = 0.8)

subset <- RunUMAP(subset, reduction = "harmony", dims = 1:30)

subset$date <- factor(subset$date, levels =c('WT_SIM4','WT_SIM6', 'WT_SIM8', 'WT_SIM10', 'myb_SIM4', 'myb_SIM6', 'myb_SIM8', 'myb_SIM10'))

# ---- Code cell 12 ----
# Visualization

plot1 <- DimPlot(subset, reduction = "umap", pt.size = 0.01, label=T)+
    theme_bw()+
    theme(panel.grid =element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
pdf("01_recluster_sct_UMAP_dim30_res0.8_1.pdf", width = 6, height = 5)
print(plot1)
dev.off()

plot2 <- DimPlot(subset, reduction = "umap", pt.size = 0.01, split.by = 'date', label=F, ncol = 4)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
pdf("02_recluster_UMAP_dim30_res0.8_splitbydate_1.pdf", width = 12, height = 6)
print(plot2)
dev.off()

plot3 <- DimPlot(subset, reduction = "umap", pt.size = 0.01, split.by = 'batch', label=T, ncol = 2)+
    theme_bw()+
    theme(panel.grid =element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
pdf("02_recluster_UMAP_dim30_res0.8_splitbybatch_1.pdf", width = 10, height = 4)
print(plot3)
dev.off()

# ---- Code cell 13 ----
select_markers <- read.table('/data/work/fig2c_genes_selected_7_7.txt', header = F)

labels <- select_markers$V2

gene <- select_markers$V1  # 修改这里：使用第二列作为基因名称

#绘制气泡图
dp <- DotPlot(
    subset,
    features = gene,  # 现在使用第二列的基因符号
    assay = 'SCT',
    cols = c("#E8F396", "#9E0142"),
    col.min = 0,
    col.max = 1,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    group.by = "seurat_clusters",
    split.by = NULL,
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = labels)  # 取消注释并添加这行来设置横坐标标签

# 保存气泡图到 PDF 文件
pdf("subset_markers_dotplot_cluster_7.pdf", width = 7, height = 5)
print(dp)
dev.off()

# ---- Code cell 14 ----
saveRDS(object = subset, file = "01_all_subset_sct_dim30_res0.2.rds")

# ---- Code cell 15 ----
Idents(object = subset) <- "batch"

# ---- Code cell 16 ----
WT <- subset(x=subset, idents = 'WT', invert = FALSE)
saveRDS(object = WT, file = "01_WT_sct_dim30_res0.3.rds")

# ---- Code cell 17 ----
mutant <- subset(x=subset, idents = 'Mutant', invert = FALSE)
saveRDS(object = mutant, file = "01_mutant_sct_dim30_res0.3.rds")

# ---- Code cell 18 ----
# Visualization

plot1 <- DimPlot(subset, reduction = "umap", pt.size = 0.01, label=T)+
    theme_bw()+
    theme(panel.grid =element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
pdf("01_recluster_sct_UMAP_dim29_res0.3_1.pdf", width = 6, height = 5)
print(plot1)
dev.off()

plot2 <- DimPlot(subset, reduction = "umap", pt.size = 0.01, split.by = 'date', label=F, ncol = 4)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
pdf("02_recluster_UMAP_dim29_res0.3_splitbydate_1.pdf", width = 12, height = 6)
print(plot2)
dev.off()

plot3 <- DimPlot(subset, reduction = "umap", pt.size = 0.01, split.by = 'batch', label=T, ncol = 2)+
    theme_bw()+
    theme(panel.grid =element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
pdf("02_recluster_UMAP_dim29_res0.3_splitbybatch_1.pdf", width = 10, height = 4)
print(plot3)
dev.off()

# ---- Code cell 19 ----
select_markers <- read.table('/data/work/fig2c_genes_selected_7_7.txt', header = F)

labels <- select_markers$V2

gene <- select_markers$V1  # 修改这里：使用第二列作为基因名称

#绘制气泡图
dp <- DotPlot(
    subset,
    features = gene,  # 现在使用第二列的基因符号
    assay = 'SCT',
    cols = c("#E8F396", "#9E0142"),
    col.min = 0,
    col.max = 1,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    group.by = "seurat_clusters",
    split.by = NULL,
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = labels)  # 取消注释并添加这行来设置横坐标标签

# 保存气泡图到 PDF 文件
pdf("subset_markers_dotplot_cluster_7.pdf", width = 7, height = 5)
print(dp)
dev.off()

# ---- Code cell 20 ----
select_markers <- read.table('/data/work/fig2c_genes_selected_7_7.txt', header = F)

labels <- select_markers$V2

gene <- select_markers$V1  # 修改这里：使用第二列作为基因名称

#绘制气泡图
dp <- DotPlot(
    WT,
    features = gene,  # 现在使用第二列的基因符号
    assay = 'SCT',
    cols = c("#E8F396", "#9E0142"),
    col.min = 0,
    col.max = 1,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    group.by = "seurat_clusters",
    split.by = NULL,
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = labels)  # 取消注释并添加这行来设置横坐标标签

# 保存气泡图到 PDF 文件
pdf("WT_markers_dotplot_cluster_7.pdf", width = 7, height = 5)
print(dp)
dev.off()

# ---- Code cell 21 ----
select_markers <- read.table('/data/work/fig2c_genes_selected_7_1.txt', header = F)

labels <- select_markers$V2

gene <- select_markers$V1  # 修改这里：使用第二列作为基因名称

#绘制气泡图
dp <- DotPlot(
    mutant,
    features = gene,  # 现在使用第二列的基因符号
    assay = 'SCT',
    cols = c("#E8F396", "#9E0142"),
    col.min = 0,
    col.max = 1,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    group.by = "seurat_clusters",
    split.by = NULL,
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = labels)  # 取消注释并添加这行来设置横坐标标签

# 保存气泡图到 PDF 文件
pdf("mutant_markers_dotplot_cluster_1.pdf", width = 7, height = 5)
print(dp)
dev.off()

# ---- Code cell 22 ----
Idents(object = subset) <- "seurat_clusters"

# How does cluster membership vary by replicate?
write.csv(x, file = "cells_eachCluster_per_date.csv")

x <- prop.table(table(Idents(subset), subset$date), margin = 2)
write.csv(x, file = "prop_cells_eachCluster_per_date.csv")

# 加载必要的包
library(ggplot2)
library(reshape2)
library(dplyr)
library(viridis)

# 读取数据
data <- read.csv("prop_cells_eachCluster_per_date.csv", header = TRUE)

# 查看数据结构
str(data)

# 重命名第一列为细胞类群
colnames(data)[1] <- "CellType"

# 转换数据为长格式
data_long <- melt(data, id.vars = "CellType",
                  variable.name = "Sample",
                  value.name = "Percentage")

# 确保细胞类群是因子类型，保持原始顺序
data_long$CellType <- factor(data_long$CellType, levels = unique(data_long$CellType))

# 绘制热图风格柱状图
# 使用白色到深蓝的渐变颜色方案
plot <- ggplot(data_long, aes(x = Sample, y = CellType, fill = Percentage)) +
  geom_tile(color = "white", size = 1) +
  geom_text(aes(label = sprintf("%.1f%%", Percentage*100)),
            color = "black", size = 2.5, fontface = "bold") +
  scale_fill_gradient(
    low = "white",
    high = "darkblue",
    name = "Percentage",
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    title = "Cell proportion",
    subtitle = "",
    x = "",
    y = "Cluster"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
    plot.subtitle = element_text(hjust = 0.5, size = 8, color = "gray40"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    panel.grid = element_blank()
  )

# 保存图像
pdf("cell_number_change_date_1.pdf", width = 4, height = 3.5)
print(plot)
dev.off()
# 如果需要保留原始dims=1:30的结果，可以重新运行一次
cat("All plots generated successfully!\n")

# ---- Code cell 23 ----
dims_list <- list(
  c(1:15),
  c(1:17),
  c(1:19),
  c(1:22),
  c(1:25),
  c(1:35)
)

# 设置因子水平（只需要设置一次，因为subset对象会被复用）

subset$date <- factor(subset$date, levels =c('WT_SIM4','WT_SIM6', 'WT_SIM8', 'WT_SIM10', 'myb_SIM4', 'myb_SIM6', 'myb_SIM8', 'myb_SIM10'))

# 循环针对不同dims生成UMAP图
for (dims in dims_list) {
  dim_max <- max(dims)
  dim_str <- paste0("dim", dim_max)

  cat("Processing", dim_str, "...\n")

  # 使用临时变量存储结果，避免覆盖原始对象
  temp_subset <- subset

  # 执行降维和聚类
  temp_subset <- FindNeighbors(temp_subset, dims = dims, reduction = "harmony", verbose = FALSE)
  temp_subset <- FindClusters(temp_subset, verbose = FALSE, resolution = 0.5)
  temp_subset <- RunUMAP(temp_subset, reduction = "harmony", dims = dims)

  # 创建主题设置函数，避免重复代码
  create_theme <- function() {
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    )
  }

  # 生成第一个图（聚类图）
  plot1 <- DimPlot(temp_subset, reduction = "umap", pt.size = 0.01, label = TRUE, ncol = 3) +
    create_theme()

  pdf_filename1 <- paste0("01_sct_UMAP_", dim_str, "_res0.5_1.pdf")
  pdf(pdf_filename1, width = 4, height = 3.5)
  print(plot1)
  dev.off()

  # 生成第二个图（按date分面）
  plot2 <- DimPlot(temp_subset, reduction = "umap", pt.size = 0.01,
                   split.by = 'date', label = FALSE, ncol = 4) +
    create_theme()

  pdf_filename2 <- paste0("02_UMAP_", dim_str, "_res0.5_splitdate_1.pdf")
  pdf(pdf_filename2, width = 8, height = 4)
  print(plot2)
  dev.off()

  # 生成第三个图（按batch分面）
  plot3 <- DimPlot(temp_subset, reduction = "umap", pt.size = 0.01,
                   split.by = 'batch', label = FALSE, ncol = 4) +
    create_theme()

  pdf_filename3 <- paste0("03_UMAP_", dim_str, "_res0.5_splitbatch_1.pdf")
  pdf(pdf_filename3, width = 8, height = 4)
  print(plot3)
  dev.off()

  cat("Saved plots for", dim_str, "\n\n")
}

# ---- Code cell 24 ----
# subcluster 14
snRNA2 <- readRDS("/data/work/05_cluster/01_recluster/01_all_subset_sct_dim30_res0.3.rds")

# ---- Code cell 25 ----
subset <- snRNA2

# ---- Code cell 26 ----
#subset <- FindNeighbors(subset, dims = 1:30, reduction = "harmony", verbose = FALSE)
subset <- FindClusters(subset, verbose = FALSE, resolution = 0.3)

subset <- RunUMAP(subset, reduction = "harmony", dims = 1:30)

subset$date <- factor(subset$date, levels =c('WT_SIM4','WT_SIM6', 'WT_SIM8', 'WT_SIM10', 'myb_SIM4', 'myb_SIM6', 'myb_SIM8', 'myb_SIM10'))

# Visualization

plot1 <- DimPlot(subset, reduction = "umap", pt.size = 0.05, label=T)+
    theme_bw()+
    theme(panel.grid =element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
pdf("01_recluster_sct_UMAP_theta2_lambda2_dim30_res0.3_1.pdf", width = 6, height = 5)
print(plot1)
dev.off()

plot2 <- DimPlot(subset, reduction = "umap", pt.size = 0.05, split.by = 'date', label=F, ncol = 4)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
pdf("02_recluster_sct_UMAP_theta2_lambda2_dim30_res0.3_splitbydate_1.pdf", width = 12, height = 6)
print(plot2)
dev.off()

plot3 <- DimPlot(subset, reduction = "umap", pt.size = 0.05, split.by = 'batch', label=T, ncol = 2)+
    theme_bw()+
    theme(panel.grid =element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
pdf("02_recluster_sct_UMAP_theta2_lambda2_dim30_res0.3_splitbybatch_1.pdf", width = 10, height = 4)
print(plot3)
dev.off()

# ---- Code cell 27 ----
seurat_obj <- snRNA2

# ---- Code cell 28 ----
# 查看当前分群
DimPlot(seurat_obj, label = TRUE, repel = TRUE)
table(seurat_obj$seurat_clusters)

# ---- Code cell 29 ----
# 提取目标类群
target_clusters <- c("14")  # 假设要对0,2,4簇进行亚群分析

sub_obj <- subset(seurat_obj, idents = target_clusters)

# 重置聚类标识
Idents(sub_obj) <- "original_clusters"  # 保存原始聚类信息

# ---- Code cell 30 ----
sub_obj <- SCTransform(sub_obj, vars.to.regress = "percent.mt", verbose = FALSE)
sub_obj <- RunPCA(sub_obj, npcs = 50, verbose = FALSE)

ep <- ElbowPlot(sub_obj, ndims = 50)
pdf("elbowplot_14.pdf", width = 4, height = 3)
print(ep)
dev.off()

# ---- Code cell 31 ----
sub_obj <- RunPCA(sub_obj, npcs = 50, verbose = FALSE)
# Harmony整合（注意指定assay）
sub_obj <- RunHarmony(sub_obj,
                           group.by.vars = "batch",  # 批次变量名
                           theta = 2,                # 调整批次校正强度（默认2）
                           lambda = 2,
                           assay.use = "SCT")

# ---- Code cell 32 ----
sub_obj <- FindNeighbors(sub_obj, dims = 1:15, reduction = "harmony", verbose = FALSE)
sub_obj <- FindClusters(sub_obj, verbose = FALSE, resolution = 0.2)

sub_obj <- RunUMAP(sub_obj, reduction = "harmony", dims = 1:15)

sub_obj$date <- factor(sub_obj$date, levels =c('WT_SIM4','WT_SIM6', 'WT_SIM8', 'WT_SIM10', 'myb_SIM4', 'myb_SIM6', 'myb_SIM8', 'myb_SIM10'))

# ---- Code cell 33 ----
# Visualization

plot1 <- DimPlot(sub_obj, reduction = "umap", pt.size = 1, label=T)+
    theme_bw()+
    theme(panel.grid =element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
pdf("01_sub_obj_sct_UMAP_dim15_res0.2_1.pdf", width = 4, height = 4)
print(plot1)
dev.off()

plot2 <- DimPlot(sub_obj, reduction = "umap", pt.size = 1, split.by = 'date', label=F, ncol = 4)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
pdf("02_sub_obj_UMAP_dim15_res0.2_splitbydate_1.pdf", width = 8, height = 4)
print(plot2)
dev.off()

plot3 <- DimPlot(sub_obj, reduction = "umap", pt.size = 1, split.by = 'batch', label=T, ncol = 2)+
    theme_bw()+
    theme(panel.grid =element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
pdf("02_sub_obj_UMAP_dim15_res0.2_splitbybatch_1.pdf", width = 8, height = 4)
print(plot3)
dev.off()

# ---- Code cell 34 ----
# 将亚群信息合并回原对象
# 方法A：创建包含原始和亚群信息的综合标签
# 首先，为原对象创建一个新列，初始值为原始聚类
seurat_obj$combined_clusters <- as.character(Idents(seurat_obj))

# ---- Code cell 35 ----
# 为子集对象创建新标签
sub_obj$subcluster <- paste0("C", Idents(sub_obj))  # 亚群标签，例如"C0", "C1"

# ---- Code cell 36 ----
# 创建亚群的完整标签（包含原始聚类信息）
sub_obj$combined_label <- paste0(
  "Main_", sub_obj$SCT_snn_res.0.3,  # 原始聚类
  "_Sub_", sub_obj$subcluster          # 亚群信息
)

# ---- Code cell 37 ----
# 提取子集细胞的标签映射
cell_subcluster_map <- data.frame(
  cell = colnames(sub_obj),
  main_cluster = sub_obj$SCT_snn_res.0.3,
  subcluster = sub_obj$subcluster,
  combined_label = sub_obj$combined_label,
  stringsAsFactors = FALSE
)

# ---- Code cell 38 ----
# 将映射应用到原对象
seurat_obj$subcluster <- "Not_analyzed"  # 默认值
seurat_obj$combined_label <- as.character(Idents(seurat_obj))

# 更新原对象中的亚群信息
for (cell in cell_subcluster_map$cell) {
  seurat_obj$subcluster[colnames(seurat_obj) == cell] <-
    cell_subcluster_map$subcluster[cell_subcluster_map$cell == cell]

  seurat_obj$combined_label[colnames(seurat_obj) == cell] <-
    cell_subcluster_map$combined_label[cell_subcluster_map$cell == cell]
}

# 验证
table(seurat_obj$subcluster)
table(seurat_obj$combined_label)

# ---- Code cell 39 ----
# 查看当前分群
DimPlot(seurat_obj, label = TRUE,  group.by = 'combined_label', repel = TRUE)
table(seurat_obj$combined_label)

# ---- Code cell 40 ----
# Visualization

plot1 <- DimPlot(seurat_obj, reduction = "umap", group.by = 'combined_label', pt.size = 0.05, label=T)+
    theme_bw()+
    theme(panel.grid =element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
pdf("01_combined_sct_UMAP_1.pdf", width = 6, height = 5)
print(plot1)
dev.off()

plot2 <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.05, group.by = 'combined_label', split.by = 'date', label=F, ncol = 4)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
pdf("02_combined_UMAP_splitbydate_1.pdf", width = 12, height = 6)
print(plot2)
dev.off()

plot3 <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.05, group.by = 'combined_label', split.by = 'batch', label=T, ncol = 2)+
    theme_bw()+
    theme(panel.grid =element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
pdf("02_combined_UMAP_splitbybatch_1.pdf", width = 10, height = 4)
print(plot3)
dev.off()
