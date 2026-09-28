#!/usr/bin/env Rscript

# Converted from Figure_S9B_cluster.ipynb.
# Notebook cell order is preserved; Markdown cells are retained as comments.

# ---- Code cell 1 ----
# Portable path configuration (override with environment variables).
project_dir <- Sys.getenv("SHOOT_REGEN_PROJECT_DIR", unset = ".")
output_root <- Sys.getenv("SHOOT_REGEN_OUTPUT_DIR", unset = file.path(project_dir, "02_results"))
external_data_dir <- Sys.getenv("SHOOT_REGEN_EXTERNAL_DATA_DIR", unset = file.path(project_dir, "external_data"))
work_dir <- Sys.getenv("SHOOT_REGEN_WORK_DIR", unset = ".")

# ---- Code cell 2 ----
library(Seurat)
packageVersion("Seurat")
# ‘4.4.0’
library(patchwork)
set.seed(128)  #设置随机数种子，使结果可重复
library(reshape2)
#library(SeuratData)
library(dplyr)
library(SeuratDisk)
library(harmony)
library(matrixStats)
packageVersion("matrixStats")

# ---- Code cell 3 ----
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)

# ---- Code cell 4 ----
dir <- c(
        '/data/work/04_WT_matrix/WT_SIM4_1',
        '/data/work/04_WT_matrix/WT_SIM4_2',
        '/data/work/04_WT_matrix/WT_SIM4_3',
        '/data/work/04_WT_matrix/WT_SIM6_1',
        '/data/work/04_WT_matrix/WT_SIM6_2',
        '/data/work/04_WT_matrix/WT_SIM6_3',
        '/data/work/04_WT_matrix/WT_SIM8_1',
        '/data/work/04_WT_matrix/WT_SIM8_2',
        '/data/work/04_WT_matrix/WT_SIM8_3',
        '/data/work/04_WT_matrix/WT_SIM10_1',
        '/data/work/04_WT_matrix/WT_SIM10_2',
        '/data/work/04_WT_matrix/WT_SIM10_3',
        '/data/work/03_matrix/myb_SIM4_1',
        '/data/work/03_matrix/myb_SIM4_2',
        '/data/work/03_matrix/myb_SIM6_1',
        '/data/work/03_matrix/myb_SIM6_2',
        '/data/work/03_matrix/myb_SIM8_1',
        '/data/work/03_matrix/myb_SIM8_2',
        '/data/work/03_matrix/myb_SIM10_1',
        '/data/work/03_matrix/myb_SIM10_2'
        )

# ---- Code cell 5 ----
sample_name <- c('WT_SIM4_1', 'WT_SIM4_2', 'WT_SIM4_3', 'WT_SIM6_1', 'WT_SIM6_2', 'WT_SIM6_3', 'WT_SIM8_1', 'WT_SIM8_2', 'WT_SIM8_3', 'WT_SIM10_1', 'WT_SIM10_2', 'WT_SIM10_3',
                  'myb_SIM4_1', 'myb_SIM4_2', 'myb_SIM6_1', 'myb_SIM6_2',  'myb_SIM8_1', 'myb_SIM8_2', 'myb_SIM10_1', 'myb_SIM10_2')

date <- c('WT_SIM4', 'WT_SIM4', 'WT_SIM4', 'WT_SIM6', 'WT_SIM6', 'WT_SIM6', 'WT_SIM8', 'WT_SIM8', 'WT_SIM8', 'WT_SIM10', 'WT_SIM10', 'WT_SIM10',
          'myb_SIM4', 'myb_SIM4',  'myb_SIM6', 'myb_SIM6',  'myb_SIM8', 'myb_SIM8', 'myb_SIM10', 'myb_SIM10')

batch <- c('WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT',
          'Mutant', 'Mutant', 'Mutant', 'Mutant', 'Mutant', 'Mutant', 'Mutant', 'Mutant')

#以下代码会把每个样本的数据创建一个seurat对象，并存放到列表scRNAlist里
snRNAlist <- list()
for(i in 1:length(dir)){
counts <- Read10X(data.dir = dir[i], gene.column = 1)
snRNAlist[[i]] <- CreateSeuratObject(counts, project = sample_name[i], min.cells = 3, min.features = 200)
snRNAlist[[i]][['date']] <- date[i]
snRNAlist[[i]][['batch']] <- batch[i]
snRNAlist[[i]][['percent.mt']] <- PercentageFeatureSet(snRNAlist[[i]], pattern = '^ATM')
}

################################################################################################################################3
#使用merge函数将11个seurat对象合并成一个seurat对象
# scRNA1 <- readRDS("00_data_raw_20220412.rds")

snRNA1 <- merge(snRNAlist[[1]], y=c(snRNAlist[[2]], snRNAlist[[3]],snRNAlist[[4]], snRNAlist[[5]], snRNAlist[[6]],
                                    snRNAlist[[7]], snRNAlist[[8]], snRNAlist[[9]], snRNAlist[[10]], snRNAlist[[11]], snRNAlist[[12]],
                                    snRNAlist[[13]],snRNAlist[[14]], snRNAlist[[15]], snRNAlist[[16]], snRNAlist[[17]], snRNAlist[[18]],
                                    snRNAlist[[19]], snRNAlist[[20]]),
                add.cell.ids= c('WT_SIM4_1', 'WT_SIM4_2', 'WT_SIM4_3', 'WT_SIM6_1', 'WT_SIM6_2', 'WT_SIM6_3',
                                'WT_SIM8_1', 'WT_SIM8_2', 'WT_SIM8_3', 'WT_SIM10_1', 'WT_SIM10_2', 'WT_SIM10_3',
                  'myb_SIM4_1', 'myb_SIM4_2', 'myb_SIM6_1', 'myb_SIM6_2',  'myb_SIM8_1', 'myb_SIM8_2', 'myb_SIM10_1', 'myb_SIM10_2'))

#查看每个样本的细胞数
write.csv(tb1,'01_cell_number_per_sample_raw.csv', row.names = T)

# ---- Code cell 6 ----
p1 <- VlnPlot(snRNA1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, raster=FALSE, ncol = 3)
pdf("quality_control_1.pdf", width = 20, height = 5)
print(p1)
dev.off()

# ---- Code cell 7 ----
# filtering

snRNA2 <- subset(x = snRNA1,
                         subset= (nFeature_RNA > 200) &
                           (nFeature_RNA < 8000)&
                           (percent.mt < 10))

# ---- Code cell 8 ----
# 绘制小提琴图
p1 <- VlnPlot(snRNA2,
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              pt.size = 0,
              raster = FALSE,
              ncol = 3)

# 保存为PDF
pdf("quality_control_filtered.pdf", width = 20, height = 5)
print(p1)
dev.off()

# ---- Code cell 9 ----
#查看每个样本的细胞数
write.csv(tb2,'02_cell_number_per_sample_filtered.csv',row.names = T)

write.csv(tb3,'03_cell_number_per_timePoint_filtered.csv',row.names = T)

# ---- Code cell 10 ----
#################################################################################################################################3
snRNA2 <- SCTransform(snRNA2, vars.to.regress = "percent.mt", verbose = FALSE)
snRNA2 <- RunPCA(snRNA2, verbose = FALSE)

# ---- Code cell 11 ----
ep <- ElbowPlot(snRNA2, ndims = 50)
pdf("elbowplot.pdf", width = 4, height = 3)
print(ep)
dev.off()

# ---- Code cell 12 ----
saveRDS(object = snRNA2, file = "01_WT_myb31_filtered_sct.rds")

# ---- Code cell 13 ----
DefaultAssay(snRNA2) <- "SCT"

# Harmony整合（注意指定assay）
snRNA2 <- RunHarmony(snRNA2,
                           group.by.vars = "batch",  # 批次变量名
                           theta = 2,                # 调整批次校正强度（默认2）
                           lambda = 2,
                           assay.use = "SCT")

# ---- Code cell 14 ----
saveRDS(object = snRNA2, file = "01_WT_myb31_filtered_sct_harmony.rds")

# ---- Code cell 15 ----
snRNA2 <- FindNeighbors(snRNA2, dims = 1:40, reduction = "harmony", verbose = FALSE)
snRNA2 <- FindClusters(snRNA2, verbose = FALSE, resolution = 0.5)

snRNA2 <- RunUMAP(snRNA2, reduction = "harmony", dims = 1:40)

snRNA2$date <- factor(snRNA2$date, levels =c('WT_SIM4','WT_SIM6', 'WT_SIM8', 'WT_SIM10', 'myb_SIM4', 'myb_SIM6', 'myb_SIM8', 'myb_SIM10'))

# ---- Code cell 16 ----
snRNA2$orig.ident <- factor(snRNA2$orig.ident, levels = c('WT_SIM4_1', 'WT_SIM4_2', 'WT_SIM4_3', 'WT_SIM6_1', 'WT_SIM6_2', 'WT_SIM6_3', 'WT_SIM8_1', 'WT_SIM8_2', 'WT_SIM8_3', 'WT_SIM10_1', 'WT_SIM10_2', 'WT_SIM10_3',
                  'myb_SIM4_1', 'myb_SIM4_2', 'myb_SIM6_1', 'myb_SIM6_2', 'myb_SIM8_1', 'myb_SIM8_2', 'myb_SIM10_1', 'myb_SIM10_2'))

# ---- Code cell 17 ----
# Idents(snRNA2) <- factor(x = Idents(snRNA2), levels = sort(levels(snRNA2)))

# Visualization

plot1 <- DimPlot(snRNA2, reduction = "umap", pt.size = 0.01, label=T, ncol = 3)+
    theme_bw()+
    theme(panel.grid =element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
pdf("01_sct_UMAP_dim40_res0.5_1.pdf", width = 6, height = 4)
print(plot1)
dev.off()

# ---- Code cell 18 ----
plot2 <- DimPlot(snRNA2, reduction = "umap", pt.size = 0.01, split.by = 'date', label=F, ncol = 4)+
    theme_bw()+
    theme(panel.grid =element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
pdf("02_UMAP_dim40_res0.5_splitbydate_1.pdf", width = 12, height = 7)
print(plot2)
dev.off()

# ---- Code cell 19 ----
plot3 <- DimPlot(snRNA2, reduction = "umap", pt.size = 0.01, split.by = 'batch', label=T, ncol = 2)+
    theme_bw()+
    theme(panel.grid =element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
pdf("02_UMAP_dim40_res0.5_splitbybatch_1.pdf", width = 10, height = 4)
print(plot3)
dev.off()

# ---- Code cell 20 ----
Idents(object = snRNA2) <- "batch"

# ---- Code cell 21 ----
WT <- subset(x=snRNA2, idents = 'WT', invert = FALSE)

# ---- Code cell 22 ----
mutant <- subset(x=snRNA2, idents = 'Mutant', invert = FALSE)

# ---- Code cell 23 ----
saveRDS(object = WT, file = "01_WT_sct_dim40_res0.5.rds")

# ---- Code cell 24 ----
saveRDS(object = mutant, file = "01_mutant_sct_dim40_res0.5.rds")

# ---- Code cell 25 ----
saveRDS(object = snRNA2, file = "01_all_sct_dim40_res0.5.rds")

# ---- Code cell 26 ----
# How does cluster membership vary by replicate?
write.csv(x, file = "cells_eachCluster_per_date.csv")

x <- prop.table(table(Idents(snRNA2), snRNA2$date), margin = 2)
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
pdf("cell_number_change_by_date_1.pdf", width = 4, height = 3.5)
print(plot)
dev.off()

# ---- Code cell 27 ----
# How does cluster membership vary by replicate?
write.csv(x, file = "cells_eachCluster_per_batch.csv")

x <- prop.table(table(Idents(snRNA2), snRNA2$batch), margin = 2)
write.csv(x, file = "prop_cells_eachCluster_per_batch.csv")

# 加载必要的包
library(ggplot2)
library(reshape2)
library(dplyr)
library(viridis)

# 读取数据
data <- read.csv("prop_cells_eachCluster_per_batch.csv", header = TRUE)

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
pdf("cell_number_change_by_batch_1.pdf", width = 2, height = 3.5)
print(plot)
dev.off()

# ---- Code cell 28 ----
DefaultAssay(snRNA2) <- "RNA"
subset <- NormalizeData(snRNA2)

# ---- Code cell 29 ----
# dotplot
Idents(object = snRNA2) <- "date"

DefaultAssay(snRNA2) <- "RNA"
snRNA2 <- NormalizeData(snRNA2, normalization.method = "LogNormalize", scale.factor = 10000)

# ---- Code cell 30 ----
Idents(object = WT) <- "seurat_clusters"

# ---- Code cell 31 ----
DefaultAssay(WT) <- "RNA"
WT <- NormalizeData(WT, normalization.method = "LogNormalize", scale.factor = 10000)

# ---- Code cell 32 ----
select_markers <- read.table('shoot_related_genes.txt', header = F)
WT <- ScaleData(WT, features = as.character(select_markers$V1))

labels <- select_markers$V2

gene <- select_markers$V1  # 修改这里：使用第二列作为基因名称

#绘制气泡图
dp <- DotPlot(
    WT,
    features = gene,  # 现在使用第二列的基因符号
    assay = NULL,
    cols = c("#E8F396", "#9E0142"),
    col.min = 0,
    col.max = 1,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    group.by = NULL,
    split.by = NULL,
    cluster.idents = FALSE,
    scale = TRUE,
    scale.by = "radius",
    scale.min = NA,
    scale.max = NA
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = labels)  # 取消注释并添加这行来设置横坐标标签

# 保存气泡图到 PDF 文件
pdf("WT_shoot_related_genes_cluster_1.pdf", width = 12, height = 4)
print(dp)
dev.off()

# ---- Code cell 33 ----
# dotplot
Idents(object = snRNA2) <- "seurat_clusters"

DefaultAssay(snRNA2) <- "RNA"
snRNA2 <- NormalizeData(snRNA2, normalization.method = "LogNormalize", scale.factor = 10000)

# ---- Code cell 34 ----
# 仅读取和缩放热图需要的基因
select_markers <- read.table('/data/work/fig2c_genes_selected_7_1.txt', header = F)
WT <- ScaleData(WT, features = as.character(select_markers$V1))

gene <- select_markers$V1
labels <- select_markers$V2

#绘制气泡图
dp <- DotPlot(
    WT,
    features = gene,
    assay = NULL,
    cols = c(c("#E8F396", "#9E0142")),
    col.min = 0,
    col.max = 1.2,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    group.by = NULL,
    split.by = NULL,
    cluster.idents = FALSE,
    scale = TRUE,
    scale.by = "radius",
    scale.min = NA,
    scale.max = NA
  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(labels = labels)  # 取消注释并添加这行来设置横坐标标签

# 保存气泡图到 PDF 文件
pdf("WT_markers_dotplot_1.pdf", width = 7, height = 6)
print(dp)
dev.off()

# ---- Code cell 35 ----
# 仅读取和缩放热图需要的基因
select_markers <- read.table('/data/work/fig2c_genes_selected_7_2.txt', header = F)
WT <- ScaleData(WT, features = as.character(select_markers$V1))

gene <- select_markers$V1
labels <- select_markers$V2

#绘制气泡图
dp <- DotPlot(
    WT,
    features = gene,
    assay = NULL,
    cols = c(c("#E8F396", "#9E0142")),
    col.min = 0,
    col.max = 1.2,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    group.by = NULL,
    split.by = NULL,
    cluster.idents = FALSE,
    scale = TRUE,
    scale.by = "radius",
    scale.min = NA,
    scale.max = NA
  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(labels = labels)  # 取消注释并添加这行来设置横坐标标签

# 保存气泡图到 PDF 文件
pdf("WT_markers_dotplot_2.pdf", width = 12, height = 6)
print(dp)
dev.off()

# ---- Code cell 36 ----
# 仅读取和缩放热图需要的基因
select_markers <- read.table('/data/work/fig2c_genes_selected_7_3.txt', header = F)
WT <- ScaleData(WT, features = as.character(select_markers$V1))

gene <- select_markers$V1
labels <- select_markers$V2

#绘制气泡图
dp <- DotPlot(
    WT,
    features = gene,
    assay = NULL,
    cols = c(c("#E8F396", "#9E0142")),
    col.min = 0,
    col.max = 1.2,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    group.by = NULL,
    split.by = NULL,
    cluster.idents = FALSE,
    scale = TRUE,
    scale.by = "radius",
    scale.min = NA,
    scale.max = NA
  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(labels = labels)  # 取消注释并添加这行来设置横坐标标签

# 保存气泡图到 PDF 文件
pdf("WT_markers_dotplot_3.pdf", width = 5, height = 6)
print(dp)
dev.off()

# ---- Code cell 37 ----
# 仅读取和缩放热图需要的基因
select_markers <- read.table('/data/work/fig2c_genes_selected_7_4.txt', header = F)
WT <- ScaleData(WT, features = as.character(select_markers$V1))

gene <- select_markers$V1
labels <- select_markers$V2

#绘制气泡图
dp <- DotPlot(
    WT,
    features = gene,
    assay = NULL,
    cols = c(c("#E8F396", "#9E0142")),
    col.min = 0,
    col.max = 1.2,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    group.by = NULL,
    split.by = NULL,
    cluster.idents = FALSE,
    scale = TRUE,
    scale.by = "radius",
    scale.min = NA,
    scale.max = NA
  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(labels = labels)  # 取消注释并添加这行来设置横坐标标签

# 保存气泡图到 PDF 文件
pdf("WT_markers_dotplot_4.pdf", width = 6, height = 6)
print(dp)
dev.off()

# ---- Code cell 38 ----
snRNA <- readRDS('/data/work/01_sct_dim20_res0.2.rds')

# ---- Code cell 39 ----
ep <- ElbowPlot(snRNA, ndims = 100)
pdf("elbowplot.pdf", width = 4, height = 3)
print(ep)
dev.off()

# ---- Code cell 40 ----
snRNA

# ---- Code cell 41 ----
snRNA <- FindNeighbors(snRNA, dims = 1:30, reduction = "pca", verbose = FALSE)
snRNA <- FindClusters(snRNA, verbose = FALSE, resolution = 0.3)

snRNA <- RunUMAP(snRNA, reduction = "pca", dims = 1:30)

snRNA$date <- factor(snRNA$date, levels =c('WT_SIM4','WT_SIM6', 'WT_SIM8', 'WT_SIM10', 'myb_SIM4', 'myb_SIM6', 'myb_SIM8', 'myb_SIM10'))

# ---- Code cell 42 ----
snRNA$orig.ident <- factor(snRNA$orig.ident, levels = c('WT_SIM4_1', 'WT_SIM4_2', 'WT_SIM4_3', 'WT_SIM6_1', 'WT_SIM6_2', 'WT_SIM6_3', 'WT_SIM8_1', 'WT_SIM8_2', 'WT_SIM8_3', 'WT_SIM10_1', 'WT_SIM10_2', 'WT_SIM10_3',
                  'myb_SIM4_1', 'myb_SIM4_2', 'myb_SIM4_3', 'myb_SIM6_1', 'myb_SIM6_2', 'myb_SIM6_3', 'myb_SIM8_1', 'myb_SIM8_2', 'myb_SIM8_3', 'myb_SIM10_1', 'myb_SIM10_2', 'myb_SIM10_3'))

# ---- Code cell 43 ----
plot1 <- DimPlot(snRNA, reduction = "umap", pt.size = 0.01, label=T, ncol = 3)+
    theme_bw()+
    theme(panel.grid =element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
pdf("01_recluster_sct_UMAP_dim30_res0.3_1.pdf", width = 6, height = 4)
print(plot1)
dev.off()

# ---- Code cell 44 ----
plot2 <- DimPlot(snRNA, reduction = "umap", pt.size = 0.01, split.by = 'date', label=F, ncol = 4)+
    theme_bw()+
    theme(panel.grid =element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
pdf("02_recluster_UMAP_dim30_res0.3_splitdate_1.pdf", width = 10, height = 4)
print(plot2)
dev.off()

# ---- Code cell 45 ----
# How does cluster membership vary by replicate?
write.csv(x, file = "snRNA_dim30_res0.3_cells_eachCluster_per_date.csv")

x <- prop.table(table(Idents(snRNA), snRNA$date), margin = 2)
write.csv(x, file = "snRNA_dim30_res0.3_prop_cells_eachCluster_per_date.csv")

# 加载必要的包
library(ggplot2)
library(reshape2)
library(dplyr)
library(viridis)

# 读取数据
data <- read.csv("snRNA_dim30_res0.3_prop_cells_eachCluster_per_date.csv", header = TRUE)

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
pdf("snRNA_dim30_res0.3_cell_number_change_1.pdf", width = 4, height = 3.5)
print(plot)
dev.off()

# ---- Code cell 46 ----
saveRDS(object = snRNA, file = "01_sct_dim30_res0.3.rds")

# ---- Code cell 47 ----
DefaultAssay(snRNA) <- "RNA"
snRNA <- NormalizeData(snRNA)

# ---- Code cell 48 ----
# 仅读取和缩放热图需要的基因
select_markers <- read.table('/data/work/fig2c_genes_selected_7_1.txt', header = F)
snRNA2 <- ScaleData(snRNA, features = as.character(select_markers$V1))

gene <- select_markers$V1
labels <- select_markers$V2

#绘制气泡图
dp <- DotPlot(
    snRNA,
    features = gene,
    assay = NULL,
    cols = c(c("#E8F396", "#9E0142")),
    col.min = 0,
    col.max = 1.2,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    group.by = NULL,
    split.by = NULL,
    cluster.idents = FALSE,
    scale = TRUE,
    scale.by = "radius",
    scale.min = NA,
    scale.max = NA
  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(labels = labels)  # 取消注释并添加这行来设置横坐标标签

# 保存气泡图到 PDF 文件
pdf("dim30_res0.3_markers_dotplot_1.pdf", width = 7, height = 6)
print(dp)
dev.off()

# ---- Code cell 49 ----
# 仅读取和缩放热图需要的基因
select_markers <- read.table('/data/work/fig2c_genes_selected_7_2.txt', header = F)
snRNA2 <- ScaleData(snRNA, features = as.character(select_markers$V1))

gene <- select_markers$V1
labels <- select_markers$V2

#绘制气泡图
dp <- DotPlot(
    snRNA,
    features = gene,
    assay = NULL,
    cols = c(c("#E8F396", "#9E0142")),
    col.min = 0,
    col.max = 1.2,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    group.by = NULL,
    split.by = NULL,
    cluster.idents = FALSE,
    scale = TRUE,
    scale.by = "radius",
    scale.min = NA,
    scale.max = NA
  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(labels = labels)  # 取消注释并添加这行来设置横坐标标签

# 保存气泡图到PDF文件
pdf("dim30_res0.3_markers_dotplot_2.pdf", width = 5, height = 6)
print(dp)
dev.off()

# ---- Code cell 50 ----
# 仅读取和缩放热图需要的基因
select_markers <- read.table('/data/work/fig2c_genes_selected_7_3.txt', header = F)
snRNA2 <- ScaleData(snRNA, features = as.character(select_markers$V1))

gene <- select_markers$V1
labels <- select_markers$V2

#绘制气泡图
dp <- DotPlot(
    snRNA,
    features = gene,
    assay = NULL,
    cols = c(c("#E8F396", "#9E0142")),
    col.min = 0,
    col.max = 1.2,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    group.by = NULL,
    split.by = NULL,
    cluster.idents = FALSE,
    scale = TRUE,
    scale.by = "radius",
    scale.min = NA,
    scale.max = NA
  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(labels = labels)  # 取消注释并添加这行来设置横坐标标签

# 保存气泡图到 PDF 文件
pdf("dim30_res0.3_markers_dotplot_3.pdf", width = 4, height = 6)
print(dp)
dev.off()

# ---- Code cell 51 ----
# 仅读取和缩放热图需要的基因
select_markers <- read.table('/data/work/fig2c_genes_selected_7_4.txt', header = F)
snRNA2 <- ScaleData(snRNA, features = as.character(select_markers$V1))

gene <- select_markers$V1
labels <- select_markers$V2

#绘制气泡图
dp <- DotPlot(
    snRNA,
    features = gene,
    assay = NULL,
    cols = c(c("#E8F396", "#9E0142")),
    col.min = 0,
    col.max = 1.2,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    group.by = NULL,
    split.by = NULL,
    cluster.idents = FALSE,
    scale = TRUE,
    scale.by = "radius",
    scale.min = NA,
    scale.max = NA
  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(labels = labels)  # 取消注释并添加这行来设置横坐标标签

# 保存气泡图到 PDF 文件
pdf("dim30_res0.3_markers_dotplot_4.pdf", width = 5, height = 6)
print(dp)
dev.off()

# ---- Code cell 52 ----
# 仅读取和缩放热图需要的基因
select_markers <- read.table('/data/work/fig2c_genes_selected_7_5.txt', header = F)
snRNA2 <- ScaleData(snRNA, features = as.character(select_markers$V1))

gene <- select_markers$V1
labels <- select_markers$V2

#绘制气泡图
dp <- DotPlot(
    snRNA,
    features = gene,
    assay = NULL,
    cols = c(c("#E8F396", "#9E0142")),
    col.min = 0,
    col.max = 1.2,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    group.by = NULL,
    split.by = NULL,
    cluster.idents = FALSE,
    scale = TRUE,
    scale.by = "radius",
    scale.min = NA,
    scale.max = NA
  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(labels = labels)  # 取消注释并添加这行来设置横坐标标签

# 保存气泡图到 PDF 文件
pdf("dim30_res0.3_markers_dotplot_5.pdf", width = 5, height = 6)
print(dp)
dev.off()

# ---- Code cell 53 ----
Idents(snRNA) <- 'seurat_clusters'

# ---- Code cell 54 ----
# 仅读取和缩放热图需要的基因
select_markers <- read.table('/data/work/fig2c_genes_selected_7_6.txt', header = F)
snRNA2 <- ScaleData(snRNA, features = as.character(select_markers$V1))

gene <- select_markers$V1
labels <- select_markers$V2

#绘制气泡图
dp <- DotPlot(
    snRNA,
    features = gene,
    assay = NULL,
    cols = c(c("#E8F396", "#9E0142")),
    col.min = 0,
    col.max = 1.2,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    group.by = NULL,
    split.by = NULL,
    cluster.idents = FALSE,
    scale = TRUE,
    scale.by = "radius",
    scale.min = NA,
    scale.max = NA
  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(labels = labels)  # 取消注释并添加这行来设置横坐标标签

# 保存气泡图到 PDF 文件
pdf("dim35_res0.3_markers_dotplot_5.pdf", width = 7, height = 5)
print(dp)
dev.off()

# ---- Code cell 55 ----
pdf("dim35_res0.3_markers_date_dotplot_2.pdf", width = 5.5, height = 5)
print(dp)
dev.off()

# ---- Code cell 56 ----
Idents(snRNA) <- snRNA$batch

# ---- Code cell 57 ----
mutant <- subset(x=snRNA, idents = 'Mutant', invert = FALSE)

# ---- Code cell 58 ----
Idents(mutant) <- mutant$seurat_clusters

# ---- Code cell 59 ----
plot1 <- DimPlot(mutant, reduction = "umap", pt.size = 0.01, label=T)+
    theme_bw()+
    theme(panel.grid =element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
pdf("01_sct_mutant_UMAP_dim30_res0.3_1.pdf", width = 6, height = 4)
print(plot1)
dev.off()

# ---- Code cell 60 ----
plot2 <- DimPlot(mutant, reduction = "umap", pt.size = 0.01, split.by = 'date', label=T, ncol = 4)+
    theme_bw()+
    theme(panel.grid =element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
pdf("02_mutant_UMAP_dim30_res0.3_splitdate_1.pdf", width = 12, height = 3)
print(plot2)
dev.off()

# ---- Code cell 61 ----
DefaultAssay(mutant) <- "RNA"
mutant <- NormalizeData(mutant, normalization.method = "LogNormalize", scale.factor = 10000)

# ---- Code cell 62 ----
# 仅读取和缩放热图需要的基因
select_markers <- read.table('/data/work/fig2c_genes_selected_7_6.txt', header = F)
mutant <- ScaleData(mutant, features = as.character(select_markers$V1))

gene <- select_markers$V1
labels <- select_markers$V2

#绘制气泡图
dp <- DotPlot(
    mutant,
    features = gene,
    assay = NULL,
    cols = c(c("#E8F396", "#9E0142")),
    col.min = 0,
    col.max = 1.2,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    group.by = NULL,
    split.by = NULL,
    cluster.idents = FALSE,
    scale = TRUE,
    scale.by = "radius",
    scale.min = NA,
    scale.max = NA
  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(labels = labels)  # 取消注释并添加这行来设置横坐标标签

# 保存气泡图到 PDF 文件
pdf("mutant_markers_dotplot_6.pdf", width = 7, height = 4.5)
print(dp)
dev.off()
