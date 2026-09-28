# Portable path configuration (override with environment variables).
project_dir <- Sys.getenv("SHOOT_REGEN_PROJECT_DIR", unset = ".")
output_root <- Sys.getenv("SHOOT_REGEN_OUTPUT_DIR", unset = file.path(project_dir, "02_results"))
external_data_dir <- Sys.getenv("SHOOT_REGEN_EXTERNAL_DATA_DIR", unset = file.path(project_dir, "external_data"))
work_dir <- Sys.getenv("SHOOT_REGEN_WORK_DIR", unset = ".")

# 2022-09-06
# root-to-calli

dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)

library(Seurat)
library(patchwork)
library(tidyverse)
library(reshape2)
library(SeuratData)
library(dplyr)
library(SeuratDisk)
library(pheatmap)
library(harmony)
library(matrixStats)
packageVersion("matrixStats")
library(glmGamPoi)
set.seed(128)  #设置随机数种子，使结果可重复
options(future.globals.maxSize = 100 * 1024^3)

##UMAP recolor
snRNA <- readRDS('01_all_sct_dim30_res0.5.rds')
Idents(snRNA) <- 'seurat_clusters'

# gene and umi counts
#新分群数据。为每一date不同样品进行颜色配色测试
Idents(object = snRNA) <- "batch"

subset <- subset(x=snRNA, idents = c('Mutant'), invert = FALSE)

Idents(object = subset) <- "orig.ident"

#最后选了这个配色！
pal<-c("#e6bb56","#e9cb95","#3a3d7a","#614c97","#a75293","#c66fad","#2480b8","#afd9ec")
y <- DimPlot(subset, reduction = "umap", label = T, pt.size = 0.1,cols =pal)+
theme_bw()+
theme(panel.grid =element_blank(),
panel.border = element_blank(),
axis.text.y = element_blank(),
axis.ticks.y = element_blank(),
axis.ticks.x = element_blank(),
axis.text.x = element_blank())
pdf("WT_mutant_UMAP_date_1.pdf", width = 8, height = 6)
print(y)
dev.off()

#小提琴质控图
v1 <- VlnPlot(subset, features = c("nFeature_RNA", "nCount_RNA"), log = TRUE, pt.size = 0, ncol = 2, cols = pal)
pdf("mutant_violin_plot.pdf", width=20, height =4)
v1
dev.off()

#箱线质控图
nCount_RNA <- subset@meta.data$nCount_RNA
nFeature_RNA <- subset@meta.data$nFeature_RNA
qc_all <- data.frame(nCount_RNA,nFeature_RNA)
qc_all$all <- "all"
qc_all$all = as.factor(qc_all$all)

qc_ncbox <- ggplot(qc_all, aes(x = all, y =nCount_RNA)) +
    geom_boxplot(fill = '#FF6666', color = "black") +
    theme_classic()
ggsave("mutant_qc_nCount_boxplot.pdf",plot = qc_ncbox, width = 4, height = 6)

qc_nfbox <- ggplot(qc_all, aes(x = all, y =nFeature_RNA)) +
  geom_boxplot(fill = '#99CCFF', color = "black") +
  theme_classic()
ggsave("mutant_qc_nFeature_boxplot.pdf",plot = qc_nfbox, width = 4, height = 6)
