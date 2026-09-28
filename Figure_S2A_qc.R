# Portable path configuration (override with environment variables).
project_dir <- Sys.getenv("SHOOT_REGEN_PROJECT_DIR", unset = ".")
output_root <- Sys.getenv("SHOOT_REGEN_OUTPUT_DIR", unset = file.path(project_dir, "02_results"))
external_data_dir <- Sys.getenv("SHOOT_REGEN_EXTERNAL_DATA_DIR", unset = file.path(project_dir, "external_data"))
work_dir <- Sys.getenv("SHOOT_REGEN_WORK_DIR", unset = ".")

#新分群数据。为每一date不同样品进行颜色配色测试

dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)
library(tidyverse)
library(Seurat)
library(patchwork)
library(reshape2)
set.seed(128)  #设置随机数种子，使结果可重复
snRNA2 <- readRDS(file.path(project_dir, "01_obj/01_WT_snRNA/01_snRNA_res1.2_SIM4-SIM16.rds"))

Idents(object = snRNA2) <- "orig.ident"

#最后选了这个配色！
pal<-c("#e6bb56","#e9cb95","#e1ebb5","#3a3d7a","#614c97","#9f9dd0","#a75293","#c66fad","#d7a0cc","#2480b8","#afd9ec","#a8d3ee","#34a356","#60b568","#a2d49b","#e55828","#d5626a","#e7989e")
y <- DimPlot(snRNA2, reduction = "umap", label = T, pt.size = 0.1,cols =pal)+
theme_bw()+
theme(panel.grid =element_blank(),
panel.border = element_blank(),
axis.text.y = element_blank(),
axis.ticks.y = element_blank(),
axis.ticks.x = element_blank(),
axis.text.x = element_blank())
pdf("root_umap_date_20230801_5.pdf", width = 8, height = 6)
print(y)
dev.off()

#小提琴质控图
v1 <- VlnPlot(snRNA2, features = c("nFeature_RNA", "nCount_RNA"), log = TRUE, pt.size = 0, ncol = 2, cols = pal)
pdf("root_snRNA_violin_plot.pdf", width=20, height =4)
v1
dev.off()

#箱线质控图
nCount_RNA <- snRNA2@meta.data$nCount_RNA
nFeature_RNA <- snRNA2@meta.data$nFeature_RNA
qc_all <- data.frame(nCount_RNA,nFeature_RNA)
qc_all$all <- "all"
qc_all$all = as.factor(qc_all$all)

qc_ncbox <- ggplot(qc_all, aes(x = all, y =nCount_RNA)) +
    geom_boxplot(fill = '#FF6666', color = "black") +
    theme_classic()
ggsave("root_snRNA_qc_nCount_boxplot.pdf",plot = qc_ncbox, width = 4, height = 6)

qc_nfbox <- ggplot(qc_all, aes(x = all, y =nFeature_RNA)) +
  geom_boxplot(fill = '#99CCFF', color = "black") +
  theme_classic()
ggsave("root_snRNA_qc_nFeature_boxplot.pdf",plot = qc_nfbox, width = 4, height = 6)
