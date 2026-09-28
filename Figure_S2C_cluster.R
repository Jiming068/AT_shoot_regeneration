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
set.seed(128)  #设置随机数种子，使结果可重复

# options(future.globals.maxSize = 5 * 1024^3)

##==合并数据集==##
#使用merge函数合并seurat对象

dir_1 <- c(
        file.path(external_data_dir, "04.Matrix"),
        file.path(external_data_dir, "04.Matrix"),
        file.path(external_data_dir, "04.Matrix"),
        file.path(external_data_dir, "04.Matrix"),
        file.path(external_data_dir, "04.Matrix"),
        file.path(external_data_dir, "04.Matrix"),
        file.path(external_data_dir, "04.Matrix"),
        file.path(external_data_dir, "04.Matrix"),
        file.path(external_data_dir, "04.Matrix"),
        file.path(external_data_dir, "04.Matrix"),
        file.path(external_data_dir, "04.Matrix"),
        file.path(external_data_dir, "04.Matrix"))

dir_2 <- c(
        file.path(external_data_dir, "04.Matrix"),
        file.path(external_data_dir, "04.Matrix"),
        file.path(external_data_dir, "04.Matrix"),
        file.path(external_data_dir, "04.Matrix"),
        file.path(external_data_dir, "04.Matrix"),
        file.path(external_data_dir, "04.Matrix"))

dir <- c(dir_1, dir_2)

sample_name <- c('SIM4_1', 'SIM4_2', 'SIM4_3', 'SIM6_1', 'SIM6_2', 'SIM6_3',
                 'SIM8_1', 'SIM8_2', 'SIM8_3', 'SIM10_1', 'SIM10_2', 'SIM10_3','SIM14_1', 'SIM14_2', 'SIM14_3',
                 'SIM16_1', 'SIM16_2', 'SIM16_3')
date <- c('SIM4', 'SIM4', 'SIM4', 'SIM6', 'SIM6', 'SIM6',
          'SIM8', 'SIM8', 'SIM8', 'SIM10', 'SIM10', 'SIM10', 'SIM14', 'SIM14', 'SIM14',
          'SIM16', 'SIM16', 'SIM16')

#以下代码会把每个样本的数据创建一个seurat对象，并存放到列表scRNAlist里
snRNAlist <- list()
for(i in 1:length(dir)){
counts <- Read10X(data.dir = dir[i], gene.column = 1)
snRNAlist[[i]] <- CreateSeuratObject(counts, project = sample_name[i], min.cells = 3, min.features = 100)
# include protoplasting-induced genes
snRNAlist[[i]][['date']] <- date[i]
snRNAlist[[i]][['percent.mt']] <- PercentageFeatureSet(snRNAlist[[i]], pattern = '^ATM')
snRNAlist[[i]][["percent.cp"]] <- PercentageFeatureSet(snRNAlist[[i]], pattern = '^ATC')
}

#saveRDS(object = snRNAlist, file = "00_data_raw_20220706.rds")

################################################################################################################################3
#使用merge函数将11个seurat对象合并成一个seurat对象
# scRNA1 <- readRDS("00_data_raw_20220412.rds")

snRNA1 <- merge(snRNAlist[[1]], y=c(snRNAlist[[2]], snRNAlist[[3]],
                snRNAlist[[4]], snRNAlist[[5]], snRNAlist[[6]], snRNAlist[[7]], snRNAlist[[8]], snRNAlist[[9]], snRNAlist[[10]], snRNAlist[[11]], snRNAlist[[12]], snRNAlist[[13]],
                snRNAlist[[14]], snRNAlist[[15]], snRNAlist[[16]], snRNAlist[[17]], snRNAlist[[18]]),
                add.cell.ids= c('SIM4_1', 'SIM4_2', 'SIM4_3', 'SIM6_1', 'SIM6_2', 'SIM6_3',
                 'SIM8_1', 'SIM8_2', 'SIM8_3', 'SIM10_1', 'SIM10_2', 'SIM10_3', 'SIM14_1', 'SIM14_2', 'SIM14_3',
                 'SIM16_1', 'SIM16_2', 'SIM16_3'))

#查看每个样本的细胞数
write.csv(tb1,'01_cell_number_per_sample_raw_20230207.csv', row.names = T)

# Visualize QC metrics as a violin plot
v1 <- VlnPlot(snRNA1, features = c("nFeature_RNA", "nCount_RNA", 'percent.mt', 'percent.cp'), slot = "data", pt.size = 0, ncol = 4)
pdf("QC_metrics_20230207.pdf", width=20, height =4)
v1
dev.off()

#################################################################################################################################3

# filtering

snRNA2 <- subset(x = snRNA1,
                         subset= (nFeature_RNA > 100) &
                           (nFeature_RNA < 6000)&
                           (percent.mt < 10))

v2 <- VlnPlot(snRNA2, features = c("nFeature_RNA", "nCount_RNA", 'percent.mt', 'percent.cp'), slot = "data", pt.size = 0, ncol = 4)
pdf("QC_metrics_filtered_20230207.pdf", width=20, height =4)
v2
dev.off()

#查看每个样本的细胞数
write.csv(tb2,'02_cell_number_per_sample_filtered_20230207.csv',row.names = T)

write.csv(tb3,'03_cell_number_per_timePoint_filtered_20230207.csv',row.names = T)

snRNA2 <- SCTransform(snRNA2, vars.to.regress = "percent.mt", verbose = FALSE)

snRNA2 <- RunPCA(snRNA2, verbose = FALSE)

snRNA2 <- FindNeighbors(snRNA2, dims = 1:30, verbose = FALSE)
snRNA2 <- FindClusters(snRNA2, verbose = FALSE, resolution = 1.2)

snRNA2 <- RunUMAP(snRNA2, reduction = "pca", dims = 1:30)

snRNA2$date <- factor(snRNA2$date, levels =c('SIM4', 'SIM6', 'SIM8', 'SIM10', 'SIM14', 'SIM16'))

snRNA2$orig.ident <- factor(snRNA2$orig.ident, levels = c('SIM4_1', 'SIM4_2', 'SIM4_3', 'SIM6_1', 'SIM6_2', 'SIM6_3',
                 'SIM8_1', 'SIM8_2', 'SIM8_3', 'SIM10_1', 'SIM10_2', 'SIM10_3', 'SIM14_1', 'SIM14_2', 'SIM14_3',
                 'SIM16_1', 'SIM16_2', 'SIM16_3'))

# Idents(snRNA2) <- factor(x = Idents(snRNA2), levels = sort(levels(snRNA2)))

# Visualization

plot1 <- DimPlot(snRNA2, reduction = "umap", label=T)
plot2 <- DimPlot(snRNA2, reduction = "umap", group.by='orig.ident', label=T)
plot3 <- DimPlot(snRNA2, reduction = "umap", group.by='date', label=T)
plot4 <- DimPlot(snRNA2, reduction = "umap", split.by='date', label=T)
plote <- plot1+plot2+plot3

pdf("01_standard_UMAP_sample_dim30_res1.2_d4-16_rm12_20230207.pdf", width = 25, height = 5)
print(plote)
dev.off()

pdf("02_standard_UMAP_date_dim30_res1.2_d4-16_rm12_20230207.pdf", width = 35, height = 5)
print(plot4)
dev.off()

plot1 <- DimPlot(snRNA2, reduction = "umap", label=T)
plot2 <- DimPlot(snRNA2, reduction = "umap", group.by='orig.ident', label=T)
plot3 <- DimPlot(snRNA2, reduction = "umap", group.by='date', label=T)
plot4 <- DimPlot(snRNA2, reduction = "umap", split.by='date', label=T, ncol = 3)
plote <- plot1+plot2+plot3

pdf("02_sct_UMAP_sample_dim30_res1.2_d4-16_rm12_20230207_1.pdf", width = 20, height = 5)
print(plote)
dev.off()

pdf("02_sct_UMAP_date_dim30_res1.2_d4-16_rm12_20230207_1.pdf", width = 15, height = 8)
print(plot4)
dev.off()

saveRDS(object = snRNA2, file = "01_sct_dim30_res1.2_d4-16_rm12_20230207.rds")

Idents(object = snRNA2) <- "seurat_clusters"

data_markers <- FindAllMarkers(snRNA2, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2)
write.table(data_markers, file='root_DEGs_dim30_res1.2_d4-16_rm12_20230207.txt', col.names = TRUE, sep = '\t')

write.csv(x, file = "root_derived_cells_per_dated4-16_rm12_20230207.csv", row.names = FALSE)

write.csv(x, file = "root_derived_cells_per_sampled4-16_rm12_20230207.csv", row.names = FALSE)

# What proportion of cells are in each cluster?
prop.table(table(Idents(snRNA2)))

# How does cluster membership vary by replicate?
write.csv(x, file = "root_derived_cells_eachCluster_per_dated4-16_rm12_20230207.csv")

x <- prop.table(table(Idents(snRNA2), snRNA2$date), margin = 2)
write.csv(x, file = "root_derived_prop_cells_eachCluster_per_dated4-16_rm12_20230207.csv")
