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

subset <- readRDS(file.path(external_data_dir, "snRNA_subset_sct_dim19_res0.2.rds"))
Idents(subset) <- 'seurat_clusters'

# ---------------------- 0b. 原始细胞类群配色方案 ----------------------
cell_colors <- c(
  # 组1: 0和8群使用对比配色（暖色系中的对比）
  "0" = "#FFB74D",  # 柔和的橙黄色
  "8" = "#EF5350",

  # 组2: 2和3群使用同一色系但不要太相似（柔和绿色系）
  "2" = "#81C784",  # 柔和的浅绿色
  "3" = "#66BB6A"   # 柔和的绿色 - 与2同系但更饱和
)

# ---------------------- 1. 设置要提取的cluster ----------------------
target_clusters <- c('0', '2', '3')

# 如果目前的Idents不是cluster编号，先设置为cluster列
# Idents(seurat_obj) <- "seurat_clusters"   # 改成你的cluster元数据列名

# ---------------------- 2. Subset目标cluster ----------------------
sub_obj <- subset(subset, idents = target_clusters)

# 确认subset后的细胞数和原始cluster分布
table(Idents(sub_obj))

# ---------------------- 3. 重新计算高变基因（可选但推荐） ----------------------
# subset后的细胞群体组成变了，原来基于全体细胞算出的高变基因
# 未必最适合这个子集，重新计算能让subclustering更准确地捕捉这几个
# cluster内部的异质性

DefaultAssay(sub_obj) <- "RNA"

sub_obj <- SCTransform(sub_obj, vars.to.regress =c("percent.mt", 'percent.cp'), verbose = FALSE)

# ---------------------- 5. 重新计算PCA（基于subset后的高变基因） ----------------------
sub_obj <- RunPCA(sub_obj, features = VariableFeatures(sub_obj), npcs = 30)

ep <- ElbowPlot(subset, ndims = 50)
pdf("elbowplot.pdf", width = 4, height = 3)
print(ep)
dev.off()

# 可以用ElbowPlot检查合适的PC数量，这里默认用前30个PC
# ElbowPlot(sub_obj, ndims = 30)

# ---------------------- 6. 重新计算邻居图（基于新的PCA） ----------------------
n_pcs <- 20   # 根据ElbowPlot结果调整，默认用前20个PC
sub_obj <- FindNeighbors(sub_obj, dims = 1:n_pcs)

# ---------------------- 7. 提高分辨率重新聚类（subclustering核心步骤） ----------------------
# 分辨率越高，划分出的cluster越多、越细。原来的5个cluster在这里
# 会被进一步拆分成更多的子cluster
new_resolution <- 0.4   # 比原来聚类时用的分辨率更高，可根据效果调整
sub_obj <- FindClusters(sub_obj, resolution = new_resolution)

# 新的subcluster编号会存储在 sub_obj$seurat_clusters
# 也可以给它另存一列，避免覆盖原始cluster信息
sub_obj$subcluster <- sub_obj$seurat_clusters
sub_obj$original_cluster <- sub_obj$seurat_clusters  # 如需保留原cluster可在subset前先备份

table(sub_obj$subcluster)

# ---------------------- 7b. 重新聚类后0-11号细胞类群的浅暖色系配色 ----------------------
# 12个子群使用浅暖色系（米黄、浅橙、浅杏、浅珊瑚、浅粉、浅棕等），
# 整体色调偏浅、偏暖，彼此又有足够区分度
subcluster_colors <- c(
  "0"  = "#EBCC50",  # Golden Yellow  柔和金黄
  "1"  = "#EBAB50",  # Marigold       万寿菊橙
  "2"  = "#EB8950",  # Tangerine      橘橙
  "3"  = "#EB6750",  # Vermilion      朱红
  "4"  = "#EB505A",  # Poppy Red      罂粟红
  "5"  = "#EB507C",  # Rose           玫红
  "6"  = "#86711C",  # Olive Gold     深橄榄金（偏暖褐）
  "7"  = "#865A1C",  # Amber Brown    琥珀棕
  "8"  = "#86431C"  # Rust           铁锈橙棕

)
# 如果RC state数量或编号与上面不一致，请根据 rc_order 的实际结果调整
# 如果实际subcluster数量与上面颜色数量不一致（resolution不同会导致
# cluster数变化），请根据 levels(sub_obj$subcluster) 的实际结果调整颜色数量：
# levels(sub_obj$subcluster)

# ---------------------- 8. 保持UMAP坐标不变（关键要求） ----------------------
# 不要重新运行 RunUMAP！直接复用原始Seurat对象里已经算好的UMAP坐标，
# 只按细胞barcode取出subset对应的那部分坐标
original_umap <- Embeddings(subset, reduction = "umap")
sub_umap <- original_umap[colnames(sub_obj), ]

# 将原始UMAP坐标重新写回subset对象的reduction中，
# 这样sub_obj的UMAP图坐标会和原图完全一致，不会因为重新计算而偏移
sub_obj[["umap"]] <- CreateDimReducObject(embeddings = sub_umap,
                                           key = "UMAP_",
                                           assay = DefaultAssay(sub_obj))

# ---------------------- 9. 可视化验证：坐标应与原图一致 ----------------------
p_original <- DimPlot(subset, reduction = "umap",
                       cells = colnames(sub_obj),
                       group.by = "seurat_clusters",
                       label = TRUE) +
  scale_color_manual(values = cell_colors) +
  ggtitle("original")

p_subcluster <- DimPlot(sub_obj, reduction = "umap",
                         group.by = "subcluster",
                         label = TRUE) +
  scale_color_manual(values = subcluster_colors) +
  ggtitle(paste0("(subcluster, resolution=", new_resolution, ")"))

print(p_original)
print(p_subcluster)

# ---------------------- 9b. 保存UMAP图至PDF文件 ----------------------
ggsave("umap_original_cluster_labels_2.pdf", plot = p_original, width = 6, height = 5)
ggsave("umap_subcluster_labels_2.pdf", plot = p_subcluster, width = 6, height = 5)

# 如果想把两张图合并保存到同一个PDF文件（各占一页），可以用下面这种方式代替：
pdf("umap_original_vs_subcluster_2.pdf", width = 6, height = 5)
print(p_original)
print(p_subcluster)
dev.off()

# 如果想把两张图并排放在同一页展示，可以用patchwork拼图：
library(patchwork)
combined_umap <- p_original + p_subcluster
ggsave("umap_original_vs_subcluster_sidebyside_2.pdf", plot = combined_umap,
        width = 11, height = 5)

# How does cluster membership vary by replicate?
write.csv(x, file = "subset_cells_eachCluster_per_date_1.csv")

x <- prop.table(table(Idents(sub_obj), sub_obj$date), margin = 2)
write.csv(x, file = "subset_prop_cells_eachCluster_per_date_2.csv")

# ---------------------- 10. 保存subset对象 ----------------------
saveRDS(sub_obj, "seurat_subset_C0-2-3-reclustered_dim20_res0.4.rds")

# ============================================================
# 【可选】如果想尝试多个分辨率、快速对比subcluster数量再决定用哪个
# ============================================================
# library(clustree)
# for (res in c(0.4, 0.6, 0.8, 1.0, 1.2, 1.5)) {
#   sub_obj <- FindClusters(sub_obj, resolution = res)
# }
# clustree(sub_obj, prefix = "RNA_snn_res.")   # 根据assay名调整prefix

# ============================================================
# 【可选】如果只想对某几个原cluster内部做subclustering，
# 且希望最终结果标注为"原cluster_子编号"（如 "0_0", "0_1", "2_0"...）
# 而不是从0开始重新编号，可用Seurat自带的subcluster功能：
# ============================================================
# seurat_obj <- FindSubCluster(seurat_obj,
#                               cluster = "0",              # 每次只能对一个原cluster操作
#                               graph.name = "RNA_snn",
#                               subcluster.name = "sub.cluster",
#                               resolution = 1.2)
# table(seurat_obj$sub.cluster)
