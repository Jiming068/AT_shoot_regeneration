# Portable path configuration (override with environment variables).
project_dir <- Sys.getenv("SHOOT_REGEN_PROJECT_DIR", unset = ".")
output_root <- Sys.getenv("SHOOT_REGEN_OUTPUT_DIR", unset = file.path(project_dir, "02_results"))
external_data_dir <- Sys.getenv("SHOOT_REGEN_EXTERNAL_DATA_DIR", unset = file.path(project_dir, "external_data"))
work_dir <- Sys.getenv("SHOOT_REGEN_WORK_DIR", unset = ".")

# ============================================================
# 单细胞转录组 - 细胞类群（cluster）相关性网络图
# 节点 = cluster，边 = 相关性高于阈值的cluster对
# ============================================================

library(Seurat)
library(igraph)
library(ggraph)
library(tidygraph)

# 【方式一】如果你已经有 Seurat 对象 seurat_obj，
seurat_obj <- readRDS(file.path(external_data_dir, "snRNA_subset_sct_dim19_res0.2.rds"))
Idents(seurat_obj) <- 'seurat_clusters'

# 数据处理
DefaultAssay(seurat_obj) <- "SCT"

# ============================================================
# 【第1步】计算cluster平均表达谱与相关性矩阵
# 与热图脚本一致，如已运行过热图脚本得到 cor_mat，可直接跳到第3步
# ============================================================

# ---------------------- 1a. 从Seurat对象计算（方式一）----------------------
# Idents(seurat_obj) <- "seurat_clusters"   # 改成你的cluster元数据列名
avg_exp <- AverageExpression(seurat_obj,
                              group.by = "seurat_clusters",
                              assays = "SCT",
                              slot = "data")
avg_exp_mat <- as.matrix(avg_exp$SCT)

if ("VariableFeatures" %in% ls("package:Seurat")) {
  hvg <- VariableFeatures(seurat_obj)
  hvg <- intersect(hvg, rownames(avg_exp_mat))
  if (length(hvg) > 0) avg_exp_mat <- avg_exp_mat[hvg, ]
}

# ---------------------- 2. 计算相关性矩阵 ----------------------
cor_mat <- cor(avg_exp_mat, method = "pearson")

cluster_order <- as.character(0:11)
cluster_order <- cluster_order[cluster_order %in% colnames(cor_mat)]
cor_mat <- cor_mat[cluster_order, cluster_order]

# ============================================================
# 【第3步】沿用之前统一的cluster配色方案
# ============================================================
cell_colors <- c(
  "0" = "#FFB74D", "8" = "#EF5350",
  "2" = "#81C784", "3" = "#66BB6A",
  "4" = "#BA68C8", "5" = "#CE93D8", "6" = "#F48FB1", "10" = "#EC407A",
  "1" = "#A1887F", "7" = "#8D6E63", "9" = "#795548", "11" = "#6D4C41"
)

# ---------------------- 4. 设置相关性阈值，构建边列表 ----------------------

# 只保留相关性高于阈值的cluster对作为连线，阈值可根据数据分布调整
threshold <- 0.9

edge_list <- data.frame()
for (i in 1:(nrow(cor_mat) - 1)) {
  for (j in (i + 1):ncol(cor_mat)) {
    r <- cor_mat[i, j]
    if (r >= threshold) {
      edge_list <- rbind(edge_list,
                          data.frame(from = rownames(cor_mat)[i],
                                     to = colnames(cor_mat)[j],
                                     weight = r))
    }
  }
}

# ---------------------- 5. 构建igraph/tidygraph网络对象 ----------------------
# 只保留在边列表中出现过的cluster（即至少有一条相关性≥阈值的连线），
# 完全没有达到阈值的孤立cluster将不会作为节点显示
connected_clusters <- unique(c(edge_list$from, edge_list$to))
connected_clusters <- cluster_order[cluster_order %in% connected_clusters]  # 保持原始顺序

if (length(connected_clusters) == 0) {
  stop("没有任何cluster对的相关性达到阈值 ", threshold, "，请适当降低threshold后重试")
}

nodes <- data.frame(name = connected_clusters)

graph_obj <- tbl_graph(nodes = nodes, edges = edge_list, directed = FALSE)

# ---------------------- 6. 绘制网络图 ----------------------
set.seed(42)  # 固定布局，保证每次运行图形一致

p <- ggraph(graph_obj, layout = "fr") +   # fr = Fruchterman-Reingold力导向布局
  geom_edge_link(aes(width = weight, alpha = weight),
                 color = "#D9432C") +
  scale_edge_width(range = c(0.3, 2.5), name = "correlation") +
  scale_edge_alpha(range = c(0.15, 0.8), guide = "none") +
  geom_node_point(aes(color = name), size = 12) +
  geom_node_text(aes(label = name), color = "white",
                 fontface = "bold", size = 4) +
  scale_color_manual(values = cell_colors, guide = "none") +
  labs(title = paste0("Cluster correlation (correlation \u2265 ", threshold, ")")) +
  theme_void(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right")

print(p)

# ---------------------- 7. 保存图片 ----------------------
ggsave("cluster_correlation_network.pdf", plot = p, width = 8, height = 7)
ggsave("cluster_correlation_network.png", plot = p, width = 8, height = 7, dpi = 300)

# 使用原始UMAP图的坐标

# ============================================================
# 单细胞转录组 - 细胞类群（cluster）相关性网络图
# 节点 = cluster，边 = 相关性高于阈值的cluster对
# ============================================================

library(Seurat)
library(igraph)
library(ggraph)
library(tidygraph)
library(dplyr)

# ============================================================
# 【第1步】计算cluster平均表达谱与相关性矩阵
# 与热图脚本一致，如已运行过热图脚本得到 cor_mat，可直接跳到第3步
# ============================================================

# ---------------------- 1a. 从Seurat对象计算（方式一）----------------------
# Idents(seurat_obj) <- "seurat_clusters"   # 改成你的cluster元数据列名
avg_exp <- AverageExpression(seurat_obj,
                              group.by = "seurat_clusters",
                              assays = "SCT",
                              slot = "data")
avg_exp_mat <- as.matrix(avg_exp$SCT)

if ("VariableFeatures" %in% ls("package:Seurat")) {
  hvg <- VariableFeatures(seurat_obj)
  hvg <- intersect(hvg, rownames(avg_exp_mat))
  if (length(hvg) > 0) avg_exp_mat <- avg_exp_mat[hvg, ]
}

# ---------------------- 1b. 或从CSV读取平均表达矩阵（方式二）----------------------
# avg_exp_mat <- read.csv("average_expression_matrix.csv",
#                          row.names = 1, check.names = FALSE)
# avg_exp_mat <- as.matrix(avg_exp_mat)

# ---------------------- 2b. 计算每个cluster在UMAP上的质心坐标 ----------------------
# 用于让网络图节点位置与UMAP图上的实际空间分布保持一致，
# 而不是用力导向算法重新布局（那样节点位置和UMAP图对不上）
umap_coords <- Embeddings(seurat_obj, reduction = "umap")  # 如果reduction名称不同，改成对应名字，如"UMAP"
cluster_labels <- as.character(Idents(seurat_obj))         # 需与上面AverageExpression用的分组一致

centroid_df <- data.frame(umap_coords, cluster = cluster_labels) %>%
  group_by(cluster) %>%
  summarise(x = mean(.data[[colnames(umap_coords)[1]]]),
            y = mean(.data[[colnames(umap_coords)[2]]]),
            .groups = "drop")
centroid_df <- as.data.frame(centroid_df)
rownames(centroid_df) <- centroid_df$cluster
# ---------------------- 2a. 计算相关性矩阵 ----------------------
cor_mat <- cor(avg_exp_mat, method = "spearman")

cluster_order <- as.character(0:11)
cluster_order <- cluster_order[cluster_order %in% colnames(cor_mat)]
cor_mat <- cor_mat[cluster_order, cluster_order]

# ============================================================
# 沿用之前统一的cluster配色方案
# ============================================================
cell_colors <- c(
  "0" = "#FFB74D", "8" = "#EF5350",
  "2" = "#81C784", "3" = "#66BB6A",
  "4" = "#BA68C8", "5" = "#CE93D8", "6" = "#F48FB1", "10" = "#EC407A",
  "1" = "#A1887F", "7" = "#8D6E63", "9" = "#795548", "11" = "#6D4C41"
)

# ---------------------- 4. 设置相关性阈值，构建边列表 ----------------------
# 只保留相关性高于阈值的cluster对作为连线，阈值可根据数据分布调整
threshold <- 0.88

edge_list <- data.frame()
for (i in 1:(nrow(cor_mat) - 1)) {
  for (j in (i + 1):ncol(cor_mat)) {
    r <- cor_mat[i, j]
    if (r >= threshold) {
      edge_list <- rbind(edge_list,
                          data.frame(from = rownames(cor_mat)[i],
                                     to = colnames(cor_mat)[j],
                                     weight = r))
    }
  }
}

# ---------------------- 5. 构建igraph/tidygraph网络对象 ----------------------
# 只保留在边列表中出现过的cluster（即至少有一条相关性≥阈值的连线），
# 完全没有达到阈值的孤立cluster将不会作为节点显示
connected_clusters <- unique(c(edge_list$from, edge_list$to))
connected_clusters <- cluster_order[cluster_order %in% connected_clusters]  # 保持原始顺序

if (length(connected_clusters) == 0) {
  stop("没有任何cluster对的相关性达到阈值 ", threshold, "，请适当降低threshold后重试")
}

# 节点坐标使用UMAP质心，保证网络图中的相对位置与UMAP图一致
nodes <- centroid_df[connected_clusters, c("cluster", "x", "y")]
colnames(nodes)[1] <- "name"

graph_obj <- tbl_graph(nodes = nodes, edges = edge_list, directed = FALSE)

# ---------------------- 6. 绘制网络图（使用UMAP质心作为节点坐标的manual布局）----------------------
layout_obj <- create_layout(graph_obj, layout = "manual",
                             x = nodes$x, y = nodes$y)

p <- ggraph(layout_obj) +
  geom_edge_link(aes(width = weight, alpha = weight),
                 color = "#D9432C") +
  scale_edge_width(range = c(0.3, 2.5), name = "correlation") +
  scale_edge_alpha(range = c(0.15, 0.8), guide = "none") +
  geom_node_point(aes(color = name), size = 12) +
  geom_node_text(aes(label = name), color = "white",
                 fontface = "bold", size = 4) +
  scale_color_manual(values = cell_colors, guide = "none") +
  labs(title = paste0("Cluster correlation (correlation \u2265 ", threshold, ")"),
       x = "UMAP_1", y = "UMAP_2") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right",
        panel.grid = element_blank())

print(p)

# ---------------------- 7. 保存图片 ----------------------
ggsave("cluster_correlation_network_spearman_2.pdf", plot = p, width = 8, height = 7)
ggsave("cluster_correlation_network_spearman_2.png", plot = p, width = 8, height = 7, dpi = 300)
