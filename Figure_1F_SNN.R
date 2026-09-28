# Portable path configuration (override with environment variables).
project_dir <- Sys.getenv("SHOOT_REGEN_PROJECT_DIR", unset = ".")
output_root <- Sys.getenv("SHOOT_REGEN_OUTPUT_DIR", unset = file.path(project_dir, "02_results"))
external_data_dir <- Sys.getenv("SHOOT_REGEN_EXTERNAL_DATA_DIR", unset = file.path(project_dir, "external_data"))
work_dir <- Sys.getenv("SHOOT_REGEN_WORK_DIR", unset = ".")

# ============================================================
# 基于局部Shared-Nearest-Neighbor (SNN) 图计算RC state间连接性
# 并可视化为RC-level连接图
#
# 思路（类似PAGA的cluster间连接性抽象）：
# 1. 取Seurat对象中FindNeighbors()已经算好的SNN图（细胞x细胞的邻接矩阵）
# 2. 统计不同RC state（subcluster）之间的细胞在SNN图上的连接强度
# 3. 按cluster大小做归一化，得到RC state间的"连接性得分"
# 4. 只保留连接性得分超过阈值的RC state对，作为图上的边
# 5. 节点位置使用UMAP质心坐标，节点颜色沿用subcluster_colors
#
# 本脚本承接前一步 subset_and_subcluster_recolor.R 得到的 sub_obj
# ============================================================
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)

library(Seurat)
library(igraph)
library(ggraph)
library(tidygraph)
library(Matrix)
library(dplyr)
library(ggplot2)

# ============================================================
# 【第1步】确认RC state（subcluster）标签与SNN图
# ============================================================

sub_obj <- readRDS(file.path(external_data_dir, "snRNA_subset_sct_dim19_res0.2.rds"))
Idents(sub_obj) <- 'seurat_clusters'

rc_labels <- as.character(sub_obj$seurat_clusters)

rc_order <- sort(unique(rc_labels))   # RC state顺序，如"0","1","2"...

# ---------------------- 1a. 提取SNN图（细胞间局部邻接关系） ----------------------
# FindNeighbors()默认会同时生成knn图和snn图，存储在sub_obj@graphs里
# 图名称通常为 "<assay>_snn"，SCTransform后assay默认是"SCT"
print(names(sub_obj@graphs))   # 运行后确认实际的图名称

snn_graph_name <- paste0("SCT_snn")   # 例如 "SCT_snn"
snn_mat <- sub_obj@graphs[[snn_graph_name]]                # 稀疏矩阵，细胞x细胞

# 确保细胞顺序一致
snn_mat <- snn_mat[colnames(sub_obj), colnames(sub_obj)]

# ============================================================
# 【第2步】统计RC state间的SNN连接强度
# ============================================================

n_rc <- length(rc_order)
conn_sum <- matrix(0, nrow = n_rc, ncol = n_rc,
                    dimnames = list(rc_order, rc_order))
n_cells_rc <- table(rc_labels)[rc_order]

# 将稀疏SNN矩阵转换为便于按块求和的triplet格式，逐个RC state对累加边权重
snn_triplet <- summary(as(snn_mat, "TsparseMatrix"))   # i, j, x (权重)

# 将细胞索引映射为RC state标签，便于分组求和
snn_triplet$rc_i <- rc_labels[cell_names[snn_triplet$i]]
snn_triplet$rc_j <- rc_labels[cell_names[snn_triplet$j]]

# 按RC state对累加边权重之和（包含state内部与state间的连接）
pair_sum <- snn_triplet %>%
  group_by(rc_i, rc_j) %>%
  summarise(total_weight = sum(x), .groups = "drop")

for (k in 1:nrow(pair_sum)) {
  ri <- pair_sum$rc_i[k]
  rj <- pair_sum$rc_j[k]
  conn_sum[ri, rj] <- conn_sum[ri, rj] + pair_sum$total_weight[k]
}

# ---------------------- 2a. 按cluster大小归一化，得到连接性得分 ----------------------
# 归一化方式：观测到的连接权重之和 / 两个state可能的细胞对数量(n_i * n_j)
# 这样可以避免大cluster因为细胞多而"看起来"连接性更强的偏差
connectivity <- matrix(0, nrow = n_rc, ncol = n_rc,
                        dimnames = list(rc_order, rc_order))

for (ri in rc_order) {
  for (rj in rc_order) {
    if (ri == rj) next
    possible_pairs <- as.numeric(n_cells_rc[ri]) * as.numeric(n_cells_rc[rj])
    connectivity[ri, rj] <- conn_sum[ri, rj] / possible_pairs
  }
}

# ============================================================
# 【第3步】设置连接性阈值，构建RC state间的边列表
# ============================================================
# 阈值需要结合connectivity矩阵的实际分布来定，先看一下分布再决定
summary(connectivity[upper.tri(connectivity)])

threshold <- quantile(connectivity[upper.tri(connectivity)], 0.6)  # 示例：取上四分位数作为阈值
# 也可以直接设置固定数值，例如：
# threshold <- 0.02

edge_list <- data.frame()
for (i in 1:(n_rc - 1)) {
  for (j in (i + 1):n_rc) {
    ri <- rc_order[i]; rj <- rc_order[j]
    w <- max(connectivity[ri, rj], connectivity[rj, ri])  # 矩阵理论上对称，取较大值保险
    if (w >= threshold) {
      edge_list <- rbind(edge_list, data.frame(from = ri, to = rj, weight = w))
    }
  }
}

# 只保留在边列表中出现过的RC state（没有任何达到阈值连接的孤立RC state不显示）
connected_rc <- unique(c(edge_list$from, edge_list$to))
connected_rc <- rc_order[rc_order %in% connected_rc]

if (length(connected_rc) == 0) {
  stop("没有任何RC state对的连接性达到阈值 ", round(threshold, 4), "，请适当降低threshold后重试")
}

# ============================================================
# 【第4步】节点坐标使用RC state在UMAP上的质心，保持空间位置关系
# ============================================================
umap_coords <- Embeddings(sub_obj, reduction = "umap")
centroid_df <- data.frame(umap_coords, rc = rc_labels[rownames(umap_coords)]) %>%
  group_by(rc) %>%
  summarise(x = mean(.data[[colnames(umap_coords)[1]]]),
            y = mean(.data[[colnames(umap_coords)[2]]]),
            .groups = "drop")
centroid_df <- as.data.frame(centroid_df)
rownames(centroid_df) <- centroid_df$rc

nodes <- centroid_df[connected_rc, c("rc", "x", "y")]
colnames(nodes)[1] <- "name"

# ============================================================
# 【第5步】RC state配色方案（沿用之前浅暖色系subcluster配色）
# ============================================================
subcluster_colors <- c(
 # 组1: 0和8群使用对比配色（暖色系中的对比）
  "0" = "#FFB74D",  # 柔和的橙黄色
  "8" = "#EF5350",  # 柔和的珊瑚红 - 与0形成暖色系内对比

  # 组2: 2和3群使用同一色系但不要太相似（柔和绿色系）
  "2" = "#81C784",  # 柔和的浅绿色
  "3" = "#66BB6A",  # 柔和的绿色 - 与2同系但更饱和

  # 组3: 4,5,6,10群使用补充配色（柔和紫色/粉色系）
  "4" = "#BA68C8",  # 柔和的淡紫色
  "5" = "#CE93D8",  # 柔和的薰衣草紫
  "6" = "#F48FB1",  # 柔和的粉红色
  "10" = "#EC407A",  # 柔和的玫瑰红

  # 组4: 1,7,9,11群使用相近色系但不要太相似（柔和棕色系）
  "1" = "#A1887F",  # 柔和的灰棕色
  "7" = "#8D6E63",  # 柔和的棕色
  "9" = "#795548",  # 柔和的深棕色
  "11" = "#6D4C41" # 柔和的巧克力棕色
)
# 如果RC state数量或编号与上面不一致，请根据 rc_order 的实际结果调整

# ============================================================
# 【第6步】构建网络对象并绘制RC-level连接图
# ============================================================
graph_obj <- tbl_graph(nodes = nodes, edges = edge_list, directed = FALSE)

layout_obj <- create_layout(graph_obj, layout = "manual",
                             x = nodes$x, y = nodes$y)

p_rc_network <- ggraph(layout_obj) +
  geom_edge_link(aes(width = weight, alpha = weight),
                 color = "#D9432C") +
  scale_edge_width(range = c(0.3, 3), name = "Connectivity") +
  scale_edge_alpha(range = c(0.2, 0.85), guide = "none") +
  geom_node_point(aes(color = name), size = 12) +
  geom_node_text(aes(label = name), color = "black",
                 fontface = "bold", size = 4) +
  scale_color_manual(values = subcluster_colors, guide = "none") +
  labs(title = paste0("SNN connectivity"),
       subtitle = "",
       x = "UMAP_1", y = "UMAP_2") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "right",
        panel.grid = element_blank())

print(p_rc_network)

# ---------------------- 保存图片 ----------------------
ggsave("RC_state_connectivity_network.pdf", plot = p_rc_network, width = 8, height = 7)
ggsave("RC_state_connectivity_network.png", plot = p_rc_network, width = 8, height = 7, dpi = 300)

# ============================================================
# 【可选】导出RC state间连接性矩阵，便于查看具体数值或用于其他分析
# ============================================================
write.csv(connectivity, "RC_state_connectivity_matrix.csv")

# ============================================================
# 【可选】叠加在原始UMAP细胞散点图之上，兼顾细胞分布和RC状态连接关系
# ============================================================
umap_cells_df <- data.frame(umap_coords, rc = rc_labels[rownames(umap_coords)])
 umap_cells_df$rc <- factor(umap_cells_df$rc, levels = rc_order)

 p_overlay <- ggplot() +
   geom_point(data = umap_cells_df, aes(x = UMAP_1, y = UMAP_2, color = rc),
              size = 0.3, alpha = 0.3) +
   geom_segment(data = edge_list %>%
                  left_join(nodes, by = c("from" = "name")) %>%
                  left_join(nodes, by = c("to" = "name"), suffix = c("", "_end")),
                aes(x = x, y = y, xend = x_end, yend = y_end,
                    linewidth = weight, alpha = weight),
                color = "#7F0000") +
   geom_point(data = nodes, aes(x = x, y = y, fill = name),
              shape = 21, size = 8, color = "white", stroke = 1) +
   geom_text(data = nodes, aes(x = x, y = y, label = name),
             color = "black", fontface = "bold", size = 3.5) +
   scale_color_manual(values = subcluster_colors, guide = "none") +
   scale_fill_manual(values = subcluster_colors, name = "RC state") +
   scale_linewidth(range = c(0.3, 2.5), guide = "none") +
   scale_alpha(range = c(0.2, 0.85), guide = "none") +
   labs(title = "SNN connectivity") +
   theme_minimal(base_size = 14) +
   theme(plot.title = element_text(hjust = 0.5, face = "bold"),
         panel.grid = element_blank())
 print(p_overlay)
 ggsave("RC_state_connectivity_over_umap.pdf", plot = p_overlay, width = 8, height = 7)

# ============================================================
# 【可选】阈值敏感性分析：不同阈值下保留多少条RC state连接边
# ============================================================
for (q in c(0.5, 0.6, 0.7, 0.75, 0.8, 0.9)) {
   t <- quantile(connectivity[upper.tri(connectivity)], q)
  n_edges <- sum(connectivity[upper.tri(connectivity)] >= t)
   cat("分位数", q, "(阈值=", round(t, 4), ") -> 保留边数:", n_edges, "\n")
 }
