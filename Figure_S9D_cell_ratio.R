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
Idents(snRNA) <- "date"

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
            color = "black", size = 3.5, fontface = "bold") +
  scale_fill_gradient(
    low = "white",
    high = "#BCAAA4",
    name = "Percentage",
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    title = "Cell proportion",
    subtitle = "",
    x = "",
    y = "Cluster"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    panel.grid = element_blank()
  )

# 保存图像
pdf("cell_number_change_1.pdf", width = 6, height = 5)
print(plot)
dev.off()
