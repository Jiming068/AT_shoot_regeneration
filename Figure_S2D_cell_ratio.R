# Portable path configuration (override with environment variables).
project_dir <- Sys.getenv("SHOOT_REGEN_PROJECT_DIR", unset = ".")
output_root <- Sys.getenv("SHOOT_REGEN_OUTPUT_DIR", unset = file.path(project_dir, "02_results"))
external_data_dir <- Sys.getenv("SHOOT_REGEN_EXTERNAL_DATA_DIR", unset = file.path(project_dir, "external_data"))
work_dir <- Sys.getenv("SHOOT_REGEN_WORK_DIR", unset = ".")

library(tidyverse)
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
set.seed(128)

dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)
snRNA2 <- readRDS(file.path(project_dir, "01_obj/01_WT_snRNA/01_snRNA_res1.2_SIM4-SIM16.rds"))
Idents(snRNA2) <- 'merged_cluster'

mylevels =  c(
"C1","C2","C3","C4","C5","C6","C7",
"C8","C9","C10","C11","C12","C13","C14","C15",
"C16","C17","C18","C19","C20","C21","C22","C23")

Idents(snRNA2) <- factor(snRNA2@meta.data$merged_cluster, levels = mylevels)

pal <- c('#40b5c4','#647a39','#636363','#7c4273','#8d6d31','#9f9dd0','#3a3d7a','#e7989e',
  '#3382bc','#e1ebb5','#6baed6','#a2d49b','#f88d41','#e9cb95','#e6bb56','#bf9f38','#5684c3','#969696','#d7a0cc','#afd9ec', '#b6cd6c', '#c66fad','#e55828')

# FIG 2B  cell ratio change

write.csv(x, file = "root_derived_cells_per_dated4-16.csv", row.names = FALSE)

write.csv(x, file = "root_derived_cells_per_sampled4-16.csv", row.names = FALSE)

# What proportion of cells are in each cluster?
prop.table(table(Idents(snRNA2)))

# How does cluster membership vary by replicate?
write.csv(x, file = "root_derived_cells_eachCluster_per_dated4-16.csv")

x <- prop.table(table(Idents(snRNA2), snRNA2$date), margin = 2)
write.csv(x, file = "root_derived_prop_cells_eachCluster_per_dated4-16.csv")

# 提取细胞类型和时间信息
cell_data <- snRNA2@meta.data %>%
  select(celltype = "merged_cluster", time = "date")

# 计算每个时间点的细胞类型比率
x <- prop.table(table(Idents(snRNA2), snRNA2$date), margin = 2)

# 转换为数据框
df <- as.data.frame(x)

# 自定义细胞类型颜色

celltype_colors <- c(
  "C1" = "#40b5c4",
  "C2" = "#647a39",
  "C3" = "#636363",
  "C4" = "#7c4273",
  "C5" = "#8d6d31",
  "C6" = "#9f9dd0",
  "C7" = "#3a3d7a",
  "C8" = "#e7989e",
  "C9" = "#3382bc",
  "C10" = "#e1ebb5",
  "C11" = "#6baed6",
  "C12" = "#a2d49b",
  "C13" = "#f88d41",
  "C14" = "#e9cb95",
  "C15" = "#e6bb56",
  "C16" = "#bf9f38",
  "C17" = "#5684c3",
  "C18" = "#969696",
  "C19" = "#d7a0cc",
  "C20" = "#afd9ec",
  "C21" = "#b6cd6c",
  "C22" = "#c66fad",
  "C23" = "#e55828"
)

p <- ggplot(df, aes(x = Time, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "", y = "Fold change", fill = "Cluster") +
  scale_fill_manual(values = celltype_colors) +  # 自定义颜色
  theme_minimal()

pdf("Cell_ratio_change.pdf", width = 6, height = 5)
print(p)
dev.off()
