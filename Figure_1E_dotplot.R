# Portable path configuration (override with environment variables).
project_dir <- Sys.getenv("SHOOT_REGEN_PROJECT_DIR", unset = ".")
output_root <- Sys.getenv("SHOOT_REGEN_OUTPUT_DIR", unset = file.path(project_dir, "02_results"))
external_data_dir <- Sys.getenv("SHOOT_REGEN_EXTERNAL_DATA_DIR", unset = file.path(project_dir, "external_data"))
work_dir <- Sys.getenv("SHOOT_REGEN_WORK_DIR", unset = ".")

library(tidyverse)
library(Seurat)
library(patchwork)
library(pheatmap)
set.seed(128)

dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)

sub_obj <- readRDS(file.path(external_data_dir, "snRNA_subset_sct_dim19_res0.2.rds"))
Idents(sub_obj) <- 'seurat_clusters'

# 数据处理
DefaultAssay(sub_obj) <- "SCT"

subset_1 <- subset(x = sub_obj, idents = c( '0', '3', '2'), invert = FALSE)

DefaultAssay(subset_1) <- "SCT"

cell_order <- c('0', '3', '2')

original_idents <- Idents(subset_1)
Idents(subset_1) <- factor(Idents(subset_1), levels = cell_order)

# 读取标记基因数据（含第三列分类信息）
select_markers <- read.table('fig2c_genes_selected_7_12.txt', header = FALSE, sep = "\t",
                              stringsAsFactors = FALSE, blank.lines.skip = TRUE)
gene <- select_markers$V1
labels <- select_markers$V2
category <- select_markers$V3

# ---------------------- 检查基因是否都能在对象中找到 ----------------------
missing_genes <- gene[!gene %in% available_genes]

if (length(missing_genes) > 0) {
  warning("以下基因未在subset_1中找到，将被DotPlot自动跳过，请检查基因ID是否一致：\n",
          paste(missing_genes, collapse = ", "))
}

keep <- gene %in% available_genes
gene <- gene[keep]
labels <- labels[keep]
category <- category[keep]

# 分类分组顺序整体倒转
category_order <- rev(unique(category))

# ---------------------- 绘制DotPlot ----------------------
dp_default <- DotPlot(
    subset_1,
    features = gene,
    cols = c("#FFEBEE", "#E53935"),
    col.min = 0,
    col.max = 0.3,
    dot.min = 0,
    dot.scale = 6,
    cluster.idents = FALSE
  )

# ---------------------- 关键修正：直接在数据层面设定因子顺序 ----------------------
# 不再使用 scale_x_discrete(limits = ...) 来控制顺序，
# 而是直接把 features.plot 的因子水平设成想要的（反转后的）基因顺序，
# 这样facet分组和坐标轴显示用的是同一套顺序，不会互相冲突导致错位
dp_default$data$features.plot <- factor(dp_default$data$features.plot,
                                          levels = rev(gene))

gene_category_map <- setNames(category, gene)
dp_default$data$category <- gene_category_map[as.character(dp_default$data$features.plot)]
dp_default$data$category <- factor(dp_default$data$category, levels = category_order)

label_map <- setNames(labels, gene)

dp_ordered <- dp_default +
  scale_y_discrete(limits = cell_order) +
  scale_x_discrete(labels = function(x) label_map[x]) +   # 不再传limits，顺序已在数据里设定好
  coord_flip() +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  theme(
    axis.text.x = element_text(angle = 0),
    axis.text.y = element_text(size = 8, face = "italic"),
    strip.text.y = element_text(angle = 0, face = "bold", size = 8),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.spacing = unit(0.3, "lines")
  )

pdf("subset_markers_dotplot_ordered_category_reversed_3.pdf", width = 5.3, height = 4)
print(dp_ordered)
dev.off()
