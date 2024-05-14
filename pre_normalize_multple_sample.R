# 准备----------------------------------------------------------------

## library
library(limma)
library(dplyr)

fs::dir_create("2-GEO校正")
rm(list = ls())

# 预处理-------------------------------------------------

# 获取文件列表，包含全路径名
file_list <- list.files(
  path = "./2-GEO下载/",
  pattern = "^(exprs|group)_[0-9]+\\.csv$",
  full.names = TRUE
)

# 初始化列表以存储数据框
exprs_list <- list()
group_list <- list()

# 循环读取文件
for (file_name in file_list) {
  code <- gsub(".*_(\\d+)\\.csv", "\\1", basename(file_name))
  code <- paste0("GSE", code)
  # 根据文件名决定是表达矩阵文件还是分组文件
  if (grepl("exprs_", basename(file_name))) {
    exprs_list[[code]] <- data.table::fread(file_name)
  } else if (grepl("group_", basename(file_name))) {
    group_list[[code]] <- read.csv(file_name)
  }
}

# 处理基因名重复的行
# TODO:前面已经处理了，但是还是有重复的
exprs_list <- lapply(exprs_list, function(x) {
  x <- aggregate(. ~ V1, data = x, mean) # 强大的分组汇总函数
  return(x)
})


# 处理列名V1转换为行名
exprs_list <- lapply(exprs_list, function(x) {
  x <- tibble::column_to_rownames(x, var = "V1")
  return(x)
})

# 校正----------------------------------------------

## 是否log2（value+1）转换
log2_transform <- function(matrix) { # 来自GEO官网
  print(paste0(c("最小值：", "最大值："), range(matrix)))
  qx <- as.numeric(quantile(matrix,
    c(0., 0.25, 0.5, 0.75, 0.99, 1.0),
    na.rm = TRUE
  ))
  log_c <- (qx[5] > 100) ||
    (qx[6] - qx[1] > 50 && qx[2] > 0)
  if (log_c) {
    matrix[which(matrix <= 0)] <- NaN
    matrix_log <- log2(matrix + 1)
    print("转换完成")
  } else {
    matrix_log <- matrix
    print("不需要转换")
  }
  return(matrix_log)
}
## log化
matrix_list <- lapply(exprs_list, log2_transform)

# 去批次--------------------------------------------

## 去批次前先合并
matrix_list <- lapply(matrix_list, function(x) {
  x <- tibble::rownames_to_column(x, var = "rownames")
  return(x)
})

combined_matrix <- matrix_list %>%
  purrr::reduce(inner_join, by = "rownames") %>%
  tibble::column_to_rownames(var = "rownames")

combined_group <- dplyr::bind_rows(group_list, .id = "project") %>%
  dplyr::arrange(desc(group)) # Normal在前，Disorder在后

## TODO:删除离谱的样本
# combined_group <- combined_group[!grepl("GSM2114232", combined_group$sample), ]

## 确保列顺序一致
combined_matrix <- combined_matrix[, combined_group$sample]


## sva中batch指定批次，mod指定保护的差异

combined_group$group <- factor(
  combined_group$group,
  levels = c("Normal", "Disorder")
)

mod <- model.matrix(~ 1 + group, data = combined_group)
batch <- combined_group$project
sva_batch_matrix <- sva::ComBat(combined_matrix, batch = batch, mod = mod)
# sva_batch_matrix <- sva::ComBat(combined_matrix, batch = batch)

# limma_batch_matrix <- limma::removeBatchEffect(combined_matrix, batch = batch) # nolint

## normalization ----
## limma包normalizeBetweenArrays
## 先去批次再归一化
batch_normal_matrix <- normalizeBetweenArrays(sva_batch_matrix) %>%
  data.frame()


# 保存--------------------------------------------
write.csv(batch_normal_matrix, "2-GEO校正/batch_normal_matrix.csv")
write.csv(combined_group, "2-GEO校正/combined_group.csv", row.names = FALSE)

## 可视化 ----
## 箱图
boxplot(combined_matrix, outline = FALSE, las = 2) ## 看看样本间是否有批次
boxplot(sva_batch_matrix, outline = FALSE, las = 2)
boxplot(batch_normal_matrix, outline = FALSE, las = 2)
dev.off()
## 密度图
limma::plotDensities(combined_matrix, legend = FALSE)
limma::plotDensities(sva_batch_matrix, legend = FALSE)
limma::plotDensities(batch_normal_matrix, legend = FALSE)

# PCA图  fviz_pca_ind
pca_before <- prcomp(t(combined_matrix),
  scale = TRUE
)
pca_after_batch <- prcomp(t(sva_batch_matrix),
  scale = TRUE
)
pca_after_normalize <- prcomp(t(batch_normal_matrix))

group_col <- c("#5E86C1", "#f59292")

library(plotly)
library(factoextra)
## 这是ggplot对象的
pca_before_p <- fviz_pca_ind(pca_before,
  title = "Before",
  label = "none",
  axes = c(1, 2),
  habillage = combined_group$project, # 按分组标注
  palette = group_col,
  addEllipses = TRUE, # 添加椭圆圈圈
  ellipse.level = 0.95 # 95％置信区间
)
pca_before_p
## 仅去批次
fviz_pca_ind(pca_after_batch,
  title = "After",
  label = "none",
  axes = c(1, 2),
  habillage = combined_group$project, # 按分组标注
  palette = group_col,
  addEllipses = TRUE, # 添加椭圆圈圈
  ellipse.level = 0.95 # 95％置信区间
)
## 去批次后再归一化的
pca_after_p <- fviz_pca_ind(pca_after_normalize,
  title = "After",
  label = "none",
  axes = c(1, 2),
  habillage = combined_group$project, # 按分组标注
  palette = group_col,
  addEllipses = TRUE, # 添加椭圆圈圈
  ellipse.level = 0.95 # 95％置信区间
)
pca_after_p
## 样本PCA
pca_sample <- fviz_pca_ind(pca_after_normalize,
  title = "After",
  label = "none",
  axes = c(1, 2),
  habillage = combined_group$sample, # 按样本标注
) +
  theme(legend.position = "none") # 不显示图注
ggplotly(pca_sample)
## 分组PCA
pca_group <- fviz_pca_ind(pca_after_normalize,
  title = "After",
  label = "none",
  axes = c(1, 2),
  habillage = combined_group$group, # 按分组标注
  palette = group_col,
  addEllipses = TRUE, # 添加椭圆圈圈
  ellipse.level = 0.95 # 95％置信区间
)
pca_group

# 美化 --------------------------------------------
library(ggplot2)
mytheme <-
  theme(
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"), # 图像
    plot.title = element_text(size = 7)
  ) +
  theme(
    panel.background = element_blank(), # 面板
    panel.grid = element_blank(),
    # panel.borde = element_rect(fill = NA, linewidth = 0.75 * 0.47)
  ) + # 添加外框
  theme(
    axis.line = element_line(size = 0.75 * 0.47), # 坐标轴
    axis.text = element_text(size = 6, color = "black"), # 坐标文字为黑色
    axis.title = element_text(size = 6),
    axis.ticks = element_line(size = 0.75 * 0.47)
  ) +
  theme(
    legend.key = element_rect(fill = "white"),
    legend.key.size = unit(c(0.3, 0.3), "cm"), # 图注
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.margin = margin(),
    legend.box.margin = margin(),
    legend.box.spacing = unit(0, "cm"),
    legend.background = element_blank(), legend.spacing = unit(0, "cm"),
    legend.box.background = element_blank()
  ) +
  theme(
    panel.grid = element_line(
      colour = "grey90",
      size = 0.75 * 0.47, linetype = 1
    ),
    plot.title = element_text(
      size = 7, hjust = 0, vjust = 0.5,
      margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
    ),
    legend.position = "top", # 图注位置（上）
    legend.direction = "horizontal",
    legend.margin = margin(0, 0, 0, 0, "cm"),
    legend.title = element_blank(),
    legend.key.size = unit(c(0.15, 0.15), "cm")
  )

pca_before_p + mytheme
ggsave("2-GEO校正/C-PCA校正前.pdf", width = 6, height = 6, units = "cm")

pca_after_p + mytheme
ggsave("2-GEO校正/D-PCA校正后.pdf", width = 6, height = 6, units = "cm")

pca_sample
ggsave("2-GEO校正/E-PCA样本.pdf")

pca_group + mytheme
ggsave("2-GEO校正/F-PCA分组.pdf", width = 6, height = 6, units = "cm")
