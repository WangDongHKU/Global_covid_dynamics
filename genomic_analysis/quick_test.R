# =============================================================================
# VOC数据快速测试脚本
# 用于测试单个文件的处理逻辑
# =============================================================================

# 加载必要的包
library(readr)
library(dplyr)
library(stringr)

# 设置工作目录
setwd("/home/Dong/Global_COVID/genomic_analysis")

# 选择一个测试文件（中国）
test_file <- "China_VOC_timeseries_wide.csv"

cat("测试文件:", test_file, "\n")

# 读取数据
data <- read_csv(test_file, show_col_types = FALSE)
cat("原始数据维度:", nrow(data), "行 x", ncol(data), "列\n")
cat("VOC类型:", paste(data$VOC, collapse = ", "), "\n")

# 获取日期列
date_cols <- colnames(data)[-1]
cat("日期范围:", date_cols[1], "到", date_cols[length(date_cols)], "\n")
cat("总日期数:", length(date_cols), "\n")

# 统计全零列
zero_count <- 0
for (col in date_cols) {
  col_sum <- sum(data[[col]], na.rm = TRUE)
  if (col_sum == 0) {
    zero_count <- zero_count + 1
  }
}
cat("全零列数量:", zero_count, "\n")

# 创建VOC Original行
voc_original_values <- c("VOC Original")
for (col in date_cols) {
  col_sum <- sum(data[[col]], na.rm = TRUE)
  voc_original_values <- c(voc_original_values, ifelse(col_sum == 0, 1, 0))
}

# 创建完整数据
voc_original_row <- data.frame(matrix(voc_original_values, nrow = 1))
colnames(voc_original_row) <- colnames(data)
data_with_original <- rbind(data, voc_original_row)

cat("添加VOC Original后维度:", nrow(data_with_original), "行 x", ncol(data_with_original), "列\n")

# 计算比例
prop_data <- data_with_original
for (col in date_cols) {
  prop_data[[col]] <- as.numeric(prop_data[[col]])
  col_total <- sum(prop_data[[col]], na.rm = TRUE)
  if (col_total > 0) {
    prop_data[[col]] <- prop_data[[col]] / col_total
  } else {
    prop_data[[col]] <- 0
  }
}

# 验证几列的比例和
cat("\n验证比例计算:\n")
sample_cols <- c(date_cols[1], date_cols[50], date_cols[100], date_cols[length(date_cols)])
for (col in sample_cols) {
  if (col %in% colnames(prop_data)) {
    col_sum <- sum(prop_data[[col]], na.rm = TRUE)
    cat(sprintf("%s: 比例和 = %.6f\n", col, col_sum))
  }
}

# 保存文件
count_file <- "China_VOC_with_original.csv"
prop_file <- "China_VOC_proportions.csv"

write_csv(data_with_original, count_file)
write_csv(prop_data, prop_file)

cat("\n文件已保存:\n")
cat("- 计数文件:", count_file, "\n")
cat("- 比例文件:", prop_file, "\n")

cat("\n测试完成！\n")
