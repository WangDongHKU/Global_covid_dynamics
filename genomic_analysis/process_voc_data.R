# =============================================================================
# VOC时间序列数据处理脚本
# 功能：为18个国家的VOC数据添加VOC Original行并计算比例
# 作者：AI Assistant
# 日期：2025年
# =============================================================================

# 加载必要的包
required_packages <- c("readr", "dplyr", "stringr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 设置工作目录
setwd("/home/Dong/Global_COVID/genomic_analysis")

# 获取所有VOC时间序列文件
voc_files <- list.files(pattern = "*_VOC_timeseries_wide.csv", full.names = TRUE)
cat(paste(rep("=", 60), collapse=""), "\n")
cat("VOC时间序列数据处理开始\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("找到", length(voc_files), "个VOC文件:\n")
for (i in seq_along(voc_files)) {
  cat(sprintf("%2d. %s\n", i, basename(voc_files[i])))
}

# 处理单个VOC文件的函数
process_voc_file <- function(file_path) {
  country_name <- str_replace(basename(file_path), "_VOC_timeseries_wide.csv", "")
  cat("\n", paste(rep("=", 50), collapse=""), "\n")
  cat("处理国家:", country_name, "\n")
  cat("文件:", basename(file_path), "\n")
  
  # 读取数据
  tryCatch({
    data <- read_csv(file_path, show_col_types = FALSE)
  }, error = function(e) {
    cat("错误：无法读取文件", file_path, "\n")
    stop(e)
  })
  
  # 验证数据结构
  if (ncol(data) < 2) {
    stop("数据文件格式错误：列数不足")
  }
  
  if (nrow(data) != 5) {
    cat("警告：期望5个VOC类型，实际找到", nrow(data), "个\n")
  }
  
  # 获取日期列（除了第一列VOC列）
  date_cols <- colnames(data)[-1]
  cat("日期范围:", date_cols[1], "到", date_cols[length(date_cols)], "\n")
  cat("总日期数:", length(date_cols), "\n")
  
  # 统计全零列数量
  zero_columns <- 0
  voc_original_values <- c("VOC Original")
  
  # 对每个日期列检查是否所有VOC都为0
  for (col in date_cols) {
    # 计算该列的总和（排除VOC列名）
    col_sum <- sum(data[[col]], na.rm = TRUE)
    
    # 如果该列所有数据都为0
    if (col_sum == 0) {
      zero_columns <- zero_columns + 1
      
      # 判断是否为2023年之后（使用正确的日期比较）
      is_after_2023 <- grepl("^2022-|^2023-|^2024-|^2025-", col)
      
      if (is_after_2023) {
        # 2022年之后：将VOC Omicron设为1，VOC Original为0
        data[data$VOC == "VOC Omicron", col] <- 1
        voc_original_values <- c(voc_original_values, 0)
        cat("调整", col, ": VOC Omicron = 1 (2022年后零列)\n")
      } else {
        # 2022年前：VOC Original = 1
        voc_original_values <- c(voc_original_values, 1)
      }
    } else {
      voc_original_values <- c(voc_original_values, 0)
    }
  }
  
  
  # 创建VOC Original行的数据框
  voc_original_row <- data.frame(matrix(voc_original_values, nrow = 1))
  colnames(voc_original_row) <- colnames(data)
  
  # 将VOC Original行添加到数据末尾
  data_with_original <- rbind(data, voc_original_row)
  
  # 计算比例数据
  prop_data <- data_with_original
  
  # 对每个日期列计算比例
  for (col in date_cols) {
    # 将列转换为数值型
    prop_data[[col]] <- as.numeric(prop_data[[col]])
    
    # 计算该列的总和
    col_total <- sum(prop_data[[col]], na.rm = TRUE)
    
    # 如果总和不为0，计算比例；否则保持为0
    if (col_total > 0) {
      prop_data[[col]] <- prop_data[[col]] / col_total
    } else {
      prop_data[[col]] <- 0  # 确保全零列的比例也为0
    }
  }
  
  # 保存添加了VOC Original的原始计数数据
  count_output_file <- str_replace(file_path, "_VOC_timeseries_wide.csv", "_VOC_with_original.csv")
  write_csv(data_with_original, count_output_file)
  cat("✓ 保存计数数据:", basename(count_output_file), "\n")
  
  # 保存比例数据
  prop_output_file <- str_replace(file_path, "_VOC_timeseries_wide.csv", "_VOC_proportions.csv")
  write_csv(prop_data, prop_output_file)
  cat("✓ 保存比例数据:", basename(prop_output_file), "\n")
  
  # 验证比例数据（随机检查几列的比例和是否为1）
  sample_cols <- sample(date_cols, min(5, length(date_cols)))
  for (col in sample_cols) {
    col_sum <- sum(prop_data[[col]], na.rm = TRUE)
    if (abs(col_sum - 1.0) > 0.001 && col_sum > 0) {
      cat("警告：", col, "列的比例和不等于1，实际为", col_sum, "\n")
    }
  }
  
  return(list(
    country = country_name,
    original_file = basename(file_path),
    count_file = basename(count_output_file),
    proportion_file = basename(prop_output_file),
    voc_count = nrow(data_with_original),
    date_count = length(date_cols),
    zero_columns = zero_columns
  ))
}

# =============================================================================
# 主处理循环
# =============================================================================

cat("\n", "="*60, "\n")
cat("开始批量处理所有文件...\n")
cat(paste(rep("=", 60), collapse=""), "\n")

results <- list()
failed_files <- character()

for (file in voc_files) {
  tryCatch({
    result <- process_voc_file(file)
    results[[length(results) + 1]] <- result
  }, error = function(e) {
    cat("❌ 处理文件", basename(file), "时出错:", e$message, "\n")
    failed_files <<- c(failed_files, basename(file))
  })
}

# =============================================================================
# 处理结果摘要
# =============================================================================

cat("\n", "="*60, "\n")
cat("处理结果摘要\n")
cat(paste(rep("=", 60), collapse=""), "\n")

if (length(results) > 0) {
  cat("✅ 成功处理的文件 (", length(results), "个):\n\n")
  
  for (i in seq_along(results)) {
    result <- results[[i]]
    cat(sprintf("%2d. %s\n", i, result$country))
    cat(sprintf("    原始文件: %s\n", result$original_file))
    cat(sprintf("    VOC数量: %d (包括VOC Original)\n", result$voc_count))
    cat(sprintf("    日期数量: %d\n", result$date_count))
    cat(sprintf("    全零列数: %d\n", result$zero_columns))
    cat(sprintf("    计数文件: %s\n", result$count_file))
    cat(sprintf("    比例文件: %s\n", result$proportion_file))
    cat("\n")
  }
  
  # 统计信息
  total_dates <- results[[1]]$date_count
  total_zero_cols <- sum(sapply(results, function(x) x$zero_columns))
  
  cat("📊 统计信息:\n")
  cat(sprintf("   - 处理的国家数量: %d\n", length(results)))
  cat(sprintf("   - 时间序列长度: %d 个时间点\n", total_dates))
  cat(sprintf("   - 总全零列数: %d\n", total_zero_cols))
  cat(sprintf("   - 生成的文件总数: %d (计数文件 + 比例文件)\n", length(results) * 2))
}

if (length(failed_files) > 0) {
  cat("❌ 处理失败的文件 (", length(failed_files), "个):\n")
  for (file in failed_files) {
    cat("   -", file, "\n")
  }
}

cat("\n", "="*60, "\n")
cat("处理完成！\n")
cat("成功:", length(results), "个文件")
if (length(failed_files) > 0) {
  cat("，失败:", length(failed_files), "个文件")
}
cat("\n")
cat(paste(rep("=", 60), collapse=""), "\n")

# =============================================================================
# 数据验证示例
# =============================================================================

if (length(results) > 0) {
  cat("\n📋 数据验证示例（以第一个成功处理的文件为例）:\n")
  
  first_result <- results[[1]]
  
  # 读取处理后的文件进行验证
  count_file_path <- paste0(str_replace(first_result$original_file, "_VOC_timeseries_wide.csv", ""), "_VOC_with_original.csv")
  prop_file_path <- paste0(str_replace(first_result$original_file, "_VOC_timeseries_wide.csv", ""), "_VOC_proportions.csv")
  
  tryCatch({
    count_data <- read_csv(count_file_path, show_col_types = FALSE)
    prop_data <- read_csv(prop_file_path, show_col_types = FALSE)
    
    cat("国家:", first_result$country, "\n")
    cat("计数文件形状:", nrow(count_data), "行 x", ncol(count_data), "列\n")
    cat("比例文件形状:", nrow(prop_data), "行 x", ncol(prop_data), "列\n")
    
    # 显示VOC类型
    cat("VOC类型:", paste(count_data$VOC, collapse = ", "), "\n")
    
    # 验证几个日期列的比例和
    date_cols <- colnames(prop_data)[-1]
    sample_dates <- sample(date_cols, min(3, length(date_cols)))
    
    cat("比例验证（随机选择的日期）:\n")
    for (date in sample_dates) {
      col_sum <- sum(prop_data[[date]], na.rm = TRUE)
      cat(sprintf("  %s: 比例和 = %.6f\n", date, col_sum))
    }
    
  }, error = function(e) {
    cat("验证时出错:", e$message, "\n")
  })
}



