# 读取并整理18个国家的多样性数据 - 每日版本
library(dplyr)
library(readr)
library(tidyr)
library(lubridate)
library(zoo)

# 定义文件路径
data_dir <- "/home/yuki/drift/outcome"

# 定义要分析的国家列表
countries <- c("China", "Japan", "South Korea", "Singapore", "India", 
               "USA", "Canada", "Mexico", "Australia", "New Zealand", 
               "Brazil", "Argentina", "United Kingdom", "Germany", 
               "France", "Italy", "Russia", "South Africa")

# 初始化数据列表
drift_data_list <- list()

# 循环读取每个国家的数据
for (country in countries) {
  file_path <- file.path(data_dir, paste0(country, "_RBD_p_distance.csv"))
  
  if (file.exists(file_path)) {
    # 读取CSV文件
    data <- read_csv(file_path, show_col_types = FALSE)
    
    # 添加国家列
    data$country <- country
    
    # 转换month列为日期格式
    data$month <- lubridate::ym(data$month)
    
    # 存储到列表中
    drift_data_list[[country]] <- data
    
    cat("已读取", country, "的数据：", nrow(data), "行\n")
  } else {
    cat("警告：找不到文件", file_path, "\n")
  }
}

# 合并所有数据
combined_drift_data <- bind_rows(drift_data_list)

# 保存合并后的数据
output_file <- "/home/Dong/Global_COVID/genomic_analysis/combined_drift_data_daily.csv"
write_csv(combined_drift_data, output_file)
cat("\n合并后的数据已保存到：", output_file, "\n")

# 将月度多样性数据转换为每日数据
cat("\n开始将月度数据转换为每日数据...\n")

# 使用线性插值转换为每日数据
cat("使用线性插值转换为每日数据\n")

# 为每个国家单独处理
daily_data_list <- list()

for (country_name in unique(combined_drift_data$country)) {
  cat("处理国家:", country_name, "\n")
  
  # 提取该国家的数据
  country_data <- combined_drift_data %>%
    filter(country == country_name) %>%
    arrange(month)
  
  if (nrow(country_data) < 2) {
    cat("  跳过", country_name, "- 数据点不足\n")
    next
  }
  
  # 创建完整的时间序列（从最早到最晚，按天）
  start_date <- floor_date(min(country_data$month), "day")
  end_date <- ceiling_date(max(country_data$month), "month") - days(1)
  
  # 生成所有天的日期
  all_days <- seq(from = start_date, to = end_date, by = "day")
  
  # 创建完整的时间框架
  daily_frame <- data.frame(
    date = all_days,
    country = country_name
  )
  
  # 线性插值核苷酸多样性
  drift_ts <- zoo(country_data$p_distance, country_data$month)
  drift_daily <- na.approx(drift_ts, xout = all_days, na.rm = FALSE, rule = 2)
  
  # 线性插值样本量（四舍五入为整数）
  sample_ts <- zoo(country_data$sample_size, country_data$month)
  sample_daily <- round(na.approx(sample_ts, xout = all_days, na.rm = FALSE, rule = 2))
  
  # 合并数据
  daily_frame$p_distance <- as.numeric(drift_daily)
  daily_frame$sample_size <- as.integer(sample_daily)
  daily_frame$year_month_day <- format(daily_frame$date, "%Y-%m-%d")
  
  daily_data_list[[country_name]] <- daily_frame
  
  cat("  完成", country_name, ":", length(all_days), "天的数据\n")
}

# 合并所有国家的每日数据
daily_combined <- bind_rows(daily_data_list)

# 保存长格式的每日数据
write_csv(daily_combined, "/home/Dong/Global_COVID/genomic_analysis/daily_drift_long.csv")
cat("保存长格式每日数据到：daily_drift_long.csv\n")

# 对核苷酸多样性进行归一化处理
cat("\n对核苷酸多样性进行归一化处理\n")
divmax <- max(daily_combined$p_distance, na.rm = TRUE)
daily_combined$p_distance_normalized <- daily_combined$p_distance / divmax

# 转换为宽格式 - 核苷酸多样性（归一化后）
cat("转换为宽格式...\n")
drift_daily_wide <- daily_combined %>%
  select(country, year_month_day, p_distance_normalized) %>%
  pivot_wider(
    names_from = year_month_day,
    values_from = p_distance_normalized,
    values_fill = 0
  ) %>%
  arrange(country)

# 转换为宽格式 - 样本量
sample_daily_wide <- daily_combined %>%
  select(country, year_month_day, sample_size) %>%
  pivot_wider(
    names_from = year_month_day,
    values_from = sample_size,
    values_fill = 0
  ) %>%
  arrange(country)

# 按照 countries 顺序重新排列行
drift_ordered <- drift_daily_wide[match(countries, drift_daily_wide$country), ]
drift_daily_wide <- drift_ordered

sample_ordered <- sample_daily_wide[match(countries, sample_daily_wide$country), ]
sample_daily_wide <- sample_ordered

# 保存宽格式数据
write_csv(drift_daily_wide, "/home/Dong/Global_COVID/genomic_analysis/drift_daily_wide.csv")
write_csv(sample_daily_wide, "/home/Dong/Global_COVID/genomic_analysis/drift_sample_size_daily_wide.csv")

cat("\n每日数据处理完成！\n")
cat("生成的文件：\n")
cat("- 原始合并数据：combined_drift_data_daily.csv\n")
cat("- 长格式每日数据：daily_drift_long.csv\n")
cat("- 宽格式每日多样性数据（归一化）：drift_daily_wide.csv\n")
cat("- 宽格式每日样本量数据：drift_sample_size_daily_wide.csv\n")

# 显示处理统计
cat("\n处理统计：\n")
cat("国家数量：", length(unique(daily_combined$country)), "\n")
cat("总天数范围：", min(daily_combined$date), "到", max(daily_combined$date), "\n")
cat("总数据点：", nrow(daily_combined), "\n")
cat("多样性数值范围（归一化前）：", round(min(daily_combined$p_distance, na.rm=TRUE), 4), 
    "到", round(max(daily_combined$p_distance, na.rm=TRUE), 4), "\n")
cat("多样性数值范围（归一化后）：", round(min(daily_combined$p_distance_normalized, na.rm=TRUE), 4), 
    "到", round(max(daily_combined$p_distance_normalized, na.rm=TRUE), 4), "\n")
