# 读取并整理18个国家的多样性数据
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
#countries <- c("China")

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
output_file <- "/home/Dong/Global_COVID/genomic_analysis/combined_drift_data1.csv"
write_csv(combined_drift_data, output_file)
cat("\n合并后的数据已保存到：", output_file, "\n")

# 将月度多样性数据转换为周度数据
cat("\n开始将月度数据转换为周度数据...\n")

# 读取现有的长格式数据
combined_data <- read_csv("/home/Dong/Global_COVID/genomic_analysis/combined_drift_data1.csv", show_col_types = FALSE)

# 使用线性插值转换为周度数据
cat("使用线性插值转换为周度数据\n")

# 为每个国家单独处理
weekly_data_list <- list()

for (country_name in unique(combined_data$country)) {
  cat("处理国家:", country_name, "\n")
  
  # 提取该国家的数据
  country_data <- combined_data %>%
    filter(country == country_name) %>%
    arrange(month)
  
  if (nrow(country_data) < 2) {
    cat("  跳过", country_name, "- 数据点不足\n")
    next
  }
  
  # 创建完整的时间序列（从最早到最晚，按周）
  start_date <- floor_date(min(country_data$month), "week")
  end_date <- ceiling_date(max(country_data$month), "week") - days(1)
  
  # 生成所有周的日期（每周一）
  all_weeks <- seq(from = start_date, to = end_date, by = "week")
  
  # 创建完整的时间框架
  weekly_frame <- data.frame(
    week = all_weeks,
    country = country_name
  )
  
  # 线性插值核苷酸多样性
  drift_ts <- zoo(country_data$p_distance, country_data$month)
  drift_weekly <- na.approx(drift_ts, xout = all_weeks, na.rm = FALSE, rule = 2)
  
  # 线性插值样本量（四舍五入为整数）
  sample_ts <- zoo(country_data$p_distance, country_data$month)
  sample_weekly <- round(na.approx(sample_ts, xout = all_weeks, na.rm = FALSE, rule = 2))
  
  # 合并数据
  weekly_frame$p_distance <- as.numeric(drift_weekly)
  weekly_frame$sample_size <- as.integer(sample_weekly)
  weekly_frame$year_week <- format(weekly_frame$week, "%Y-W%U")
  
  weekly_data_list[[country_name]] <- weekly_frame
}

# 合并所有国家的周度数据
weekly_combined <- bind_rows(weekly_data_list)

# 保存长格式的周度数据
write_csv(weekly_combined, "/home/Dong/Global_COVID/genomic_analysis/weekly_drift_long1.csv")


# 对核苷酸多样性进行归一化处理
cat("对核苷酸多样性进行归一化处理\n")
divmax <- max(weekly_combined$p_distance)
weekly_combined$p_distance <- weekly_combined$p_distance / divmax

# 转换为宽格式 - 核苷酸多样性
cat("\n转换为宽格式...\n")
drift_weekly_wide <- weekly_combined %>%
  select(country, year_week, p_distance) %>%
  pivot_wider(
    names_from = year_week,
    values_from = p_distance,
    values_fill = 0
  ) %>%
  arrange(country)

# 转换为宽格式 - 样本量
sample_weekly_wide <- weekly_combined %>%
  select(country, year_week, p_distance) %>%
  pivot_wider(
    names_from = year_week,
    values_from = p_distance,
    values_fill = 0
  ) %>%
  arrange(country)

# 按照 countries 顺序重新排列行
drift_ordered <- drift_weekly_wide[match(countries, drift_weekly_wide$country), ]
drift_weekly_wide <- drift_ordered

# 保存宽格式数据
write_csv(drift_weekly_wide, "/home/Dong/Global_COVID/genomic_analysis/drift_weekly_wide1.csv")
write_csv(sample_weekly_wide, "/home/Dong/Global_COVID/genomic_analysis/drift_sample_size_weekly_wide1.csv")






# 为每个国家单独处理
daily_data_list <- list()

for (country_name in unique(combined_data$country)) {
  cat("处理国家:", country_name, "\n")
  
  # 提取该国家的数据
  country_data <- combined_data %>%
    filter(country == country_name) %>%
    arrange(month)
  
  if (nrow(country_data) < 2) {
    cat("  跳过", country_name, "- 数据点不足\n")
    next
  }
  
  # 创建完整的时间序列（从最早到最晚，按天）
  start_date <- floor_date(min(country_data$month), "day")
  end_date <- ceiling_date(max(country_data$month), "day") - days(1)
  
  # 生成所有天的日期（每天）
  all_days <- seq(from = start_date, to = end_date, by = "day")
  
  # 创建完整的时间框架
  daily_frame <- data.frame(
    day = all_days,
    country = country_name
  )
  
  # 线性插值核苷酸多样性
  drift_ts <- zoo(country_data$p_distance, country_data$month)
  drift_daily <- na.approx(drift_ts, xout = all_days, na.rm = FALSE, rule = 2)
  
  # 线性插值样本量（四舍五入为整数）
  sample_ts <- zoo(country_data$p_distance, country_data$month)
  sample_daily <- round(na.approx(sample_ts, xout = all_days, na.rm = FALSE, rule = 2))
  
  # 合并数据
  daily_frame$p_distance <- as.numeric(drift_daily)
  daily_frame$sample_size <- as.integer(sample_daily)
  daily_frame$year_day <- format(daily_frame$day, "%Y-%m-%d")
  
  daily_data_list[[country_name]] <- daily_frame
}

# 合并所有国家的日度数据
daily_combined <- bind_rows(daily_data_list)

# 保存长格式的日度数据
write_csv(daily_combined, "/home/Dong/Global_COVID/genomic_analysis/daily_drift_long1.csv")

# 对核苷酸多样性进行归一化处理
cat("对核苷酸多样性进行归一化处理\n")
divmax <- max(daily_combined$p_distance, na.rm = TRUE)
daily_combined$p_distance_normalized <- daily_combined$p_distance / divmax

# 转换为宽格式 - 核苷酸多样性（归一化后）
cat("\n转换为宽格式...\n")
drift_daily_wide <- daily_combined %>%
  select(country, year_day, p_distance_normalized) %>%
  pivot_wider(
    names_from = year_day,
    values_from = p_distance_normalized,
    values_fill = 0
  ) %>%
  arrange(country)

# 按照 countries 顺序重新排列行
drift_ordered <- drift_daily_wide[match(countries, drift_daily_wide$country), ]
drift_daily_wide <- drift_ordered

# 保存宽格式数据
write_csv(drift_daily_wide, "/home/Dong/Global_COVID/genomic_analysis/drift_daily_wide1.csv")
write_csv(sample_daily_wide, "/home/Dong/Global_COVID/genomic_analysis/drift_sample_size_daily_wide1.csv")
