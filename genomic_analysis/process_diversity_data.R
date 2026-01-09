# 读取并整理18个国家的多样性数据
library(dplyr)
library(readr)
library(tidyr)
library(readr)
library(tidyr)
# 定义文件路径
data_dir <- "/home/yuki/diversity/outcome"

# 定义要分析的国家列表
countries <- c("China", "Japan", "South Korea", "Singapore", "India", 
               "USA", "Canada", "Mexico", "Australia", "New Zealand", 
               "Brazil", "Argentina", "United Kingdom", "Germany", 
               "France", "Italy", "Russia", "South Africa")

# 初始化数据列表
diversity_data_list <- list()

# 循环读取每个国家的数据
for (country in countries) {
  file_path <- file.path(data_dir, paste0(country, "_diversity_results.csv"))
  
  if (file.exists(file_path)) {
    # 读取CSV文件
    data <- read_csv(file_path, show_col_types = FALSE)
    
    # 添加国家列
    data$country <- country
    
    # 转换month列为日期格式
    data$month <- as.Date(data$month)
    
    # 存储到列表中
    diversity_data_list[[country]] <- data
    
    cat("已读取", country, "的数据：", nrow(data), "行\n")
  } else {
    cat("警告：找不到文件", file_path, "\n")
  }
}

# 合并所有数据
combined_diversity_data <- bind_rows(diversity_data_list)

# 保存合并后的数据
output_file <- "/home/Dong/Global_COVID/genomic_analysis/combined_diversity_data.csv"
write_csv(combined_diversity_data, output_file)
cat("\n合并后的数据已保存到：", output_file, "\n")

# 显示前几行数据
cat("\n前10行数据预览：\n")
print(head(combined_diversity_data, 10))

# 转换为宽格式：18行（国家）× 月份列
cat("\n转换数据为宽格式（国家 × 月份）...\n")

# 为了便于识别，创建年月格式的列名
combined_diversity_data <- combined_diversity_data %>%
  mutate(year_month = format(month, "%Y-%m"))

# 转换为宽格式 - 核苷酸多样性
diversity_wide <- combined_diversity_data %>%
  select(country, year_month, nucleotide_diversity) %>%
  pivot_wider(
    names_from = year_month,
    values_from = nucleotide_diversity,
    values_fill = NA
  ) %>%
  arrange(country)

# 转换为宽格式 - 样本量
sample_size_wide <- combined_diversity_data %>%
  select(country, year_month, sample_size) %>%
  pivot_wider(
    names_from = year_month,
    values_from = sample_size,
    values_fill = NA
  ) %>%
  arrange(country)

# 保存宽格式数据
diversity_wide_file <- "/home/Dong/Global_COVID/genomic_analysis/diversity_wide_format.csv"
sample_size_wide_file <- "/home/Dong/Global_COVID/genomic_analysis/sample_size_wide_format.csv"

write_csv(diversity_wide, diversity_wide_file)
write_csv(sample_size_wide, sample_size_wide_file)




# 将月度多样性数据转换为周度数据
library(dplyr)
library(readr)
library(tidyr)
library(lubridate)
library(zoo)

cat("开始将月度数据转换为周度数据...\n")

# 读取现有的长格式数据
combined_data <- read_csv("/home/Dong/Global_COVID/genomic_analysis/combined_diversity_data.csv", show_col_types = FALSE)

# 方法1: 线性插值转换为周度数据
cat("方法1: 使用线性插值转换为周度数据\n")

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
  diversity_ts <- zoo(country_data$nucleotide_diversity, country_data$month)
  diversity_weekly <- na.approx(diversity_ts, xout = all_weeks, na.rm = FALSE, rule = 2)
  
  # 线性插值样本量（四舍五入为整数）
  sample_ts <- zoo(country_data$sample_size, country_data$month)
  sample_weekly <- round(na.approx(sample_ts, xout = all_weeks, na.rm = FALSE, rule = 2))
  
  # 合并数据
  weekly_frame$nucleotide_diversity <- as.numeric(diversity_weekly)
  weekly_frame$sample_size <- as.integer(sample_weekly)
  weekly_frame$year_week <- format(weekly_frame$week, "%Y-W%U")
  
  weekly_data_list[[country_name]] <- weekly_frame
}

# 合并所有国家的周度数据
weekly_combined <- bind_rows(weekly_data_list)

cat("生成的周度数据概况:\n")
cat("总行数:", nrow(weekly_combined), "\n")
cat("国家数:", n_distinct(weekly_combined$country), "\n")
cat("时间范围:", min(weekly_combined$week), "到", max(weekly_combined$week), "\n")
cat("总周数:", n_distinct(weekly_combined$week), "\n")

# 保存长格式的周度数据
write_csv(weekly_combined, "/home/Dong/Global_COVID/genomic_analysis/weekly_diversity_long.csv")
cat("长格式周度数据已保存到: weekly_diversity_long.csv\n")

# 转换为宽格式 - 核苷酸多样性
cat("\n转换为宽格式...\n")
diversity_weekly_wide <- weekly_combined %>%
  select(country, year_week, nucleotide_diversity) %>%
  pivot_wider(
    names_from = year_week,
    values_from = nucleotide_diversity,
    values_fill = NA
  ) %>%
  arrange(country)

# 转换为宽格式 - 样本量
sample_weekly_wide <- weekly_combined %>%
  select(country, year_week, sample_size) %>%
  pivot_wider(
    names_from = year_week,
    values_from = sample_size,
    values_fill = NA
  ) %>%
  arrange(country)

# 保存宽格式数据
write_csv(diversity_weekly_wide, "/home/Dong/Global_COVID/genomic_analysis/diversity_weekly_wide.csv")
write_csv(sample_weekly_wide, "/home/Dong/Global_COVID/genomic_analysis/sample_size_weekly_wide.csv")

cat("宽格式周度数据已保存:\n")
cat("- 核苷酸多样性: diversity_weekly_wide.csv\n")
cat("- 样本量: sample_size_weekly_wide.csv\n")

# 显示数据概况
cat("\n宽格式数据维度:\n")
cat("核苷酸多样性:", nrow(diversity_weekly_wide), "行 ×", ncol(diversity_weekly_wide)-1, "周\n")
cat("样本量:", nrow(sample_weekly_wide), "行 ×", ncol(sample_weekly_wide)-1, "周\n")

# 方法2: 每月数据重复填充到该月的所有周
cat("\n方法2: 月度数据重复填充到各周\n")

weekly_repeat_list <- list()

for (country_name in unique(combined_data$country)) {
  country_data <- combined_data %>%
    filter(country == country_name) %>%
    arrange(month)
  
  if (nrow(country_data) == 0) next
  
  weekly_repeat_data <- data.frame()
  
  for (i in 1:nrow(country_data)) {
    month_date <- country_data$month[i]
    # 获取该月的所有周一日期
    month_start <- floor_date(month_date, "month")
    month_end <- ceiling_date(month_date, "month") - days(1)
    month_weeks <- seq(from = floor_date(month_start, "week"), 
                       to = ceiling_date(month_end, "week") - days(1), 
                       by = "week")
    
    # 筛选真正属于该月的周
    month_weeks <- month_weeks[month(month_weeks) == month(month_date) | 
                              (month_weeks >= month_start & month_weeks <= month_end)]
    
    # 为该月的每一周重复相同的数据
    week_data <- data.frame(
      week = month_weeks,
      country = country_name,
      nucleotide_diversity = country_data$nucleotide_diversity[i],
      sample_size = country_data$sample_size[i],
      year_week = format(month_weeks, "%Y-W%U"),
      original_month = month_date
    )
    
    weekly_repeat_data <- rbind(weekly_repeat_data, week_data)
  }
  
  weekly_repeat_list[[country_name]] <- weekly_repeat_data
}

weekly_repeat_combined <- bind_rows(weekly_repeat_list)

# 保存重复填充方法的结果
write_csv(weekly_repeat_combined, "/home/Dong/Global_COVID/genomic_analysis/weekly_diversity_repeat.csv")

cat("重复填充方法结果已保存到: weekly_diversity_repeat.csv\n")
cat("该方法生成", nrow(weekly_repeat_combined), "行数据\n")

# 比较两种方法
cat("\n两种方法对比:\n")
cat("方法1 (线性插值):", nrow(weekly_combined), "行\n")
cat("方法2 (重复填充):", nrow(weekly_repeat_combined), "行\n")

cat("\n转换完成！\n")
cat("推荐使用方法1的线性插值结果，因为它提供了更平滑的时间序列\n")# 将月度多样性数据转换为周度数据
library(dplyr)
library(readr)
library(tidyr)
library(lubridate)
library(zoo)

cat("开始将月度数据转换为周度数据...\n")

# 读取现有的长格式数据
combined_data <- read_csv("/home/Dong/Global_COVID/genomic_analysis/combined_diversity_data.csv", show_col_types = FALSE)

# 方法1: 线性插值转换为周度数据
cat("方法1: 使用线性插值转换为周度数据\n")

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
  diversity_ts <- zoo(country_data$nucleotide_diversity, country_data$month)
  diversity_weekly <- na.approx(diversity_ts, xout = all_weeks, na.rm = FALSE, rule = 2)
  
  # 线性插值样本量（四舍五入为整数）
  sample_ts <- zoo(country_data$sample_size, country_data$month)
  sample_weekly <- round(na.approx(sample_ts, xout = all_weeks, na.rm = FALSE, rule = 2))
  
  # 合并数据
  weekly_frame$nucleotide_diversity <- as.numeric(diversity_weekly)
  weekly_frame$sample_size <- as.integer(sample_weekly)
  weekly_frame$year_week <- format(weekly_frame$week, "%Y-W%U")
  
  weekly_data_list[[country_name]] <- weekly_frame
}

# 合并所有国家的周度数据
weekly_combined <- bind_rows(weekly_data_list)

cat("生成的周度数据概况:\n")
cat("总行数:", nrow(weekly_combined), "\n")
cat("国家数:", n_distinct(weekly_combined$country), "\n")
cat("时间范围:", min(weekly_combined$week), "到", max(weekly_combined$week), "\n")
cat("总周数:", n_distinct(weekly_combined$week), "\n")

# 保存长格式的周度数据
write_csv(weekly_combined, "/home/Dong/Global_COVID/genomic_analysis/weekly_diversity_long.csv")
cat("长格式周度数据已保存到: weekly_diversity_long.csv\n")
divmax=max(weekly_combined$nucleotide_diversity)
weekly_combined$nucleotide_diversity<-weekly_combined$nucleotide_diversity/divmax
# 转换为宽格式 - 核苷酸多样性
cat("\n转换为宽格式...\n")
diversity_weekly_wide <- weekly_combined %>%
  select(country, year_week, nucleotide_diversity) %>%
  pivot_wider(
    names_from = year_week,
    values_from = nucleotide_diversity,
    values_fill = 0
  ) %>%
  arrange(country)

# 转换为宽格式 - 样本量
sample_weekly_wide <- weekly_combined %>%
  select(country, year_week, sample_size) %>%
  pivot_wider(
    names_from = year_week,
    values_from = sample_size,
    values_fill = 0
  ) %>%
  arrange(country)

# 保存宽格式数据
write_csv(diversity_weekly_wide, "/home/Dong/Global_COVID/genomic_analysis/diversity_weekly_wide.csv")
write_csv(sample_weekly_wide, "/home/Dong/Global_COVID/genomic_analysis/sample_size_weekly_wide.csv")


