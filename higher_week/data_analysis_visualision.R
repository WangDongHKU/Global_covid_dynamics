# 设置编译器优化级别为最高 (-O3)，并启用本地 CPU 架构优化 (-march=native)
Sys.setenv("CXXFLAGS" = "-O3 -march=native -ffast-math")

# 可选：针对多核 CPU 启用并行编译（如 Make 的 -j 参数）
Sys.setenv("MAKEFLAGS" = "-j4")  # 根据 CPU 核心数调整（例如 4 核）

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(tidyverse)
library(ggpubr)
library(R.matlab)
#library(xlsx)
knitr::opts_chunk$set(cache = T, echo = T, message = F, warning = F,include = T)
theme_set(theme_bw())
# Colourblind friendly colours
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
scale_colour_discrete <- function(...)
  scale_colour_manual(..., values = cbPalette)
scale_fill_discrete <- function(...)
  scale_fill_manual(..., values = cbPalette)

setwd("/home/Dong/Global_COVID/higher_week")
  
load("/home/Dong/Global_COVID/higher_week/higher_model_100.Rdata")

nuts_fit_4_summary <- summary(stan_fit, pars = c("S0","E0","I0","p","v","sp","npi","phi","thetap","rate","diversity_rate","drift_rate"))$summary

print(nuts_fit_4_summary,scientific=FALSE,digits=3)
  
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(tidyverse)
library(ggpubr)
library(R.matlab)
#library(xlsx)
knitr::opts_chunk$set(cache = T, echo = T, message = F, warning = F,include = T)
theme_set(theme_bw())
# Colourblind friendly colours
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
scale_colour_discrete <- function(...)
  scale_colour_manual(..., values = cbPalette)
scale_fill_discrete <- function(...)
  scale_fill_manual(..., values = cbPalette)

custom_labeller = custom_labeller <- function(variable, value) {
  return(c("CN","JP","KR","SG","IN","US","CA","MX","AU","NZ","BR","AR","GB","DE","FR","IT","RU","ZA","SA"))}

library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
# 1. 整理Stan预测病例数据
pred_cases <- tail(multistrain_fit$pred_cases[,seq(1:18),seq(1:244)], 200)
# Reshape the data
data_frame <- reshape2::melt(pred_cases)
colnames(data_frame) <- c("Iteration", "Variable", "Length", "Value")
# Calculate credible intervals
credible_intervals <- data_frame %>%
  group_by(Variable, Length) %>%
  summarize(
    lower = quantile(Value, 0.025),
    upper = quantile(Value, 0.975),
    lower1 = quantile(Value, 0.25),
    upper1 = quantile(Value, 0.75),
    mean = mean(Value)#quantile(Value, 0.5)
  ) %>%
  mutate(Variable = factor(Variable))


# Prepare ILI data
#data_ILI <- sqrt(1+data_lst$cases) %>% as.matrix()%>% t()%>% tibble::as_data_frame() %>% mutate(week = row_number())

data_ILI <-data_lst$cases
data_ILI <- t(data_ILI) %>% as.data.frame()
data_ILI <- data_ILI[-1,]
colnames(data_ILI) <- seq_len(ncol(data_ILI))
data_ILI$id <- seq_len(nrow(data_ILI))


#data_ILI$Variable <- seq_len(col(data_ILI))

data_ILI_long <- reshape2::melt(data_ILI, id.vars = "id")
colnames(data_ILI_long) <- c("Length", "Variable", "Value")
# 转换为数值型
data_ILI_long$Value <- as.numeric(data_ILI_long$Value)

# 再次检查类型
class(data_ILI_long$Value)  # 应返回 "numeric"

# Combine credible intervals and ILI data
credible_intervals <- credible_intervals %>%
  left_join(data_ILI_long, by = c("Variable", "Length"))

# 5. 绘图
plot1 <- ggplot() +
  geom_point(
    data = data_ILI_long, 
    aes(x = Length, y = Value, color = "Data"), 
    size = 0.3, alpha = 0.8
  ) + 
  geom_ribbon(
    data = credible_intervals, 
    aes(x = Length, ymin = lower, ymax = upper, fill = "95% CI"), 
    alpha = 0.15
  ) +
  geom_line(
    data = credible_intervals, 
    aes(x = Length, y = mean, color = "Fitting")
  ) +
  theme_classic() +
  scale_fill_manual(name = "Interval", 
                    values = c("95% CI" = "red")) +
  scale_color_manual(name = "Legend", 
                     values = c("Data" = "black", "Fitting" = "red")) +
  facet_wrap(~ Variable, ncol = 5, scales = "free", labeller = custom_labeller) #, strip.position = "top"

# ggsave("Fig3_social_impact_fitting2_hig15.pdf", plot1, width = 19, height = 11)
# pdf("Fig3_social_impact_fitting2_copy.pdf", width = 19, height = 11, useDingbats = FALSE)
plot1

# dev.off()

ggsave("cases_fit_by_country_hig.pdf", plot1, width = 16, height = 10, dpi = 300)
