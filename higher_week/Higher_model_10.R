library(rstan)
options(mc.cores = 4)  # 明确设置为4核
rstan_options(auto_write = TRUE)
library(dplyr)

.libPaths(c("/home/dongw21/R_libs/4.4.3", .libPaths()))
if (!require(stringr, quietly = TRUE)) {
  # 尝试多个 CRAN 镜像
  cran_mirrors <- c(
    "https://mirrors.tuna.tsinghua.edu.cn/CRAN/",  # 清华大学镜像
    "https://mirror.lzu.edu.cn/CRAN/",             # 兰州大学镜像
    "https://cloud.r-project.org/"                 # 官方镜像
  )
  installed <- FALSE
  for (mirror in cran_mirrors) {
    tryCatch({
      install.packages("stringr", 
                       lib = "/home/dongw21/R_libs/4.4.3",
                       repos = mirror,
                       dependencies = TRUE)
      library(stringr, lib.loc = "/home/dongw21/R_libs/4.4.3")
      installed <- TRUE
      break
    }, error = function(e) {
      cat("Failed to install from", mirror, "\n")
    })
  }
  if (!installed) {
    stop("无法安装 stringr 包，请检查网络连接或手动安装")
  }
}
# 设置R库路径
# 加载RcppParallel以提供TBB支持
if (!require(RcppParallel, quietly = TRUE)) {
  install.packages("RcppParallel", 
                   lib = "/home/dongw21/R_libs/4.4.3",
                   repos = "https://cloud.r-project.org/")
  library(RcppParallel, lib.loc = "/home/dongw21/R_libs/4.4.3")
}
# 设置TBB线程数
RcppParallel::setThreadOptions(numThreads = 4)

library(ggplot2)
library(R.matlab)
#library(xlsx)
# Colourblind friendly colours
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
scale_colour_discrete <- function(...)
  scale_colour_manual(..., values = cbPalette)
scale_fill_discrete <- function(...)
  scale_fill_manual(..., values = cbPalette)

#load("/scr/u/dongw21/Global_COVID/standard_week/standard_model_mainland4000.Rdata")
#a_raw <- matrix(apply(multistrain_fit$a_raw, c(2,3),  function(x) quantile(x, probs = 0.5, na.rm = TRUE)), ncol = 22, nrow  = 18)
#a_raw

load("/scr/u/dongw21/Global_COVID/higher_week/a_raw.RData")
a_raw

#load("/scr/u/dongw21/Global_COVID/higher_week/higher_model_40.Rdata")
#posterior_samples <- rstan::extract(fit)

m <- rstan::stan_model('/scr/u/dongw21//Global_COVID/higher_week/global_covid_hig_final1.stan')
m
options(max.prints=999)
subtypes <- c("CN","JP","KR","SG","IN","US","CA","MX","AU","NZ","BR","AR","GB","DE","FR","IT","RU","ZA")
row_names <-subtypes
#load("weekly_cases.RData")
#load("weekly_Stringency.RData")
#load("/scr/u/dongw21/Global_COVID/mobilityweek_cov.RData")
load("/scr/u/dongw21/Global_COVID/mobilityday_cov.RData")

#mobility <- aperm(weekly_mobility, perm = c(3, 1, 2))
#mobility <- mobility[,1:18, 1:18]
mobility <- aperm(daily_mobility, perm = c(3, 1, 2))

countries <- c("China", "Japan", "South Korea", "Singapore", "India", 
               "USA", "Canada", "Mexico", "Australia", "New Zealand", 
               "Brazil", "Argentina", "United Kingdom", "Germany", 
               "France", "Italy", "Russia", "South Africa")
drift <- read.csv("/scr/u/dongw21/Global_COVID/genomic_analysis/drift_weekly_wide1.csv",skip = 0, nrows = 18)
drift <- drift[,1:270] 
drift_ordered <- drift[match(countries, drift$country), ]
# 更新diversity变量
drift <- drift_ordered[,4:247] / max(drift_ordered[,4:247])

cases <- read.csv("/scr/u/dongw21/Global_COVID/weekly_cases.csv",skip = 0, nrows = 18, row.names = row_names)
cases <- cases[,1:244]
NPI <- read.csv("/scr/u/dongw21/Global_COVID/Stringency_daily.csv",skip = 0, nrows = 18, row.names = row_names)

population <- c(1408000000, 124400000, 51800000, 5910000, 1415000000, 
                333000000, 36990000, 128000000, 26640000, 5330000, 
                214000000, 46000000, 67800000, 83100000, 67300000, 
                59000000, 143000000, 60000000)

cases_numeric <- cases[, -1]  # 排除第1列非数值列
cases_numeric <- as.matrix(cases_numeric)  # 强制转换为数值矩阵
cases_per_million <- sweep(cases_numeric, MARGIN = 1, STATS = population, FUN = "/") * 1e6


K=18;error <- 10^(-18)

X <-  seq(from=1, to=244*7, by=1);  num_knots  <-  20; spline_degree <- 3;

#num_knots  <-  20; 
#interpolate_row <- function(row) {
#  x_original <- 1:17
#  x_new <- seq(1, 17, length.out = 22) spline(x_original, row, xout = x_new, method = "natural")$y}
#a_raw <- t(apply(a_raw, 1, interpolate_row))  

# 数据验证和清理
cases[is.na(cases)] <- 0  # 将NA替换为0
cases[cases < 0] <- 0     # 确保没有负值
NPI[is.na(NPI)] <- 0      # NPI数据清理
drift[is.na(drift)] <- 0  # drift数据清理

# Mobility数据清理
for(i in 1:dim(mobility)[1]) {
  mobility[i,,][is.na(mobility[i,,])] <- 0
}

1:18

data_lst <- list(
  K = K, 
  W = 244,
  N = population,
  cases = cases, 
  NPI1 = NPI,
  Mobility = mobility,
  error = error,
  num_data = 244*7,
  num_knots = num_knots,
  knots = unname(quantile(X,probs=seq(from=0, to=1, length.out = num_knots))),
  spline_degree = spline_degree,
  X = X,
  a_raw1 = a_raw
)

num_chains <- 4
theta_p <- c(0.019512733, 0.010913988, 0.020111787, 0.009285352, 0.014031503, 0.018381520, 0.018327190, 0.010055737, 0.012994192, 0.011080714,  0.010028904, 0.015409260, 0.009179656, 0.017541051, 0.015997826, 0.014274534, 0.014799658, 0.019314975)

init_lst <- purrr::map(1:num_chains, function(i) {
  list(
    thetap = theta_p%>% as.array(),
    S0 =  c(0.980,0.992),
    E0 =  c(0.0005,0.0001),
    I0 =  c(0.0005,0.0003),
    npi = 0.0009142214,
    phi   = runif(K,0.8,1.2),
    a_raw = a_raw, 
    rate=runif(1,0.028,0.030),
   p= c(0.2102387, 0.2123927, 0.2188851, 0.2167698, 0.2179786, 0.2125815, 0.2194682, 0.2140872, 0.2017209, 0.2027180, 0.2155867, 0.2166512, 0.2084634, 0.2178874, 0.2049716, 0.2020677, 0.2031399, 0.2007355), 
   v = runif(3,0.70,0.73),
   sp= c(1.970110, 1.951526, 1.992947, 2.388402, 2.007994, 1.871765, 1.636906, 1.594359, 2.249022, 2.315858, 1.590475, 1.555430, 1.790928, 1.826716, 1.905530, 1.576002, 2.486666, 1.560680)
   )
})

parameters_stoch = c("p","v","sp","pred_cases","phi","a_raw","S0","E0","I0","npi","S","R","I","E","rate","log_lik","thetap" )

#fit1 <- stan(model_code = m, data = data_lst, iter = 10, chains = 1, warmup = 0) "log_lik"

stan_fit <- rstan::sampling(
  m, 
  data = data_lst, 
  pars = parameters_stoch,
  chains = num_chains, 
  iter = 2000,  
  thin = 1, 
  verbose = TRUE,
  show_messages = TRUE,
  seed = sample.int(.Machine$integer.max, 1),
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  init = init_lst
)
#
rstan::check_hmc_diagnostics(stan_fit)

multistrain_fit <- rstan::extract(stan_fit)

save(list = ls(), file = "/scr/u/dongw21/Global_COVID/higher_week/higher_model_2000.Rdata")

nuts_fit_4_summary <- summary(stan_fit, pars = c("S0","E0","I0","p","v","sp","npi","phi","thetap","rate"))$summary

print(nuts_fit_4_summary,scientific=FALSE,digits=3,probs=c(0.025,0.975))


library(rstan)
# 确保RStan使用4个核心（与Slurm配置一致）
options(mc.cores = 4)
rstan_options(auto_write = TRUE)
#library(tidyverse)
library(ggpubr)
#library(R.matlab)
#library(xlsx)
# Colourblind friendly colours
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
scale_colour_discrete <- function(...)
  scale_colour_manual(..., values = cbPalette)
scale_fill_discrete <- function(...)
  scale_fill_manual(..., values = cbPalette)
custom_labeller = custom_labeller <- function(variable, value) {
  return(c("CN","JP","KR","SG","IN","US","CA","MX","AU","NZ","BR","AR","GB","DE","FR","IT","RU","ZA"))}
setwd("/scr/u/dongw21/Global_COVID/higher_week")

# Define the mapping as a named vector (levels -> labels)
country_map <- c(
  "1" = "CN", "2" = "JP", "3" = "KR", "4" = "SG", "5" = "IN",
  "6" = "US", "7" = "CA", "8" = "MX", "9" = "AU", "10" = "NZ",
  "11" = "BR", "12" = "AR", "13" = "GB", "14" = "DE", "15" = "FR",
  "16" = "IT", "17" = "RU", "18" = "ZA"
)

# Wrap it as a labeller
custom_labeller <- as_labeller(country_map)

library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)

# Define the mapping as a named vector (levels -> labels)
country_map <- c(
  "1" = "CN", "2" = "JP", "3" = "KR", "4" = "SG", "5" = "IN",
  "6" = "US", "7" = "CA", "8" = "MX", "9" = "AU", "10" = "NZ",
  "11" = "BR", "12" = "AR", "13" = "GB", "14" = "DE", "15" = "FR",
  "16" = "IT", "17" = "RU", "18" = "ZA"
)

# Wrap it as a labeller
custom_labeller <- as_labeller(country_map)

pred_cases <- tail(multistrain_fit$pred_cases[,seq(1:18),seq(1:244)], 1000)
# Reshape the data
data_frame <- reshape2::melt(pred_cases)

colnames(data_frame) <- c("Iteration", "Variable", "Length", "Value")
# ... (your credible_intervals code)
credible_intervals <- data_frame %>%
  group_by(Variable, Length) %>%
  summarize(
    lower = quantile(Value, 0.025),
    upper = quantile(Value, 0.975),
    lower1 = quantile(Value, 0.25),
    upper1 = quantile(Value, 0.75),
    mean = mean(Value)
  ) %>%
  mutate(Variable = factor(Variable, levels = 1:18, labels = names(country_map)))  # Explicit levels/labels


data_ILI <-data_lst$cases
data_ILI <- t(data_ILI) %>% as.data.frame()
data_ILI <- data_ILI[-1,]
colnames(data_ILI) <- seq_len(ncol(data_ILI))
data_ILI$id <- seq_len(nrow(data_ILI))

# After melting data_ILI
data_ILI_long <- reshape2::melt(data_ILI, id.vars = "id")
colnames(data_ILI_long) <- c("Length", "Variable", "Value")
data_ILI_long$Value <- as.numeric(data_ILI_long$Value)
# 正确写法（Higher_model_10.R 已是这样）
data_ILI_long$Variable <- factor(
  data_ILI_long$Variable,
  levels = 1:18,
  labels = names(country_map)
)


# Combine (as before)
credible_intervals <- credible_intervals %>%
  left_join(data_ILI_long %>% select(Variable, Length, Value), by = c("Variable", "Length"))  # Select only needed cols

# Plot
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
  scale_y_continuous(trans = "sqrt") +  # Use this for sqrt scale (scale_y_sqrt() is deprecated in some versions)
  scale_fill_manual(name = "Interval", 
                    values = c("95% CI" = "red")) +
  scale_color_manual(name = "Legend", 
                     values = c("Data" = "black", "Fitting" = "red")) +
  facet_wrap(~ Variable, ncol = 5, scales = "free", labeller = custom_labeller)  # Now uses the new labeller

plot1  # Or ggsave() as commented
# dev.off()
ggsave("cases_fit_by_country_hig4400.pdf", plot1, width = 16, height = 10, dpi = 300)
