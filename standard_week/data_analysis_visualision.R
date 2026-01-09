# 设置编译器优化级别为最高 (-O3)，并启用本地 CPU 架构优化 (-march=native)
Sys.setenv("CXXFLAGS" = "-O3 -march=native -ffast-math")
# 可选：针对多核 CPU 启用并行编译（如 Make 的 -j 参数）
Sys.setenv("MAKEFLAGS" = "-j4")  # 根据 CPU 核心数调整（例如 4 核）
options(bitmapType = "cairo")

X11(type = "cairo")
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
#library(tidyverse)
library(ggpubr)
# library(R.matlab)
#library(xlsx)
knitr::opts_chunk$set(cache = T, echo = T, message = F, warning = F,include = T)
theme_set(theme_bw())
# Colourblind friendly colours
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
scale_colour_discrete <- function(...)
  scale_colour_manual(..., values = cbPalette)
scale_fill_discrete <- function(...)
  scale_fill_manual(..., values = cbPalette)

setwd ("home/dongw21/Global_COVID/standard_week")
  
load("/home/dongw21/Global_COVID/standard_week/standard_model_mainland1000.Rdata")
#
nuts_fit_4_summary <- summary(stan_fit, pars = c("S0","E0","I0","npi","phi","thetap"))$summary

print(nuts_fit_4_summary,scientific=FALSE,digits=3)


library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
# library(ggpubr)
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


# After melting pred_cases
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

# After melting data_ILI
data_ILI_long <- reshape2::melt(data_ILI, id.vars = "id")
colnames(data_ILI_long) <- c("Length", "Variable", "Value")
data_ILI_long$Value <- as.numeric(data_ILI_long$Value)
data_ILI_long$Variable <- factor(data_ILI_long$Variable, levels = 1:18, labels = names(country_map))  # Match above



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

ggsave("cases_fit_by_countrystd_2000.pdf", plot1, width = 16, height = 10, dpi = 300)




