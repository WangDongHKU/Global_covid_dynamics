library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(tidyverse)
library(ggplot2)
library(R.matlab)
library(dplyr)
library(ggpubr)
library(reshape2)
library(tidyr)

knitr::opts_chunk$set(cache = T, echo = T, message = F, warning = F, include = T)
theme_set(theme_bw())

# Colourblind friendly colours
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
scale_colour_discrete <- function(...)
  scale_colour_manual(..., values = cbPalette)
scale_fill_discrete <- function(...)
  scale_fill_manual(..., values = cbPalette)

# Load Stan model
m <- rstan::stan_model('/scr/u/dongw21/Global_COVID/higher_week/global_covid_hig_final.stan')
m
options(max.prints=999)

# Define countries and data
subtypes <- c("CN","JP","KR","SG","IN","US","CA","MX","AU","NZ","BR","AR","GB","DE","FR","IT","RU","ZA")
row_names <- subtypes

# Load mobility data
load("/scr/u/dongw21/Global_COVID/mobilityday_cov.RData")
mobility <- aperm(daily_mobility, perm = c(3, 1, 2))

countries <- c("China", "Japan", "South Korea", "Singapore", "India", 
               "USA", "Canada", "Mexico", "Australia", "New Zealand", 
               "Brazil", "Argentina", "United Kingdom", "Germany", 
               "France", "Italy", "Russia", "South Africa")

# Load and process drift data
drift <- read.csv("/scr/u/dongw21/Global_COVID/genomic_analysis/drift_weekly_wide1.csv", skip = 0, nrows = 18)
drift <- drift[,1:270] 
drift_ordered <- drift[match(countries, drift$country), ]
drift <- drift_ordered[,4:247] / max(drift_ordered[,4:247])

# Load cases and NPI data
cases <- read.csv("/scr/u/dongw21/Global_COVID/weekly_cases.csv", skip = 0, nrows = 18, row.names = row_names)
cases <- cases[,1:244]
NPI <- read.csv("/scr/u/dongw21/Global_COVID/Stringency_daily.csv", skip = 0, nrows = 18, row.names = row_names)

# Population data
population <- c(1408000000, 124400000, 51800000, 5910000, 1415000000, 
                333000000, 36990000, 128000000, 26640000, 5330000, 
                214000000, 46000000, 67800000, 83100000, 67300000, 
                59000000, 143000000, 60000000)

# Process cases data
cases_numeric <- cases[, -1]  # Remove first non-numeric column
cases_numeric <- as.matrix(cases_numeric)  # Convert to numeric matrix
cases_per_million <- sweep(cases_numeric, MARGIN = 1, STATS = population, FUN = "/") * 1e6

# Model parameters
K <- 18
error <- 10^(-18)
X <- seq(from=1, to=244*7, by=1)
num_knots <- 20
spline_degree <- 3

# Load a_raw data
load("/scr/u/dongw21/Global_COVID/higher_week/a_raw.RData")
a_raw

# Data validation and cleaning
cases[is.na(cases)] <- 0  # Replace NA with 0
cases[cases < 0] <- 0     # Ensure no negative values
NPI[is.na(NPI)] <- 0      # Clean NPI data
drift[is.na(drift)] <- 0  # Clean drift data

# Clean Mobility data
for(i in 1:dim(mobility)[1]) {
  mobility[i,,][is.na(mobility[i,,])] <- 0
}

# Prepare data list for Stan
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
  knots = unname(quantile(X, probs=seq(from=0, to=1, length.out = num_knots))),
  spline_degree = spline_degree,
  X = X,
  a_raw1 = a_raw
)

# Set up chains and initial values
num_chains <- 4
init_lst <- purrr::map(1:num_chains, function(i) {
  list(
    thetap = runif(K, 0.009, 0.021),
    S0 = c(0.982, 0.992),
    E0 = c(0.0005, 0.0001),
    I0 = c(0.0005, 0.0003),
    npi = runif(1, 0.0009, 0.0011),
    phi = runif(K, 0.8, 1.2),
    a_raw = a_raw, 
    rate = runif(1, 0.02, 0.035),
    p = runif(K, 0.2, 0.22), 
    v = runif(3, 0.70, 0.79),
    sp = runif(K, 1.5, 2.5)
  )
})

# Parameters to track
parameters_stoch <- c("p", "v", "sp", "pred_cases", "phi", "a_raw", "S0", "E0", "I0", "npi", "S", "R", "I", "E", "rate", "log_lik", "thetap")

# Run Stan sampling
stan_fit <- rstan::sampling(
  m, 
  data = data_lst, 
  pars = parameters_stoch,
  chains = num_chains, 
  iter = 400,  
  thin = 1, 
  verbose = TRUE,
  show_messages = TRUE,
  seed = sample.int(.Machine$integer.max, 1),
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  init = init_lst
)

# Check diagnostics
rstan::check_hmc_diagnostics(stan_fit)

# Extract results
multistrain_fit <- rstan::extract(stan_fit)

# Save results
save(list = ls(), file = "/scr/u/dongw21/Global_COVID/higher_week/higher_model_40.Rdata")

# Print summary
nuts_fit_4_summary <- summary(stan_fit, pars = c("S0", "E0", "I0", "p", "v", "sp", "npi", "phi", "thetap", "rate"))$summary
print(nuts_fit_4_summary, scientific=FALSE, digits=3, probs=c(0.025, 0.975))

# Set working directory for plotting
setwd("/scr/u/dongw21/Global_COVID/higher_week")

# Define country mapping for plotting
country_map <- c(
  "1" = "CN", "2" = "JP", "3" = "KR", "4" = "SG", "5" = "IN",
  "6" = "US", "7" = "CA", "8" = "MX", "9" = "AU", "10" = "NZ",
  "11" = "BR", "12" = "AR", "13" = "GB", "14" = "DE", "15" = "FR",
  "16" = "IT", "17" = "RU", "18" = "ZA"
)

# Create custom labeller
custom_labeller <- as_labeller(country_map)

# Prepare prediction data for plotting
pred_cases <- tail(multistrain_fit$pred_cases[, seq(1:18), seq(1:244)], 1000)
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
    mean = mean(Value),
    .groups = 'drop'
  ) %>%
  mutate(Variable = factor(Variable, levels = 1:18, labels = names(country_map)))

# Prepare observed data for plotting
data_ILI <- data_lst$cases
data_ILI <- t(data_ILI) %>% as.data.frame()
data_ILI <- data_ILI[-1,]
colnames(data_ILI) <- seq_len(ncol(data_ILI))
data_ILI$id <- seq_len(nrow(data_ILI))

# Melt observed data
data_ILI_long <- reshape2::melt(data_ILI, id.vars = "id")
colnames(data_ILI_long) <- c("Length", "Variable", "Value")
data_ILI_long$Value <- as.numeric(data_ILI_long$Value)
data_ILI_long$Variable <- factor(data_ILI_long$Variable, levels = 1:18, labels = names(country_map))

# Combine data for plotting
credible_intervals <- credible_intervals %>%
  left_join(data_ILI_long %>% select(Variable, Length, Value), by = c("Variable", "Length"))

# Create plot
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
  scale_y_continuous(trans = "sqrt") +
  scale_fill_manual(name = "Interval", 
                    values = c("95% CI" = "red")) +
  scale_color_manual(name = "Legend", 
                     values = c("Data" = "black", "Fitting" = "red")) +
  facet_wrap(~ Variable, ncol = 5, scales = "free", labeller = custom_labeller)

# Display and save plot
plot1
ggsave("cases_fit_by_country_hig4400.pdf", plot1, width = 16, height = 10, dpi = 300)
