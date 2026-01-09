# Global_covid_dynamics
# Global COVID-19 Transmission Dynamics: A Comprehensive Study

## Project Overview

This project presents a comprehensive study of global COVID-19 transmission dynamics, integrating epidemiological modeling, genomic analysis, and multi-factor impact assessment. The research covers 18 countries over a time span from early 2020 to early 2025, employing Bayesian statistical modeling methods (Stan) and phylogenetic analysis approaches to investigate COVID-19 transmission mechanisms, viral evolution dynamics, and the effectiveness of various intervention measures.

### Research Objectives

1. **Transmission Dynamics Modeling**: Using the SEIR (Susceptible-Exposed-Infectious-Recovered) model framework combined with B-spline functions to perform time-series modeling of COVID-19 transmission across 18 countries
2. **Multi-Factor Impact Analysis**: Assessing the impact of non-pharmaceutical interventions (NPI), human mobility, and variants of concern (VOC) on transmission
3. **Genomic Evolution Analysis**: Tracking viral evolution trajectories through phylogenetic analysis, calculating antigenic drift and phylogenetic diversity
4. **Model Comparison**: Comparing the fitting performance of standard and higher-order models to evaluate different modeling strategies

## Project Structure

```
Global_COVID/
├── standard_week/          # Standard SEIR model
│   ├── global_covid_std_final.stan    # Stan model definition
│   ├── standard_model_10.R           # Main analysis script
│   ├── data_analysis_visualision.R   # Data analysis and visualization
│   └── README_进度跟踪.md             # Progress tracking documentation
│
├── higher_week/            # Higher-order SEIR model (with additional parameters)
│   ├── global_covid_hig_final.stan    # Higher-order Stan model
│   ├── Higher_model_10_clean.R        # Higher-order model analysis script
│   └── a_raw.RData                    # Pre-trained parameters
│
├── genomic_analysis/        # Genomic analysis module
│   ├── genomic_analysis.R              # Main phylogenetic analysis script
│   ├── README_phylogenetic_analysis.md # Phylogenetic analysis documentation
│   ├── drift_weekly_wide1.csv          # Antigenic drift data
│   └── *_VOC_proportions.csv           # Country-specific VOC proportion data
│
├── semi_mechinistic/        # Semi-mechanistic model
│   └── semi_mechinistic.stan           # Semi-mechanistic Stan model
│
├── weekly_cases.csv         # Weekly case data for 18 countries (2020-2025)
├── Stringency_daily.csv     # Policy stringency index (daily)
├── Stringency_weekly.csv    # Policy stringency index (weekly)
├── global_mobility_OAG.CSV  # Global mobility data
├── mobilityday_cov.RData    # Daily mobility data (R format)
└── mobilityweek_cov.RData   # Weekly mobility data (R format)
```

## Study Countries

The research covers the following 18 countries:

| Code | Country | Code | Country |
|------|---------|------|---------|
| CN | China | US | United States |
| JP | Japan | CA | Canada |
| KR | South Korea | MX | Mexico |
| SG | Singapore | AU | Australia |
| IN | India | NZ | New Zealand |
| BR | Brazil | GB | United Kingdom |
| AR | Argentina | DE | Germany |
| FR | France | IT | Italy |
| RU | Russia | ZA | South Africa |

## Data Description

### 1. Case Data (`weekly_cases.csv`)
- **Format**: 18 rows (countries) × 244 columns (weeks)
- **Time Range**: December 30, 2019 to January 27, 2025
- **Unit**: Weekly new confirmed cases
- **Data Source**: Official health departments of respective countries

### 2. Policy Stringency Data (`Stringency_daily.csv`, `Stringency_weekly.csv`)
- **Indicator**: Oxford COVID-19 Government Response Tracker (OxCGRT) Stringency Index
- **Range**: 0-100, higher values indicate stricter policies
- **Dimensions**: 18 countries × time series

### 3. Mobility Data (`mobilityday_cov.RData`, `mobilityweek_cov.RData`)
- **Format**: Three-dimensional array [time points × 18 countries × 18 countries]
- **Meaning**: Inter-country human mobility flow matrix
- **Data Source**: OAG (Official Airline Guide) aviation data

### 4. VOC Proportion Data (`genomic_analysis/*_VOC_proportions.csv`)
- **VOC Types**:
  - VOC Original (Original strain)
  - VOC Alpha (Alpha variant)
  - VOC Beta (Beta variant)
  - VOC Gamma (Gamma variant)
  - VOC Delta (Delta variant)
  - VOC Omicron (Omicron variant)
- **Format**: One CSV file per country, containing proportions of each VOC at different time points

### 5. Antigenic Drift Data (`genomic_analysis/drift_weekly_wide1.csv`)
- **Definition**: Average genetic distance between adjacent time points
- **Purpose**: Quantifying the rate of change in viral antigenic properties
- **Format**: 18 countries × time series

## Model Description

### Standard Model (`standard_week/`)

**Model Framework**: SEIR model + B-spline time-varying parameters

**Main Parameters**:
- `thetap[K]`: Initial transmission rate for each country
- `S0[2]`, `E0[2]`, `I0[2]`: Initial proportions of susceptible, exposed, and infectious populations
- `npi`: Non-pharmaceutical intervention effectiveness parameter
- `phi[K]`: Overdispersion parameter for each country
- `a_raw[K, num_basis]`: B-spline coefficients (time-varying transmission rate)
- `rate`: Recovery rate

**Key Features**:
- Uses B-spline functions to model time-varying transmission rates
- Integrates mobility matrices and NPI data
- Accounts for overdispersion using negative binomial distribution

### Higher-Order Model (`higher_week/`)

**Extended Parameters** (compared to standard model):
- `p[K]`: Country-specific vaccination-related parameters
- `v[3]`: Vaccine effectiveness parameters
- `sp[K]`: Country-specific parameters

**Application Scenarios**:
- Requires more fine-grained parameterization
- Considers vaccination impact
- Requires higher model flexibility

### Genomic Analysis (`genomic_analysis/`)

**Main Functions**:

1. **Phylogenetic Tree Construction**
   - Methods: Maximum Likelihood (ML), Neighbor-Joining (NJ), UPGMA
   - Model Selection: Automatic selection of best-fit nucleotide substitution models (JC, F81, K80, HKY, TrN, GTR, etc.)

2. **Phylogenetic Diversity Calculation**
   - Definition: Sum of all branch lengths in the phylogenetic tree
   - Application: Tracking viral evolution speed and differentiation degree

3. **Antigenic Drift Analysis**
   - Calculates average genetic distance between adjacent time points
   - Identifies significant drift events (Z-score > 2.0)

4. **VOC Time Series Analysis**
   - Tracks proportion changes of each VOC across different countries
   - Calculates VOC-specific case numbers

## Usage Instructions

### System Requirements

**R Package Dependencies**:
```r
# Stan-related
library(rstan)
library(reshape2)
library(dplyr)
library(tidyr)

# Visualization
library(ggplot2)

# Genomic analysis
library(ape)
library(msa)
library(Biostrings)
library(ggtree)
library(treeio)
```

**System Requirements**:
- R >= 4.0
- Stan >= 2.26
- MAFFT (for sequence alignment)
- Recommended memory >= 32GB (for large-scale genomic analysis)

### Running the Standard Model

```r
# 1. Set working directory
setwd("/scr/u/dongw21/Global_COVID/standard_week")

# 2. Run main analysis script
source("standard_model_10.R")

# Or submit job using SLURM
sbatch standard_model_10.sh
```

**Model Running Parameters**:
- Number of chains: 4
- Iterations: 4000 (adjustable as needed)
- Adaptation period: Auto-adjusted
- Control parameters: `adapt_delta = 0.95`, `max_treedepth = 15`

### Running the Higher-Order Model

```r
setwd("/scr/u/dongw21/Global_COVID/higher_week")
source("Higher_model_10_clean.R")
```

### Running Genomic Analysis

```r
setwd("/scr/u/dongw21/Global_COVID/genomic_analysis")

# Install dependencies (first run)
source("install_phylo_packages_simple.R")

# Run complete analysis
source("genomic_analysis.R")
```

**Analysis Levels**:
1. Full sample analysis (stratified sampling, maximum 30 sequences per month)
2. Quarterly analysis (grouped by quarter)
3. Representative month detailed analysis (months with sample size ≥ 50)
4. Phylogenetic diversity statistics

### Monitoring Model Running Progress

```bash
# Real-time monitoring (recommended)
cd /scr/u/dongw21/Global_COVID/standard_week
./monitor_progress.sh

# Simple check
./check_progress.sh

# Check SLURM job status
squeue -u $USER
```

## Output Files Description

### Model Fitting Results

**Standard Model Outputs**:
- `standard_model_global4000.Rdata`: Complete model results (R format)
- `cases_fit_by_countrystd5000.pdf`: Case fitting plots by country
- `voc_specific_cases_by_country_line.pdf`: VOC-specific case time series
- `voc_model_fitting_results5000.pdf`: VOC-specific model fitting performance

**Higher-Order Model Outputs**:
- `higher_model_40.Rdata`: Higher-order model results
- `cases_fit_by_country_hig4400.pdf`: Higher-order model fitting plots

### Genomic Analysis Outputs

**Phylogenetic Trees**:
- `global_phylo_basic_tree.pdf`: Full basic phylogenetic tree
- `global_phylo_circular_tree.pdf`: Full circular phylogenetic tree
- `global_phylo_annotated_tree.pdf`: Time-annotated tree
- `quarterly_phylo_YYYY-QX_*.pdf`: Quarterly phylogenetic trees
- `phylo_YYYY-MM_*.pdf`: Monthly phylogenetic trees

**Data Files**:
- `global_phylo_tree.newick`: Phylogenetic tree in Newick format
- `monthly_phylogenetic_diversity.csv`: Monthly phylogenetic diversity
- `monthly_antigenic_drift.csv`: Monthly antigenic drift analysis
- `significant_antigenic_drift_events.csv`: Significant drift events

**VOC Analysis**:
- `VOC_*_timeseries.csv`: Time series data for each VOC
- `VOC_timeseries_summary.csv`: VOC time series summary statistics

## Key Results Interpretation

### 1. Model Fitting Quality

**Evaluation Metrics**:
- R-hat < 1.01: Good chain convergence
- ESS (Effective Sample Size) > 400: Sufficient sampling
- 95% credible intervals cover observed data: Good model fit

**Visual Inspection**:
- Predicted mean line should be close to observed data points
- 95% credible intervals should contain most observed values
- Residuals should be randomly distributed without systematic bias

### 2. Parameter Estimation

**Key Parameter Meanings**:
- `thetap`: Initial transmission rate, reflecting transmission intensity in the early stage of the epidemic in each country
- `npi`: NPI effectiveness, higher values indicate more effective interventions
- `phi`: Overdispersion parameter, reflecting data variability
- `a_raw`: B-spline coefficients, describing temporal patterns of transmission rate changes

### 3. VOC Dynamics

**VOC Proportion Changes**:
- Track the emergence, growth, and decline of each VOC across different countries
- Identify VOC replacement events (e.g., Delta replaced by Omicron)

**VOC-Specific Cases**:
- Decompose total cases by VOC proportions
- Assess the contribution of each VOC to total case numbers

### 4. Phylogenetic Analysis

**Phylogenetic Diversity (PD)**:
- High PD values: Active viral evolution, high genetic diversity
- Low PD values: Relatively stable virus, low genetic diversity

**Antigenic Drift**:
- Z-score > 2.0: Significant drift event
- May correspond to emergence of new VOC or important mutations

## Model Diagnostics

### Stan Model Diagnostics

```r
# Check HMC diagnostics
rstan::check_hmc_diagnostics(stan_fit)

# View parameter summary
summary(stan_fit, pars = c("npi", "phi", "thetap", "rate"))

# Check chain convergence
traceplot(stan_fit, pars = "npi")
```

### Common Issues

1. **Chain Non-Convergence**
   - Increase number of iterations
   - Adjust `adapt_delta` (e.g., 0.99)
   - Check data quality

2. **Low Sampling Efficiency**
   - Check `max_treedepth` settings
   - Consider reparameterization
   - Use more appropriate prior distributions

3. **Insufficient Memory**
   - Reduce number of iterations
   - Use fewer chains
   - Process data in batches

## Performance Optimization Recommendations

### Computational Resources

**Standard Model**:
- Recommended CPU: 8-16 cores
- Memory: 16-32GB
- Runtime: 4-8 hours (4000 iterations)

**Higher-Order Model**:
- Recommended CPU: 16-32 cores
- Memory: 32-64GB
- Runtime: 8-16 hours

**Genomic Analysis**:
- Recommended CPU: 16-32 cores
- Memory: 64-128GB
- Runtime: Several hours to days (depending on number of sequences)

### Optimization Strategies

1. **Parallel Computing**:
   ```r
   options(mc.cores = parallel::detectCores())
   ```

2. **Compiler Optimization**:
   ```r
   Sys.setenv("CXXFLAGS" = "-O3 -march=native -ffast-math")
   ```

3. **Data Preprocessing**:
   - Clean and validate data in advance
   - Save intermediate results in RData format

## Citations and Acknowledgments

### Software and Packages Used

- **Stan**: Bayesian statistical modeling
- **RStan**: R interface to Stan
- **ape**: Analyses of Phylogenetics and Evolution
- **msa**: Multiple Sequence Alignment
- **MAFFT**: Multiple sequence alignment program

### Data Sources

- COVID-19 case data: Official health departments of respective countries
- Policy stringency data: Oxford COVID-19 Government Response Tracker
- Mobility data: OAG (Official Airline Guide)
- Genomic sequences: GISAID

## Contact

For questions or suggestions, please:
- Check README files in subdirectories for detailed documentation
- Review SLURM output files (`.out` and `.err`) for error troubleshooting

## Changelog

### 2024-2025
- Added VOC-specific case analysis
- Enhanced genomic analysis pipeline
- Optimized model visualization
- Added progress tracking system

---

**Note**: This project involves large-scale computation and complex statistical modeling. It is recommended to run on high-performance computing clusters with sufficient computational resources. Please carefully read the README files in each subdirectory before first use.
