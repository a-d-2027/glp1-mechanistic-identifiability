# Mechanistic identifiability with correlated mediators

Reproducibility repository for:

**Mechanistic identifiability with correlated mediators: a simulation-based study motivated by GLP-1 signaling and cognition**

## Overview

This repository contains the complete R workflow used for the computational analyses in the manuscript. The study uses **synthetic data only**; no participant-level or identifiable human data are required.

The analysis reproduces:

- literature-calibrated synthetic cohort generation;
- parallel SEM estimation of neural and metabolic downstream pathway products;
- the contrast `Delta = psi_H - psi_M`;
- replicate-specific analysis-scale ground truth;
- bias, MAE, RMSE, 95% confidence-interval coverage, type I error, power, and pathway discrimination;
- mediator-correlation stress testing;
- sample-size stress testing;
- classical mediator-measurement-error stress testing;
- mediator-permutation falsification;
- individual-level 2SLS genetic-instrument simulation;
- GWAS-summary IVW and MR-Egger simulation;
- all computational manuscript and supplementary figures;
- all computational result tables.

## Repository structure

```text
glp1-mechanistic-identifiability/
├── README.md
└── analysis/
    └── glp1_mechanistic_identifiability_analysis.R
```

Running the script creates:

```text
outputs/
├── data/
├── figures/
└── session_info.txt
```

## Required R packages

- MASS
- dplyr
- tidyr
- ggplot2
- lavaan
- AER
- readr
- tibble
- patchwork

The script checks for missing packages and stops with installation instructions. It does not automatically install software.

## Running the complete analysis

From the repository root:

```bash
Rscript analysis/glp1_mechanistic_identifiability_analysis.R
```

Default settings reproduce the manuscript:

```r
PRIMARY_MC_REPS <- 500
SECONDARY_MC_REPS <- 50
```

For a short code check, temporarily reduce both values to 5.

## Reproducibility

Fixed seeds are defined in `run_manuscript_analysis()`. Each Monte Carlo replicate receives a deterministic derived seed.

A completed run writes:

```text
outputs/session_info.txt
```

which records the R version, platform, and package versions.

## Analysis-scale ground truth

The SEM is fitted to standardized signaling, mediators, and outcome. Consequently, raw data-generating coefficients are transformed to the standardized analysis scale within each replicate:

```text
a_true = beta_g / SD(H_latent)
b_true = gamma_g / SD(M_latent)
c_true = delta_h / SD(Y_raw)
d_true = delta_m / SD(Y_raw)

psi_H,true = a_true * c_true
psi_M,true = b_true * d_true
Delta_true = psi_H,true - psi_M,true
```

Bias, RMSE, and confidence-interval coverage are evaluated against these analysis-scale targets.


## Complete stochastic specification

The R code and `table_s2_simulation_parameters.csv` provide the exact data-generating distributions used in the simulation. In particular:

- `epsilon_G ~ N(0,1)` and `epsilon_Y ~ N(0,1)`;
- `(epsilon_H, epsilon_M)` is bivariate normal with mean `(0,0)` and covariance matrix `[[1, rho_HM], [rho_HM, 1]]`, so both mediator residual variances equal 1;
- receptor sensitivity is `N(1.3, 0.1^2)` for treated observations and `N(1.0, 0.1^2)` for controls, followed by within-cohort standardization;
- age is drawn from `N(71.8, 7.1^2)`, clipped to `[55, 90]`, and standardized;
- BMI is drawn from `N(26.5, 5.2^2)`, clipped to `[16, 55]`, and standardized;
- the HOMA-like residual is `N(0, 0.5^2)`.

The age and BMI implementation uses clipping (`pmax`/`pmin`), not sampling from a mathematically truncated normal distribution.

## Pathway discrimination

Correct discrimination is defined only for scenarios with a prespecified dominant pathway:

- neural-dominant: `Delta_hat > 0`;
- metabolic-dominant: `Delta_hat < 0`.

Classification accuracy is not assigned to the shared or null scenarios. The heterogeneous scenario is summarized using the operating characteristics of `Delta`.

## Manuscript figure mapping

The current R script generates the full computational figure set.

| Manuscript item | R function | Output file |
|---|---|---|
| Figure 2. Estimated pathway contrast across mediator correlation | `plot_figure_2_pathway_contrast()` | `outputs/figures/figure_2_pathway_contrast.png` |
| Figure 3. GWAS-summary genetic-instrument simulation | `plot_figure_3_gwas_summary()` | `outputs/figures/figure_3_gwas_summary.png` |
| Figure 4. Sample-size stress test | `plot_figure_4_sample_size_stress_test()` | `outputs/figures/figure_4_sample_size_stress_test.png` |
| Figure 5. Mediator-measurement-error stress test | `plot_figure_5_measurement_error_stress_test()` | `outputs/figures/figure_5_measurement_error_stress_test.png` |
| Supplementary Figure S1. SEM fit diagnostics | `plot_figure_s1_sem_fit_diagnostics()` | `outputs/figures/figure_s1_sem_fit_diagnostics.png` |
| Supplementary Figure S2. Individual-level 2SLS | `plot_figure_s2_individual_2sls()` | `outputs/figures/figure_s2_individual_2sls.png` |
| Supplementary Figure S3. SEM/2SLS concordance | `plot_figure_s3_sem_2sls_concordance()` | `outputs/figures/figure_s3_sem_2sls_concordance.png` |

### Current Figure 2

`plot_figure_2_pathway_contrast()` reproduces the **current manuscript layout**: one plotting panel with the five data-generating scenarios on the x-axis and three side-by-side boxplots within each scenario for `rho_HM = 0.2, 0.6, 0.8`.

The plotted quantity is:

```text
Delta_hat = psi_H_hat - psi_M_hat
```

Positive values favor the neural pathway and negative values favor the metabolic pathway.


### Current Figure 3

`plot_figure_3_gwas_summary()` generates the **current two-panel Figure 3 directly in R** from the synthetic GWAS-summary results returned by `run_gwas_summary_analysis()`.

The upper panel, titled **Genetic Validation**, plots the simulated SNP-outcome associations (`by`) against SNP-GLP-1 signaling associations (`bx`) for cognition, education, and Alzheimer disease, with outcome-specific regression lines through the origin.

The lower panel, titled **IVW MR Estimates**, plots the IVW point estimate and 95% confidence interval for each outcome.

The final combined figure is written to:

```text
outputs/figures/figure_3_gwas_summary.png
```

No pre-existing image file is read or copied by this function.


## Manuscript table mapping

| Manuscript item | Output file |
|---|---|
| Main Table 1 | `outputs/data/table_1_delta_operating_characteristics_summary.csv` |
| Main Table 2 | `outputs/data/table_2_individual_2sls_summary.csv` |
| Supplementary Table S2 | `outputs/data/table_s2_simulation_parameters.csv` |
| Supplementary Table S3 | `outputs/data/table_s3_primary_delta_operating_characteristics.csv` |
| Supplementary Table S4 | `outputs/data/table_s4_permutation_falsification.csv` |
| Supplementary Table S5 | `outputs/data/table_s5_pathway_recovery_analysis_scale.csv` |
| Supplementary Table S6 | `outputs/data/table_s6_sample_size_stress_test.csv` |
| Supplementary Table S7 | `outputs/data/table_s7_measurement_error_stress_test.csv` |
| Supplementary Table S8 | `outputs/data/table_s8_individual_2sls.csv` |

**Supplementary Table S1** is a literature-context table and is not a computational simulation output; it is therefore not generated by the analysis script.

## Additional exported files

The workflow also writes:

- `primary_mc_runs.csv`;
- `operating_characteristics_all.csv`;
- `secondary_mc_runs.csv`;
- `gwas_summary_simulation.csv`;
- `gwas_mr_summary.csv`;
- `discrimination_interpretation.csv`.

These files allow independent reconstruction and checking of the manuscript summaries.

## Data availability

All analytic data are generated by the R script. No external participant-level dataset is required.

## Citation

Please cite the accompanying manuscript when using or adapting this code.
