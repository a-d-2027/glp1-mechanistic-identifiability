# ==============================================================================
# Mechanistic identifiability with correlated mediators
# Reproducibility analysis for the GLP-1 signaling and cognition simulation study
# ==============================================================================
#
# Purpose
# -------
# This script reproduces the simulation analyses reported in the manuscript:
#
#   "Mechanistic identifiability with correlated mediators:
#    a simulation-based study motivated by GLP-1 signaling and cognition"
#
# The analysis uses synthetic data only. No participant-level or identifiable
# human data are required.
#
# Main analyses
# -------------
# 1. Generate literature-calibrated synthetic cohorts with two correlated
#    candidate mediators.
# 2. Fit the prespecified parallel structural equation model (SEM).
# 3. Evaluate the pathway contrast
#         Delta = psi_H - psi_M
#    using replicate-specific truth on the standardized analysis scale.
# 4. Estimate bias, MAE, RMSE, 95% confidence-interval coverage, type I error,
#    power, and probability of correct discrimination where a dominant pathway
#    is prespecified.
# 5. Stress-test sample size and classical mediator measurement error.
# 6. Run mediator-permutation falsification.
# 7. Run secondary individual-level 2SLS and GWAS-summary genetic-instrument
#    simulations.
# 8. Export manuscript-ready tables, figures, replicate-level results, and
#    session information.
#
# Reproducibility
# ---------------
# Primary operating-characteristic analyses use 500 Monte Carlo replicates per
# design cell. Secondary 2SLS/permutation summaries use 50 replicates per
# scenario-by-correlation cell, matching the manuscript.
#
# The primary random-number seed is set explicitly in run_manuscript_analysis().
# Each Monte Carlo replicate then receives a deterministic derived seed.
#
# Suggested repository layout
# ---------------------------
#   analysis/
#     glp1_mechanistic_identifiability_analysis.R
#   outputs/
#     data/
#     figures/
#     session_info.txt
#
# Run
# ---
# From the repository root:
#
#   Rscript analysis/glp1_mechanistic_identifiability_analysis.R
#
# For a short code check, change PRIMARY_MC_REPS and SECONDARY_MC_REPS in the
# configuration block near the end of this file.
#
# ==============================================================================


# ==============================================================================
# 1. Packages and global configuration
# ==============================================================================

required_packages <- c(
  "MASS",
  "dplyr",
  "tidyr",
  "ggplot2",
  "lavaan",
  "AER",
  "readr",
  "tibble",
  "patchwork"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required R packages: ",
      paste(missing_packages, collapse = ", "),
      "\nInstall them before running this script, for example:\n",
      "install.packages(c(",
      paste(sprintf('"%s"', missing_packages), collapse = ", "),
      "))"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lavaan)
  library(readr)
  library(tibble)
  library(patchwork)
})

options(stringsAsFactors = FALSE)

SCENARIOS <- c(
  "neural_dominant",
  "metabolic_dominant",
  "shared",
  "heterogeneous",
  "null"
)

PRIMARY_RHOS <- c(0.2, 0.6, 0.8)
SAMPLE_SIZES <- c(500, 1000, 2500, 5000)
MEDIATOR_RELIABILITIES <- c(1.00, 0.80, 0.60)
STRESS_TEST_RHO <- 0.60


# ==============================================================================
# 2. General utilities
# ==============================================================================

#' Create a directory if it does not already exist.
ensure_directory <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}


#' Standardize a numeric vector to mean 0 and standard deviation 1.
standardize <- function(x) {
  as.numeric(scale(x))
}


#' Format a numeric vector as "mean (SD)".
format_mean_sd <- function(x, digits = 3) {
  sprintf(
    paste0("%.", digits, "f (%.", digits, "f)"),
    mean(x, na.rm = TRUE),
    stats::sd(x, na.rm = TRUE)
  )
}


#' Convert internal scenario names to manuscript labels.
scenario_label <- function(x) {
  dplyr::recode(
    x,
    neural_dominant    = "Neural dominant",
    metabolic_dominant = "Metabolic dominant",
    shared              = "Shared",
    heterogeneous       = "Heterogeneous",
    null                = "Null",
    .default = x
  )
}


#' Publication-style plotting theme with an explicit white background.
publication_theme <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.key = ggplot2::element_rect(fill = "white", color = NA)
    )
}


# ==============================================================================
# 3. Simulation parameters
# ==============================================================================

#' Literature-based marginal anchors used to construct the synthetic cohort.
#'
#' These values define plausible marginal scales and GWAS-scale sample sizes.
#' They are not participant-level data.
publication_anchors <- function() {
  list(
    age_mean = 71.8,
    age_sd = 7.1,
    age_min = 55,
    age_max = 90,

    bmi_mean = 26.5,
    bmi_sd = 5.2,
    bmi_min = 16,
    bmi_max = 55,

    hrv_mean = 42,
    hrv_sd = 15,

    cognition_mean = 50,
    cognition_sd = 10,

    homa_mean = 2.6,
    homa_sd = 1.8,

    gwas_sample_sizes = c(
      Cognition = 269867,
      Education = 1100000,
      Alzheimers = 94437
    )
  )
}


#' Return structural parameters for one data-generating scenario.
scenario_parameters <- function(
    scenario = SCENARIOS,
    rho_hm = 0.60
) {
  scenario <- match.arg(scenario)

  alpha_t <- 0.90
  beta_g <- 0.60
  gamma_g <- 0.60

  delta_h <- 0.45
  delta_m <- 0.15
  delta_hs <- 0.20

  if (scenario == "metabolic_dominant") {
    delta_h <- 0.08
    delta_m <- 0.35
    delta_hs <- 0
  } else if (scenario == "shared") {
    delta_h <- 0.25
    delta_m <- 0.25
    delta_hs <- 0
  } else if (scenario == "heterogeneous") {
    delta_h <- 0.22
    delta_m <- 0.18
    delta_hs <- 0.22
  } else if (scenario == "null") {
    delta_h <- 0
    delta_m <- 0
    delta_hs <- 0
  }

  list(
    alpha_t = alpha_t,
    beta_g = beta_g,
    gamma_g = gamma_g,
    delta_h = delta_h,
    delta_m = delta_m,
    rho_hm = rho_hm,
    delta_hs = delta_hs,
    gwas_true_betas = c(
      Cognition = 0.18,
      Education = 0.12,
      Alzheimers = -0.22
    )
  )
}


# ==============================================================================
# 4. Synthetic cohort generation
# ==============================================================================

#' Add classical measurement error to a standardized mediator.
#'
#' For X_obs = X_true + e and Var(X_true) approximately equal to 1,
#'
#'   reliability = Var(X_true) / {Var(X_true) + Var(e)}.
#'
#' The observed mediator is re-standardized before SEM fitting, matching the
#' analysis convention used in the manuscript.
apply_mediator_measurement_error <- function(x_true, reliability = 1) {
  if (
    !is.numeric(reliability) ||
    length(reliability) != 1 ||
    reliability <= 0 ||
    reliability > 1
  ) {
    stop("reliability must be in (0, 1].", call. = FALSE)
  }

  if (isTRUE(all.equal(reliability, 1))) {
    return(as.numeric(x_true))
  }

  error_sd <- sqrt((1 - reliability) / reliability)
  x_observed <- x_true + stats::rnorm(
    length(x_true),
    mean = 0,
    sd = error_sd
  )

  standardize(x_observed)
}


#' Simulate one synthetic cohort and attach replicate-specific analysis truth.
#'
#' The outcome is generated from the latent true mediators. Measurement error,
#' when requested, is added only after outcome generation. This makes the
#' measurement-error stress test a direct assessment of mediator attenuation.
#'
#' The fitted SEM uses standardized G, H, M, and Y. Therefore raw generating
#' coefficients are transformed to the corresponding standardized analysis
#' scale before bias, RMSE, or coverage are calculated.
simulate_digital_twin_cohort <- function(
    n = 5000,
    scenario = "neural_dominant",
    rho_hm = 0.60,
    mediator_reliability = 1.00,
    seed = 1
) {
  set.seed(seed)

  pars <- scenario_parameters(scenario, rho_hm)
  anchors <- publication_anchors()

  age_years <- pmax(
    anchors$age_min,
    pmin(
      stats::rnorm(n, anchors$age_mean, anchors$age_sd),
      anchors$age_max
    )
  )

  bmi_value <- pmax(
    anchors$bmi_min,
    pmin(
      stats::rnorm(n, anchors$bmi_mean, anchors$bmi_sd),
      anchors$bmi_max
    )
  )

  age <- standardize(age_years)
  bmi <- standardize(bmi_value)

  treatment <- stats::rbinom(n, 1, 0.5)
  pgs <- stats::rnorm(n, 0, 1)

  homa_ir <- pmax(
    0.1,
    anchors$homa_mean +
      anchors$homa_sd * bmi +
      stats::rnorm(n, 0, 0.5)
  )
  homa_std <- standardize(homa_ir)

  g_latent <- (
    pars$alpha_t * treatment -
      0.20 * bmi -
      0.15 * homa_std +
      0.15 * pgs +
      stats::rnorm(n, 0, 1)
  )
  g <- standardize(g_latent)

  residual_covariance <- matrix(
    c(1, pars$rho_hm, pars$rho_hm, 1),
    nrow = 2,
    byrow = TRUE
  )

  mediator_errors <- MASS::mvrnorm(
    n,
    mu = c(0, 0),
    Sigma = residual_covariance
  )

  h_latent <- pars$beta_g * g - 0.10 * age + mediator_errors[, 1]
  m_latent <- (
    pars$gamma_g * g -
      0.30 * bmi -
      0.20 * homa_std +
      mediator_errors[, 2]
  )

  sd_h_latent <- stats::sd(h_latent)
  sd_m_latent <- stats::sd(m_latent)

  h_true <- standardize(h_latent)
  m_true <- standardize(m_latent)

  receptor_sensitivity <- ifelse(
    treatment == 1,
    stats::rnorm(n, 1.3, 0.1),
    stats::rnorm(n, 1.0, 0.1)
  )
  sensitivity_std <- standardize(receptor_sensitivity)
  h_by_sensitivity_true <- h_true * sensitivity_std

  y_raw <- (
    pars$delta_h * h_true +
      pars$delta_m * m_true +
      ifelse(
        scenario == "heterogeneous",
        pars$delta_hs * h_by_sensitivity_true,
        0
      ) +
      0.10 * age +
      stats::rnorm(n, 0, 1)
  )

  sd_y_raw <- stats::sd(y_raw)
  y <- standardize(y_raw)

  h_observed <- apply_mediator_measurement_error(
    h_true,
    mediator_reliability
  )
  m_observed <- apply_mediator_measurement_error(
    m_true,
    mediator_reliability
  )

  # Analysis-scale truth.
  a_true <- pars$beta_g / sd_h_latent
  b_true <- pars$gamma_g / sd_m_latent
  c_true <- pars$delta_h / sd_y_raw
  d_true <- pars$delta_m / sd_y_raw

  psi_h_true <- a_true * c_true
  psi_m_true <- b_true * d_true
  delta_true <- psi_h_true - psi_m_true

  realized_reliability_h <- suppressWarnings(
    stats::cor(h_true, h_observed, use = "complete.obs")^2
  )
  realized_reliability_m <- suppressWarnings(
    stats::cor(m_true, m_observed, use = "complete.obs")^2
  )

  cohort <- tibble::tibble(
    id = seq_len(n),
    treatment = treatment,
    age_years = age_years,
    bmi_value = bmi_value,
    homa_ir = homa_ir,
    homa_std = homa_std,
    pgs = pgs,
    G = g,
    H_true = h_true,
    M_true = m_true,
    H = h_observed,
    M = m_observed,
    Y = y,
    hrv_ms = anchors$hrv_mean + anchors$hrv_sd * h_observed,
    cognition_t = anchors$cognition_mean + anchors$cognition_sd * y,
    Age = age,
    BMI = bmi,
    scenario = scenario,
    rho_hm = rho_hm,
    mediator_reliability = mediator_reliability,
    realized_reliability_h = realized_reliability_h,
    realized_reliability_m = realized_reliability_m
  )

  attr(cohort, "analysis_truth") <- tibble::tibble(
    a_true = a_true,
    b_true = b_true,
    c_true = c_true,
    d_true = d_true,
    psi_h_true = psi_h_true,
    psi_m_true = psi_m_true,
    delta_true = delta_true,
    sd_h_latent = sd_h_latent,
    sd_m_latent = sd_m_latent,
    sd_y_raw = sd_y_raw,
    realized_reliability_h = realized_reliability_h,
    realized_reliability_m = realized_reliability_m
  )

  cohort
}


# ==============================================================================
# 5. Structural equation model
# ==============================================================================

#' Return the prespecified parallel-pathway SEM.
pathway_sem_model <- function() {
  "
  G ~ treatment + BMI + pgs
  H ~ a*G + Age
  M ~ b*G + BMI
  Y ~ c*H + d*M + Age

  H ~~ cov_hm*M

  psi_H := a*c
  psi_M := b*d
  pathway_sum := psi_H + psi_M
  Delta := psi_H - psi_M
  "
}


#' Fit the pathway SEM and extract estimates used in the Monte Carlo analysis.
fit_pathway_sem <- function(cohort, ci_level = 0.95) {
  fit <- lavaan::sem(
    pathway_sem_model(),
    data = cohort,
    fixed.x = TRUE,
    meanstructure = FALSE
  )

  converged <- isTRUE(lavaan::lavInspect(fit, "converged"))

  fit_indices <- tryCatch(
    lavaan::fitMeasures(fit, c("cfi", "tli", "rmsea", "srmr")),
    error = function(e) {
      c(cfi = NA_real_, tli = NA_real_, rmsea = NA_real_, srmr = NA_real_)
    }
  )

  estimates <- lavaan::parameterEstimates(
    fit,
    standardized = TRUE,
    ci = TRUE,
    level = ci_level
  )

  extract_path <- function(lhs_name, rhs_name) {
    row <- estimates |>
      dplyr::filter(
        op == "~",
        lhs == lhs_name,
        rhs == rhs_name
      ) |>
      dplyr::slice(1)

    if (nrow(row) == 0) {
      return(NA_real_)
    }

    row$est[1]
  }

  extract_covariance <- function(lhs_name, rhs_name) {
    row <- estimates |>
      dplyr::filter(
        op == "~~",
        lhs == lhs_name,
        rhs == rhs_name
      ) |>
      dplyr::slice(1)

    if (nrow(row) == 0) {
      return(NA_real_)
    }

    row$est[1]
  }

  extract_defined_parameter <- function(parameter_name) {
    row <- estimates |>
      dplyr::filter(
        op == ":=",
        lhs == parameter_name
      ) |>
      dplyr::slice(1)

    if (nrow(row) == 0) {
      return(
        tibble::tibble(
          estimate = NA_real_,
          se = NA_real_,
          p_value = NA_real_,
          lower = NA_real_,
          upper = NA_real_
        )
      )
    }

    row |>
      dplyr::transmute(
        estimate = est,
        se = se,
        p_value = pvalue,
        lower = ci.lower,
        upper = ci.upper
      )
  }

  neural <- extract_defined_parameter("psi_H")
  metabolic <- extract_defined_parameter("psi_M")
  pathway_sum <- extract_defined_parameter("pathway_sum")
  delta <- extract_defined_parameter("Delta")

  tibble::tibble(
    converged = converged,
    cfi = unname(fit_indices["cfi"]),
    tli = unname(fit_indices["tli"]),
    rmsea = unname(fit_indices["rmsea"]),
    srmr = unname(fit_indices["srmr"]),

    a_hat = extract_path("H", "G"),
    b_hat = extract_path("M", "G"),
    c_hat = extract_path("Y", "H"),
    d_hat = extract_path("Y", "M"),
    covariance_hm_hat = extract_covariance("H", "M"),

    psi_h_hat = neural$estimate,
    psi_h_se = neural$se,
    psi_h_p = neural$p_value,
    psi_h_lower = neural$lower,
    psi_h_upper = neural$upper,

    psi_m_hat = metabolic$estimate,
    psi_m_se = metabolic$se,
    psi_m_p = metabolic$p_value,
    psi_m_lower = metabolic$lower,
    psi_m_upper = metabolic$upper,

    pathway_sum_hat = pathway_sum$estimate,

    delta_hat = delta$estimate,
    delta_se = delta$se,
    delta_p = delta$p_value,
    delta_lower = delta$lower,
    delta_upper = delta$upper
  )
}


# ==============================================================================
# 6. Primary Monte Carlo operating-characteristic analysis
# ==============================================================================

#' Construct the primary, sample-size, and measurement-error design grids.
build_simulation_design <- function(
    scenarios = SCENARIOS,
    primary_rhos = PRIMARY_RHOS,
    sample_sizes = SAMPLE_SIZES,
    reliability_levels = MEDIATOR_RELIABILITIES,
    stress_test_rho = STRESS_TEST_RHO
) {
  primary <- tidyr::crossing(
    analysis = "Primary correlation stress test",
    scenario = scenarios,
    rho_hm = primary_rhos,
    n_subjects = 5000L,
    mediator_reliability = 1.00
  )

  sample_size <- tidyr::crossing(
    analysis = "Sample-size stress test",
    scenario = scenarios,
    rho_hm = stress_test_rho,
    n_subjects = sample_sizes,
    mediator_reliability = 1.00
  )

  measurement_error <- tidyr::crossing(
    analysis = "Mediator-measurement-error stress test",
    scenario = scenarios,
    rho_hm = stress_test_rho,
    n_subjects = 5000L,
    mediator_reliability = reliability_levels
  )

  list(
    primary = primary,
    sample_size = sample_size,
    measurement_error = measurement_error,
    all = dplyr::bind_rows(
      primary,
      sample_size,
      measurement_error
    )
  )
}


#' Run the SEM Monte Carlo simulation over a supplied design grid.
run_operating_characteristic_simulation <- function(
    design,
    n_mc = 500,
    seed_base = 91001,
    ci_level = 0.95,
    progress_every = 25
) {
  required_columns <- c(
    "analysis",
    "scenario",
    "rho_hm",
    "n_subjects",
    "mediator_reliability"
  )

  if (!all(required_columns %in% names(design))) {
    stop(
      "design is missing one or more required columns: ",
      paste(required_columns, collapse = ", "),
      call. = FALSE
    )
  }

  results <- vector("list", nrow(design) * n_mc)
  result_index <- 0L

  for (design_index in seq_len(nrow(design))) {
    current <- design[design_index, ]

    for (replicate in seq_len(n_mc)) {
      result_index <- result_index + 1L
      replicate_seed <- (
        seed_base +
          100000L * design_index +
          replicate
      )

      if (
        progress_every > 0 &&
        (replicate == 1 || replicate %% progress_every == 0)
      ) {
        message(
          sprintf(
            paste0(
              "Monte Carlo grid %d/%d | replicate %d/%d | ",
              "%s | scenario=%s | rho=%.2f | N=%d | reliability=%.2f"
            ),
            design_index,
            nrow(design),
            replicate,
            n_mc,
            current$analysis,
            current$scenario,
            current$rho_hm,
            current$n_subjects,
            current$mediator_reliability
          )
        )
      }

      cohort <- simulate_digital_twin_cohort(
        n = current$n_subjects,
        scenario = current$scenario,
        rho_hm = current$rho_hm,
        mediator_reliability = current$mediator_reliability,
        seed = replicate_seed
      )

      truth <- attr(cohort, "analysis_truth")
      sem_result <- fit_pathway_sem(
        cohort,
        ci_level = ci_level
      )

      results[[result_index]] <- dplyr::bind_cols(
        tibble::tibble(
          analysis = current$analysis,
          scenario = current$scenario,
          scenario_label = scenario_label(current$scenario),
          rho_hm = current$rho_hm,
          n_subjects = current$n_subjects,
          mediator_reliability = current$mediator_reliability,
          replicate = replicate,
          seed = replicate_seed
        ),
        truth,
        sem_result
      )
    }
  }

  dplyr::bind_rows(results)
}


#' Summarize pathway-specific and Delta operating characteristics.
summarize_operating_characteristics <- function(
    simulation_runs,
    alpha = 0.05
) {
  simulation_runs |>
    dplyr::mutate(
      error_h = psi_h_hat - psi_h_true,
      error_m = psi_m_hat - psi_m_true,
      error_delta = delta_hat - delta_true,

      cover_h = (
        psi_h_lower <= psi_h_true &
          psi_h_upper >= psi_h_true
      ),
      cover_m = (
        psi_m_lower <= psi_m_true &
          psi_m_upper >= psi_m_true
      ),
      cover_delta = (
        delta_lower <= delta_true &
          delta_upper >= delta_true
      ),

      correct_discrimination = dplyr::case_when(
        scenario == "neural_dominant" ~ delta_hat > 0,
        scenario == "metabolic_dominant" ~ delta_hat < 0,
        TRUE ~ NA
      ),

      significant_correct_discrimination = dplyr::case_when(
        scenario == "neural_dominant" ~ delta_lower > 0,
        scenario == "metabolic_dominant" ~ delta_upper < 0,
        TRUE ~ NA
      ),

      reject_delta_zero = (
        !is.na(delta_p) &
          delta_p < alpha
      )
    ) |>
    dplyr::group_by(
      analysis,
      scenario,
      scenario_label,
      rho_hm,
      n_subjects,
      mediator_reliability
    ) |>
    dplyr::summarise(
      n_runs = dplyr::n(),
      convergence = mean(converged, na.rm = TRUE),

      psi_h_truth = mean(psi_h_true, na.rm = TRUE),
      psi_m_truth = mean(psi_m_true, na.rm = TRUE),
      delta_truth = mean(delta_true, na.rm = TRUE),

      psi_h_bias = mean(error_h, na.rm = TRUE),
      psi_h_mae = mean(abs(error_h), na.rm = TRUE),
      psi_h_rmse = sqrt(mean(error_h^2, na.rm = TRUE)),
      psi_h_coverage = mean(cover_h, na.rm = TRUE),

      psi_m_bias = mean(error_m, na.rm = TRUE),
      psi_m_mae = mean(abs(error_m), na.rm = TRUE),
      psi_m_rmse = sqrt(mean(error_m^2, na.rm = TRUE)),
      psi_m_coverage = mean(cover_m, na.rm = TRUE),

      delta_bias = mean(error_delta, na.rm = TRUE),
      delta_mae = mean(abs(error_delta), na.rm = TRUE),
      delta_rmse = sqrt(mean(error_delta^2, na.rm = TRUE)),
      delta_coverage = mean(cover_delta, na.rm = TRUE),

      delta_type1_error = dplyr::if_else(
        dplyr::first(scenario) == "null",
        mean(reject_delta_zero, na.rm = TRUE),
        NA_real_
      ),

      delta_power = dplyr::if_else(
        dplyr::first(scenario) %in% c(
          "neural_dominant",
          "metabolic_dominant",
          "heterogeneous"
        ),
        mean(reject_delta_zero, na.rm = TRUE),
        NA_real_
      ),

      probability_correct_discrimination = dplyr::if_else(
        dplyr::first(scenario) %in% c(
          "neural_dominant",
          "metabolic_dominant"
        ),
        mean(correct_discrimination, na.rm = TRUE),
        NA_real_
      ),

      probability_significant_correct_discrimination = dplyr::if_else(
        dplyr::first(scenario) %in% c(
          "neural_dominant",
          "metabolic_dominant"
        ),
        mean(significant_correct_discrimination, na.rm = TRUE),
        NA_real_
      ),

      probability_neural_gt_metabolic = mean(
        delta_hat > 0,
        na.rm = TRUE
      ),

      realized_reliability_h = mean(
        realized_reliability_h,
        na.rm = TRUE
      ),
      realized_reliability_m = mean(
        realized_reliability_m,
        na.rm = TRUE
      ),

      .groups = "drop"
    )
}


# ==============================================================================
# 7. Secondary validation analyses
# ==============================================================================

#' Fit the individual-level 2SLS genetic-instrument model.
fit_genetic_2sls <- function(cohort) {
  model <- AER::ivreg(
    Y ~ G + Age + BMI | pgs + Age + BMI,
    data = cohort
  )

  model_summary <- summary(
    model,
    diagnostics = TRUE
  )

  tibble::tibble(
    iv_estimate = unname(
      model_summary$coefficients["G", "Estimate"]
    ),
    iv_se = unname(
      model_summary$coefficients["G", "Std. Error"]
    ),
    iv_p = unname(
      model_summary$coefficients["G", "Pr(>|t|)"]
    ),
    first_stage_f = as.numeric(
      model_summary$diagnostics["Weak instruments", "statistic"]
    )
  )
}


#' Permute the observed mediators and refit the SEM as a falsification test.
run_permutation_falsification <- function(cohort) {
  permuted <- cohort
  permuted$H <- cohort$H[sample(nrow(cohort))]
  permuted$M <- cohort$M[sample(nrow(cohort))]

  result <- fit_pathway_sem(permuted)

  tibble::tibble(
    psi_h_permuted = result$psi_h_hat,
    psi_m_permuted = result$psi_m_hat,
    pathway_sum_permuted = result$pathway_sum_hat
  )
}


#' Run the secondary SEM, 2SLS, and permutation simulations.
run_secondary_simulation <- function(
    n_mc = 50,
    n_subjects = 5000,
    scenarios = SCENARIOS,
    rhos = PRIMARY_RHOS,
    seed_base = 7001,
    progress = TRUE
) {
  total <- length(scenarios) * length(rhos) * n_mc
  results <- vector("list", total)
  index <- 0L

  for (scenario in scenarios) {
    for (rho_hm in rhos) {
      for (replicate in seq_len(n_mc)) {
        index <- index + 1L
        seed <- seed_base + index

        if (isTRUE(progress) && (replicate == 1 || replicate %% 10 == 0)) {
          message(
            sprintf(
              "Secondary analysis %d/%d | %s | rho=%.2f | replicate=%d",
              index,
              total,
              scenario,
              rho_hm,
              replicate
            )
          )
        }

        cohort <- simulate_digital_twin_cohort(
          n = n_subjects,
          scenario = scenario,
          rho_hm = rho_hm,
          mediator_reliability = 1.00,
          seed = seed
        )

        sem_result <- fit_pathway_sem(cohort)
        iv_result <- fit_genetic_2sls(cohort)
        permutation_result <- run_permutation_falsification(cohort)

        results[[index]] <- dplyr::bind_cols(
          tibble::tibble(
            scenario = scenario,
            scenario_label = scenario_label(scenario),
            rho_hm = rho_hm,
            replicate = replicate,
            seed = seed
          ),
          sem_result,
          iv_result,
          permutation_result
        )
      }
    }
  }

  dplyr::bind_rows(results)
}


#' Simulate GWAS-summary SNP associations for one outcome.
simulate_gwas_summary_data <- function(
    n_snps = 50,
    true_beta = 0.18,
    sample_size = 269867,
    pleiotropy_eta = 0,
    seed = 4001,
    outcome = "Cognition"
) {
  set.seed(seed)

  bx <- abs(stats::rnorm(
    n_snps,
    mean = 0.10,
    sd = 0.02
  ))

  pleiotropy_term <- stats::rnorm(
    n_snps,
    mean = 0,
    sd = 0.10
  )

  se_scale <- sqrt(250000 / sample_size)
  se <- stats::runif(
    n_snps,
    min = 0.005,
    max = 0.015
  ) * se_scale

  by <- (
    bx * true_beta +
      pleiotropy_eta * pleiotropy_term +
      stats::rnorm(n_snps, 0, mean(se))
  )

  tibble::tibble(
    snp = paste0("rs", seq_len(n_snps)),
    outcome = outcome,
    sample_size = sample_size,
    bx = bx,
    by = by,
    se = se,
    weight = 1 / se^2
  )
}


#' Fit inverse-variance weighted regression through the origin.
fit_ivw <- function(gwas_data) {
  model <- stats::lm(
    by ~ 0 + bx,
    data = gwas_data,
    weights = weight
  )
  model_summary <- summary(model)

  estimate <- model_summary$coefficients["bx", "Estimate"]
  standard_error <- model_summary$coefficients["bx", "Std. Error"]
  p_value <- model_summary$coefficients["bx", "Pr(>|t|)"]
  residual <- gwas_data$by - estimate * gwas_data$bx
  q_statistic <- sum(gwas_data$weight * residual^2)

  tibble::tibble(
    method = "IVW",
    beta = estimate,
    se = standard_error,
    p = p_value,
    Q = q_statistic,
    egger_intercept = NA_real_,
    egger_intercept_se = NA_real_,
    egger_intercept_p = NA_real_
  )
}


#' Fit weighted MR-Egger regression.
fit_mr_egger <- function(gwas_data) {
  model <- stats::lm(
    by ~ bx,
    data = gwas_data,
    weights = weight
  )
  model_summary <- summary(model)

  tibble::tibble(
    method = "MR-Egger",
    beta = model_summary$coefficients["bx", "Estimate"],
    se = model_summary$coefficients["bx", "Std. Error"],
    p = model_summary$coefficients["bx", "Pr(>|t|)"],
    Q = NA_real_,
    egger_intercept = model_summary$coefficients[
      "(Intercept)",
      "Estimate"
    ],
    egger_intercept_se = model_summary$coefficients[
      "(Intercept)",
      "Std. Error"
    ],
    egger_intercept_p = model_summary$coefficients[
      "(Intercept)",
      "Pr(>|t|)"
    ]
  )
}


#' Run the complete GWAS-summary genetic-instrument analysis.
run_gwas_summary_analysis <- function(
    seed_base = 4001,
    pleiotropy_eta = 0
) {
  anchors <- publication_anchors()
  parameters <- scenario_parameters(
    scenario = "neural_dominant",
    rho_hm = 0.60
  )

  outcomes <- names(parameters$gwas_true_betas)

  gwas_data <- dplyr::bind_rows(
    lapply(seq_along(outcomes), function(i) {
      outcome <- outcomes[i]

      simulate_gwas_summary_data(
        n_snps = 50,
        true_beta = parameters$gwas_true_betas[outcome],
        sample_size = unname(
          anchors$gwas_sample_sizes[outcome]
        ),
        pleiotropy_eta = pleiotropy_eta,
        seed = seed_base + i,
        outcome = outcome
      )
    })
  )

  mr_summary <- gwas_data |>
    dplyr::group_by(outcome, sample_size) |>
    dplyr::group_modify(
      ~ dplyr::bind_rows(
        fit_ivw(.x),
        fit_mr_egger(.x)
      )
    ) |>
    dplyr::ungroup()

  list(
    gwas_data = gwas_data,
    mr_summary = mr_summary
  )
}


# ==============================================================================
# 8. Manuscript tables
# ==============================================================================

#' Export the simulation parameter table corresponding to Supplementary Table S2.
export_simulation_parameter_table <- function(output_dir) {
  ensure_directory(output_dir)

  parameter_table <- tibble::tribble(
    ~component, ~parameter, ~description, ~implemented_value,

    "Baseline covariate", "Age",
    "Normal draw, clipped to the specified age range, then standardized",
    "Z_Age ~ N(71.8, 7.1^2); Age_Real = max(55, min(Z_Age, 90)); Age = scale(Age_Real)",

    "Baseline covariate", "BMI",
    "Normal draw, clipped to the specified BMI range, then standardized",
    "Z_BMI ~ N(26.5, 5.2^2); BMI_Real = max(16, min(Z_BMI, 55)); BMI = scale(BMI_Real)",

    "Treatment", "treatment",
    "Randomized treatment assignment",
    "Bernoulli(0.5)",

    "Genetic instrument", "pgs",
    "Simulated polygenic score",
    "N(0,1)",

    "Metabolic burden", "epsilon_HOMA",
    "Residual used in the HOMA-like variable",
    "N(0, 0.5^2); HOMA_IR = max(0.1, 2.6 + 1.8*BMI + epsilon_HOMA)",

    "GLP-1 signaling", "alpha_t",
    "Treatment to GLP-1 signaling coefficient",
    "0.90",

    "GLP-1 signaling", "epsilon_G",
    "GLP-1 signaling residual",
    "N(0,1); variance = 1",

    "Neural pathway", "beta_g",
    "GLP-1 signaling to neural mediator coefficient",
    "0.60",

    "Metabolic pathway", "gamma_g",
    "GLP-1 signaling to metabolic mediator coefficient",
    "0.60",

    "Mediator residuals", "epsilon_H, epsilon_M",
    "Bivariate-normal neural and metabolic mediator residuals",
    "N2((0,0)', Sigma); Sigma = [[1, rho_HM], [rho_HM, 1]]; marginal residual variances = 1",

    "Mediator overlap", "rho_hm",
    "Residual neural-metabolic mediator correlation",
    "0.2, 0.6, 0.8",

    "Receptor sensitivity", "Sensitivity | treatment=1",
    "Receptor sensitivity in the treated group before standardization",
    "N(1.3, 0.1^2)",

    "Receptor sensitivity", "Sensitivity | treatment=0",
    "Receptor sensitivity in the control group before standardization",
    "N(1.0, 0.1^2)",

    "Receptor sensitivity", "S_std",
    "Standardized receptor sensitivity",
    "scale(Sensitivity)",

    "Outcome", "epsilon_Y",
    "Cognition residual",
    "N(0,1); variance = 1",

    "Neural-dominant", "delta_h / delta_m",
    "Mediator-to-cognition coefficients",
    "0.45 / 0.15",

    "Metabolic-dominant", "delta_h / delta_m",
    "Mediator-to-cognition coefficients",
    "0.08 / 0.35",

    "Shared", "delta_h / delta_m",
    "Mediator-to-cognition coefficients",
    "0.25 / 0.25",

    "Heterogeneous", "delta_h / delta_m",
    "Mediator-to-cognition coefficients",
    "0.22 / 0.18",

    "Heterogeneous", "delta_hs",
    "Neural mediator by receptor-sensitivity coefficient",
    "0.22",

    "Null", "delta_h / delta_m",
    "Mediator-to-cognition coefficients",
    "0 / 0",

    "Primary design", "N",
    "Simulated cohort size",
    "5000",

    "Primary design", "N_MC",
    "Monte Carlo replicates per scenario x rho_hm condition",
    "500",

    "Sample-size stress test", "N",
    "Simulated cohort sizes",
    "500, 1000, 2500, 5000",

    "Measurement-error stress test", "R",
    "Mediator reliability",
    "1.00, 0.80, 0.60",

    "GWAS summary", "beta_true",
    "Programmed causal slopes",
    "Cognition 0.18; Education 0.12; Alzheimer -0.22"
  )

  readr::write_csv(
    parameter_table,
    file.path(
      output_dir,
      "table_s2_simulation_parameters.csv"
    )
  )

  invisible(parameter_table)
}

#' Create and export tables used in the manuscript and supplement.
export_manuscript_tables <- function(
    primary_runs,
    operating_characteristics,
    secondary_runs,
    output_dir
) {
  data_dir <- file.path(output_dir, "data")
  ensure_directory(data_dir)

  parameter_table <- export_simulation_parameter_table(
    data_dir
  )

  # Replicate-level and complete operating-characteristic files.
  readr::write_csv(
    primary_runs,
    file.path(data_dir, "primary_mc_runs.csv")
  )

  readr::write_csv(
    operating_characteristics,
    file.path(data_dir, "operating_characteristics_all.csv")
  )

  readr::write_csv(
    secondary_runs,
    file.path(data_dir, "secondary_mc_runs.csv")
  )

  primary_summary <- operating_characteristics |>
    dplyr::filter(
      analysis == "Primary correlation stress test"
    )

  sample_size_summary <- operating_characteristics |>
    dplyr::filter(
      analysis == "Sample-size stress test"
    )

  measurement_error_summary <- operating_characteristics |>
    dplyr::filter(
      analysis == "Mediator-measurement-error stress test"
    )

  # Main Table 1: compact ranges across mediator-correlation levels.
  table_1 <- primary_summary |>
    dplyr::group_by(
      scenario,
      scenario_label
    ) |>
    dplyr::summarise(
      delta_truth_min = min(delta_truth, na.rm = TRUE),
      delta_truth_max = max(delta_truth, na.rm = TRUE),
      delta_bias_min = min(delta_bias, na.rm = TRUE),
      delta_bias_max = max(delta_bias, na.rm = TRUE),
      delta_rmse_min = min(delta_rmse, na.rm = TRUE),
      delta_rmse_max = max(delta_rmse, na.rm = TRUE),
      coverage_min = min(delta_coverage, na.rm = TRUE),
      coverage_max = max(delta_coverage, na.rm = TRUE),
      power_min = if (
        all(is.na(delta_power))
      ) NA_real_ else min(delta_power, na.rm = TRUE),
      power_max = if (
        all(is.na(delta_power))
      ) NA_real_ else max(delta_power, na.rm = TRUE),
      type1_error_min = if (
        all(is.na(delta_type1_error))
      ) NA_real_ else min(delta_type1_error, na.rm = TRUE),
      type1_error_max = if (
        all(is.na(delta_type1_error))
      ) NA_real_ else max(delta_type1_error, na.rm = TRUE),
      discrimination_min = if (
        all(is.na(probability_correct_discrimination))
      ) NA_real_ else min(
        probability_correct_discrimination,
        na.rm = TRUE
      ),
      discrimination_max = if (
        all(is.na(probability_correct_discrimination))
      ) NA_real_ else max(
        probability_correct_discrimination,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  readr::write_csv(
    table_1,
    file.path(
      data_dir,
      "table_1_delta_operating_characteristics_summary.csv"
    )
  )

  # Supplementary Table S3: complete primary correlation-stress results.
  readr::write_csv(
    primary_summary,
    file.path(
      data_dir,
      "table_s3_primary_delta_operating_characteristics.csv"
    )
  )

  # Supplementary Table S5.
  pathway_recovery <- operating_characteristics |>
    dplyr::select(
      analysis,
      scenario,
      scenario_label,
      rho_hm,
      n_subjects,
      mediator_reliability,
      n_runs,
      psi_h_truth,
      psi_h_bias,
      psi_h_mae,
      psi_h_rmse,
      psi_h_coverage,
      psi_m_truth,
      psi_m_bias,
      psi_m_mae,
      psi_m_rmse,
      psi_m_coverage
    )

  readr::write_csv(
    pathway_recovery,
    file.path(
      data_dir,
      "table_s5_pathway_recovery_analysis_scale.csv"
    )
  )

  # Supplementary Tables S6 and S7.
  readr::write_csv(
    sample_size_summary,
    file.path(
      data_dir,
      "table_s6_sample_size_stress_test.csv"
    )
  )

  readr::write_csv(
    measurement_error_summary,
    file.path(
      data_dir,
      "table_s7_measurement_error_stress_test.csv"
    )
  )

  # Secondary individual-level 2SLS output: main Table 2 / Table S8.
  table_2sls <- secondary_runs |>
    dplyr::group_by(
      scenario,
      scenario_label,
      rho_hm
    ) |>
    dplyr::summarise(
      n_runs = dplyr::n(),
      iv_estimate = format_mean_sd(iv_estimate),
      iv_se = format_mean_sd(iv_se),
      first_stage_f = format_mean_sd(first_stage_f),
      .groups = "drop"
    )

  readr::write_csv(
    table_2sls,
    file.path(
      data_dir,
      "table_s8_individual_2sls.csv"
    )
  )

  # Main Table 2: compact ranges across mediator-correlation levels.
  table_2 <- secondary_runs |>
    dplyr::group_by(
      scenario,
      scenario_label
    ) |>
    dplyr::summarise(
      iv_estimate_min = min(iv_estimate, na.rm = TRUE),
      iv_estimate_max = max(iv_estimate, na.rm = TRUE),
      iv_se_min = min(iv_se, na.rm = TRUE),
      iv_se_max = max(iv_se, na.rm = TRUE),
      first_stage_f_min = min(first_stage_f, na.rm = TRUE),
      first_stage_f_max = max(first_stage_f, na.rm = TRUE),
      .groups = "drop"
    )

  readr::write_csv(
    table_2,
    file.path(
      data_dir,
      "table_2_individual_2sls_summary.csv"
    )
  )

  # Supplementary Table S4.
  permutation_summary <- secondary_runs |>
    dplyr::group_by(
      scenario,
      scenario_label,
      rho_hm
    ) |>
    dplyr::summarise(
      neural_permuted = format_mean_sd(psi_h_permuted),
      metabolic_permuted = format_mean_sd(psi_m_permuted),
      sum_permuted = format_mean_sd(pathway_sum_permuted),
      .groups = "drop"
    )

  readr::write_csv(
    permutation_summary,
    file.path(
      data_dir,
      "table_s4_permutation_falsification.csv"
    )
  )

  # Explicit interpretation table for discrimination.
  discrimination_interpretation <- primary_summary |>
    dplyr::transmute(
      scenario_code = scenario,
      scenario = scenario_label,
      rho_hm = rho_hm,
      delta_truth = delta_truth,
      delta_bias = delta_bias,
      delta_rmse = delta_rmse,
      delta_coverage = delta_coverage,
      type1_error = delta_type1_error,
      power = delta_power,
      probability_correct_discrimination =
        probability_correct_discrimination,
      interpretation = dplyr::case_when(
        scenario_code == "neural_dominant" ~
          "Correct discrimination: estimated Delta > 0",
        scenario_code == "metabolic_dominant" ~
          "Correct discrimination: estimated Delta < 0",
        scenario_code == "shared" ~
          "No prespecified dominant pathway; classification not defined",
        scenario_code == "null" ~
          "No prespecified dominant pathway; report type I error",
        TRUE ~
          "Heterogeneous scenario; summarize Delta operating characteristics"
      )
    )

  readr::write_csv(
    discrimination_interpretation,
    file.path(
      data_dir,
      "discrimination_interpretation.csv"
    )
  )

  invisible(
    list(
      table_1 = table_1,
      table_2 = table_2,
      table_s2 = parameter_table,
      primary = primary_summary,
      sample_size = sample_size_summary,
      measurement_error = measurement_error_summary,
      pathway_recovery = pathway_recovery,
      two_stage_least_squares = table_2sls,
      permutation = permutation_summary
    )
  )
}


# ==============================================================================
# 9. Manuscript figures
# ==============================================================================

#' Main Figure 1: prespecified causal architecture.
#'
#' Generated with ggplot2 so the complete computational figure set is
#' reproducible directly from this repository.
plot_figure_1_causal_architecture <- function(
    output_file
) {
  ensure_directory(dirname(output_file))

  nodes <- tibble::tribble(
    ~node, ~x, ~y, ~label,
    "T",   0.0, 1.0, "Treatment\n(T)",
    "G",   1.5, 1.0, "GLP-1 signaling\n(G)",
    "H",   3.0, 1.6, "Neural mediator\n(H: HRV)",
    "M",   3.0, 0.4, "Metabolic mediator\n(M)",
    "Y",   4.6, 1.0, "Cognition\n(Y)",
    "PGS", 1.5, 2.0, "Polygenic score\n(PGS)"
  )

  edges <- tibble::tribble(
    ~x, ~y, ~xend, ~yend,
    0.3, 1.0, 1.15, 1.0,
    1.85, 1.0, 2.65, 1.52,
    1.85, 1.0, 2.65, 0.48,
    3.35, 1.55, 4.25, 1.05,
    3.35, 0.45, 4.25, 0.95,
    1.5, 1.72, 1.5, 1.28
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(
        x = x, y = y, xend = xend, yend = yend
      ),
      linewidth = 0.8,
      arrow = grid::arrow(
        length = grid::unit(0.18, "inches"),
        type = "closed"
      )
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 2.70, y = 1.30,
        xend = 2.70, yend = 0.70
      ),
      linetype = 2,
      linewidth = 0.8
    ) +
    ggplot2::annotate(
      "text",
      x = 2.83,
      y = 1.00,
      label = expression(rho[HM]),
      hjust = 0,
      size = 4.3
    ) +
    ggplot2::geom_label(
      data = nodes,
      ggplot2::aes(x = x, y = y, label = label),
      size = 4.0,
      label.size = 0.5,
      label.padding = grid::unit(0.20, "lines"),
      fill = "white"
    ) +
    ggplot2::coord_cartesian(
      xlim = c(-0.4, 5.1),
      ylim = c(0.0, 2.35),
      clip = "off"
    ) +
    ggplot2::labs(
      title = "Prespecified causal architecture for pathway analysis"
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      )
    )

  ggplot2::ggsave(
    output_file,
    p,
    width = 10,
    height = 5,
    dpi = 300,
    bg = "white"
  )

  invisible(p)
}


#' Main Figure 2: estimated pathway contrast across mediator correlation.
#'
#' Reproduces the current manuscript layout: five data-generating scenarios
#' on the x-axis with side-by-side boxplots for rho_HM = 0.2, 0.6, and 0.8.
#' Positive Delta values favor the neural pathway; negative values favor the
#' metabolic pathway. Shared and null scenarios are descriptive and are not
#' treated as classification problems.
plot_figure_2_pathway_contrast <- function(
    primary_runs,
    output_file
) {
  ensure_directory(dirname(output_file))

  scenario_levels <- c(
    "Neural dominant",
    "Metabolic dominant",
    "Heterogeneous",
    "Shared",
    "Null"
  )

  plot_data <- primary_runs |>
    dplyr::filter(
      analysis == "Primary correlation stress test"
    ) |>
    dplyr::mutate(
      scenario_label = factor(
        scenario_label,
        levels = scenario_levels
      ),
      rho_group = factor(
        rho_hm,
        levels = c(0.2, 0.6, 0.8),
        labels = c("0.2", "0.6", "0.8")
      )
    )

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = scenario_label,
      y = delta_hat,
      fill = rho_group,
      group = interaction(
        scenario_label,
        rho_group,
        drop = TRUE
      )
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = 2,
      linewidth = 0.7
    ) +
    ggplot2::geom_boxplot(
      position = ggplot2::position_dodge2(
        width = 0.80,
        preserve = "single",
        padding = 0.15
      ),
      width = 0.24,
      alpha = 0.35,
      outlier.alpha = 0.12,
      outlier.size = 1.0
    ) +
    ggplot2::labs(
      title = "Estimated pathway contrast across mediator correlation",
      x = "Data-generating scenario",
      y = expression(
        hat(Delta) == hat(psi)[H] - hat(psi)[M]
      ),
      fill = expression(rho[HM])
    ) +
    publication_theme(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(
        angle = 0,
        hjust = 0.5
      )
    )

  ggplot2::ggsave(
    output_file,
    p,
    width = 12,
    height = 6.8,
    dpi = 300,
    bg = "white"
  )

  invisible(p)
}


#' Main Figure 4: sample-size stress test.
plot_figure_4_sample_size_stress_test <- function(
    operating_characteristics,
    output_file
) {
  ensure_directory(dirname(output_file))

  plot_data <- operating_characteristics |>
    dplyr::filter(
      analysis == "Sample-size stress test",
      scenario %in% c(
        "neural_dominant",
        "metabolic_dominant",
        "heterogeneous",
        "null"
      )
    ) |>
    dplyr::select(
      scenario_label,
      n_subjects,
      delta_rmse,
      delta_coverage,
      delta_type1_error,
      delta_power,
      probability_correct_discrimination
    ) |>
    tidyr::pivot_longer(
      cols = c(
        delta_rmse,
        delta_coverage,
        delta_type1_error,
        delta_power,
        probability_correct_discrimination
      ),
      names_to = "metric",
      values_to = "value"
    ) |>
    dplyr::filter(!is.na(value)) |>
    dplyr::mutate(
      metric = dplyr::recode(
        metric,
        delta_rmse = "RMSE of Delta",
        delta_coverage = "95% CI coverage",
        delta_type1_error = "Type I error",
        delta_power = "Power",
        probability_correct_discrimination =
          "Pr(correct discrimination)"
      )
    )

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = n_subjects,
      y = value
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::facet_grid(
      metric ~ scenario_label,
      scales = "free_y"
    ) +
    ggplot2::scale_x_continuous(
      breaks = sort(unique(plot_data$n_subjects))
    ) +
    ggplot2::labs(
      title = "Sample-size stress test for pathway contrast Delta",
      subtitle = "rho_HM = 0.60; mediator reliability = 1.00",
      x = "Simulated cohort size",
      y = NULL
    ) +
    publication_theme(base_size = 12)

  ggplot2::ggsave(
    output_file,
    p,
    width = 14,
    height = 10,
    dpi = 300,
    bg = "white"
  )

  invisible(p)
}


#' Main Figure 5: classical mediator-measurement-error stress test.
plot_figure_5_measurement_error_stress_test <- function(
    operating_characteristics,
    output_file
) {
  ensure_directory(dirname(output_file))

  plot_data <- operating_characteristics |>
    dplyr::filter(
      analysis == "Mediator-measurement-error stress test",
      scenario %in% c(
        "neural_dominant",
        "metabolic_dominant",
        "heterogeneous",
        "null"
      )
    ) |>
    dplyr::mutate(
      measurement_error_percent =
        100 * (1 - mediator_reliability)
    ) |>
    dplyr::select(
      scenario_label,
      measurement_error_percent,
      delta_rmse,
      delta_coverage,
      delta_type1_error,
      delta_power,
      probability_correct_discrimination
    ) |>
    tidyr::pivot_longer(
      cols = c(
        delta_rmse,
        delta_coverage,
        delta_type1_error,
        delta_power,
        probability_correct_discrimination
      ),
      names_to = "metric",
      values_to = "value"
    ) |>
    dplyr::filter(!is.na(value)) |>
    dplyr::mutate(
      metric = dplyr::recode(
        metric,
        delta_rmse = "RMSE of Delta",
        delta_coverage = "95% CI coverage",
        delta_type1_error = "Type I error",
        delta_power = "Power",
        probability_correct_discrimination =
          "Pr(correct discrimination)"
      )
    )

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = measurement_error_percent,
      y = value
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::facet_grid(
      metric ~ scenario_label,
      scales = "free_y"
    ) +
    ggplot2::scale_x_continuous(
      breaks = sort(
        unique(plot_data$measurement_error_percent)
      )
    ) +
    ggplot2::labs(
      title = "Mediator-measurement-error stress test for pathway contrast Delta",
      subtitle = "N = 5000; rho_HM = 0.60",
      x = "Mediator measurement error (%)",
      y = NULL
    ) +
    publication_theme(base_size = 12)

  ggplot2::ggsave(
    output_file,
    p,
    width = 14,
    height = 10,
    dpi = 300,
    bg = "white"
  )

  invisible(p)
}


#' Supplementary Figure S1: SEM fit diagnostics.
plot_figure_s1_sem_fit_diagnostics <- function(
    primary_runs,
    output_file
) {
  ensure_directory(dirname(output_file))

  plot_data <- primary_runs |>
    dplyr::filter(
      analysis == "Primary correlation stress test"
    ) |>
    dplyr::select(
      scenario_label,
      rho_hm,
      cfi,
      tli,
      rmsea,
      srmr
    ) |>
    tidyr::pivot_longer(
      cols = c(cfi, tli, rmsea, srmr),
      names_to = "fit_index",
      values_to = "value"
    )

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = factor(rho_hm),
      y = value
    )
  ) +
    ggplot2::geom_boxplot(
      outlier.alpha = 0.10
    ) +
    ggplot2::facet_grid(
      fit_index ~ scenario_label,
      scales = "free_y"
    ) +
    ggplot2::labs(
      title = "SEM fit diagnostics",
      x = expression(rho[HM]),
      y = NULL
    ) +
    publication_theme(base_size = 11)

  ggplot2::ggsave(
    output_file,
    p,
    width = 13,
    height = 8,
    dpi = 300,
    bg = "white"
  )

  invisible(p)
}


#' Supplementary Figure S2: individual-level 2SLS estimates.
plot_figure_s2_individual_2sls <- function(
    secondary_runs,
    output_file
) {
  ensure_directory(dirname(output_file))

  summary_data <- secondary_runs |>
    dplyr::group_by(
      scenario_label,
      rho_hm
    ) |>
    dplyr::summarise(
      estimate = mean(iv_estimate, na.rm = TRUE),
      se = mean(iv_se, na.rm = TRUE),
      lower = estimate - 1.96 * se,
      upper = estimate + 1.96 * se,
      .groups = "drop"
    )

  p <- ggplot2::ggplot(
    summary_data,
    ggplot2::aes(
      x = estimate,
      y = scenario_label
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = 2
    ) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        xmin = lower,
        xmax = upper
      ),
      orientation = "y",
      width = 0.18
    ) +
    ggplot2::facet_wrap(
      ~ rho_hm,
      ncol = 1
    ) +
    ggplot2::labs(
      title = "Individual-level 2SLS estimates",
      subtitle = "Mean estimate with 95% interval based on mean standard error",
      x = "Estimated signaling-to-cognition effect",
      y = NULL
    ) +
    publication_theme(base_size = 11)

  ggplot2::ggsave(
    output_file,
    p,
    width = 9,
    height = 10,
    dpi = 300,
    bg = "white"
  )

  invisible(p)
}


#' Main Figure 3: GWAS-summary genetic-instrument simulation.
plot_figure_3_gwas_summary <- function(
    gwas_results,
    output_file
) {
  ensure_directory(dirname(output_file))

  gwas_data <- gwas_results$gwas_data
  mr_summary <- gwas_results$mr_summary

  scatter_plot <- ggplot2::ggplot(
    gwas_data,
    ggplot2::aes(
      x = bx,
      y = by,
      color = outcome
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.55,
      size = 1.6
    ) +
    ggplot2::geom_smooth(
      method = "lm",
      formula = y ~ x - 1,
      se = TRUE,
      linewidth = 0.8
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = 2
    ) +
    ggplot2::labs(
      title = "GWAS-summary simulation",
      x = "SNP effect on GLP-1 signaling",
      y = "SNP effect on outcome",
      color = "Outcome"
    ) +
    publication_theme(base_size = 11)

  forest_data <- mr_summary |>
    dplyr::filter(method == "IVW") |>
    dplyr::mutate(
      lower = beta - 1.96 * se,
      upper = beta + 1.96 * se,
      outcome_label = paste0(
        outcome,
        " (N=",
        format(sample_size, big.mark = ",", scientific = FALSE),
        ")"
      )
    )

  forest_plot <- ggplot2::ggplot(
    forest_data,
    ggplot2::aes(
      x = beta,
      y = outcome_label
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = 2
    ) +
    ggplot2::geom_point(size = 2.7) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        xmin = lower,
        xmax = upper
      ),
      orientation = "y",
      width = 0.18
    ) +
    ggplot2::labs(
      title = "IVW estimates",
      x = "Estimated causal effect",
      y = NULL
    ) +
    publication_theme(base_size = 11)

  combined <- scatter_plot / forest_plot +
    patchwork::plot_layout(
      heights = c(1.1, 0.9)
    )

  ggplot2::ggsave(
    output_file,
    combined,
    width = 12,
    height = 11,
    dpi = 300,
    bg = "white"
  )

  invisible(combined)
}


#' Supplementary Figure S3: SEM and 2SLS directional concordance.
plot_figure_s3_sem_2sls_concordance <- function(
    secondary_runs,
    output_file
) {
  ensure_directory(dirname(output_file))

  p <- ggplot2::ggplot(
    secondary_runs,
    ggplot2::aes(
      x = psi_h_hat,
      y = iv_estimate,
      color = scenario_label,
      shape = factor(rho_hm)
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = 2,
      alpha = 0.5
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = 2,
      alpha = 0.5
    ) +
    ggplot2::geom_point(
      alpha = 0.70,
      size = 2.5
    ) +
    ggplot2::geom_smooth(
      method = "lm",
      se = FALSE,
      linewidth = 0.8
    ) +
    ggplot2::labs(
      title = "SEM and 2SLS concordance",
      subtitle = "The two analyses target different causal quantities",
      x = "Neural downstream pathway product (SEM)",
      y = "Signaling-to-cognition estimate (2SLS)",
      color = "Scenario",
      shape = expression(rho[HM])
    ) +
    publication_theme(base_size = 11)

  ggplot2::ggsave(
    output_file,
    p,
    width = 11,
    height = 8,
    dpi = 300,
    bg = "white"
  )

  invisible(p)
}


# ==============================================================================
# 10. Complete analysis workflow
# ==============================================================================

#' Run all analyses needed to reproduce the manuscript results.
run_manuscript_analysis <- function(
    output_dir = "outputs",
    primary_mc_reps = 500,
    secondary_mc_reps = 50,
    primary_seed = 91001,
    secondary_seed = 7001,
    gwas_seed = 4001
) {
  ensure_directory(output_dir)
  data_dir <- file.path(output_dir, "data")
  figure_dir <- file.path(output_dir, "figures")
  ensure_directory(data_dir)
  ensure_directory(figure_dir)

  message("1/5 Building simulation design...")
  design <- build_simulation_design()

  message("2/5 Running primary SEM operating-characteristic simulations...")
  primary_runs <- run_operating_characteristic_simulation(
    design = design$all,
    n_mc = primary_mc_reps,
    seed_base = primary_seed
  )

  operating_characteristics <- summarize_operating_characteristics(
    primary_runs
  )

  message("3/5 Running secondary 2SLS and permutation simulations...")
  secondary_runs <- run_secondary_simulation(
    n_mc = secondary_mc_reps,
    n_subjects = 5000,
    seed_base = secondary_seed
  )

  message("4/5 Running GWAS-summary genetic-instrument simulation...")
  gwas_results <- run_gwas_summary_analysis(
    seed_base = gwas_seed,
    pleiotropy_eta = 0
  )

  readr::write_csv(
    gwas_results$gwas_data,
    file.path(data_dir, "gwas_summary_simulation.csv")
  )
  readr::write_csv(
    gwas_results$mr_summary,
    file.path(data_dir, "gwas_mr_summary.csv")
  )

  tables <- export_manuscript_tables(
    primary_runs = primary_runs,
    operating_characteristics = operating_characteristics,
    secondary_runs = secondary_runs,
    output_dir = output_dir
  )

  message("5/5 Creating manuscript figures...")

  plot_figure_1_causal_architecture(
    file.path(
      figure_dir,
      "figure_1_causal_architecture.png"
    )
  )

  plot_figure_2_pathway_contrast(
    primary_runs,
    file.path(
      figure_dir,
      "figure_2_pathway_contrast.png"
    )
  )

  plot_figure_3_gwas_summary(
    gwas_results,
    file.path(
      figure_dir,
      "figure_3_gwas_summary.png"
    )
  )

  plot_figure_4_sample_size_stress_test(
    operating_characteristics,
    file.path(
      figure_dir,
      "figure_4_sample_size_stress_test.png"
    )
  )

  plot_figure_5_measurement_error_stress_test(
    operating_characteristics,
    file.path(
      figure_dir,
      "figure_5_measurement_error_stress_test.png"
    )
  )

  plot_figure_s1_sem_fit_diagnostics(
    primary_runs,
    file.path(
      figure_dir,
      "figure_s1_sem_fit_diagnostics.png"
    )
  )

  plot_figure_s2_individual_2sls(
    secondary_runs,
    file.path(
      figure_dir,
      "figure_s2_individual_2sls.png"
    )
  )

  plot_figure_s3_sem_2sls_concordance(
    secondary_runs,
    file.path(
      figure_dir,
      "figure_s3_sem_2sls_concordance.png"
    )
  )

  # Record the R environment used for the analysis.
  session_file <- file.path(
    output_dir,
    "session_info.txt"
  )
  capture.output(
    utils::sessionInfo(),
    file = session_file
  )

  message("Analysis complete.")
  message("Outputs written to: ", normalizePath(output_dir))

  invisible(
    list(
      design = design,
      primary_runs = primary_runs,
      operating_characteristics = operating_characteristics,
      secondary_runs = secondary_runs,
      gwas_results = gwas_results,
      tables = tables
    )
  )
}


# ==============================================================================
# 11. Execution settings
# ==============================================================================
#
# These values reproduce the manuscript analysis.
#
# Expected computational burden is substantial because the primary design
# contains 50 design cells and 500 SEM fits per cell. For a short code check,
# set PRIMARY_MC_REPS = 5 and SECONDARY_MC_REPS = 5.
#
# ==============================================================================

OUTPUT_DIR <- "outputs"
PRIMARY_MC_REPS <- 500
SECONDARY_MC_REPS <- 50

analysis_results <- run_manuscript_analysis(
  output_dir = OUTPUT_DIR,
  primary_mc_reps = PRIMARY_MC_REPS,
  secondary_mc_reps = SECONDARY_MC_REPS,
  primary_seed = 91001,
  secondary_seed = 7001,
  gwas_seed = 4001
)
