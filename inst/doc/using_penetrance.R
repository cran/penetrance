## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(penetrance)

## ----dat, include=FALSE-------------------------------------------------------
dat <- simulated_families

## -----------------------------------------------------------------------------
# The default prior 
prior_params_default <- list(
    asymptote = list(g1 = 1, g2 = 1),
    threshold = list(min = 15, max = 35),
    median = list(m1 = 2, m2 = 2),
    first_quartile = list(q1 = 6, q2 = 3)
)

distribution_data_test <- data.frame(
  row.names = c("min", "first_quartile", "median", "max"),
  age = c(NA, NA, NA, NA),
  at_risk = c(NA, NA, NA, NA)
)

## -----------------------------------------------------------------------------
#' @param twins A list specifying identical twins or triplets in the family. For example, to indicate that "ora024" and "ora027" are identical twins, and "aey063" and "aey064" are identical twins, use the following format: `twins <- list(c("ora024", "ora027"), c("aey063", "aey064"))`.
#' @param n_chains Integer, the number of chains for parallel computation. Default is 1.
#' @param n_iter_per_chain Integer, the number of iterations for each chain. Default is 10000.
#' @param ncores Integer, the number of cores for parallel computation. Default is 6.
#' @param baseline_data Data for the baseline risk estimates (probability of developing cancer), such as population-level risk from a cancer registry. Default data, for exemplary purposes, is for Colorectal cancer from the SEER database.
#' @param max_age Integer, the maximum age considered for analysis. Default is 94.
#' @param remove_proband Logical, indicating whether to remove probands from the analysis. Default is FALSE.
#' @param age_imputation Logical, indicating whether to perform age imputation. Default is FALSE.
#' @param median_max Logical, indicating whether to use the baseline median age or `max_age` as an upper bound for the median proposal. Default is TRUE.
#' @param BaselineNC Logical, indicating that the non-carrier penetrance is assumed to be the baseline penetrance. Default is TRUE.
#' @param var Numeric vector, variances for the proposal distribution in the Metropolis-Hastings algorithm. Default is `c(0.1, 0.1, 2, 2, 5, 5, 5, 5)`.
#' @param burn_in Numeric, the fraction of results to discard as burn-in (0 to 1). Default is 0 (no burn-in).
#' @param thinning_factor Integer, the factor by which to thin the results. Default is 1 (no thinning).
#' @param imp_interval Integer, the interval at which age imputation should be performed when age_imputation = TRUE.
#' @param distribution_data Data for generating prior distributions.
#' @param prev Numeric, prevalence of the carrier status. Default is 0.0001.
#' @param sample_size Optional numeric, sample size for distribution generation.
#' @param ratio Optional numeric, ratio parameter for distribution generation.
#' @param prior_params List, parameters for prior distributions.
#' @param risk_proportion Numeric, proportion of risk for distribution generation.
#' @param summary_stats Logical, indicating whether to include summary statistics in the output. Default is TRUE.
#' @param rejection_rates Logical, indicating whether to include rejection rates in the output. Default is TRUE.
#' @param density_plots Logical, indicating whether to include density plots in the output. Default is TRUE.
#' @param plot_trace Logical, indicating whether to include trace plots in the output. Default is TRUE.
#' @param penetrance_plot Logical, indicating whether to include penetrance plots in the output. Default is TRUE.
#' @param penetrance_plot_pdf Logical, indicating whether to include PDF plots in the output. Default is TRUE.
#' @param plot_loglikelihood Logical, indicating whether to include log-likelihood plots in the output. Default is TRUE.
#' @param plot_acf Logical, indicating whether to include autocorrelation function (ACF) plots for posterior samples. Default is TRUE.
#' @param probCI Numeric, probability level for credible intervals in penetrance plots. Must be between 0 and 1. Default is 0.95.
#' @param sex_specific Logical, indicating whether to use sex-specific parameters in the analysis. Default is TRUE.

## -----------------------------------------------------------------------------
# This is exemplary SEER penetrance data for Colorectal Cancer. 
baseline_data_default <- data.frame(
  Age = 1:94,
  Female = c(
    2.8e-07, 9e-08, 1e-08, 3e-08, 5e-08, 7e-08, 9e-08, 7e-08, 5e-08, 4e-08, 2e-08, 4e-08, 3.2e-07, 6.3e-07, 9.5e-07, 1.27e-06,
    1.94e-06, 4.73e-06, 7.88e-06, 1.102e-05, 1.417e-05, 1.916e-05, 3.517e-05, 5.303e-05, 7.088e-05, 8.873e-05, 0.00010961,
    1.486e-04, 1.906e-04, 0.00023261, 0.00027461, 0.00032101, 0.00039385, 0.00047109, 0.00054832, 0.00062554, 7.166e-04,
    0.00089081, 0.00107885, 0.00126684, 0.00145479, 0.00163827, 0.00179526, 0.00194779, 0.00210026, 0.00225264, 0.00239638,
    0.00248859, 0.00257215, 0.00265562, 0.00273901, 0.00281841, 0.00287439, 0.00292639, 0.00297831, 0.00303014, 0.00309486,
    0.00323735, 0.00339265, 0.00354775, 0.00370266, 0.00385852, 0.00402115, 0.0041847, 0.00434799, 4.511e-03, 0.00466219,
    0.00474384, 0.00481372, 0.00488337, 0.0049528, 0.00500362, 0.00494413, 0.00486620, 0.00478821, 0.00471016, 0.00462597,
    0.00450511, 0.00437818, 0.00425130, 0.00412450, 0.00399836, 0.00387577, 0.00375380, 0.00363187, 0.00351002, 0.00338186,
    0.00321541, 0.00304280, 0.00287050, 0.00269853, 0.00253244, 0.00239963, 0.00227262
  ),
  Male = c(
    2.8e-07, 9e-08, 1e-08, 3e-08, 5e-08, 7e-08, 9e-08, 7e-08, 5e-08, 4e-08, 2e-08, 4e-08, 3.2e-07, 6.3e-07, 9.5e-07, 1.27e-06,
    1.94e-06, 4.73e-06, 7.88e-06, 1.102e-05, 1.417e-05, 1.916e-05, 3.517e-05, 5.303e-05, 7.088e-05, 8.873e-05, 0.00010961,
    1.486e-04, 1.906e-04, 0.00023261, 0.00027461, 0.00032101, 0.00039385, 0.00047109, 0.00054832, 0.00062554, 7.166e-04,
    0.00089081, 0.00107885, 0.00126684, 0.00145479, 0.00163827, 0.00179526, 0.00194779, 0.00210026, 0.00225264, 0.00239638,
    0.00248859, 0.00257215, 0.00265562, 0.00273901, 0.00281841, 0.00287439, 0.00292639, 0.00297831, 0.00303014, 0.00309486,
    0.00323735, 0.00339265, 0.00354775, 0.00370266, 0.00385852, 0.00402115, 0.0041847, 0.00434799, 4.511e-03, 0.00466219,
    0.00474384, 0.00481372, 0.00488337, 0.0049528, 0.00500362, 0.00494413, 0.00486620, 0.00478821, 0.00471016, 0.00462597,
    0.00450511, 0.00437818, 0.00425130, 0.00412450, 0.00399836, 0.00387577, 0.00375380, 0.00363187, 0.00351002, 0.00338186,
    0.00321541, 0.00304280, 0.00287050, 0.00269853, 0.00253244, 0.00239963, 0.00227262
  )
)

## ----eval=FALSE---------------------------------------------------------------
# output <- penetrance(
#     pedigree  = dat, twins = NULL, n_chains = 1, n_iter_per_chain = 20000,
#    ncores = 2, baseline_data = baseline_data_default, prev  = 0.0001,
#     prior_params = prior_params_default, burn_in = 0.1, median_max = TRUE,
#     age_imputation = TRUE, sex_specific = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# plot_penetrance(data = output$combined_chains, prob = 0.95)
# plot_pdf(data = output$combined_chains, prob = 0.95)

