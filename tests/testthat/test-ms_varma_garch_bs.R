## == == == == == == == == == == == == == == == == == == == == == == ==
## Unit Tests for MS-VARMA-GARCH Functionality
## == == == == == == == == == == == == == == == == == == == == == == == 

## ---- Test Setup ----
## Create minimal data and specifications for testing.


## UNIVARIATE SETUPS ===========================================================

set.seed(123)
y_test <- as.matrix(arima.sim(n = 100, list(ar = 0.5)))
colnames(y_test) <- "series_1"

## Basic univariate spec with Normal distribution
spec_test_uni_norm <- list(
  ## State 1
  list(
    arma_order = c(1, 0),
    garch_model = "garch",
    garch_order = c(1, 1),
    distribution = "norm",
    start_pars = list(
      arma_pars = c(ar1 = 0.1),
      garch_pars = list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
      dist_pars = NULL  ## Normal has no additional parameters
    )
  ),
  ## State 2
  list(
    arma_order = c(1, 0),
    garch_model = "garch",
    garch_order = c(1, 1),
    distribution = "norm",
    start_pars = list(
      arma_pars = c(ar1 = 0.8),
      garch_pars = list(omega = 0.2, alpha1 = 0.2, beta1 = 0.7),
      dist_pars = NULL
    )
  )
)

## Univariate spec with Student-t distribution
spec_test_uni_std <- list(
  ## State 1
  list(
    arma_order = c(1, 0),
    garch_model = "garch",
    garch_order = c(1, 1),
    distribution = "std",
    start_pars = list(
      arma_pars = c(ar1 = 0.1),
      garch_pars = list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
      dist_pars = list(shape = 8.0)  ## Degrees of freedom
    )
  ),
  ## State 2
  list(
    arma_order = c(1, 0),
    garch_model = "garch",
    garch_order = c(1, 1),
    distribution = "std",
    start_pars = list(
      arma_pars = c(ar1 = 0.8),
      garch_pars = list(omega = 0.2, alpha1 = 0.2, beta1 = 0.7),
      dist_pars = list(shape = 10.0)
    )
  )
)

## Backward compatibility alias
spec_test_uni <- spec_test_uni_norm


## MULTIVARIATE SETUPS =========================================================

y_test_mv <- matrix(rnorm(200), ncol = 2)
colnames(y_test_mv) <- c("series_1", "series_2")

## ---- DCC-MVN: Fast Smoke Test Spec ----
## Simple 2-series DCC with MVN distribution (no shape parameter)
## This is the fastest multivariate spec for smoke tests
spec_test_mv_smoke <- list(
  list(
    var_order = 1, 
    garch_spec_fun = "dcc_modelspec",
    distribution = "mvn",
    garch_spec_args = list(
      dcc_order = c(1, 1),
      dynamics = "dcc",
      garch_model = list(univariate = list(
        list(model = "garch", garch_order = c(1, 1), distribution = "norm"),
        list(model = "garch", garch_order = c(1, 1), distribution = "norm")
      ))
    ),
    start_pars = list(
      var_pars = rep(0.1, 6),  ## Intercept + 2 lags for 2 series
      garch_pars = list(
        list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
        list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8)
      ),
      dcc_pars = list(alpha_1 = 0.05, beta_1 = 0.90),
      dist_pars = NULL  ## MVN has no shape parameter
    )
  ),
  list(
    var_order = 1, 
    garch_spec_fun = "dcc_modelspec",
    distribution = "mvn",
    garch_spec_args = list(
      dcc_order = c(1, 1),
      dynamics = "dcc",
      garch_model = list(univariate = list(
        list(model = "garch", garch_order = c(1, 1), distribution = "norm"),
        list(model = "garch", garch_order = c(1, 1), distribution = "norm")
      ))
    ),
    start_pars = list(
      var_pars = rep(0.1, 6),
      garch_pars = list(
        list(omega = 0.2, alpha1 = 0.15, beta1 = 0.75),
        list(omega = 0.2, alpha1 = 0.15, beta1 = 0.75)
      ),
      dcc_pars = list(alpha_1 = 0.1, beta_1 = 0.85),
      dist_pars = NULL
    )
  )
)

## ---- DCC-MVT: Integration Test Spec ----
## DCC with Multivariate Student-t distribution
## More complex due to shape parameter estimation
spec_test_mv_dcc_mvt <- list(
  list(
    var_order = 1,
    garch_spec_fun = "dcc_modelspec",
    distribution = "mvt",
    garch_spec_args = list(
      dcc_order = c(1, 1),
      dynamics = "dcc",
      garch_model = list(univariate = list(
        list(model = "garch", garch_order = c(1, 1), distribution = "norm"),
        list(model = "garch", garch_order = c(1, 1), distribution = "norm")
      ))
    ),
    start_pars = list(
      var_pars = rep(0.1, 6),
      garch_pars = list(
        list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
        list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8)
      ),
      dcc_pars = list(alpha_1 = 0.05, beta_1 = 0.90),
      dist_pars = list(shape = 8.0)  ## Degrees of freedom for MVT
    )
  ),
  list(
    var_order = 1,
    garch_spec_fun = "dcc_modelspec",
    distribution = "mvt",
    garch_spec_args = list(
      dcc_order = c(1, 1),
      dynamics = "dcc",
      garch_model = list(univariate = list(
        list(model = "garch", garch_order = c(1, 1), distribution = "norm"),
        list(model = "garch", garch_order = c(1, 1), distribution = "norm")
      ))
    ),
    start_pars = list(
      var_pars = rep(0.1, 6),
      garch_pars = list(
        list(omega = 0.2, alpha1 = 0.15, beta1 = 0.75),
        list(omega = 0.2, alpha1 = 0.15, beta1 = 0.75)
      ),
      dcc_pars = list(alpha_1 = 0.1, beta_1 = 0.85),
      dist_pars = list(shape = 10.0)
    )
  )
)

## ---- Constant Correlation Spec (DCC with order 0,0) ----
## Useful for faster testing when correlation dynamics aren't important
spec_test_mv_constant_corr <- list(
  list(
    var_order = 1,
    garch_spec_fun = "dcc_modelspec",
    distribution = "mvn",
    garch_spec_args = list(
      dcc_order = c(0, 0),  ## Constant correlation
      garch_model = list(univariate = list(
        list(model = "garch", garch_order = c(1, 1), distribution = "norm"),
        list(model = "garch", garch_order = c(1, 1), distribution = "norm")
      ))
    ),
    start_pars = list(
      var_pars = rep(0.1, 6),
      garch_pars = list(
        list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
        list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8)
      ),
      dcc_pars = NULL,  ## No DCC parameters for constant correlation
      dist_pars = NULL
    )
  ),
  list(
    var_order = 1,
    garch_spec_fun = "dcc_modelspec",
    distribution = "mvn",
    garch_spec_args = list(
      dcc_order = c(0, 0),
      garch_model = list(univariate = list(
        list(model = "garch", garch_order = c(1, 1), distribution = "norm"),
        list(model = "garch", garch_order = c(1, 1), distribution = "norm")
      ))
    ),
    start_pars = list(
      var_pars = rep(0.1, 6),
      garch_pars = list(
        list(omega = 0.2, alpha1 = 0.15, beta1 = 0.75),
        list(omega = 0.2, alpha1 = 0.15, beta1 = 0.75)
      ),
      dcc_pars = NULL,
      dist_pars = NULL
    )
  )
)

## ---- Backward Compatibility Alias ----
## For existing tests that use spec_mv_dcc
spec_uni_garch_dcc <- list(
  model = "garch", 
  garch_order = c(1, 1), 
  distribution = "norm"
)

dcc_spec_args <- list(
  dcc_order = c(1, 1), 
  dynamics = "dcc",
  distribution = "mvn", 
  garch_model = list(univariate = list(spec_uni_garch_dcc, spec_uni_garch_dcc))
)

spec_mv_dcc <- list(
  list(
    var_order = 1, 
    garch_spec_fun = "dcc_modelspec",
    distribution = "mvn",
    garch_spec_args = dcc_spec_args, 
    start_pars = list(
      var_pars = rep(0.1, 6), 
      garch_pars = list(
        list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
        list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8)
      ),
      dcc_pars = list(alpha_1 = 0.05, beta_1 = 0.90),
      dist_pars = NULL
    )
  ),
  list(
    var_order = 1, 
    garch_spec_fun = "dcc_modelspec",
    distribution = "mvn",
    garch_spec_args = dcc_spec_args, 
    start_pars = list(
      var_pars = rep(0.1, 6), 
      garch_pars = list(
        list(omega = 0.2, alpha1 = 0.2, beta1 = 0.7),
        list(omega = 0.2, alpha1 = 0.2, beta1 = 0.7)
      ),
      dcc_pars = list(alpha_1 = 0.1, beta_1 = 0.85),
      dist_pars = NULL
    )
  )
)


## HELPER: Generate Larger Test Data ===========================================

#' Generate test data for more realistic convergence tests
#' @param n Number of observations
#' @param k Number of series (1 for univariate, >1 for multivariate)
#' @param seed Random seed
#' @return Matrix of simulated data
generate_test_data <- function(n = 200, k = 1, seed = 123) {
  set.seed(seed)
  
  if (k == 1) {
    y <- as.matrix(arima.sim(n = n, list(ar = 0.5)))
    colnames(y) <- "series_1"
  } else {
    y <- matrix(rnorm(n * k, sd = 0.5), ncol = k)
    colnames(y) <- paste0("series_", 1:k)
  }
  
  return(y)
}


## == == == == == == == == == == == == == == == == == == == == == == ==
## SUMMARY OF AVAILABLE SPECS
## == == == == == == == == == == == == == == == == == == == == == == ==

## UNIVARIATE:
## - spec_test_uni_norm  : Basic ARMA-GARCH with Normal distribution
## - spec_test_uni_std   : ARMA-GARCH with Student-t distribution
## - spec_test_uni       : Alias for spec_test_uni_norm (backward compatibility)

## MULTIVARIATE:
## - spec_test_mv_smoke          : DCC-MVN (fastest, for smoke tests)
## - spec_test_mv_dcc_mvt        : DCC-MVT (with shape parameter)
## - spec_test_mv_constant_corr  : Constant correlation (no DCC dynamics)
## - spec_mv_dcc                 : Alias for backward compatibility

## TEST DATA:
## - y_test              : Univariate (100 obs)
## - y_test_mv           : Bivariate (100 obs)
## - generate_test_data(): Function to create custom test data

## USAGE EXAMPLES:
## 
## # Univariate smoke test
## fit <- fit_ms_varma_garch(y = y_test, M = 2, spec = spec_test_uni_norm,
##                           control = list(max_iter = 1))
##
## # Multivariate smoke test (fast)
## fit <- fit_ms_varma_garch(y = y_test_mv, M = 2, spec = spec_test_mv_smoke,
##                           model_type = "multivariate", 
##                           control = list(max_iter = 1))
##
## # Multivariate with Student-t
## fit <- fit_ms_varma_garch(y = y_test_mv, M = 2, spec = spec_test_mv_dcc_mvt,
##                           model_type = "multivariate",
##                           control = list(max_iter = 5))
##
## # Generate larger test data
## y_large <- generate_test_data(n = 500, k = 2, seed = 456)


## PART 1: Fast Tests (Always Run) =============================================

context("MS-VARMA-GARCH: Fast Input Validation and Smoke Tests")

test_that("Input validation for fit_ms_varma_garch() works correctly", {
  ## Test 'y' argument
  expect_error(
    fit_ms_varma_garch(y = "not_a_matrix", M = 2, spec = spec_test_uni),
    "Input 'y' must be a numeric matrix or data frame."
  )
  
  y_with_na <- y_test; y_with_na[5] <- NA
  expect_error(
    fit_ms_varma_garch(y = y_with_na, M = 2, spec = spec_test_uni),
    "Input matrix 'y' contains non-finite values"
  )
  
  expect_error(
    fit_ms_varma_garch(
      y = y_test[1, , drop = FALSE], 
      M = 2, 
      spec = spec_test_uni, 
      d = 1
    ),
    "The number of observations must be greater than the differencing order"
  )
  
  ## Test 'M' argument
  expect_error(fit_ms_varma_garch(y = y_test, M = 1, spec = spec_test_uni),
               "'M' must be an integer >= 2.")
  
  ## Test 'spec' argument
  expect_error(fit_ms_varma_garch(y = y_test, M = 2, spec = "not_a_list"),
               "'spec' must be a list of length M.")
  
  expect_error(fit_ms_varma_garch(y = y_test, M = 2, spec = spec_test_uni[1]),
               "'spec' must be a list of length M.")
})

test_that("Input validation for ms_varma_garch_bs() works correctly", {
  ## Test that it catches its own validation errors
  expect_error(
    ms_varma_garch_bs(x = y_test, M = 2, spec = spec_test_uni, num_boots = -1),
    "num_boots must be a positive integer."
  )
  
  expect_error(
    ms_varma_garch_bs(
      x = y_test, 
      M = 2, 
      spec = spec_test_uni, 
      n_boot = NULL, 
      num_blocks = NULL
    ),
    "Must provide a valid value for either n_boot or num_blocks"
  )
  
  ## Test that it propagates errors from the underlying fitter
  expect_error(ms_varma_garch_bs(x = y_test, M = 2, spec = "not_a_list"),
               "'spec' must be a list of length M.")
})


## PART 2: Slow Integration Tests (Skip on CRAN) ===============================

context("MS-VARMA-GARCH: Slow Integration and Convergence Tests")

test_that("Smoke test: Fitter runs for 1 iteration and returns correct structure (univariate)", {
  skip_on_cran()
  
  fit <- fit_ms_varma_garch(y = y_test, M = 2, spec = spec_test_uni,
                            control = list(max_iter = 1))
  
  expect_s3_class(fit, "msm.fit")
  
  expect_named(
    fit, 
    c("model_fits", "P", "log_likelihood", "smoothed_probabilities", "aic", 
      "bic", "d", "y", "call", "convergence", "warnings")
  )
  
  expect_equal(dim(fit$P), c(2, 2))
  
  expect_equal(nrow(fit$smoothed_probabilities), nrow(y_test))
})


test_that("Smoke test: Fitter runs for 1 iteration and returns correct structure (multivariate)", {
  skip_on_cran()
  
  fit <- fit_ms_varma_garch(
    y = y_test_mv, 
    M = 2, 
    spec = spec_test_mv_smoke, 
    model_type = "multivariate",
    control = list(max_iter = 1)
  )
  
  expect_s3_class(fit, "msm.fit")
  
  expect_equal(dim(fit$P), c(2, 2))
  
  expect_equal(nrow(fit$smoothed_probabilities), nrow(y_test_mv))
})


## WARNING: This test may take a long time
test_that("Full estimation converges (univariate)", {
  skip_on_cran()
  
  ## A very short run to check for convergence without taking hours
  y_sim_short <- as.matrix(arima.sim(n = 200, list(ar = 0.5)))
  colnames(y_sim_short) <- "series_1"
  fit <- fit_ms_varma_garch(
    y = y_sim_short, 
    M = 2, 
    spec = spec_test_uni,
    control = list(max_iter = 10, tol = 0.1)
  )
  
  expect_true(is.finite(fit$log_likelihood))
  
  ## A basic sanity check on one parameter
  ar1_est_s1 <- fit$model_fits[[1]]$arma_pars[["ar1"]]
  expect_true(ar1_est_s1 > -1 && ar1_est_s1 < 1)
})


## WARNING: This test may take a long time
test_that("tsbs() with ms_varma_garch runs without error (univariate)", {
  skip_on_cran()
  result <- tsbs(
    x = y_test,
    bs_type = "ms_varma_garch",
    num_boots = 2,
    num_blocks = 10, ## Needs num_blocks for resampling
    ## Arguments for the fitter
    M = 2,
    spec = spec_test_uni,
    control = list(max_iter = 2) ## Keep the fitting part very short
  )
  
  expect_named(result, c("bootstrap_series", "func_outs", "func_out_means"))
  
  expect_length(result$bootstrap_series, 2)
  
  expect_true(is.matrix(result$bootstrap_series[[1]]))
})


## WARNING: This test may take a long time
test_that("Full estimation converges (multivariate)", {
  skip_on_cran()
  
  set.seed(123)
  
  ## Generate bivariate data with GARCH(1,1) effects
  n <- 200
  k <- 2
  
  ## True GARCH parameters
  omega_true <- c(0.1, 0.15)
  alpha_true <- c(0.1, 0.12)
  beta_true <- c(0.8, 0.75)
  
  ## Simulate GARCH series
  y_sim_mv_short <- matrix(0, n, k)
  h <- matrix(0, n, k)
  
  for (i in 1:k) {
    h[1, i] <- omega_true[i] / (1 - alpha_true[i] - beta_true[i])
    y_sim_mv_short[1, i] <- rnorm(1) * sqrt(h[1, i])
    
    for (t in 2:n) {
      h[t, i] <- omega_true[i] + alpha_true[i] * y_sim_mv_short[t-1, i]^2 + 
        beta_true[i] * h[t-1, i]
      y_sim_mv_short[t, i] <- rnorm(1) * sqrt(h[t, i])
    }
  }
  
  colnames(y_sim_mv_short) <- c("s1", "s2")
  
  ## Now fit the model
  suppressWarnings({
    fit <- fit_ms_varma_garch(
      y = y_sim_mv_short, 
      M = 2, 
      spec = spec_mv_dcc, 
      model_type = "multivariate",
      control = list(max_iter = 10, tol = 0.01)
    )
  })
  
  expect_true(is.finite(fit$log_likelihood))
  
  ## Test that alpha1 is positive (since we know the true data has ARCH effects)
  dcc_alpha_s1 <- fit$model_fits[[1]]$garch_pars[[1]]$alpha1
  expect_true(dcc_alpha_s1 > 0 && dcc_alpha_s1 < 1,
              info = paste("alpha1 =", dcc_alpha_s1, "should be strictly positive"))
  
  ## Optionally: Check that at least one state has substantial ARCH effects
  ## (though which state gets which regime is not guaranteed)
  alpha_values <- c(
    fit$model_fits[[1]]$garch_pars[[1]]$alpha1,
    fit$model_fits[[2]]$garch_pars[[1]]$alpha1
  )
  expect_true(any(alpha_values > 0.05),
              info = "At least one state should have alpha1 > 0.05")
})


test_that("tsbs() with ms_varma_garch runs without error (multivariate)", {
  skip_on_cran()
  
  ## A very small run to ensure the full user-facing pipeline connects without
  ## errors. This is an integration test, not a statistical validity test.
  result <- tsbs(
    x = y_test_mv,
    bs_type = "ms_varma_garch",
    num_boots = 2,
    num_blocks = 10,
    ## Arguments for the fitter
    M = 2,
    spec = spec_mv_dcc,
    model_type = "multivariate",
    control = list(max_iter = 2)
  )
  
  expect_named(result, c("bootstrap_series", "func_outs", "func_out_means"))
  
  expect_length(result$bootstrap_series, 2)
  
  expect_true(is.matrix(result$bootstrap_series[[1]]))
  
  expect_equal(ncol(result$bootstrap_series[[1]]), 2)
})


## PART 3: Generalized Univariate Log-Likelihood Calculation ===================

context("Generalized Univariate Log-Likelihood Calculation")

test_that("calculate_loglik_vector_r works for 'norm' distribution", {
  set.seed(123)
  y_test <- arima.sim(model = list(ar = 0.5), n = 100)
  
  spec_norm <- list(
    arma_order = c(1,0),
    garch_model = "garch", garch_order = c(1,1),
    distribution = "norm",
    start_pars = list(dist_pars = NULL)
  )
  current_pars_norm <- list(
    arma_pars = list(ar1 = 0.5),
    garch_pars = list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
    dist_pars = NULL
  )
  
  ll_vec_calculated <- calculate_loglik_vector_r(
    y_test, 
    current_pars_norm, 
    spec_norm, 
    "univariate"
  )
  
  # Ground truth
  residuals <- stats::arima(
    y_test, 
    order = c(1,0,0), 
    fixed = 0.5, 
    include.mean = FALSE
  )$residuals
  garch_spec_obj <- tsgarch::garch_modelspec(
    xts::xts(residuals, Sys.Date()-(length(residuals):1)),
    model="garch", 
    garch_order = c(1,1), 
    distribution = "norm"
  )
  garch_spec_obj$parmatrix[estimate == 1, value := c(0.1, 0.1, 0.8)]
  fit_truth <- tsmethods::tsfilter(garch_spec_obj)
  
  # Use the 'residuals' variable we already have, not fit_truth$residuals
  ll_vec_truth <- dnorm(residuals, 0, fit_truth$sigma, log = TRUE)
  
  T_eff <- length(ll_vec_truth)
  T_orig <- length(ll_vec_calculated)
  expect_equal(
    ll_vec_calculated[(T_orig - T_eff + 1):T_orig], 
    as.numeric(ll_vec_truth), 
    tolerance = 1e-6
  )
})


test_that("calculate_loglik_vector_r works for 'std' distribution", {
  set.seed(123)
  y_test <- arima.sim(model = list(ar = 0.5), n = 100)
  
  spec_std <- list(
    arma_order = c(1,0),
    garch_model = "garch", garch_order = c(1,1),
    distribution = "std",
    start_pars = list(dist_pars = list(shape = 8))
  )
  current_pars_std <- list(
    arma_pars = list(ar1 = 0.5),
    garch_pars = list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
    dist_pars = list(shape = 8.0)
  )
  
  ll_vec_calculated <- calculate_loglik_vector_r(
    y_test, 
    current_pars_std, 
    spec_std, 
    "univariate"
  )
  
  ## Ground truth
  residuals <- stats::arima(
    y_test, 
    order = c(1,0,0), 
    fixed = 0.5, 
    include.mean = FALSE
  )$residuals
  garch_spec_obj <- tsgarch::garch_modelspec(
    xts::xts(
      residuals, 
      Sys.Date()-(length(residuals):1)
    ),
    model="garch", 
    garch_order = c(1,1), 
    distribution = "std"
  )
  garch_spec_obj$parmatrix[estimate == 1, value := c(0.1, 0.1, 0.8, 8.0)]
  fit_truth <- tsmethods::tsfilter(garch_spec_obj)
  
  ## Use the 'residuals' variable we already have, not fit_truth$residuals
  ll_vec_truth <- tsdistributions::dstd(
    residuals, 
    0, 
    fit_truth$sigma, 
    shape = 8.0, 
    log = TRUE
  )
  
  T_eff <- length(ll_vec_truth)
  T_orig <- length(ll_vec_calculated)
  expect_equal(
    ll_vec_calculated[(T_orig - T_eff + 1):T_orig], 
    as.numeric(ll_vec_truth), 
    tolerance = 1e-6
  )
})


## PART 4: Generalize the Univariate Parameter Estimation ======================

context("Generalized Univariate Parameter Estimation")

test_that("estimate_garch_weighted_r recovers known GARCH-std parameters", {
  skip_on_cran() # This test can be a bit slow for CRAN checks
  
  ## 1. Define TRUE parameters and simulate data from that process
  true_pars <- list(omega = 0.2, alpha1 = 0.15, beta1 = 0.75, shape = 7.0)
  
  ## Simulate data
  my_vector <- rnorm(2000)
  dates <- seq(Sys.Date(), by = "day", length.out = length(my_vector))
  my_xts <- xts(my_vector, order.by = dates)
  
  ## Create a spec with the true parameters to simulate from
  sim_spec <- tsgarch::garch_modelspec(
    y = my_xts, ## dummy data
    model = "garch", 
    garch_order = c(1,1),
    distribution = "std"
  )
  ## The := operator in data.table is used for assignment by reference - it 
  ## modifies the data.table in place without making copies, making it very 
  ## memory efficient.
  sim_spec$parmatrix[estimate == 1, value := unlist(true_pars)]
  
  set.seed(42)
  sim_path <- simulate(sim_spec, n.sim = 2000, n.start = 100)
  sim_residuals <- as.numeric(sim_path$y) ## Use the simulated series as our residuals
  
  ## 2. Set up the inputs for our function
  ## Use slightly perturbed starting values to make sure the optimizer is working
  start_pars_test <- list(
    garch_pars = list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
    dist_pars = list(shape = 6.0)
  )
  
  spec_test <- list(
    garch_model = "garch", garch_order = c(1,1),
    distribution = "std",
    start_pars = start_pars_test
  )
  
  ## Use weights of 1, which makes this a standard (unweighted) MLE problem
  test_weights <- rep(1, length(sim_residuals))
  
  ## 3. Call our refactored estimation function
  fit_result <- estimate_garch_weighted_r(
    residuals = sim_residuals,
    weights = test_weights,
    spec = spec_test,
    model_type = "univariate"
  )
  
  ## 4. Check if the estimated parameters are close to the true ones
  ## We use a reasonably large tolerance because MLE on simulated data has noise
  est_coeffs <- fit_result$coefficients
  
  expect_equal(est_coeffs$omega,  true_pars$omega,  tolerance = 0.15)
  expect_equal(est_coeffs$alpha1, true_pars$alpha1, tolerance = 0.15)
  expect_equal(est_coeffs$beta1,  true_pars$beta1,  tolerance = 0.15)
  ## Shape is harder to estimate
  expect_equal(est_coeffs$shape,  true_pars$shape,  tolerance = 1.5) 
})


## PART 5: Test the Univariate Main EM Orchestrator ============================

## Testing an orchestrator function like this is best done with mocking. 
## We will replace the complex estimation functions (estimate_arma_weighted_r 
## and estimate_garch_weighted_r) with simple "mock" versions that return 
## predictable results. This allows us to test only the logic of 
## perform_m_step_parallel_r itself—specifically, whether it correctly 
## separates the parameters.

context("Generalized Univariate M-Step Orchestrator")

test_that("perform_m_step_parallel_r correctly structures the returned parameters", {
  ## 1. Define dummy inputs
  dummy_y <- rnorm(100)
  dummy_weights <- matrix(0.5, nrow = 100, ncol = 2)
  
  ## Define a spec for two states: one normal, one student-t
  spec_test <- list(
    state1 = list( # Normal distribution
      start_pars = list(
        arma_pars = list(ar1 = 0), 
        garch_pars = list(omega = 0.1), 
        dist_pars = NULL
      )
    ),
    state2 = list( # Student-t distribution
      start_pars = list(
        arma_pars = list(ar1 = 0), 
        garch_pars = list(omega = 0.2), 
        dist_pars = list(shape = 5)
      )
    )
  )
  
  ## 2. Define the expected output structure
  expected_output <- list(
    list( ## State 1
      arma_pars = list(ar1 = 0.123), ## From mock_arma_fit
      garch_pars = list(omega = 0.456), ## Separated from mock_garch_fit
      dist_pars = list() ## Separated (empty)
    ),
    list( ## State 2
      arma_pars = list(ar1 = 0.123), ## From mock_arma_fit
      garch_pars = list(omega = 0.789), ## Separated from mock_garch_fit
      dist_pars = list(shape = 9.99) ## Correctly separated
    )
  )
  
  ## 3. Create mock versions of the estimation functions
  ## This mock function will be used for both states
  mock_arma_fit <- function(...) {
    list(coefficients = list(ar1 = 0.123), residuals = rnorm(100))
  }
  
  ## This mock needs to know which state it's being called for
  mock_garch_fit <- function(spec, ...) {
    if (is.null(spec$start_pars$dist_pars)) { ## State 1 (norm)
      return(list(coefficients = list(omega = 0.456)))
    } else { ## State 2 (std)
      return(list(coefficients = list(omega = 0.789, shape = 9.99)))
    }
  }
  
  ## 4. Run the test using testthat::with_mock
  ## This temporarily replaces the real functions with our mocks for this test
  future::plan(future::sequential)
  actual_output <- testthat::with_mocked_bindings(
    estimate_arma_weighted_r = mock_arma_fit,
    estimate_garch_weighted_r = mock_garch_fit,
    {
      perform_m_step_parallel_r(dummy_y, dummy_weights, spec_test, "univariate")
    }
  )
  
  # 5. Compare the actual result to our expected structure
  expect_identical(actual_output, expected_output)
})


## PART 6: Test the Multivariate DCC Log-Likelihood Calculation ================
context("Multivariate DCC Log-Likelihood Calculation")

test_that("calculate_loglik_vector_r produces valid output for DCC-MVN model", {
  skip_on_cran()
  
  ## 1. Setup test data
  set.seed(123)
  y_matrix <- matrix(rnorm(200, sd = 0.5), 100, 2)
  colnames(y_matrix) <- c("series_1", "series_2")
  
  ## Convert to xts
  y_test <- xts::xts(y_matrix, order.by = Sys.Date() - (nrow(y_matrix):1))
  
  ## 2. Specification for DCC-MVN (Normal marginals + MVN copula)
  spec_mv <- list(
    var_order = 1,
    garch_spec_fun = "dcc_modelspec",
    distribution = "mvn",
    garch_spec_args = list(
      dcc_order = c(1, 1),
      garch_model = list(univariate = list(
        list(model = "garch", garch_order = c(1, 1), distribution = "norm"),
        list(model = "garch", garch_order = c(1, 1), distribution = "norm")
      ))
    ),
    start_pars = list(
      dist_pars = NULL
    )
  )
  
  ## 3. Parameter values
  current_pars_mv <- list(
    var_pars = c(0.1, 0.5, 0.1, 0.1, 0.2, 0.4),
    garch_pars = list(
      list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
      list(omega = 0.2, alpha1 = 0.15, beta1 = 0.7)
    ),
    alpha_1 = 0.05,
    beta_1 = 0.90,
    dist_pars = NULL
  )
  
  ## 4. Call our function
  ll_vec_calculated <- calculate_loglik_vector_r(
    y = y_test,
    current_pars = current_pars_mv,
    spec = spec_mv,
    model_type = "multivariate"
  )
  
  ## 5. Verify output properties
  T_obs <- nrow(y_test)
  var_order <- spec_mv$var_order
  
  ## Check basic properties
  expect_true(is.numeric(ll_vec_calculated))
  expect_equal(length(ll_vec_calculated), T_obs)
  expect_true(all(is.finite(ll_vec_calculated)))
  
  ## The first var_order observations should be padded
  expect_true(abs(ll_vec_calculated[1]) < 1e-6 || 
                abs(ll_vec_calculated[1]) < abs(mean(ll_vec_calculated[(var_order+1):T_obs])))
  
  ## Effective observations (after padding) should be reasonable
  effective_ll <- ll_vec_calculated[(var_order+1):T_obs]
  
  ## Should not be too extreme (catching numerical issues)
  expect_true(all(effective_ll > -100))
  expect_true(all(effective_ll < 100))
  
  ## The effective observations should have reasonable values
  effective_ll <- ll_vec_calculated[(var_order+1):T_obs]
  expect_true(sd(effective_ll) > 0.1)
  expect_true(mean(effective_ll) < 10)
  expect_true(mean(effective_ll) > -10)
})


test_that("calculate_loglik_vector_r is consistent across parameter changes", {
  skip_on_cran()
  
  ## Test that changing parameters changes the log-likelihood
  ## Uses SAME parameters as diagnostic script for consistency
  
  set.seed(456)
  y_matrix <- matrix(rnorm(200, sd = 0.5), 100, 2)
  colnames(y_matrix) <- c("series_1", "series_2")
  y_test <- xts::xts(y_matrix, order.by = Sys.Date() - (nrow(y_matrix):1))
  
  spec_mv <- list(
    var_order = 1,
    garch_spec_fun = "dcc_modelspec",
    distribution = "mvn",
    garch_spec_args = list(
      dcc_order = c(1, 1),
      garch_model = list(univariate = list(
        list(model = "garch", garch_order = c(1, 1), distribution = "norm"),
        list(model = "garch", garch_order = c(1, 1), distribution = "norm")
      ))
    ),
    start_pars = list(dist_pars = NULL)
  )
  
  ## Base parameters
  pars_base <- list(
    var_pars = c(0.1, 0.5, 0.1, 0.1, 0.2, 0.4),
    garch_pars = list(
      # list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
      # list(omega = 0.2, alpha1 = 0.15, beta1 = 0.7)
      list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
      list(omega = 0.2, alpha1 = 0.15, beta1 = 0.7)
    ),
    alpha_1 = 0.05,
    beta_1 = 0.90,
    dist_pars = NULL
  )
  
  ## Very different parameters
  pars_different <- list(
    var_pars = c(0.5, 0.8, 0.3, 0.4, 0.5, 0.6),  ## Very different
    garch_pars = list(
      list(omega = 0.5, alpha1 = 0.2, beta1 = 0.7),  ## Much higher
      list(omega = 0.6, alpha1 = 0.3, beta1 = 0.6)   ## Much higher
    ),
    alpha_1 = 0.20,  ## 4x higher (was 0.05)
    beta_1 = 0.70,   ## Much lower (was 0.90)
    dist_pars = NULL
  )
  
  ## Calculate log-likelihoods
  ll_base <- calculate_loglik_vector_r(y_test, pars_base, spec_mv, "multivariate")
  ll_diff <- calculate_loglik_vector_r(y_test, pars_different, spec_mv, "multivariate")
  
  ## All should be valid
  expect_true(all(is.finite(ll_base)))
  expect_true(all(is.finite(ll_diff)))
  
  ## Vectors should NOT be identical
  expect_false(identical(ll_base, ll_diff), 
               info = "Parameter changes should affect log-likelihood")
  
  ## Total log-likelihoods should differ substantially
  sum_base <- sum(ll_base)
  sum_diff <- sum(ll_diff)
  diff_amount <- abs(sum_diff - sum_base)
  
  ## Diagnostic showed difference of ~55, so expect at least 10
  expect_true(diff_amount > 10,
              info = paste0("Total LL should differ substantially. ",
                            "Base: ", round(sum_base, 4), 
                            ", Different: ", round(sum_diff, 4),
                            ", Diff: ", round(diff_amount, 4)))
  
  ## Test individual component changes
  
  ## 1. DCC parameters only (large changes)
  pars_dcc_only <- pars_base
  pars_dcc_only$alpha_1 <- 0.20  ## 4x increase from 0.05
  pars_dcc_only$beta_1 <- 0.70   ## Large decrease from 0.90
  
  ll_dcc <- calculate_loglik_vector_r(y_test, pars_dcc_only, spec_mv, "multivariate")
  
  expect_false(identical(ll_base, ll_dcc),
               info = "DCC parameter changes should affect log-likelihood")
  
  ## 2. GARCH parameters only (large change)
  pars_garch_only <- pars_base
  pars_garch_only$garch_pars[[1]]$omega <- 0.5  ## 5x increase from 0.1
  
  ll_garch <- calculate_loglik_vector_r(
    y = y_test, 
    current_pars = pars_garch_only, 
    spec = spec_mv, 
    model_type = "multivariate"
  )
  
  expect_false(identical(ll_base, ll_garch),
               info = "GARCH parameter changes should affect log-likelihood")
  
  ## 3. VAR parameters only
  pars_var_only <- pars_base
  pars_var_only$var_pars[2] <- 0.8  ## Increase from 0.5
  
  ll_var <- calculate_loglik_vector_r(y_test, pars_var_only, spec_mv, "multivariate")
  
  expect_false(identical(ll_base, ll_var),
               info = "VAR parameter changes should affect log-likelihood")
})


test_that("calculate_loglik_vector_r works with MVT distribution", {
  skip_on_cran()
  
  set.seed(789)
  y_matrix <- matrix(rnorm(200, sd = 0.5), 100, 2)
  colnames(y_matrix) <- c("series_1", "series_2")
  y_test <- xts::xts(y_matrix, order.by = Sys.Date() - (nrow(y_matrix):1))
  
  spec_mvt <- list(
    var_order = 1,
    garch_spec_fun = "dcc_modelspec",
    distribution = "mvt",
    garch_spec_args = list(
      dcc_order = c(1, 1),
      garch_model = list(univariate = list(
        list(model = "garch", garch_order = c(1, 1), distribution = "norm"),
        list(model = "garch", garch_order = c(1, 1), distribution = "norm")
      ))
    ),
    start_pars = list(dist_pars = list(shape = 8))
  )
  
  pars_mvt <- list(
    var_pars = c(0.1, 0.5, 0.1, 0.1, 0.2, 0.4),
    garch_pars = list(
      list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
      list(omega = 0.2, alpha1 = 0.15, beta1 = 0.7)
    ),
    alpha_1 = 0.05,
    beta_1 = 0.90,
    shape = 8.0
  )
  
  ll_vec <- calculate_loglik_vector_r(y_test, pars_mvt, spec_mvt, "multivariate")
  
  ## Verify output
  expect_true(is.numeric(ll_vec))
  expect_true(all(is.finite(ll_vec)))
  expect_equal(length(ll_vec), nrow(y_test))
  
  ## Should not be too extreme
  expect_true(all(ll_vec > -100))
  expect_true(all(ll_vec < 100))
  
  ## Different shape parameter should give different likelihood
  pars_mvt2 <- pars_mvt
  pars_mvt2$shape <- 5.0
  
  ll_vec2 <- calculate_loglik_vector_r(y_test, pars_mvt2, spec_mvt, "multivariate")
  
  expect_false(identical(ll_vec, ll_vec2))
})


test_that("calculate_loglik_vector_r matches tsmarch calculation", {
  skip_on_cran()
  
  ## This test verifies consistency with tsmarch
  
  set.seed(999)
  y_matrix <- matrix(rnorm(200, sd = 0.5), 100, 2)
  colnames(y_matrix) <- c("series_1", "series_2")
  y_test <- xts::xts(y_matrix, order.by = Sys.Date() - (nrow(y_matrix):1))
  
  spec_mv <- list(
    var_order = 1,
    garch_spec_fun = "dcc_modelspec",
    distribution = "mvn",
    garch_spec_args = list(
      dcc_order = c(1, 1),
      dynamics = "dcc",
      garch_model = list(univariate = list(
        list(model = "garch", garch_order = c(1, 1), distribution = "norm"),
        list(model = "garch", garch_order = c(1, 1), distribution = "norm")
      ))
    ),
    start_pars = list(dist_pars = NULL)
  )
  
  ## Use tsmarch to estimate and get parameters
  var_order <- 1
  k <- ncol(y_test)
  T_obs <- nrow(y_test)
  
  ## Compute VAR residuals
  X_lagged <- matrix(1, nrow = T_obs - var_order, ncol = 1 + k * var_order)
  X_lagged[, 2:(1+k)] <- coredata(y_test)[1:(T_obs-1), ]
  y_target <- coredata(y_test)[(var_order+1):T_obs, ]
  var_pars <- rep(0.1, 1 + k * var_order)
  beta_mat <- matrix(var_pars, nrow = 1 + k * var_order, ncol = k)
  residuals_data <- y_target - X_lagged %*% beta_mat
  
  ## Create DCC spec and estimate
  residuals_xts <- xts::xts(residuals_data, order.by = Sys.Date() - (nrow(residuals_data):1))
  
  ## Estimate univariate GARCH models
  uni_models <- lapply(1:k, function(i) {
    suppressWarnings(
      estimate(
        tsgarch::garch_modelspec(
          y = residuals_xts[,i],
          model = "garch",
          garch_order = c(1, 1),
          distribution = "norm"
        ),
        keep_tmb = TRUE
      )
    )
  })
  
  multi_est <- tsgarch::to_multi_estimate(uni_models)
  names(multi_est) <- paste0("series_", 1:k)
  
  ## Create DCC spec
  dcc_spec <- tsmarch::dcc_modelspec(
    object = multi_est,
    dcc_order = c(1, 1),
    dynamics = "dcc",
    distribution = "mvn"
  )
  
  ## Estimate DCC
  dcc_fit <- suppressWarnings(estimate(dcc_spec))
  
  ## Extract parameters
  estimated_pars <- list(
    var_pars = var_pars,
    garch_pars = list(
      list(omega = uni_models[[1]]$parmatrix[parameter == "omega"]$value,
           alpha1 = uni_models[[1]]$parmatrix[parameter == "alpha1"]$value,
           beta1 = uni_models[[1]]$parmatrix[parameter == "beta1"]$value),
      list(omega = uni_models[[2]]$parmatrix[parameter == "omega"]$value,
           alpha1 = uni_models[[2]]$parmatrix[parameter == "alpha1"]$value,
           beta1 = uni_models[[2]]$parmatrix[parameter == "beta1"]$value)
    ),
    alpha_1 = dcc_fit$parmatrix[parameter == "alpha_1"]$value,
    beta_1 = dcc_fit$parmatrix[parameter == "beta_1"]$value,
    dist_pars = NULL
  )
  
  ## Get log-likelihood from our function
  ll_ours <- calculate_loglik_vector_r(
    y = y_test,
    current_pars = estimated_pars,
    spec = spec_mv,
    model_type = "multivariate"
  )
  
  ## Get total log-likelihood from tsmarch
  ll_tsmarch_total <- dcc_fit$loglik
  ll_ours_total <- sum(ll_ours)
  
  ## Determine which sign convention tsmarch uses
  diff_direct <- abs(ll_ours_total - ll_tsmarch_total)
  diff_negated <- abs(ll_ours_total - (-ll_tsmarch_total))
  
  if (diff_direct < diff_negated) {
    ## tsmarch returns log-likelihood (same sign)
    expect_equal(ll_ours_total, ll_tsmarch_total, tolerance = 1e-3,
                 info = paste0("Total log-likelihood should match tsmarch (same sign).\n",
                               "Our LL: ", round(ll_ours_total, 4), "\n",
                               "tsmarch LL: ", round(ll_tsmarch_total, 4)))
  } else {
    ## tsmarch returns negative log-likelihood
    expect_equal(ll_ours_total, -ll_tsmarch_total, tolerance = 1e-3,
                 info = paste0("Total log-likelihood should match negative of tsmarch.\n",
                               "Our LL: ", round(ll_ours_total, 4), "\n",
                               "tsmarch NLL: ", round(ll_tsmarch_total, 4), "\n",
                               "tsmarch -NLL: ", round(-ll_tsmarch_total, 4)))
  }
})


## PART 7: Test Multivariate Parameter Estimation Structure ====================

context("Multivariate DCC Parameter Estimation")

test_that("estimate_garch_weighted_r returns correct structure for DCC", {
  skip_on_cran()
  
  set.seed(789)
  residuals_matrix <- matrix(rnorm(200, sd = 0.5), 100, 2)
  colnames(residuals_matrix) <- c("series_1", "series_2")
  
  weights <- rep(1, 100)
  
  spec_dcc <- list(
    garch_spec_fun = "dcc_modelspec",
    distribution = "mvn",
    garch_spec_args = list(
      dcc_order = c(1, 1),
      garch_model = list(univariate = list(
        list(model = "garch", garch_order = c(1, 1), distribution = "norm"),
        list(model = "garch", garch_order = c(1, 1), distribution = "norm")
      ))
    ),
    start_pars = list(
      garch_pars = list(
        list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
        list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8)
      ),
      dcc_pars = list(alpha_1 = 0.05, beta_1 = 0.90),
      dist_pars = NULL
    )
  )
  
  result <- estimate_garch_weighted_r(
    residuals = residuals_matrix,
    weights = weights,
    spec = spec_dcc,
    model_type = "multivariate"
  )
  
  ## Verify structure
  expect_type(result, "list")
  expect_named(result, c("coefficients", "warnings"))
  
  expect_true("garch_pars" %in% names(result$coefficients))
  expect_true("dcc_pars" %in% names(result$coefficients))
  
  expect_length(result$coefficients$garch_pars, 2)
  
  expect_named(result$coefficients$garch_pars[[1]], c("omega", "alpha1", "beta1"))
  expect_named(result$coefficients$garch_pars[[2]], c("omega", "alpha1", "beta1"))
  
  expect_named(result$coefficients$dcc_pars, c("alpha_1", "beta_1"))
  
  ## Check parameter bounds
  for (i in 1:2) {
    garch_pars <- result$coefficients$garch_pars[[i]]
    expect_true(garch_pars$omega > 0)
    expect_true(garch_pars$alpha1 >= 0)
    expect_true(garch_pars$beta1 >= 0)
    expect_true((garch_pars$alpha1 + garch_pars$beta1) < 1)
  }
  
  dcc_pars <- result$coefficients$dcc_pars
  expect_true(dcc_pars$alpha_1 >= 0)
  expect_true(dcc_pars$beta_1 >= 0)
  expect_true((dcc_pars$alpha_1 + dcc_pars$beta_1) < 1)
})


test_that("estimate_garch_weighted_r handles DCC-MVT with shape parameter", {
  skip_on_cran()
  
  set.seed(101112)
  residuals_matrix <- matrix(rnorm(200, sd = 0.5), 100, 2)
  colnames(residuals_matrix) <- c("series_1", "series_2")
  
  weights <- rep(1, 100)
  
  spec_dcc_t <- list(
    garch_spec_fun = "dcc_modelspec",
    distribution = "mvt",
    garch_spec_args = list(
      dcc_order = c(1, 1),
      garch_model = list(univariate = list(
        list(model = "garch", garch_order = c(1, 1), distribution = "norm"),
        list(model = "garch", garch_order = c(1, 1), distribution = "norm")
      ))
    ),
    start_pars = list(
      garch_pars = list(
        list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
        list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8)
      ),
      dcc_pars = list(alpha_1 = 0.05, beta_1 = 0.90),
      dist_pars = list(shape = 8.0)
    )
  )
  
  result <- estimate_garch_weighted_r(
    residuals = residuals_matrix,
    weights = weights,
    spec = spec_dcc_t,
    model_type = "multivariate"
  )
  
  expect_true("dist_pars" %in% names(result$coefficients))
  expect_named(result$coefficients$dist_pars, "shape")
  expect_true(result$coefficients$dist_pars$shape > 2)
})

test_that("Sigma is computed correctly with fixed parameters in multivariate DCC", {
  skip_on_cran()
  
  set.seed(456)
  
  ## Generate simple bivariate GARCH data
  n <- 150
  k <- 2
  
  omega_true <- c(0.1, 0.12)
  alpha_true <- c(0.08, 0.10)
  beta_true <- c(0.85, 0.80)
  
  y_sim <- matrix(0, n, k)
  h <- matrix(0, n, k)
  
  for (i in 1:k) {
    h[1, i] <- omega_true[i] / (1 - alpha_true[i] - beta_true[i])
    y_sim[1, i] <- rnorm(1) * sqrt(h[1, i])
    
    for (t in 2:n) {
      h[t, i] <- omega_true[i] + alpha_true[i] * y_sim[t-1, i]^2 + 
        beta_true[i] * h[t-1, i]
      y_sim[t, i] <- rnorm(1) * sqrt(h[t, i])
    }
  }
  
  colnames(y_sim) <- c("s1", "s2")
  
  ## Compute VAR(1) residuals
  var_order <- 1
  T_obs <- nrow(y_sim)
  X_lagged <- matrix(1, nrow = T_obs - var_order, ncol = 1 + k * var_order)
  X_lagged[, 2:(1+k)] <- y_sim[1:(T_obs-1), ]
  y_target <- y_sim[(var_order+1):T_obs, ]
  
  var_pars <- c(0, 0, 0)  # Simple intercept-only for testing
  beta_mat <- matrix(var_pars, nrow = 1 + k * var_order, ncol = k)
  residuals_data <- y_target - X_lagged %*% beta_mat
  
  residuals_xts <- xts::xts(residuals_data, 
                            order.by = Sys.Date() - (nrow(residuals_data):1))
  
  ## STEP 1: Estimate univariate GARCH models and extract parameters
  uni_models_estimated <- lapply(1:k, function(i) {
    suppressWarnings(
      estimate(
        tsgarch::garch_modelspec(
          y = residuals_xts[,i],
          model = "garch",
          garch_order = c(1, 1),
          distribution = "norm"
        ),
        keep_tmb = TRUE
      )
    )
  })
  
  ## Extract the estimated parameters
  params_estimated <- list(
    garch_pars = lapply(1:k, function(i) {
      list(
        omega = uni_models_estimated[[i]]$parmatrix[parameter == "omega"]$value,
        alpha1 = uni_models_estimated[[i]]$parmatrix[parameter == "alpha1"]$value,
        beta1 = uni_models_estimated[[i]]$parmatrix[parameter == "beta1"]$value
      )
    }),
    alpha_1 = 0.03,
    beta_1 = 0.92,
    dist_pars = NULL
  )
  
  ## STEP 2: Use create_garch_spec_object_r to build spec with these parameters
  spec_for_test <- list(
    var_order = 1,
    garch_spec_fun = "dcc_modelspec",
    distribution = "mvn",
    garch_spec_args = list(
      dcc_order = c(1, 1),
      dynamics = "dcc",
      garch_model = list(univariate = list(
        list(model = "garch", garch_order = c(1, 1), distribution = "norm"),
        list(model = "garch", garch_order = c(1, 1), distribution = "norm")
      ))
    )
  )
  
  ## Create spec object with our parameters
  garch_spec_obj <- create_garch_spec_object_r(
    residuals = residuals_data,
    spec = spec_for_test,
    model_type = "multivariate",
    current_pars = params_estimated
  )
  
  ## STEP 3: Verify that sigma in the univariate models matches the parameters
  for (i in 1:k) {
    uni_model <- garch_spec_obj$univariate[[paste0("series_", i)]]
    
    ## Check parameters are set correctly
    omega_in_spec <- uni_model$parmatrix[parameter == "omega"]$value
    alpha1_in_spec <- uni_model$parmatrix[parameter == "alpha1"]$value
    beta1_in_spec <- uni_model$parmatrix[parameter == "beta1"]$value
    
    expect_equal(omega_in_spec, params_estimated$garch_pars[[i]]$omega,
                 tolerance = 1e-6,
                 info = paste("Series", i, "omega should match"))
    
    expect_equal(alpha1_in_spec, params_estimated$garch_pars[[i]]$alpha1,
                 tolerance = 1e-6,
                 info = paste("Series", i, "alpha1 should match"))
    
    expect_equal(beta1_in_spec, params_estimated$garch_pars[[i]]$beta1,
                 tolerance = 1e-6,
                 info = paste("Series", i, "beta1 should match"))
    
    ## CRITICAL TEST: Recompute sigma manually with TMB and compare
    if (!is.null(uni_model$tmb)) {
      estimate_col <- NULL
      pars_for_tmb <- uni_model$parmatrix[estimate_col == 1, "value", drop = TRUE]
      
      ## Get sigma from TMB report
      tmb_report <- uni_model$tmb$report(pars_for_tmb)
      
      ## Strip initialization period
      maxpq <- max(uni_model$spec$model$order)
      sigma_from_tmb <- if (maxpq > 0) {
        tmb_report$sigma[-(1:maxpq)]
      } else {
        tmb_report$sigma
      }
      
      ## Get sigma stored in the object
      sigma_in_object <- as.numeric(sigma(uni_model))
      
      ## These should be IDENTICAL (or very close)
      expect_equal(sigma_from_tmb, sigma_in_object,
                   tolerance = 1e-10,
                   info = paste("Series", i, 
                                "sigma from TMB should match stored sigma"))
    }
  }
  
  ## STEP 4: Verify DCC parameters are set correctly
  dcc_alpha <- garch_spec_obj$parmatrix[parameter == "alpha_1"]$value
  dcc_beta <- garch_spec_obj$parmatrix[parameter == "beta_1"]$value
  
  expect_equal(dcc_alpha, params_estimated$alpha_1, tolerance = 1e-6)
  expect_equal(dcc_beta, params_estimated$beta_1, tolerance = 1e-6)
})