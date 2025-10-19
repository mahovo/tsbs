## Test file for compute_loglik_fixed function

## Helper to get fixtures directory
.get_fixtures_dir <- function() {
  if (file.exists("tests/testthat/fixtures")) {
    return("tests/testthat/fixtures")
  } else if (file.exists("fixtures")) {
    return("fixtures")
  } else if (exists("test_path")) {
    return(test_path("fixtures"))
  } else {
    stop("Cannot find fixtures directory")
  }
}

## Setup: Create mock objects for testing
## Note: These tests will need actual tsmarch objects for integration testing

test_that("compute_loglik_fixed validates input object class", {
  
  ## Test with invalid object class
  expect_error(
    compute_loglik_fixed(
      object = list(some = "data"),
      params = list(alpha_1 = 0.05)
    ),
    "must be an estimated tsmarch model"
  )
  
  ## Test with NULL object
  expect_error(
    compute_loglik_fixed(
      object = NULL,
      params = list(alpha_1 = 0.05)
    ),
    "must be an estimated tsmarch model"
  )
})

test_that("compute_loglik_fixed validates params argument", {
  skip_if_not_installed("tsmarch")
  skip_on_cran()
  
  ## Load a fixture
  dcc_fit <- readRDS(file.path(.get_fixtures_dir(), "dcc_dynamic_fit.rds"))
  
  ## Test with non-list params
  expect_error(
    compute_loglik_fixed(
      object = dcc_fit,
      params = "not a list"
    ),
    "params must be a named list"
  )
  
  ## Test with numeric vector instead of list
  expect_error(
    compute_loglik_fixed(
      object = dcc_fit,
      params = c(0.05, 0.90)
    ),
    "params must be a named list"
  )
})

test_that(".update_parmatrix correctly updates parameter values", {
  skip_if_not_installed("data.table")
  
  ## Create a simple parmatrix
  parmatrix <- data.table::data.table(
    parameter = c("alpha_1", "beta_1", "shape"),
    value = c(0.05, 0.90, 5.0),
    lower = c(0, 0, 2.01),
    upper = c(1, 1, 50),
    estimate = c(1, 1, 1)
  )
  
  ## Update parameters
  new_params <- list(alpha_1 = 0.03, beta_1 = 0.95)
  updated <- .update_parmatrix(parmatrix, new_params)
  
  ## Check updates
  expect_equal(updated[parameter == "alpha_1"]$value, 0.03)
  expect_equal(updated[parameter == "beta_1"]$value, 0.95)
  expect_equal(updated[parameter == "shape"]$value, 5.0)  ## Unchanged
})

test_that(".update_parmatrix warns about unknown parameters", {
  skip_if_not_installed("data.table")
  
  parmatrix <- data.table::data.table(
    parameter = c("alpha_1", "beta_1"),
    value = c(0.05, 0.90),
    estimate = c(1, 1)
  )
  
  ## Try to update non-existent parameter
  expect_warning(
    .update_parmatrix(parmatrix, list(gamma_1 = 0.02)),
    "Parameter 'gamma_1' not found"
  )
})

test_that(".update_parmatrix handles multiple matching parameters", {
  skip_if_not_installed("data.table")
  
  ## Edge case: duplicate parameter names (shouldn't happen in practice)
  parmatrix <- data.table::data.table(
    parameter = c("alpha_1", "alpha_1", "beta_1"),
    value = c(0.05, 0.05, 0.90),
    estimate = c(1, 1, 1)
  )
  
  expect_warning(
    .update_parmatrix(parmatrix, list(alpha_1 = 0.03)),
    "Multiple parameters named 'alpha_1'"
  )
})

## == == == == == == == == == == == == == == == == == == == == == == == == ==
## Integration tests with actual tsmarch objects
## == == == == == == == == == == == == == == == == == == == == == == == == ==

## Test DCC Dynamic Model ====================================================

test_that("compute_loglik_fixed works with DCC dynamic (mvn)", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  skip_if_not_installed("xts")
  skip_on_cran()
  
  #3 Get or create test fixture
  dcc_fit <- readRDS(file.path(.get_fixtures_dir(), "dcc_dynamic_fit.rds"))
  
  ## Get estimated parameters
  est_params <- coef(dcc_fit)
  
  ## Convert to named list
  params_list <- as.list(est_params)
  names(params_list) <- names(est_params)
  
  ## Compute log-likelihood at estimated parameters
  ll_fixed <- compute_loglik_fixed(dcc_fit, params = params_list)
  
  ## Get the estimated log-likelihood
  series1_nll <- dcc_fit$spec$univariate$series1$loglik
  series2_nll <- dcc_fit$spec$univariate$series2$loglik
  mv_nll <- dcc_fit$loglik
  nll_estimated <- series1_nll + series2_nll + mv_nll
  ll_estimated <- -nll_estimated
  
  ## Should match within numerical tolerance
  expect_equal(ll_fixed, ll_estimated, tolerance = 1e-4,
               info = "Log-likelihood at estimated params should match estimation")
})

test_that("compute_loglik_fixed detects parameter changes in DCC dynamic", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  skip_if_not_installed("xts")
  skip_on_cran()
  
  dcc_fit <- readRDS(file.path(.get_fixtures_dir(), "dcc_dynamic_fit.rds"))
  
  ## Get estimated parameters
  est_params <- coef(dcc_fit)
  
  ## Compute at estimated parameters
  ll_at_est <- compute_loglik_fixed(dcc_fit, 
                                    params = list(alpha_1 = est_params["alpha_1"],
                                                  beta_1 = est_params["beta_1"]))
  
  ## Compute at different parameters (should be lower)
  ll_at_other <- compute_loglik_fixed(dcc_fit,
                                      params = list(alpha_1 = 0.03,
                                                    beta_1 = 0.95))
  
  ## Likelihood at estimated parameters should be higher
  expect_true(ll_at_est >= ll_at_other,
              info = "Likelihood should be maximized at estimated parameters")
})

test_that("compute_loglik_fixed return_components works for DCC dynamic", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  skip_if_not_installed("xts")
  skip_on_cran()
  
  dcc_fit <- readRDS(file.path(.get_fixtures_dir(), "dcc_dynamic_fit.rds"))
  est_params <- coef(dcc_fit)
  
  ## Get components
  ll_comp <- compute_loglik_fixed(
    dcc_fit, 
    params = as.list(est_params),
    return_components = TRUE
  )
  
  ## Check structure
  expect_type(ll_comp, "list")
  expect_named(ll_comp, c("loglik", "garch_loglik", "multivariate_loglik"))
  
  ## Check that components sum to total
  sum_comp <- ll_comp$garch_loglik + ll_comp$multivariate_loglik
  expect_equal(ll_comp$loglik, sum_comp, tolerance = 1e-10)
  
  ## Check against simple call
  ll_simple <- compute_loglik_fixed(dcc_fit, params = as.list(est_params))
  expect_equal(ll_comp$loglik, ll_simple, tolerance = 1e-10)
})

## Test DCC Constant Model ==================================================

test_that("compute_loglik_fixed works with DCC constant (mvn)", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  skip_if_not_installed("xts")
  skip_on_cran()
  
  dcc_fit <- readRDS(file.path(.get_fixtures_dir(), "dcc_constant_fit.rds"))
  
  ## For constant MVN, there are no parameters to estimate
  ## Should still compute likelihood correctly
  ll_fixed <- compute_loglik_fixed(dcc_fit, params = list())
  ll_estimated <- as.numeric(logLik(dcc_fit))
  
  expect_equal(ll_fixed, ll_estimated, tolerance = 1e-4)
})

## Test Student-t DCC Model =================================================

test_that("compute_loglik_fixed works with DCC Student-t", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  skip_if_not_installed("xts")
  skip_on_cran()
  
  dcc_fit <- readRDS(file.path(.get_fixtures_dir(), "dcc_student_fit.rds"))
  
  ## Get estimated parameters (should include shape parameter)
  est_params <- coef(dcc_fit)
  
  ## Should have alpha, beta, and shape
  expect_true("shape" %in% names(est_params))
  
  ## Compute log-likelihood
  ll_fixed <- compute_loglik_fixed(dcc_fit, params = as.list(est_params))
  series1_nll <- dcc_fit$spec$univariate$series1$loglik
  series2_nll <- dcc_fit$spec$univariate$series2$loglik
  mv_nll <- dcc_fit$loglik
  nll_estimated <- series1_nll + series2_nll + mv_nll
  ll_estimated <- -nll_estimated
  
  expect_equal(ll_fixed, ll_estimated, tolerance = 1e-4)
})

test_that("compute_loglik_fixed handles different shape values", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  skip_if_not_installed("xts")
  skip_on_cran()
  
  dcc_fit <- readRDS(file.path(.get_fixtures_dir(), "dcc_student_fit.rds"))
  est_params <- coef(dcc_fit)
  
  ## Try different degrees of freedom
  shapes <- c(3, 5, 10, 30)
  lls <- sapply(shapes, function(nu) {
    compute_loglik_fixed(dcc_fit, 
                         params = list(
                           alpha_1 = est_params["alpha_1"],
                           beta_1 = est_params["beta_1"],
                           shape = nu
                         ))
  })
  
  ## Should all be finite
  expect_true(all(is.finite(lls)))
  
  ## Likelihood at estimated shape should be highest
  est_shape_idx <- which.min(abs(shapes - est_params["shape"]))
  if (est_shape_idx > 1 && est_shape_idx < length(shapes)) {
    ## If estimated shape is in the interior of our grid, 
    ## it should have higher likelihood than neighbors
    expect_true(lls[est_shape_idx] >= max(lls[est_shape_idx - 1], 
                                          lls[est_shape_idx + 1]))
  }
})

## Test Copula-GARCH Model ==================================================

test_that("compute_loglik_fixed works with Copula-GARCH dynamic", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  skip_if_not_installed("xts")
  skip_on_cran()
  
  cgarch_fit <- readRDS(file.path(.get_fixtures_dir(), "copula_dynamic_fit.rds"))
  
  ## Get estimated parameters
  est_params <- coef(cgarch_fit)
  
  ## Compute log-likelihood
  ll_fixed <- compute_loglik_fixed(cgarch_fit, params = as.list(est_params))
  series1_nll <- cgarch_fit$spec$univariate$series1$loglik
  series2_nll <- cgarch_fit$spec$univariate$series2$loglik
  copula_nll <- cgarch_fit$copula_nll
  model_nll <- series1_nll + series2_nll + copula_nll
  ll_expected <- -model_nll
  
  expect_equal(ll_fixed, ll_expected, tolerance = 1e-4)
})

## Profile Likelihood Tests =================================================

test_that("profile likelihood is smooth around optimum", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  skip_if_not_installed("xts")
  skip_on_cran()
  
  dcc_fit <- readRDS(file.path(.get_fixtures_dir(), "dcc_dynamic_fit.rds"))
  est_params <- coef(dcc_fit)
  
  ## Create grid around estimated alpha
  alpha_est <- est_params["alpha_1"]
  beta_est <- est_params["beta_1"]
  alpha_grid <- seq(alpha_est - 0.02, alpha_est + 0.02, length.out = 5)
  
  ## Compute profile likelihood
  profile_ll <- sapply(alpha_grid, function(a) {
    compute_loglik_fixed(dcc_fit, 
                         params = list(alpha_1 = a, beta_1 = beta_est))
  })
  
  ## Should all be finite
  expect_true(all(is.finite(profile_ll)),
              info = "All profile likelihoods should be finite")
  
  ## The likelihood at the estimated parameter should be among the highest
  ## (allowing for numerical tolerance and that we're fixing beta)
  max_ll <- max(profile_ll)
  ll_at_est <- profile_ll[3]  # Middle point is the estimated value
  
  ## Check that estimated is within a reasonable range of maximum
  expect_true(ll_at_est >= max_ll - 1,
              info = paste("LL at estimated alpha should be near maximum.",
                           "Estimated:", ll_at_est, "Max:", max_ll))
})

## Likelihood Ratio Test Example ============================================

test_that("likelihood ratio test framework works", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  skip_if_not_installed("xts")
  skip_on_cran()
  
  dcc_fit <- readRDS(file.path(.get_fixtures_dir(), "dcc_dynamic_fit.rds"))
  
  ## H0: alpha = 0.05, beta = 0.90 (restricted)
  ll_restricted <- compute_loglik_fixed(dcc_fit,
                                        params = list(alpha_1 = 0.05, beta_1 = 0.90))
  
  ## H1: free parameters (unrestricted)
  ll_unrestricted <- as.numeric(logLik(dcc_fit))
  
  ## LR statistic
  lr_stat <- 2 * (ll_unrestricted - ll_restricted)
  
  ## Should be non-negative
  expect_true(lr_stat >= 0,
              info = "LR statistic should be non-negative")
  
  ## Should be finite
  expect_true(is.finite(lr_stat))
})

## Error Handling Tests =====================================================

test_that("compute_loglik_fixed handles invalid parameters gracefully", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  skip_if_not_installed("xts")
  skip_on_cran()
  
  dcc_fit <- readRDS(file.path(.get_fixtures_dir(), "dcc_dynamic_fit.rds"))
  
  ## Try with unknown parameter name
  expect_warning(
    compute_loglik_fixed(dcc_fit, params = list(unknown_param = 0.5)),
    "not found in model parmatrix"
  )
})

test_that("compute_loglik_fixed preserves original object", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  skip_if_not_installed("xts")
  skip_on_cran()
  
  dcc_fit <- readRDS(file.path(.get_fixtures_dir(), "dcc_dynamic_fit.rds"))
  
  ## Store original parameters
  original_params <- coef(dcc_fit)
  
  ## Compute with different parameters
  ll <- compute_loglik_fixed(dcc_fit, params = list(alpha_1 = 0.03, beta_1 = 0.95))
  
  ## Check original not modified
  new_params <- coef(dcc_fit)
  expect_equal(original_params, new_params)
})

## == == == == == == == == == == == == == == == == == == == == == == == == ==
## Tests with ll_vec = TRUE                                             =====
## == == == == == == == == == == == == == == == == == == == == == == == == ==

test_that("compute_loglik_fixed returns correct dimensions with ll_vec = TRUE", {
  skip_on_cran()
  
  ## Suppress warnings from GARCH estimation (known issue with small samples)
  suppressWarnings({
    ## Generate small sample data
    set.seed(123)
    n <- 100
    returns <- matrix(rnorm(n * 2), ncol = 2)
    returns <- xts(returns, order.by = seq.Date(Sys.Date() - n + 1, 
                                                Sys.Date(), by = "day"))
    colnames(returns) <- c("series1", "series2")
    
    ## Estimate univariate GARCH models
    spec1 <- garch_modelspec(returns[,1], model = "garch", order = c(1,1))
    spec2 <- garch_modelspec(returns[,2], model = "garch", order = c(1,1))
    fit1 <- estimate(spec1, keep_tmb = TRUE)
    fit2 <- estimate(spec2, keep_tmb = TRUE)
    
    ## Combine into multivariate
    garch_fits <- to_multi_estimate(list(fit1, fit2))
    
    ## Estimate DCC model
    dcc_spec <- dcc_modelspec(garch_fits, dynamics = "dcc", dcc_order = c(1,1))
    dcc_fit <- estimate(dcc_spec)
    
    ## Get estimated parameters
    est_params <- coef(dcc_fit)
    
    ## Test ll_vec = TRUE
    ll_vec_result <- compute_loglik_fixed(
      dcc_fit, 
      params = as.list(est_params), 
      ll_vec = TRUE
    )
  })
  
  ## Should return a vector
  expect_true(is.numeric(ll_vec_result))
  expect_true(is.vector(ll_vec_result))
  expect_equal(length(ll_vec_result), n)
  expect_false(any(is.na(ll_vec_result)))
})

test_that("compute_loglik_fixed ll_vec sums to total loglik for DCC dynamic", {
  skip_on_cran()
  
  suppressWarnings({
    ## Generate small sample data
    set.seed(456)
    n <- 100
    returns <- matrix(rnorm(n * 2), ncol = 2)
    returns <- xts(returns, order.by = seq.Date(Sys.Date() - n + 1, 
                                                Sys.Date(), by = "day"))
    colnames(returns) <- c("asset1", "asset2")
    
    ## Estimate univariate GARCH models
    spec1 <- garch_modelspec(returns[,1], model = "garch", order = c(1,1))
    spec2 <- garch_modelspec(returns[,2], model = "garch", order = c(1,1))
    fit1 <- estimate(spec1, keep_tmb = TRUE)
    fit2 <- estimate(spec2, keep_tmb = TRUE)
    
    ## Combine into multivariate
    garch_fits <- to_multi_estimate(list(fit1, fit2))
    
    ## Estimate DCC model
    dcc_spec <- dcc_modelspec(garch_fits, dynamics = "dcc", dcc_order = c(1,1))
    dcc_fit <- estimate(dcc_spec)
    
    ## Get estimated parameters
    est_params <- coef(dcc_fit)
    
    ## Get per-observation log-likelihoods
    ll_vec_result <- compute_loglik_fixed(dcc_fit, params = as.list(est_params), 
                                          ll_vec = TRUE)
    
    ## Get total log-likelihood
    ll_total <- compute_loglik_fixed(dcc_fit, params = as.list(est_params))
  })
  
  ## Sum of per-observation should equal total
  expect_equal(sum(ll_vec_result), ll_total, tolerance = 1e-6)
})

test_that("compute_loglik_fixed ll_vec works for DCC constant", {
  skip_on_cran()
  
  suppressWarnings({
    ## Generate small sample data
    set.seed(789)
    n <- 100
    returns <- matrix(rnorm(n * 2), ncol = 2)
    returns <- xts(returns, order.by = seq.Date(Sys.Date() - n + 1, 
                                                Sys.Date(), by = "day"))
    colnames(returns) <- c("asset1", "asset2")
    
    ## Estimate univariate GARCH models
    spec1 <- garch_modelspec(returns[,1], model = "garch", order = c(1,1))
    spec2 <- garch_modelspec(returns[,2], model = "garch", order = c(1,1))
    fit1 <- estimate(spec1, keep_tmb = TRUE)
    fit2 <- estimate(spec2, keep_tmb = TRUE)
    
    ## Combine into multivariate
    garch_fits <- to_multi_estimate(list(fit1, fit2))
    
    ## Estimate constant correlation DCC model
    dcc_spec <- dcc_modelspec(garch_fits, dynamics = "constant")
    dcc_fit <- estimate(dcc_spec)
    
    ## Get estimated parameters (may be empty for mvn constant)
    est_params <- coef(dcc_fit)
    
    ## Get per-observation log-likelihoods
    ll_vec_result <- compute_loglik_fixed(dcc_fit, params = as.list(est_params), 
                                          ll_vec = TRUE)
  })
    
  ## Should return a vector of correct length
  expect_true(is.numeric(ll_vec_result))
  expect_equal(length(ll_vec_result), n)
  
  suppressWarnings({
    ## Get total log-likelihood
    ll_total <- compute_loglik_fixed(dcc_fit, params = as.list(est_params))
  })
    
  ## Sum should equal total
  expect_equal(sum(ll_vec_result), ll_total, tolerance = 1e-6)
})

test_that("compute_loglik_fixed ll_vec works for Copula-GARCH dynamic", {
  skip_on_cran()
  
  ## Generate small sample data
  set.seed(111)
  n <- 100
  returns <- matrix(rnorm(n * 2), ncol = 2)
  returns <- xts(returns, order.by = seq.Date(Sys.Date() - n + 1, 
                                              Sys.Date(), by = "day"))
  colnames(returns) <- c("asset1", "asset2")
  
  ## Estimate univariate GARCH models
  spec1 <- garch_modelspec(returns[,1], model = "garch", order = c(1,1))
  spec2 <- garch_modelspec(returns[,2], model = "garch", order = c(1,1))
  fit1 <- estimate(spec1, keep_tmb = TRUE)
  fit2 <- estimate(spec2, keep_tmb = TRUE)
  
  ## Combine into multivariate
  garch_fits <- to_multi_estimate(list(fit1, fit2))
  
  ## Estimate Copula-GARCH model
  cgarch_spec <- cgarch_modelspec(garch_fits, dynamics = "dcc", dcc_order = c(1,1))
  cgarch_fit <- estimate(cgarch_spec)
  
  ## Get estimated parameters
  est_params <- coef(cgarch_fit)
  
  ## Get per-observation log-likelihoods
  ll_vec_result <- compute_loglik_fixed(cgarch_fit, params = as.list(est_params), 
                                        ll_vec = TRUE)
  
  ## Should return a vector
  expect_true(is.numeric(ll_vec_result))
  expect_equal(length(ll_vec_result), n)
  
  ## Get total log-likelihood
  ll_total <- compute_loglik_fixed(cgarch_fit, params = as.list(est_params))
  
  ## Sum should equal total
  expect_equal(sum(ll_vec_result), ll_total, tolerance = 1e-6)
})

test_that("compute_loglik_fixed ll_vec works for Copula-GARCH constant", {
  skip_on_cran()
  
  ## Generate small sample data
  set.seed(222)
  n <- 100
  returns <- matrix(rnorm(n * 2), ncol = 2)
  returns <- xts(returns, order.by = seq.Date(Sys.Date() - n + 1, 
                                              Sys.Date(), by = "day"))
  colnames(returns) <- c("asset1", "asset2")
  
  ## Estimate univariate GARCH models
  spec1 <- garch_modelspec(returns[,1], model = "garch", order = c(1,1))
  spec2 <- garch_modelspec(returns[,2], model = "garch", order = c(1,1))
  fit1 <- estimate(spec1, keep_tmb = TRUE)
  fit2 <- estimate(spec2, keep_tmb = TRUE)
  
  ## Combine into multivariate
  garch_fits <- to_multi_estimate(list(fit1, fit2))
  
  ## Estimate constant Copula-GARCH model
  cgarch_spec <- cgarch_modelspec(garch_fits, dynamics = "constant")
  cgarch_fit <- estimate(cgarch_spec)
  
  ## Get estimated parameters
  est_params <- coef(cgarch_fit)
  
  ## Get per-observation log-likelihoods
  ll_vec_result <- compute_loglik_fixed(cgarch_fit, params = as.list(est_params), 
                                        ll_vec = TRUE)
  
  ## Should return a vector
  expect_true(is.numeric(ll_vec_result))
  expect_equal(length(ll_vec_result), n)
  
  ## Get total log-likelihood
  ll_total <- compute_loglik_fixed(cgarch_fit, params = as.list(est_params))
  
  ## Sum should equal total
  expect_equal(sum(ll_vec_result), ll_total, tolerance = 1e-6)
})

test_that("compute_loglik_fixed errors when ll_vec and return_components both TRUE", {
  skip_on_cran()
  
  ## Generate small sample data
  set.seed(333)
  n <- 50
  returns <- matrix(rnorm(n * 2), ncol = 2)
  returns <- xts(returns, order.by = seq.Date(Sys.Date() - n + 1, 
                                              Sys.Date(), by = "day"))
  colnames(returns) <- c("series1", "series2")
  
  ## Estimate univariate GARCH models
  spec1 <- garch_modelspec(returns[,1], model = "garch", order = c(1,1))
  spec2 <- garch_modelspec(returns[,2], model = "garch", order = c(1,1))
  fit1 <- estimate(spec1, keep_tmb = TRUE)
  fit2 <- estimate(spec2, keep_tmb = TRUE)
  
  ## Combine into multivariate
  garch_fits <- to_multi_estimate(list(fit1, fit2))
  
  ## Estimate DCC model
  dcc_spec <- dcc_modelspec(garch_fits, dynamics = "dcc", dcc_order = c(1,1))
  dcc_fit <- estimate(dcc_spec)
  
  ## Get estimated parameters
  est_params <- coef(dcc_fit)
  
  ## Should error when both ll_vec and return_components are TRUE
  expect_error(
    compute_loglik_fixed(dcc_fit, params = as.list(est_params), 
                         ll_vec = TRUE, return_components = TRUE),
    "Cannot use both return_components = TRUE and ll_vec = TRUE"
  )
})

test_that("compute_loglik_fixed ll_vec values are reasonable", {
  skip_on_cran()
  
  ## Generate small sample data
  set.seed(444)
  n <- 100
  returns <- matrix(rnorm(n * 2), ncol = 2)
  returns <- xts(returns, order.by = seq.Date(Sys.Date() - n + 1, 
                                              Sys.Date(), by = "day"))
  colnames(returns) <- c("asset1", "asset2")
  
  ## Estimate univariate GARCH models
  spec1 <- garch_modelspec(returns[,1], model = "garch", order = c(1,1))
  spec2 <- garch_modelspec(returns[,2], model = "garch", order = c(1,1))
  fit1 <- estimate(spec1, keep_tmb = TRUE)
  fit2 <- estimate(spec2, keep_tmb = TRUE)
  
  ## Combine into multivariate
  garch_fits <- to_multi_estimate(list(fit1, fit2))
  
  ## Estimate DCC model
  dcc_spec <- dcc_modelspec(garch_fits, dynamics = "dcc", dcc_order = c(1,1))
  dcc_fit <- estimate(dcc_spec)
  
  ## Get estimated parameters
  est_params <- coef(dcc_fit)
  
  ## Get per-observation log-likelihoods
  ll_vec_result <- compute_loglik_fixed(dcc_fit, params = as.list(est_params), 
                                        ll_vec = TRUE)
  
  ## All values should be finite (not Inf, -Inf, or NaN)
  expect_true(all(is.finite(ll_vec_result)))
  
  ## Log-likelihoods should be negative (since they're probabilities < 1)
  expect_true(all(ll_vec_result < 0))
  
  ## Should not have any extreme outliers (within reasonable range)
  ## For normal distributions, log-likelihoods typically range from -10 to 0
  expect_true(all(ll_vec_result > -100))
})

test_that("compute_loglik_fixed ll_vec works with alternative parameters", {
  skip_on_cran()
  
  ## Generate small sample data
  set.seed(555)
  n <- 100
  returns <- matrix(rnorm(n * 2), ncol = 2)
  returns <- xts(returns, order.by = seq.Date(Sys.Date() - n + 1, 
                                              Sys.Date(), by = "day"))
  colnames(returns) <- c("asset1", "asset2")
  
  ## Estimate univariate GARCH models
  spec1 <- garch_modelspec(returns[,1], model = "garch", order = c(1,1))
  spec2 <- garch_modelspec(returns[,2], model = "garch", order = c(1,1))
  fit1 <- estimate(spec1, keep_tmb = TRUE)
  fit2 <- estimate(spec2, keep_tmb = TRUE)
  
  ## Combine into multivariate
  garch_fits <- to_multi_estimate(list(fit1, fit2))
  
  ## Estimate DCC model
  dcc_spec <- dcc_modelspec(garch_fits, dynamics = "dcc", dcc_order = c(1,1))
  dcc_fit <- estimate(dcc_spec)
  
  ## Test with alternative parameters
  alt_params <- list(alpha_1 = 0.05, beta_1 = 0.90)
  
  ll_vec_alt <- compute_loglik_fixed(dcc_fit, params = alt_params, ll_vec = TRUE)
  ll_total_alt <- compute_loglik_fixed(dcc_fit, params = alt_params)
  
  ## Should return a vector
  expect_true(is.numeric(ll_vec_alt))
  expect_equal(length(ll_vec_alt), n)
  
  ## Sum should equal total
  expect_equal(sum(ll_vec_alt), ll_total_alt, tolerance = 1e-6)
  
  ## Get estimated parameter results
  est_params <- coef(dcc_fit)
  ll_vec_est <- compute_loglik_fixed(dcc_fit, params = as.list(est_params), 
                                     ll_vec = TRUE)
  
  ## Per-observation likelihoods should differ
  expect_false(isTRUE(all.equal(ll_vec_alt, ll_vec_est)))
})
