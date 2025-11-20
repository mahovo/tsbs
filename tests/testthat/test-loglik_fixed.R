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
  ll_at_est <- compute_loglik_fixed(
    dcc_fit,
    params = list(alpha_1 = est_params["alpha_1"],
    beta_1 = est_params["beta_1"])
  )
  
  ## Compute at different parameters (should be lower)
  ll_at_other <- compute_loglik_fixed(
    dcc_fit,
    params = list(alpha_1 = 0.03, beta_1 = 0.95)
  )
  
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

## Test GOGARCH Model ==================================================

## Unit tests for .compute_gogarch_loglik()
# Setup: Create a GOGARCH fit object for testing
setup_gogarch_test <- function() {
  suppressWarnings({
    set.seed(123)
    n <- 100
    returns <- matrix(rnorm(n * 2), ncol = 2)
    returns <- xts::xts(returns, order.by = seq.Date(Sys.Date() - n + 1, Sys.Date(), by = "day"))
    colnames(returns) <- c("asset1", "asset2")
    
    gogarch_spec <- tsmarch::gogarch_modelspec(returns, distribution = "norm", 
                                               model = "garch", order = c(1, 1), 
                                               components = 2)
    gogarch_fit <- estimate(gogarch_spec)
    
    return(gogarch_fit)
  })
}

test_that("GOGARCH: scalar likelihood at estimated parameters matches stored value", {
  gogarch_fit <- setup_gogarch_test()
  
  ll_computed <- suppressWarnings(compute_loglik_fixed(gogarch_fit, params = list()))
  ll_stored <- gogarch_fit$loglik
  
  expect_equal(ll_computed, ll_stored, tolerance = 1e-10)
})

test_that("GOGARCH: vector likelihood sums to scalar likelihood at estimated parameters", {
  gogarch_fit <- setup_gogarch_test()
  
  ll_scalar <- suppressWarnings(compute_loglik_fixed(gogarch_fit, params = list()))
  ll_vector <- suppressWarnings(compute_loglik_fixed(gogarch_fit, params = list(), ll_vec = TRUE))
  
  expect_equal(sum(ll_vector), ll_scalar, tolerance = 1e-10)
  expect_equal(length(ll_vector), 100)
})

test_that("GOGARCH: likelihood decreases with non-optimal parameters", {
  gogarch_fit <- setup_gogarch_test()
  
  ll_estimated <- suppressWarnings(compute_loglik_fixed(gogarch_fit, params = list()))
  
  ## Test with different parameters
  test_params <- list(
    omega_1 = 0.03, alpha_1 = 0.04, beta_1 = 0.95,
    omega_2 = 0.35, alpha_2 = 0.00, beta_2 = 0.65
  )
  ll_test <- suppressWarnings(compute_loglik_fixed(gogarch_fit, params = test_params))
  
  expect_true(ll_estimated > ll_test)
})

test_that("GOGARCH: vector mode with changed parameters is consistent", {
  gogarch_fit <- setup_gogarch_test()
  
  test_params <- list(
    omega_1 = 0.03, alpha_1 = 0.04, beta_1 = 0.95,
    omega_2 = 0.35, alpha_2 = 0.00, beta_2 = 0.65
  )
  
  ll_scalar <- suppressWarnings(compute_loglik_fixed(gogarch_fit, params = test_params))
  ll_vector <- suppressWarnings(compute_loglik_fixed(gogarch_fit, params = test_params, ll_vec = TRUE))
  
  expect_equal(sum(ll_vector), ll_scalar, tolerance = 1e-10)
  expect_equal(length(ll_vector), 100)
})

test_that("GOGARCH: partial parameter update works correctly", {
  gogarch_fit <- setup_gogarch_test()
  
  ll_estimated <- suppressWarnings(compute_loglik_fixed(gogarch_fit, params = list()))
  
  ## Update only component 1 parameters
  partial_params <- list(
    omega_1 = 0.05, alpha_1 = 0.03, beta_1 = 0.92
  )
  ll_partial <- suppressWarnings(compute_loglik_fixed(gogarch_fit, params = partial_params))
  
  ## Should be different from estimated
  expect_false(isTRUE(all.equal(ll_partial, ll_estimated, tolerance = 1e-10)))
  
  ## Should be numeric and finite
  expect_true(is.numeric(ll_partial))
  expect_true(is.finite(ll_partial))
})

test_that("GOGARCH: return_components gives correct decomposition", {
  gogarch_fit <- setup_gogarch_test()
  
  result <- suppressWarnings(compute_loglik_fixed(gogarch_fit, params = list(), 
                                                  return_components = TRUE))
  
  ## Check structure
  expect_true(is.list(result))
  expect_named(result, c("loglik", "component_logliks", "jacobian_adjustment"))
  
  ## Check that manual sum equals total
  manual_sum <- sum(result$component_logliks) + result$jacobian_adjustment
  expect_equal(result$loglik, manual_sum, tolerance = 1e-10)
  
  ## Check that total matches stored value
  expect_equal(result$loglik, gogarch_fit$loglik, tolerance = 1e-10)
  
  ## Check component dimensions
  expect_equal(length(result$component_logliks), 2)
  expect_true(is.numeric(result$jacobian_adjustment))
})

test_that("GOGARCH: vector + components mode has consistent structure", {
  gogarch_fit <- setup_gogarch_test()
  
  result <- suppressWarnings(compute_loglik_fixed(gogarch_fit, params = list(), 
                                                  return_components = TRUE, ll_vec = TRUE))
  
  ## Check structure
  expect_true(is.list(result))
  expect_named(result, c("loglik", "lik_vector", "component_logliks", 
                         "jacobian_adjustment"))
  
  ## Check consistency
  expect_equal(result$loglik, sum(result$lik_vector), tolerance = 1e-10)
  expect_equal(result$loglik, gogarch_fit$loglik, tolerance = 1e-10)
  
  ## Check dimensions
  expect_equal(length(result$lik_vector), 100)
  expect_equal(length(result$component_logliks), 2)
})

test_that("GOGARCH: Jacobian adjustment is computed correctly", {
  gogarch_fit <- setup_gogarch_test()
  
  result <- suppressWarnings(compute_loglik_fixed(gogarch_fit, params = list(), 
                                                  return_components = TRUE))
  
  ## Manually compute Jacobian
  K <- gogarch_fit$ica$K
  if (nrow(K) == ncol(K)) {
    expected_jacobian <- log(abs(det(K)))
  } else {
    expected_jacobian <- log(abs(det(K %*% t(K))))
  }
  
  expect_equal(result$jacobian_adjustment, expected_jacobian, tolerance = 1e-10)
})

test_that("GOGARCH: component log-likelihoods sum correctly with Jacobian", {
  gogarch_fit <- setup_gogarch_test()
  
  result <- suppressWarnings(compute_loglik_fixed(gogarch_fit, params = list(), 
                                                  return_components = TRUE))
  
  ## Check that component logliks plus Jacobian equals total
  manual_total <- sum(result$component_logliks) + result$jacobian_adjustment
  
  expect_equal(result$loglik, manual_total, tolerance = 1e-10)
  
  ## Check that all components are finite
  expect_true(all(is.finite(result$component_logliks)))
  expect_equal(length(result$component_logliks), 2)
  
  ## Check that components are negative (log-likelihoods typically negative)
  ## This is a sanity check for proper data
  expect_true(all(result$component_logliks < 0))
})

test_that("GOGARCH: observation-wise likelihoods are all finite", {
  gogarch_fit <- setup_gogarch_test()
  
  ll_vector <- suppressWarnings(compute_loglik_fixed(gogarch_fit, params = list(), ll_vec = TRUE))
  
  expect_true(all(is.finite(ll_vector)))
  expect_true(is.numeric(ll_vector))
  expect_equal(length(ll_vector), 100)
})

test_that("GOGARCH: parameter names are correctly matched", {
  gogarch_fit <- setup_gogarch_test()
  
  ## Test with explicit parameter names for both components
  params_explicit <- list(
    omega_1 = 0.04, alpha_1 = 0.05, beta_1 = 0.90,
    omega_2 = 0.30, alpha_2 = 0.05, beta_2 = 0.70
  )
  
  ll_explicit <- suppressWarnings(compute_loglik_fixed(gogarch_fit, params = params_explicit))
  
  ## Should be different from estimated
  ll_estimated <- suppressWarnings(compute_loglik_fixed(gogarch_fit, params = list()))
  expect_false(isTRUE(all.equal(ll_explicit, ll_estimated, tolerance = 1e-6)))
  
  ## Should be finite and numeric
  expect_true(is.finite(ll_explicit))
  expect_true(is.numeric(ll_explicit))
})

test_that("GOGARCH: empty params list uses estimated parameters", {
  gogarch_fit <- setup_gogarch_test()
  
  ll_empty <- suppressWarnings(compute_loglik_fixed(gogarch_fit, params = list()))
  ll_stored <- gogarch_fit$loglik
  
  ## Should exactly match stored value
  expect_equal(ll_empty, ll_stored, tolerance = 1e-10)
})

test_that("GOGARCH: vector mode distributes Jacobian adjustment correctly", {
  gogarch_fit <- setup_gogarch_test()
  
  result <- suppressWarnings(compute_loglik_fixed(gogarch_fit, params = list(), 
                                                  return_components = TRUE, ll_vec = TRUE))
  
  n_obs <- length(result$lik_vector)
  
  ## The Jacobian adjustment should be distributed equally across observations
  ## So the mean contribution per observation should be jacobian/n_obs
  expected_jacobian_per_obs <- result$jacobian_adjustment / n_obs
  
  ## Get component contributions without Jacobian
  suppressWarnings({
    comp1_vec <- log(gogarch_fit$univariate[[1]]$tmb$report(
      get("last.par.best", envir = gogarch_fit$univariate[[1]]$tmb$env)
    )$ll_vector)
    comp2_vec <- log(gogarch_fit$univariate[[2]]$tmb$report(
      get("last.par.best", envir = gogarch_fit$univariate[[2]]$tmb$env)
    )$ll_vector)
  })
  
  ## Manual calculation: sum of components plus distributed Jacobian
  manual_vec <- comp1_vec + comp2_vec + expected_jacobian_per_obs
  
  expect_equal(result$lik_vector, manual_vec, tolerance = 1e-10)
})

test_that("GOGARCH: handles different parameter combinations correctly", {
  gogarch_fit <- setup_gogarch_test()
  
  ## Test 1: Only omega parameters
  params1 <- list(omega_1 = 0.05, omega_2 = 0.30)
  ll1 <- suppressWarnings(compute_loglik_fixed(gogarch_fit, params = params1))
  
  ## Test 2: Only alpha parameters
  params2 <- list(alpha_1 = 0.03, alpha_2 = 0.02)
  ll2 <- suppressWarnings(compute_loglik_fixed(gogarch_fit, params = params2))
  
  ## Test 3: Only beta parameters
  params3 <- list(beta_1 = 0.90, beta_2 = 0.85)
  ll3 <- suppressWarnings(compute_loglik_fixed(gogarch_fit, params = params3))
  
  ## All should be different
  expect_false(isTRUE(all.equal(ll1, ll2, tolerance = 1e-6)))
  expect_false(isTRUE(all.equal(ll2, ll3, tolerance = 1e-6)))
  expect_false(isTRUE(all.equal(ll1, ll3, tolerance = 1e-6)))
  
  ## All should be finite
  expect_true(all(is.finite(c(ll1, ll2, ll3))))
})

test_that("GOGARCH: vector mode with partial parameters is consistent", {
  gogarch_fit <- setup_gogarch_test()
  
  partial_params <- list(omega_1 = 0.04, beta_2 = 0.70)
  
  ll_scalar <- suppressWarnings(compute_loglik_fixed(gogarch_fit, params = partial_params))
  ll_vector <- suppressWarnings(compute_loglik_fixed(gogarch_fit, params = partial_params, 
                                                     ll_vec = TRUE))
  
  expect_equal(sum(ll_vector), ll_scalar, tolerance = 1e-10)
})

test_that("GOGARCH: works with different model specifications", {
  ## This test ensures the implementation is robust to different GOGARCH setups
  suppressWarnings({
    set.seed(456)
    n <- 50
    returns <- matrix(rnorm(n * 3), ncol = 3)  # 3 assets instead of 2
    returns <- xts::xts(returns, order.by = seq.Date(Sys.Date() - n + 1, Sys.Date(), by = "day"))
    colnames(returns) <- c("asset1", "asset2", "asset3")
    
    gogarch_spec <- tsmarch::gogarch_modelspec(returns, distribution = "norm", 
                                               model = "garch", order = c(1, 1), 
                                               components = 3)
    gogarch_fit <- estimate(gogarch_spec)
    
    ll_computed <- compute_loglik_fixed(gogarch_fit, params = list())
    ll_stored <- gogarch_fit$loglik
  })
  
  expect_equal(ll_computed, ll_stored, tolerance = 1e-10)
  
  ## Test with parameters
  test_params <- list(
    omega_1 = 0.03, alpha_1 = 0.04, beta_1 = 0.90,
    omega_2 = 0.03, alpha_2 = 0.04, beta_2 = 0.90,
    omega_3 = 0.03, alpha_3 = 0.04, beta_3 = 0.90
  )
  ll_test <- suppressWarnings(compute_loglik_fixed(gogarch_fit, params = test_params))
  expect_true(is.finite(ll_test))
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
    compute_loglik_fixed(
      dcc_fit,
      params = list(alpha_1 = a, beta_1 = beta_est)
    )
  })
  
  ## Should all be finite
  expect_true(
    all(is.finite(profile_ll)),
    info = "All profile likelihoods should be finite"
  )
  
  ## The likelihood at the estimated parameter should be among the highest
  ## (allowing for numerical tolerance and that we're fixing beta)
  max_ll <- max(profile_ll)
  ll_at_est <- profile_ll[3]  # Middle point is the estimated value
  
  ## Check that estimated is within a reasonable range of maximum
  expect_true(
    ll_at_est >= max_ll - 1,
    info = paste("LL at estimated alpha should be near maximum.",
      "Estimated:", ll_at_est, "Max:", max_ll)
    )
})

## Likelihood Ratio Test Example ============================================

test_that("likelihood ratio test framework works", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  skip_if_not_installed("xts")
  skip_on_cran()
  
  dcc_fit <- readRDS(file.path(.get_fixtures_dir(), "dcc_dynamic_fit.rds"))
  
  ## H0: alpha = 0.05, beta = 0.90 (restricted)
  ll_restricted <- compute_loglik_fixed(
    dcc_fit,
    params = list(alpha_1 = 0.05, beta_1 = 0.90)
  )
  
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
  ll <- compute_loglik_fixed(
    dcc_fit, 
    params = list(alpha_1 = 0.03, beta_1 = 0.95)
  )
  
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
    returns <- xts(
      returns, 
      order.by = seq.Date(Sys.Date() - n + 1, Sys.Date(), by = "day")
    )
    colnames(returns) <- c("series1", "series2")
    
    ## Estimate univariate GARCH models
    spec1 <- garch_modelspec(returns[ ,1], model = "garch", order = c(1, 1))
    spec2 <- garch_modelspec(returns[ ,2], model = "garch", order = c(1, 1))
    fit1 <- estimate(spec1, keep_tmb = TRUE)
    fit2 <- estimate(spec2, keep_tmb = TRUE)
    
    ## Combine into multivariate
    garch_fits <- to_multi_estimate(list(fit1, fit2))
    
    ## Estimate DCC model
    dcc_spec <- dcc_modelspec(garch_fits, dynamics = "dcc", dcc_order = c(1, 1))
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
    returns <- xts(
      returns, 
      order.by = seq.Date(Sys.Date() - n + 1, Sys.Date(), by = "day")
    )
    colnames(returns) <- c("asset1", "asset2")
    
    ## Estimate univariate GARCH models
    spec1 <- garch_modelspec(returns[ ,1], model = "garch", order = c(1, 1))
    spec2 <- garch_modelspec(returns[ ,2], model = "garch", order = c(1, 1))
    fit1 <- estimate(spec1, keep_tmb = TRUE)
    fit2 <- estimate(spec2, keep_tmb = TRUE)
    
    ## Combine into multivariate
    garch_fits <- to_multi_estimate(list(fit1, fit2))
    
    ## Estimate DCC model
    dcc_spec <- dcc_modelspec(garch_fits, dynamics = "dcc", dcc_order = c(1, 1))
    dcc_fit <- estimate(dcc_spec)
    
    ## Get estimated parameters
    est_params <- coef(dcc_fit)
    
    ## Get per-observation log-likelihoods
    ll_vec_result <- compute_loglik_fixed(
      dcc_fit, 
      params = as.list(est_params),
      ll_vec = TRUE
    )
    
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
    returns <- xts(
      returns, 
      order.by = seq.Date(Sys.Date() - n + 1, Sys.Date(), by = "day")
    )
    colnames(returns) <- c("asset1", "asset2")
    
    ## Estimate univariate GARCH models
    spec1 <- garch_modelspec(returns[ ,1], model = "garch", order = c(1, 1))
    spec2 <- garch_modelspec(returns[ ,2], model = "garch", order = c(1, 1))
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
    ll_vec_result <- compute_loglik_fixed(
      dcc_fit, 
      params = as.list(est_params),
      ll_vec = TRUE
    )
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
  
  suppressWarnings({
    ## Generate small sample data
    set.seed(111)
    n <- 100
    returns <- matrix(rnorm(n * 2), ncol = 2)
    returns <- xts(
      returns, 
      order.by = seq.Date(Sys.Date() - n + 1, Sys.Date(), by = "day")
    )
    colnames(returns) <- c("asset1", "asset2")
    
    ## Estimate univariate GARCH models
    spec1 <- garch_modelspec(returns[ ,1], model = "garch", order = c(1, 1))
    spec2 <- garch_modelspec(returns[ ,2], model = "garch", order = c(1, 1))
    fit1 <- estimate(spec1, keep_tmb = TRUE)
    fit2 <- estimate(spec2, keep_tmb = TRUE)
    
    ## Combine into multivariate
    garch_fits <- to_multi_estimate(list(fit1, fit2))
    
    ## Estimate Copula-GARCH model
    cgarch_spec <- cgarch_modelspec(
      garch_fits, 
      dynamics = "dcc", 
      dcc_order = c(1, 1)
    )
    cgarch_fit <- estimate(cgarch_spec)
    
    ## Get estimated parameters
    est_params <- coef(cgarch_fit)
    
    ## Get per-observation log-likelihoods
    ll_vec_result <- compute_loglik_fixed(
      cgarch_fit, 
      params = as.list(est_params),
      ll_vec = TRUE
    )
  })
    
  ## Should return a vector
  expect_true(is.numeric(ll_vec_result))
  expect_equal(length(ll_vec_result), n)
  
  suppressWarnings({
    ## Get total log-likelihood
    ll_total <- compute_loglik_fixed(cgarch_fit, params = as.list(est_params))
  })
  
  ## Sum should equal total
  expect_equal(sum(ll_vec_result), ll_total, tolerance = 1e-6)
})

test_that("compute_loglik_fixed ll_vec works for Copula-GARCH constant", {
  skip_on_cran()
  
  suppressWarnings({
    ## Generate small sample data
    set.seed(222)
    n <- 100
    returns <- matrix(rnorm(n * 2), ncol = 2)
    returns <- xts(
      returns, 
      order.by = seq.Date(Sys.Date() - n + 1, Sys.Date(), by = "day")
    )
    colnames(returns) <- c("asset1", "asset2")
    
    ## Estimate univariate GARCH models
    spec1 <- garch_modelspec(returns[ ,1], model = "garch", order = c(1, 1))
    spec2 <- garch_modelspec(returns[ ,2], model = "garch", order = c(1, 1))
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
    ll_vec_result <- compute_loglik_fixed(
      cgarch_fit, 
      params = as.list(est_params),
      ll_vec = TRUE
    )
  })
    
  ## Should return a vector
  expect_true(is.numeric(ll_vec_result))
  expect_equal(length(ll_vec_result), n)
  
  suppressWarnings({
    ## Get total log-likelihood
    ll_total <- compute_loglik_fixed(cgarch_fit, params = as.list(est_params))
  })
    
  ## Sum should equal total
  expect_equal(sum(ll_vec_result), ll_total, tolerance = 1e-6)
})

test_that("compute_loglik_fixed errors when ll_vec and return_components both TRUE", {
  skip_on_cran()
  
  suppressWarnings({
    ## Generate small sample data
    set.seed(333)
    n <- 50
    returns <- matrix(rnorm(n * 2), ncol = 2)
    returns <- xts(
      returns, 
      order.by = seq.Date(Sys.Date() - n + 1, Sys.Date(), by = "day")
    )
    colnames(returns) <- c("series1", "series2")
    
    ## Estimate univariate GARCH models
    spec1 <- garch_modelspec(returns[ ,1], model = "garch", order = c(1, 1))
    spec2 <- garch_modelspec(returns[ ,2], model = "garch", order = c(1, 1))
    fit1 <- estimate(spec1, keep_tmb = TRUE)
    fit2 <- estimate(spec2, keep_tmb = TRUE)
    
    ## Combine into multivariate
    garch_fits <- to_multi_estimate(list(fit1, fit2))
    
    ## Estimate DCC model
    dcc_spec <- dcc_modelspec(garch_fits, dynamics = "dcc", dcc_order = c(1, 1))
    dcc_fit <- estimate(dcc_spec)
    
    ## Get estimated parameters
    est_params <- coef(dcc_fit)
  })
  
  ## Should error when both ll_vec and return_components are TRUE
  expect_error(
    compute_loglik_fixed(
      dcc_fit, 
      params = as.list(est_params), 
      ll_vec = TRUE, 
      return_components = TRUE
    ),
    "Cannot use both return_components = TRUE and ll_vec = TRUE"
  )
})

test_that("compute_loglik_fixed ll_vec values are reasonable", {
  skip_on_cran()
  
  suppressWarnings({
    ## Generate small sample data
    set.seed(444)
    n <- 100
    returns <- matrix(rnorm(n * 2), ncol = 2)
    returns <- xts(
      returns, 
      order.by = seq.Date(Sys.Date() - n + 1, Sys.Date(), by = "day")
    )
    colnames(returns) <- c("asset1", "asset2")
    
    ## Estimate univariate GARCH models
    spec1 <- garch_modelspec(returns[ ,1], model = "garch", order = c(1, 1))
    spec2 <- garch_modelspec(returns[ ,2], model = "garch", order = c(1, 1))
    fit1 <- estimate(spec1, keep_tmb = TRUE)
    fit2 <- estimate(spec2, keep_tmb = TRUE)
    
    ## Combine into multivariate
    garch_fits <- to_multi_estimate(list(fit1, fit2))
    
    ## Estimate DCC model
    dcc_spec <- dcc_modelspec(garch_fits, dynamics = "dcc", dcc_order = c(1, 1))
    dcc_fit <- estimate(dcc_spec)
    
    ## Get estimated parameters
    est_params <- coef(dcc_fit)
    
    ## Get per-observation log-likelihoods
    ll_vec_result <- compute_loglik_fixed(
      dcc_fit, 
      params = as.list(est_params),
      ll_vec = TRUE
    )
  })
  
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
  
  suppressWarnings({
    ## Generate small sample data
    set.seed(555)
    n <- 100
    returns <- matrix(rnorm(n * 2), ncol = 2)
    returns <- xts(
      returns, 
      order.by = seq.Date(Sys.Date() - n + 1, Sys.Date(), by = "day")
    )
    colnames(returns) <- c("asset1", "asset2")
    
    ## Estimate univariate GARCH models
    spec1 <- garch_modelspec(returns[ ,1], model = "garch", order = c(1, 1))
    spec2 <- garch_modelspec(returns[ ,2], model = "garch", order = c(1, 1))
    fit1 <- estimate(spec1, keep_tmb = TRUE)
    fit2 <- estimate(spec2, keep_tmb = TRUE)
    
    ## Combine into multivariate
    garch_fits <- to_multi_estimate(list(fit1, fit2))
    
    ## Estimate DCC model
    dcc_spec <- dcc_modelspec(garch_fits, dynamics = "dcc", dcc_order = c(1, 1))
    dcc_fit <- estimate(dcc_spec)
    
    ## Test with alternative parameters
    alt_params <- list(alpha_1 = 0.05, beta_1 = 0.90)
    
    ll_vec_alt <- compute_loglik_fixed(
      dcc_fit, 
      params = alt_params, 
      ll_vec = TRUE
    )
    ll_total_alt <- compute_loglik_fixed(dcc_fit, params = alt_params)
  })
  
  ## Should return a vector
  expect_true(is.numeric(ll_vec_alt))
  expect_equal(length(ll_vec_alt), n)
  
  ## Sum should equal total
  expect_equal(sum(ll_vec_alt), ll_total_alt, tolerance = 1e-6)
  
  suppressWarnings({
    ## Get estimated parameter results
    est_params <- coef(dcc_fit)
    ll_vec_est <- compute_loglik_fixed(
      dcc_fit, 
      params = as.list(est_params),
      ll_vec = TRUE
    )
  })
  
  ## Per-observation likelihoods should differ
  expect_false(isTRUE(all.equal(ll_vec_alt, ll_vec_est)))
})
