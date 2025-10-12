# Test file for compute_loglik_fixed function

# Helper to get fixtures directory
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

# Setup: Create mock objects for testing
# Note: These tests will need actual tsmarch objects for integration testing

test_that("compute_loglik_fixed validates input object class", {
  
  # Test with invalid object class
  expect_error(
    compute_loglik_fixed(
      object = list(some = "data"),
      params = list(alpha_1 = 0.05)
    ),
    "must be an estimated tsmarch model"
  )
  
  # Test with NULL object
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
  
  # Load a fixture
  dcc_fit <- readRDS(file.path(.get_fixtures_dir(), "dcc_dynamic_fit.rds"))
  
  # Test with non-list params
  expect_error(
    compute_loglik_fixed(
      object = dcc_fit,
      params = "not a list"
    ),
    "params must be a named list"
  )
  
  # Test with numeric vector instead of list
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
  
  # Create a simple parmatrix
  parmatrix <- data.table::data.table(
    parameter = c("alpha_1", "beta_1", "shape"),
    value = c(0.05, 0.90, 5.0),
    lower = c(0, 0, 2.01),
    upper = c(1, 1, 50),
    estimate = c(1, 1, 1)
  )
  
  # Update parameters
  new_params <- list(alpha_1 = 0.03, beta_1 = 0.95)
  updated <- .update_parmatrix(parmatrix, new_params)
  
  # Check updates
  expect_equal(updated[parameter == "alpha_1"]$value, 0.03)
  expect_equal(updated[parameter == "beta_1"]$value, 0.95)
  expect_equal(updated[parameter == "shape"]$value, 5.0)  # Unchanged
})

test_that(".update_parmatrix warns about unknown parameters", {
  skip_if_not_installed("data.table")
  
  parmatrix <- data.table::data.table(
    parameter = c("alpha_1", "beta_1"),
    value = c(0.05, 0.90),
    estimate = c(1, 1)
  )
  
  # Try to update non-existent parameter
  expect_warning(
    .update_parmatrix(parmatrix, list(gamma_1 = 0.02)),
    "Parameter 'gamma_1' not found"
  )
})

test_that(".update_parmatrix handles multiple matching parameters", {
  skip_if_not_installed("data.table")
  
  # Edge case: duplicate parameter names (shouldn't happen in practice)
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

# ============================================================================
# Integration tests with actual tsmarch objects
# ============================================================================

# Test DCC Dynamic Model ====================================================

test_that("compute_loglik_fixed works with DCC dynamic (mvn)", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  skip_if_not_installed("xts")
  skip_on_cran()
  
  # Get or create test fixture
  dcc_fit <- readRDS(file.path(.get_fixtures_dir(), "dcc_dynamic_fit.rds"))
  
  # Get estimated parameters
  est_params <- coef(dcc_fit)
  
  # Convert to named list
  params_list <- as.list(est_params)
  names(params_list) <- names(est_params)
  
  # Compute log-likelihood at estimated parameters
  ll_fixed <- compute_loglik_fixed(dcc_fit, params = params_list)
  
  # Get the estimated log-likelihood
  #ll_estimated <- as.numeric(logLik(dcc_fit))
  series1_nll <- dcc_fit$spec$univariate$series1$loglik
  series2_nll <- dcc_fit$spec$univariate$series2$loglik
  mv_nll <- dcc_fit$loglik
  nll_estimated <- series1_nll + series2_nll + mv_nll
  ll_estimated <- -nll_estimated
  
  # Should match within numerical tolerance
  expect_equal(ll_fixed, ll_estimated, tolerance = 1e-4,
               info = "Log-likelihood at estimated params should match estimation")
})

test_that("compute_loglik_fixed detects parameter changes in DCC dynamic", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  skip_if_not_installed("xts")
  skip_on_cran()
  
  dcc_fit <- readRDS(file.path(.get_fixtures_dir(), "dcc_dynamic_fit.rds"))
  
  # Get estimated parameters
  est_params <- coef(dcc_fit)
  
  # Compute at estimated parameters
  ll_at_est <- compute_loglik_fixed(dcc_fit, 
                                    params = list(alpha_1 = est_params["alpha_1"],
                                                  beta_1 = est_params["beta_1"]))
  
  # Compute at different parameters (should be lower)
  ll_at_other <- compute_loglik_fixed(dcc_fit,
                                      params = list(alpha_1 = 0.03,
                                                    beta_1 = 0.95))
  
  # Likelihood at estimated parameters should be higher
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
  
  # Get components
  ll_comp <- compute_loglik_fixed(
    dcc_fit, 
    params = as.list(est_params),
    return_components = TRUE
  )
  
  # Check structure
  expect_type(ll_comp, "list")
  expect_named(ll_comp, c("loglik", "garch_loglik", "multivariate_loglik"))
  
  # Check that components sum to total
  sum_comp <- ll_comp$garch_loglik + ll_comp$multivariate_loglik
  expect_equal(ll_comp$loglik, sum_comp, tolerance = 1e-10)
  
  # Check against simple call
  ll_simple <- compute_loglik_fixed(dcc_fit, params = as.list(est_params))
  expect_equal(ll_comp$loglik, ll_simple, tolerance = 1e-10)
})

# Test DCC Constant Model ==================================================

test_that("compute_loglik_fixed works with DCC constant (mvn)", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  skip_if_not_installed("xts")
  skip_on_cran()
  
  dcc_fit <- readRDS(file.path(.get_fixtures_dir(), "dcc_constant_fit.rds"))
  
  # For constant MVN, there are no parameters to estimate
  # Should still compute likelihood correctly
  ll_fixed <- compute_loglik_fixed(dcc_fit, params = list())
  ll_estimated <- as.numeric(logLik(dcc_fit))
  
  expect_equal(ll_fixed, ll_estimated, tolerance = 1e-4)
})

# Test Student-t DCC Model =================================================

test_that("compute_loglik_fixed works with DCC Student-t", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  skip_if_not_installed("xts")
  skip_on_cran()
  
  dcc_fit <- readRDS(file.path(.get_fixtures_dir(), "dcc_student_fit.rds"))
  
  # Get estimated parameters (should include shape parameter)
  est_params <- coef(dcc_fit)
  
  # Should have alpha, beta, and shape
  expect_true("shape" %in% names(est_params))
  
  # Compute log-likelihood
  ll_fixed <- compute_loglik_fixed(dcc_fit, params = as.list(est_params))
  # ll_estimated <- as.numeric(logLik(dcc_fit))
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
  
  # Try different degrees of freedom
  shapes <- c(3, 5, 10, 30)
  lls <- sapply(shapes, function(nu) {
    compute_loglik_fixed(dcc_fit, 
                         params = list(
                           alpha_1 = est_params["alpha_1"],
                           beta_1 = est_params["beta_1"],
                           shape = nu
                         ))
  })
  
  # Should all be finite
  expect_true(all(is.finite(lls)))
  
  # Likelihood at estimated shape should be highest
  est_shape_idx <- which.min(abs(shapes - est_params["shape"]))
  if (est_shape_idx > 1 && est_shape_idx < length(shapes)) {
    # If estimated shape is in the interior of our grid, 
    # it should have higher likelihood than neighbors
    expect_true(lls[est_shape_idx] >= max(lls[est_shape_idx - 1], 
                                          lls[est_shape_idx + 1]))
  }
})

# Test Copula-GARCH Model ==================================================

test_that("compute_loglik_fixed works with Copula-GARCH dynamic", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  skip_if_not_installed("xts")
  skip_on_cran()
  
  cgarch_fit <- readRDS(file.path(.get_fixtures_dir(), "copula_dynamic_fit.rds"))
  
  # Get estimated parameters
  est_params <- coef(cgarch_fit)
  
  # Compute log-likelihood
  ll_fixed <- compute_loglik_fixed(cgarch_fit, params = as.list(est_params))
  # ll_estimated <- as.numeric(logLik(cgarch_fit))
  series1_nll <- cgarch_fit$spec$univariate$series1$loglik
  series2_nll <- cgarch_fit$spec$univariate$series2$loglik
  copula_nll <- cgarch_fit$copula_nll
  model_nll <- series1_nll + series2_nll + copula_nll
  ll_expected <- -model_nll
  
  expect_equal(ll_fixed, ll_expected, tolerance = 1e-4)
})

# Profile Likelihood Tests =================================================

test_that("profile likelihood is smooth around optimum", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  skip_if_not_installed("xts")
  skip_on_cran()
  
  dcc_fit <- readRDS(file.path(.get_fixtures_dir(), "dcc_dynamic_fit.rds"))
  est_params <- coef(dcc_fit)
  
  # Create grid around estimated alpha
  alpha_est <- est_params["alpha_1"]
  beta_est <- est_params["beta_1"]
  alpha_grid <- seq(alpha_est - 0.02, alpha_est + 0.02, length.out = 5)
  
  # Compute profile likelihood
  profile_ll <- sapply(alpha_grid, function(a) {
    compute_loglik_fixed(dcc_fit, 
                         params = list(alpha_1 = a, beta_1 = beta_est))
  })
  
  # Should all be finite
  expect_true(all(is.finite(profile_ll)))
  
  # Maximum should be at or near the estimated value
  max_idx <- which.max(profile_ll)
  expect_true(max_idx >= 2 && max_idx <= 4,
              info = "Maximum should be near middle of grid")
})

# Likelihood Ratio Test Example ============================================

test_that("likelihood ratio test framework works", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  skip_if_not_installed("xts")
  skip_on_cran()
  
  dcc_fit <- readRDS(file.path(.get_fixtures_dir(), "dcc_dynamic_fit.rds"))
  
  # H0: alpha = 0.05, beta = 0.90 (restricted)
  ll_restricted <- compute_loglik_fixed(dcc_fit,
                                        params = list(alpha_1 = 0.05, beta_1 = 0.90))
  
  # H1: free parameters (unrestricted)
  ll_unrestricted <- as.numeric(logLik(dcc_fit))
  
  # LR statistic
  lr_stat <- 2 * (ll_unrestricted - ll_restricted)
  
  # Should be non-negative
  expect_true(lr_stat >= 0,
              info = "LR statistic should be non-negative")
  
  # Should be finite
  expect_true(is.finite(lr_stat))
})

# Error Handling Tests =====================================================

test_that("compute_loglik_fixed handles invalid parameters gracefully", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  skip_if_not_installed("xts")
  skip_on_cran()
  
  dcc_fit <- readRDS(file.path(.get_fixtures_dir(), "dcc_dynamic_fit.rds"))
  
  # Try with unknown parameter name
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
  
  # Store original parameters
  original_params <- coef(dcc_fit)
  
  # Compute with different parameters
  ll <- compute_loglik_fixed(dcc_fit, params = list(alpha_1 = 0.03, beta_1 = 0.95))
  
  # Check original not modified
  new_params <- coef(dcc_fit)
  expect_equal(original_params, new_params)
})
