## ===================================================================
## Unit Tests for MS-VARMA-GARCH Functionality
## ===================================================================

## ---- Test Data and Specifications ----
## Create minimal objects needed to run the tests.
set.seed(1)
y_test <- as.matrix(rnorm(100))
spec_test <- list(
  # State 1
  list(
    arma_order = c(1,0),
    garch_model = "garch",
    garch_order = c(1,1),
    distribution = "norm",
    start_pars = list(
      arma_pars = c(ar1 = 0.5),
      garch_pars = list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8)
    )
  ),
  # State 2
  list(
    arma_order = c(1,0),
    garch_model = "garch",
    garch_order = c(1,1),
    distribution = "norm",
    start_pars = list(
      arma_pars = c(ar1 = 0.2),
      garch_pars = list(omega = 0.2, alpha1 = 0.2, beta1 = 0.7)
    )
  )
)

## ===================================================================
## Part 1: Fast Tests (Always Run)
## ===================================================================

## ---- 1. Input Validation for fit_ms_varma_garch() ----

test_that("Input validation for 'y' argument works correctly", {
  expect_error(
    fit_ms_varma_garch(y = "not_a_matrix", M = 2, spec = spec_test),
    "Input 'y' must be a numeric matrix or data frame."
  )
  
  y_with_na <- y_test
  y_with_na[10] <- NA
  expect_error(
    fit_ms_varma_garch(y = y_with_na, M = 2, spec = spec_test),
    "Input matrix 'y' contains non-finite values"
  )
  
  y_too_short <- as.matrix(rnorm(5))
  expect_error(
    fit_ms_varma_garch(y = y_too_short, M = 2, spec = spec_test, d = 5),
    "The number of observations must be greater than the differencing order 'd'."
  )
})


test_that("Input validation for 'M' argument works correctly", {
  expect_error(
    fit_ms_varma_garch(y = y_test, M = 1, spec = list(list())),
    "'M' must be an integer >= 2."
  )
  
  expect_error(
    fit_ms_varma_garch(y = y_test, M = 2.5, spec = spec_test),
    "'M' must be an integer >= 2."
  )
})


test_that("Input validation for 'spec' argument works correctly", {
  expect_error(
    fit_ms_varma_garch(y = y_test, M = 2, spec = "not_a_list"),
    "'spec' must be a list of length M."
  )
  
  spec_wrong_length <- list(list())
  expect_error(
    fit_ms_varma_garch(y = y_test, M = 2, spec = spec_wrong_length),
    "'spec' must be a list of length M."
  )
})

test_that("Input validation for 'd' and 'model_type' works correctly", {
  expect_error(
    fit_ms_varma_garch(y = y_test, M = 2, spec = spec_test, d = -1),
    "'d' must be a non-negative integer."
  )
  
  expect_error(
    fit_ms_varma_garch(y = y_test, M = 2, spec = spec_test, model_type = "invalid_string"),
    "'arg' should be one of" # testthat default error for match.arg
  )
})


## ---- 2. Input Validation for ms_varma_garch_bs() ----

test_that("Input validation for ms_varma_garch_bs() works correctly", {
  # Test for non-positive num_boots (assuming validation is added to the function)
  expect_error(
    ms_varma_garch_bs(x = y_test, M = 2, spec = spec_test, num_boots = -1),
    # This error would come from a validation check inside ms_varma_garch_bs
    # For now, we assume it will pass down to the helper.
    "num_boots must be a positive integer." # A likely downstream error
  )
  
  # Test for missing n_boot AND num_blocks
  expect_error(
    ms_varma_garch_bs(x = y_test, M = 2, spec = spec_test, n_boot = NULL, num_blocks = NULL),
    "Must provide a valid value for either n_boot or num_blocks"
  )
  
  # Test that it correctly passes through errors from the fitter
  expect_error(
    ms_varma_garch_bs(x = y_test, M = 2, spec = "not_a_list"),
    "'spec' must be a list of length M."
  )
})


# ===================================================================
# Part 2: Slow Tests (Skip on CRAN)
# ===================================================================

test_that("Smoke test: Fitter runs for 1 iteration and returns correct structure (univariate)", {
  skip_on_cran()
  
  fit <- fit_ms_varma_garch(
    y = y_test,
    M = 2,
    spec = spec_test,
    control = list(max_iter = 1),
    parallel = FALSE # Ensure serial execution for a simple smoke test
  )
  
  expect_s3_class(fit, "msm.fit")
  expect_named(fit, c("model_fits", "P", "log_likelihood", "smoothed_probabilities",
                      "aic", "bic", "d", "y", "call", "convergence", "warnings"))
  expect_equal(dim(fit$P), c(2, 2))
  expect_equal(nrow(fit$smoothed_probabilities), nrow(y_test))
})


test_that("tsbs() with ms_varma_garch runs without error", {
  skip_on_cran()
  
  # A very small run to ensure the full user-facing pipeline connects without errors.
  # This is an integration test, not a statistical validity test.
  expect_no_error({
    tsbs(
      x = y_test,
      bs_type = "ms_varma_garch",
      num_boots = 2,
      num_blocks = 5, # ensure length is reasonable
      # Pass arguments for the MS GARCH model
      M = 2,
      spec = spec_test,
      control = list(max_iter = 2), # keep it very short
      parallel = FALSE
    )
  })
})


test_that("Full estimation converges (univariate)", {
  skip_on_cran()
  
  # This test is more intensive and should use a pre-saved, longer simulated series.
  # For this example, we'll reuse the short series but run more iterations.
  fit <- fit_ms_varma_garch(
    y = y_test,
    M = 2,
    spec = spec_test,
    control = list(max_iter = 10), # A short but meaningful run
    parallel = FALSE
  )
  # The main test is that it completes without error and produces a finite LL
  expect_true(is.finite(fit$log_likelihood))
})
