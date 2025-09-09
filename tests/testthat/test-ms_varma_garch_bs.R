## == == == == == == == == == == == == == == == == == == == == == == ==
## Unit Tests for MS-VARMA-GARCH Functionality
## == == == == == == == == == == == == == == == == == == == == == == == 

## ---- Test Setup ----
## Create minimal data and specifications for testing.

## Univariate setup
set.seed(123)
y_test <- as.matrix(arima.sim(n = 100, list(ar = 0.5)))
colnames(y_test) <- "series_1"
spec_test_uni <- list(
  ## State 1
  list(arma_order = c(1,0),
     garch_model = "garch",
     garch_order = c(1,1),
     distribution = "norm",
     start_pars = list(
       arma_pars = c(ar1 = 0.1),
      garch_pars = list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8)
    )
  ),
  ## State 2
  list(arma_order = c(1,0),
   garch_model = "garch",
   garch_order = c(1,1),
   distribution = "norm",
   start_pars = list(
     arma_pars = c(ar1 = 0.8),
     garch_pars = list(omega = 0.2, alpha1 = 0.2, beta1 = 0.7)
    )
  )
)


## Multivariate setup
y_test_mv <- matrix(rnorm(200), ncol = 2)
colnames(y_test_mv) <- c("series_1", "series_2")

## A simple GOGARCH spec for fast smoke tests
spec_test_mv_smoke <- list(
  list(
    var_order = 1, 
    garch_spec_fun = "gogarch_modelspec", 
    garch_spec_args = list(model = "garch"), 
    start_pars = list(
      var_pars = rep(0.1, 6), 
      garch_pars = list()
    )
  ),
  list(
    var_order = 1, 
    garch_spec_fun = "gogarch_modelspec", 
    garch_spec_args = list(model = "garch"), 
    start_pars = list(
      var_pars = rep(0.1, 6), 
      garch_pars = list()
    )
  )
)

## A more complex DCC spec for slower integration tests
spec_uni_garch_dcc <- list(model = "garch", garch_order = c(1,1), distribution = "norm")

## Valid values of `distribution`: `c("mvn", "mvt")`
dcc_spec_args <- list(dcc_order = c(1,1), distribution = "mvn", garch_model = list(univariate = list(spec_uni_garch_dcc, spec_uni_garch_dcc)))
spec_mv_dcc <- list(
  list(var_order = 1, garch_spec_fun = "dcc_modelspec", garch_spec_args = dcc_spec_args, start_pars = list(var_pars=rep(0.1, 6), garch_pars=list(dcc_alpha=0.1, dcc_beta=0.8))),
  list(var_order = 1, garch_spec_fun = "dcc_modelspec", garch_spec_args = dcc_spec_args, start_pars = list(var_pars=rep(0.1, 6), garch_pars=list(dcc_alpha=0.2, dcc_beta=0.7)))
)

## == == == == == == == == == == == == == == == == == == == == == == ==
## PART 1: Fast Tests (Always Run)
## == == == == == == == == == == == == == == == == == == == == == == ==

context("MS-VARMA-GARCH: Fast Input Validation and Smoke Tests")

test_that("Input validation for fit_ms_varma_garch() works correctly", {
  ## Test 'y' argument
  expect_error(fit_ms_varma_garch(y = "not_a_matrix", M = 2, spec = spec_test_uni),
               "Input 'y' must be a numeric matrix or data frame.")
  
  y_with_na <- y_test; y_with_na[5] <- NA
  expect_error(fit_ms_varma_garch(y = y_with_na, M = 2, spec = spec_test_uni),
               "Input matrix 'y' contains non-finite values")
  
  expect_error(fit_ms_varma_garch(y = y_test[1, , drop = FALSE], M = 2, spec = spec_test_uni, d = 1),
               "The number of observations must be greater than the differencing order")
  
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
  expect_error(ms_varma_garch_bs(x = y_test, M = 2, spec = spec_test_uni, num_boots = -1),
               "num_boots must be a positive integer.")
  
  expect_error(ms_varma_garch_bs(x = y_test, M = 2, spec = spec_test_uni, n_boot = NULL, num_blocks = NULL),
               "Must provide a valid value for either n_boot or num_blocks")
  
  ## Test that it propagates errors from the underlying fitter
  expect_error(ms_varma_garch_bs(x = y_test, M = 2, spec = "not_a_list"),
               "'spec' must be a list of length M.")
})


## == == == == == == == == == == == == == == == == == == == == == == ==
## PART 2: Slow Integration Tests (Skip on CRAN)
## == == == == == == == == == == == == == == == == == == == == == == ==

context("MS-VARMA-GARCH: Slow Integration and Convergence Tests")

test_that("Smoke test: Fitter runs for 1 iteration and returns correct structure (univariate)", {
  skip_on_cran()
  
  fit <- fit_ms_varma_garch(y = y_test, M = 2, spec = spec_test_uni,
                            control = list(max_iter = 1))
  
  expect_s3_class(fit, "msm.fit")
  
  expect_named(fit, c("model_fits", "P", "log_likelihood", "smoothed_probabilities",
                      "aic", "bic", "d", "y", "call", "convergence", "warnings"))
  
  expect_equal(dim(fit$P), c(2, 2))
  
  expect_equal(nrow(fit$smoothed_probabilities), nrow(y_test))
})


test_that("Smoke test: Fitter runs for 1 iteration and returns correct structure (multivariate)", {
  skip_on_cran()
  
  fit <- fit_ms_varma_garch(y = y_test_mv, M = 2, spec = spec_test_mv_smoke, model_type = "multivariate",
                            control = list(max_iter = 1))
  
  expect_s3_class(fit, "msm.fit")
  
  expect_equal(dim(fit$P), c(2, 2))
  
  expect_equal(nrow(fit$smoothed_probabilities), nrow(y_test_mv))
})


test_that("Full estimation converges (univariate)", {
  skip_on_cran()
  
  ## A very short run to check for convergence without taking hours
  y_sim_short <- as.matrix(arima.sim(n = 200, list(ar = 0.5)))
  colnames(y_sim_short) <- "series_1"
  fit <- fit_ms_varma_garch(y = y_sim_short, M = 2, spec = spec_test_uni,
                            control = list(max_iter = 10, tol = 1e-4))
  
  expect_true(is.finite(fit$log_likelihood))
  
  ## A basic sanity check on one parameter
  ar1_est_s1 <- fit$model_fits[[1]]$arma_pars[["ar1"]]
  expect_true(ar1_est_s1 > -1 && ar1_est_s1 < 1)
})


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

test_that("Full estimation converges (multivariate)", {
  skip_on_cran()
  
  y_sim_mv_short <- matrix(rnorm(400), ncol = 2)
  colnames(y_sim_mv_short) <- c("s1", "s2")
  fit <- fit_ms_varma_garch(y = y_sim_mv_short, M = 2, spec = spec_mv_dcc, model_type = "multivariate",
                            control = list(max_iter = 5, tol = 1e-3)) ## Very short run
  
  expect_true(is.finite(fit$log_likelihood))
  
  dcc_alpha_s1 <- fit$model_fits[[1]]$garch_pars$dcc_alpha
  expect_true(dcc_alpha_s1 > 0 && dcc_alpha_s1 < 1)
})


test_that("tsbs() with ms_varma_garch runs without error (multivariate)", {
  skip_on_cran()

  ## A very small run to ensure the full user-facing pipeline connects without errors.
  ## This is an integration test, not a statistical validity test.  
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
