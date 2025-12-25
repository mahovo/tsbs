## == == == == == == == == == == == == == == == == == == == == == == == == == ==
## Unit Tests for calculate_loglik_vector_r() - Multivariate Models
## == == == == == == == == == == == == == == == == == == == == == == == == == ==
## These tests verify that the refactored calculate_loglik_vector_r() function
## correctly computes per-observation log-likelihoods for multivariate models
## using the compute_loglik_fixed() approach.

## Test Setup: Create Specifications -------------------------------------------

## Univariate GARCH spec for building multivariate models
spec_uni_garch <- list(
  model = "garch", 
  garch_order = c(1, 1), 
  distribution = "norm"
)

## DCC(1,1) specification with MVN distribution
dcc_spec_args_mvn <- list(
  dcc_order = c(1, 1), 
  distribution = "mvn", 
  garch_model = list(
    univariate = list(spec_uni_garch, spec_uni_garch)
  )
)

## DCC(1,1) specification with Student-t distribution
dcc_spec_args_mvt <- list(
  dcc_order = c(1, 1), 
  distribution = "mvt", 
  garch_model = list(
    univariate = list(spec_uni_garch, spec_uni_garch)
  )
)

## DCC Constant (0,0) specification with MVN
dcc_spec_args_constant <- list(
  dcc_order = c(0, 0), 
  distribution = "mvn", 
  garch_model = list(
    univariate = list(spec_uni_garch, spec_uni_garch)
  )
)

## Copula-GARCH Dynamic specification
copula_spec_args_dynamic <- list(
  dcc_order = c(1, 1), 
  copula = "mvn", 
  garch_model = list(
    univariate = list(spec_uni_garch, spec_uni_garch)
  )
)

## Copula-GARCH Constant specification
copula_spec_args_constant <- list(
  dcc_order = c(0, 0), 
  copula = "mvn", 
  garch_model = list(
    univariate = list(spec_uni_garch, spec_uni_garch)
  )
)


## Test 1: DCC Dynamic MVN - Basic Functionality -------------------------------

test_that("calculate_loglik_vector_r works for DCC dynamic MVN", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  
  set.seed(123)
  y_sim <- matrix(rnorm(400), ncol = 2)
  colnames(y_sim) <- c("s1", "s2")
  
  spec <- list(
    var_order = 1,
    garch_spec_fun = "dcc_modelspec",
    garch_spec_args = dcc_spec_args_mvn,
    distribution = "mvn",
    start_pars = list(
      var_pars = rep(0.1, 6),
      garch_pars = list(
        list(omega = 0.01, alpha = 0.05, beta = 0.90),
        list(omega = 0.01, alpha = 0.05, beta = 0.90)
      ),
      alpha_1 = 0.05,
      beta_1 = 0.90
    )
  )
  
  current_pars <- list(
    var_pars = rep(0.1, 6),
    garch_pars = list(
      list(omega = 0.01, alpha = 0.05, beta = 0.90),
      list(omega = 0.01, alpha = 0.05, beta = 0.90)
    ),
    alpha_1 = 0.05,
    beta_1 = 0.90
  )
  
  ## Call the function
  ll_vec <- calculate_loglik_vector_r(
    y = y_sim,
    current_pars = current_pars,
    spec = spec,
    model_type = "multivariate"
  )
  
  ## Check output
  expect_type(ll_vec, "double")
  expect_length(ll_vec, nrow(y_sim))
  expect_true(all(is.finite(ll_vec)))
  expect_true(all(ll_vec <= 0))  ## Log-likelihoods should be negative or zero
})


## Test 2: DCC Dynamic MVT - With Shape Parameter ------------------------------

test_that("calculate_loglik_vector_r works for DCC dynamic MVT", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  
  set.seed(456)
  y_sim <- matrix(rt(400, df = 5), ncol = 2)
  colnames(y_sim) <- c("s1", "s2")
  
  spec <- list(
    var_order = 1,
    garch_spec_fun = "dcc_modelspec",
    garch_spec_args = dcc_spec_args_mvt,
    distribution = "mvt",
    start_pars = list(
      var_pars = rep(0.1, 6),
      garch_pars = list(
        list(omega = 0.01, alpha = 0.05, beta = 0.90),
        list(omega = 0.01, alpha = 0.05, beta = 0.90)
      ),
      dcc_alpha = 0.05,
      dcc_beta = 0.90,
      shape = 5
    )
  )
  
  current_pars <- list(
    var_pars = rep(0.1, 6),
    garch_pars = list(
      list(omega = 0.01, alpha = 0.05, beta = 0.90),
      list(omega = 0.01, alpha = 0.05, beta = 0.90)
    ),
    dcc_alpha = 0.05,
    dcc_beta = 0.90,
    shape = 5
  )
  
  ll_vec <- calculate_loglik_vector_r(
    y = y_sim,
    current_pars = current_pars,
    spec = spec,
    model_type = "multivariate"
  )
  
  expect_type(ll_vec, "double")
  expect_length(ll_vec, nrow(y_sim))
  expect_true(all(is.finite(ll_vec)))
})


## Test 3: DCC Constant - No Dynamic Parameters --------------------------------

test_that("calculate_loglik_vector_r works for DCC constant correlation", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  
  set.seed(789)
  y_sim <- matrix(rnorm(400), ncol = 2)
  colnames(y_sim) <- c("s1", "s2")
  
  spec <- list(
    var_order = 1,
    garch_spec_fun = "dcc_modelspec",
    garch_spec_args = dcc_spec_args_constant,
    distribution = "mvn",
    start_pars = list(
      var_pars = rep(0.1, 6),
      garch_pars = list(
        list(omega = 0.01, alpha = 0.05, beta = 0.90),
        list(omega = 0.01, alpha = 0.05, beta = 0.90)
      )
      ## No dcc_alpha, dcc_beta for constant correlation
    )
  )
  
  current_pars <- list(
    var_pars = rep(0.1, 6),
    garch_pars = list(
      list(omega = 0.01, alpha = 0.05, beta = 0.90),
      list(omega = 0.01, alpha = 0.05, beta = 0.90)
    )
  )
  
  ll_vec <- calculate_loglik_vector_r(
    y = y_sim,
    current_pars = current_pars,
    spec = spec,
    model_type = "multivariate"
  )
  
  expect_type(ll_vec, "double")
  expect_length(ll_vec, nrow(y_sim))
  expect_true(all(is.finite(ll_vec)))
})


## Test 4: Copula-GARCH Dynamic ------------------------------------------------

test_that("calculate_loglik_vector_r works for Copula-GARCH dynamic", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  
  set.seed(321)
  y_sim <- matrix(rnorm(400), ncol = 2)
  colnames(y_sim) <- c("s1", "s2")
  
  spec <- list(
    var_order = 1,
    garch_spec_fun = "cgarch_modelspec",
    garch_spec_args = copula_spec_args_dynamic,
    distribution = "mvn",
    start_pars = list(
      var_pars = rep(0.1, 6),
      garch_pars = list(
        list(omega = 0.01, alpha = 0.05, beta = 0.90),
        list(omega = 0.01, alpha = 0.05, beta = 0.90)
      ),
      alpha_1 = 0.05,
      beta_1 = 0.90
    )
  )
  
  current_pars <- list(
    var_pars = rep(0.1, 6),
    garch_pars = list(
      list(omega = 0.01, alpha = 0.05, beta = 0.90),
      list(omega = 0.01, alpha = 0.05, beta = 0.90)
    ),
    alpha_1 = 0.05,
    beta_1 = 0.90
  )
  
  ll_vec <- calculate_loglik_vector_r(
    y = y_sim,
    current_pars = current_pars,
    spec = spec,
    model_type = "multivariate"
  )
  
  expect_type(ll_vec, "double")
  expect_length(ll_vec, nrow(y_sim))
  expect_true(all(is.finite(ll_vec)))
})


## Test 5: Copula-GARCH Constant -----------------------------------------------

test_that("calculate_loglik_vector_r works for Copula-GARCH constant", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  
  set.seed(654)
  y_sim <- matrix(rnorm(400), ncol = 2)
  colnames(y_sim) <- c("s1", "s2")
  
  spec <- list(
    var_order = 1,
    garch_spec_fun = "cgarch_modelspec",
    garch_spec_args = copula_spec_args_constant,
    distribution = "mvn",
    start_pars = list(
      var_pars = rep(0.1, 6),
      garch_pars = list(
        list(omega = 0.01, alpha = 0.05, beta = 0.90),
        list(omega = 0.01, alpha = 0.05, beta = 0.90)
      )
    )
  )
  
  current_pars <- list(
    var_pars = rep(0.1, 6),
    garch_pars = list(
      list(omega = 0.01, alpha = 0.05, beta = 0.90),
      list(omega = 0.01, alpha = 0.05, beta = 0.90)
    )
  )
  
  ll_vec <- calculate_loglik_vector_r(
    y = y_sim,
    current_pars = current_pars,
    spec = spec,
    model_type = "multivariate"
  )
  
  expect_type(ll_vec, "double")
  expect_length(ll_vec, nrow(y_sim))
  expect_true(all(is.finite(ll_vec)))
})


## Test 6: Consistency with compute_loglik_fixed() - Direct Comparison ---------

test_that("calculate_loglik_vector_r is consistent with compute_loglik_fixed()", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  skip_if_not_installed("xts")
  
  set.seed(999)
  n <- 300
  y_sim <- matrix(rnorm(n * 2), ncol = 2)
  colnames(y_sim) <- c("s1", "s2")
  
  residuals_xts <- xts::xts(y_sim, order.by = Sys.Date() - (n:1))
  
  ## Estimate univariate GARCH models
  spec1 <- tsgarch::garch_modelspec(residuals_xts[,1], model = "garch", order = c(1,1))
  spec2 <- tsgarch::garch_modelspec(residuals_xts[,2], model = "garch", order = c(1,1))
  fit1 <- suppressWarnings(estimate(spec1, keep_tmb = TRUE))
  fit2 <- suppressWarnings(estimate(spec2, keep_tmb = TRUE))
  
  garch_fits <- tsgarch::to_multi_estimate(list(fit1, fit2))
  
  ## Estimate DCC model
  dcc_spec <- tsmarch::dcc_modelspec(garch_fits, dynamics = "dcc", dcc_order = c(1,1))
  dcc_fit <- suppressWarnings(estimate(dcc_spec))
  
  ## Reference: compute_loglik_fixed with ll_vec = TRUE
  ## Returns n - maxpq - 1 observations (burn-in + placeholder removed)
  ll_vec_reference <- compute_loglik_fixed(
    object = dcc_fit,
    params = list(),
    ll_vec = TRUE
  )
  
  maxpq <- max(dcc_spec$dynamics$order)
  
  ## compute_loglik_fixed returns n - maxpq observations (burn-in removed, no padding)
  expect_equal(length(ll_vec_reference), n - maxpq,
               info = "compute_loglik_fixed should return n - maxpq observations")
  
  ## Test calculate_loglik_vector_r with var_order = 0 (no VAR mean)
  spec_for_calc <- list(
    var_order = 0,
    garch_spec_fun = "dcc_modelspec",
    garch_spec_args = dcc_spec_args_mvn,
    distribution = "mvn"
  )
  
  current_pars <- list(
    var_pars = numeric(0),
    garch_pars = list(
      list(omega = coef(fit1)["omega"], 
           alpha = coef(fit1)["alpha1"],
           beta = coef(fit1)["beta1"]),
      list(omega = coef(fit2)["omega"], 
           alpha = coef(fit2)["alpha1"],
           beta = coef(fit2)["beta1"])
    ),
    dcc_alpha = coef(dcc_fit)["alpha_1"],
    dcc_beta = coef(dcc_fit)["beta_1"]
  )
  
  ll_vec_calc <- calculate_loglik_vector_r(
    y = y_sim,
    current_pars = current_pars,
    spec = spec_for_calc,
    model_type = "multivariate"
  )
  
  ## calculate_loglik_vector_r returns n observations (padded)
  expect_equal(length(ll_vec_calc), n,
               info = "calculate_loglik_vector_r should return n observations (with padding)")
  
  ## The actual log-likelihood values should match (after accounting for padding)
  ## ll_vec_calc has maxpq leading zeros, then n - maxpq actual values
  expect_equal(ll_vec_calc[(maxpq + 1):n], ll_vec_reference, tolerance = 1e-10,
               info = "Log-likelihood values should match after removing padding")
})


## Test 7: No Warnings About TMB Objects ---------------------------------------

test_that("calculate_loglik_vector_r produces no TMB warnings", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  
  set.seed(111)
  y_sim <- matrix(rnorm(400), ncol = 2)
  colnames(y_sim) <- c("s1", "s2")
  
  spec <- list(
    var_order = 1,
    garch_spec_fun = "dcc_modelspec",
    garch_spec_args = dcc_spec_args_mvn,
    distribution = "mvn",
    start_pars = list(
      var_pars = rep(0.1, 6),
      garch_pars = list(
        list(omega = 0.01, alpha = 0.05, beta = 0.90),
        list(omega = 0.01, alpha = 0.05, beta = 0.90)
      ),
      dcc_alpha = 0.05,
      dcc_beta = 0.90
    )
  )
  
  current_pars <- list(
    var_pars = rep(0.1, 6),
    garch_pars = list(
      list(omega = 0.01, alpha = 0.05, beta = 0.90),
      list(omega = 0.01, alpha = 0.05, beta = 0.90)
    ),
    dcc_alpha = 0.05,
    dcc_beta = 0.90
  )
  
  ## Capture warnings
  warnings_captured <- character()
  withCallingHandlers(
    {
      ll_vec <- calculate_loglik_vector_r(
        y = y_sim,
        current_pars = current_pars,
        spec = spec,
        model_type = "multivariate"
      )
    },
    warning = function(w) {
      warnings_captured <<- c(warnings_captured, conditionMessage(w))
    }
  )
  
  ## Check that there are no warnings about TMB objects
  tmb_warnings <- grep("TMB object not found", warnings_captured, value = TRUE)
  expect_equal(length(tmb_warnings), 0, 
                info = "Should not produce warnings about TMB objects")
})


## Test 8: Performance - Faster Than Old Approach ------------------------------

test_that("calculate_loglik_vector_r is reasonably fast", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  
  set.seed(222)
  y_sim <- matrix(rnorm(1000), ncol = 2)
  colnames(y_sim) <- c("s1", "s2")
  
  spec <- list(
    var_order = 1,
    garch_spec_fun = "dcc_modelspec",
    garch_spec_args = dcc_spec_args_mvn,
    distribution = "mvn",
    start_pars = list(
      var_pars = rep(0.1, 6),
      garch_pars = list(
        list(omega = 0.01, alpha = 0.05, beta = 0.90),
        list(omega = 0.01, alpha = 0.05, beta = 0.90)
      ),
      dcc_alpha = 0.05,
      dcc_beta = 0.90
    )
  )
  
  current_pars <- list(
    var_pars = rep(0.1, 6),
    garch_pars = list(
      list(omega = 0.01, alpha = 0.05, beta = 0.90),
      list(omega = 0.01, alpha = 0.05, beta = 0.90)
    ),
    dcc_alpha = 0.05,
    dcc_beta = 0.90
  )
  
  ## Time the execution
  timing <- system.time({
    ll_vec <- calculate_loglik_vector_r(
      y = y_sim,
      current_pars = current_pars,
      spec = spec,
      model_type = "multivariate"
    )
  })
  
  ## Should complete in reasonable time (adjust threshold as needed)
  expect_true(timing["elapsed"] < 5, 
              info = "Function should complete in less than 5 seconds")
})


## Test 9: Edge Case - Very Small Sample ---------------------------------------

test_that("calculate_loglik_vector_r handles small samples", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  
  set.seed(333)
  y_sim <- matrix(rnorm(100), ncol = 2)  # Only 50 observations
  colnames(y_sim) <- c("s1", "s2")
  
  spec <- list(
    var_order = 1,
    garch_spec_fun = "dcc_modelspec",
    garch_spec_args = dcc_spec_args_mvn,
    distribution = "mvn",
    start_pars = list(
      var_pars = rep(0.1, 6),
      garch_pars = list(
        list(omega = 0.01, alpha = 0.05, beta = 0.90),
        list(omega = 0.01, alpha = 0.05, beta = 0.90)
      ),
      dcc_alpha = 0.05,
      dcc_beta = 0.90
    )
  )
  
  current_pars <- list(
    var_pars = rep(0.1, 6),
    garch_pars = list(
      list(omega = 0.01, alpha = 0.05, beta = 0.90),
      list(omega = 0.01, alpha = 0.05, beta = 0.90)
    ),
    dcc_alpha = 0.05,
    dcc_beta = 0.90
  )
  
  ## Should not error on small samples
  expect_silent({
    ll_vec <- calculate_loglik_vector_r(
      y = y_sim,
      current_pars = current_pars,
      spec = spec,
      model_type = "multivariate"
    )
  })
  
  expect_type(ll_vec, "double")
  expect_length(ll_vec, nrow(y_sim))
})


## Test 10: Integration Test - Full MS-VARMA-GARCH Estimation ------------------

test_that("Full MS-VARMA-GARCH estimation converges without TMB warnings", {
  skip_if_not_installed("tsmarch")
  skip_if_not_installed("tsgarch")
  
  set.seed(444)
  y_sim_mv_short <- matrix(rnorm(400), ncol = 2)
  colnames(y_sim_mv_short) <- c("s1", "s2")
  
  spec_mv_dcc <- list(
    list(
      var_order = 1,
      garch_spec_fun = "dcc_modelspec",
      garch_spec_args = dcc_spec_args_mvn,
      distribution = "mvn",
      start_pars = list(
        var_pars = rep(0.1, 6),
        garch_pars = list(
          list(omega = 0.01, alpha = 0.05, beta = 0.90),
          list(omega = 0.01, alpha = 0.05, beta = 0.90)
        ),
        dcc_alpha = 0.1,
        dcc_beta = 0.8
      )
    ),
    list(
      var_order = 1,
      garch_spec_fun = "dcc_modelspec",
      garch_spec_args = dcc_spec_args_mvn,
      distribution = "mvn",
      start_pars = list(
        var_pars = rep(0.1, 6),
        garch_pars = list(
          list(omega = 0.02, alpha = 0.08, beta = 0.85),
          list(omega = 0.02, alpha = 0.08, beta = 0.85)
        ),
        dcc_alpha = 0.2,
        dcc_beta = 0.7
      )
    )
  )
  
  ## Capture warnings during estimation
  warnings_captured <- character()
  fit <- withCallingHandlers(
    {
      fit_ms_varma_garch(
        y = y_sim_mv_short,
        M = 2,
        spec = spec_mv_dcc,
        model_type = "multivariate",
        control = list(max_iter = 5, tol = 1e-3)
      )
    },
    warning = function(w) {
      warnings_captured <<- c(warnings_captured, conditionMessage(w))
    }
  )
  
  ## Check results
  expect_true(is.finite(fit$log_likelihood))
  
  ## Check that there are no TMB-related warnings
  tmb_warnings <- grep("TMB object not found|fallback", warnings_captured, 
                       value = TRUE, ignore.case = TRUE)
  expect_length(tmb_warnings, 0,
                info = "Full estimation should not produce TMB warnings")
})