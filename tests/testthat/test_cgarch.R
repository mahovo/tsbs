## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## Unit Tests for Copula GARCH (cgarch_modelspec) Implementation
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
##
## Test Organization:
## - - - - - - - - - 
## PART 1:  Copula GARCH Specification Tests
## PART 2:  PIT Transformation Tests
## PART 3:  Copula Residual Computation Tests
## PART 4:  Copula Log-Likelihood Tests
## PART 5:  Copula Parameter Estimation Tests
## PART 6:  Integration with MS-VARMA-GARCH Framework
## PART 7:  Comparison with DCC Implementation
##
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


#### ______________________________________________________________________ ####
#### TEST SETUP                                                             ####

## Required packages
# library(testthat)
# library(tsmarch)
# library(tsgarch)
# library(tsmethods)
# library(tsdistributions)
# library(xts)
# library(data.table)


## ---- Test Data Setups ----

set.seed(123)
n_obs <- 500

## Bivariate test data with correlation structure
rho_true <- 0.6
Sigma_true <- matrix(c(1, rho_true, rho_true, 1), 2, 2)
L <- chol(Sigma_true)
z_raw <- matrix(rnorm(n_obs * 2), ncol = 2)
y_test_mv_corr <- z_raw %*% L
colnames(y_test_mv_corr) <- c("series_1", "series_2")


## ---- Copula GARCH Specifications ----

## CGARCH-MVN: Gaussian Copula with parametric transformation
spec_cgarch_mvn_parametric <- list(
  list(
    var_order = 1, 
    garch_spec_fun = "cgarch_modelspec",
    distribution = "mvn",
    garch_spec_args = list(
      dcc_order = c(1, 1),
      dynamics = "dcc",
      transformation = "parametric",
      copula = "mvn",
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
      dist_pars = NULL
    )
  ),
  list(
    var_order = 1, 
    garch_spec_fun = "cgarch_modelspec",
    distribution = "mvn",
    garch_spec_args = list(
      dcc_order = c(1, 1),
      dynamics = "dcc",
      transformation = "parametric",
      copula = "mvn",
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

## CGARCH-MVT: Student-t Copula with parametric transformation
spec_cgarch_mvt_parametric <- list(
  list(
    var_order = 1,
    garch_spec_fun = "cgarch_modelspec",
    distribution = "mvt",
    garch_spec_args = list(
      dcc_order = c(1, 1),
      dynamics = "dcc",
      transformation = "parametric",
      copula = "mvt",
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
      dist_pars = list(shape = 8.0)
    )
  ),
  list(
    var_order = 1,
    garch_spec_fun = "cgarch_modelspec",
    distribution = "mvt",
    garch_spec_args = list(
      dcc_order = c(1, 1),
      dynamics = "dcc",
      transformation = "parametric",
      copula = "mvt",
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

## CGARCH with empirical transformation
spec_cgarch_mvn_empirical <- list(
  list(
    var_order = 1, 
    garch_spec_fun = "cgarch_modelspec",
    distribution = "mvn",
    garch_spec_args = list(
      dcc_order = c(1, 1),
      dynamics = "dcc",
      transformation = "empirical",
      copula = "mvn",
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
      dist_pars = NULL
    )
  )
)

## Constant correlation copula (no DCC dynamics)
spec_cgarch_constant <- list(
  list(
    var_order = 1, 
    garch_spec_fun = "cgarch_modelspec",
    distribution = "mvn",
    garch_spec_args = list(
      dcc_order = c(0, 0),
      dynamics = "constant",
      transformation = "parametric",
      copula = "mvn",
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
      dcc_pars = NULL,
      dist_pars = NULL
    )
  )
)


#### ______________________________________________________________________ ####
#### PART 1: Copula GARCH Specification Tests                               ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


test_that("cgarch_modelspec can be created via tsmarch", {
  ## Create univariate GARCH fits first
  residuals_xts <- xts::xts(y_test_mv_corr, 
                            order.by = Sys.Date() - (n_obs:1))
  
  spec1 <- tsgarch::garch_modelspec(
    y = residuals_xts[, 1],
    model = "garch",
    garch_order = c(1, 1),
    distribution = "norm"
  )
  spec2 <- tsgarch::garch_modelspec(
    y = residuals_xts[, 2],
    model = "garch",
    garch_order = c(1, 1),
    distribution = "norm"
  )
  
  fit1 <- estimate(spec1, keep_tmb = TRUE)
  fit2 <- estimate(spec2, keep_tmb = TRUE)
  
  garch_fits <- to_multi_estimate(list(fit1, fit2))
  
  ## Create cgarch specification
  cgarch_spec <- cgarch_modelspec(
    garch_fits,
    dynamics = "dcc",
    dcc_order = c(1, 1),
    transformation = "parametric",
    copula = "mvn"
  )
  
  expect_s3_class(cgarch_spec, "cgarch.spec")
  expect_equal(cgarch_spec$copula, "mvn")
  expect_equal(cgarch_spec$transformation, "parametric")
  expect_equal(cgarch_spec$dynamics$model, "dcc")
  expect_equal(cgarch_spec$dynamics$order, c(1, 1))
})


test_that("cgarch_modelspec supports mvt copula", {
  residuals_xts <- xts::xts(y_test_mv_corr, 
                            order.by = Sys.Date() - (n_obs:1))
  
  spec1 <- tsgarch::garch_modelspec(y = residuals_xts[, 1], model = "garch",
                                    garch_order = c(1, 1), distribution = "norm")
  spec2 <- tsgarch::garch_modelspec(y = residuals_xts[, 2], model = "garch",
                                    garch_order = c(1, 1), distribution = "norm")
  
  fit1 <- estimate(spec1, keep_tmb = TRUE)
  fit2 <- estimate(spec2, keep_tmb = TRUE)
  
  garch_fits <- to_multi_estimate(list(fit1, fit2))
  
  cgarch_spec <- cgarch_modelspec(
    garch_fits,
    dynamics = "dcc",
    dcc_order = c(1, 1),
    transformation = "parametric",
    copula = "mvt"
  )
  
  expect_s3_class(cgarch_spec, "cgarch.spec")
  expect_equal(cgarch_spec$copula, "mvt")
  
  ## Check that shape parameter is in parmatrix
  expect_true("shape" %in% cgarch_spec$parmatrix$parameter)
})


test_that("cgarch_modelspec supports constant dynamics", {
  residuals_xts <- xts::xts(y_test_mv_corr, 
                            order.by = Sys.Date() - (n_obs:1))
  
  spec1 <- tsgarch::garch_modelspec(y = residuals_xts[, 1], model = "garch",
                                    garch_order = c(1, 1), distribution = "norm")
  spec2 <- tsgarch::garch_modelspec(y = residuals_xts[, 2], model = "garch",
                                    garch_order = c(1, 1), distribution = "norm")
  
  fit1 <- estimate(spec1, keep_tmb = TRUE)
  fit2 <- estimate(spec2, keep_tmb = TRUE)
  
  garch_fits <- to_multi_estimate(list(fit1, fit2))
  
  cgarch_spec <- cgarch_modelspec(
    garch_fits,
    dynamics = "constant",
    transformation = "parametric",
    copula = "mvn"
  )
  
  expect_s3_class(cgarch_spec, "cgarch.spec")
  expect_equal(cgarch_spec$dynamics$model, "constant")
})


#### ______________________________________________________________________ ####
#### PART 2: PIT Transformation Tests                                       ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


test_that("compute_pit_transform produces uniform margins for parametric transform", {
  ## Create standardized normal residuals
  std_resid <- matrix(rnorm(500 * 2), ncol = 2)
  
  ## Create mock univariate fit list with normal distribution
  uni_fit_list <- list()
  for (i in 1:2) {
    uni_fit_list[[i]] <- list(
      spec = list(distribution = "norm"),
      parmatrix = data.table(
        parameter = c("skew", "shape", "lambda"),
        value = c(1, 4, 0)
      )
    )
  }
  
  ## Compute PIT
  u_matrix <- compute_pit_transform(
    std_residuals = std_resid,
    uni_fit_list = uni_fit_list,
    transformation = "parametric",
    copula_dist = "mvn"
  )
  
  ## Check output dimensions
  expect_equal(dim(u_matrix), c(500, 2))
  
  ## Check that values are in (0, 1)
  expect_true(all(u_matrix > 0 & u_matrix < 1))
  
  ## For normal distribution, the result should be approximately uniform
  ## Use KS test (p-value > 0.01 for uniformity)
  ks_result1 <- ks.test(u_matrix[, 1], "punif")
  ks_result2 <- ks.test(u_matrix[, 2], "punif")
  
  expect_true(ks_result1$p.value > 0.01)
  expect_true(ks_result2$p.value > 0.01)
})


test_that("compute_pit_transform works with empirical transformation", {
  std_resid <- matrix(rnorm(500 * 2), ncol = 2)
  
  uni_fit_list <- list(list(), list())  # Not used for empirical
  
  u_matrix <- compute_pit_transform(
    std_residuals = std_resid,
    uni_fit_list = uni_fit_list,
    transformation = "empirical",
    copula_dist = "mvn"
  )
  
  expect_equal(dim(u_matrix), c(500, 2))
  expect_true(all(u_matrix > 0 & u_matrix < 1))
  
  ## Empirical should produce exactly uniform marginals (up to discretization)
  expect_true(max(u_matrix) < 1)
  expect_true(min(u_matrix) > 0)
})


#### ______________________________________________________________________ ####
#### PART 3: Copula Residual Computation Tests                              ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


test_that("compute_copula_residuals produces standard normal for mvn copula", {
  ## Create uniform margins
  u_matrix <- matrix(runif(500 * 2), ncol = 2)
  
  z_matrix <- compute_copula_residuals(
    u_matrix = u_matrix,
    copula_dist = "mvn"
  )
  
  expect_equal(dim(z_matrix), c(500, 2))
  
  ## For Gaussian copula, Z should be approximately standard normal
  ## Mean should be close to 0
  expect_true(abs(mean(z_matrix[, 1])) < 0.15)
  expect_true(abs(mean(z_matrix[, 2])) < 0.15)
  
  ## Variance should be close to 1
  expect_true(abs(var(z_matrix[, 1]) - 1) < 0.2)
  expect_true(abs(var(z_matrix[, 2]) - 1) < 0.2)
})


test_that("compute_copula_residuals handles mvt copula correctly", {
  u_matrix <- matrix(runif(500 * 2), ncol = 2)
  
  z_matrix <- compute_copula_residuals(
    u_matrix = u_matrix,
    copula_dist = "mvt",
    dist_pars = list(shape = 5)
  )
  
  expect_equal(dim(z_matrix), c(500, 2))
  
  ## For Student-t copula with variance scaling, should have unit variance
  expect_true(abs(var(z_matrix[, 1]) - 1) < 0.3)
  expect_true(abs(var(z_matrix[, 2]) - 1) < 0.3)
})


#### ______________________________________________________________________ ####
#### PART 4: Copula Log-Likelihood Tests                                    ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


test_that("copula_nll returns finite values for valid parameters", {
  ## Create copula residuals
  z_matrix <- matrix(rnorm(500 * 2), ncol = 2)
  weights <- rep(1, 500)
  Qbar <- cov(z_matrix)
  
  ## Use reparameterized parameters
  params <- c(psi = 0, phi = 0)  # Corresponds to alpha = beta = 0.5 * (1 - pers)
  
  nll <- copula_nll(
    params = params,
    z_matrix = z_matrix,
    weights = weights,
    Qbar = Qbar,
    copula_dist = "mvn",
    use_reparam = TRUE
  )
  
  expect_true(is.finite(nll))
  expect_true(nll > 0)  # NLL should be positive
})


test_that("copula_nll handles mvt copula with shape parameter", {
  z_matrix <- matrix(rnorm(500 * 2), ncol = 2)
  weights <- rep(1, 500)
  Qbar <- cov(z_matrix)
  
  params <- c(psi = 0, phi = 0, shape = 8)
  
  nll <- copula_nll(
    params = params,
    z_matrix = z_matrix,
    weights = weights,
    Qbar = Qbar,
    copula_dist = "mvt",
    use_reparam = TRUE
  )
  
  expect_true(is.finite(nll))
})


test_that("copula_nll returns large values for invalid parameters", {
  z_matrix <- matrix(rnorm(500 * 2), ncol = 2)
  weights <- rep(1, 500)
  Qbar <- cov(z_matrix)
  
  ## Parameters that violate stationarity (persistence > 1)
  params <- c(psi = 10, phi = 0)  # Very high persistence
  
  nll <- copula_nll(
    params = params,
    z_matrix = z_matrix,
    weights = weights,
    Qbar = Qbar,
    copula_dist = "mvn",
    use_reparam = FALSE
  )
  
  expect_true(nll >= 1e10 || is.infinite(nll))
})


#### ______________________________________________________________________ ####
#### PART 5: Copula Parameter Estimation Tests                              ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


test_that("estimate_copula_parameters_weighted returns valid structure", {
  ## Test the copula parameter estimation function directly
  
  skip_if_not(exists("estimate_copula_parameters_weighted"),
              "estimate_copula_parameters_weighted not available")
  
  set.seed(456)
  T_obs <- 300
  k <- 2
  
  ## Generate test residuals
  residuals <- matrix(rnorm(T_obs * k), ncol = k)
  colnames(residuals) <- c("series_1", "series_2")
  
  ## Uniform weights
  weights <- rep(1, T_obs)
  
  ## GARCH parameters (Stage 1 results)
  garch_pars <- list(
    list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
    list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8)
  )
  
  ## DCC starting parameters
  dcc_start_pars <- list(alpha_1 = 0.05, beta_1 = 0.90)
  
  ## Spec
  spec <- list(
    var_order = 1,
    garch_spec_fun = "cgarch_modelspec",
    distribution = "mvn",
    garch_spec_args = list(
      dcc_order = c(1, 1),
      dynamics = "dcc",
      transformation = "parametric",
      copula = "mvn",
      garch_model = list(univariate = list(
        list(model = "garch", garch_order = c(1, 1), distribution = "norm"),
        list(model = "garch", garch_order = c(1, 1), distribution = "norm")
      ))
    )
  )
  
  result <- estimate_copula_parameters_weighted(
    residuals = residuals,
    weights = weights,
    garch_pars = garch_pars,
    dcc_start_pars = dcc_start_pars,
    dist_start_pars = NULL,
    spec = spec,
    transformation = "parametric",
    copula_dist = "mvn",
    verbose = FALSE
  )
  
  ## Check result structure
  expect_true(is.list(result))
  expect_true("dcc_pars" %in% names(result) || "alpha" %in% names(result))
})


test_that("estimate_garch_weighted_cgarch returns valid coefficient structure", {
  ## Test the main CGARCH estimation function
  
  skip_if_not(exists("estimate_garch_weighted_cgarch"),
              "estimate_garch_weighted_cgarch not available")
  
  ## Create test residuals
  set.seed(789)
  T_obs <- 300
  residuals <- matrix(rnorm(T_obs * 2), ncol = 2)
  colnames(residuals) <- c("series_1", "series_2")
  
  ## Uniform weights
  weights <- rep(1, T_obs)
  
  ## CGARCH specification
  spec <- list(
    var_order = 1,
    garch_spec_fun = "cgarch_modelspec",
    distribution = "mvn",
    garch_spec_args = list(
      dcc_order = c(1, 1),
      dynamics = "dcc",
      transformation = "parametric",
      copula = "mvn",
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
      dist_pars = NULL
    )
  )
  
  result <- estimate_garch_weighted_cgarch(
    residuals = residuals,
    weights = weights,
    spec = spec,
    verbose = FALSE
  )
  
  ## Check result structure
  expect_true(is.list(result))
  expect_true("coefficients" %in% names(result))
  
  coefs <- result$coefficients
  expect_true("garch_pars" %in% names(coefs))
  expect_true("correlation_type" %in% names(coefs))
  
  ## Check GARCH parameters exist for each series
  expect_equal(length(coefs$garch_pars), 2)
  
  ## Check correlation type is valid
  
  expect_true(coefs$correlation_type %in% c("constant", "dynamic"))
})


test_that("estimate_garch_weighted_cgarch handles MVT copula", {
  ## Test with Student-t copula (has shape parameter)
  
  skip_if_not(exists("estimate_garch_weighted_cgarch"),
              "estimate_garch_weighted_cgarch not available")
  
  set.seed(321)
  T_obs <- 300
  residuals <- matrix(rnorm(T_obs * 2), ncol = 2)
  colnames(residuals) <- c("series_1", "series_2")
  weights <- rep(1, T_obs)
  
  spec <- list(
    var_order = 1,
    garch_spec_fun = "cgarch_modelspec",
    distribution = "mvt",
    garch_spec_args = list(
      dcc_order = c(1, 1),
      dynamics = "dcc",
      transformation = "parametric",
      copula = "mvt",
      garch_model = list(univariate = list(
        list(model = "garch", garch_order = c(1, 1), distribution = "std"),
        list(model = "garch", garch_order = c(1, 1), distribution = "std")
      ))
    ),
    start_pars = list(
      var_pars = rep(0.1, 6),
      garch_pars = list(
        list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8, shape = 8),
        list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8, shape = 8)
      ),
      dcc_pars = list(alpha_1 = 0.05, beta_1 = 0.90),
      dist_pars = list(shape = 8)
    )
  )
  
  result <- estimate_garch_weighted_cgarch(
    residuals = residuals,
    weights = weights,
    spec = spec,
    verbose = FALSE
  )
  
  expect_true(is.list(result))
  expect_true("coefficients" %in% names(result))
  
  ## MVT should have dist_pars with shape
  if (!is.null(result$coefficients$dist_pars)) {
    expect_true("shape" %in% names(result$coefficients$dist_pars) ||
                  length(result$coefficients$dist_pars) > 0)
  }
})


test_that("estimate_garch_weighted_cgarch handles constant correlation", {
  ## Test that constant correlation model works
  
  skip_if_not(exists("estimate_garch_weighted_cgarch"),
              "estimate_garch_weighted_cgarch not available")
  
  set.seed(654)
  T_obs <- 300
  residuals <- matrix(rnorm(T_obs * 2), ncol = 2)
  colnames(residuals) <- c("series_1", "series_2")
  weights <- rep(1, T_obs)
  
  spec <- list(
    var_order = 1,
    garch_spec_fun = "cgarch_modelspec",
    distribution = "mvn",
    garch_spec_args = list(
      dcc_order = c(1, 1),
      dynamics = "constant",  ## Force constant correlation
      transformation = "parametric",
      copula = "mvn",
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
      dcc_pars = NULL,  ## No DCC parameters for constant
      dist_pars = NULL
    )
  )
  
  result <- estimate_garch_weighted_cgarch(
    residuals = residuals,
    weights = weights,
    spec = spec,
    force_constant = TRUE,
    verbose = FALSE
  )
  
  expect_true(is.list(result))
  expect_equal(result$coefficients$correlation_type, "constant")
})


#### ______________________________________________________________________ ####
#### PART 6: Integration with MS-VARMA-GARCH Framework                      ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


test_that("create_garch_spec_object_r creates valid cgarch.spec", {
  ## Test that the spec creation function handles cgarch_modelspec
  
  skip_if_not(exists("create_garch_spec_object_r"),
              "create_garch_spec_object_r not available")
  
  set.seed(111)
  T_obs <- 200
  residuals <- matrix(rnorm(T_obs * 2), ncol = 2)
  colnames(residuals) <- c("series_1", "series_2")
  
  spec <- list(
    var_order = 1,
    garch_spec_fun = "cgarch_modelspec",
    distribution = "mvn",
    garch_spec_args = list(
      dcc_order = c(1, 1),
      dynamics = "dcc",
      transformation = "parametric",
      copula = "mvn",
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
      dist_pars = NULL
    )
  )
  
  garch_spec_obj <- create_garch_spec_object_r(
    residuals = residuals,
    spec = spec,
    model_type = "multivariate",
    current_pars = spec$start_pars,
    verbose = FALSE
  )
  
  ## Check it's a valid cgarch.spec object
  expect_true(inherits(garch_spec_obj, "cgarch.spec"))
  
  ## Check it has the expected components
  expect_true(!is.null(garch_spec_obj$parmatrix))
  expect_true(!is.null(garch_spec_obj$dynamics))
  
  ## Check transformation and copula are set correctly
  expect_equal(garch_spec_obj$transformation, "parametric")
  expect_equal(garch_spec_obj$copula, "mvn")
})


test_that("calculate_loglik_vector_r computes valid log-likelihoods for cgarch", {
  ## Test log-likelihood computation with cgarch_modelspec
  
  skip_if_not(exists("calculate_loglik_vector_r"),
              "calculate_loglik_vector_r not available")
  
  set.seed(222)
  T_obs <- 200
  k <- 2
  
  ## Generate test data
  y <- matrix(rnorm(T_obs * k), ncol = k)
  colnames(y) <- c("series_1", "series_2")
  
  spec <- list(
    var_order = 1,
    garch_spec_fun = "cgarch_modelspec",
    distribution = "mvn",
    garch_spec_args = list(
      dcc_order = c(1, 1),
      dynamics = "dcc",
      transformation = "parametric",
      copula = "mvn",
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
      dist_pars = NULL
    )
  )
  
  ll_vec <- calculate_loglik_vector_r(
    y = y,
    current_pars = spec$start_pars,
    spec = spec,
    model_type = "multivariate"
  )
  
  ## Check basic properties
  expect_true(is.numeric(ll_vec))
  expect_true(length(ll_vec) > 0)
  expect_true(all(is.finite(ll_vec)))
  
  ## Log-likelihoods should be negative (for most observations)
  ## Allow some tolerance since copula LL can be positive for well-correlated data
  expect_true(mean(ll_vec) < 10)  # Reasonable upper bound
})


test_that("CGARCH routes correctly in estimate_garch_weighted_multivariate", {
  ## Verify that cgarch_modelspec is handled by the multivariate estimator
  
  skip_if_not(exists("estimate_garch_weighted_multivariate"),
              "estimate_garch_weighted_multivariate not available")
  
  set.seed(333)
  T_obs <- 300
  residuals <- matrix(rnorm(T_obs * 2), ncol = 2)
  colnames(residuals) <- c("series_1", "series_2")
  weights <- rep(1, T_obs)
  
  spec <- list(
    var_order = 1,
    garch_spec_fun = "cgarch_modelspec",
    distribution = "mvn",
    garch_spec_args = list(
      dcc_order = c(1, 1),
      dynamics = "dcc",
      transformation = "parametric",
      copula = "mvn",
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
      dist_pars = NULL
    )
  )
  
  ## This should not error - it should route to appropriate estimator
  result <- estimate_garch_weighted_multivariate(
    residuals = residuals,
    weights = weights,
    spec = spec,
    verbose = FALSE
  )
  
  expect_true(is.list(result))
  expect_true("coefficients" %in% names(result))
})


test_that("fit_ms_varma_garch runs with cgarch_modelspec (smoke test)", {
  ## Smoke test: verify fit_ms_varma_garch accepts cgarch specs
  
  skip_if_not(exists("fit_ms_varma_garch"),
              "fit_ms_varma_garch not available")
  
  set.seed(444)
  T_obs <- 200
  y <- matrix(rnorm(T_obs * 2), ncol = 2)
  colnames(y) <- c("series_1", "series_2")
  
  ## CGARCH specification for 2 states
  spec_cgarch <- list(
    ## State 1
    list(
      var_order = 1,
      garch_spec_fun = "cgarch_modelspec",
      distribution = "mvn",
      garch_spec_args = list(
        dcc_order = c(1, 1),
        dynamics = "dcc",
        transformation = "parametric",
        copula = "mvn",
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
        dist_pars = NULL
      )
    ),
    ## State 2
    list(
      var_order = 1,
      garch_spec_fun = "cgarch_modelspec",
      distribution = "mvn",
      garch_spec_args = list(
        dcc_order = c(1, 1),
        dynamics = "dcc",
        transformation = "parametric",
        copula = "mvn",
        garch_model = list(univariate = list(
          list(model = "garch", garch_order = c(1, 1), distribution = "norm"),
          list(model = "garch", garch_order = c(1, 1), distribution = "norm")
        ))
      ),
      start_pars = list(
        var_pars = rep(0.2, 6),
        garch_pars = list(
          list(omega = 0.2, alpha1 = 0.15, beta1 = 0.75),
          list(omega = 0.2, alpha1 = 0.15, beta1 = 0.75)
        ),
        dcc_pars = list(alpha_1 = 0.08, beta_1 = 0.85),
        dist_pars = NULL
      )
    )
  )
  
  ## Run with minimal iterations (smoke test)
  fit <- fit_ms_varma_garch(
    y = y,
    M = 2,
    spec = spec_cgarch,
    model_type = "multivariate",
    control = list(max_iter = 2, tol = 1),
    verbose = FALSE
  )
  
  ## Basic structure checks
  expect_true(is.list(fit))
  expect_true("smoothed_probabilities" %in% names(fit))  # State probabilities
  expect_true("model_fits" %in% names(fit))              # Parameters per state
  expect_true("log_likelihood" %in% names(fit))          # Log-likelihood
  expect_true("P" %in% names(fit))                       # Transition matrix
})


#### ______________________________________________________________________ ####
#### PART 7: Comparison with DCC Implementation                             ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


test_that("CGARCH and DCC produce consistent results for same data", {
  ## For Gaussian copula with parametric transformation,
  ## results should be similar to DCC-MVN
  
  residuals_xts <- xts::xts(y_test_mv_corr, 
                            order.by = Sys.Date() - (n_obs:1))
  
  spec1 <- tsgarch::garch_modelspec(y = residuals_xts[, 1], model = "garch",
                                    garch_order = c(1, 1), distribution = "norm")
  spec2 <- tsgarch::garch_modelspec(y = residuals_xts[, 2], model = "garch",
                                    garch_order = c(1, 1), distribution = "norm")
  
  fit1 <- estimate(spec1, keep_tmb = TRUE)
  fit2 <- estimate(spec2, keep_tmb = TRUE)
  
  garch_fits <- to_multi_estimate(list(fit1, fit2))
  
  ## Create DCC specification
  dcc_spec <- dcc_modelspec(
    garch_fits,
    dynamics = "dcc",
    dcc_order = c(1, 1),
    distribution = "mvn"
  )
  
  ## Create CGARCH specification
  cgarch_spec <- cgarch_modelspec(
    garch_fits,
    dynamics = "dcc",
    dcc_order = c(1, 1),
    transformation = "parametric",
    copula = "mvn"
  )
  
  ## Estimate both
  dcc_fit <- estimate(dcc_spec)
  cgarch_fit <- estimate(cgarch_spec)
  
  ## Compare estimated parameters (should be similar but not identical
  ## due to different likelihood formulations)
  dcc_alpha <- coef(dcc_fit)["alpha_1"]
  dcc_beta <- coef(dcc_fit)["beta_1"]
  
  cgarch_alpha <- coef(cgarch_fit)["alpha_1"]
  cgarch_beta <- coef(cgarch_fit)["beta_1"]
  
  ## Parameters should be in the same ballpark (within 0.2)
  expect_true(abs(dcc_alpha - cgarch_alpha) < 0.2)
  expect_true(abs(dcc_beta - cgarch_beta) < 0.2)
})


test_that("CGARCH captures different transformation effects", {
  ## Compare parametric vs empirical transformation
  
  residuals_xts <- xts::xts(y_test_mv_corr, 
                            order.by = Sys.Date() - (n_obs:1))
  
  spec1 <- tsgarch::garch_modelspec(y = residuals_xts[, 1], model = "garch",
                                    garch_order = c(1, 1), distribution = "norm")
  spec2 <- tsgarch::garch_modelspec(y = residuals_xts[, 2], model = "garch",
                                    garch_order = c(1, 1), distribution = "norm")
  
  fit1 <- estimate(spec1, keep_tmb = TRUE)
  fit2 <- estimate(spec2, keep_tmb = TRUE)
  
  garch_fits <- to_multi_estimate(list(fit1, fit2))
  
  ## Parametric transformation
  cgarch_param <- cgarch_modelspec(
    garch_fits,
    dynamics = "constant",
    transformation = "parametric",
    copula = "mvn"
  )
  
  ## Empirical transformation
  cgarch_empir <- cgarch_modelspec(
    garch_fits,
    dynamics = "constant",
    transformation = "empirical",
    copula = "mvn"
  )
  
  fit_param <- estimate(cgarch_param)
  fit_empir <- estimate(cgarch_empir)
  
  ## Both should produce valid estimates
  expect_s3_class(fit_param, c("cgarch.estimate", "cgarch.constant"))
  expect_s3_class(fit_empir, c("cgarch.estimate", "cgarch.constant"))
  
  ## Log-likelihoods should be finite
  expect_true(is.finite(as.numeric(logLik(fit_param))))
  expect_true(is.finite(as.numeric(logLik(fit_empir))))
})


## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## END OF COPULA GARCH TESTS
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =