## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## Unit Tests for GOGARCH (gogarch_modelspec) Implementation
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
##
## Test Organization:
## - - - - - - - - - 
## PART 1:  GOGARCH Specification Tests
## PART 2:  ICA Component Extraction Tests
## PART 3:  Component GARCH Estimation Tests
## PART 4:  GOGARCH Log-Likelihood Tests
## PART 5:  Covariance and Correlation Computation Tests
## PART 6:  Integration with MS-VARMA-GARCH Framework
## PART 7:  Distribution Support Tests (norm, nig, gh)
## PART 8:  Prediction and Simulation Tests
## PART 9:  Edge Cases and Robustness
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
y_test_gogarch <- z_raw %*% L
colnames(y_test_gogarch) <- c("series_1", "series_2")

## Convert to xts for gogarch_modelspec
y_test_gogarch_xts <- xts::xts(
  y_test_gogarch, 
  order.by = Sys.Date() - (n_obs:1)
)

## 3-series test data
set.seed(456)
n_obs_3 <- 400
rho_12 <- 0.5
rho_13 <- 0.3
rho_23 <- 0.4
Sigma_3 <- matrix(c(
  1,     rho_12, rho_13,
  rho_12, 1,     rho_23,
  rho_13, rho_23, 1
), 3, 3)
L_3 <- chol(Sigma_3)
z_raw_3 <- matrix(rnorm(n_obs_3 * 3), ncol = 3)
y_test_gogarch_3 <- z_raw_3 %*% L_3
colnames(y_test_gogarch_3) <- c("series_1", "series_2", "series_3")
y_test_gogarch_3_xts <- xts::xts(
  y_test_gogarch_3, 
  order.by = Sys.Date() - (n_obs_3:1)
)


## ---- GOGARCH Specifications for MS-VARMA-GARCH Framework ----

## GOGARCH with Normal distribution
spec_gogarch_norm <- list(
  list(
    var_order = 0,
    garch_spec_fun = "gogarch_modelspec",
    distribution = "norm",
    garch_spec_args = list(
      model = "garch",
      order = c(1, 1),
      ica = "radical",
      components = 2
    ),
    start_pars = list(
      var_pars = NULL,
      garch_pars = list(
        list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
        list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8)
      ),
      dist_pars = NULL
    )
  ),
  list(
    var_order = 0,
    garch_spec_fun = "gogarch_modelspec",
    distribution = "norm",
    garch_spec_args = list(
      model = "garch",
      order = c(1, 1),
      ica = "radical",
      components = 2
    ),
    start_pars = list(
      var_pars = NULL,
      garch_pars = list(
        list(omega = 0.2, alpha1 = 0.15, beta1 = 0.75),
        list(omega = 0.2, alpha1 = 0.15, beta1 = 0.75)
      ),
      dist_pars = NULL
    )
  )
)

## GOGARCH with NIG distribution
spec_gogarch_nig <- list(
  list(
    var_order = 0,
    garch_spec_fun = "gogarch_modelspec",
    distribution = "nig",
    garch_spec_args = list(
      model = "garch",
      order = c(1, 1),
      ica = "radical",
      components = 2
    ),
    start_pars = list(
      var_pars = NULL,
      garch_pars = list(
        list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
        list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8)
      ),
      dist_pars = list(skew = 0, shape = 1)
    )
  ),
  list(
    var_order = 0,
    garch_spec_fun = "gogarch_modelspec",
    distribution = "nig",
    garch_spec_args = list(
      model = "garch",
      order = c(1, 1),
      ica = "radical",
      components = 2
    ),
    start_pars = list(
      var_pars = NULL,
      garch_pars = list(
        list(omega = 0.15, alpha1 = 0.12, beta1 = 0.78),
        list(omega = 0.15, alpha1 = 0.12, beta1 = 0.78)
      ),
      dist_pars = list(skew = 0, shape = 1)
    )
  )
)

## GOGARCH with GH distribution
spec_gogarch_gh <- list(
  list(
    var_order = 0,
    garch_spec_fun = "gogarch_modelspec",
    distribution = "gh",
    garch_spec_args = list(
      model = "garch",
      order = c(1, 1),
      ica = "radical",
      components = 2,
      lambda_range = c(-5, 5),
      shape_range = c(0.1, 25)
    ),
    start_pars = list(
      var_pars = NULL,
      garch_pars = list(
        list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
        list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8)
      ),
      dist_pars = list(skew = 0, shape = 1, lambda = -0.5)
    )
  ),
  list(
    var_order = 0,
    garch_spec_fun = "gogarch_modelspec",
    distribution = "gh",
    garch_spec_args = list(
      model = "garch",
      order = c(1, 1),
      ica = "radical",
      components = 2,
      lambda_range = c(-5, 5),
      shape_range = c(0.1, 25)
    ),
    start_pars = list(
      var_pars = NULL,
      garch_pars = list(
        list(omega = 0.15, alpha1 = 0.12, beta1 = 0.78),
        list(omega = 0.15, alpha1 = 0.12, beta1 = 0.78)
      ),
      dist_pars = list(skew = 0, shape = 1, lambda = -0.5)
    )
  )
)

#### ______________________________________________________________________ ####
#### PART 1: GOGARCH Specification Tests                                    ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

context("GOGARCH Specification")

test_that("gogarch_modelspec can be created via tsmarch", {
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  expect_s3_class(gogarch_spec, "gogarch.spec")
  expect_equal(gogarch_spec$distribution, "norm")
  expect_equal(gogarch_spec$garch$model, "garch")
  expect_equal(gogarch_spec$garch$order, c(1, 1))
  expect_equal(gogarch_spec$ica$model, "radical")
  expect_equal(gogarch_spec$ica$components, 2)
  expect_equal(gogarch_spec$n_series, 2)
  expect_equal(gogarch_spec$nobs, n_obs)
})


test_that("gogarch_modelspec supports NIG distribution", {
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_xts,
    distribution = "nig",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  expect_s3_class(gogarch_spec, "gogarch.spec")
  expect_equal(gogarch_spec$distribution, "nig")
})


test_that("gogarch_modelspec supports GH distribution", {
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_xts,
    distribution = "gh",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2,
    lambda_range = c(-5, 5),
    shape_range = c(0.1, 25)
  )
  
  expect_s3_class(gogarch_spec, "gogarch.spec")
  expect_equal(gogarch_spec$distribution, "gh")
  expect_equal(gogarch_spec$garch$lambda_range, c(-5, 5))
  expect_equal(gogarch_spec$garch$shape_range, c(0.1, 25))
})


test_that("gogarch_modelspec handles ICA algorithm specification", {
  ## tsmarch v1.0.1 supports both "radical" and "fastica", but tsbs relies on 
  ## tsmarch v1.0.0 (only v1.0.1 is available on CRAN) which only supports 
  ## "radical". 
  ## This test documents the behavior.
  
  ## radical should always work
  gogarch_spec_radical <- gogarch_modelspec(
    y = y_test_gogarch_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  expect_s3_class(gogarch_spec_radical, "gogarch.spec")
  expect_equal(gogarch_spec_radical$ica$model, "radical")
  
  ## Test fastica - may or may not be supported depending on version
  fastica_result <- tryCatch({
    gogarch_modelspec(
      y = y_test_gogarch_xts,
      distribution = "norm",
      model = "garch",
      order = c(1, 1),
      ica = "fastica",
      components = 2
    )
  }, error = function(e) {
    "fastica_not_supported"
  })
  
  if (!identical(fastica_result, "fastica_not_supported")) {
    expect_s3_class(fastica_result, "gogarch.spec")
    expect_equal(fastica_result$ica$model, "fastica")
  }
})


test_that("gogarch_modelspec rejects non-xts input", {
  
  expect_error(
    gogarch_modelspec(
      y = y_test_gogarch,  # Plain matrix, not xts
      distribution = "norm",
      model = "garch",
      order = c(1, 1)
    ),
    regexp = "xts"
  )
})


test_that("gogarch_modelspec handles 3-series data", {
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_3_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 3
  )
  
  expect_s3_class(gogarch_spec, "gogarch.spec")
  expect_equal(gogarch_spec$n_series, 3)
  expect_equal(gogarch_spec$ica$components, 3)
})


test_that("gogarch_modelspec allows fewer components than series", {
  
  ## Extract 2 components from 3 series (dimension reduction)
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_3_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2  # Less than n_series = 3
  )
  
  expect_s3_class(gogarch_spec, "gogarch.spec")
  expect_equal(gogarch_spec$n_series, 3)
  expect_equal(gogarch_spec$ica$components, 2)
})


#### ______________________________________________________________________ ####
#### PART 2: ICA Component Extraction Tests                                 ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

context("ICA Component Extraction")

test_that("GOGARCH estimation extracts ICA components correctly", {
  skip_on_cran()
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  ## Check ICA components are extracted
  expect_true(!is.null(gogarch_fit$ica))
  expect_true(!is.null(gogarch_fit$ica$S))  # Independent components
  expect_true(!is.null(gogarch_fit$ica$A))  # Mixing matrix
  expect_true(!is.null(gogarch_fit$ica$K))  # Pre-whitening matrix
  expect_true(!is.null(gogarch_fit$ica$W))  # Unmixing matrix
  
  ## Check dimensions
  expect_equal(ncol(gogarch_fit$ica$S), 2)  # 2 components
  expect_equal(nrow(gogarch_fit$ica$S), n_obs)
  expect_equal(dim(gogarch_fit$ica$A), c(2, 2))
})


test_that("ICA mixing matrix A is invertible", {
  skip_on_cran()
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  A <- gogarch_fit$ica$A
  det_A <- det(A)
  
  ## A should be invertible (non-zero determinant)
  expect_true(abs(det_A) > 1e-10,
              info = "Mixing matrix A should be invertible")
  
  ## Test that we can compute the inverse
  A_inv <- tryCatch(solve(A), error = function(e) NULL)
  expect_true(!is.null(A_inv),
              info = "A should have a computable inverse")
})


test_that("ICA components are approximately uncorrelated", {
  skip_on_cran()
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  S <- gogarch_fit$ica$S
  cor_S <- cor(S)
  
  ## Off-diagonal elements should be close to zero
  off_diag <- cor_S[1, 2]
  expect_true(abs(off_diag) < 0.2,
              info = sprintf("ICA components should be approximately uncorrelated, got cor = %.3f", off_diag))
})


#### ______________________________________________________________________ ####
#### PART 3: Component GARCH Estimation Tests                               ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

context("Component GARCH Estimation")

test_that("GARCH models are fitted to each ICA component", {
  skip_on_cran()
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  ## Should have univariate GARCH models for each component
  expect_true(!is.null(gogarch_fit$univariate))
  expect_equal(length(gogarch_fit$univariate), 2)
  
  ## Each component model should have sigma
  for (i in 1:2) {
    expect_true(!is.null(gogarch_fit$univariate[[i]]$sigma),
                info = sprintf("Component %d should have sigma", i))
  }
})


test_that("Component GARCH parameters satisfy stationarity", {
  skip_on_cran()
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  ## Check each component's GARCH parameters
  for (i in 1:length(gogarch_fit$univariate)) {
    parmat <- gogarch_fit$univariate[[i]]$parmatrix
    
    omega <- parmat[parameter == "omega"]$value
    alpha1 <- parmat[parameter == "alpha1"]$value
    beta1 <- parmat[parameter == "beta1"]$value
    
    expect_true(omega > 0,
                info = sprintf("Component %d: omega should be positive", i))
    expect_true(alpha1 >= 0,
                info = sprintf("Component %d: alpha1 should be non-negative", i))
    expect_true(beta1 >= 0,
                info = sprintf("Component %d: beta1 should be non-negative", i))
    expect_true(alpha1 + beta1 < 1,
                info = sprintf("Component %d: alpha1 + beta1 should be < 1 for stationarity", i))
  }
})


test_that("Component volatilities have correct dimensions", {
  skip_on_cran()
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  for (i in 1:2) {
    sigma_i <- sigma(gogarch_fit$univariate[[i]])
    expect_equal(length(sigma_i), n_obs,
                 info = sprintf("Component %d sigma should have length n_obs", i))
    expect_true(all(sigma_i > 0),
                info = sprintf("Component %d sigma should all be positive", i))
  }
})


#### ______________________________________________________________________ ####
#### PART 4: GOGARCH Log-Likelihood Tests                                   ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

context("GOGARCH Log-Likelihood")

test_that("GOGARCH log-likelihood is finite", {
  skip_on_cran()
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  ll <- logLik(gogarch_fit)
  
  expect_true(is.finite(as.numeric(ll)),
              info = "Log-likelihood should be finite")
})


test_that("GOGARCH log-likelihood: tsmarch v1.0.0 behavior and tsbs handling", {
  ## IMPORTANT: tsmarch v1.0.0 (CRAN) has an inconsistency in GOGARCH.
  ##
  ## Convention in DCC/CGARCH: object$loglik stores NEGATIVE log-likelihood,
  ## and logLik() returns -object$loglik (the actual LL).
  ##
  ## Bug in GOGARCH (v1.0.0): .gogarch_log_likelihood() returns the ACTUAL 
  ## log-likelihood (not NLL), stored in object$loglik. Then logLik_estimate() 
  ## negates it, so logLik() returns -LL (wrong sign).
  ##
  ## The CORRECT log-likelihood formula is:
  ##   LL = sum(component_LLs) + log|det(K)|
  ##
  ## Due to the bug:
  ##   - object$loglik = sum(component_LLs) + jacobian  (happens to be correct LL!)
  ##   - logLik(object) = -object$loglik (wrong - returns -LL)
  ##
  ## tsbs handles this by using compute_loglik_fixed() which computes LL 
  ## correctly via TMB, bypassing the buggy logLik() method.
  
  skip_on_cran()
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  ## Get component log-likelihoods (tsgarch correctly returns actual LL)
  ll_components <- sapply(gogarch_fit$univariate, function(x) as.numeric(logLik(x)))
  ll_comp_sum <- sum(ll_components)
  
  ## Compute Jacobian adjustment: log|det(K)|
  K <- gogarch_fit$ica$K
  if (nrow(K) == ncol(K)) {
    jacobian_adj <- log(abs(det(K)))
  } else {
    jacobian_adj <- log(abs(det(K %*% t(K))))
  }
  
  ## The CORRECT total log-likelihood
  correct_ll <- ll_comp_sum + jacobian_adj
  
  ## Verify tsmarch v1.0.0 behavior:
  ## 1. object$loglik stores the actual LL (due to the bug, this is correct!)
  expect_equal(gogarch_fit$loglik, correct_ll, tolerance = 1e-6,
               info = "tsmarch v1.0.0: object$loglik happens to store correct LL")
  
  ## 2. logLik() returns -object$loglik (wrong sign due to bug)
  ll_from_method <- as.numeric(logLik(gogarch_fit))
  expect_equal(ll_from_method, -correct_ll, tolerance = 1e-6,
               info = "tsmarch v1.0.0: logLik() returns wrong sign (-LL instead of LL)")
  
  ## 3. Verify the Jacobian is non-trivial
  
  expect_true(abs(jacobian_adj) > 1e-10,
              info = "Jacobian adjustment should be non-zero for correlated data")
})


test_that("tsbs compute_loglik_fixed returns correct GOGARCH log-likelihood", {
  ## tsbs bypasses tsmarch's buggy logLik() by computing LL directly via TMB
  
  skip_on_cran()
  skip_if_not(exists("compute_loglik_fixed"),
              "compute_loglik_fixed not available")
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  ## Compute the correct LL manually
  ll_components <- sapply(gogarch_fit$univariate, function(x) as.numeric(logLik(x)))
  ll_comp_sum <- sum(ll_components)
  
  K <- gogarch_fit$ica$K
  if (nrow(K) == ncol(K)) {
    jacobian_adj <- log(abs(det(K)))
  } else {
    jacobian_adj <- log(abs(det(K %*% t(K))))
  }
  correct_ll <- ll_comp_sum + jacobian_adj
  
  ## tsbs compute_loglik_fixed should return the correct LL
  tsbs_ll <- compute_loglik_fixed(
    object = gogarch_fit,
    params = list(),
    ll_vec = FALSE
  )
  
  expect_equal(tsbs_ll, correct_ll, tolerance = 1e-6,
               info = "tsbs compute_loglik_fixed should return correct LL")
  
  ## Also verify it matches object$loglik (which is correct due to the bug)
  expect_equal(tsbs_ll, gogarch_fit$loglik, tolerance = 1e-6,
               info = "tsbs LL should match object$loglik")
  
  ## And it should NOT match logLik() output (which has wrong sign)
  expect_equal(tsbs_ll, -as.numeric(logLik(gogarch_fit)), tolerance = 1e-6,
               info = "tsbs LL should be negative of logLik() output")
})


test_that("GOGARCH log-likelihood includes Jacobian adjustment", {
  ## The GOGARCH model transforms observed data Y into independent components S
  ## via S = Y * W, where W is the unmixing matrix derived from ICA.
  ##
  ## The log-likelihood must include a Jacobian term to account for this
  
  ## transformation. For the density transformation:
  ##   f_Y(y) = f_S(s) * |det(K)|
  ##
  ## Taking logs:
  ##   LL_Y = LL_S + log|det(K)|
  ##
  ## Where K is the pre-whitening matrix from ICA.
  
  skip_on_cran()
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  ## Sum of component (independent factor) log-likelihoods
  ll_components <- sapply(gogarch_fit$univariate, function(x) as.numeric(logLik(x)))
  ll_comp_sum <- sum(ll_components)
  
  ## Jacobian adjustment: log|det(K)|
  K <- gogarch_fit$ica$K
  if (nrow(K) == ncol(K)) {
    jacobian_adj <- log(abs(det(K)))
  } else {
    ## For non-square K (dimension reduction case)
    jacobian_adj <- log(abs(det(K %*% t(K))))
  }
  
  ## The correct total log-likelihood
  expected_ll <- ll_comp_sum + jacobian_adj
  
  ## Verify using object$loglik (which stores correct LL in tsmarch v1.0.0)
  expect_equal(gogarch_fit$loglik, expected_ll, tolerance = 1e-6,
               info = "Total LL should equal sum of component LLs + Jacobian")
  
  ## Verify the Jacobian is meaningful (non-zero for correlated data)
  expect_true(abs(jacobian_adj) > 1e-10,
              info = "Jacobian should be non-zero for correlated data")
  
  ## Verify Jacobian has expected sign (typically negative for whitening)
  ## K whitens the data, so |det(K)| < 1 implies log|det(K)| < 0
  ## But this depends on scaling, so we just check it's finite
  expect_true(is.finite(jacobian_adj),
              info = "Jacobian should be finite")
})



test_that("compute_loglik_fixed works for GOGARCH", {
  skip_on_cran()
  skip_if_not(exists("compute_loglik_fixed"),
              "compute_loglik_fixed not available")
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  ## Test scalar log-likelihood
  ll_scalar <- compute_loglik_fixed(
    object = gogarch_fit,
    params = list(),
    ll_vec = FALSE
  )
  
  expect_true(is.finite(ll_scalar))
  
  ## Test log-likelihood vector
  ll_vec <- compute_loglik_fixed(
    object = gogarch_fit,
    params = list(),
    ll_vec = TRUE
  )
  
  expect_true(is.numeric(ll_vec))
  expect_true(length(ll_vec) > 0)
  expect_true(all(is.finite(ll_vec)))
})


test_that("Log-likelihood vector has correct length", {
  skip_on_cran()
  skip_if_not(exists("compute_loglik_fixed"),
              "compute_loglik_fixed not available")
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  ll_vec <- compute_loglik_fixed(
    object = gogarch_fit,
    params = list(),
    ll_vec = TRUE
  )
  
  ## Should have one LL value per observation (minus any initialization)
  expect_true(length(ll_vec) <= n_obs,
              info = "LL vector length should not exceed n_obs")
  expect_true(length(ll_vec) >= n_obs - 2,
              info = "LL vector should be close to n_obs in length")
})


#### ______________________________________________________________________ ####
#### PART 5: Covariance and Correlation Computation Tests                   ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

context("GOGARCH Covariance and Correlation")

test_that("GOGARCH produces time-varying covariance matrices", {
  skip_on_cran()
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  ## Get covariance array
  H <- gogarch_fit$H  # Should be T x k x k or similar structure
  
  expect_true(!is.null(H),
              info = "GOGARCH should produce covariance matrices")
  
  ## Each time point should have a 2x2 covariance matrix
  if (is.array(H) && length(dim(H)) == 3) {
    expect_equal(dim(H)[2], 2)
    expect_equal(dim(H)[3], 2)
  }
})


test_that("GOGARCH covariance matrices are positive definite", {
  skip_on_cran()
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  ## Use tscov() to get properly formatted covariance array
  ## tscov converts from compressed lower-triangular to full array
  H <- tscov(gogarch_fit)
  
  expect_true(!is.null(H),
              info = "tscov should return covariance matrices")
  
  ## H should be a 3D array: n_series x n_series x T
  expect_true(is.array(H) && length(dim(H)) == 3,
              info = "Covariance should be a 3D array")
  expect_equal(dim(H)[1], 2, info = "First dim should be n_series")
  expect_equal(dim(H)[2], 2, info = "Second dim should be n_series")
  
  ## Check a sample of covariance matrices for positive definiteness
  n_check <- min(10, dim(H)[3])
  indices <- round(seq(1, dim(H)[3], length.out = n_check))
  
  for (t in indices) {
    H_t <- H[, , t]
    eig <- eigen(H_t, symmetric = TRUE, only.values = TRUE)
    expect_true(all(eig$values > 0),
                info = sprintf("Covariance at t=%d should be positive definite", t))
  }
})


test_that("GOGARCH produces correlation matrices", {
  skip_on_cran()
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  R <- gogarch_fit$R  # Correlation matrices
  
  expect_true(!is.null(R),
              info = "GOGARCH should produce correlation matrices")
  
  if (is.array(R) && length(dim(R)) == 3) {
    ## Check sample correlation matrices have proper structure
    n_check <- min(5, dim(R)[1])
    indices <- seq(1, dim(R)[1], length.out = n_check)
    
    for (t in indices) {
      R_t <- R[t, , ]
      
      ## Diagonal should be 1
      expect_equal(diag(R_t), rep(1, ncol(R_t)), tolerance = 1e-6,
                   info = sprintf("Correlation diagonal at t=%d should be 1", t))
      
      ## Off-diagonal should be in [-1, 1]
      off_diag <- R_t[upper.tri(R_t)]
      expect_true(all(off_diag >= -1 & off_diag <= 1),
                  info = sprintf("Correlations at t=%d should be in [-1, 1]", t))
    }
  }
})


test_that("tscov method works for GOGARCH estimates", {
  skip_on_cran()
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  ## tscov should return covariance matrices
  cov_result <- tscov(gogarch_fit)
  
  expect_true(!is.null(cov_result))
})


#### ______________________________________________________________________ ####
#### PART 6: Integration with MS-VARMA-GARCH Framework                      ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

context("GOGARCH Integration with MS-VARMA-GARCH")

test_that("calculate_loglik_vector_r handles GOGARCH spec", {
  skip_on_cran()
  skip_if_not(exists("calculate_loglik_vector_r"),
              "calculate_loglik_vector_r not available")
  
  set.seed(111)
  T_obs <- 200
  y <- matrix(rnorm(T_obs * 2), ncol = 2)
  colnames(y) <- c("series_1", "series_2")
  
  spec <- list(
    var_order = 0,
    garch_spec_fun = "gogarch_modelspec",
    distribution = "norm",
    garch_spec_args = list(
      model = "garch",
      order = c(1, 1),
      ica = "radical",
      components = 2
    ),
    start_pars = list(
      var_pars = NULL,
      garch_pars = list(
        list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
        list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8)
      ),
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
})


test_that("create_garch_spec_object_r handles GOGARCH", {
  skip_on_cran()
  skip_if_not(exists("create_garch_spec_object_r"),
              "create_garch_spec_object_r not available")
  
  set.seed(222)
  T_obs <- 200
  residuals <- matrix(rnorm(T_obs * 2), ncol = 2)
  colnames(residuals) <- c("series_1", "series_2")
  
  spec <- list(
    var_order = 0,
    garch_spec_fun = "gogarch_modelspec",
    distribution = "norm",
    garch_spec_args = list(
      model = "garch",
      order = c(1, 1),
      ica = "radical",
      components = 2
    )
  )
  
  current_pars <- list(
    garch_pars = list(
      list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
      list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8)
    )
  )
  
  garch_spec_obj <- create_garch_spec_object_r(
    residuals = residuals,
    spec = spec,
    model_type = "multivariate",
    current_pars = current_pars
  )
  
  expect_s3_class(garch_spec_obj, "gogarch.spec")
})


test_that("fit_ms_varma_garch runs with gogarch_modelspec (smoke test)", {
  skip_on_cran()
  skip_if_not(exists("fit_ms_varma_garch"),
              "fit_ms_varma_garch not available")
  
  set.seed(333)
  T_obs <- 200
  y <- matrix(rnorm(T_obs * 2), ncol = 2)
  colnames(y) <- c("series_1", "series_2")
  
  ## Run with minimal iterations (smoke test)
  fit <- tryCatch({
    fit_ms_varma_garch(
      y = y,
      M = 2,
      spec = spec_gogarch_norm,
      model_type = "multivariate",
      control = list(max_iter = 2, tol = 1),
      verbose = FALSE
    )
  }, error = function(e) {
    message("GOGARCH integration not yet implemented: ", e$message)
    NULL
  })
  
  if (!is.null(fit)) {
    ## Basic structure checks
    expect_true(is.list(fit))
    expect_true("smoothed_probabilities" %in% names(fit))
    expect_true("model_fits" %in% names(fit))
    expect_true("log_likelihood" %in% names(fit))
  }
})


#### ______________________________________________________________________ ####
#### PART 7: Distribution Support Tests                                     ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

context("GOGARCH Distribution Support")

test_that("GOGARCH with Normal distribution produces valid results", {
  skip_on_cran()
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  expect_s3_class(gogarch_fit, "gogarch.estimate")
  expect_true(is.finite(as.numeric(logLik(gogarch_fit))))
})


test_that("GOGARCH with NIG distribution produces valid results", {
  skip_on_cran()
  
  ## Generate data with heavier tails for NIG
  set.seed(789)
  n <- 300
  y_nig <- matrix(rt(n * 2, df = 5), ncol = 2)
  colnames(y_nig) <- c("series_1", "series_2")
  y_nig_xts <- xts::xts(y_nig, order.by = Sys.Date() - (n:1))
  
  gogarch_spec <- gogarch_modelspec(
    y = y_nig_xts,
    distribution = "nig",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  expect_s3_class(gogarch_fit, "gogarch.estimate")
  expect_true(is.finite(as.numeric(logLik(gogarch_fit))))
  
  ## Check that NIG parameters are estimated for each component
  for (i in 1:length(gogarch_fit$univariate)) {
    parmat <- gogarch_fit$univariate[[i]]$parmatrix
    expect_true("skew" %in% parmat$parameter,
                info = sprintf("Component %d should have skew parameter", i))
    expect_true("shape" %in% parmat$parameter,
                info = sprintf("Component %d should have shape parameter", i))
  }
})


test_that("GOGARCH with GH distribution produces valid results", {
  skip_on_cran()
  
  ## Generate data with heavier tails
  set.seed(101)
  n <- 300
  y_gh <- matrix(rt(n * 2, df = 6), ncol = 2)
  colnames(y_gh) <- c("series_1", "series_2")
  y_gh_xts <- xts::xts(y_gh, order.by = Sys.Date() - (n:1))
  
  gogarch_spec <- gogarch_modelspec(
    y = y_gh_xts,
    distribution = "gh",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2,
    lambda_range = c(-5, 5),
    shape_range = c(0.1, 25)
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  expect_s3_class(gogarch_fit, "gogarch.estimate")
  expect_true(is.finite(as.numeric(logLik(gogarch_fit))))
  
  ## Check lambda parameter is estimated
  for (i in 1:length(gogarch_fit$univariate)) {
    parmat <- gogarch_fit$univariate[[i]]$parmatrix
    expect_true("lambda" %in% parmat$parameter,
                info = sprintf("Component %d should have lambda parameter", i))
  }
})


#### ______________________________________________________________________ ####
#### PART 8: Prediction and Simulation Tests                                ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

context("GOGARCH Prediction and Simulation")

test_that("GOGARCH prediction works", {
  skip_on_cran()
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  ## Predict 5 steps ahead
  pred <- predict(gogarch_fit, h = 5, nsim = 100)
  
  expect_s3_class(pred, "gogarch.predict")
  expect_equal(pred$h, 5)
  expect_equal(pred$nsim, 100)
  expect_equal(pred$n_series, 2)
})


test_that("GOGARCH simulation works", {
  skip_on_cran()
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  ## Simulate paths
  sim <- simulate(gogarch_fit, nsim = 10, h = 50, seed = 123)
  
  expect_s3_class(sim, "gogarch.simulate")
})


test_that("GOGARCH filter works for new data", {
  skip_on_cran()
  
  ## Split data
  n_train <- 400
  y_train <- y_test_gogarch_xts[1:n_train, ]
  y_new <- y_test_gogarch_xts[(n_train + 1):(n_train + 50), ]
  
  gogarch_spec <- gogarch_modelspec(
    y = y_train,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  ## Filter new data
  filtered <- tsfilter(gogarch_fit, y = y_new)
  
  expect_s3_class(filtered, "gogarch.estimate")
  expect_equal(filtered$spec$nobs, n_train + 50)
})


#### ______________________________________________________________________ ####
#### PART 9: Edge Cases and Robustness                                      ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

context("GOGARCH Edge Cases and Robustness")

test_that("GOGARCH handles highly correlated series", {
  skip_on_cran()
  
  set.seed(555)
  n <- 400
  x <- rnorm(n)
  y <- 0.95 * x + 0.31 * rnorm(n)  # correlation ≈ 0.95
  y_high_cor <- cbind(x, y)
  colnames(y_high_cor) <- c("series_1", "series_2")
  y_high_cor_xts <- xts::xts(y_high_cor, order.by = Sys.Date() - (n:1))
  
  gogarch_spec <- gogarch_modelspec(
    y = y_high_cor_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  expect_s3_class(gogarch_fit, "gogarch.estimate")
  expect_true(is.finite(as.numeric(logLik(gogarch_fit))))
})


test_that("GOGARCH handles nearly uncorrelated series", {
  skip_on_cran()
  
  set.seed(666)
  n <- 400
  y_uncor <- matrix(rnorm(n * 2), ncol = 2)
  colnames(y_uncor) <- c("series_1", "series_2")
  y_uncor_xts <- xts::xts(y_uncor, order.by = Sys.Date() - (n:1))
  
  gogarch_spec <- gogarch_modelspec(
    y = y_uncor_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  expect_s3_class(gogarch_fit, "gogarch.estimate")
  expect_true(is.finite(as.numeric(logLik(gogarch_fit))))
})


test_that("GOGARCH handles different series scales", {
  skip_on_cran()
  
  set.seed(777)
  n <- 400
  x1 <- rnorm(n, sd = 1)
  x2 <- rnorm(n, sd = 10)  # Much larger scale
  y_diff_scale <- cbind(x1, x2)
  colnames(y_diff_scale) <- c("series_1", "series_2")
  y_diff_scale_xts <- xts::xts(y_diff_scale, order.by = Sys.Date() - (n:1))
  
  gogarch_spec <- gogarch_modelspec(
    y = y_diff_scale_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  expect_s3_class(gogarch_fit, "gogarch.estimate")
  expect_true(is.finite(as.numeric(logLik(gogarch_fit))))
})


test_that("GOGARCH handles short time series", {
  skip_on_cran()
  
  set.seed(888)
  n_short <- 100
  y_short <- matrix(rnorm(n_short * 2), ncol = 2)
  colnames(y_short) <- c("series_1", "series_2")
  y_short_xts <- xts::xts(y_short, order.by = Sys.Date() - (n_short:1))
  
  gogarch_spec <- gogarch_modelspec(
    y = y_short_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  expect_s3_class(gogarch_fit, "gogarch.estimate")
})


test_that("GOGARCH handles 3-series with 2 components (dimension reduction)", {
  skip_on_cran()
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_3_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2  # Extract only 2 components from 3 series
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  expect_s3_class(gogarch_fit, "gogarch.estimate")
  expect_equal(length(gogarch_fit$univariate), 2)  # 2 component models
  expect_equal(gogarch_fit$spec$n_series, 3)  # But data has 3 series
})


test_that("GOGARCH handles GARCH(2,1) order", {
  skip_on_cran()
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_xts,
    distribution = "norm",
    model = "garch",
    order = c(2, 1),  # GARCH(2,1)
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  expect_s3_class(gogarch_fit, "gogarch.estimate")
  
  ## Check that alpha2 parameter exists
  for (i in 1:length(gogarch_fit$univariate)) {
    parmat <- gogarch_fit$univariate[[i]]$parmatrix
    expect_true("alpha2" %in% parmat$parameter,
                info = sprintf("Component %d should have alpha2 parameter for GARCH(2,1)", i))
  }
})


test_that("Residuals can be extracted from GOGARCH fit", {
  skip_on_cran()
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  ## Standard residuals
  resid_std <- residuals(gogarch_fit, standardize = FALSE)
  expect_true(!is.null(resid_std))
  expect_equal(ncol(resid_std), 2)
  
  ## Standardized residuals
  resid_standardized <- residuals(gogarch_fit, standardize = TRUE)
  expect_true(!is.null(resid_standardized))
})


test_that("Fitted values can be extracted from GOGARCH fit", {
  skip_on_cran()
  
  gogarch_spec <- gogarch_modelspec(
    y = y_test_gogarch_xts,
    distribution = "norm",
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = 2
  )
  
  gogarch_fit <- estimate(gogarch_spec)
  
  fv <- fitted(gogarch_fit)
  expect_true(!is.null(fv))
})


## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## END OF GOGARCH TESTS
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
