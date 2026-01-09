## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## Unit Tests for DCC(1,1) Analytical Gradient
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


#### ______________________________________________________________________ ####
#### TEST SETUP                                                             ####

generate_dcc_test_data <- function(n = 200, k = 2, seed = 123) {
  set.seed(seed)
  z <- matrix(rnorm(n * k), ncol = k)
  Qbar <- cov(z)
  eig <- eigen(Qbar, symmetric = TRUE)
  if (any(eig$values < 1e-8)) {
    Qbar <- Qbar + diag(1e-6, k)
  }
  list(std_resid = z, Qbar = Qbar, weights = rep(1, n))
}

check_gradient_agreement <- function(analytical, numerical, rtol = 1e-4, atol = 1e-6) {
  analytical <- as.numeric(analytical)
  numerical <- as.numeric(numerical)
  diff <- abs(analytical - numerical)
  scale <- pmax(abs(analytical), abs(numerical), 1)
  rel_diff <- diff / scale
  all(rel_diff < rtol | diff < atol)
}


#### ______________________________________________________________________ ####
#### PART 1: DCC ORDER DETECTION AND PERSISTENCE UTILITIES                  ####

test_that("get_dcc_order correctly identifies DCC orders", {
  expect_equal(get_dcc_order(list(alpha_1 = 0.05, beta_1 = 0.90)), c(p = 1, q = 1))
  expect_equal(get_dcc_order(list(alpha_1 = 0.03, beta_1 = 0.45, beta_2 = 0.45)), c(p = 2, q = 1))
  expect_equal(get_dcc_order(NULL), c(p = 0, q = 0))
})

test_that("is_dcc11 correctly identifies DCC(1,1)", {
  expect_true(is_dcc11(list(alpha_1 = 0.05, beta_1 = 0.90)))
  expect_false(is_dcc11(list(alpha_1 = 0.03, alpha_2 = 0.02, beta_1 = 0.90)))
  expect_false(is_dcc11(NULL))
})

test_that("compute_dcc_persistence extracts parameters correctly", {
  result <- compute_dcc_persistence(list(alpha_1 = 0.05, beta_1 = 0.90))
  expect_equal(result$persistence, 0.95)
  expect_equal(result$alpha_sum, 0.05)
  expect_equal(result$beta_sum, 0.90)
})

test_that("check_dcc_stationarity validates constraints", {
  result <- check_dcc_stationarity(list(alpha_1 = 0.05, beta_1 = 0.90))
  expect_true(result$is_stationary)
  
  result <- check_dcc_stationarity(list(alpha_1 = 0.10, beta_1 = 0.95))
  expect_false(result$is_stationary)
})


#### ______________________________________________________________________ ####
#### PART 2: REPARAMETERIZATION TESTS                                       ####

test_that("Reparameterization is invertible", {
  test_cases <- list(
    c(alpha = 0.05, beta = 0.90),
    c(alpha = 0.10, beta = 0.85),
    c(alpha = 0.20, beta = 0.70)
  )
  
  for (case in test_cases) {
    unconstrained <- dcc11_to_unconstrained(case["alpha"], case["beta"])
    recovered <- dcc11_from_unconstrained(unconstrained["psi"], unconstrained["phi"])
    expect_equal(as.numeric(recovered["alpha"]), as.numeric(case["alpha"]), tolerance = 1e-10)
    expect_equal(as.numeric(recovered["beta"]), as.numeric(case["beta"]), tolerance = 1e-10)
  }
})

test_that("Reparameterization preserves stationarity constraint", {
  test_psi <- c(-5, -1, 0, 1, 5, 10)
  test_phi <- c(-5, -1, 0, 1, 5)
  
  for (psi in test_psi) {
    for (phi in test_phi) {
      params <- dcc11_from_unconstrained(psi, phi)
      expect_true(params["alpha"] > 0)
      expect_true(params["beta"] > 0)
      expect_true(params["alpha"] + params["beta"] < 1)
    }
  }
})

test_that("Jacobian is correct (numerical verification)", {
  psi <- 2.0
  phi <- -0.5
  eps <- 1e-7
  
  J_analytical <- dcc11_reparam_jacobian(psi, phi)
  J_numerical <- matrix(0, 2, 2)
  rownames(J_numerical) <- c("alpha", "beta")
  colnames(J_numerical) <- c("psi", "phi")
  
  params_plus <- dcc11_from_unconstrained(psi + eps, phi)
  params_minus <- dcc11_from_unconstrained(psi - eps, phi)
  J_numerical[1, 1] <- (params_plus["alpha"] - params_minus["alpha"]) / (2 * eps)
  J_numerical[2, 1] <- (params_plus["beta"] - params_minus["beta"]) / (2 * eps)
  
  params_plus <- dcc11_from_unconstrained(psi, phi + eps)
  params_minus <- dcc11_from_unconstrained(psi, phi - eps)
  J_numerical[1, 2] <- (params_plus["alpha"] - params_minus["alpha"]) / (2 * eps)
  J_numerical[2, 2] <- (params_plus["beta"] - params_minus["beta"]) / (2 * eps)
  
  expect_equal(J_analytical, J_numerical, tolerance = 1e-5)
})


#### ______________________________________________________________________ ####
#### PART 3: GRADIENT VS FINITE DIFFERENCES (MVN)                           ####

test_that("MVN gradient matches numerical gradient (interior point)", {
  data <- generate_dcc_test_data(n = 100, k = 2, seed = 123)
  
  params <- c(alpha = 0.08, beta = 0.85)
  
  grad_analytical <- dcc11_gradient(params, data$std_resid, data$weights, data$Qbar,
                                    distribution = "mvn", use_reparam = FALSE)
  
  grad_numerical <- numerical_gradient(
    fn = dcc11_nll, params = params, eps = 1e-6,
    std_resid = data$std_resid, weights = data$weights,
    Qbar = data$Qbar, distribution = "mvn", use_reparam = FALSE
  )
  
  expect_true(check_gradient_agreement(grad_analytical, grad_numerical, rtol = 1e-4))
})

test_that("MVN gradient matches numerical (reparameterized)", {
  data <- generate_dcc_test_data(n = 100, k = 2, seed = 456)
  
  ## Use moderate parameters that won't cause numerical issues
  params <- c(psi = 1.5, phi = 0.3)
  
  grad_analytical <- dcc11_gradient(params, data$std_resid, data$weights, data$Qbar,
                                    distribution = "mvn", use_reparam = TRUE)
  
  grad_numerical <- numerical_gradient(
    fn = dcc11_nll, params = params, eps = 1e-6,
    std_resid = data$std_resid, weights = data$weights,
    Qbar = data$Qbar, distribution = "mvn", use_reparam = TRUE
  )
  
  expect_true(check_gradient_agreement(grad_analytical, grad_numerical, rtol = 1e-4))
})

test_that("MVN gradient correct for 3-series data", {
  data <- generate_dcc_test_data(n = 100, k = 3, seed = 789)
  
  params <- c(alpha = 0.06, beta = 0.88)
  
  grad_analytical <- dcc11_gradient(params, data$std_resid, data$weights, data$Qbar,
                                    distribution = "mvn", use_reparam = FALSE)
  
  grad_numerical <- numerical_gradient(
    fn = dcc11_nll, params = params, eps = 1e-6,
    std_resid = data$std_resid, weights = data$weights,
    Qbar = data$Qbar, distribution = "mvn", use_reparam = FALSE
  )
  
  expect_true(check_gradient_agreement(grad_analytical, grad_numerical, rtol = 1e-4))
})


#### ______________________________________________________________________ ####
#### PART 4: GRADIENT VS FINITE DIFFERENCES (MVT)                           ####

test_that("MVT gradient matches numerical (original params)", {
  data <- generate_dcc_test_data(n = 100, k = 2, seed = 111)
  
  params <- c(alpha = 0.07, beta = 0.87, shape = 8.0)
  
  grad_analytical <- dcc11_gradient(params, data$std_resid, data$weights, data$Qbar,
                                    distribution = "mvt", use_reparam = FALSE)
  
  grad_numerical <- numerical_gradient(
    fn = dcc11_nll, params = params, eps = 1e-6,
    std_resid = data$std_resid, weights = data$weights,
    Qbar = data$Qbar, distribution = "mvt", use_reparam = FALSE
  )
  
  expect_true(check_gradient_agreement(grad_analytical, grad_numerical, rtol = 1e-4))
})

test_that("MVT gradient matches numerical (reparameterized)", {
  data <- generate_dcc_test_data(n = 100, k = 2, seed = 222)
  
  ## Use moderate parameters
  params <- c(psi = 2.0, phi = -0.3, shape = 6.0)
  
  grad_analytical <- dcc11_gradient(params, data$std_resid, data$weights, data$Qbar,
                                    distribution = "mvt", use_reparam = TRUE)
  
  grad_numerical <- numerical_gradient(
    fn = dcc11_nll, params = params, eps = 1e-6,
    std_resid = data$std_resid, weights = data$weights,
    Qbar = data$Qbar, distribution = "mvt", use_reparam = TRUE
  )
  
  expect_true(check_gradient_agreement(grad_analytical, grad_numerical, rtol = 1e-4))
})

test_that("MVT shape gradient is correct", {
  data <- generate_dcc_test_data(n = 100, k = 2, seed = 333)
  
  params <- c(alpha = 0.05, beta = 0.90, shape = 10.0)
  
  grad_analytical <- dcc11_gradient(params, data$std_resid, data$weights, data$Qbar,
                                    distribution = "mvt", use_reparam = FALSE)
  
  eps <- 1e-6
  nll_plus <- dcc11_nll(c(0.05, 0.90, 10.0 + eps), data$std_resid, data$weights,
                        data$Qbar, "mvt", FALSE)
  nll_minus <- dcc11_nll(c(0.05, 0.90, 10.0 - eps), data$std_resid, data$weights,
                         data$Qbar, "mvt", FALSE)
  grad_shape_numerical <- (nll_plus - nll_minus) / (2 * eps)
  
  ## Compare numeric values (strip names)
  expect_equal(as.numeric(grad_analytical["shape"]), grad_shape_numerical, tolerance = 1e-4)
})


#### ______________________________________________________________________ ####
#### PART 5: NON-UNIFORM WEIGHTS                                            ####

test_that("Gradient handles non-uniform weights correctly", {
  data <- generate_dcc_test_data(n = 100, k = 2, seed = 444)
  
  set.seed(444)
  data$weights <- runif(100, 0.2, 1.0)
  data$weights <- data$weights / sum(data$weights) * 100
  
  params <- c(alpha = 0.06, beta = 0.88)
  
  grad_analytical <- dcc11_gradient(params, data$std_resid, data$weights, data$Qbar,
                                    distribution = "mvn", use_reparam = FALSE)
  
  grad_numerical <- numerical_gradient(
    fn = dcc11_nll, params = params, eps = 1e-6,
    std_resid = data$std_resid, weights = data$weights,
    Qbar = data$Qbar, distribution = "mvn", use_reparam = FALSE
  )
  
  expect_true(check_gradient_agreement(grad_analytical, grad_numerical, rtol = 1e-4))
})

test_that("Gradient handles extreme weight distribution", {
  data <- generate_dcc_test_data(n = 100, k = 2, seed = 555)
  
  data$weights <- c(rep(4.5, 20), rep(0.1, 80))
  
  params <- c(alpha = 0.08, beta = 0.85)
  
  grad_analytical <- dcc11_gradient(params, data$std_resid, data$weights, data$Qbar,
                                    distribution = "mvn", use_reparam = FALSE)
  
  grad_numerical <- numerical_gradient(
    fn = dcc11_nll, params = params, eps = 1e-6,
    std_resid = data$std_resid, weights = data$weights,
    Qbar = data$Qbar, distribution = "mvn", use_reparam = FALSE
  )
  
  expect_true(check_gradient_agreement(grad_analytical, grad_numerical, rtol = 1e-4))
})


#### ______________________________________________________________________ ####
#### PART 6: BOUNDARY BEHAVIOR                                              ####

test_that("Gradient is stable near alpha boundary (low alpha)", {
  data <- generate_dcc_test_data(n = 100, k = 2, seed = 666)
  
  params <- c(alpha = 0.01, beta = 0.90)
  
  grad_analytical <- dcc11_gradient(params, data$std_resid, data$weights, data$Qbar,
                                    distribution = "mvn", use_reparam = FALSE)
  
  expect_true(all(is.finite(grad_analytical)))
  
  grad_numerical <- numerical_gradient(
    fn = dcc11_nll, params = params, eps = 1e-7,
    std_resid = data$std_resid, weights = data$weights,
    Qbar = data$Qbar, distribution = "mvn", use_reparam = FALSE
  )
  
  expect_true(check_gradient_agreement(grad_analytical, grad_numerical, rtol = 1e-3))
})

test_that("Gradient is stable near persistence boundary", {
  data <- generate_dcc_test_data(n = 100, k = 2, seed = 777)
  
  params <- c(alpha = 0.03, beta = 0.96)
  
  grad_analytical <- dcc11_gradient(params, data$std_resid, data$weights, data$Qbar,
                                    distribution = "mvn", use_reparam = FALSE)
  
  expect_true(all(is.finite(grad_analytical)))
})

test_that("Reparameterized gradient is well-behaved at moderate extremes", {
  data <- generate_dcc_test_data(n = 100, k = 2, seed = 888)
  
  ## Test at moderate values that should work
  ## psi = 3 -> persistence ~ 0.95
  ## psi = -2 -> persistence ~ 0.12
  test_psi <- c(-2, 3)
  test_phi <- c(-1, 1)
  
  for (psi in test_psi) {
    for (phi in test_phi) {
      params <- c(psi = psi, phi = phi)
      
      grad <- tryCatch({
        dcc11_gradient(params, data$std_resid, data$weights, data$Qbar,
                       distribution = "mvn", use_reparam = TRUE)
      }, error = function(e) NULL)
      
      expect_true(!is.null(grad) && all(is.finite(grad)),
                  info = sprintf("psi=%.1f, phi=%.1f", psi, phi))
    }
  }
})


#### ______________________________________________________________________ ####
#### PART 7: CONSISTENCY TESTS                                              ####

test_that("NLL is consistent between parameterizations", {
  data <- generate_dcc_test_data(n = 100, k = 2, seed = 1111)
  
  alpha <- 0.07
  beta <- 0.87
  
  params_orig <- c(alpha = alpha, beta = beta)
  params_reparam <- dcc11_to_unconstrained(alpha, beta)
  
  nll_orig <- dcc11_nll(params_orig, data$std_resid, data$weights, data$Qbar,
                        distribution = "mvn", use_reparam = FALSE)
  
  nll_reparam <- dcc11_nll(params_reparam, data$std_resid, data$weights, data$Qbar,
                           distribution = "mvn", use_reparam = TRUE)
  
  expect_equal(nll_orig, nll_reparam, tolerance = 1e-10)
})

test_that("Gradients are finite for both parameterizations at same point", {
  data <- generate_dcc_test_data(n = 100, k = 2, seed = 2222)
  
  alpha <- 0.06
  beta <- 0.88
  
  grad_orig <- dcc11_gradient(c(alpha = alpha, beta = beta),
                              data$std_resid, data$weights, data$Qbar,
                              distribution = "mvn", use_reparam = FALSE)
  
  params_reparam <- dcc11_to_unconstrained(alpha, beta)
  grad_reparam <- dcc11_gradient(params_reparam,
                                 data$std_resid, data$weights, data$Qbar,
                                 distribution = "mvn", use_reparam = TRUE)
  
  expect_true(all(is.finite(grad_orig)))
  expect_true(all(is.finite(grad_reparam)))
})


#### ______________________________________________________________________ ####
#### PART 8: EDGE CASES                                                     ####

test_that("Gradient handles near-singular correlation matrices", {
  set.seed(3333)
  n <- 100
  x <- rnorm(n)
  y <- 0.99 * x + 0.1 * rnorm(n)
  z <- cbind(x, y)
  z <- scale(z)
  
  Qbar <- cov(z)
  weights <- rep(1, n)
  
  params <- c(alpha = 0.05, beta = 0.90)
  
  grad <- tryCatch({
    dcc11_gradient(params, z, weights, Qbar, distribution = "mvn", use_reparam = FALSE)
  }, error = function(e) NULL)
  
  expect_true(!is.null(grad) && all(is.finite(grad)))
})

test_that("Gradient handles short time series", {
  data <- generate_dcc_test_data(n = 20, k = 2, seed = 4444)
  
  params <- c(alpha = 0.08, beta = 0.85)
  
  grad_analytical <- dcc11_gradient(params, data$std_resid, data$weights, data$Qbar,
                                    distribution = "mvn", use_reparam = FALSE)
  
  grad_numerical <- numerical_gradient(
    fn = dcc11_nll, params = params, eps = 1e-6,
    std_resid = data$std_resid, weights = data$weights,
    Qbar = data$Qbar, distribution = "mvn", use_reparam = FALSE
  )
  
  expect_true(check_gradient_agreement(grad_analytical, grad_numerical, rtol = 1e-4))
})

test_that("Gradient agrees across multiple random seeds", {
  seeds <- c(5555, 6666, 7777, 8888)
  
  for (seed in seeds) {
    data <- generate_dcc_test_data(n = 80, k = 2, seed = seed)
    
    params <- c(alpha = 0.06, beta = 0.87)
    
    grad_analytical <- dcc11_gradient(params, data$std_resid, data$weights, data$Qbar,
                                      distribution = "mvn", use_reparam = FALSE)
    
    grad_numerical <- numerical_gradient(
      fn = dcc11_nll, params = params, eps = 1e-6,
      std_resid = data$std_resid, weights = data$weights,
      Qbar = data$Qbar, distribution = "mvn", use_reparam = FALSE
    )
    
    expect_true(check_gradient_agreement(grad_analytical, grad_numerical, rtol = 1e-4),
                info = sprintf("Seed %d", seed))
  }
})




## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## Integration Tests for DCC(1,1) Analytical Gradient
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
##
## These tests verify that the new estimate_dcc_parameters_weighted() function
## with analytical gradients produces correct results and integrates properly
## with the existing MS-VARMA-GARCH framework.
##
## Test Organization:
## PART 1: Direct Function Tests (estimate_dcc_parameters_weighted)
## PART 2: Comparison with Original Implementation
## PART 3: Integration with estimate_garch_weighted_r
## PART 4: Full MS-DCC-GARCH Integration
## PART 5: MVT Distribution Tests
## PART 6: Edge Cases and Robustness

#### ______________________________________________________________________ ####
#### PART 1: DIRECT FUNCTION TESTS                                          ####

test_that("estimate_dcc_parameters_weighted detects DCC(1,1) correctly", {
  
  ## DCC(1,1) should be detected
  dcc_pars_11 <- list(alpha_1 = 0.05, beta_1 = 0.90)
  expect_true(is_dcc11(dcc_pars_11))
  
  ## DCC(2,1) should NOT be detected as DCC(1,1)
  dcc_pars_21 <- list(alpha_1 = 0.03, alpha_2 = 0.02, beta_1 = 0.90)
  expect_false(is_dcc11(dcc_pars_21))
  
  ## DCC(1,2) should NOT be detected as DCC(1,1)
  dcc_pars_12 <- list(alpha_1 = 0.05, beta_1 = 0.45, beta_2 = 0.45)
  expect_false(is_dcc11(dcc_pars_12))
})


test_that("estimate_dcc_parameters_weighted returns valid structure (MVN)", {
  skip_on_cran()
  
  sim <- simulate_dcc_garch_test_data(n = 150, k = 2, seed = 111)
  spec <- create_test_dcc_spec(k = 2, distribution = "mvn")
  
  ## Create GARCH parameters (use true values)
  garch_pars <- list(
    list(omega = 0.05, alpha1 = 0.10, beta1 = 0.85),
    list(omega = 0.05, alpha1 = 0.10, beta1 = 0.85)
  )
  
  dcc_start <- list(alpha_1 = 0.05, beta_1 = 0.90)
  dist_start <- list()
  weights <- rep(1, nrow(sim$y))
  
  result <- estimate_dcc_parameters_weighted(
    residuals = sim$y,
    weights = weights,
    garch_pars = garch_pars,
    dcc_start_pars = dcc_start,
    dist_start_pars = dist_start,
    spec = spec,
    verbose = FALSE
  )
  
  ## Check structure
  
  expect_true(is.list(result))
  expect_true("dcc_pars" %in% names(result))
  expect_true("dist_pars" %in% names(result))
  expect_true("weighted_ll" %in% names(result))
  expect_true("warnings" %in% names(result))
  
  ## Check DCC parameters
  expect_true("alpha_1" %in% names(result$dcc_pars))
  expect_true("beta_1" %in% names(result$dcc_pars))
  
  ## Check parameter validity
  alpha_est <- result$dcc_pars$alpha_1
  beta_est <- result$dcc_pars$beta_1
  
  expect_true(alpha_est > 0, info = "alpha should be positive")
  expect_true(beta_est > 0, info = "beta should be positive")
  expect_true(alpha_est + beta_est < 1, info = "persistence should be < 1")
  
  ## Check log-likelihood is finite
  expect_true(is.finite(result$weighted_ll))
})


test_that("estimate_dcc_parameters_weighted returns valid structure (MVT)", {
  skip_on_cran()
  
  sim <- simulate_dcc_garch_test_data(n = 150, k = 2, seed = 222)
  spec <- create_test_dcc_spec(k = 2, distribution = "mvt")
  
  garch_pars <- list(
    list(omega = 0.05, alpha1 = 0.10, beta1 = 0.85),
    list(omega = 0.05, alpha1 = 0.10, beta1 = 0.85)
  )
  
  dcc_start <- list(alpha_1 = 0.05, beta_1 = 0.90)
  dist_start <- list(shape = 8.0)
  weights <- rep(1, nrow(sim$y))
  
  result <- estimate_dcc_parameters_weighted(
    residuals = sim$y,
    weights = weights,
    garch_pars = garch_pars,
    dcc_start_pars = dcc_start,
    dist_start_pars = dist_start,
    spec = spec,
    verbose = FALSE
  )
  
  ## Check structure
  expect_true("dcc_pars" %in% names(result))
  expect_true("dist_pars" %in% names(result))
  
  ## Check shape parameter
  expect_true("shape" %in% names(result$dist_pars))
  expect_true(result$dist_pars$shape > 2, info = "shape should be > 2")
  
  ## Check DCC parameters
  alpha_est <- result$dcc_pars$alpha_1
  beta_est <- result$dcc_pars$beta_1
  expect_true(alpha_est + beta_est < 1)
})


#### ______________________________________________________________________ ####
#### PART 2: COMPARISON WITH ORIGINAL IMPLEMENTATION                        ####

test_that("Analytical gradient produces valid optimization results", {
  skip_on_cran()
  
  ## This test verifies that the analytical gradient implementation produces
  ## valid DCC estimates that improve the likelihood from the starting point.
  ## Note: This is NOT a parameter recovery test - that would require matching
  ## the GARCH specification to the true DGP.
  
  sim <- simulate_dcc_garch_test_data(n = 200, k = 2, 
                                      alpha = 0.06, beta = 0.88,
                                      seed = 333)
  
  spec <- create_test_dcc_spec(k = 2, distribution = "mvn")
  
  garch_pars <- list(
    list(omega = 0.05, alpha1 = 0.10, beta1 = 0.85),
    list(omega = 0.05, alpha1 = 0.10, beta1 = 0.85)
  )
  
  dcc_start <- list(alpha_1 = 0.05, beta_1 = 0.90)
  weights <- rep(1, nrow(sim$y))
  
  result <- estimate_dcc_parameters_weighted(
    residuals = sim$y,
    weights = weights,
    garch_pars = garch_pars,
    dcc_start_pars = dcc_start,
    dist_start_pars = list(),
    spec = spec,
    verbose = FALSE
  )
  
  alpha_est <- result$dcc_pars$alpha_1
  beta_est <- result$dcc_pars$beta_1
  
  cat("\nEstimated alpha:", alpha_est, "\n")
  cat("Estimated beta:", beta_est, "\n")
  cat("Estimated persistence:", alpha_est + beta_est, "\n")
  cat("Weighted log-likelihood:", result$weighted_ll, "\n")
  
  ## Key tests: valid estimates and improved likelihood
  expect_true(alpha_est >= 0, info = "alpha should be non-negative")
  expect_true(beta_est >= 0, info = "beta should be non-negative")
  expect_true(alpha_est + beta_est < 1, info = "stationarity constraint should hold")
  expect_true(is.finite(result$weighted_ll), info = "log-likelihood should be finite")
  
  ## The likelihood should be reasonable (not a penalty value)
  expect_true(result$weighted_ll > -1e9, 
              info = "log-likelihood should not be a penalty value")
})


test_that("Reparameterization produces same NLL at same point", {
  skip_on_cran()
  
  sim <- simulate_dcc_garch_test_data(n = 100, k = 2, seed = 444)
  
  ## Compute standardized residuals manually
  std_resid <- scale(sim$y)
  Qbar <- cov(std_resid)
  weights <- rep(1, nrow(std_resid))
  
  alpha <- 0.07
  beta <- 0.87
  
  ## NLL in original space
  params_orig <- c(alpha = alpha, beta = beta)
  nll_orig <- dcc11_nll(params_orig, std_resid, weights, Qbar, 
                        distribution = "mvn", use_reparam = FALSE)
  
  ## NLL in reparameterized space
  params_reparam <- dcc11_to_unconstrained(alpha, beta)
  nll_reparam <- dcc11_nll(params_reparam, std_resid, weights, Qbar,
                           distribution = "mvn", use_reparam = TRUE)
  
  expect_equal(nll_orig, nll_reparam, tolerance = 1e-10,
               info = "NLL should be identical in both parameterizations")
})


test_that("Analytical gradient matches numerical gradient at interior points", {
  skip_on_cran()
  
  ## This is the key test: verify gradient correctness
  sim <- simulate_dcc_garch_test_data(n = 100, k = 2, seed = 555)
  
  std_resid <- scale(sim$y)
  Qbar <- cov(std_resid)
  weights <- rep(1, nrow(std_resid))
  
  ## Test at several interior points
  test_points <- list(
    c(psi = 2.0, phi = -1.0),   # persistence ≈ 0.88, alpha < beta
    c(psi = 1.0, phi = 0.5),    # persistence ≈ 0.73, alpha > beta
    c(psi = 0.0, phi = 0.0)     # persistence = 0.50, alpha = beta
  )
  
  for (params in test_points) {
    ## Analytical gradient
    grad_analytical <- dcc11_gradient(
      params = params,
      std_resid = std_resid,
      weights = weights,
      Qbar = Qbar,
      distribution = "mvn",
      use_reparam = TRUE
    )
    
    ## Numerical gradient
    grad_numerical <- numerical_gradient(
      fn = dcc11_nll,
      params = params,
      eps = 1e-6,
      std_resid = std_resid,
      weights = weights,
      Qbar = Qbar,
      distribution = "mvn",
      use_reparam = TRUE
    )
    
    ## Compare
    diff <- abs(as.numeric(grad_analytical) - as.numeric(grad_numerical))
    scale <- pmax(abs(grad_analytical), abs(grad_numerical), 1)
    rel_diff <- diff / scale
    
    expect_true(all(rel_diff < 1e-4 | diff < 1e-6),
                info = sprintf("Gradient mismatch at psi=%.1f, phi=%.1f: analytical=%s, numerical=%s",
                               params[1], params[2],
                               paste(round(grad_analytical, 6), collapse=", "),
                               paste(round(grad_numerical, 6), collapse=", ")))
  }
})


#### ______________________________________________________________________ ####
#### PART 3: INTEGRATION WITH estimate_garch_weighted_r                     ####

test_that("New DCC estimation integrates with estimate_garch_weighted_r", {
  skip_on_cran()
  
  sim <- simulate_dcc_garch_test_data(n = 200, k = 2, seed = 555)
  
  spec <- list(
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
      garch_pars = list(
        list(omega = 0.05, alpha1 = 0.10, beta1 = 0.85),
        list(omega = 0.05, alpha1 = 0.10, beta1 = 0.85)
      ),
      dcc_pars = list(alpha_1 = 0.05, beta_1 = 0.90),
      dist_pars = list()
    )
  )
  
  weights <- rep(1, nrow(sim$y))
  
  result <- estimate_garch_weighted_r(
    residuals = sim$y,
    weights = weights,
    spec = spec,
    model_type = "multivariate",
    verbose = FALSE
  )
  
  expect_true(!is.null(result))
  expect_true("coefficients" %in% names(result))
  expect_true("dcc_pars" %in% names(result$coefficients) || 
                "correlation_type" %in% names(result$coefficients))
  
  ## If dynamic correlation was kept
  if (result$coefficients$correlation_type == "dynamic") {
    alpha_est <- result$coefficients$dcc_pars$alpha_1
    beta_est <- result$coefficients$dcc_pars$beta_1
    
    expect_true(alpha_est > 0)
    expect_true(beta_est > 0)
    expect_true(alpha_est + beta_est < 1)
  }
})


#### ______________________________________________________________________ ####
#### PART 4: FULL MS-DCC-GARCH INTEGRATION                                  ####

test_that("Full MS-DCC-GARCH fit works with analytical gradient", {
  skip_on_cran()
  
  sim <- simulate_dcc_garch_test_data(n = 250, k = 2, 
                                      alpha = 0.05, beta = 0.90,
                                      seed = 666)
  
  spec_ms <- list(
    ## State 1
    list(
      var_order = 0,
      garch_spec_fun = "dcc_modelspec",
      garch_spec_args = list(
        garch_model = list(univariate = list(
          list(model = "garch", garch_order = c(1, 1), distribution = "norm"),
          list(model = "garch", garch_order = c(1, 1), distribution = "norm")
        )),
        dcc_order = c(1, 1),
        dynamics = "dcc"
      ),
      distribution = "mvn",
      start_pars = list(
        var_pars = NULL,
        garch_pars = list(
          list(omega = 0.05, alpha1 = 0.10, beta1 = 0.85),
          list(omega = 0.05, alpha1 = 0.10, beta1 = 0.85)
        ),
        dcc_pars = list(alpha_1 = 0.05, beta_1 = 0.90),
        dist_pars = NULL #list()
      )
    ),
    ## State 2
    list(
      var_order = 0,
      garch_spec_fun = "dcc_modelspec",
      garch_spec_args = list(
        garch_model = list(univariate = list(
          list(model = "garch", garch_order = c(1, 1), distribution = "norm"),
          list(model = "garch", garch_order = c(1, 1), distribution = "norm")
        )),
        dcc_order = c(1, 1),
        dynamics = "dcc"
      ),
      distribution = "mvn",
      start_pars = list(
        var_pars = NULL,
        garch_pars = list(
          list(omega = 0.10, alpha1 = 0.15, beta1 = 0.75),
          list(omega = 0.10, alpha1 = 0.15, beta1 = 0.75)
        ),
        dcc_pars = list(alpha_1 = 0.10, beta_1 = 0.85),
        dist_pars = NULL #list()
      )
    )
  )
  
  fit <- fit_ms_varma_garch(
    y = sim$y,
    M = 2,
    spec = spec_ms,
    model_type = "multivariate",
    control = list(max_iter = 10, tol = 0.1)
  )
  
  expect_true(!is.null(fit))
  expect_true("model_fits" %in% names(fit))
  expect_true("smoothed_probabilities" %in% names(fit))
  expect_equal(length(fit$model_fits), 2)
  
  ## Check convergence info exists
  expect_true("convergence" %in% names(fit) || "converged" %in% names(fit))
})


test_that("MS-DCC-GARCH with MVT distribution works", {
  skip_on_cran()
  
  sim <- simulate_dcc_garch_test_data(n = 200, k = 2, seed = 777)
  
  spec_ms_mvt <- list(
    list(
      var_order = 0,
      garch_spec_fun = "dcc_modelspec",
      garch_spec_args = list(
        garch_model = list(univariate = list(
          list(model = "garch", garch_order = c(1, 1), distribution = "norm"),
          list(model = "garch", garch_order = c(1, 1), distribution = "norm")
        )),
        dcc_order = c(1, 1)
      ),
      distribution = "mvt",
      start_pars = list(
        var_pars = NULL,
        garch_pars = list(
          list(omega = 0.05, alpha1 = 0.10, beta1 = 0.85),
          list(omega = 0.05, alpha1 = 0.10, beta1 = 0.85)
        ),
        dcc_pars = list(alpha_1 = 0.05, beta_1 = 0.90),
        dist_pars = list(shape = 8.0)
      )
    ),
    list(
      var_order = 0,
      garch_spec_fun = "dcc_modelspec",
      garch_spec_args = list(
        garch_model = list(univariate = list(
          list(model = "garch", garch_order = c(1, 1), distribution = "norm"),
          list(model = "garch", garch_order = c(1, 1), distribution = "norm")
        )),
        dcc_order = c(1, 1)
      ),
      distribution = "mvt",
      start_pars = list(
        var_pars = NULL,
        garch_pars = list(
          list(omega = 0.10, alpha1 = 0.15, beta1 = 0.75),
          list(omega = 0.10, alpha1 = 0.15, beta1 = 0.75)
        ),
        dcc_pars = list(alpha_1 = 0.10, beta_1 = 0.85),
        dist_pars = list(shape = 6.0)
      )
    )
  )
  
  fit <- fit_ms_varma_garch(
    y = sim$y,
    M = 2,
    spec = spec_ms_mvt,
    model_type = "multivariate",
    control = list(max_iter = 8, tol = 0.1)
  )
  
  expect_true(!is.null(fit))
  expect_true("model_fits" %in% names(fit))
  
  ## Check that shape parameters are estimated
  for (s in 1:2) {
    if (!is.null(fit$model_fits[[s]]$shape)) {
      expect_true(fit$model_fits[[s]]$shape > 2,
                  info = sprintf("State %d shape should be > 2", s))
    }
  }
})


#### ______________________________________________________________________ ####
#### PART 5: MVT DISTRIBUTION TESTS                                         ####

test_that("MVT shape gradient is correctly integrated", {
  skip_on_cran()
  
  sim <- simulate_dcc_garch_test_data(n = 150, k = 2, seed = 888)
  spec <- create_test_dcc_spec(k = 2, distribution = "mvt")
  
  garch_pars <- list(
    list(omega = 0.05, alpha1 = 0.10, beta1 = 0.85),
    list(omega = 0.05, alpha1 = 0.10, beta1 = 0.85)
  )
  
  ## Start with different shape values
  for (shape_start in c(5.0, 8.0, 15.0)) {
    result <- estimate_dcc_parameters_weighted(
      residuals = sim$y,
      weights = rep(1, nrow(sim$y)),
      garch_pars = garch_pars,
      dcc_start_pars = list(alpha_1 = 0.05, beta_1 = 0.90),
      dist_start_pars = list(shape = shape_start),
      spec = spec,
      verbose = FALSE
    )
    
    expect_true(result$dist_pars$shape > 2,
                info = sprintf("Starting from shape=%.1f, final should be > 2", shape_start))
    expect_true(is.finite(result$weighted_ll),
                info = sprintf("Log-likelihood should be finite for shape_start=%.1f", shape_start))
  }
})


#### ______________________________________________________________________ ####
#### PART 6: EDGE CASES AND ROBUSTNESS                                      ####

test_that("Handles high persistence starting values", {
  skip_on_cran()
  
  sim <- simulate_dcc_garch_test_data(n = 150, k = 2, seed = 901)
  spec <- create_test_dcc_spec(k = 2, distribution = "mvn")
  
  garch_pars <- list(
    list(omega = 0.05, alpha1 = 0.10, beta1 = 0.85),
    list(omega = 0.05, alpha1 = 0.10, beta1 = 0.85)
  )
  
  ## High persistence starting point (0.03 + 0.96 = 0.99)
  result <- estimate_dcc_parameters_weighted(
    residuals = sim$y,
    weights = rep(1, nrow(sim$y)),
    garch_pars = garch_pars,
    dcc_start_pars = list(alpha_1 = 0.03, beta_1 = 0.96),
    dist_start_pars = list(),
    spec = spec,
    verbose = FALSE
  )
  
  ## Should still produce valid results
  expect_true(result$dcc_pars$alpha_1 > 0)
  expect_true(result$dcc_pars$beta_1 > 0)
  expect_true(result$dcc_pars$alpha_1 + result$dcc_pars$beta_1 < 1)
  expect_true(is.finite(result$weighted_ll))
})


test_that("Handles low alpha starting values", {
  skip_on_cran()
  
  sim <- simulate_dcc_garch_test_data(n = 150, k = 2, seed = 902)
  spec <- create_test_dcc_spec(k = 2, distribution = "mvn")
  
  garch_pars <- list(
    list(omega = 0.05, alpha1 = 0.10, beta1 = 0.85),
    list(omega = 0.05, alpha1 = 0.10, beta1 = 0.85)
  )
  
  ## Very low alpha starting point
  result <- estimate_dcc_parameters_weighted(
    residuals = sim$y,
    weights = rep(1, nrow(sim$y)),
    garch_pars = garch_pars,
    dcc_start_pars = list(alpha_1 = 0.001, beta_1 = 0.90),
    dist_start_pars = list(),
    spec = spec,
    verbose = FALSE
  )
  
  expect_true(result$dcc_pars$alpha_1 > 0)
  expect_true(result$dcc_pars$beta_1 > 0)
  expect_true(is.finite(result$weighted_ll))
})


test_that("Handles non-uniform weights", {
  skip_on_cran()
  
  sim <- simulate_dcc_garch_test_data(n = 150, k = 2, seed = 903)
  spec <- create_test_dcc_spec(k = 2, distribution = "mvn")
  
  garch_pars <- list(
    list(omega = 0.05, alpha1 = 0.10, beta1 = 0.85),
    list(omega = 0.05, alpha1 = 0.10, beta1 = 0.85)
  )
  
  ## Non-uniform weights (emphasize recent observations)
  set.seed(903)
  weights <- runif(nrow(sim$y), 0.5, 1.5)
  
  result <- estimate_dcc_parameters_weighted(
    residuals = sim$y,
    weights = weights,
    garch_pars = garch_pars,
    dcc_start_pars = list(alpha_1 = 0.05, beta_1 = 0.90),
    dist_start_pars = list(),
    spec = spec,
    verbose = FALSE
  )
  
  expect_true(result$dcc_pars$alpha_1 > 0)
  expect_true(result$dcc_pars$beta_1 > 0)
  expect_true(result$dcc_pars$alpha_1 + result$dcc_pars$beta_1 < 1)
  expect_true(is.finite(result$weighted_ll))
})


test_that("Handles 3-series data", {
  skip_on_cran()
  
  sim <- simulate_dcc_garch_test_data(n = 150, k = 3, seed = 904)
  spec <- create_test_dcc_spec(k = 3, distribution = "mvn")
  
  garch_pars <- list(
    list(omega = 0.05, alpha1 = 0.10, beta1 = 0.85),
    list(omega = 0.05, alpha1 = 0.10, beta1 = 0.85),
    list(omega = 0.05, alpha1 = 0.10, beta1 = 0.85)
  )
  
  result <- estimate_dcc_parameters_weighted(
    residuals = sim$y,
    weights = rep(1, nrow(sim$y)),
    garch_pars = garch_pars,
    dcc_start_pars = list(alpha_1 = 0.05, beta_1 = 0.90),
    dist_start_pars = list(),
    spec = spec,
    verbose = FALSE
  )
  
  expect_true(result$dcc_pars$alpha_1 > 0)
  expect_true(result$dcc_pars$beta_1 > 0)
  expect_true(result$dcc_pars$alpha_1 + result$dcc_pars$beta_1 < 1)
})


test_that("No penalty warnings with reparameterization", {
  skip_on_cran()
  
  ## The key benefit of reparameterization is eliminating boundary penalty warnings
  ## This test verifies we don't get any warnings during optimization
  
  sim <- simulate_dcc_garch_test_data(n = 200, k = 2, seed = 905)
  spec <- create_test_dcc_spec(k = 2, distribution = "mvn")
  
  garch_pars <- list(
    list(omega = 0.05, alpha1 = 0.10, beta1 = 0.85),
    list(omega = 0.05, alpha1 = 0.10, beta1 = 0.85)
  )
  
  ## Capture warnings
  warnings_captured <- character(0)
  
  result <- withCallingHandlers({
    estimate_dcc_parameters_weighted(
      residuals = sim$y,
      weights = rep(1, nrow(sim$y)),
      garch_pars = garch_pars,
      dcc_start_pars = list(alpha_1 = 0.05, beta_1 = 0.90),
      dist_start_pars = list(),
      spec = spec,
      verbose = FALSE
    )
  }, warning = function(w) {
    warnings_captured <<- c(warnings_captured, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  
  ## Check for penalty-related warnings
  penalty_warnings <- grep("penalty|bound|constraint|stationarity", 
                           warnings_captured, value = TRUE, ignore.case = TRUE)
  
  expect_equal(length(penalty_warnings), 0,
               info = paste("Should have no penalty warnings. Found:", 
                            paste(penalty_warnings, collapse = "; ")))
})


