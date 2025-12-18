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
