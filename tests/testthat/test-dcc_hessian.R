## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## Unit Tests for DCC(1,1) Hessian and Standard Errors
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
##
## Prerequisites:
##   - source("dcc_gradient.R")
##   - source("dcc_hessian.R")
##
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


#### ______________________________________________________________________ ####
#### HELPER FUNCTIONS                                                       ####

#' Generate test data with true DCC dynamics for Hessian tests
#' @param n Number of observations
#' @param k Number of series
#' @param alpha_true True DCC alpha parameter
#' @param beta_true True DCC beta parameter
#' @param seed Random seed
generate_hessian_test_data <- function(n = 300, k = 2, 
                                       alpha_true = 0.05, beta_true = 0.90,
                                       seed = 123) {
  set.seed(seed)
  
  ## Unconditional correlation matrix
  Qbar <- diag(k)
  if (k >= 2) {
    for (i in 1:(k-1)) {
      for (j in (i+1):k) {
        Qbar[i, j] <- Qbar[j, i] <- 0.3  ## Lower unconditional correlation
      }
    }
  }
  
  ## Initialize arrays
  z <- matrix(0, n, k)
  Q <- array(0, dim = c(k, k, n))
  R <- array(0, dim = c(k, k, n))
  
  ## Initial values
  Q[,,1] <- Qbar
  R[,,1] <- Qbar
  
  ## First observation - draw from unconditional
  L1 <- chol(Qbar)
  z[1, ] <- as.vector(t(L1) %*% rnorm(k))
  
  ## Simulate DCC process
  for (t in 2:n) {
    ## Update Q using lagged z
    z_lag <- z[t-1, , drop = FALSE]
    Q[,,t] <- (1 - alpha_true - beta_true) * Qbar + 
      alpha_true * (t(z_lag) %*% z_lag) + 
      beta_true * Q[,,t-1]
    
    ## Normalize to correlation
    d_t <- sqrt(diag(Q[,,t]))
    d_t[d_t < 1e-8] <- 1e-8
    D_inv <- diag(1 / d_t, k)
    R[,,t] <- D_inv %*% Q[,,t] %*% D_inv
    
    ## Ensure R_t is a valid correlation matrix
    R_t <- R[,,t]
    ## Force symmetry
    R_t <- (R_t + t(R_t)) / 2
    ## Ensure PD
    eig <- eigen(R_t, symmetric = TRUE)
    if (any(eig$values < 1e-8)) {
      eig$values[eig$values < 1e-8] <- 1e-8
      R_t <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
      ## Re-normalize to have unit diagonal
      d_t <- sqrt(diag(R_t))
      D_inv <- diag(1 / d_t, k)
      R_t <- D_inv %*% R_t %*% D_inv
    }
    R[,,t] <- R_t
    
    ## Generate correlated innovations from R_t
    L <- chol(R_t)
    z[t, ] <- as.vector(t(L) %*% rnorm(k))
  }
  
  ## Compute sample correlation as Qbar estimate
  Qbar_est <- cor(z)
  
  ## Diagnostic: check that correlations actually vary
  if (k == 2) {
    rolling_cor <- sapply(30:n, function(t) cor(z[(t-29):t, 1], z[(t-29):t, 2]))
    cor_range <- range(rolling_cor)
  }
  
  list(
    std_resid = z,
    Qbar = Qbar_est,
    weights = rep(1, n),
    true_alpha = alpha_true,
    true_beta = beta_true,
    true_persistence = alpha_true + beta_true
  )
}


#### ______________________________________________________________________ ####
#### PART 1: NUMERICAL HESSIAN TESTS                                        ####

test_that("numerical_hessian computes correct Hessian for simple function", {
  ## Test on a simple quadratic function where we know the Hessian
  ## f(x, y) = x^2 + 3*y^2 + 2*x*y
  ## H = [[2, 2], [2, 6]]
  ##
  ## The error in second-order finite differences is O(ε²) for the truncation
  ## error, but there's also O(1/ε²) roundoff error, so the optimal ε is around
  ## 1e-4 to 1e-5, giving errors around 1e-8 to 1e-6 in ideal cases.
  
  f <- function(params) {
    x <- params[1]
    y <- params[2]
    x^2 + 3*y^2 + 2*x*y
  }
  
  H <- numerical_hessian(f, c(1, 1), eps = 1e-5)
  
  expected_H <- matrix(c(2, 2, 2, 6), 2, 2)
  
  ## Numerical differentiation has O(eps^2) error, so tolerance ~1e-8 for eps=1e-5
  ## But finite precision effects can make it slightly worse
  expect_equal(H, expected_H, tolerance = 1e-4)
})


test_that("numerical_hessian_richardson improves accuracy", {
  ## Rosenbrock function: f(x,y) = (1-x)^2 + 100*(y-x^2)^2
  ## More complex, tests Richardson extrapolation
  
  rosenbrock <- function(params) {
    x <- params[1]
    y <- params[2]
    (1 - x)^2 + 100 * (y - x^2)^2
  }
  
  ## At (1, 1), the minimum:
  ## H = [[802, -400], [-400, 200]]
  
  H_standard <- numerical_hessian(rosenbrock, c(1, 1), eps = 1e-4)
  H_richardson <- numerical_hessian_richardson(rosenbrock, c(1, 1), eps = 1e-4)
  
  expected_H <- matrix(c(802, -400, -400, 200), 2, 2)
  
  ## Richardson should be closer to true value
  error_standard <- max(abs(H_standard - expected_H))
  error_richardson <- max(abs(H_richardson - expected_H))
  
  expect_lt(error_richardson, error_standard)
  expect_lt(error_richardson, 1e-4)
})


test_that("numerical_hessian is symmetric", {
  data <- generate_hessian_test_data(n = 100, seed = 111)
  
  params <- c(alpha = 0.07, beta = 0.85)
  
  H <- numerical_hessian(
    fn = dcc11_nll,
    params = params,
    eps = 1e-5,
    std_resid = data$std_resid,
    weights = data$weights,
    Qbar = data$Qbar,
    distribution = "mvn",
    use_reparam = FALSE
  )
  
  expect_equal(H, t(H), tolerance = 1e-10)
})


#### ______________________________________________________________________ ####
#### PART 2: OBSERVED INFORMATION MATRIX                                    ####

test_that("observed information is positive definite at interior point", {
  data <- generate_hessian_test_data(n = 200, seed = 222)
  
  params <- c(alpha = 0.06, beta = 0.88)
  
  I_obs <- dcc11_observed_information(
    params = params,
    std_resid = data$std_resid,
    weights = data$weights,
    Qbar = data$Qbar,
    distribution = "mvn",
    use_reparam = FALSE
  )
  
  eig <- eigen(I_obs, symmetric = TRUE, only.values = TRUE)$values
  
  ## At a reasonable interior point, information should be positive definite
  expect_true(all(eig > 0), 
              info = sprintf("Eigenvalues should be positive: %s", 
                             paste(round(eig, 4), collapse = ", ")))
})


test_that("observed information works in reparameterized space", {
  data <- generate_hessian_test_data(n = 200, seed = 333)
  
  ## Transform to (psi, phi) space
  params_reparam <- dcc11_to_unconstrained(0.06, 0.88)
  
  I_obs <- dcc11_observed_information(
    params = params_reparam,
    std_resid = data$std_resid,
    weights = data$weights,
    Qbar = data$Qbar,
    distribution = "mvn",
    use_reparam = TRUE
  )
  
  expect_true(is.matrix(I_obs))
  expect_equal(dim(I_obs), c(2, 2))
  expect_true(all(is.finite(I_obs)))
})


#### ______________________________________________________________________ ####
#### PART 3: STANDARD ERRORS                                                ####

test_that("dcc11_standard_errors returns valid structure", {
  data <- generate_hessian_test_data(n = 200, seed = 444)
  
  ## First, find the actual MLE for this data
  ## (the arbitrary params c(0.06, 0.88) may not be at the minimum)
  start_params <- c(alpha = 0.06, beta = 0.88)
  
  opt_result <- optim(
    par = start_params,
    fn = dcc11_nll,
    method = "L-BFGS-B",
    lower = c(1e-6, 1e-6),
    upper = c(0.5, 0.99),
    std_resid = data$std_resid,
    weights = data$weights,
    Qbar = data$Qbar,
    distribution = "mvn",
    use_reparam = FALSE
  )
  
  ## Use the MLE for SE computation
  params_mle <- opt_result$par
  names(params_mle) <- c("alpha", "beta")
  
  se_result <- dcc11_standard_errors(
    params = params_mle,
    std_resid = data$std_resid,
    weights = data$weights,
    Qbar = data$Qbar,
    distribution = "mvn",
    use_reparam = FALSE
  )
  
  expect_true(is.list(se_result))
  expect_true("se" %in% names(se_result))
  expect_true("vcov" %in% names(se_result))
  expect_true("info" %in% names(se_result))
  
  expect_equal(length(se_result$se), 2)
  expect_equal(dim(se_result$vcov), c(2, 2))
  
  ## At the MLE, standard errors should be positive (information matrix should be PD)
  expect_true(all(se_result$se > 0), 
              info = sprintf("SEs should be positive at MLE: %s",
                             paste(round(se_result$se, 6), collapse = ", ")))
})


test_that("standard errors are reasonable magnitude", {
  ## Use data with clear DCC dynamics
  ## Higher alpha (0.08) makes dynamics more identifiable
  data <- generate_hessian_test_data(n = 500, k = 2, 
                                     alpha_true = 0.08, beta_true = 0.85,
                                     seed = 555)
  
  ## Find MLE - should be close to true values with proper DCC data
  start_params <- c(alpha = 0.08, beta = 0.85)
  
  opt_result <- optim(
    par = start_params,
    fn = dcc11_nll,
    method = "L-BFGS-B",
    lower = c(1e-4, 1e-4),  ## Slightly higher lower bound
    upper = c(0.5, 0.99),
    std_resid = data$std_resid,
    weights = data$weights,
    Qbar = data$Qbar,
    distribution = "mvn",
    use_reparam = FALSE
  )
  
  params_mle <- opt_result$par
  names(params_mle) <- c("alpha", "beta")
  
  cat("\nTrue params: alpha =", data$true_alpha, ", beta =", data$true_beta, "\n")
  cat("MLE params:  alpha =", round(params_mle["alpha"], 4), 
      ", beta =", round(params_mle["beta"], 4), "\n")
  cat("Convergence:", opt_result$convergence, "\n")
  
  se_result <- dcc11_standard_errors(
    params = params_mle,
    std_resid = data$std_resid,
    weights = data$weights,
    Qbar = data$Qbar,
    distribution = "mvn",
    use_reparam = FALSE
  )
  
  cat("SEs:         alpha =", round(se_result$se["alpha"], 4), 
      ", beta =", round(se_result$se["beta"], 4), "\n")
  
  ## Check that MLE is not at boundary
  expect_gt(params_mle["alpha"], 0.01, 
            label = "alpha MLE should not be at lower boundary")
  expect_gt(params_mle["beta"], 0.01, 
            label = "beta MLE should not be at lower boundary")
  
  ## With proper DCC data (n=500), SEs should be reasonably small
  expect_lt(se_result$se["alpha"], 0.10, 
            label = "alpha SE should be < 0.10 with n=500 DCC data")
  expect_lt(se_result$se["beta"], 0.20, 
            label = "beta SE should be < 0.20 with n=500 DCC data")
  
  ## But not extremely small (would indicate numerical issues)
  expect_gt(se_result$se["alpha"], 1e-4)
  expect_gt(se_result$se["beta"], 1e-4)
})


test_that("standard errors decrease with sample size", {
  ## Asymptotic property: SE should scale as 1/sqrt(n)
  ## Generate one long dataset and use subsets
  
  data_full <- generate_hessian_test_data(n = 600, k = 2,
                                          alpha_true = 0.08, beta_true = 0.85,
                                          seed = 666)
  
  ## Helper to find MLE and compute SE for first n observations
  get_se_for_n <- function(n, data) {
    z <- data$std_resid[1:n, , drop = FALSE]
    Qbar <- cor(z)
    weights <- rep(1, n)
    
    opt_result <- optim(
      par = c(0.08, 0.85),
      fn = dcc11_nll,
      method = "L-BFGS-B",
      lower = c(1e-4, 1e-4),
      upper = c(0.5, 0.99),
      std_resid = z,
      weights = weights,
      Qbar = Qbar,
      distribution = "mvn",
      use_reparam = FALSE
    )
    
    params_mle <- opt_result$par
    names(params_mle) <- c("alpha", "beta")
    
    se <- dcc11_standard_errors(
      params = params_mle,
      std_resid = z,
      weights = weights,
      Qbar = Qbar,
      distribution = "mvn",
      use_reparam = FALSE
    )$se
    
    list(params = params_mle, se = se)
  }
  
  result_n200 <- get_se_for_n(200, data_full)
  result_n500 <- get_se_for_n(500, data_full)
  
  cat("\nn=200: MLE =", round(result_n200$params, 4), 
      ", SE =", round(result_n200$se, 4), "\n")
  cat("n=500: MLE =", round(result_n500$params, 4), 
      ", SE =", round(result_n500$se, 4), "\n")
  
  ## Only compare if both MLEs are in the interior (not at boundary)
  ## Otherwise the comparison is meaningless
  alpha_interior_200 <- result_n200$params["alpha"] > 0.01
  
  alpha_interior_500 <- result_n500$params["alpha"] > 0.01
  beta_interior_200 <- result_n200$params["beta"] > 0.01
  beta_interior_500 <- result_n500$params["beta"] > 0.01
  
  if (alpha_interior_200 && alpha_interior_500) {
    ## SE at n=500 should be smaller than at n=200 for alpha
    expect_lt(result_n500$se["alpha"], result_n200$se["alpha"],
              label = "alpha SE should decrease with sample size")
  } else {
    cat("Skipping alpha SE comparison - MLE at boundary\n")
  }
  
  if (beta_interior_200 && beta_interior_500) {
    ## SE at n=500 should be smaller than at n=200 for beta
    expect_lt(result_n500$se["beta"], result_n200$se["beta"],
              label = "beta SE should decrease with sample size")
  } else {
    cat("Skipping beta SE comparison - MLE at boundary\n")
  }
  
  ## At minimum, SEs should be finite and positive
  expect_true(all(is.finite(result_n200$se)))
  expect_true(all(is.finite(result_n500$se)))
  expect_true(all(result_n200$se > 0))
  expect_true(all(result_n500$se > 0))
})


test_that("MVT standard errors include shape parameter", {
  ## Use properly simulated DCC data
  data <- generate_hessian_test_data(n = 300, k = 2,
                                     alpha_true = 0.08, beta_true = 0.85,
                                     seed = 777)
  
  ## Find MLE for MVT model
  start_params <- c(alpha = 0.08, beta = 0.85, shape = 8.0)
  
  opt_result <- optim(
    par = start_params,
    fn = dcc11_nll,
    method = "L-BFGS-B",
    lower = c(1e-4, 1e-4, 2.5),
    upper = c(0.5, 0.99, 50),
    std_resid = data$std_resid,
    weights = data$weights,
    Qbar = data$Qbar,
    distribution = "mvt",
    use_reparam = FALSE
  )
  
  params_mle <- opt_result$par
  names(params_mle) <- c("alpha", "beta", "shape")
  
  cat("\nMVT MLE: alpha =", round(params_mle["alpha"], 4),
      ", beta =", round(params_mle["beta"], 4),
      ", shape =", round(params_mle["shape"], 2), "\n")
  
  se_result <- dcc11_standard_errors(
    params = params_mle,
    std_resid = data$std_resid,
    weights = data$weights,
    Qbar = data$Qbar,
    distribution = "mvt",
    use_reparam = FALSE
  )
  
  cat("SEs: alpha =", round(se_result$se["alpha"], 4),
      ", beta =", round(se_result$se["beta"], 4),
      ", shape =", round(se_result$se["shape"], 2), "\n")
  
  expect_equal(length(se_result$se), 3)
  expect_true("shape" %in% names(se_result$se))
  
  ## Shape SE should be positive if MLE is in interior
  if (params_mle["shape"] > 3 && params_mle["shape"] < 49) {
    expect_true(se_result$se["shape"] > 0,
                info = sprintf("Shape SE should be > 0, got %g", se_result$se["shape"]))
  }
  
  ## At minimum, should have valid structure
  expect_equal(dim(se_result$vcov), c(3, 3))
})


#### ______________________________________________________________________ ####
#### PART 4: TRANSFORMATION BETWEEN PARAMETERIZATIONS                       ####

test_that("SE transformation via delta method is consistent", {
  ## This test verifies that the delta method transformation between
  ## parameterizations gives approximately consistent results.
  ##
  ## Strategy: Test the transformation at a known interior point.
  ##
  ## Note: Some discrepancy is expected because we're comparing:
  ## 1. SEs from numerical Hessian in (alpha, beta) space
  ## 2. SEs from numerical Hessian in (psi, phi) space, then transformed
  ## The numerical differentiation has different step size characteristics
  ## in each parameterization.
  
  ## Generate data with strong DCC dynamics
  data <- generate_hessian_test_data(n = 500, k = 2,
                                     alpha_true = 0.15, beta_true = 0.75,
                                     seed = 8888)
  
  ## Use a known interior point (not necessarily MLE)
  params_orig <- c(alpha = 0.10, beta = 0.80)
  
  cat("\nTest point: alpha =", params_orig["alpha"],
      ", beta =", params_orig["beta"], "\n")
  
  ## Compute SE in original space at this point
  se_orig <- dcc11_standard_errors(
    params = params_orig,
    std_resid = data$std_resid,
    weights = data$weights,
    Qbar = data$Qbar,
    distribution = "mvn",
    use_reparam = FALSE
  )
  
  cat("Original SEs: alpha =", round(se_orig$se["alpha"], 5),
      ", beta =", round(se_orig$se["beta"], 5), "\n")
  
  ## Transform to reparameterized space
  params_reparam <- dcc11_to_unconstrained(params_orig["alpha"], params_orig["beta"])
  cat("Reparam point: psi =", round(params_reparam["psi"], 4),
      ", phi =", round(params_reparam["phi"], 4), "\n")
  
  ## Compute SE in reparameterized space at the corresponding point
  se_reparam <- dcc11_standard_errors(
    params = params_reparam,
    std_resid = data$std_resid,
    weights = data$weights,
    Qbar = data$Qbar,
    distribution = "mvn",
    use_reparam = TRUE
  )
  
  cat("Reparam SEs: psi =", round(se_reparam$se["psi"], 5),
      ", phi =", round(se_reparam$se["phi"], 5), "\n")
  
  ## Check that both information matrices are positive definite
  eig_orig <- eigen(se_orig$info, symmetric = TRUE, only.values = TRUE)$values
  eig_reparam <- eigen(se_reparam$info, symmetric = TRUE, only.values = TRUE)$values
  
  cat("Original info eigenvalues:", round(eig_orig, 2), "\n")
  cat("Reparam info eigenvalues:", round(eig_reparam, 2), "\n")
  
  ## If either information matrix is not PD, the test is not meaningful
  if (any(eig_orig <= 0)) {
    skip("Original information matrix not positive definite")
  }
  if (any(eig_reparam <= 0)) {
    skip("Reparameterized information matrix not positive definite")
  }
  
  ## Transform reparam SE to original space using delta method
  se_transformed <- dcc11_transform_se(se_reparam, to_reparam = FALSE)
  
  cat("Transformed SEs: alpha =", round(se_transformed$se["alpha"], 5),
      ", beta =", round(se_transformed$se["beta"], 5), "\n")
  
  ## Compute relative differences
  rel_diff_alpha <- abs(se_transformed$se["alpha"] - se_orig$se["alpha"]) / se_orig$se["alpha"]
  rel_diff_beta <- abs(se_transformed$se["beta"] - se_orig$se["beta"]) / se_orig$se["beta"]
  
  cat("Relative differences: alpha =", round(rel_diff_alpha, 3),
      ", beta =", round(rel_diff_beta, 3), "\n")
  
  ## The transformed SE should be within 50% of directly computed SE
  ## This is a loose tolerance because numerical Hessians computed in
  ## different parameterizations can differ due to step size effects
  expect_lt(rel_diff_alpha, 0.50,
            label = "Alpha SE relative difference should be < 50%")
  expect_lt(rel_diff_beta, 0.50,
            label = "Beta SE relative difference should be < 50%")
  
  ## Also verify that both are positive and finite
  expect_true(all(is.finite(se_transformed$se)))
  expect_true(all(se_transformed$se > 0))
})


#### ______________________________________________________________________ ####
#### PART 5: CONFIDENCE INTERVALS                                           ####

test_that("confidence intervals have correct structure", {
  data <- generate_hessian_test_data(n = 200, seed = 999)
  
  ## Find MLE first
  start_params <- c(alpha = 0.06, beta = 0.88)
  
  opt_result <- optim(
    par = start_params,
    fn = dcc11_nll,
    method = "L-BFGS-B",
    lower = c(1e-6, 1e-6),
    upper = c(0.5, 0.99),
    std_resid = data$std_resid,
    weights = data$weights,
    Qbar = data$Qbar,
    distribution = "mvn",
    use_reparam = FALSE
  )
  
  params_mle <- opt_result$par
  names(params_mle) <- c("alpha", "beta")
  
  se_result <- dcc11_standard_errors(
    params = params_mle,
    std_resid = data$std_resid,
    weights = data$weights,
    Qbar = data$Qbar,
    distribution = "mvn",
    use_reparam = FALSE
  )
  
  ci <- dcc11_confint(se_result, level = 0.95)
  
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 2)
  expect_equal(ncol(ci), 4)
  expect_equal(colnames(ci), c("estimate", "se", "lower", "upper"))
  
  ## Lower should be less than estimate, estimate less than upper
  expect_true(all(ci[, "lower"] < ci[, "estimate"]))
  expect_true(all(ci[, "estimate"] < ci[, "upper"]))
  
  ## CI width should be approximately 2 * 1.96 * SE for 95% CI
  width <- ci[, "upper"] - ci[, "lower"]
  expected_width <- 2 * qnorm(0.975) * ci[, "se"]
  expect_equal(width, expected_width, tolerance = 1e-10)
})


test_that("99% CI is wider than 95% CI", {
  data <- generate_hessian_test_data(n = 200, seed = 1010)
  
  ## Find MLE first
  start_params <- c(alpha = 0.06, beta = 0.88)
  
  opt_result <- optim(
    par = start_params,
    fn = dcc11_nll,
    method = "L-BFGS-B",
    lower = c(1e-6, 1e-6),
    upper = c(0.5, 0.99),
    std_resid = data$std_resid,
    weights = data$weights,
    Qbar = data$Qbar,
    distribution = "mvn",
    use_reparam = FALSE
  )
  
  params_mle <- opt_result$par
  names(params_mle) <- c("alpha", "beta")
  
  se_result <- dcc11_standard_errors(
    params = params_mle,
    std_resid = data$std_resid,
    weights = data$weights,
    Qbar = data$Qbar,
    distribution = "mvn",
    use_reparam = FALSE
  )
  
  ci_95 <- dcc11_confint(se_result, level = 0.95)
  ci_99 <- dcc11_confint(se_result, level = 0.99)
  
  width_95 <- ci_95[, "upper"] - ci_95[, "lower"]
  width_99 <- ci_99[, "upper"] - ci_99[, "lower"]
  
  expect_true(all(width_99 > width_95))
})


#### ______________________________________________________________________ ####
#### PART 6: ESTIMATION SUMMARY                                             ####

test_that("dcc11_estimation_summary returns complete results", {
  data <- generate_hessian_test_data(n = 300, seed = 1111)
  
  ## Find MLE first
  start_params <- c(alpha = 0.06, beta = 0.88)
  
  opt_result <- optim(
    par = start_params,
    fn = dcc11_nll,
    method = "L-BFGS-B",
    lower = c(1e-6, 1e-6),
    upper = c(0.5, 0.99),
    std_resid = data$std_resid,
    weights = data$weights,
    Qbar = data$Qbar,
    distribution = "mvn",
    use_reparam = FALSE
  )
  
  params_mle <- opt_result$par
  names(params_mle) <- c("alpha", "beta")
  
  summary <- dcc11_estimation_summary(
    params = params_mle,
    std_resid = data$std_resid,
    weights = data$weights,
    Qbar = data$Qbar,
    distribution = "mvn",
    use_reparam = FALSE
  )
  
  expect_true(inherits(summary, "dcc11_summary"))
  
  ## Check key components
  expect_true(!is.null(summary$n_obs))
  expect_true(!is.null(summary$log_likelihood))
  expect_true(!is.null(summary$alpha))
  expect_true(!is.null(summary$beta))
  expect_true(!is.null(summary$persistence))
  expect_true(!is.null(summary$se))
  expect_true(!is.null(summary$ci))
  expect_true(!is.null(summary$aic))
  expect_true(!is.null(summary$bic))
  
  ## Persistence should be alpha + beta
  expect_equal(summary$persistence, summary$alpha + summary$beta)
})


#### ______________________________________________________________________ ####
#### PART 7: EDGE CASES                                                     ####

test_that("handles near-boundary parameters gracefully", {
  data <- generate_hessian_test_data(n = 200, seed = 1212)
  
  ## Find MLE, which may be near boundary
  ## Start with high persistence
  start_params <- c(alpha = 0.02, beta = 0.95)
  
  opt_result <- optim(
    par = start_params,
    fn = dcc11_nll,
    method = "L-BFGS-B",
    lower = c(1e-6, 1e-6),
    upper = c(0.5, 0.999),  ## Allow high persistence
    std_resid = data$std_resid,
    weights = data$weights,
    Qbar = data$Qbar,
    distribution = "mvn",
    use_reparam = FALSE
  )
  
  params_mle <- opt_result$par
  names(params_mle) <- c("alpha", "beta")
  
  se_result <- dcc11_standard_errors(
    params = params_mle,
    std_resid = data$std_resid,
    weights = data$weights,
    Qbar = data$Qbar,
    distribution = "mvn",
    use_reparam = FALSE
  )
  
  ## Should still produce results (possibly with warnings)
  expect_true(is.list(se_result))
  expect_true(length(se_result$se) == 2)
})


test_that("handles small samples", {
  data_small <- generate_hessian_test_data(n = 50, seed = 1313)
  data_large <- generate_hessian_test_data(n = 500, seed = 1313)
  
  ## Find MLE for small sample
  opt_small <- optim(
    par = c(0.05, 0.90),
    fn = dcc11_nll,
    method = "L-BFGS-B",
    lower = c(1e-6, 1e-6),
    upper = c(0.5, 0.99),
    std_resid = data_small$std_resid,
    weights = data_small$weights,
    Qbar = data_small$Qbar,
    distribution = "mvn",
    use_reparam = FALSE
  )
  
  params_small <- opt_small$par
  names(params_small) <- c("alpha", "beta")
  
  se_small <- dcc11_standard_errors(
    params = params_small,
    std_resid = data_small$std_resid,
    weights = data_small$weights,
    Qbar = data_small$Qbar,
    distribution = "mvn",
    use_reparam = FALSE
  )
  
  ## Find MLE for large sample
  opt_large <- optim(
    par = c(0.05, 0.90),
    fn = dcc11_nll,
    method = "L-BFGS-B",
    lower = c(1e-6, 1e-6),
    upper = c(0.5, 0.99),
    std_resid = data_large$std_resid,
    weights = data_large$weights,
    Qbar = data_large$Qbar,
    distribution = "mvn",
    use_reparam = FALSE
  )
  
  params_large <- opt_large$par
  names(params_large) <- c("alpha", "beta")
  
  se_large <- dcc11_standard_errors(
    params = params_large,
    std_resid = data_large$std_resid,
    weights = data_large$weights,
    Qbar = data_large$Qbar,
    distribution = "mvn",
    use_reparam = FALSE
  )
  
  ## Should produce results
  expect_true(is.list(se_small))
  
  ## SEs from small sample should be larger
  expect_gt(se_small$se["alpha"], se_large$se["alpha"])
  expect_gt(se_small$se["beta"], se_large$se["beta"])
})



#### ______________________________________________________________________ ####
#### PART 8: ROBUST SE WRAPPER TESTS                                        ####

test_that("compute_dcc_standard_errors_robust handles dynamic correlation", {
  data <- generate_hessian_test_data(n = 300, k = 2,
                                     alpha_true = 0.10, beta_true = 0.80,
                                     seed = 2001)
  
  ## Simulate a result from estimate_dcc_parameters_weighted
  dcc_result <- list(
    dcc_pars = list(alpha_1 = 0.08, beta_1 = 0.82),
    dist_pars = list(),
    correlation_type = "dynamic"
  )
  
  result <- compute_dcc_standard_errors_robust(
    dcc_result = dcc_result,
    std_resid = data$std_resid,
    weights = data$weights,
    Qbar = data$Qbar,
    distribution = "mvn"
  )
  
  expect_equal(result$correlation_type, "dynamic")
  expect_true(is.numeric(result$se))
  expect_equal(length(result$se), 2)
  expect_true("alpha" %in% names(result$se))
  expect_true("beta" %in% names(result$se))
  
  cat("\nRobust SE result (dynamic):\n")
  cat("  Valid:", result$valid, "\n")
  cat("  Reason:", result$reason, "\n")
  cat("  SEs:", round(result$se, 4), "\n")
})


test_that("compute_dcc_standard_errors_robust handles constant correlation", {
  data <- generate_hessian_test_data(n = 200, seed = 2002)
  
  ## Constant correlation result
  dcc_result <- list(
    dcc_pars = list(),
    dist_pars = list(),
    correlation_type = "constant"
  )
  
  result <- compute_dcc_standard_errors_robust(
    dcc_result = dcc_result,
    std_resid = data$std_resid,
    weights = data$weights,
    Qbar = data$Qbar,
    distribution = "mvn"
  )
  
  expect_equal(result$correlation_type, "constant")
  expect_equal(result$reason, "constant_correlation")
  expect_true(result$valid)  ## Valid because there's nothing to estimate
  expect_null(result$se)
})


test_that("compute_dcc_standard_errors_robust handles boundary estimates", {
  data <- generate_hessian_test_data(n = 200, seed = 2003)
  
  ## Near-boundary result (alpha very small)
  dcc_result <- list(
    dcc_pars = list(alpha_1 = 1e-6, beta_1 = 0.95),
    dist_pars = list(),
    correlation_type = "dynamic"
  )
  
  result <- compute_dcc_standard_errors_robust(
    dcc_result = dcc_result,
    std_resid = data$std_resid,
    weights = data$weights,
    Qbar = data$Qbar,
    distribution = "mvn",
    boundary_threshold = 1e-4
  )
  
  expect_false(result$valid)
  expect_true(grepl("boundary", result$reason))
  expect_true(all(is.na(result$se)))
  
  cat("\nRobust SE result (boundary):\n")
  cat("  Valid:", result$valid, "\n")
  cat("  Reason:", result$reason, "\n")
})


test_that("compute_dcc_standard_errors_robust handles MVT distribution", {
  data <- generate_hessian_test_data(n = 300, k = 2,
                                     alpha_true = 0.10, beta_true = 0.80,
                                     seed = 2004)
  
  dcc_result <- list(
    dcc_pars = list(alpha_1 = 0.08, beta_1 = 0.82),
    dist_pars = list(shape = 8.0),
    correlation_type = "dynamic"
  )
  
  result <- compute_dcc_standard_errors_robust(
    dcc_result = dcc_result,
    std_resid = data$std_resid,
    weights = data$weights,
    Qbar = data$Qbar,
    distribution = "mvt"
  )
  
  expect_equal(length(result$se), 3)
  expect_true("shape" %in% names(result$se))
  
  cat("\nRobust SE result (MVT):\n")
  cat("  Valid:", result$valid, "\n")
  cat("  SEs:", round(result$se, 4), "\n")
})


test_that("compute_dcc_standard_errors_robust handles nested coefficient structure", {
  data <- generate_hessian_test_data(n = 300, seed = 2005)
  
  ## Structure as returned by estimate_garch_weighted_r
  dcc_result <- list(
    coefficients = list(
      dcc_pars = list(alpha_1 = 0.08, beta_1 = 0.82),
      dist_pars = list(),
      correlation_type = "dynamic"
    )
  )
  
  result <- compute_dcc_standard_errors_robust(
    dcc_result = dcc_result,
    std_resid = data$std_resid,
    weights = data$weights,
    Qbar = data$Qbar,
    distribution = "mvn"
  )
  
  expect_equal(result$correlation_type, "dynamic")
  expect_true(is.numeric(result$se))
})
