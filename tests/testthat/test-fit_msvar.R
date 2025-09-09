#' Simulate data from a 2-state MS-VAR(1) process
#'
#' @param T_obs Number of observations to generate.
#' @param beta1, beta2 True VAR coefficient matrices for each state.
#' @param sigma1, sigma2 True error covariance matrices for each state.
#' @param P True 2x2 transition probability matrix.
#' @return A list containing the simulated data `y` and the true state sequence `s`.
simulate_msvar_data <- function(T_obs, beta1, beta2, sigma1, sigma2, P) {
  k <- ncol(beta1)
  p <- 1 ## VAR(1)
  
  ## Generate the Markov chain for the states
  states <- numeric(T_obs)
  states[1] <- sample(1:2, 1, prob = c(0.5, 0.5)) ## Start in a random state
  for (t in 2:T_obs) {
    states[t] <- sample(1:2, 1, prob = P[states[t - 1], ])
  }
  
  ## Generate the time series data
  y <- matrix(0, nrow = T_obs, ncol = k)
  y[1, ] <- MASS::mvrnorm(1, mu = rep(0, k), Sigma = diag(k)) ## Initial value
  
  for (t in (p + 1):T_obs) {
    X_t <- c(1, y[t - 1, ]) ## Intercept and lag
    
    if (states[t] == 1) {
      mean_t <- X_t %*% beta1
      y[t, ] <- MASS::mvrnorm(1, mu = mean_t, Sigma = sigma1)
    } else {
      mean_t <- X_t %*% beta2
      y[t, ] <- MASS::mvrnorm(1, mu = mean_t, Sigma = sigma2)
    }
  }
  
  ## Add column names for clarity
  colnames(y) <- paste0("y", 1:k)
  
  return(list(y = y, states = states))
}


## --- Test Data Setup ---
set.seed(42)
T_obs <- 400
y1 <- arima.sim(model = list(ar = 0.8), n = T_obs)
y2 <- 0.6 * y1 + arima.sim(model = list(ar = -0.4), n = T_obs)
test_data <- cbind(y1, y2)
univar_data <- matrix(y1, ncol = 1)


## == == == == == == == == == == == == == == == == == == == == == == ==
## Section 1: Input Validation and Error Handling
## == == == == == == == == == == == == == == == == == == == == == == ==

context("Input Validation and Error Handling")

test_that("Function rejects non-numeric or invalid inputs without crashing", {
  ## Test with a character string
  expect_error(
    fit_msvar("not a matrix"), 
    "Input 'y' must be a numeric matrix or data frame.",
    fixed = TRUE
  )
  
  ## Test with a data frame containing a character column
  expect_error(
    fit_msvar(data.frame(x = letters[1:5], y = 1:5)),
    "Input 'y' must be numeric. Character or other non-numeric data is not allowed.",
    fixed = TRUE
  )
  
  ## Test with NA values (this is the one that was failing)
  expect_error(
    fit_msvar(matrix(c(1, 2, 3, NA), 2, 2)),
    "Input matrix 'y' contains non-finite values (NA, NaN, Inf).",
    fixed = TRUE
  )
})

test_that("Function handles insufficient data", {
  ## This error message also has parentheses, so add fixed = TRUE
  expect_error(
    fit_msvar(matrix(1, 1, 2)), 
    "Input matrix 'y' must have at least 2 rows for a VAR(1) model.",
    fixed = TRUE
  )
})


## == == == == == == == == == == == == == == == == == == == == == == ==
## Section 2: Output Structure and Formatting
## == == == == == == == == == == == == == == == == == == == == == == ==

context("Output Structure and Formatting")

test_that("Output object is a correctly structured list", {
  fit <- fit_msvar(test_data, max_iter = 5) # Few iterations for speed
  
  expect_true(is.list(fit))
  expected_names <- c("beta1", "beta2", "sigma1", "sigma2", "P", "log_likelihood", "smoothed_probabilities")
  expect_named(fit, expected_names, ignore.order = TRUE)
  
  expect_true(is.numeric(fit$log_likelihood))
})

test_that("Output matrices have correct dimensions and names", {
  k <- ncol(test_data)
  p <- 1
  T_eff <- nrow(test_data) - p
  
  ## Test with named input data
  fit_named <- fit_msvar(test_data, max_iter = 5)
  
  ## Beta matrices
  expect_equal(dim(fit_named$beta1), c(1 + p * k, k))
  expect_equal(colnames(fit_named$beta1), colnames(test_data))
  expect_equal(rownames(fit_named$beta1), c("const", "y1_lag1", "y2_lag1"))
  expect_equal(dim(fit_named$beta2), dim(fit_named$beta1))
  expect_equal(dimnames(fit_named$beta2), dimnames(fit_named$beta1))
  
  ## Sigma matrices
  expect_equal(dim(fit_named$sigma1), c(k, k))
  expect_equal(dimnames(fit_named$sigma1), list(colnames(test_data), colnames(test_data)))
  expect_equal(dim(fit_named$sigma2), dim(fit_named$sigma1))
  
  ## Probabilities
  expect_equal(dim(fit_named$smoothed_probabilities), c(T_eff, 2))
  expect_equal(colnames(fit_named$smoothed_probabilities), c("State1_Prob", "State2_Prob"))
  
  ## Test with unnamed input data
  fit_unnamed <- fit_msvar(unname(test_data), max_iter = 5)
  expect_equal(colnames(fit_unnamed$beta1), c("y1", "y2"))
})


## == == == == == == == == == == == == == == == == == == == == == == ==
## Section 3: Mathematical Properties of the Output
## == == == == == == == == == == == == == == == == == == == == == == ==

context("Mathematical Properties of Output")

test_that("Output parameters satisfy mathematical constraints", {
  # Use a longer run to get more stable estimates
  fit <- fit_msvar(test_data, max_iter = 100)
  
  ## Transition matrix rows must sum to 1
  expect_equal(rowSums(fit$P), c(1, 1), tolerance = 1e-9)
  
  ## All transition probabilities must be in [0, 1]
  expect_true(all(fit$P >= 0 & fit$P <= 1))
  
  ## Covariance matrices must be symmetric
  expect_true(isSymmetric(fit$sigma1, tol = 1e-9))
  expect_true(isSymmetric(fit$sigma2, tol = 1e-9))
  
  ## Covariance matrices must be positive semi-definite (all eigenvalues >= 0)
  ## For non-degenerate data, they should be strictly positive definite
  expect_true(all(eigen(fit$sigma1, only.values = TRUE)$values > 0))
  expect_true(all(eigen(fit$sigma2, only.values = TRUE)$values > 0))
  
  ## Smoothed probabilities must sum to 1 for each observation
  expect_equal(
    rowSums(fit$smoothed_probabilities), 
    rep(1, nrow(fit$smoothed_probabilities)), 
    tolerance = 1e-9
  )
})

test_that("Log-likelihood increases or stays constant with each iteration", {
  ## This is a fundamental property of the EM algorithm
  fit_1_iter <- fit_msvar(test_data, max_iter = 1)
  fit_5_iter <- fit_msvar(test_data, max_iter = 5)
  fit_10_iter <- fit_msvar(test_data, max_iter = 10)
  
  expect_gte(fit_5_iter$log_likelihood, fit_1_iter$log_likelihood - 1e-9) # Allow for float tolerance
  expect_gte(fit_10_iter$log_likelihood, fit_5_iter$log_likelihood - 1e-9)
})


## == == == == == == == == == == == == == == == == == == == == == == ==
## Section 4: Core Algorithm Correctness (Simulation)
## == == == == == == == == == == == == == == == == == == == == == == ==

context("Core Algorithm Correctness (Simulation)")

test_that("Fitter recovers known parameters from simulated data", {
  skip_on_cran() # This test can be slow, so skip on CRAN checks
  
  ## 1. Define TRUE parameters for a 2-variable system
  true_beta1 <- matrix(c(0.1, 0.7, 0.1,  # For y1
                         0.2, 0.2, 0.3), # For y2
                       nrow = 3, byrow = FALSE)
  ## Add the correct dimension names
  dimnames(true_beta1) <- list(c("const", "y1_lag1", "y2_lag1"), c("y1", "y2"))
  
  true_beta2 <- matrix(c(-0.1, -0.5, 0.4,
                         -0.2, 0.2, -0.6),
                       nrow = 3, byrow = FALSE)
  ## Add the correct dimension names
  dimnames(true_beta2) <- list(c("const", "y1_lag1", "y2_lag1"), c("y1", "y2"))
  
  true_sigma1 <- matrix(c(1.0, 0.3,
                          0.3, 0.8), 2, 2)
  ## Add the correct dimension names
  dimnames(true_sigma1) <- list(c("y1", "y2"), c("y1", "y2"))
  
  true_sigma2 <- matrix(c(1.5, -0.5,
                          -0.5, 1.2), 2, 2)
  ## Add the correct dimension names
  dimnames(true_sigma2) <- list(c("y1", "y2"), c("y1", "y2"))
  
  true_P <- matrix(c(0.95, 0.05,
                     0.10, 0.90), 2, 2, byrow = TRUE)
  
  ## 2. Simulate data from this true model
  sim_data <- simulate_msvar_data(
    T_obs = 1000, # Use a large T for better estimates
    beta1 = true_beta1, beta2 = true_beta2,
    sigma1 = true_sigma1, sigma2 = true_sigma2,
    P = true_P
  )
  
  ## 3. Fit the model to the simulated data
  fit <- fit_msvar(sim_data$y, max_iter = 500, tol = 1e-7)
  
  ## 4. Check for Label Switching
  ## The algorithm might label our true "State 1" as its "State 2".
  ## We check which estimated state is closer to our true State 1.
  dist1 <- sum((fit$beta1 - true_beta1)^2)
  dist2 <- sum((fit$beta2 - true_beta1)^2)
  
  if (dist1 <= dist2) {
    ## No label switching
    est_beta1 <- fit$beta1; est_beta2 <- fit$beta2
    est_sigma1 <- fit$sigma1; est_sigma2 <- fit$sigma2
    est_P <- fit$P
  } else {
    ## Label switching occurred
    est_beta1 <- fit$beta2; est_beta2 <- fit$beta1
    est_sigma1 <- fit$sigma2; est_sigma2 <- fit$sigma1
    ## Swap rows and columns of P
    est_P <- fit$P[c(2, 1), c(2, 1)]
  }
  
  ## 5. Compare estimated parameters to true parameters
  ## The tolerance must be reasonably large due to statistical estimation noise.
  ## For a large T and a well-specified model, we expect to get close.
  expect_equal(est_beta1, true_beta1, tolerance = 0.2)
  expect_equal(est_beta2, true_beta2, tolerance = 0.2)
  expect_equal(est_sigma1, true_sigma1, tolerance = 0.25)
  expect_equal(est_sigma2, true_sigma2, tolerance = 0.25)
  expect_equal(est_P, true_P, tolerance = 0.1)
})


## == == == == == == == == == == == == == == == == == == == == == == ==
## Section 5: Comparison with a Reference Implementation
## == == == == == == == == == == == == == == == == == == == == == == ==

context("Comparison with a Reference Implementation (MSwM)")

test_that("Results are comparable to the MSwM package", {
  skip_on_cran()
  if (!requireNamespace("MSwM", quietly = TRUE)) {
    skip("MSwM package not available for comparison")
  }
  
  ## Use a simpler, univariate case for a clear comparison
  set.seed(123)
  univar_data <- arima.sim(list(order = c(1,0,0), ar = 0.5), n = 300)
  
  ## 1. Prepare data for MSwM (needs explicit lags)
  mswm_df <- data.frame(y = univar_data[-1], y_lag1 = univar_data[-length(univar_data)])
  
  ## 2. Fit with MSwM
  ## sw tells msmFit which coefficients to switch (intercept, lag, and variance)
  mod.mswm <- MSwM::msmFit(
    lm(y ~ y_lag1, data = mswm_df), 
    k = 2, 
    sw = c(TRUE, TRUE, TRUE) # Switch Intercept, Lag, and Variance
  )
  
  ## 3. Fit with our function
  mod.mine <- fit_msvar(matrix(univar_data, ncol = 1), max_iter = 200, tol = 1e-7)
  
  ## 4. Compare key results
  
  ## Log-likelihoods should be very close
  ## This is the best anchor point as it is invariant to label switching
  ## NOTE: mod.mine$log_likelihood is negative log-likelihood, while
  ## mod.mswm@Fit@logLikel is log-likelihood (not negative). Hence the minus.
  expect_equal(mod.mine$log_likelihood, -mod.mswm@Fit@logLikel, tolerance = 1e-3)
  
  ## Smoothed probabilities should also be highly correlated
  ## Handle label switching by checking correlation
  probs_mswm <- mod.mswm@Fit@smoProb
  cor1 <- cor(mod.mine$smoothed_probabilities[, 1], probs_mswm[-1, 1])
  ## NOTE: The state sequence is one observation shorter than the original 
  ## series because the VAR model requires an initial lag. We prepend the first 
  ## state to align the sequence with the original data's length. So we align 
  ## the two series by removing the first observation from the MSwM output
  
  if (cor1 < 0.5) { ## If correlation is negative, labels are swapped
    probs_mine_aligned <- mod.mine$smoothed_probabilities[, c(2, 1)]
  } else {
    probs_mine_aligned <- mod.mine$smoothed_probabilities
  }
  
  ## Check that the aligned probabilities are very close
  ## Use Mean Absolute Error as a metric
  expect_lt(mean(abs(probs_mine_aligned - probs_mswm[-1, ])), 0.05)
})
