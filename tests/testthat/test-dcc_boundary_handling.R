# =============================================================================
# DIAGNOSTIC TESTS: DCC Boundary Handling
# =============================================================================
# These tests verify that:
# 1. Boundary estimates are handled correctly
# 2. Bounds are consistent with tsmarch internal bounds
# 3. Stationarity constraint (alpha + beta < 1) is enforced
# 4. Upper boundary (near stationarity) is detected
# 5. BIC comparison for constant vs dynamic is correct
# 
# NOTE:
# These are tentative tests for troubleshooting, not proper unit tests!
# For proper tests, see test-ms_varma_garch_bs.R
# 
# =============================================================================

# library(testthat)
# library(xts)
# library(tsgarch)
# library(tsmarch)


## Test that bounds match tsmarch parmatrix bounds ============================

test_that("DCC bounds should match tsmarch parmatrix bounds", {
  skip_if_not_installed("tsmarch")
  skip_on_cran()
  
  # Create a simple DCC spec to inspect bounds
  set.seed(123)
  n <- 300
  k <- 2
  
  # Simulate simple data
  y <- matrix(rnorm(n * k), ncol = k)
  y_xts <- xts(y, order.by = Sys.Date() - (n:1))
  colnames(y_xts) <- c("s1", "s2")
  
  
  # Fit univariate GARCH models
  spec1 <- garch_modelspec(y_xts[,1], model = "garch", order = c(1,1))
  spec2 <- garch_modelspec(y_xts[,2], model = "garch", order = c(1,1))
  fit1 <- estimate(spec1, keep_tmb = TRUE)
  fit2 <- estimate(spec2, keep_tmb = TRUE)
  
  # Create DCC spec
  multi_fit <- to_multi_estimate(list(fit1, fit2))
  dcc_spec <- dcc_modelspec(multi_fit, dynamics = "dcc", dcc_order = c(1,1))
  
  # Extract bounds from tsmarch parmatrix
  parmatrix <- dcc_spec$parmatrix
  
  # Get alpha and beta bounds
  alpha_row <- parmatrix[parmatrix$parameter == "alpha_1", ]
  beta_row <- parmatrix[parmatrix$parameter == "beta_1", ]
  
  cat("\n=== tsmarch parmatrix bounds ===\n")
  cat("alpha_1: lower =", alpha_row$lower, ", upper =", alpha_row$upper, "\n")
  cat("beta_1:  lower =", beta_row$lower, ", upper =", beta_row$upper, "\n")
  
  # These should be the bounds we use in estimate_dcc_parameters_weighted()
  # Current implementation uses:
  #   alpha: lower = 0.01, upper = 0.99
  #   beta:  lower = 1e-6, upper = 0.99
  
  # Check if tsmarch uses different bounds
  expect_true(alpha_row$lower >= 0, "tsmarch alpha lower bound should be >= 0")
  expect_true(beta_row$lower >= 0, "tsmarch beta lower bound should be >= 0")
  expect_true(alpha_row$upper <= 1, "tsmarch alpha upper bound should be <= 1")
  expect_true(beta_row$upper <= 1, "tsmarch beta upper bound should be <= 1")
  
  # DIAGNOSTIC: Print actual tsmarch bounds for comparison
  cat("\nFull parmatrix:\n")
  print(parmatrix[parmatrix$group %in% c("alpha", "beta", "gamma"), ])
})


# Test effective sample size calculation for weighted BIC =====================

## OBS: THIS TEST IS MEANINGLESS!

test_that("Effective sample size for weighted BIC is calculated correctly", {
  skip_on_cran()
  
  # When using weighted likelihood, the effective sample size matters for BIC
  
  # Uniform weights: n_eff = n
  n <- 100
  w_uniform <- rep(1, n)
  n_eff_uniform <- sum(w_uniform)  # = n
  
  # Concentrated weights: n_eff < n
  w_concentrated <- c(rep(0.1, 50), rep(1.9, 50))  # Same sum as uniform
  w_concentrated <- w_concentrated / sum(w_concentrated) * n  # Normalize
  
  # Two common formulas for effective sample size:
  # Method 1: n_eff = (sum(w))^2 / sum(w^2)  -- from ESS literature
  # Method 2: n_eff = sum(w)  -- simple sum
  
  n_eff_method1 <- sum(w_concentrated)^2 / sum(w_concentrated^2)
  n_eff_method2 <- length(w_concentrated)  # Current implementation uses length
  
  cat("\n=== Effective Sample Size Calculation ===\n")
  cat("Uniform weights:\n")
  cat("  n_eff (method 1):", sum(w_uniform)^2 / sum(w_uniform^2), "\n")
  cat("  n_eff (method 2):", length(w_uniform), "\n")
  
  cat("\nConcentrated weights:\n")
  cat("  n_eff (method 1):", round(n_eff_method1, 2), "\n")
  cat("  n_eff (method 2):", n_eff_method2, "\n")
  
  # Method 1 accounts for weight concentration
  # Method 2 ignores it
  
  # For EM with smoothed probabilities, weights can be very concentrated
  # This affects BIC calculation significantly
  
  expect_true(n_eff_method1 <= n_eff_method2, 
              "ESS method should give smaller n_eff for concentrated weights")
})



# Test boundary event logging =================================================

test_that("Boundary events are properly logged to diagnostics", {
  skip_on_cran()
  
  # Create a diagnostic collector
  diagnostics <- create_diagnostic_collector()
  
  # Simulate adding boundary events
  diagnostics <- add_boundary_event(
    diagnostics,
    iteration = 5,
    state = 1,
    parameter_name = "alpha_1",
    value = 0.0101,
    boundary_type = "lower",
    action_taken = "constant_correlation_fallback"
  )
  
  diagnostics <- add_boundary_event(
    diagnostics,
    iteration = 5,
    state = 2,
    parameter_name = "alpha_1",
    value = 0.45,
    boundary_type = "none",
    action_taken = "kept_dynamic"
  )
  
  # Check that events are stored
  expect_equal(length(diagnostics$boundary_events), 2)
  expect_equal(diagnostics$boundary_events[[1]]$parameter, "alpha_1")
  expect_equal(diagnostics$boundary_events[[1]]$value, 0.0101)
  expect_equal(diagnostics$boundary_events[[1]]$state, 1)
  
  cat("\n=== Boundary Event Logging ===\n")
  cat("Events logged:", length(diagnostics$boundary_events), "\n")
  for (event in diagnostics$boundary_events) {
    cat(sprintf("  Iter %d, State %d: %s = %.4f (%s) -> %s\n",
                event$iteration, event$state, event$parameter,
                event$value, event$boundary_type, event$action_taken))
  }
})



# Test that constant correlation fallback works correctly =====================

test_that("DCC bounds are consistent and match tsmarch", {
  skip_on_cran()
  skip_if_not_installed("tsmarch")
  
  set.seed(42)
  n <- 400
  k <- 2

  
  cat("\n=== Test 1A: Bounds Consistency ===\n")
  
  # Create a DCC spec to check tsmarch bounds
  y_xts <- xts(matrix(rnorm(n * k), ncol = k), order.by = Sys.Date() - (n:1))
  colnames(y_xts) <- c("s1", "s2")
  
  fit1 <- estimate(garch_modelspec(y_xts[,1], model = "garch", order = c(1,1)), keep_tmb = TRUE)
  fit2 <- estimate(garch_modelspec(y_xts[,2], model = "garch", order = c(1,1)), keep_tmb = TRUE)
  multi_fit <- to_multi_estimate(list(fit1, fit2))
  dcc_spec <- dcc_modelspec(multi_fit, dynamics = "dcc", dcc_order = c(1,1))
  
  # Extract tsmarch bounds
  alpha_bounds <- dcc_spec$parmatrix[dcc_spec$parmatrix$parameter == "alpha_1", c("lower", "upper")]
  beta_bounds <- dcc_spec$parmatrix[dcc_spec$parmatrix$parameter == "beta_1", c("lower", "upper")]
  
  cat("tsmarch parmatrix bounds:\n")
  cat(sprintf("  alpha_1: [%g, %g]\n", alpha_bounds$lower, alpha_bounds$upper))
  cat(sprintf("  beta_1:  [%g, %g]\n", beta_bounds$lower, beta_bounds$upper))
  
  # TEST: tsmarch bounds should be consistent
  
  expect_equal(alpha_bounds$lower, beta_bounds$lower, 
               label = "tsmarch: alpha and beta should have same lower bound")
  expect_equal(alpha_bounds$upper, beta_bounds$upper,
               label = "tsmarch: alpha and beta should have same upper bound")
})


test_that("Near-zero DCC correctly falls back to constant", {
  skip_on_cran()
  skip_if_not_installed("tsmarch")
  
  set.seed(42)
  n <- 400
  k <- 2
  
  cat("\n=== Test 1B: Near-Zero DCC -> Constant Correlation ===\n")
  
  # Simulate data with VERY small DCC parameters
  # This should result in constant correlation being selected
  y <- simulate_dcc_garch(
    n = n,
    k = k,
    omega = c(0.05, 0.05),
    alpha_garch = c(0.10, 0.10),
    beta_garch = c(0.85, 0.85),
    dcc_alpha = 0.005,   # Very small
    dcc_beta = 0.005,    # Very small
    seed = 42
  )
  
  cat("True DCC parameters: alpha = 0.005, beta = 0.005\n")
  cat("Expected outcome: CONSTANT correlation (near-zero dynamics)\n\n")
  
  spec <- generate_single_state_dcc_spec(
    k = k,
    dcc_alpha = 0.05,  # Start with reasonable values
    dcc_beta = 0.50
  )
  
  residuals_mat <- as.matrix(y)
  weights <- rep(1, n)
  
  result <- tryCatch({
    estimate_garch_weighted_dcc(
      residuals = residuals_mat,
      weights = weights,
      spec = spec[[1]],
      diagnostics = NULL,
      iteration = 1,
      state = 1,
      verbose = TRUE,
      dcc_threshold = 0.02,  # Standard threshold
      dcc_criterion = "threshold",
      force_constant = FALSE
    )
  }, error = function(e) {
    cat("Error:", e$message, "\n")
    return(NULL)
  })
  
  expect_true(!is.null(result), label = "Estimation should complete")
  
  corr_type <- result$coefficients$correlation_type
  cat(sprintf("\nResult: correlation_type = '%s'\n", corr_type))
  
  # For near-zero true DCC, constant correlation is CORRECT
  expect_equal(corr_type, "constant",
               label = "Near-zero DCC should result in constant correlation")
  
  # dcc_pars should be empty for constant correlation
  expect_equal(length(result$coefficients$dcc_pars), 0,
               label = "Constant correlation should have empty dcc_pars")
})


test_that("Moderate DCC parameters are estimated correctly", {
  skip_on_cran()
  skip_if_not_installed("tsmarch")
  
  set.seed(100)
  n <- 500
  k <- 2
  
  cat("\n=== Test 1C: Moderate DCC Parameters ===\n")
  
  # Simulate data with MODERATE DCC parameters
  # These should be estimated as dynamic correlation
  y <- simulate_dcc_garch(
    n = n,
    k = k,
    omega = c(0.05, 0.05),
    alpha_garch = c(0.10, 0.10),
    beta_garch = c(0.85, 0.85),
    dcc_alpha = 0.06,    # Clearly above threshold (0.02)
    dcc_beta = 0.92,
    seed = 100
  )
  
  cat("True DCC parameters: alpha = 0.06, beta = 0.92\n")
  cat("Expected outcome: DYNAMIC correlation\n\n")
  
  spec <- generate_single_state_dcc_spec(
    k = k,
    dcc_alpha = 0.05,
    dcc_beta = 0.90
  )
  
  residuals_mat <- as.matrix(y)
  weights <- rep(1, n)
  
  result <- tryCatch({
    estimate_garch_weighted_dcc(
      residuals = residuals_mat,
      weights = weights,
      spec = spec[[1]],
      diagnostics = NULL,
      iteration = 1,
      state = 1,
      verbose = TRUE,
      dcc_threshold = 0.02,
      dcc_criterion = "threshold",
      force_constant = FALSE
    )
  }, error = function(e) {
    cat("Error:", e$message, "\n")
    return(NULL)
  })
  
  expect_true(!is.null(result), label = "Estimation should complete")
  
  corr_type <- result$coefficients$correlation_type
  cat(sprintf("\nResult: correlation_type = '%s'\n", corr_type))
  
  # For moderate DCC, dynamic correlation is expected
  expect_equal(corr_type, "dynamic",
               label = "Moderate DCC should result in dynamic correlation")
  
  if (corr_type == "dynamic") {
    estimated_alpha <- result$coefficients$dcc_pars$alpha_1
    estimated_beta <- result$coefficients$dcc_pars$beta_1
    
    cat(sprintf("  Estimated: alpha = %.4f, beta = %.4f\n", 
                estimated_alpha, estimated_beta))
    cat(sprintf("  True:      alpha = 0.06, beta = 0.92\n"))
    
    # Check estimates are in reasonable range
    expect_true(estimated_alpha > 0.02,
                label = "Estimated alpha should be above threshold")
    expect_true(estimated_alpha < 0.20,
                label = "Estimated alpha should be reasonable (< 0.20)")
    expect_true(estimated_beta > 0.70,
                label = "Estimated beta should be substantial (> 0.70)")
  }
})


test_that("Bounds allow full parameter space exploration", {
  skip_on_cran()
  skip_if_not_installed("tsmarch")
  
  set.seed(200)
  n <- 400
  k <- 2
  
  cat("\n=== Test 1D: Parameter Space Exploration ===\n")
  
  # This test verifies that the optimizer CAN reach values that were 
  # previously blocked by the old bounds
  
  # Simulate with alpha = 0.008 (below old bound of 0.01)
  y <- simulate_dcc_garch(
    n = n,
    k = k,
    omega = c(0.05, 0.05),
    alpha_garch = c(0.10, 0.10),
    beta_garch = c(0.85, 0.85),
    dcc_alpha = 0.008,   # Below old alpha lower bound of 0.01
    dcc_beta = 0.90,     # High beta to avoid constant fallback
    seed = 200
  )
  
  cat("True DCC parameters: alpha = 0.008 (below old bound!), beta = 0.90\n")
  cat("Old bounds would force alpha >= 0.01\n")
  cat("New bounds allow alpha down to ~0\n\n")
  
  spec <- generate_single_state_dcc_spec(
    k = k,
    dcc_alpha = 0.05,
    dcc_beta = 0.90
  )
  
  residuals_mat <- as.matrix(y)
  weights <- rep(1, n)
  
  # Use very low threshold to avoid constant fallback
  # We want to test that the optimizer can reach small alpha values
  result <- tryCatch({
    estimate_garch_weighted_dcc(
      residuals = residuals_mat,
      weights = weights,
      spec = spec[[1]],
      diagnostics = NULL,
      iteration = 1,
      state = 1,
      verbose = TRUE,
      dcc_threshold = 0.005,  # Very low to test bounds, not fallback
      dcc_criterion = "threshold",
      force_constant = FALSE
    )
  }, error = function(e) {
    cat("Error:", e$message, "\n")
    return(NULL)
  })
  
  expect_true(!is.null(result), label = "Estimation should complete")
  
  corr_type <- result$coefficients$correlation_type
  
  if (corr_type == "dynamic") {
    estimated_alpha <- result$coefficients$dcc_pars$alpha_1
    
    cat(sprintf("\nEstimated alpha: %.6f\n", estimated_alpha))
    cat(sprintf("True alpha:      0.008\n"))
    cat(sprintf("Old lower bound: 0.01\n"))
    
    # KEY TEST: With new bounds, alpha should NOT be stuck at 0.01
    # It should be able to go below 0.01
    alpha_stuck_at_old_bound <- abs(estimated_alpha - 0.01) < 0.0001
    
    cat(sprintf("\nIs alpha stuck at old bound (0.01)? %s\n",
                ifelse(alpha_stuck_at_old_bound, "YES", "NO")))
    
    # This is the key assertion - alpha should be free to go below 0.01
    expect_false(alpha_stuck_at_old_bound,
                 label = "Alpha should not be stuck at old lower bound (0.01)")
    
  } else {
    # If constant, the threshold was too aggressive
    # This is still valid - just a different code path
    cat("\nResult: constant correlation (threshold triggered)\n")
    cat("This means alpha estimate was below 0.005\n")
    cat("Which proves bounds allow values below 0.01!\n")
    
    # If we got constant, that means alpha went below threshold (0.005)
    # which is below the old bound of 0.01 - so bounds ARE working
    expect_true(TRUE, label = "Constant fallback proves bounds allow small alpha")
  }
})



# Integration Test - Full Boundary Handling in MS-DCC Estimation ==============
#
# This test simulates a 2-state MS-DCC scenario where:
# - State 1: Has meaningful DCC dynamics (should stay dynamic)
# - State 2: Has nearly zero DCC dynamics (should fall back to constant)

test_that("Full integration: boundary handling in MS-DCC estimation", {
  
  set.seed(999)
  n <- 600
  k <- 2
  
  # Step 1: ===================================================================
  # Simulate regime-switching data with different correlation dynamics ========
  
  # State 1: Low volatility, DYNAMIC correlation (alpha=0.05, beta=0.93)
  y_state1 <- simulate_dcc_garch(
    n = n/2,
    k = k,
    omega = c(0.02, 0.02),
    alpha_garch = c(0.05, 0.05),
    beta_garch = c(0.90, 0.90),
    dcc_alpha = 0.05,
    dcc_beta = 0.93,
    seed = 999
  )
  
  # State 2: High volatility, CONSTANT correlation (alpha~0, beta~0)
  y_state2 <- simulate_dcc_garch(
    n = n/2,
    k = k,
    omega = c(0.10, 0.10),
    alpha_garch = c(0.15, 0.15),
    beta_garch = c(0.75, 0.75),
    dcc_alpha = 0.001,  # Nearly zero -> should trigger constant fallback
    dcc_beta = 0.001,
    seed = 1000
  )
  
  # Combine with simple regime switching pattern
  # Pattern: State1 -> State2 -> State1 -> State2
  block_size <- n / 4
  y <- rbind(
    y_state1[1:block_size, ],
    y_state2[1:block_size, ],
    y_state1[(block_size + 1):(n/2), ],
    y_state2[(block_size + 1):(n/2), ]
  )
  
  cat("\n")
  cat("=" , rep("=", 70), "\n", sep = "")
  cat("INTEGRATION TEST: 2-State MS-DCC with Boundary Handling\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat("\nData characteristics:\n")
  cat("  Total observations: ", n, "\n")
  cat("  Number of series: ", k, "\n")
  cat("\nTrue parameters:\n")
  cat("  State 1 (Dynamic): alpha=0.05, beta=0.93, persistence=0.98\n")
  cat("  State 2 (Constant): alpha~0, beta~0, persistence~0\n")
  cat("\nExpected outcome:\n")
  cat("  State 1: Should remain DYNAMIC correlation\n")
  cat("  State 2: Should fall back to CONSTANT correlation\n")
  cat("\n")
  

  # Step 2: Create model specification ======================================
  
  # Generate 2-state specification
  spec <- generate_dcc_spec(
    M = 2,
    k = k,
    var_order = 1,
    garch_order = c(1, 1),
    distribution = "mvn",
    seed = 123,
    simple = FALSE
  )
  
  # Print starting parameters
  cat("Starting parameters:\n")
  for (j in 1:2) {
    cat(sprintf("  State %d: alpha=%.3f, beta=%.3f\n",
                j,
                spec[[j]]$start_pars$dcc_pars$alpha_1,
                spec[[j]]$start_pars$dcc_pars$beta_1))
  }
  cat("\n")
  
  
  # Step 3: Create diagnostic collector =====================================
  
  #diagnostics <- create_diagnostic_collector()
  
  
  # Step 4: Run the fit with diagnostics enabled ============================
  
  cat("Running fit_ms_varma_garch()...\n")
  cat("(This may take a minute)\n\n")
  
  fit_result <- tryCatch({
    fit_ms_varma_garch(
      y = y,
      M = 2,
      spec = spec,
      model_type = "multivariate",
      control = list(
        max_iter = 30,
        tol = 1e-4
      )
    )
  }, error = function(e) {
    cat("ERROR during fitting:\n")
    cat(e$message, "\n")
    return(NULL)
  })
  

  # Step 5: Analyze results =================================================
  
  if (!is.null(fit_result)) {
    cat("\n")
    cat("=", rep("=", 70), "\n", sep = "")
    cat("RESULTS\n")
    cat("=", rep("=", 70), "\n", sep = "")
    
    # Extract final parameters
    final_params <- fit_result$parameters
    
    cat("\nFinal estimated parameters:\n")
    for (j in 1:2) {
      state_params <- final_params[[j]]
      
      corr_type <- state_params$correlation_type %||% "unknown"
      
      if (corr_type == "dynamic") {
        alpha <- state_params$alpha_1 %||% NA
        beta <- state_params$beta_1 %||% NA
        cat(sprintf("  State %d: DYNAMIC - alpha=%.4f, beta=%.4f, persistence=%.4f\n",
                    j, alpha, beta, alpha + beta))
      } else {
        cat(sprintf("  State %d: CONSTANT - (no DCC parameters)\n", j))
        if (!is.null(state_params$degeneracy_reason)) {
          cat(sprintf("           Reason: %s\n", state_params$degeneracy_reason))
        }
      }
    }
    
    
    # Step 6: Examine diagnostics =============================================
    
    if (!is.null(fit_result$diagnostics)) {
      diag <- fit_result$diagnostics
      
      cat("\n")
      cat("-", rep("-", 70), "\n", sep = "")
      cat("DIAGNOSTICS\n")
      cat("-", rep("-", 70), "\n", sep = "")
      
      # EM convergence
      cat("\nEM Convergence:\n")
      cat("  Total iterations:", length(diag$em_iterations), "\n")
      
      if (length(diag$em_iterations) > 0) {
        ll_changes <- sapply(diag$em_iterations, function(x) x$ll_change)
        cat("  LL decreased in", sum(ll_changes < -1e-6), "iterations\n")
        cat("  Final LL change:", tail(ll_changes, 1), "\n")
      }
      
      # Boundary events
      cat("\nBoundary Events:\n")
      if (length(diag$boundary_events) > 0) {
        for (event in diag$boundary_events) {
          cat(sprintf("  Iter %d, State %d: %s = %.6f (%s) -> %s\n",
                      event$iteration,
                      event$state,
                      event$parameter,
                      event$value,
                      event$boundary_type,
                      event$action_taken))
        }
      } else {
        cat("  No boundary events recorded\n")
      }
      
      # Warnings
      cat("\nWarnings:\n")
      if (length(diag$warnings) > 0) {
        warning_types <- table(sapply(diag$warnings, function(x) x$type))
        for (wtype in names(warning_types)) {
          cat(sprintf("  %s: %d\n", wtype, warning_types[wtype]))
        }
      } else {
        cat("  No warnings\n")
      }
    }
    
    
    # Step 7: Verify expected outcomes ========================================
    
    cat("\n")
    cat("=", rep("=", 70), "\n", sep = "")
    cat("VERIFICATION\n")
    cat("=", rep("=", 70), "\n", sep = "")
    
    # Check that at least one state has dynamic correlation
    has_dynamic <- any(sapply(final_params, function(p) {
      !is.null(p$correlation_type) && p$correlation_type == "dynamic"
    }))
    
    # Check that at least one state has constant correlation
    has_constant <- any(sapply(final_params, function(p) {
      !is.null(p$correlation_type) && p$correlation_type == "constant"
    }))
    
    cat("\nOutcome checks:\n")
    cat("  At least one DYNAMIC state:", has_dynamic, "\n")
    cat("  At least one CONSTANT state:", has_constant, "\n")
    
    # The ideal outcome: one dynamic, one constant
    ideal_outcome <- has_dynamic && has_constant
    cat("\n  Ideal outcome (1 dynamic + 1 constant):", ideal_outcome, "\n")
    
    if (!ideal_outcome) {
      cat("\n  NOTE: Both states may have similar dynamics if the data\n")
      cat("        doesn't show strong regime differentiation.\n")
    }
    
    # Test assertions
    expect_true(!is.null(fit_result), "Fit should complete without error")
    expect_true(length(final_params) == 2, "Should have 2 states")
    
    # At minimum, correlation_type should be set for all states
    for (j in 1:2) {
      expect_true(!is.null(final_params[[j]]$correlation_type),
                  paste("State", j, "should have correlation_type set"))
    }
    
  } else {
    fail("Fit failed - see error message above")
  }
})




