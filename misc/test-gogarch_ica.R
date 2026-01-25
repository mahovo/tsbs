## Test script for GOGARCH ICA improvements
## Following Directions.txt: Step-by-step testing

library(tsmarch)

## ============================================================================
## TEST 1: Basic ICA on clean data
## ============================================================================

cat("\n=== TEST 1: Basic tsmarch::radical() on simulated data ===\n")

set.seed(42)
T_obs <- 252
k <- 5

# Generate clean test data similar to y_train
test_data <- matrix(rnorm(T_obs * k), nrow = T_obs, ncol = k)
colnames(test_data) <- c("SPY", "EFA", "BND", "GLD", "VNQ")

cat("Test data dimensions:", dim(test_data), "\n")
cat("Any non-finite values:", any(!is.finite(test_data)), "\n")

# Test basic radical
ica_result <- tryCatch({
  tsmarch::radical(test_data, components = k, demean = FALSE, trace = FALSE)
}, error = function(e) {
  cat("ERROR in radical():", e$message, "\n")
  NULL
})

if (!is.null(ica_result)) {
  cat("SUCCESS: radical() completed\n")
  cat("S dimensions:", dim(ica_result$S), "\n")
  cat("A dimensions:", dim(ica_result$A), "\n")
  cat("W dimensions:", dim(ica_result$W), "\n")
  cat("K dimensions:", dim(ica_result$K), "\n")
} else {
  cat("FAILED: radical() did not complete\n")
}

## ============================================================================
## TEST 2: Test improved_ica_decomposition on clean data
## ============================================================================

cat("\n=== TEST 2: improved_ica_decomposition() on clean data ===\n")

ica_improved <- tryCatch({
  improved_ica_decomposition(
    residuals = test_data,
    n_components = k,
    n_restarts = 1,  # Single restart for speed
    demean = FALSE,
    verbose = TRUE
  )
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  NULL
})

if (!is.null(ica_improved)) {
  cat("\nSUCCESS: improved_ica_decomposition() completed\n")
  cat("Method used:", ica_improved$method, "\n")
  cat("S dimensions:", dim(ica_improved$S), "\n")
  cat("Quality - Independence score:", ica_improved$quality$independence_score, "\n")
  cat("Quality - Total negentropy:", ica_improved$quality$total_negentropy, "\n")
} else {
  cat("FAILED: improved_ica_decomposition() did not complete\n")
}

## ============================================================================
## TEST 3: Test with data containing NaN
## ============================================================================

cat("\n=== TEST 3: improved_ica_decomposition() with NaN values ===\n")

test_data_nan <- test_data
test_data_nan[1, 1] <- NaN
test_data_nan[10, 3] <- Inf

cat("Introduced", sum(!is.finite(test_data_nan)), "non-finite values\n")

ica_nan <- tryCatch({
  improved_ica_decomposition(
    residuals = test_data_nan,
    n_components = k,
    n_restarts = 1,
    verbose = TRUE
  )
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  NULL
})

if (!is.null(ica_nan)) {
  cat("\nResult method:", ica_nan$method, "\n")
  cat("Is fallback:", ica_nan$quality$is_fallback, "\n")
  if (ica_nan$method == "identity_fallback_nan") {
    cat("SUCCESS: Gracefully handled NaN input with fallback\n")
  }
} else {
  cat("FAILED: Should have returned fallback, not error\n")
}

## ============================================================================
## TEST 4: Test estimate_garch_weighted_gogarch with clean data
## ============================================================================

cat("\n=== TEST 4: estimate_garch_weighted_gogarch() basic test ===\n")

# Create a simple spec similar to what's used in portfolio demo
spec <- list(
  garch_spec_args = list(
    model = "garch",
    order = c(1, 1),
    ica = "radical",
    components = k
  ),
  distribution = "norm",
  start_pars = list(
    garch_pars = lapply(1:k, function(i) list(omega = 0.02, alpha1 = 0.05, beta1 = 0.90)),
    dist_pars = NULL
  )
)

weights <- rep(1, T_obs)

gogarch_result <- tryCatch({
  estimate_garch_weighted_gogarch(
    residuals = test_data,
    weights = weights,
    spec = spec,
    diagnostics = NULL,
    verbose = TRUE
  )
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  print(traceback())
  NULL
})

if (!is.null(gogarch_result)) {
  cat("\nSUCCESS: estimate_garch_weighted_gogarch() completed\n")
  cat("Number of components:", length(gogarch_result$coefficients$garch_pars), "\n")
  cat("ICA method:", gogarch_result$coefficients$ica_info$method, "\n")
} else {
  cat("FAILED: estimate_garch_weighted_gogarch() did not complete\n")
}

## ============================================================================
## TEST 5: Test compute_gogarch_loglik_ms
## ============================================================================

cat("\n=== TEST 5: compute_gogarch_loglik_ms() ===\n")

if (!is.null(gogarch_result)) {
  loglik <- tryCatch({
    compute_gogarch_loglik_ms(
      residuals = test_data,
      garch_pars = gogarch_result$coefficients$garch_pars,
      ica_info = gogarch_result$coefficients$ica_info,
      distribution = "norm",
      return_vector = FALSE
    )
  }, error = function(e) {
    cat("ERROR:", e$message, "\n")
    NA
  })
  
  cat("Log-likelihood:", loglik, "\n")
  if (is.finite(loglik)) {
    cat("SUCCESS: Valid log-likelihood computed\n")
  } else {
    cat("FAILED: Log-likelihood is not finite\n")
  }
}

## ============================================================================
## SUMMARY
## ============================================================================

cat("\n=== TEST SUMMARY ===\n")
cat("Complete the above tests interactively and report results.\n")
cat("Key questions:\n")
cat("1. Does basic radical() work? (TEST 1)\n")
cat("2. Does improved_ica_decomposition() work on clean data? (TEST 2)\n")
cat("3. Does improved_ica_decomposition() handle NaN gracefully? (TEST 3)\n")
cat("4. Does estimate_garch_weighted_gogarch() work? (TEST 4)\n")
cat("5. Does compute_gogarch_loglik_ms() produce finite results? (TEST 5)\n")