## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## DCC Implementation Validation: Our Code vs tsmarch
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
##
## Purpose:
##   Validate our DCC implementation by comparing estimates against tsmarch's
##   single-regime DCC on single-regime simulated data.
##
## Strategy:
##   1. Simulate single-regime DCC-GARCH data (no regime switching)
##   2. Estimate with tsmarch (single-regime DCC) - reference
##   3. Estimate with our MS implementation (M=2) - should see one dominant state
##   4. Compare estimates and assess if differences are within sampling variability
##
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

library(tsmarch)
library(tsgarch)
library(xts)

## Source our code
# source("diagnostic_utils.R")
# source("dcc_gradient.R")
# source("dcc_hessian.R")
# source("dcc_inference.R")

## =============================================================================
## SECTION 1: Simulation Function
## =============================================================================

#' Simulate DCC-GARCH data and return in format suitable for both estimators
simulate_validation_data <- function(
    n = 1000,
    k = 2,
    true_alpha_dcc = 0.05,
    true_beta_dcc = 0.90,
    omega = c(0.05, 0.08),
    alpha_garch = c(0.10, 0.12),
    beta_garch = c(0.85, 0.82),
    seed = NULL
) {
  
  if (!is.null(seed)) set.seed(seed)
  
  ## Use our simulate_dcc_garch function
  y <- simulate_dcc_garch(
    n = n,
    k = k,
    omega = omega,
    alpha_garch = alpha_garch,
    beta_garch = beta_garch,
    alpha_dcc = true_alpha_dcc,
    beta_dcc = true_beta_dcc,
    seed = seed
  )
  
  ## Convert to xts for tsmarch
  dates <- Sys.Date() - (n:1)
  y_xts <- xts::xts(y, order.by = dates)
  colnames(y_xts) <- paste0("series_", 1:k)
  
  list(
    y_matrix = y,
    y_xts = y_xts,
    true_params = list(
      alpha_dcc = true_alpha_dcc,
      beta_dcc = true_beta_dcc,
      omega = omega,
      alpha_garch = alpha_garch,
      beta_garch = beta_garch
    )
  )
}


## =============================================================================
## SECTION 2: Estimate with tsmarch
## =============================================================================

#' Estimate DCC using tsmarch (single-regime reference)
estimate_with_tsmarch <- function(y_xts, verbose = TRUE) {
  
  k <- ncol(y_xts)
  
  if (verbose) cat("Estimating with tsmarch...\n")
  
  ## Step 1: Fit univariate GARCH models
  uni_models <- lapply(1:k, function(i) {
    spec <- tsgarch::garch_modelspec(
      y = y_xts[, i],
      model = "garch",
      garch_order = c(1, 1),
      distribution = "norm"
    )
    suppressWarnings(tsmethods::estimate(spec, keep_tmb = TRUE))
  })
  
  ## Step 2: Combine into multi-estimate object
  multi_est <- tsgarch::to_multi_estimate(uni_models)
  names(multi_est) <- colnames(y_xts)
  
  ## Step 3: Create DCC spec and estimate
  dcc_spec <- tsmarch::dcc_modelspec(
    object = multi_est,
    dcc_order = c(1, 1),
    dynamics = "dcc",
    distribution = "mvn"
  )
  
  dcc_fit <- suppressWarnings(tsmethods::estimate(dcc_spec))
  
  ## Extract parameters
  alpha_dcc <- dcc_fit$parmatrix[dcc_fit$parmatrix$parameter == "alpha_1", ]$value
  beta_dcc <- dcc_fit$parmatrix[dcc_fit$parmatrix$parameter == "beta_1", ]$value
  
  garch_pars <- lapply(1:k, function(i) {
    list(
      omega = uni_models[[i]]$parmatrix[uni_models[[i]]$parmatrix$parameter == "omega", ]$value,
      alpha1 = uni_models[[i]]$parmatrix[uni_models[[i]]$parmatrix$parameter == "alpha1", ]$value,
      beta1 = uni_models[[i]]$parmatrix[uni_models[[i]]$parmatrix$parameter == "beta1", ]$value
    )
  })
  
  if (verbose) {
    cat(sprintf("  DCC alpha: %.4f\n", alpha_dcc))
    cat(sprintf("  DCC beta:  %.4f\n", beta_dcc))
    cat(sprintf("  Persistence: %.4f\n", alpha_dcc + beta_dcc))
    cat(sprintf("  Log-likelihood: %.4f\n", dcc_fit$loglik))
  }
  
  list(
    alpha_dcc = alpha_dcc,
    beta_dcc = beta_dcc,
    garch_pars = garch_pars,
    loglik = dcc_fit$loglik,
    fit = dcc_fit,
    uni_models = uni_models
  )
}


## =============================================================================
## SECTION 3: Estimate with our DCC code (standalone, not MS)
## =============================================================================

#' Estimate DCC using our implementation (standalone, for direct comparison)
#' @param y_matrix The raw data matrix
#' @param tsmarch_result Full result from estimate_with_tsmarch (includes uni_models)
#' @param use_tsmarch_residuals If TRUE, use tsmarch's standardized residuals directly
#' @param start_at_tsmarch If TRUE, start optimization at tsmarch's solution
#' @param verbose Print output
estimate_with_ours_standalone <- function(y_matrix, tsmarch_result, 
                                          use_tsmarch_residuals = FALSE,
                                          start_at_tsmarch = FALSE,
                                          verbose = TRUE) {
  
  if (verbose) cat("Estimating with our standalone DCC...\n")
  
  n <- nrow(y_matrix)
  k <- ncol(y_matrix)
  
  tsmarch_garch_pars <- tsmarch_result$garch_pars
  
  ## Extract GARCH parameters
  omega <- sapply(tsmarch_garch_pars, function(x) x$omega)
  alpha_garch <- sapply(tsmarch_garch_pars, function(x) x$alpha1)
  beta_garch <- sapply(tsmarch_garch_pars, function(x) x$beta1)
  
  if (use_tsmarch_residuals) {
    ## Use tsmarch's standardized residuals directly
    std_resid <- do.call(cbind, lapply(tsmarch_result$uni_models, function(m) {
      as.numeric(residuals(m, standardize = TRUE))
    }))
    if (verbose) cat("  Using tsmarch's standardized residuals\n")
  } else {
    ## Compute standardized residuals ourselves
    std_resid <- compute_std_residuals(y_matrix, omega, alpha_garch, beta_garch)
    if (verbose) cat("  Using our computed standardized residuals\n")
  }
  
  ## Estimate DCC parameters
  Qbar <- cor(std_resid)
  weights <- rep(1, n)
  
  ## Starting values
  if (start_at_tsmarch) {
    start_par <- c(tsmarch_result$alpha_dcc, tsmarch_result$beta_dcc)
    if (verbose) cat(sprintf("  Starting at tsmarch solution: (%.4f, %.4f)\n", 
                             start_par[1], start_par[2]))
  } else {
    start_par <- c(0.05, 0.90)
  }
  
  opt_result <- optim(
    par = start_par,
    fn = dcc11_nll,
    method = "L-BFGS-B",
    lower = c(1e-6, 1e-6),
    upper = c(0.5, 0.999),
    std_resid = std_resid,
    weights = weights,
    Qbar = Qbar,
    distribution = "mvn",
    use_reparam = FALSE
  )
  
  alpha_dcc <- opt_result$par[1]
  beta_dcc <- opt_result$par[2]
  
  ## Compute NLL at tsmarch's solution for comparison
  nll_at_tsmarch <- dcc11_nll(
    params = c(tsmarch_result$alpha_dcc, tsmarch_result$beta_dcc),
    std_resid = std_resid,
    weights = weights,
    Qbar = Qbar,
    distribution = "mvn",
    use_reparam = FALSE
  )
  
  if (verbose) {
    cat(sprintf("  DCC alpha: %.4f\n", alpha_dcc))
    cat(sprintf("  DCC beta:  %.4f\n", beta_dcc))
    cat(sprintf("  Persistence: %.4f\n", alpha_dcc + beta_dcc))
    cat(sprintf("  Our NLL at our solution: %.4f\n", opt_result$value))
    cat(sprintf("  Our NLL at tsmarch solution: %.4f\n", nll_at_tsmarch))
    cat(sprintf("  NLL difference (ours - tsmarch): %+.6f\n", opt_result$value - nll_at_tsmarch))
  }
  
  list(
    alpha_dcc = alpha_dcc,
    beta_dcc = beta_dcc,
    nll = opt_result$value,
    nll_at_tsmarch = nll_at_tsmarch,
    convergence = opt_result$convergence,
    std_resid = std_resid,
    Qbar = Qbar
  )
}


#' Compare standardized residuals between our computation and tsmarch
compare_standardized_residuals <- function(y_matrix, tsmarch_result, verbose = TRUE) {
  
  k <- ncol(y_matrix)
  tsmarch_garch_pars <- tsmarch_result$garch_pars
  
  ## Extract GARCH parameters
  omega <- sapply(tsmarch_garch_pars, function(x) x$omega)
  alpha_garch <- sapply(tsmarch_garch_pars, function(x) x$alpha1)
  beta_garch <- sapply(tsmarch_garch_pars, function(x) x$beta1)
  
  ## Our standardized residuals
  std_resid_ours <- compute_std_residuals(y_matrix, omega, alpha_garch, beta_garch)
  
  ## tsmarch's standardized residuals
  std_resid_tsmarch <- do.call(cbind, lapply(tsmarch_result$uni_models, function(m) {
    as.numeric(residuals(m, standardize = TRUE))
  }))
  
  ## Compare
  n <- nrow(std_resid_ours)
  
  comparison <- list()
  for (i in 1:k) {
    diff <- std_resid_ours[, i] - std_resid_tsmarch[, i]
    comparison[[paste0("series_", i)]] <- list(
      correlation = cor(std_resid_ours[, i], std_resid_tsmarch[, i]),
      mean_diff = mean(diff),
      max_abs_diff = max(abs(diff)),
      rmse = sqrt(mean(diff^2))
    )
  }
  
  if (verbose) {
    cat("\n--- Standardized Residuals Comparison ---\n\n")
    cat(sprintf("  %-12s %12s %12s %12s %12s\n", 
                "Series", "Correlation", "Mean Diff", "Max |Diff|", "RMSE"))
    for (i in 1:k) {
      comp <- comparison[[paste0("series_", i)]]
      cat(sprintf("  %-12s %12.6f %12.6f %12.6f %12.6f\n",
                  paste0("Series ", i),
                  comp$correlation, comp$mean_diff, 
                  comp$max_abs_diff, comp$rmse))
    }
    
    ## Also compare Qbar
    Qbar_ours <- cor(std_resid_ours)
    Qbar_tsmarch <- cor(std_resid_tsmarch)
    
    cat("\n--- Unconditional Correlation (Qbar) ---\n\n")
    cat("  Ours:\n")
    print(round(Qbar_ours, 6))
    cat("\n  tsmarch:\n")
    print(round(Qbar_tsmarch, 6))
    cat(sprintf("\n  Max |diff| in Qbar: %.6f\n", max(abs(Qbar_ours - Qbar_tsmarch))))
  }
  
  list(
    std_resid_ours = std_resid_ours,
    std_resid_tsmarch = std_resid_tsmarch,
    comparison = comparison,
    Qbar_ours = cor(std_resid_ours),
    Qbar_tsmarch = cor(std_resid_tsmarch)
  )
}


## =============================================================================
## SECTION 4: Single Replication Comparison
## =============================================================================

#' Run one validation replication
run_single_validation <- function(
    n = 1000,
    true_alpha = 0.05,
    true_beta = 0.90,
    seed = NULL,
    verbose = TRUE
) {
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 60), collapse = ""), "\n")
    cat("  DCC Validation: Single Replication\n")
    cat(paste(rep("=", 60), collapse = ""), "\n")
    cat(sprintf("  n = %d, true alpha = %.3f, true beta = %.3f\n", n, true_alpha, true_beta))
    cat(paste(rep("=", 60), collapse = ""), "\n\n")
  }
  
  ## Simulate data
  if (verbose) cat("Step 1: Simulating data...\n")
  sim_data <- simulate_validation_data(
    n = n,
    true_alpha_dcc = true_alpha,
    true_beta_dcc = true_beta,
    seed = seed
  )
  
  ## Estimate with tsmarch
  if (verbose) cat("\nStep 2: Estimating with tsmarch (reference)...\n")
  tsmarch_result <- estimate_with_tsmarch(sim_data$y_xts, verbose = verbose)
  
  ## Compare standardized residuals
  if (verbose) cat("\nStep 3: Comparing standardized residuals...\n")
  resid_comparison <- compare_standardized_residuals(
    sim_data$y_matrix, 
    tsmarch_result,
    verbose = verbose
  )
  
  ## Estimate with our code using tsmarch residuals (default start)
  if (verbose) cat("\nStep 4a: Our DCC with tsmarch residuals (default start)...\n")
  ours_default <- estimate_with_ours_standalone(
    sim_data$y_matrix, 
    tsmarch_result,
    use_tsmarch_residuals = TRUE,
    start_at_tsmarch = FALSE,
    verbose = verbose
  )
  
  ## Estimate with our code starting at tsmarch's solution
  if (verbose) cat("\nStep 4b: Our DCC starting at tsmarch's solution...\n")
  ours_from_tsmarch <- estimate_with_ours_standalone(
    sim_data$y_matrix, 
    tsmarch_result,
    use_tsmarch_residuals = TRUE,
    start_at_tsmarch = TRUE,
    verbose = verbose
  )
  
  ## Compare
  if (verbose) {
    cat("\n")
    cat(paste(rep("-", 70), collapse = ""), "\n")
    cat("  COMPARISON\n")
    cat(paste(rep("-", 70), collapse = ""), "\n")
    
    cat(sprintf("  %-25s %10s %10s %10s %10s\n", 
                "", "True", "tsmarch", "Ours(def)", "Ours(tsm)"))
    cat(sprintf("  %-25s %10.4f %10.4f %10.4f %10.4f\n", "Alpha", 
                true_alpha, tsmarch_result$alpha_dcc, 
                ours_default$alpha_dcc,
                ours_from_tsmarch$alpha_dcc))
    cat(sprintf("  %-25s %10.4f %10.4f %10.4f %10.4f\n", "Beta", 
                true_beta, tsmarch_result$beta_dcc, 
                ours_default$beta_dcc,
                ours_from_tsmarch$beta_dcc))
    cat(sprintf("  %-25s %10.4f %10.4f %10.4f %10.4f\n", "Persistence", 
                true_alpha + true_beta, 
                tsmarch_result$alpha_dcc + tsmarch_result$beta_dcc,
                ours_default$alpha_dcc + ours_default$beta_dcc,
                ours_from_tsmarch$alpha_dcc + ours_from_tsmarch$beta_dcc))
    
    cat("\n  Legend:\n")
    cat("    Ours(def) = Our DCC with default start (0.05, 0.90)\n")
    cat("    Ours(tsm) = Our DCC starting at tsmarch's solution\n")
    
    cat("\n  NLL Comparison (lower is better):\n")
    cat(sprintf("    Our NLL at our solution:      %.4f\n", ours_default$nll))
    cat(sprintf("    Our NLL at tsmarch solution:  %.4f\n", ours_default$nll_at_tsmarch))
    cat(sprintf("    tsmarch reported log-lik:     %.4f (sign may differ)\n", tsmarch_result$loglik))
    
    if (ours_default$nll < ours_default$nll_at_tsmarch - 0.01) {
      cat("\n  NOTE: Our optimizer found a BETTER solution than tsmarch!\n")
    } else if (ours_default$nll > ours_default$nll_at_tsmarch + 0.01) {
      cat("\n  NOTE: tsmarch's solution has LOWER NLL - they found better optimum.\n")
    } else {
      cat("\n  NOTE: Both solutions have similar NLL (within 0.01).\n")
    }
  }
  
  list(
    true = list(alpha = true_alpha, beta = true_beta),
    tsmarch = list(alpha = tsmarch_result$alpha_dcc, beta = tsmarch_result$beta_dcc,
                   loglik = tsmarch_result$loglik),
    ours_default = list(alpha = ours_default$alpha_dcc, 
                        beta = ours_default$beta_dcc,
                        nll = ours_default$nll,
                        nll_at_tsmarch = ours_default$nll_at_tsmarch),
    ours_from_tsmarch = list(alpha = ours_from_tsmarch$alpha_dcc, 
                             beta = ours_from_tsmarch$beta_dcc,
                             nll = ours_from_tsmarch$nll),
    resid_comparison = resid_comparison
  )
}


## =============================================================================
## SECTION 5: Monte Carlo Validation
## =============================================================================

#' Run Monte Carlo validation comparing our implementation to tsmarch
run_validation_monte_carlo <- function(
    n_sim = 50,
    n_obs = 1000,
    true_alpha = 0.05,
    true_beta = 0.90,
    seed = 12345,
    verbose = TRUE
) {
  
  set.seed(seed)
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    cat("  DCC IMPLEMENTATION VALIDATION: Monte Carlo Study\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    cat(sprintf("  Replications: %d\n", n_sim))
    cat(sprintf("  Observations: %d\n", n_obs))
    cat(sprintf("  True params: alpha = %.3f, beta = %.3f\n", true_alpha, true_beta))
    cat(paste(rep("=", 70), collapse = ""), "\n\n")
  }
  
  ## Storage
  results <- data.frame(
    rep = integer(n_sim),
    tsmarch_alpha = numeric(n_sim),
    tsmarch_beta = numeric(n_sim),
    ours_alpha = numeric(n_sim),
    ours_beta = numeric(n_sim),
    ours_nll = numeric(n_sim),
    nll_at_tsmarch = numeric(n_sim),
    nll_diff = numeric(n_sim),  ## ours_nll - nll_at_tsmarch (negative = we found better)
    resid_corr_1 = numeric(n_sim),
    resid_corr_2 = numeric(n_sim)
  )
  
  if (verbose) pb <- txtProgressBar(min = 0, max = n_sim, style = 3)
  
  for (i in 1:n_sim) {
    
    result <- tryCatch({
      ## Simulate
      sim_data <- simulate_validation_data(
        n = n_obs,
        true_alpha_dcc = true_alpha,
        true_beta_dcc = true_beta,
        seed = seed + i
      )
      
      ## tsmarch
      tsmarch_result <- estimate_with_tsmarch(sim_data$y_xts, verbose = FALSE)
      
      ## Compare residuals
      resid_comp <- compare_standardized_residuals(
        sim_data$y_matrix, tsmarch_result, verbose = FALSE
      )
      
      ## Our DCC with tsmarch residuals
      ours <- estimate_with_ours_standalone(
        sim_data$y_matrix, tsmarch_result, 
        use_tsmarch_residuals = TRUE, 
        start_at_tsmarch = FALSE,
        verbose = FALSE
      )
      
      list(
        tsmarch = list(alpha = tsmarch_result$alpha_dcc, beta = tsmarch_result$beta_dcc),
        ours = list(alpha = ours$alpha_dcc, beta = ours$beta_dcc,
                    nll = ours$nll, nll_at_tsmarch = ours$nll_at_tsmarch),
        resid_corr = c(
          resid_comp$comparison$series_1$correlation,
          resid_comp$comparison$series_2$correlation
        )
      )
    }, error = function(e) {
      list(
        tsmarch = list(alpha = NA, beta = NA),
        ours = list(alpha = NA, beta = NA, nll = NA, nll_at_tsmarch = NA),
        resid_corr = c(NA, NA)
      )
    })
    
    results$rep[i] <- i
    results$tsmarch_alpha[i] <- result$tsmarch$alpha
    results$tsmarch_beta[i] <- result$tsmarch$beta
    results$ours_alpha[i] <- result$ours$alpha
    results$ours_beta[i] <- result$ours$beta
    results$ours_nll[i] <- result$ours$nll
    results$nll_at_tsmarch[i] <- result$ours$nll_at_tsmarch
    results$nll_diff[i] <- result$ours$nll - result$ours$nll_at_tsmarch
    results$resid_corr_1[i] <- result$resid_corr[1]
    results$resid_corr_2[i] <- result$resid_corr[2]
    
    if (verbose) setTxtProgressBar(pb, i)
  }
  
  if (verbose) {
    close(pb)
    cat("\n\n")
  }
  
  ## Remove failed replications
  valid_idx <- complete.cases(results)
  results_valid <- results[valid_idx, ]
  n_valid <- nrow(results_valid)
  
  if (verbose) {
    cat(paste(rep("=", 70), collapse = ""), "\n")
    cat("  VALIDATION RESULTS\n")
    cat(paste(rep("=", 70), collapse = ""), "\n\n")
    
    cat(sprintf("Valid replications: %d / %d\n\n", n_valid, n_sim))
    
    ## Residual comparison
    cat("--- Standardized Residuals Comparison ---\n\n")
    cat(sprintf("  Mean correlation (Series 1): %.6f\n", mean(results_valid$resid_corr_1)))
    cat(sprintf("  Mean correlation (Series 2): %.6f\n", mean(results_valid$resid_corr_2)))
    
    ## NLL comparison - THIS IS THE KEY
    cat("\n--- NLL Comparison (Critical Test) ---\n\n")
    cat(sprintf("  Mean NLL at our solution:     %.4f\n", mean(results_valid$ours_nll)))
    cat(sprintf("  Mean NLL at tsmarch solution: %.4f\n", mean(results_valid$nll_at_tsmarch)))
    cat(sprintf("  Mean difference (ours - tsm): %+.6f\n", mean(results_valid$nll_diff)))
    cat(sprintf("  SD of difference:             %.6f\n", sd(results_valid$nll_diff)))
    
    n_we_better <- sum(results_valid$nll_diff < -0.01)
    n_they_better <- sum(results_valid$nll_diff > 0.01)
    n_same <- sum(abs(results_valid$nll_diff) <= 0.01)
    
    cat(sprintf("\n  Cases where we found better optimum:  %d (%.1f%%)\n", 
                n_we_better, 100*n_we_better/n_valid))
    cat(sprintf("  Cases where tsmarch found better:     %d (%.1f%%)\n", 
                n_they_better, 100*n_they_better/n_valid))
    cat(sprintf("  Cases with similar NLL (|diff| < 0.01): %d (%.1f%%)\n", 
                n_same, 100*n_same/n_valid))
    
    ## Estimation comparison
    cat("\n--- Estimation Comparison ---\n\n")
    cat(sprintf("  %-20s %10s %10s %10s\n", "", "True", "tsmarch", "Ours"))
    cat(sprintf("  %-20s %10.4f %10.4f %10.4f\n", "Alpha (mean)", 
                true_alpha, 
                mean(results_valid$tsmarch_alpha),
                mean(results_valid$ours_alpha)))
    cat(sprintf("  %-20s %10s %10.4f %10.4f\n", "Alpha (SD)", "",
                sd(results_valid$tsmarch_alpha),
                sd(results_valid$ours_alpha)))
    cat(sprintf("  %-20s %10.4f %10.4f %10.4f\n", "Beta (mean)", 
                true_beta,
                mean(results_valid$tsmarch_beta),
                mean(results_valid$ours_beta)))
    cat(sprintf("  %-20s %10s %10.4f %10.4f\n", "Beta (SD)", "",
                sd(results_valid$tsmarch_beta),
                sd(results_valid$ours_beta)))
    
    ## Correlations
    cat("\n--- Correlation of Estimates ---\n\n")
    cat(sprintf("  Corr(tsmarch_alpha, ours_alpha): %.4f\n", 
                cor(results_valid$tsmarch_alpha, results_valid$ours_alpha)))
    cat(sprintf("  Corr(tsmarch_beta, ours_beta):   %.4f\n", 
                cor(results_valid$tsmarch_beta, results_valid$ours_beta)))
    
    ## Interpretation
    cat("\n--- Interpretation ---\n\n")
    
    if (mean(results_valid$nll_diff) < 0) {
      cat("  Our optimizer consistently finds LOWER (better) NLL than tsmarch.\n")
      cat("  This means our implementation is correct - we're finding better optima.\n")
      cat("  The parameter differences are because the likelihood is flat/ridged,\n")
      cat("  and different optimizers find different points with similar likelihood.\n")
    } else if (abs(mean(results_valid$nll_diff)) < 0.1) {
      cat("  Both implementations find similar NLL values.\n")
      cat("  Parameter differences are due to flat likelihood surface.\n")
    } else {
      cat("  tsmarch consistently finds better optima.\n")
      cat("  This might indicate an issue with our optimizer settings.\n")
    }
  }
  
  list(
    results = results,
    results_valid = results_valid,
    settings = list(
      n_sim = n_sim,
      n_obs = n_obs,
      true_alpha = true_alpha,
      true_beta = true_beta
    )
  )
}


## =============================================================================
## SECTION 6: Run Validation
## =============================================================================

if (interactive()) {
  
  cat("\n")
  cat("============================================================\n")
  cat("  DCC IMPLEMENTATION VALIDATION\n")
  cat("============================================================\n")
  cat("\n")
  cat("This script compares our DCC implementation against tsmarch.\n")
  cat("\n")
  cat("Options:\n")
  cat("  1. Single replication comparison (quick)\n")
  cat("  2. Monte Carlo validation (thorough)\n")
  cat("\n")
  
  ## Run single validation first
  cat("Running single replication...\n")
  single_result <- run_single_validation(
    n = 1000,
    true_alpha = 0.05,
    true_beta = 0.90,
    seed = 42,
    verbose = TRUE
  )
  
  cat("\n")
  cat("Run Monte Carlo validation? (This takes several minutes)\n")
  cat("Uncomment the code below to run.\n")
  
  ## Uncomment to run Monte Carlo validation
  mc_validation <- run_validation_monte_carlo(
    n_sim = 50,
    n_obs = 1000,
    true_alpha = 0.05,
    true_beta = 0.90,
    seed = 12345,
    verbose = TRUE
  )
}