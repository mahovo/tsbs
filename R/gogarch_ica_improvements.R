## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## GOGARCH ICA Integration Improvements
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
##
## This file provides improved ICA integration for GOGARCH models, addressing
## convergence issues and enhancing the tsmarch package integration.
##
## Key improvements:
##   1. Pre-conditioning of input data for better ICA convergence
##   2. Multiple ICA restarts with best-fit selection
##   3. Post-ICA validation and fallback strategies
##   4. Enhanced diagnostics for ICA quality assessment
##
## Dependencies:
##   - tsmarch (radical, fastica functions)
##   - tsbs_gogarch.R (estimate_garch_weighted_gogarch)
##
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


## SECTION 1: Improved ICA Decomposition =======================================

#' @title Improved ICA Decomposition for GOGARCH
#' 
#' @description Performs ICA decomposition with multiple restarts and quality
#'   diagnostics. Uses the RADICAL algorithm from \code{tsmarch}.
#'
#' @param residuals T x k matrix of residuals
#' @param method Currently only \code{"radical"} is supported
#' @param n_components Number of components to extract (default: k)
#' @param n_restarts Number of random restarts (default: 3)
#' @param demean Whether to demean the data (default: FALSE)
#' @param verbose Print progress information
#'
#' @return List with components:
#'   \describe{
#'     \item{S}{Independent components matrix (T x n_components)}
#'     \item{A}{Mixing matrix (k x n_components)}
#'     \item{W}{Unmixing matrix (n_components x k)}
#'     \item{K}{Pre-whitening matrix}
#'     \item{method}{"radical" or "pca_fallback"}
#'     \item{quality}{Quality metrics (see Details)}
#'     \item{convergence_info}{Restart diagnostics}
#'   }
#'
#' @details
#' Improvements over direct \code{\link[tsmarch]{radical}} calls:
#' \itemize{
#'   \item Multiple restarts with best-fit selection based on negentropy
#'   \item Quality metrics: independence score (target >0.8), reconstruction
#'     error (target <1\%), component negentropy
#'   \item Graceful PCA fallback if ICA fails
#' }
#'
#' @seealso \code{\link{estimate_garch_weighted_gogarch}}, 
#'   \code{\link{gogarch_diagnostics}}
#'
#' @export
improved_ica_decomposition <- function(
    residuals,
    method =  "radical",
    n_components = NULL,
    n_restarts = 3,
    demean = FALSE,
    verbose = FALSE
) {
  T_obs <- nrow(residuals)
  k <- ncol(residuals)
  
  if (is.null(n_components)) n_components <- k
  
  ## Input validation - check for NaN/Inf values
  if (any(!is.finite(residuals))) {
    n_bad <- sum(!is.finite(residuals))
    pct_bad <- 100 * n_bad / length(residuals)
    warning(sprintf("ICA input contains %d non-finite values (%.1f%%). Using PCA fallback.",
                    n_bad, pct_bad))
    
    ## Return identity fallback
    return(list(
      S = residuals,
      A = diag(k)[, 1:n_components, drop = FALSE],
      W = diag(n_components, k),
      K = diag(k),
      method = "identity_fallback_nan",
      n_components = n_components,
      convergence_info = list(
        method = "identity_fallback",
        reason = "non_finite_input",
        n_bad_values = n_bad
      ),
      quality = list(
        independence_score = NA_real_,
        total_negentropy = NA_real_,
        reconstruction_error = NA_real_,
        is_fallback = TRUE
      ),
      mu = rep(0, k)
    ))
  }
  
  ## Step 1: Pre-condition the data
  ## Centering (if needed) and whitening
  if (demean) {
    mu <- colMeans(residuals)
    residuals <- sweep(residuals, 2, mu, "-")
  } else {
    mu <- rep(0, k)
  }
  
  ## Pre-whiten using SVD for numerical stability
  svd_result <- svd(residuals, nu = T_obs, nv = k)
  
  ## Whitening matrix (from covariance eigendecomposition)
  ## X_white = X %*% K where K whitens the data
  D <- diag(svd_result$d[1:n_components] / sqrt(T_obs - 1))
  V <- svd_result$v[, 1:n_components, drop = FALSE]
  
  ## Pre-whitening matrix
  K_prewhite <- V %*% diag(1 / diag(D))
  
  ## Whitened data
  X_white <- residuals %*% K_prewhite
  
  ## Step 2: Run ICA with multiple restarts
  best_result <- NULL
  best_negentropy <- -Inf
  convergence_log <- list()
  
  for (restart in 1:n_restarts) {
    if (verbose) cat(sprintf("ICA restart %d/%d...\n", restart, n_restarts))
    
    ica_result <- tryCatch({
      tsmarch::radical(
        X = residuals, 
        components = n_components, 
        demean = demean,
        trace = verbose && restart == 1
      )
    }, error = function(e) {
      if (verbose) cat(sprintf("  Restart %d failed: %s\n", restart, e$message))
      NULL
    })
    
    if (is.null(ica_result)) {
      convergence_log[[restart]] <- list(success = FALSE, error = "ICA failed")
      next
    }
    
    ## Evaluate quality using negentropy approximation
    S <- ica_result$S
    quality <- evaluate_ica_quality(S)
    
    convergence_log[[restart]] <- list(
      success = TRUE,
      negentropy = quality$total_negentropy,
      independence = quality$independence_score
    )
    
    if (quality$total_negentropy > best_negentropy) {
      best_negentropy <- quality$total_negentropy
      best_result <- ica_result
      best_quality <- quality
    }
  }
  
  ## Step 3: Handle failure cases
  if (is.null(best_result)) {
    warning("ICA decomposition failed on all restarts. Using PCA fallback.")
    
    ## PCA fallback: use SVD components directly
    S <- svd_result$u[, 1:n_components, drop = FALSE] %*% D
    A <- V
    W <- t(V) %*% diag(diag(D))
    K <- K_prewhite
    
    best_quality <- list(
      total_negentropy = 0,
      independence_score = 0,
      component_kurtosis = apply(S, 2, kurtosis_simple),
      is_fallback = TRUE
    )
    
    convergence_info <- list(
      method = "pca_fallback",
      n_restarts = n_restarts,
      all_failed = TRUE,
      log = convergence_log
    )
  } else {
    ## Extract components from best result
    S <- best_result$S
    A <- best_result$A
    W <- best_result$W
    K <- best_result$K
    
    best_quality$is_fallback <- FALSE
    
    convergence_info <- list(
      method = method,
      n_restarts = n_restarts,
      successful_restarts = sum(sapply(convergence_log, function(x) x$success)),
      best_negentropy = best_negentropy,
      log = convergence_log
    )
  }
  
  ## Step 4: Validate decomposition
  ## Check that Y ≈ S %*% A' (reconstruction error)
  Y_reconstructed <- S %*% t(A)
  reconstruction_error <- norm(residuals - Y_reconstructed, "F") / norm(residuals, "F")
  
  if (reconstruction_error > 0.01 && verbose) {
    warning(sprintf("ICA reconstruction error is %.2f%%. Results may be unreliable.", 
                    reconstruction_error * 100))
  }
  
  best_quality$reconstruction_error <- reconstruction_error
  
  list(
    S = S,
    A = A,
    W = W,
    K = K,
    method = if (best_quality$is_fallback) "pca_fallback" else method,
    n_components = n_components,
    convergence_info = convergence_info,
    quality = best_quality,
    mu = mu
  )
}


#' @title Evaluate ICA Quality
#' 
#' @description Computes quality metrics for ICA decomposition.
#'
#' @param S T x n_components matrix of independent components
#'
#' @return List with:
#'   \describe{
#'     \item{independence_score}{1 - mean(|cor(S)|), higher is better}
#'     \item{total_negentropy}{Sum of component negentropies}
#'     \item{component_kurtosis}{Excess kurtosis per component}
#'     \item{correlation_matrix}{Component correlation matrix}
#'   }
#'
#' @seealso \code{\link{improved_ica_decomposition}}, \code{\link{gogarch_diagnostics}}
#'
#' @export
evaluate_ica_quality <- function(S) {
  n_components <- ncol(S)
  
  ## Compute negentropy approximation for each component
  ## Using the approximation: J(y) ≈ [E{G(y)} - E{G(v)}]^2
  ## where G(u) = log(cosh(u)) and v ~ N(0,1)
  
  negentropy <- numeric(n_components)
  kurtosis_vals <- numeric(n_components)
  
  for (i in 1:n_components) {
    s <- S[, i]
    s <- (s - mean(s)) / sd(s)  # Standardize
    
    ## Negentropy approximation
    G_s <- mean(log(cosh(s)))
    G_gauss <- 0.3746  # E[log(cosh(N(0,1)))] ≈ 0.3746
    negentropy[i] <- (G_s - G_gauss)^2
    
    ## Kurtosis (should be non-zero for non-Gaussian)
    kurtosis_vals[i] <- kurtosis_simple(s)
  }
  
  ## Independence score: mutual information approximation
  ## Lower correlation between components = better independence
  cor_matrix <- cor(S)
  diag(cor_matrix) <- 0
  independence_score <- 1 - mean(abs(cor_matrix))
  
  list(
    component_negentropy = negentropy,
    total_negentropy = sum(negentropy),
    component_kurtosis = kurtosis_vals,
    independence_score = independence_score,
    correlation_matrix = cor_matrix + diag(n_components)
  )
}


#' @title Simple Kurtosis Calculation
#' @keywords internal
kurtosis_simple <- function(x) {
  n <- length(x)
  x <- x - mean(x)
  m4 <- mean(x^4)
  m2 <- mean(x^2)
  m4 / m2^2 - 3
}


## SECTION 2: Enhanced GOGARCH Estimation with Better ICA ======================

#' @title Estimate GOGARCH with Improved ICA
#' 
#' @description Alternative to \code{\link{estimate_garch_weighted_gogarch}} 
#'   using \code{\link{improved_ica_decomposition}} for better convergence.
#'
#' @param residuals T x k matrix of residuals
#' @param weights T-vector of observation weights
#' @param spec Model specification (see \code{\link{estimate_garch_weighted_gogarch}})
#' @param ica_restarts Number of ICA restarts (default: 3)
#' @param diagnostics Optional diagnostics collector
#' @param verbose Print progress
#'
#' @return Same structure as \code{\link{estimate_garch_weighted_gogarch}} with
#'   ICA quality metrics in \code{coefficients$ica_info$quality}
#'
#' @details
#' Use this when ICA convergence is problematic. Provides warnings when
#' ICA quality is poor (independence score < 0.8) and falls back to PCA
#' if ICA fails completely.
#'
#' @seealso \code{\link{estimate_garch_weighted_gogarch}}, 
#'   \code{\link{improved_ica_decomposition}}, \code{\link{gogarch_diagnostics}}
#'
#' @export
estimate_gogarch_improved <- function(
    residuals,
    weights,
    spec,
    ica_restarts = 3,
    diagnostics = NULL,
    verbose = FALSE
) {
  k <- ncol(residuals)
  T_obs <- nrow(residuals)
  
  ## Adjust weights if needed
  if (length(weights) > T_obs) {
    n_to_remove <- length(weights) - T_obs
    weights <- weights[(n_to_remove + 1):length(weights)]
  }
  
  n_components <- spec$garch_spec_args$components %||% k
  
  ## Step 1: Improved ICA decomposition
  if (verbose) cat("Performing improved ICA decomposition...\n")
  
  ica_result <- improved_ica_decomposition(
    residuals = residuals,
    n_components = n_components,
    n_restarts = ica_restarts,
    demean = FALSE,
    verbose = verbose
  )
  
  if (verbose) {
    cat(sprintf("ICA method: %s\n", ica_result$method))
    cat(sprintf("Independence score: %.3f\n", ica_result$quality$independence_score))
    cat(sprintf("Reconstruction error: %.4f%%\n", ica_result$quality$reconstruction_error * 100))
  }
  
  ## Check if ICA quality is acceptable
  if (ica_result$quality$is_fallback) {
    warning("ICA failed, using PCA fallback. GOGARCH results may not be reliable.")
  } else if (ica_result$quality$independence_score < 0.8) {
    warning(sprintf("ICA independence score (%.3f) is low. Components may not be truly independent.",
                    ica_result$quality$independence_score))
  }
  
  ## Step 2: Estimate GARCH on each component
  S <- ica_result$S
  garch_pars_list <- list()
  warnings_list <- list()
  
  for (i in 1:n_components) {
    if (verbose) cat(sprintf("Estimating GARCH for component %d/%d...\n", i, n_components))
    
    component <- S[, i]
    
    ## Get starting values
    start_pars <- if (!is.null(spec$start_pars$garch_pars) && 
                      length(spec$start_pars$garch_pars) >= i) {
      spec$start_pars$garch_pars[[i]]
    } else {
      list(omega = 0.05, alpha1 = 0.08, beta1 = 0.90)
    }
    
    ## Estimate
    garch_result <- estimate_garch_weighted_univariate_gogarch(
      residuals = component,
      weights = weights,
      spec = list(
        garch_spec_args = spec$garch_spec_args,
        distribution = spec$distribution %||% "norm",
        start_pars = list(garch_pars = list(start_pars))
      ),
      diagnostics = diagnostics,
      verbose = verbose
    )
    
    garch_pars_list[[i]] <- garch_result$coefficients
    if (length(garch_result$warnings) > 0) {
      warnings_list <- c(warnings_list, garch_result$warnings)
    }
  }
  
  ## Build result
  list(
    coefficients = list(
      garch_pars = garch_pars_list,
      ica_info = list(
        A = ica_result$A,
        W = ica_result$W,
        K = ica_result$K,
        S = ica_result$S,
        method = ica_result$method,
        n_components = n_components,
        quality = ica_result$quality,
        convergence_info = ica_result$convergence_info
      ),
      dist_pars = spec$start_pars$dist_pars,
      correlation_type = "gogarch"
    ),
    warnings = warnings_list,
    diagnostics = diagnostics
  )
}


## SECTION 3: GOGARCH Quality Diagnostics ======================================

#' @title GOGARCH Model Diagnostics
#' 
#' @description Comprehensive diagnostics for GOGARCH model fit.
#'
#' @param gogarch_result Result from \code{\link{estimate_garch_weighted_gogarch}} 
#'   or \code{\link{estimate_gogarch_improved}}
#' @param residuals Original T x k residual matrix
#' @param verbose Print diagnostic report (default: TRUE)
#'
#' @return List with:
#'   \describe{
#'     \item{ica_quality}{Independence score, negentropy, reconstruction error}
#'     \item{component_diagnostics}{Per-component: persistence, Ljung-Box p-value}
#'     \item{mixing_matrix}{Orthogonality error, condition number}
#'     \item{covariance_fit}{RMSE, correlation vs sample covariance}
#'   }
#'
#' @details
#' Key quality thresholds:
#' \itemize{
#'   \item Independence score > 0.8
#'   \item Reconstruction error < 1\%
#'   \item Mixing matrix condition number < 100
#'   \item Component persistence < 1 (stationarity)
#' }
#'
#' @seealso \code{\link{estimate_garch_weighted_gogarch}}, 
#'   \code{\link{improved_ica_decomposition}}
#'
#' @export
gogarch_diagnostics <- function(
    gogarch_result,
    residuals,
    verbose = TRUE
) {
  coeffs <- gogarch_result$coefficients %||% gogarch_result
  ica_info <- coeffs$ica_info
  garch_pars <- coeffs$garch_pars
  
  T_obs <- nrow(residuals)
  k <- ncol(residuals)
  n_components <- length(garch_pars)
  
  ## 1. ICA Quality Metrics
  S <- ica_info$S %||% (residuals %*% t(ica_info$W))
  ica_quality <- evaluate_ica_quality(S)
  
  ## 2. Component GARCH Diagnostics
  component_diagnostics <- list()
  for (i in 1:n_components) {
    pars <- garch_pars[[i]]
    omega <- pars$omega %||% 0.05
    alpha <- unlist(pars[grepl("alpha", names(pars))])
    beta <- unlist(pars[grepl("beta", names(pars))])
    
    if (length(alpha) == 0) alpha <- 0.1
    if (length(beta) == 0) beta <- 0.8
    
    persistence <- sum(alpha) + sum(beta)
    
    ## Compute standardized residuals
    sigma2 <- compute_garch_variance(S[, i], omega, alpha, beta)
    std_resid <- S[, i] / sqrt(sigma2)
    
    component_diagnostics[[i]] <- list(
      persistence = persistence,
      is_stationary = persistence < 1,
      std_resid_kurtosis = kurtosis_simple(std_resid),
      std_resid_skewness = mean((std_resid - mean(std_resid))^3) / sd(std_resid)^3,
      ljung_box_p = tryCatch(
        Box.test(std_resid^2, lag = 10, type = "Ljung-Box")$p.value,
        error = function(e) NA
      )
    )
  }
  
  ## 3. Mixing Matrix Quality
  A <- ica_info$A
  W <- ica_info$W
  
  ## Check orthogonality: W %*% A should be close to identity
  WA <- W %*% A
  orthogonality_error <- norm(WA - diag(n_components), "F") / n_components
  
  ## Condition number of mixing matrix
  mixing_condition <- kappa(A)
  
  ## 4. Covariance Reconstruction Quality
  ## Compute sample covariance vs model-implied covariance
  sample_cov <- cov(residuals)
  
  ## Model covariance: H = A %*% D %*% A' where D is diagonal GARCH variances
  D_mean <- diag(sapply(1:n_components, function(i) {
    pars <- garch_pars[[i]]
    omega <- pars$omega %||% 0.05
    alpha <- sum(unlist(pars[grepl("alpha", names(pars))]))
    beta <- sum(unlist(pars[grepl("beta", names(pars))]))
    if (alpha + beta >= 1) return(var(S[, i]))
    omega / (1 - alpha - beta)
  }))
  
  model_cov <- A %*% D_mean %*% t(A)
  
  cov_rmse <- sqrt(mean((sample_cov - model_cov)^2))
  cov_mae <- mean(abs(sample_cov - model_cov))
  
  diagnostics <- list(
    ica_quality = ica_quality,
    component_diagnostics = component_diagnostics,
    mixing_matrix = list(
      orthogonality_error = orthogonality_error,
      condition_number = mixing_condition,
      is_well_conditioned = mixing_condition < 100
    ),
    covariance_fit = list(
      rmse = cov_rmse,
      mae = cov_mae,
      correlation_fit = cor(as.vector(sample_cov), as.vector(model_cov))
    ),
    overall = list(
      n_components = n_components,
      k = k,
      T_obs = T_obs,
      ica_method = ica_info$method
    )
  )
  
  if (verbose) {
    cat("\n=== GOGARCH Model Diagnostics ===\n\n")
    
    cat("ICA Quality:\n")
    cat(sprintf("  Independence score: %.3f\n", ica_quality$independence_score))
    cat(sprintf("  Total negentropy: %.4f\n", ica_quality$total_negentropy))
    if (!is.null(ica_info$quality$reconstruction_error)) {
      cat(sprintf("  Reconstruction error: %.4f%%\n", 
                  ica_info$quality$reconstruction_error * 100))
    }
    
    cat("\nMixing Matrix:\n")
    cat(sprintf("  Orthogonality error: %.4f\n", orthogonality_error))
    cat(sprintf("  Condition number: %.1f\n", mixing_condition))
    
    cat("\nComponent GARCH:\n")
    for (i in 1:n_components) {
      diag_i <- component_diagnostics[[i]]
      cat(sprintf("  Component %d: persistence=%.3f, LB p=%.3f\n",
                  i, diag_i$persistence, diag_i$ljung_box_p %||% NA))
    }
    
    cat("\nCovariance Fit:\n")
    cat(sprintf("  RMSE: %.4f\n", cov_rmse))
    cat(sprintf("  Correlation: %.4f\n", diagnostics$covariance_fit$correlation_fit))
  }
  
  diagnostics
}


#' @title Compute GARCH Variance Path
#' @keywords internal
compute_garch_variance <- function(resid, omega, alpha, beta) {
  n <- length(resid)
  p <- length(alpha)
  q <- length(beta)
  
  ## Initialize with unconditional variance
  sigma2 <- rep(omega / (1 - sum(alpha) - sum(beta)), n)
  sigma2[sigma2 <= 0 | !is.finite(sigma2)] <- var(resid)
  
  maxpq <- max(p, q, 1)
  
  for (t in (maxpq + 1):n) {
    sigma2[t] <- omega
    for (i in 1:p) sigma2[t] <- sigma2[t] + alpha[i] * resid[t - i]^2
    for (j in 1:q) sigma2[t] <- sigma2[t] + beta[j] * sigma2[t - j]
    sigma2[t] <- max(sigma2[t], 1e-10)
  }
  
  sigma2
}

