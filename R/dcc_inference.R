## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## DCC(1,1) Inference: Bootstrap SEs and Profile Likelihood CIs
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
##
## This file provides alternative inference methods for DCC parameters that
## do not rely on Hessian-based standard errors:
##
##   1. Bootstrap Standard Errors
##      - Parametric bootstrap (resample from fitted model)
##      - Residual bootstrap (resample standardized residuals)
##      - BCa intervals for bias correction
##
##   2. Profile Likelihood Confidence Intervals
##      - Does not assume quadratic likelihood
##      - Based on likelihood ratio test inversion
##      - Handles asymmetric uncertainty correctly
##
##   3. Comprehensive Inference
##      - Computes all three methods (Hessian, Bootstrap, Profile)
##      - Provides comparison and recommendations
##
## WHY THESE METHODS ARE NEEDED
## ────────────────────────────
## Validation studies (see validate_dcc_vs_tsmarch.R) demonstrate that:
##
##   1. The DCC likelihood is ~8x flatter in beta than alpha direction
##   2. Hessian-based SE for beta is only ~16% of bootstrap SE
##   3. This is a fundamental property of DCC, not an implementation error
##   4. Our implementation matches tsmarch to 6 decimal places
##
## Quantitative findings (n=1000, α=0.05, β=0.90):
##   ┌─────────────────────────────────────────────────────┐
##   │ Parameter │ Hessian SE │ Bootstrap SE │ Ratio      │
##   ├─────────────────────────────────────────────────────┤
##   │ Alpha     │ 0.019      │ 0.025        │ 0.78       │
##   │ Beta      │ 0.055      │ 0.343        │ 0.16       │
##   └─────────────────────────────────────────────────────┘
##
## RECOMMENDATIONS
## ───────────────
##   - Alpha: Hessian-based SE is acceptable
##   - Beta:  ALWAYS use bootstrap or profile likelihood
##   - For publication: Report bootstrap SE and/or profile CI for beta
##
## See DCC_INFERENCE_DOCUMENTATION.R for complete findings and guidelines.
##
## Dependencies:
##   - dcc_gradient.R (NLL computation)
##   - dcc_hessian.R (Hessian-based SE for comparison)
##   - diagnostic_utils.R (simulate_dcc_garch for parametric bootstrap)
##
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


## SECTION 1: Bootstrap Standard Errors ========================================

#' @title Bootstrap Standard Errors for DCC(1,1) Parameters
#' @description Compute bootstrap standard errors using either parametric or
#'   residual resampling.
#' @param std_resid T x k matrix of standardized residuals
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional correlation matrix
#' @param mle_params MLE estimates c(alpha, beta)
#' @param n_boot Number of bootstrap replications (default 200)
#' @param method Bootstrap method: "parametric" or "residual"
#' @param distribution "mvn" or "mvt"
#' @param shape Shape parameter for MVT (required if distribution="mvt")
#' @param verbose Print progress
#' @param seed Random seed for reproducibility
#' @return List with:
#'   \item{se}{Bootstrap standard errors}
#'   \item{boot_estimates}{Matrix of bootstrap estimates (n_boot x 2)}
#'   \item{bias}{Bootstrap estimate of bias}
#'   \item{ci_percentile}{Percentile confidence intervals}
#'   \item{ci_basic}{Basic bootstrap confidence intervals}
#'   \item{ci_bca}{BCa confidence intervals (if available)}
#' @export
dcc11_bootstrap_se <- function(
    std_resid,
    weights,
    Qbar,
    mle_params,
    n_boot = 200,
    method = c("residual", "parametric"),
    distribution = "mvn",
    shape = NULL,
    verbose = TRUE,
    seed = NULL
) {
  
  method <- match.arg(method)
  
  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(std_resid)
  k <- ncol(std_resid)
  
  alpha_mle <- mle_params[1]
  beta_mle <- mle_params[2]
  
  ## Storage for bootstrap estimates
  boot_estimates <- matrix(NA, n_boot, 2)
  colnames(boot_estimates) <- c("alpha", "beta")
  
  if (verbose) {
    cat(sprintf("Running %s bootstrap with %d replications...\n", method, n_boot))
    pb <- txtProgressBar(min = 0, max = n_boot, style = 3)
  }
  
  for (b in 1:n_boot) {
    
    if (method == "residual") {
      ## Residual bootstrap: resample rows of standardized residuals
      boot_idx <- sample(1:n, n, replace = TRUE)
      boot_resid <- std_resid[boot_idx, , drop = FALSE]
      boot_weights <- weights[boot_idx]
      boot_Qbar <- cor(boot_resid)
      
    } else {
      ## Parametric bootstrap: simulate from fitted DCC model
      boot_resid <- simulate_dcc_residuals(
        n = n,
        k = k,
        alpha = alpha_mle,
        beta = beta_mle,
        Qbar = Qbar,
        distribution = distribution,
        shape = shape
      )
      boot_weights <- weights
      boot_Qbar <- Qbar  ## Use original Qbar, not sample correlation
    }
    
    ## Re-estimate DCC parameters
    opt_result <- tryCatch({
      optim(
        par = mle_params,
        fn = dcc11_nll,
        method = "L-BFGS-B",
        lower = c(1e-6, 1e-6),
        upper = c(0.5, 0.999),
        std_resid = boot_resid,
        weights = boot_weights,
        Qbar = boot_Qbar,
        distribution = distribution,
        use_reparam = FALSE
      )
    }, error = function(e) NULL)
    
    if (!is.null(opt_result) && opt_result$convergence == 0) {
      boot_estimates[b, ] <- opt_result$par
    }
    
    if (verbose) setTxtProgressBar(pb, b)
  }
  
  if (verbose) {
    close(pb)
    cat("\n")
  }
  
  ## Remove failed bootstrap replications
  valid_idx <- complete.cases(boot_estimates)
  n_valid <- sum(valid_idx)
  boot_estimates_valid <- boot_estimates[valid_idx, , drop = FALSE]
  
  if (n_valid < n_boot * 0.8) {
    warning(sprintf("Only %d/%d bootstrap replications succeeded", n_valid, n_boot))
  }
  
  ## Compute bootstrap statistics
  boot_mean <- colMeans(boot_estimates_valid)
  boot_se <- apply(boot_estimates_valid, 2, sd)
  boot_bias <- boot_mean - mle_params
  
  ## Percentile confidence intervals (95%)
  ci_percentile <- apply(boot_estimates_valid, 2, quantile, probs = c(0.025, 0.975))
  
  ## Basic bootstrap CI: 2*theta_hat - quantile
  ci_basic <- matrix(NA, 2, 2)
  rownames(ci_basic) <- c("2.5%", "97.5%")
  colnames(ci_basic) <- c("alpha", "beta")
  ci_basic[1, ] <- 2 * mle_params - ci_percentile[2, ]
  ci_basic[2, ] <- 2 * mle_params - ci_percentile[1, ]
  
  ## BCa confidence intervals (bias-corrected and accelerated)
  ci_bca <- tryCatch({
    compute_bca_ci(boot_estimates_valid, mle_params, std_resid, weights, Qbar, distribution)
  }, error = function(e) {
    if (verbose) cat("BCa CI computation failed, skipping\n")
    NULL
  })
  
  if (verbose) {
    cat(sprintf("Valid bootstrap replications: %d/%d\n", n_valid, n_boot))
    cat(sprintf("Bootstrap SE: alpha=%.4f, beta=%.4f\n", boot_se[1], boot_se[2]))
  }
  
  list(
    se = boot_se,
    boot_estimates = boot_estimates,
    n_valid = n_valid,
    mean = boot_mean,
    bias = boot_bias,
    ci_percentile = ci_percentile,
    ci_basic = ci_basic,
    ci_bca = ci_bca,
    method = method,
    n_boot = n_boot
  )
}


#' @title Simulate DCC Standardized Residuals
#' @description Simulate standardized residuals from a DCC(1,1) model.
#'   Used for parametric bootstrap.
#' @keywords internal
simulate_dcc_residuals <- function(
    n,
    k,
    alpha,
    beta,
    Qbar,
    distribution = "mvn",
    shape = NULL
) {
  
  z <- matrix(0, n, k)
  Q <- Qbar
  
  for (t in 1:n) {
    ## Normalize Q to correlation R
    d <- sqrt(diag(Q))
    d[d < 1e-8] <- 1e-8
    D_inv <- diag(1 / d, k)
    R <- D_inv %*% Q %*% D_inv
    
    ## Ensure PD
    R <- (R + t(R)) / 2
    eig <- eigen(R, symmetric = TRUE)
    if (any(eig$values < 1e-8)) {
      eig$values[eig$values < 1e-8] <- 1e-8
      R <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
      d <- sqrt(diag(R))
      D_inv <- diag(1 / d, k)
      R <- D_inv %*% R %*% D_inv
    }
    
    ## Draw from multivariate distribution
    if (distribution == "mvn") {
      L <- chol(R)
      eps <- rnorm(k)
      z[t, ] <- as.vector(t(L) %*% eps)
    } else {
      ## MVT
      if (is.null(shape)) shape <- 8
      L <- chol(R)
      eps <- rnorm(k)
      chi_sq <- rchisq(1, df = shape)
      z[t, ] <- as.vector(t(L) %*% eps) * sqrt(shape / chi_sq)
    }
    
    ## Update Q for next period
    if (t < n) {
      z_t <- z[t, , drop = FALSE]
      Q <- (1 - alpha - beta) * Qbar + alpha * (t(z_t) %*% z_t) + beta * Q
    }
  }
  
  z
}


#' @title Compute BCa Confidence Intervals
#' @description Bias-corrected and accelerated bootstrap confidence intervals.
#' @keywords internal
compute_bca_ci <- function(
    boot_estimates,
    mle_params,
    std_resid,
    weights,
    Qbar,
    distribution,
    alpha_level = 0.05
) {
  
  n_boot <- nrow(boot_estimates)
  n <- nrow(std_resid)
  
  ci_bca <- matrix(NA, 2, 2)
  rownames(ci_bca) <- c(paste0(alpha_level/2 * 100, "%"), 
                        paste0((1 - alpha_level/2) * 100, "%"))
  colnames(ci_bca) <- c("alpha", "beta")
  
  for (j in 1:2) {
    theta_hat <- mle_params[j]
    theta_boot <- boot_estimates[, j]
    
    ## Bias correction factor z0
    z0 <- qnorm(mean(theta_boot < theta_hat))
    
    ## Acceleration factor (jackknife)
    theta_jack <- numeric(n)
    for (i in 1:n) {
      jack_resid <- std_resid[-i, , drop = FALSE]
      jack_weights <- weights[-i]
      jack_Qbar <- cor(jack_resid)
      
      jack_opt <- tryCatch({
        optim(
          par = mle_params,
          fn = dcc11_nll,
          method = "L-BFGS-B",
          lower = c(1e-6, 1e-6),
          upper = c(0.5, 0.999),
          std_resid = jack_resid,
          weights = jack_weights,
          Qbar = jack_Qbar,
          distribution = distribution,
          use_reparam = FALSE
        )
      }, error = function(e) NULL)
      
      if (!is.null(jack_opt) && jack_opt$convergence == 0) {
        theta_jack[i] <- jack_opt$par[j]
      } else {
        theta_jack[i] <- NA
      }
    }
    
    theta_jack <- theta_jack[!is.na(theta_jack)]
    theta_jack_mean <- mean(theta_jack)
    
    ## Acceleration
    num <- sum((theta_jack_mean - theta_jack)^3)
    denom <- 6 * sum((theta_jack_mean - theta_jack)^2)^1.5
    a <- if (abs(denom) > 1e-10) num / denom else 0
    
    ## Adjusted quantiles
    z_alpha_lower <- qnorm(alpha_level / 2)
    z_alpha_upper <- qnorm(1 - alpha_level / 2)
    
    alpha1 <- pnorm(z0 + (z0 + z_alpha_lower) / (1 - a * (z0 + z_alpha_lower)))
    alpha2 <- pnorm(z0 + (z0 + z_alpha_upper) / (1 - a * (z0 + z_alpha_upper)))
    
    ci_bca[1, j] <- quantile(theta_boot, alpha1, na.rm = TRUE)
    ci_bca[2, j] <- quantile(theta_boot, alpha2, na.rm = TRUE)
  }
  
  ci_bca
}


## SECTION 2: Profile Likelihood Confidence Intervals ==========================

#' @title Profile Likelihood Confidence Intervals for DCC(1,1)
#' @description Compute profile likelihood confidence intervals by inverting
#'   the likelihood ratio test. Does not assume quadratic likelihood.
#' @param std_resid T x k matrix of standardized residuals
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional correlation matrix
#' @param mle_params MLE estimates c(alpha, beta)
#' @param mle_nll NLL at MLE (optional, computed if not provided)
#' @param conf_level Confidence level (default 0.95)
#' @param param Which parameter to profile: "alpha", "beta", or "both"
#' @param n_points Number of points for profile (default 50)
#' @param distribution "mvn" or "mvt"
#' @param verbose Print progress
#' @return List with profile likelihood results and confidence intervals
#' @export
dcc11_profile_likelihood_ci <- function(
    std_resid,
    weights,
    Qbar,
    mle_params,
    mle_nll = NULL,
    conf_level = 0.95,
    param = c("both", "alpha", "beta"),
    n_points = 50,
    distribution = "mvn",
    verbose = TRUE
) {
  
  param <- match.arg(param)
  
  alpha_mle <- mle_params[1]
  beta_mle <- mle_params[2]
  
  ## Compute MLE NLL if not provided
  if (is.null(mle_nll)) {
    mle_nll <- dcc11_nll(
      params = mle_params,
      std_resid = std_resid,
      weights = weights,
      Qbar = Qbar,
      distribution = distribution,
      use_reparam = FALSE
    )
  }
  
  ## Critical value for LR test (chi-squared with 1 df)
  chi_sq_crit <- qchisq(conf_level, df = 1)
  nll_threshold <- mle_nll + chi_sq_crit / 2
  
  results <- list()
  
  ## Profile alpha
  if (param %in% c("both", "alpha")) {
    if (verbose) cat("Profiling alpha...\n")
    
    profile_alpha <- profile_one_parameter(
      param_idx = 1,
      param_name = "alpha",
      mle_params = mle_params,
      mle_nll = mle_nll,
      nll_threshold = nll_threshold,
      std_resid = std_resid,
      weights = weights,
      Qbar = Qbar,
      distribution = distribution,
      n_points = n_points
    )
    
    results$alpha <- profile_alpha
  }
  
  ## Profile beta
  if (param %in% c("both", "beta")) {
    if (verbose) cat("Profiling beta...\n")
    
    profile_beta <- profile_one_parameter(
      param_idx = 2,
      param_name = "beta",
      mle_params = mle_params,
      mle_nll = mle_nll,
      nll_threshold = nll_threshold,
      std_resid = std_resid,
      weights = weights,
      Qbar = Qbar,
      distribution = distribution,
      n_points = n_points
    )
    
    results$beta <- profile_beta
  }
  
  results$mle_params <- mle_params
  results$mle_nll <- mle_nll
  results$conf_level <- conf_level
  results$nll_threshold <- nll_threshold
  
  if (verbose) {
    cat("\nProfile Likelihood CIs:\n")
    if (!is.null(results$alpha)) {
      cat(sprintf("  alpha: [%.4f, %.4f]\n", 
                  results$alpha$ci[1], results$alpha$ci[2]))
    }
    if (!is.null(results$beta)) {
      cat(sprintf("  beta:  [%.4f, %.4f]\n", 
                  results$beta$ci[1], results$beta$ci[2]))
    }
  }
  
  results
}


#' @title Profile One Parameter
#' @description Internal function to compute profile likelihood for one parameter.
#' @keywords internal
profile_one_parameter <- function(
    param_idx,
    param_name,
    mle_params,
    mle_nll,
    nll_threshold,
    std_resid,
    weights,
    Qbar,
    distribution,
    n_points
) {
  
  other_idx <- if (param_idx == 1) 2 else 1
  mle_value <- mle_params[param_idx]
  
  ## Determine range to search
  if (param_name == "alpha") {
    ## Alpha typically small, search from near 0 to ~0.3
    search_min <- 1e-6
    search_max <- min(0.4, 0.999 - mle_params[2] - 0.01)  ## Stationarity
  } else {
    ## Beta typically large, search from ~0.5 to near 1
    search_min <- max(0.01, 1e-6)
    search_max <- min(0.999, 0.999 - mle_params[1] - 0.01)  ## Stationarity
  }
  
  ## Create grid of values to profile over
  ## Denser near MLE, sparser at edges
  grid_lower <- seq(search_min, mle_value, length.out = n_points / 2)
  grid_upper <- seq(mle_value, search_max, length.out = n_points / 2)
  param_grid <- unique(c(grid_lower, grid_upper))
  param_grid <- sort(param_grid)
  
  ## Compute profile likelihood at each grid point
  profile_nll <- numeric(length(param_grid))
  profile_other <- numeric(length(param_grid))
  
  for (i in seq_along(param_grid)) {
    fixed_value <- param_grid[i]
    
    ## Optimize over the other parameter
    opt_result <- tryCatch({
      if (param_idx == 1) {
        ## Fixed alpha, optimize beta
        optimize(
          f = function(beta) {
            if (fixed_value + beta >= 0.999) return(1e10)
            dcc11_nll(
              params = c(fixed_value, beta),
              std_resid = std_resid,
              weights = weights,
              Qbar = Qbar,
              distribution = distribution,
              use_reparam = FALSE
            )
          },
          interval = c(1e-6, 0.999 - fixed_value - 1e-6),
          maximum = FALSE
        )
      } else {
        ## Fixed beta, optimize alpha
        optimize(
          f = function(alpha) {
            if (alpha + fixed_value >= 0.999) return(1e10)
            dcc11_nll(
              params = c(alpha, fixed_value),
              std_resid = std_resid,
              weights = weights,
              Qbar = Qbar,
              distribution = distribution,
              use_reparam = FALSE
            )
          },
          interval = c(1e-6, 0.999 - fixed_value - 1e-6),
          maximum = FALSE
        )
      }
    }, error = function(e) list(objective = NA, minimum = NA))
    
    profile_nll[i] <- opt_result$objective
    profile_other[i] <- opt_result$minimum
  }
  
  ## Find CI bounds by interpolation
  valid_idx <- !is.na(profile_nll)
  param_valid <- param_grid[valid_idx]
  nll_valid <- profile_nll[valid_idx]
  
  ## Lower bound: largest value below MLE where NLL crosses threshold
  lower_idx <- which(param_valid < mle_value)
  if (length(lower_idx) > 0) {
    lower_nll <- nll_valid[lower_idx]
    lower_param <- param_valid[lower_idx]
    
    ## Find where NLL crosses threshold
    crosses_lower <- which(lower_nll > nll_threshold)
    if (length(crosses_lower) > 0) {
      cross_idx <- max(crosses_lower)
      if (cross_idx < length(lower_idx)) {
        ## Linear interpolation
        x1 <- lower_param[cross_idx]
        x2 <- lower_param[cross_idx + 1]
        y1 <- lower_nll[cross_idx]
        y2 <- lower_nll[cross_idx + 1]
        ci_lower <- x1 + (nll_threshold - y1) * (x2 - x1) / (y2 - y1)
      } else {
        ci_lower <- lower_param[cross_idx]
      }
    } else {
      ci_lower <- min(lower_param)  ## Didn't cross, use grid min
    }
  } else {
    ci_lower <- search_min
  }
  
  ## Upper bound: smallest value above MLE where NLL crosses threshold
  upper_idx <- which(param_valid > mle_value)
  if (length(upper_idx) > 0) {
    upper_nll <- nll_valid[upper_idx]
    upper_param <- param_valid[upper_idx]
    
    crosses_upper <- which(upper_nll > nll_threshold)
    if (length(crosses_upper) > 0) {
      cross_idx <- min(crosses_upper)
      if (cross_idx > 1) {
        ## Linear interpolation
        x1 <- upper_param[cross_idx - 1]
        x2 <- upper_param[cross_idx]
        y1 <- upper_nll[cross_idx - 1]
        y2 <- upper_nll[cross_idx]
        ci_upper <- x1 + (nll_threshold - y1) * (x2 - x1) / (y2 - y1)
      } else {
        ci_upper <- upper_param[cross_idx]
      }
    } else {
      ci_upper <- max(upper_param)  ## Didn't cross, use grid max
    }
  } else {
    ci_upper <- search_max
  }
  
  list(
    param_name = param_name,
    mle = mle_value,
    ci = c(lower = ci_lower, upper = ci_upper),
    grid = param_grid,
    profile_nll = profile_nll,
    profile_other = profile_other,
    nll_threshold = nll_threshold
  )
}


#' @title Plot Profile Likelihood
#' @description Create a visualization of the profile likelihood with CI.
#' @param profile_result Result from profile_one_parameter or dcc11_profile_likelihood_ci
#' @param param_name Which parameter to plot ("alpha" or "beta")
#' @return plotly object
#' @export
plot_profile_likelihood <- function(profile_result, param_name = "alpha") {
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required")
  }
  
  if (!is.null(profile_result[[param_name]])) {
    prof <- profile_result[[param_name]]
  } else {
    prof <- profile_result
  }
  
  valid_idx <- !is.na(prof$profile_nll)
  
  p <- plotly::plot_ly() %>%
    ## Profile likelihood curve
    plotly::add_trace(
      x = prof$grid[valid_idx],
      y = prof$profile_nll[valid_idx],
      type = "scatter",
      mode = "lines",
      line = list(color = "blue", width = 2),
      name = "Profile NLL"
    ) %>%
    ## Threshold line
    plotly::add_trace(
      x = range(prof$grid[valid_idx]),
      y = c(prof$nll_threshold, prof$nll_threshold),
      type = "scatter",
      mode = "lines",
      line = list(color = "red", width = 2, dash = "dash"),
      name = sprintf("%.0f%% threshold", profile_result$conf_level * 100)
    ) %>%
    ## MLE point
    plotly::add_trace(
      x = prof$mle,
      y = profile_result$mle_nll,
      type = "scatter",
      mode = "markers",
      marker = list(color = "green", size = 12, symbol = "diamond"),
      name = sprintf("MLE = %.4f", prof$mle)
    ) %>%
    ## CI bounds
    plotly::add_trace(
      x = c(prof$ci["lower"], prof$ci["lower"]),
      y = c(min(prof$profile_nll, na.rm = TRUE), prof$nll_threshold),
      type = "scatter",
      mode = "lines",
      line = list(color = "orange", width = 1, dash = "dot"),
      name = sprintf("Lower = %.4f", prof$ci["lower"])
    ) %>%
    plotly::add_trace(
      x = c(prof$ci["upper"], prof$ci["upper"]),
      y = c(min(prof$profile_nll, na.rm = TRUE), prof$nll_threshold),
      type = "scatter",
      mode = "lines",
      line = list(color = "orange", width = 1, dash = "dot"),
      name = sprintf("Upper = %.4f", prof$ci["upper"])
    ) %>%
    plotly::layout(
      title = sprintf("Profile Likelihood: %s", prof$param_name),
      xaxis = list(title = prof$param_name),
      yaxis = list(title = "Profile NLL"),
      showlegend = TRUE
    )
  
  p
}


## SECTION 3: Comprehensive Inference Summary ==================================

#' @title Comprehensive DCC Inference
#' @description Compute all three types of inference for DCC parameters:
#'   Hessian-based SEs, bootstrap SEs, and profile likelihood CIs.
#' @param std_resid T x k matrix of standardized residuals
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional correlation matrix
#' @param mle_params MLE estimates c(alpha, beta) (optional, computed if NULL)
#' @param distribution "mvn" or "mvt"
#' @param n_boot Number of bootstrap replications
#' @param boot_method Bootstrap method: "residual" or "parametric"
#' @param n_profile_points Points for profile likelihood
#' @param conf_level Confidence level
#' @param verbose Print progress
#' @param seed Random seed
#' @return List with all inference results
#' @export
dcc11_comprehensive_inference <- function(
    std_resid,
    weights,
    Qbar,
    mle_params = NULL,
    distribution = "mvn",
    n_boot = 200,
    boot_method = "residual",
    n_profile_points = 50,
    conf_level = 0.95,
    verbose = TRUE,
    seed = NULL
) {
  
  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(std_resid)
  
  ## Step 1: Find MLE if not provided
  if (is.null(mle_params)) {
    if (verbose) cat("Finding MLE...\n")
    opt_result <- optim(
      par = c(0.05, 0.90),
      fn = dcc11_nll,
      method = "L-BFGS-B",
      lower = c(1e-6, 1e-6),
      upper = c(0.5, 0.999),
      std_resid = std_resid,
      weights = weights,
      Qbar = Qbar,
      distribution = distribution,
      use_reparam = FALSE
    )
    mle_params <- opt_result$par
    mle_nll <- opt_result$value
  } else {
    mle_nll <- dcc11_nll(
      params = mle_params,
      std_resid = std_resid,
      weights = weights,
      Qbar = Qbar,
      distribution = distribution,
      use_reparam = FALSE
    )
  }
  
  names(mle_params) <- c("alpha", "beta")
  
  if (verbose) {
    cat(sprintf("MLE: alpha=%.4f, beta=%.4f\n", mle_params[1], mle_params[2]))
    cat("\n")
  }
  
  ## Step 2: Hessian-based SEs
  if (verbose) cat("Computing Hessian-based SEs...\n")
  hessian_result <- tryCatch(
    withCallingHandlers({
      dcc11_standard_errors(
        params = mle_params,
        std_resid = std_resid,
        weights = weights,
        Qbar = Qbar,
        distribution = distribution,
        use_reparam = FALSE
      )
    }, warning = function(w) {
      invokeRestart("muffleWarning")
    }),
    error = function(e) list(se = c(NA, NA), warning = e$message)
  )
  
  z_crit <- qnorm(1 - (1 - conf_level) / 2)
  hessian_ci <- matrix(NA, 2, 2)
  rownames(hessian_ci) <- c("lower", "upper")
  colnames(hessian_ci) <- c("alpha", "beta")
  if (all(is.finite(hessian_result$se))) {
    hessian_ci[1, ] <- mle_params - z_crit * hessian_result$se
    hessian_ci[2, ] <- mle_params + z_crit * hessian_result$se
  }
  
  ## Step 3: Bootstrap SEs
  if (verbose) cat("\nComputing bootstrap SEs...\n")
  boot_result <- dcc11_bootstrap_se(
    std_resid = std_resid,
    weights = weights,
    Qbar = Qbar,
    mle_params = mle_params,
    n_boot = n_boot,
    method = boot_method,
    distribution = distribution,
    verbose = verbose,
    seed = seed
  )
  
  ## Step 4: Profile likelihood CIs
  if (verbose) cat("\nComputing profile likelihood CIs...\n")
  profile_result <- dcc11_profile_likelihood_ci(
    std_resid = std_resid,
    weights = weights,
    Qbar = Qbar,
    mle_params = mle_params,
    mle_nll = mle_nll,
    conf_level = conf_level,
    param = "both",
    n_points = n_profile_points,
    distribution = distribution,
    verbose = verbose
  )
  
  ## Compile summary
  summary_table <- data.frame(
    Parameter = c("alpha", "beta"),
    MLE = mle_params,
    Hessian_SE = hessian_result$se,
    Hessian_CI_lower = hessian_ci[1, ],
    Hessian_CI_upper = hessian_ci[2, ],
    Boot_SE = boot_result$se,
    Boot_CI_lower = boot_result$ci_percentile[1, ],
    Boot_CI_upper = boot_result$ci_percentile[2, ],
    Profile_CI_lower = c(profile_result$alpha$ci[1], profile_result$beta$ci[1]),
    Profile_CI_upper = c(profile_result$alpha$ci[2], profile_result$beta$ci[2]),
    row.names = NULL
  )
  
  if (verbose) {
    cat("\n")
    cat("=== Inference Summary ===\n")
    cat(sprintf("Confidence level: %.0f%%\n\n", conf_level * 100))
    print(summary_table, digits = 4)
    
    cat("\nCI Width Comparison:\n")
    for (param in c("alpha", "beta")) {
      hess_width <- hessian_ci[2, param] - hessian_ci[1, param]
      boot_width <- boot_result$ci_percentile[2, param] - boot_result$ci_percentile[1, param]
      prof_ci <- if (param == "alpha") profile_result$alpha$ci else profile_result$beta$ci
      prof_width <- prof_ci[2] - prof_ci[1]
      
      cat(sprintf("  %s: Hessian=%.4f, Bootstrap=%.4f, Profile=%.4f\n",
                  param, hess_width, boot_width, prof_width))
    }
  }
  
  list(
    mle = mle_params,
    mle_nll = mle_nll,
    hessian = list(
      se = hessian_result$se,
      ci = hessian_ci,
      info_matrix = hessian_result$info_matrix
    ),
    bootstrap = boot_result,
    profile = profile_result,
    summary = summary_table,
    conf_level = conf_level,
    settings = list(
      n = n,
      n_boot = n_boot,
      boot_method = boot_method,
      n_profile_points = n_profile_points,
      distribution = distribution
    )
  )
}


#' @title Print Inference Comparison
#' @description Print a formatted comparison of inference methods.
#' @param inf_result Result from dcc11_comprehensive_inference
#' @export
print_inference_comparison <- function(inf_result) {
  
  cat("\n")
  cat("================================================================\n")
  cat("  DCC(1,1) Inference Comparison\n")
  cat("================================================================\n\n")
  
  cat(sprintf("Sample size: %d\n", inf_result$settings$n))
  cat(sprintf("MLE: alpha = %.4f, beta = %.4f\n", 
              inf_result$mle[1], inf_result$mle[2]))
  cat(sprintf("Confidence level: %.0f%%\n\n", inf_result$conf_level * 100))
  
  cat("Standard Errors:\n")
  cat(sprintf("  %-10s %12s %12s\n", "Method", "Alpha", "Beta"))
  cat(sprintf("  %-10s %12.4f %12.4f\n", "Hessian", 
              inf_result$hessian$se[1], inf_result$hessian$se[2]))
  cat(sprintf("  %-10s %12.4f %12.4f\n", "Bootstrap", 
              inf_result$bootstrap$se[1], inf_result$bootstrap$se[2]))
  cat(sprintf("  %-10s %12.2f %12.2f\n", "Ratio", 
              inf_result$hessian$se[1] / inf_result$bootstrap$se[1],
              inf_result$hessian$se[2] / inf_result$bootstrap$se[2]))
  cat("\n")
  
  cat("Confidence Intervals:\n")
  cat(sprintf("  %-10s %20s %20s\n", "Method", "Alpha", "Beta"))
  cat(sprintf("  %-10s [%8.4f, %8.4f] [%8.4f, %8.4f]\n", "Hessian",
              inf_result$hessian$ci[1, 1], inf_result$hessian$ci[2, 1],
              inf_result$hessian$ci[1, 2], inf_result$hessian$ci[2, 2]))
  cat(sprintf("  %-10s [%8.4f, %8.4f] [%8.4f, %8.4f]\n", "Bootstrap",
              inf_result$bootstrap$ci_percentile[1, 1], 
              inf_result$bootstrap$ci_percentile[2, 1],
              inf_result$bootstrap$ci_percentile[1, 2], 
              inf_result$bootstrap$ci_percentile[2, 2]))
  cat(sprintf("  %-10s [%8.4f, %8.4f] [%8.4f, %8.4f]\n", "Profile",
              inf_result$profile$alpha$ci[1], inf_result$profile$alpha$ci[2],
              inf_result$profile$beta$ci[1], inf_result$profile$beta$ci[2]))
  cat("\n")
  
  ## Note about Hessian reliability
  hess_boot_ratio <- inf_result$hessian$se / inf_result$bootstrap$se
  if (any(hess_boot_ratio < 0.7, na.rm = TRUE)) {
    cat("NOTE: Hessian SEs appear underestimated (ratio < 0.7).\n")
    cat("      Bootstrap or Profile CIs recommended for inference.\n")
  }
  
  invisible(inf_result)
}