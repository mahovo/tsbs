## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## DCC(1,1) Estimation Diagnostics
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
##
## This file provides diagnostic tools for assessing DCC parameter estimation:
##
##   1. Monte Carlo Simulation Study
##      - Bias, RMSE, coverage probability
##      - Distribution of estimates across replications
##
##   2. Likelihood Surface Visualization
##      - Contour plots of NLL surface
##      - Parameter trace plots
##
##   3. Statistical Tests
##      - Coverage probability of confidence intervals
##      - Normality of standardized estimates
##
## Dependencies:
##   - dcc_gradient.R (NLL computation, transformations)
##   - dcc_hessian.R (standard errors)
##   - diagnostic_utils.R (simulate_dcc_garch)
##   - plotly (for interactive visualizations)
##
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


## SECTION 1: Helper Functions =================================================

#' @title Compute Standardized Residuals from Simulated DCC-GARCH Data
#' @description Given simulated returns and known GARCH parameters, compute
#'   the standardized residuals needed for DCC estimation.
#' @param y T x k matrix of simulated returns
#' @param omega Vector of GARCH omega parameters (length k)
#' @param alpha_garch Vector of GARCH alpha parameters (length k)
#' @param beta_garch Vector of GARCH beta parameters (length k)
#' @return T x k matrix of standardized residuals
#' @keywords internal
compute_std_residuals <- function(y, omega, alpha_garch, beta_garch) {
  n <- nrow(y)
  k <- ncol(y)
  
  ## Initialize
  h <- matrix(0, n, k)
  z <- matrix(0, n, k)
  
  ## Unconditional variance as starting value
  for (i in 1:k) {
    h[1, i] <- omega[i] / (1 - alpha_garch[i] - beta_garch[i])
  }
  
  ## Compute conditional variances and standardized residuals
  for (t in 1:n) {
    for (i in 1:k) {
      z[t, i] <- y[t, i] / sqrt(h[t, i])
      
      if (t < n) {
        h[t+1, i] <- omega[i] + alpha_garch[i] * y[t, i]^2 + beta_garch[i] * h[t, i]
      }
    }
  }
  
  z
}


## SECTION 2: Monte Carlo Simulation Study =====================================

#' @title Run Monte Carlo Study for DCC(1,1) Estimation
#' @description Performs a Monte Carlo simulation study to assess the accuracy
#'   of DCC parameter estimation. Computes bias, RMSE, and coverage probabilities.
#' @param n_sim Number of simulation replications
#' @param n_obs Number of observations per replication
#' @param k Number of series (default 2)
#' @param true_alpha True DCC alpha parameter
#' @param true_beta True DCC beta parameter
#' @param omega Vector of GARCH omega parameters (default: rep(0.05, k))
#' @param alpha_garch Vector of GARCH alpha parameters (default: rep(0.10, k))
#' @param beta_garch Vector of GARCH beta parameters (default: rep(0.85, k))
#' @param confidence_level Confidence level for coverage (default 0.95)
#' @param verbose Print progress
#' @param seed Base seed for reproducibility
#' @return List with:
#'   \item{estimates}{Matrix of estimates (n_sim x 2)}
#'   \item{std_errors}{Matrix of standard errors (n_sim x 2)}
#'   \item{bias}{Bias for each parameter}
#'   \item{rmse}{RMSE for each parameter}
#'   \item{coverage}{Coverage probability for each parameter}
#'   \item{summary}{Summary data frame}
#' @export
run_dcc_monte_carlo <- function(
    n_sim = 100,
    n_obs = 500,
    k = 2,
    true_alpha = 0.05,
    true_beta = 0.90,
    omega = NULL,
    alpha_garch = NULL,
    beta_garch = NULL,
    confidence_level = 0.95,
    verbose = TRUE,
    seed = 12345
) {
  
  set.seed(seed)
  
  ## Default GARCH parameters
  if (is.null(omega)) omega <- rep(0.05, k)
  if (is.null(alpha_garch)) alpha_garch <- rep(0.10, k)
  if (is.null(beta_garch)) beta_garch <- rep(0.85, k)
  
  ## Storage
  estimates <- matrix(NA, n_sim, 2)
  std_errors <- matrix(NA, n_sim, 2)
  ci_lower <- matrix(NA, n_sim, 2)
  ci_upper <- matrix(NA, n_sim, 2)
  convergence <- logical(n_sim)
  valid_se <- logical(n_sim)
  singular_info_count <- 0  ## Count singular information matrix warnings
  
  colnames(estimates) <- colnames(std_errors) <- c("alpha", "beta")
  colnames(ci_lower) <- colnames(ci_upper) <- c("alpha", "beta")
  
  true_params <- c(alpha = true_alpha, beta = true_beta)
  z_crit <- qnorm(1 - (1 - confidence_level) / 2)
  
  if (verbose) {
    cat("Running Monte Carlo study:\n")
    cat(sprintf("  Replications: %d\n", n_sim))
    cat(sprintf("  Observations: %d\n", n_obs))
    cat(sprintf("  True params: alpha=%.3f, beta=%.3f\n", true_alpha, true_beta))
    cat("\n")
  }
  
  pb <- if (verbose) txtProgressBar(min = 0, max = n_sim, style = 3) else NULL
  
  for (i in 1:n_sim) {
    ## Simulate data using simulate_dcc_garch from diagnostic_utils.R
    y <- simulate_dcc_garch(
      n = n_obs,
      k = k,
      omega = omega,
      alpha_garch = alpha_garch,
      beta_garch = beta_garch,
      alpha_dcc = true_alpha,
      beta_dcc = true_beta,
      seed = seed + i
    )
    
    ## Compute standardized residuals (using known GARCH parameters)
    std_resid <- compute_std_residuals(y, omega, alpha_garch, beta_garch)
    
    ## Estimate DCC parameters
    Qbar <- cor(std_resid)
    weights <- rep(1, n_obs)
    
    opt_result <- tryCatch({
      optim(
        par = c(true_alpha, true_beta),  ## Start at true values
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
    }, error = function(e) list(convergence = 99))
    
    convergence[i] <- opt_result$convergence == 0
    
    if (convergence[i]) {
      estimates[i, ] <- opt_result$par
      
      ## Compute standard errors (suppress warnings, collect count)
      se_result <- tryCatch(
        withCallingHandlers({
          dcc11_standard_errors(
            params = opt_result$par,
            std_resid = std_resid,
            weights = weights,
            Qbar = Qbar,
            distribution = "mvn",
            use_reparam = FALSE
          )
        }, warning = function(w) {
          if (grepl("singular", w$message, ignore.case = TRUE)) {
            singular_info_count <<- singular_info_count + 1
          }
          invokeRestart("muffleWarning")
        }),
        error = function(e) list(se = c(NA, NA))
      )
      
      std_errors[i, ] <- se_result$se
      valid_se[i] <- all(is.finite(se_result$se)) && all(se_result$se > 0)
      
      if (valid_se[i]) {
        ci_lower[i, ] <- estimates[i, ] - z_crit * std_errors[i, ]
        ci_upper[i, ] <- estimates[i, ] + z_crit * std_errors[i, ]
      }
    }
    
    if (verbose) setTxtProgressBar(pb, i)
  }
  
  if (verbose) {
    close(pb)
    cat("\n")
    if (singular_info_count > 0) {
      cat(sprintf("  Note: %d replications had singular information matrices (SEs set to NA)\n", 
                  singular_info_count))
    }
  }
  
  ## Compute summary statistics
  valid_idx <- convergence & !is.na(estimates[, 1])
  n_valid <- sum(valid_idx)
  
  ## Bias
  bias <- colMeans(estimates[valid_idx, , drop = FALSE], na.rm = TRUE) - true_params
  
  ## RMSE
  mse <- colMeans((estimates[valid_idx, , drop = FALSE] - 
                     matrix(true_params, n_valid, 2, byrow = TRUE))^2, na.rm = TRUE)
  rmse <- sqrt(mse)
  
  ## Coverage
  coverage <- rep(NA, 2)
  names(coverage) <- c("alpha", "beta")
  
  valid_ci_idx <- valid_idx & valid_se
  if (sum(valid_ci_idx) > 0) {
    for (j in 1:2) {
      covered <- ci_lower[valid_ci_idx, j] <= true_params[j] & 
        ci_upper[valid_ci_idx, j] >= true_params[j]
      coverage[j] <- mean(covered, na.rm = TRUE)
    }
  }
  
  ## Mean SE vs empirical SD
  mean_se <- colMeans(std_errors[valid_ci_idx, , drop = FALSE], na.rm = TRUE)
  empirical_sd <- apply(estimates[valid_idx, , drop = FALSE], 2, sd, na.rm = TRUE)
  
  ## Summary table
  summary_df <- data.frame(
    parameter = c("alpha", "beta"),
    true_value = true_params,
    mean_estimate = colMeans(estimates[valid_idx, , drop = FALSE], na.rm = TRUE),
    bias = bias,
    rmse = rmse,
    empirical_sd = empirical_sd,
    mean_se = mean_se,
    se_ratio = mean_se / empirical_sd,
    coverage = coverage,
    stringsAsFactors = FALSE
  )
  
  if (verbose) {
    cat("\n=== Monte Carlo Results ===\n")
    cat(sprintf("Valid replications: %d / %d (%.1f%%)\n", 
                n_valid, n_sim, 100 * n_valid / n_sim))
    cat(sprintf("Valid SEs: %d / %d (%.1f%%)\n",
                sum(valid_ci_idx), n_sim, 100 * sum(valid_ci_idx) / n_sim))
    cat("\n")
    print(summary_df, digits = 4)
  }
  
  list(
    estimates = estimates,
    std_errors = std_errors,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    convergence = convergence,
    valid_se = valid_se,
    true_params = true_params,
    bias = bias,
    rmse = rmse,
    coverage = coverage,
    mean_se = mean_se,
    empirical_sd = empirical_sd,
    summary = summary_df,
    settings = list(
      n_sim = n_sim,
      n_obs = n_obs,
      k = k,
      true_alpha = true_alpha,
      true_beta = true_beta,
      confidence_level = confidence_level
    )
  )
}


#' @title Print Monte Carlo Results
#' @export
print_mc_results <- function(mc_result) {
  cat("\n")
  cat("========================================\n")
  cat("  DCC(1,1) Monte Carlo Simulation Study\n")
  cat("========================================\n\n")
  
  s <- mc_result$settings
  cat(sprintf("Settings:\n"))
  cat(sprintf("  Replications: %d\n", s$n_sim))
  cat(sprintf("  Observations: %d\n", s$n_obs))
  cat(sprintf("  True alpha: %.4f\n", s$true_alpha))
  cat(sprintf("  True beta:  %.4f\n", s$true_beta))
  cat(sprintf("  True persistence: %.4f\n", s$true_alpha + s$true_beta))
  cat("\n")
  
  cat("Results:\n")
  cat(sprintf("  Converged: %d / %d (%.1f%%)\n",
              sum(mc_result$convergence), s$n_sim,
              100 * mean(mc_result$convergence)))
  cat(sprintf("  Valid SEs: %d / %d (%.1f%%)\n",
              sum(mc_result$valid_se), s$n_sim,
              100 * mean(mc_result$valid_se)))
  cat("\n")
  
  cat("Parameter Estimates:\n")
  cat(sprintf("  Alpha: bias=%.4f, RMSE=%.4f, coverage=%.1f%%\n",
              mc_result$bias["alpha"], mc_result$rmse["alpha"],
              100 * mc_result$coverage["alpha"]))
  cat(sprintf("  Beta:  bias=%.4f, RMSE=%.4f, coverage=%.1f%%\n",
              mc_result$bias["beta"], mc_result$rmse["beta"],
              100 * mc_result$coverage["beta"]))
  cat("\n")
  
  cat("Standard Error Calibration:\n")
  cat(sprintf("  Alpha: mean_SE=%.4f, empirical_SD=%.4f, ratio=%.2f\n",
              mc_result$mean_se["alpha"], mc_result$empirical_sd["alpha"],
              mc_result$mean_se["alpha"] / mc_result$empirical_sd["alpha"]))
  cat(sprintf("  Beta:  mean_SE=%.4f, empirical_SD=%.4f, ratio=%.2f\n",
              mc_result$mean_se["beta"], mc_result$empirical_sd["beta"],
              mc_result$mean_se["beta"] / mc_result$empirical_sd["beta"]))
  cat("\n")
  
  invisible(mc_result)
}


## SECTION 3: Likelihood Surface Visualization =================================

#' @title Compute NLL Surface for DCC(1,1)
#' @description Compute the negative log-likelihood on a grid of (alpha, beta) values.
#' @param std_resid T x k matrix of standardized residuals
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional correlation matrix
#' @param alpha_range Range for alpha (default c(0.01, 0.20))
#' @param beta_range Range for beta (default c(0.70, 0.98))
#' @param n_grid Number of grid points per dimension (default 50)
#' @param distribution "mvn" or "mvt"
#' @return List with alpha_grid, beta_grid, nll_surface, and mle location
#' @export
compute_nll_surface <- function(
    std_resid,
    weights,
    Qbar,
    alpha_range = c(0.01, 0.20),
    beta_range = c(0.70, 0.98),
    n_grid = 50,
    distribution = "mvn"
) {
  
  alpha_grid <- seq(alpha_range[1], alpha_range[2], length.out = n_grid)
  beta_grid <- seq(beta_range[1], beta_range[2], length.out = n_grid)
  
  nll_surface <- matrix(NA, n_grid, n_grid)
  
  for (i in 1:n_grid) {
    for (j in 1:n_grid) {
      ## Skip non-stationary region
      if (alpha_grid[i] + beta_grid[j] >= 1) {
        nll_surface[i, j] <- NA
        next
      }
      
      params <- c(alpha_grid[i], beta_grid[j])
      
      nll_surface[i, j] <- tryCatch({
        dcc11_nll(
          params = params,
          std_resid = std_resid,
          weights = weights,
          Qbar = Qbar,
          distribution = distribution,
          use_reparam = FALSE
        )
      }, error = function(e) NA)
    }
  }
  
  ## Find MLE on grid
  valid_idx <- which(!is.na(nll_surface), arr.ind = TRUE)
  if (nrow(valid_idx) > 0) {
    min_idx <- which.min(nll_surface)
    min_ij <- arrayInd(min_idx, dim(nll_surface))
    mle_grid <- c(alpha = alpha_grid[min_ij[1]], beta = beta_grid[min_ij[2]])
    mle_nll <- nll_surface[min_idx]
  } else {
    mle_grid <- c(alpha = NA, beta = NA)
    mle_nll <- NA
  }
  
  list(
    alpha_grid = alpha_grid,
    beta_grid = beta_grid,
    nll_surface = nll_surface,
    mle_grid = mle_grid,
    mle_nll = mle_nll
  )
}


#' @title Plot NLL Contours with plotly
#' @description Create an interactive contour plot of the DCC likelihood surface.
#' @param surface Result from compute_nll_surface()
#' @param true_params Optional vector c(alpha, beta) of true parameters
#' @param mle_params Optional vector c(alpha, beta) of MLE estimates
#' @param title Plot title
#' @return plotly object
#' @export
plot_nll_contours <- function(
    surface,
    true_params = NULL,
    mle_params = NULL,
    title = "DCC(1,1) Negative Log-Likelihood Surface"
) {
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for this function")
  }
  
  ## Create base contour plot
  p <- plotly::plot_ly(
    x = surface$alpha_grid,
    y = surface$beta_grid,
    z = t(surface$nll_surface),  ## Transpose for correct orientation
    type = "contour",
    contours = list(
      showlabels = TRUE,
      labelfont = list(size = 10, color = "white")
    ),
    colorscale = "Viridis",
    colorbar = list(title = "NLL")
  )
  
  ## Add stationarity boundary
  alpha_boundary <- seq(0, 1, length.out = 100)
  beta_boundary <- 1 - alpha_boundary
  
  p <- p %>% plotly::add_trace(
    x = alpha_boundary,
    y = beta_boundary,
    type = "scatter",
    mode = "lines",
    line = list(color = "red", width = 2, dash = "dash"),
    name = "α + β = 1",
    showlegend = TRUE
  )
  
  ## Add true parameters if provided
  if (!is.null(true_params)) {
    p <- p %>% plotly::add_trace(
      x = true_params[1],
      y = true_params[2],
      type = "scatter",
      mode = "markers",
      marker = list(color = "green", size = 12, symbol = "star"),
      name = sprintf("True (%.3f, %.3f)", true_params[1], true_params[2]),
      showlegend = TRUE
    )
  }
  
  ## Add MLE if provided
  if (!is.null(mle_params)) {
    p <- p %>% plotly::add_trace(
      x = mle_params[1],
      y = mle_params[2],
      type = "scatter",
      mode = "markers",
      marker = list(color = "red", size = 12, symbol = "x"),
      name = sprintf("MLE (%.3f, %.3f)", mle_params[1], mle_params[2]),
      showlegend = TRUE
    )
  }
  
  ## Layout
  p <- p %>% plotly::layout(
    title = title,
    xaxis = list(title = "α (alpha)"),
    yaxis = list(title = "β (beta)"),
    showlegend = TRUE
  )
  
  p
}


#' @title Plot 3D NLL Surface with plotly
#' @description Create an interactive 3D surface plot of the DCC likelihood.
#' @param surface Result from compute_nll_surface()
#' @param true_params Optional vector c(alpha, beta) of true parameters
#' @param mle_params Optional vector c(alpha, beta) of MLE estimates
#' @param title Plot title
#' @return plotly object
#' @export
plot_nll_surface_3d <- function(
    surface,
    true_params = NULL,
    mle_params = NULL,
    title = "DCC(1,1) Negative Log-Likelihood Surface"
) {
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for this function")
  }
  
  ## Create 3D surface
  p <- plotly::plot_ly(
    x = surface$alpha_grid,
    y = surface$beta_grid,
    z = t(surface$nll_surface),
    type = "surface",
    colorscale = "Viridis",
    colorbar = list(title = "NLL")
  )
  
  ## Add true parameters marker
  if (!is.null(true_params)) {
    ## Get NLL at true params
    i <- which.min(abs(surface$alpha_grid - true_params[1]))
    j <- which.min(abs(surface$beta_grid - true_params[2]))
    z_true <- surface$nll_surface[i, j]
    
    p <- p %>% plotly::add_trace(
      x = true_params[1],
      y = true_params[2],
      z = z_true,
      type = "scatter3d",
      mode = "markers",
      marker = list(color = "green", size = 8, symbol = "diamond"),
      name = "True"
    )
  }
  
  ## Add MLE marker
  if (!is.null(mle_params)) {
    i <- which.min(abs(surface$alpha_grid - mle_params[1]))
    j <- which.min(abs(surface$beta_grid - mle_params[2]))
    z_mle <- surface$nll_surface[i, j]
    
    p <- p %>% plotly::add_trace(
      x = mle_params[1],
      y = mle_params[2],
      z = z_mle,
      type = "scatter3d",
      mode = "markers",
      marker = list(color = "red", size = 8, symbol = "x"),
      name = "MLE"
    )
  }
  
  ## Layout
  p <- p %>% plotly::layout(
    title = title,
    scene = list(
      xaxis = list(title = "α"),
      yaxis = list(title = "β"),
      zaxis = list(title = "NLL")
    )
  )
  
  p
}


## SECTION 4: Monte Carlo Visualization ========================================

#' @title Plot Monte Carlo Estimate Distribution
#' @description Create histograms and scatter plots of MC estimates.
#' @param mc_result Result from run_dcc_monte_carlo()
#' @return plotly object with subplots
#' @export
plot_mc_estimates <- function(mc_result) {
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for this function")
  }
  
  valid_idx <- mc_result$convergence & !is.na(mc_result$estimates[, 1])
  alpha_est <- mc_result$estimates[valid_idx, "alpha"]
  beta_est <- mc_result$estimates[valid_idx, "beta"]
  
  true_alpha <- mc_result$true_params["alpha"]
  true_beta <- mc_result$true_params["beta"]
  
  ## Create scatter plot of estimates
  p1 <- plotly::plot_ly(
    x = alpha_est,
    y = beta_est,
    type = "scatter",
    mode = "markers",
    marker = list(color = "blue", opacity = 0.5, size = 6),
    name = "Estimates"
  ) %>%
    plotly::add_trace(
      x = true_alpha,
      y = true_beta,
      type = "scatter",
      mode = "markers",
      marker = list(color = "red", size = 15, symbol = "star"),
      name = "True"
    ) %>%
    plotly::add_trace(
      x = mean(alpha_est),
      y = mean(beta_est),
      type = "scatter",
      mode = "markers",
      marker = list(color = "green", size = 12, symbol = "diamond"),
      name = "Mean"
    ) %>%
    plotly::layout(
      title = "DCC Parameter Estimates",
      xaxis = list(title = "α"),
      yaxis = list(title = "β")
    )
  
  ## Alpha histogram
  p2 <- plotly::plot_ly(
    x = alpha_est,
    type = "histogram",
    name = "α estimates",
    marker = list(color = "steelblue")
  ) %>%
    plotly::add_trace(
      x = c(true_alpha, true_alpha),
      y = c(0, max(table(cut(alpha_est, 30)))),
      type = "scatter",
      mode = "lines",
      line = list(color = "red", width = 2, dash = "dash"),
      name = "True α"
    ) %>%
    plotly::layout(
      title = "Alpha Distribution",
      xaxis = list(title = "α"),
      yaxis = list(title = "Count")
    )
  
  ## Beta histogram
  p3 <- plotly::plot_ly(
    x = beta_est,
    type = "histogram",
    name = "β estimates",
    marker = list(color = "steelblue")
  ) %>%
    plotly::add_trace(
      x = c(true_beta, true_beta),
      y = c(0, max(table(cut(beta_est, 30)))),
      type = "scatter",
      mode = "lines",
      line = list(color = "red", width = 2, dash = "dash"),
      name = "True β"
    ) %>%
    plotly::layout(
      title = "Beta Distribution",
      xaxis = list(title = "β"),
      yaxis = list(title = "Count")
    )
  
  list(scatter = p1, alpha_hist = p2, beta_hist = p3)
}


#' @title Plot Coverage Diagnostic
#' @description Visualize confidence interval coverage for each replication.
#' @param mc_result Result from run_dcc_monte_carlo()
#' @param param Which parameter ("alpha" or "beta")
#' @param max_show Maximum number of replications to show
#' @return plotly object
#' @export
plot_coverage_diagnostic <- function(mc_result, param = "alpha", max_show = 50) {
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for this function")
  }
  
  valid_idx <- which(mc_result$valid_se)
  if (length(valid_idx) > max_show) {
    valid_idx <- valid_idx[1:max_show]
  }
  
  j <- if (param == "alpha") 1 else 2
  true_val <- mc_result$true_params[j]
  
  est <- mc_result$estimates[valid_idx, j]
  lower <- mc_result$ci_lower[valid_idx, j]
  upper <- mc_result$ci_upper[valid_idx, j]
  
  ## Determine if CI covers true value
  covers <- lower <= true_val & upper >= true_val
  colors <- ifelse(covers, "blue", "red")
  
  ## Create plot
  n <- length(valid_idx)
  
  p <- plotly::plot_ly() %>%
    ## Add CIs as segments
    plotly::add_trace(
      x = rep(1:n, each = 2),
      y = as.vector(rbind(lower, upper)),
      type = "scatter",
      mode = "lines",
      line = list(color = "gray", width = 1),
      showlegend = FALSE
    )
  
  ## Add points
  for (i in 1:n) {
    p <- p %>% plotly::add_trace(
      x = i,
      y = est[i],
      type = "scatter",
      mode = "markers",
      marker = list(color = colors[i], size = 6),
      showlegend = FALSE
    )
  }
  
  ## Add true value line
  p <- p %>% plotly::add_trace(
    x = c(0, n + 1),
    y = c(true_val, true_val),
    type = "scatter",
    mode = "lines",
    line = list(color = "green", width = 2),
    name = sprintf("True %s = %.3f", param, true_val)
  )
  
  ## Layout
  p <- p %>% plotly::layout(
    title = sprintf("%s: CI Coverage (%.1f%% empirical)", 
                    toupper(param), 100 * mean(covers)),
    xaxis = list(title = "Replication"),
    yaxis = list(title = param),
    showlegend = TRUE
  )
  
  p
}


## SECTION 5: Parameter Trace Visualization (for EM/optimization) ==============

#' @title Extract DCC Parameter Trajectory from MS Diagnostics
#' @description Extract the evolution of DCC parameters (alpha, beta) across EM 
#'   iterations from an ms_diagnostics object.
#' @param diagnostics An object of class \code{ms_diagnostics} from fit_ms_varma_garch()
#' @param state Integer state index (default 1)
#' @return Data frame with columns: iteration, alpha, beta
#'   Returns NULL if DCC parameters not found.
#' @export
extract_dcc_trajectory <- function(diagnostics, state = 1) {
  
  if (!inherits(diagnostics, "ms_diagnostics")) {
    stop("diagnostics must be an object of class 'ms_diagnostics'")
  }
  
  ## Extract alpha and beta trajectories
  alpha_traj <- extract_param_trajectory(diagnostics, state = state, param_name = "alpha_1")
  beta_traj <- extract_param_trajectory(diagnostics, state = state, param_name = "beta_1")
  
  if (is.null(alpha_traj) || is.null(beta_traj)) {
    warning("Could not extract DCC parameters for state ", state)
    return(NULL)
  }
  
  ## Combine into single data frame
  data.frame(
    iteration = alpha_traj$iteration,
    alpha = alpha_traj$value,
    beta = beta_traj$value
  )
}


#' @title Extract DCC Trajectory with Log-Likelihood
#' @description Extract DCC parameter evolution along with log-likelihood values.
#' @param diagnostics An object of class \code{ms_diagnostics}
#' @param state Integer state index (default 1)
#' @return Data frame with columns: iteration, alpha, beta, log_lik
#' @export
extract_dcc_trajectory_with_ll <- function(diagnostics, state = 1) {
  
  ## Get parameter trajectory
  traj <- extract_dcc_trajectory(diagnostics, state = state)
  
  if (is.null(traj)) return(NULL)
  
  ## Get log-likelihood trajectory
  ll_after <- extract_ll_trajectory(diagnostics, type = "after")
  
  ## Match iterations (LL may have different length due to initialization)
  n_iter <- nrow(traj)
  if (length(ll_after) >= n_iter) {
    traj$log_lik <- ll_after[1:n_iter]
  } else {
    ## Pad with NA if needed
    traj$log_lik <- c(ll_after, rep(NA, n_iter - length(ll_after)))
  }
  
  traj
}


#' @title Plot Parameter Trace from Optimization
#' @description Visualize the path of parameter estimates during optimization
#'   overlaid on the likelihood contours.
#' @param trace_data Data frame with columns: iteration, alpha, beta, nll
#' @param surface Result from compute_nll_surface() (optional)
#' @param true_params True parameter values (optional)
#' @param title Plot title
#' @return plotly object
#' @export
plot_optimization_trace <- function(
    trace_data,
    surface = NULL,
    true_params = NULL,
    title = "Optimization Trace"
) {
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for this function")
  }
  
  ## Start with contours if surface provided
  if (!is.null(surface)) {
    p <- plotly::plot_ly(
      x = surface$alpha_grid,
      y = surface$beta_grid,
      z = t(surface$nll_surface),
      type = "contour",
      contours = list(showlabels = FALSE),
      colorscale = "Viridis",
      showscale = FALSE
    )
  } else {
    p <- plotly::plot_ly()
  }
  
  ## Add trace path
  p <- p %>% plotly::add_trace(
    x = trace_data$alpha,
    y = trace_data$beta,
    type = "scatter",
    mode = "lines+markers",
    line = list(color = "orange", width = 2),
    marker = list(
      color = trace_data$iteration,
      colorscale = "Viridis",
      size = 8,
      showscale = TRUE,
      colorbar = list(title = "Iteration", len = 0.5)
    ),
    name = "Trace"
  )
  
  ## Add start point
  p <- p %>% plotly::add_trace(
    x = trace_data$alpha[1],
    y = trace_data$beta[1],
    type = "scatter",
    mode = "markers",
    marker = list(color = "green", size = 15, symbol = "circle"),
    name = "Start"
  )
  
  ## Add end point
  n <- nrow(trace_data)
  p <- p %>% plotly::add_trace(
    x = trace_data$alpha[n],
    y = trace_data$beta[n],
    type = "scatter",
    mode = "markers",
    marker = list(color = "red", size = 15, symbol = "x"),
    name = "End"
  )
  
  ## Add true params if provided
  if (!is.null(true_params)) {
    p <- p %>% plotly::add_trace(
      x = true_params[1],
      y = true_params[2],
      type = "scatter",
      mode = "markers",
      marker = list(color = "gold", size = 15, symbol = "star"),
      name = "True"
    )
  }
  
  ## Layout
  p <- p %>% plotly::layout(
    title = title,
    xaxis = list(title = "α"),
    yaxis = list(title = "β")
  )
  
  p
}


#' @title Visualize DCC Parameter Evolution from EM Fit
#' @description Plot the evolution of DCC parameters across EM iterations,
#'   optionally overlaid on the likelihood surface.
#' @param fit A fitted model from fit_ms_varma_garch() with diagnostics=TRUE
#' @param state Integer state index (default 1)
#' @param show_surface Logical: compute and show likelihood contours? (default TRUE)
#' @param true_params Optional true parameter values c(alpha, beta) for simulation studies
#' @param n_grid Grid size for surface computation (default 50)
#' @return List with trace data, surface (if computed), and plotly visualization
#' @export
visualize_dcc_evolution <- function(
    fit,
    state = 1,
    show_surface = TRUE,
    true_params = NULL,
    n_grid = 50
) {
  
  ## Check for diagnostics
  
  if (is.null(fit$diagnostics)) {
    stop("Fit object does not contain diagnostics. Re-run with diagnostics=TRUE")
  }
  
  ## Extract DCC trajectory
  traj <- extract_dcc_trajectory_with_ll(fit$diagnostics, state = state)
  
  if (is.null(traj)) {
    stop("Could not extract DCC trajectory for state ", state)
  }
  
  ## Get final estimates
  final_alpha <- traj$alpha[nrow(traj)]
  final_beta <- traj$beta[nrow(traj)]
  
  ## Compute likelihood surface if requested
  surface <- NULL
  if (show_surface) {
    ## Need standardized residuals - extract from fit or recompute
    ## This assumes fit contains the necessary information
    if (!is.null(fit$std_resid)) {
      std_resid <- fit$std_resid
    } else if (!is.null(fit$data) && !is.null(fit$garch_fit)) {
      ## Would need to extract from GARCH fit
      warning("Cannot extract standardized residuals from fit. Skipping surface.")
      show_surface <- FALSE
    } else {
      warning("Cannot compute likelihood surface without standardized residuals")
      show_surface <- FALSE
    }
    
    if (show_surface) {
      Qbar <- cor(std_resid)
      weights <- rep(1, nrow(std_resid))
      
      ## Determine grid range based on trajectory
      alpha_range <- range(traj$alpha, na.rm = TRUE)
      beta_range <- range(traj$beta, na.rm = TRUE)
      
      ## Expand slightly
      alpha_range <- c(max(0.001, alpha_range[1] - 0.02), 
                       min(0.30, alpha_range[2] + 0.05))
      beta_range <- c(max(0.50, beta_range[1] - 0.05), 
                      min(0.99, beta_range[2] + 0.02))
      
      surface <- compute_nll_surface(
        std_resid = std_resid,
        weights = weights,
        Qbar = Qbar,
        alpha_range = alpha_range,
        beta_range = beta_range,
        n_grid = n_grid
      )
    }
  }
  
  ## Create plot
  p <- plot_optimization_trace(
    trace_data = traj,
    surface = surface,
    true_params = true_params,
    title = sprintf("DCC Parameter Evolution (State %d): %d iterations → (%.3f, %.3f)",
                    state, nrow(traj), final_alpha, final_beta)
  )
  
  list(
    trace = traj,
    surface = surface,
    plot = p,
    final_params = c(alpha = final_alpha, beta = final_beta)
  )
}


#' @title Quick Visualization of DCC Optimization Path (Standalone)
#' @description Convenience function for visualizing optimization on simulated data
#'   without needing a full MS-VARMA-GARCH fit. Useful for testing/diagnostics.
#' @param n Number of observations
#' @param true_alpha True DCC alpha
#' @param true_beta True DCC beta
#' @param start_alpha Starting value for alpha
#' @param start_beta Starting value for beta
#' @param omega Vector of GARCH omega parameters
#' @param alpha_garch Vector of GARCH alpha parameters
#' @param beta_garch Vector of GARCH beta parameters
#' @param n_grid Grid size for surface (default 50)
#' @param seed Random seed
#' @return List with trace data, surface, and plot
#' @export
visualize_standalone_optimization <- function(
    n = 500,
    true_alpha = 0.05,
    true_beta = 0.90,
    start_alpha = 0.10,
    start_beta = 0.80,
    omega = c(0.05, 0.08),
    alpha_garch = c(0.10, 0.12),
    beta_garch = c(0.85, 0.82),
    n_grid = 50,
    seed = 42
) {
  
  k <- length(omega)
  
  ## Simulate data
  y <- simulate_dcc_garch(
    n = n,
    k = k,
    omega = omega,
    alpha_garch = alpha_garch,
    beta_garch = beta_garch,
    alpha_dcc = true_alpha,
    beta_dcc = true_beta,
    seed = seed
  )
  
  ## Compute standardized residuals
  std_resid <- compute_std_residuals(y, omega, alpha_garch, beta_garch)
  Qbar <- cor(std_resid)
  weights <- rep(1, n)
  
  ## Compute likelihood surface
  surface <- compute_nll_surface(
    std_resid = std_resid,
    weights = weights,
    Qbar = Qbar,
    # alpha_range = c(max(0.001, min(true_alpha, start_alpha) - 0.05), 
    #                 min(0.30, max(true_alpha, start_alpha) + 0.10)),
    # beta_range = c(max(0.50, min(true_beta, start_beta) - 0.15), 
    #                min(0.99, max(true_beta, start_beta) + 0.05)),
    alpha_range = c(0, 1),
    beta_range = c(0, 1),
    n_grid = n_grid
  )
  
  ## Run optimization and record path manually
  trace_list <- list()
  iter_count <- 0
  
  nll_with_trace <- function(params) {
    nll <- dcc11_nll(
      params = params,
      std_resid = std_resid,
      weights = weights,
      Qbar = Qbar,
      distribution = "mvn",
      use_reparam = FALSE
    )
    iter_count <<- iter_count + 1
    trace_list[[iter_count]] <<- data.frame(
      iteration = iter_count,
      alpha = params[1],
      beta = params[2],
      nll = nll
    )
    nll
  }
  
  opt_result <- optim(
    par = c(start_alpha, start_beta),
    fn = nll_with_trace,
    method = "L-BFGS-B",
    lower = c(1e-6, 1e-6),
    upper = c(0.5, 0.999)
  )
  
  ## Combine trace
  trace_df <- do.call(rbind, trace_list)
  
  ## Remove duplicates
  if (nrow(trace_df) > 1) {
    keep <- c(TRUE, diff(trace_df$alpha) != 0 | diff(trace_df$beta) != 0)
    trace_df <- trace_df[keep, ]
    trace_df$iteration <- seq_len(nrow(trace_df))
  }
  
  ## Create plot
  p <- plot_optimization_trace(
    trace_data = trace_df,
    surface = surface,
    true_params = c(true_alpha, true_beta),
    title = sprintf("Optimization Path: Start (%.2f, %.2f) → MLE (%.3f, %.3f)",
                    start_alpha, start_beta, 
                    opt_result$par[1], opt_result$par[2])
  )
  
  list(
    trace = trace_df,
    mle = opt_result$par,
    surface = surface,
    plot = p,
    true_params = c(alpha = true_alpha, beta = true_beta),
    start_params = c(alpha = start_alpha, beta = start_beta)
  )
}


## SECTION 6: Statistical Tests ================================================

#' @title Test Normality of Standardized Estimates
#' @description Test whether (theta_hat - theta_0) / SE follows N(0,1).
#' @param mc_result Result from run_dcc_monte_carlo()
#' @return List with test statistics for each parameter
#' @export
test_estimate_normality <- function(mc_result) {
  
  valid_idx <- mc_result$convergence & mc_result$valid_se
  
  results <- list()
  
  for (param in c("alpha", "beta")) {
    j <- if (param == "alpha") 1 else 2
    
    est <- mc_result$estimates[valid_idx, j]
    se <- mc_result$std_errors[valid_idx, j]
    true_val <- mc_result$true_params[j]
    
    ## Standardized estimates
    z_scores <- (est - true_val) / se
    
    ## Shapiro-Wilk test (if n <= 5000)
    n <- length(z_scores)
    if (n >= 3 && n <= 5000) {
      sw_test <- shapiro.test(z_scores)
    } else {
      sw_test <- list(statistic = NA, p.value = NA)
    }
    
    ## Basic summary
    results[[param]] <- list(
      n = n,
      mean_z = mean(z_scores),
      sd_z = sd(z_scores),
      skewness = mean((z_scores - mean(z_scores))^3) / sd(z_scores)^3,
      kurtosis = mean((z_scores - mean(z_scores))^4) / sd(z_scores)^4 - 3,
      shapiro_W = as.numeric(sw_test$statistic),
      shapiro_p = sw_test$p.value,
      z_scores = z_scores
    )
  }
  
  ## Print summary
  cat("\n=== Normality Test for Standardized Estimates ===\n\n")
  cat("Under correct specification: (θ̂ - θ₀) / SE ~ N(0, 1)\n\n")
  
  for (param in c("alpha", "beta")) {
    r <- results[[param]]
    cat(sprintf("%s:\n", toupper(param)))
    cat(sprintf("  N valid: %d\n", r$n))
    cat(sprintf("  Mean z-score: %.3f (should be ≈ 0)\n", r$mean_z))
    cat(sprintf("  SD z-score: %.3f (should be ≈ 1)\n", r$sd_z))
    cat(sprintf("  Skewness: %.3f (should be ≈ 0)\n", r$skewness))
    cat(sprintf("  Excess kurtosis: %.3f (should be ≈ 0)\n", r$kurtosis))
    if (!is.na(r$shapiro_p)) {
      cat(sprintf("  Shapiro-Wilk p-value: %.4f\n", r$shapiro_p))
    }
    cat("\n")
  }
  
  invisible(results)
}


#' @title Summary Function for Complete Diagnostic Analysis
#' @description Run a complete diagnostic analysis including MC study,
#'   likelihood surface, and statistical tests.
#' @param n_sim Number of MC replications
#' @param n_obs Observations per replication
#' @param true_alpha True DCC alpha
#' @param true_beta True DCC beta
#' @param omega Vector of GARCH omega parameters
#' @param alpha_garch Vector of GARCH alpha parameters
#' @param beta_garch Vector of GARCH beta parameters
#' @param seed Random seed
#' @param plot Logical: create plots?
#' @return List with all diagnostic results
#' @export
run_complete_diagnostics <- function(
    n_sim = 50,
    n_obs = 500,
    true_alpha = 0.05,
    true_beta = 0.90,
    omega = c(0.05, 0.08),
    alpha_garch = c(0.10, 0.12),
    beta_garch = c(0.85, 0.82),
    seed = 42,
    plot = TRUE
) {
  
  k <- length(omega)
  
  cat("\n")
  cat("================================================================\n")
  cat("  Complete DCC(1,1) Estimation Diagnostics\n")
  cat("================================================================\n\n")
  
  ## 1. Run Monte Carlo study
  cat("Step 1: Running Monte Carlo simulation...\n")
  mc_result <- run_dcc_monte_carlo(
    n_sim = n_sim,
    n_obs = n_obs,
    k = k,
    true_alpha = true_alpha,
    true_beta = true_beta,
    omega = omega,
    alpha_garch = alpha_garch,
    beta_garch = beta_garch,
    seed = seed,
    verbose = TRUE
  )
  
  ## 2. Generate one dataset for surface visualization
  cat("\nStep 2: Computing likelihood surface...\n")
  y <- simulate_dcc_garch(
    n = n_obs,
    k = k,
    omega = omega,
    alpha_garch = alpha_garch,
    beta_garch = beta_garch,
    alpha_dcc = true_alpha,
    beta_dcc = true_beta,
    seed = seed
  )
  
  std_resid <- compute_std_residuals(y, omega, alpha_garch, beta_garch)
  
  surface <- compute_nll_surface(
    std_resid = std_resid,
    weights = rep(1, n_obs),
    Qbar = cor(std_resid),
    alpha_range = c(max(0.001, true_alpha - 0.10), min(0.30, true_alpha + 0.15)),
    beta_range = c(max(0.50, true_beta - 0.20), min(0.99, true_beta + 0.05)),
    n_grid = 40
  )
  
  ## 3. Statistical tests
  cat("\nStep 3: Running statistical tests...\n")
  normality_tests <- test_estimate_normality(mc_result)
  
  ## 4. Create plots
  plots <- NULL
  if (plot && requireNamespace("plotly", quietly = TRUE)) {
    cat("\nStep 4: Creating visualizations...\n")
    
    plots <- list(
      contour = plot_nll_contours(
        surface,
        true_params = c(true_alpha, true_beta),
        title = "DCC(1,1) Likelihood Surface"
      ),
      mc_scatter = plot_mc_estimates(mc_result)$scatter,
      alpha_hist = plot_mc_estimates(mc_result)$alpha_hist,
      beta_hist = plot_mc_estimates(mc_result)$beta_hist,
      coverage_alpha = plot_coverage_diagnostic(mc_result, "alpha"),
      coverage_beta = plot_coverage_diagnostic(mc_result, "beta")
    )
    
    cat("  Plots created. Access via result$plots\n")
  }
  
  cat("\n================================================================\n")
  cat("  Diagnostics Complete\n")
  cat("================================================================\n")
  
  list(
    mc_result = mc_result,
    surface = surface,
    normality_tests = normality_tests,
    plots = plots,
    settings = list(
      n_sim = n_sim,
      n_obs = n_obs,
      true_alpha = true_alpha,
      true_beta = true_beta,
      seed = seed
    )
  )
}