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


#### ______________________________________________________________________ ####
## SECTION 1: Helper Functions =================================================

#' @title Compute Standardized Residuals from Simulated DCC-GARCH Data
#' @description Given simulated returns and known GARCH parameters, compute
#'   the standardized residuals needed for DCC estimation.
#' @param y T x k matrix of simulated returns
#' @param omega Vector of GARCH omega parameters (length k)
#' @param alpha_garch Vector of GARCH alpha parameters (length k)
#' @param beta_garch Vector of GARCH beta parameters (length k)
#' @return T x k matrix of standardized residuals
#' 
#' @export
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


#### ______________________________________________________________________ ####
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


#### ______________________________________________________________________ ####
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


#### ______________________________________________________________________ ####
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


#### ______________________________________________________________________ ####
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

#### ______________________________________________________________________ ####
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


#### ______________________________________________________________________ ####
#### SECTION 1: CGARCH Helper Functions                                     ####

#' @title Compute Copula Residuals from Simulated CGARCH Data
#' @description Given simulated returns and known GARCH parameters, compute
#'   the standardized residuals, apply PIT transformation, and convert to
#'   copula residuals needed for CGARCH estimation.
#' @param y T x k matrix of simulated returns
#' @param omega Vector of GARCH omega parameters (length k)
#' @param alpha_garch Vector of GARCH alpha parameters (length k)
#' @param beta_garch Vector of GARCH beta parameters (length k)
#' @param transformation PIT transformation type ("parametric", "empirical", "spd")
#' @param copula Copula distribution ("mvn" or "mvt")
#' @param shape Degrees of freedom for MVT copula (default 8)
#' @return List with:
#'   \item{std_resid}{T x k matrix of standardized residuals}
#'   \item{u_matrix}{T x k matrix of PIT-transformed uniforms}
#'   \item{z_matrix}{T x k matrix of copula residuals}
#' @export
compute_cgarch_residuals <- function(y, omega, alpha_garch, beta_garch,
                                     transformation = "parametric",
                                     copula = "mvn",
                                     shape = 8) {
  n <- nrow(y)
  k <- ncol(y)
  
  ## Initialize
  h <- matrix(0, n, k)
  std_resid <- matrix(0, n, k)
  
  ## Compute conditional variances and standardized residuals
  for (i in 1:k) {
    ## Unconditional variance as starting value
    h[1, i] <- omega[i] / (1 - alpha_garch[i] - beta_garch[i])
    
    for (t in 1:n) {
      std_resid[t, i] <- y[t, i] / sqrt(h[t, i])
      
      if (t < n) {
        h[t+1, i] <- omega[i] + alpha_garch[i] * y[t, i]^2 + beta_garch[i] * h[t, i]
      }
    }
  }
  
  ## PIT transform to uniforms
  u_matrix <- matrix(0, n, k)
  
  if (transformation == "empirical") {
    for (i in 1:k) {
      u_matrix[, i] <- rank(std_resid[, i]) / (n + 1)
    }
  } else {
    ## Parametric: assume normal marginals (most common)
    for (i in 1:k) {
      u_matrix[, i] <- pnorm(std_resid[, i])
    }
  }
  
  ## Bound away from extremes
  u_matrix[u_matrix < 1e-10] <- 1e-10
  u_matrix[u_matrix > 1 - 1e-10] <- 1 - 1e-10
  
  ## Convert to copula residuals
  z_matrix <- matrix(0, n, k)
  
  if (copula == "mvn") {
    z_matrix <- qnorm(u_matrix)
  } else {
    ## MVT: use t-quantiles
    z_matrix <- qt(u_matrix, df = shape) * sqrt(shape / (shape - 2))
  }
  
  list(
    std_resid = std_resid,
    u_matrix = u_matrix,
    z_matrix = z_matrix
  )
}


#### ______________________________________________________________________ ####
#### SECTION 2: CGARCH Monte Carlo Simulation Study                         ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


#' @title Run Monte Carlo Study for CGARCH(1,1) Estimation with ADCC Support
#' @description Performs a Monte Carlo simulation study to assess the accuracy
#'   of Copula GARCH parameter estimation. Computes bias, RMSE, and coverage
#'   probabilities for DCC/ADCC correlation dynamics parameters.
#' @param n_sim Number of simulation replications
#' @param n_obs Number of observations per replication
#' @param k Number of series (default 2)
#' @param true_alpha True DCC alpha parameter
#' @param true_beta True DCC beta parameter
#' @param true_gamma True ADCC gamma parameter (default NULL for standard DCC)
#' @param omega Vector of GARCH omega parameters (default: rep(0.05, k))
#' @param alpha_garch Vector of GARCH alpha parameters (default: rep(0.10, k))
#' @param beta_garch Vector of GARCH beta parameters (default: rep(0.85, k))
#' @param copula Copula distribution ("mvn" or "mvt")
#' @param true_shape True shape parameter for MVT copula (default 8)
#' @param shape Alias for true_shape (for backward compatibility)
#' @param transformation PIT transformation type ("parametric", "empirical", "spd")
#' @param confidence_level Confidence level for coverage (default 0.95)
#' @param verbose Print progress
#' @param seed Base seed for reproducibility
#' @return List with:
#'   \item{estimates}{Matrix of estimates (n_sim x n_params)}
#'   \item{std_errors}{Matrix of standard errors (n_sim x n_params)}
#'   \item{bias}{Bias for each parameter}
#'   \item{rmse}{RMSE for each parameter}
#'   \item{coverage}{Coverage probability for each parameter}
#'   \item{persistence}{Vector of persistence values per replication}
#'   \item{summary}{Summary data frame}
#' @export
run_cgarch_monte_carlo <- function(
    n_sim = 100,
    n_obs = 500,
    k = 2,
    true_alpha = 0.05,
    true_beta = 0.90,
    true_gamma = NULL,
    omega = NULL,
    alpha_garch = NULL,
    beta_garch = NULL,
    copula = "mvn",
    true_shape = 8,
    shape = NULL,
    transformation = "parametric",
    confidence_level = 0.95,
    verbose = TRUE,
    seed = 12345
) {
  
  set.seed(seed)
  
  ## Handle shape/true_shape alias
  if (!is.null(shape) && is.null(true_shape)) {
    true_shape <- shape
  } else if (!is.null(shape)) {
    true_shape <- shape  ## shape takes precedence if both provided
  }
  
  ## Determine if using ADCC (has gamma parameter)
  use_adcc <- !is.null(true_gamma) && true_gamma != 0
  
  ## Default GARCH parameters
  if (is.null(omega)) omega <- rep(0.05, k)
  if (is.null(alpha_garch)) alpha_garch <- rep(0.10, k)
  if (is.null(beta_garch)) beta_garch <- rep(0.85, k)
  
  ## Number of parameters depends on copula type and ADCC
  ## Base: alpha, beta
  ## + gamma if ADCC
  ## + shape if MVT copula
  n_base_params <- 2
  if (use_adcc) n_base_params <- 3
  n_params <- n_base_params + if (copula == "mvt") 1 else 0
  
  ## Parameter names
  if (use_adcc) {
    if (copula == "mvt") {
      param_names <- c("alpha", "beta", "gamma", "shape")
    } else {
      param_names <- c("alpha", "beta", "gamma")
    }
  } else {
    if (copula == "mvt") {
      param_names <- c("alpha", "beta", "shape")
    } else {
      param_names <- c("alpha", "beta")
    }
  }
  
  ## Storage
  estimates <- matrix(NA, n_sim, n_params)
  std_errors <- matrix(NA, n_sim, n_params)
  ci_lower <- matrix(NA, n_sim, n_params)
  ci_upper <- matrix(NA, n_sim, n_params)
  convergence <- logical(n_sim)
  valid_se <- logical(n_sim)
  persistence <- numeric(n_sim)
  singular_info_count <- 0
  
  colnames(estimates) <- colnames(std_errors) <- param_names
  colnames(ci_lower) <- colnames(ci_upper) <- param_names
  
  ## True parameters vector
  if (use_adcc) {
    if (copula == "mvt") {
      true_params <- c(alpha = true_alpha, beta = true_beta, 
                       gamma = true_gamma, shape = true_shape)
    } else {
      true_params <- c(alpha = true_alpha, beta = true_beta, gamma = true_gamma)
    }
    true_persistence <- true_alpha + true_beta + 0.5 * true_gamma
  } else {
    if (copula == "mvt") {
      true_params <- c(alpha = true_alpha, beta = true_beta, shape = true_shape)
    } else {
      true_params <- c(alpha = true_alpha, beta = true_beta)
    }
    true_persistence <- true_alpha + true_beta
  }
  
  z_crit <- qnorm(1 - (1 - confidence_level) / 2)
  
  if (verbose) {
    cat("Running CGARCH Monte Carlo study:\n")
    cat(sprintf("  Replications: %d\n", n_sim))
    cat(sprintf("  Observations: %d\n", n_obs))
    cat(sprintf("  Copula: %s\n", copula))
    cat(sprintf("  Dynamics: %s\n", if (use_adcc) "ADCC" else "DCC"))
    cat(sprintf("  True alpha: %.3f\n", true_alpha))
    cat(sprintf("  True beta: %.3f\n", true_beta))
    if (use_adcc) {
      cat(sprintf("  True gamma: %.3f\n", true_gamma))
    }
    if (copula == "mvt") {
      cat(sprintf("  True shape: %.1f\n", true_shape))
    }
    cat(sprintf("  True persistence: %.3f\n", true_persistence))
    cat("\n")
  }
  
  pb <- if (verbose) txtProgressBar(min = 0, max = n_sim, style = 3) else NULL
  
  for (i in 1:n_sim) {
    if (verbose) setTxtProgressBar(pb, i)
    
    tryCatch({
      ## Simulate data with ADCC support
      y <- simulate_cgarch(
        n = n_obs,
        k = k,
        omega = omega,
        alpha_garch = alpha_garch,
        beta_garch = beta_garch,
        alpha_dcc = true_alpha,
        beta_dcc = true_beta,
        gamma_dcc = if (use_adcc) true_gamma else NULL,
        copula = copula,
        shape = true_shape,
        seed = seed + i
      )
      
      ## Compute standardized residuals and PIT transform
      resid_result <- compute_cgarch_residuals(
        y, omega, alpha_garch, beta_garch,
        transformation = transformation,
        copula = copula,
        shape = true_shape
      )
      
      z_matrix <- resid_result$z_matrix
      Qbar <- cor(z_matrix)
      weights <- rep(1, n_obs)
      
      ## Fit model
      if (use_adcc) {
        ## Fit ADCC model
        fit_result <- fit_adcc_copula(
          z_matrix = z_matrix,
          weights = weights,
          Qbar = Qbar,
          copula_dist = copula,
          shape = if (copula == "mvt") NULL else true_shape  ## Estimate shape for MVT
        )
      } else {
        ## Fit standard DCC model
        fit_result <- fit_cgarch_copula(
          z_matrix = z_matrix,
          weights = weights,
          Qbar = Qbar,
          copula_dist = copula
        )
      }
      
      if (!is.null(fit_result) && fit_result$convergence) {
        convergence[i] <- TRUE
        
        ## Extract estimates
        if (use_adcc) {
          estimates[i, "alpha"] <- fit_result$alpha
          estimates[i, "beta"] <- fit_result$beta
          estimates[i, "gamma"] <- fit_result$gamma
          persistence[i] <- fit_result$alpha + fit_result$beta + 0.5 * fit_result$gamma
          if (copula == "mvt" && !is.null(fit_result$shape)) {
            estimates[i, "shape"] <- fit_result$shape
          }
        } else {
          estimates[i, "alpha"] <- fit_result$alpha
          estimates[i, "beta"] <- fit_result$beta
          persistence[i] <- fit_result$alpha + fit_result$beta
          if (copula == "mvt" && !is.null(fit_result$shape)) {
            estimates[i, "shape"] <- fit_result$shape
          }
        }
        
        ## Compute standard errors using cgarch_standard_errors from hessian_se.R
        se_result <- tryCatch({
          if (use_adcc) {
            ## For ADCC, construct result and use compute_cgarch_standard_errors_robust
            cgarch_result <- list(
              dcc_pars = list(alpha_1 = fit_result$alpha, 
                              gamma_1 = fit_result$gamma,
                              beta_1 = fit_result$beta),
              dist_pars = if (copula == "mvt") list(shape = fit_result$shape) else list(),
              correlation_type = "dynamic"
            )
            ## Note: ADCC SE computation may need special handling
            ## For now, use standard CGARCH SE as approximation (ignoring gamma)
            params <- if (copula == "mvt") {
              c(fit_result$alpha, fit_result$beta, fit_result$shape)
            } else {
              c(fit_result$alpha, fit_result$beta)
            }
            cgarch_standard_errors(
              params = params,
              z_matrix = z_matrix,
              weights = weights,
              Qbar = Qbar,
              copula_dist = copula,
              use_reparam = FALSE,
              method = "hessian"
            )
          } else {
            ## For standard DCC dynamics, use cgarch_standard_errors directly
            params <- if (copula == "mvt") {
              c(fit_result$alpha, fit_result$beta, fit_result$shape)
            } else {
              c(fit_result$alpha, fit_result$beta)
            }
            cgarch_standard_errors(
              params = params,
              z_matrix = z_matrix,
              weights = weights,
              Qbar = Qbar,
              copula_dist = copula,
              use_reparam = FALSE,
              method = "hessian"
            )
          }
        }, error = function(e) NULL)
        
        ## Extract SEs from result
        if (!is.null(se_result) && !is.null(se_result$se) && 
            all(is.finite(se_result$se))) {
          valid_se[i] <- TRUE
          
          std_errors[i, "alpha"] <- se_result$se["alpha"]
          std_errors[i, "beta"] <- se_result$se["beta"]
          
          if (use_adcc) {
            ## For ADCC, we don't have SE for gamma from standard method
            ## Set to NA or compute separately
            std_errors[i, "gamma"] <- NA
          }
          
          if (copula == "mvt" && "shape" %in% names(se_result$se)) {
            std_errors[i, "shape"] <- se_result$se["shape"]
          }
          
          ## Compute confidence intervals for parameters with valid SEs
          for (j in 1:n_params) {
            if (!is.na(std_errors[i, j])) {
              ci_lower[i, j] <- estimates[i, j] - z_crit * std_errors[i, j]
              ci_upper[i, j] <- estimates[i, j] + z_crit * std_errors[i, j]
            }
          }
        }
      }
      
    }, error = function(e) {
      ## Silently continue on error
    })
  }
  
  if (verbose) {
    close(pb)
    cat("\n")
  }
  
  ## Compute summary statistics
  valid_idx <- convergence
  valid_ci_idx <- convergence & valid_se
  n_valid <- sum(valid_idx)
  
  ## Bias and RMSE
  bias <- colMeans(estimates[valid_idx, , drop = FALSE], na.rm = TRUE) - true_params
  rmse <- sqrt(colMeans((estimates[valid_idx, , drop = FALSE] - 
                           matrix(true_params, nrow = n_valid, ncol = n_params, 
                                  byrow = TRUE))^2, na.rm = TRUE))
  
  ## Coverage
  coverage <- rep(NA, n_params)
  names(coverage) <- param_names
  
  if (sum(valid_ci_idx) > 0) {
    for (j in 1:n_params) {
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
    parameter = param_names,
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
    cat("\n=== CGARCH Monte Carlo Results ===\n")
    cat(sprintf("Model: %s with %s copula\n", 
                if (use_adcc) "ADCC" else "DCC", toupper(copula)))
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
    persistence = persistence,
    true_params = true_params,
    true_persistence = true_persistence,
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
      true_gamma = true_gamma,
      true_shape = if (copula == "mvt") true_shape else NULL,
      true_persistence = true_persistence,
      copula = copula,
      transformation = transformation,
      confidence_level = confidence_level,
      use_adcc = use_adcc
    )
  )
}


#' @title Print CGARCH Monte Carlo Results
#' @export
print_cgarch_mc_results <- function(mc_result) {
  cat("\n")
  cat("============================================\n")
  cat("  CGARCH(1,1) Monte Carlo Simulation Study\n")
  cat("============================================\n\n")
  
  s <- mc_result$settings
  cat(sprintf("Settings:\n"))
  cat(sprintf("  Replications: %d\n", s$n_sim))
  cat(sprintf("  Observations: %d\n", s$n_obs))
  cat(sprintf("  Copula: %s\n", toupper(s$copula)))
  cat(sprintf("  Transformation: %s\n", s$transformation))
  cat(sprintf("  True alpha: %.4f\n", s$true_alpha))
  cat(sprintf("  True beta:  %.4f\n", s$true_beta))
  cat(sprintf("  True persistence: %.4f\n", s$true_alpha + s$true_beta))
  if (!is.null(s$true_shape)) {
    cat(sprintf("  True shape: %.1f\n", s$true_shape))
  }
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
  for (param in names(mc_result$true_params)) {
    cat(sprintf("  %s: bias=%.4f, RMSE=%.4f, coverage=%.1f%%\n",
                toupper(param),
                mc_result$bias[param], mc_result$rmse[param],
                100 * mc_result$coverage[param]))
  }
  cat("\n")
  
  cat("Standard Error Calibration:\n")
  for (param in names(mc_result$true_params)) {
    ratio <- mc_result$mean_se[param] / mc_result$empirical_sd[param]
    cat(sprintf("  %s: mean_SE=%.4f, empirical_SD=%.4f, ratio=%.2f\n",
                toupper(param),
                mc_result$mean_se[param], mc_result$empirical_sd[param],
                ratio))
  }
  cat("\n")
  
  invisible(mc_result)
}


#### ______________________________________________________________________ ####
#### SECTION 3: CGARCH Likelihood Surface Visualization                     ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


#' @title Compute NLL Surface for CGARCH(1,1)
#' @description Compute the negative log-likelihood on a grid of (alpha, beta) values
#'   for a Copula GARCH model.
#' @param z_matrix T x k matrix of copula residuals
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional correlation matrix
#' @param alpha_range Range for alpha (default c(0.01, 0.20))
#' @param beta_range Range for beta (default c(0.70, 0.98))
#' @param n_grid Number of grid points per dimension (default 50)
#' @param copula Copula distribution ("mvn" or "mvt")
#' @param shape Degrees of freedom for MVT (only used if copula="mvt")
#' @return List with alpha_grid, beta_grid, nll_surface, and mle location
#' @export
compute_cgarch_nll_surface <- function(
    z_matrix,
    weights,
    Qbar,
    alpha_range = c(0.01, 0.20),
    beta_range = c(0.70, 0.98),
    n_grid = 50,
    copula = "mvn",
    shape = 8
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
        if (copula == "mvn") {
          copula_nll_mvn(
            params = params,
            z_matrix = z_matrix,
            weights = weights,
            Qbar = Qbar,
            copula = "mvn"
          )
        } else {
          copula_nll_mvt(
            params = c(params, shape),
            z_matrix = z_matrix,
            weights = weights,
            Qbar = Qbar
          )
        }
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
    mle_nll = mle_nll,
    copula = copula,
    shape = shape
  )
}


#' @title Plot CGARCH NLL Contours
#' @description Create an interactive contour plot of the CGARCH likelihood surface.
#' @param surface Result from compute_cgarch_nll_surface()
#' @param true_params Optional vector c(alpha, beta) of true parameters
#' @param mle_params Optional vector c(alpha, beta) of MLE estimates
#' @param title Plot title
#' @return plotly object
#' @export
plot_cgarch_nll_contours <- function(
    surface,
    true_params = NULL,
    mle_params = NULL,
    title = NULL
) {
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for this function")
  }
  
  if (is.null(title)) {
    title <- sprintf("CGARCH(1,1) %s Copula Likelihood Surface",
                     toupper(surface$copula))
  }
  
  ## Create base contour plot
  p <- plotly::plot_ly(
    x = surface$alpha_grid,
    y = surface$beta_grid,
    z = t(surface$nll_surface),
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
    name = "\u03b1 + \u03b2 = 1",
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
    xaxis = list(title = "\u03b1 (alpha)"),
    yaxis = list(title = "\u03b2 (beta)"),
    showlegend = TRUE
  )
  
  p
}


#### ______________________________________________________________________ ####
#### SECTION 4: CGARCH Parameter Trajectory Extraction                      ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


#' @title Extract CGARCH Parameter Trajectory from MS Diagnostics
#' @description Extract the evolution of CGARCH parameters (alpha, beta, gamma, shape)
#'   across EM iterations from an ms_diagnostics object.
#' @param diagnostics An object of class \code{ms_diagnostics} from fit_ms_varma_garch()
#' @param state Integer state index (default 1)
#' @return Data frame with columns: iteration, alpha, beta, gamma (if ADCC), shape (if MVT)
#'   Returns NULL if CGARCH parameters not found.
#' @export
extract_cgarch_trajectory <- function(diagnostics, state = 1) {
  
  if (!inherits(diagnostics, "ms_diagnostics")) {
    stop("diagnostics must be an object of class 'ms_diagnostics'")
  }
  
  state_key <- paste0("state_", state)
  state_data <- diagnostics$parameter_evolution[[state_key]]
  
  if (is.null(state_data) || length(state_data) == 0) {
    warning("No parameter evolution data found for state ", state)
    return(NULL)
  }
  
  ## Check model type
  first_params <- state_data[[1]]$parameters
  model_type <- .detect_model_type_from_params(first_params)
  
  if (model_type != "cgarch") {
    warning("State ", state, " does not appear to be a CGARCH model")
    return(NULL)
  }
  
  ## Extract trajectories
  n_iter <- length(state_data)
  
  result <- data.frame(
    iteration = sapply(state_data, function(x) x$iteration),
    alpha = sapply(state_data, function(x) {
      p <- x$parameters
      p$alpha_1 %||% p$dcc_pars$alpha_1 %||% NA
    }),
    beta = sapply(state_data, function(x) {
      p <- x$parameters
      p$beta_1 %||% p$dcc_pars$beta_1 %||% NA
    })
  )
  
  ## Add persistence
  result$persistence <- result$alpha + result$beta
  
  ## Check for gamma (ADCC)
  gamma_vals <- sapply(state_data, function(x) {
    p <- x$parameters
    p$gamma_1 %||% p$dcc_pars$gamma_1 %||% p$adcc_pars$gamma_1 %||% NA
  })
  if (!all(is.na(gamma_vals))) {
    result$gamma <- gamma_vals
    result$persistence <- result$alpha + result$beta + 0.5 * result$gamma
  }
  
  ## Check for shape (MVT copula)
  shape_vals <- sapply(state_data, function(x) {
    p <- x$parameters
    p$shape %||% p$dist_pars$shape %||% NA
  })
  if (!all(is.na(shape_vals))) {
    result$shape <- shape_vals
  }
  
  ## Check for copula and transformation info
  result$copula <- sapply(state_data, function(x) {
    x$parameters$copula %||% NA
  })
  result$transformation <- sapply(state_data, function(x) {
    x$parameters$transformation %||% x$parameters$pit_method %||% NA
  })
  
  result
}


#' @title Extract CGARCH Trajectory with Log-Likelihood
#' @description Extract CGARCH parameter evolution along with log-likelihood values.
#' @param diagnostics An object of class \code{ms_diagnostics}
#' @param state Integer state index (default 1)
#' @return Data frame with columns: iteration, alpha, beta, gamma, shape, log_lik
#' @export
extract_cgarch_trajectory_with_ll <- function(diagnostics, state = 1) {
  
  ## Get parameter trajectory
  traj <- extract_cgarch_trajectory(diagnostics, state = state)
  
  if (is.null(traj)) return(NULL)
  
  ## Get log-likelihood trajectory
  ll_after <- extract_ll_trajectory(diagnostics, type = "after")
  
  ## Match iterations
  n_iter <- nrow(traj)
  if (length(ll_after) >= n_iter) {
    traj$log_lik <- ll_after[1:n_iter]
  } else {
    traj$log_lik <- c(ll_after, rep(NA, n_iter - length(ll_after)))
  }
  
  traj
}


#### ______________________________________________________________________ ####
#### SECTION 5: CGARCH Complete Diagnostics                                 ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


#' @title Summary Function for Complete CGARCH Diagnostic Analysis
#' @description Run a complete diagnostic analysis including MC study,
#'   likelihood surface, and statistical tests for CGARCH models.
#' @param n_sim Number of MC replications
#' @param n_obs Observations per replication
#' @param true_alpha True DCC alpha
#' @param true_beta True DCC beta
#' @param omega Vector of GARCH omega parameters
#' @param alpha_garch Vector of GARCH alpha parameters
#' @param beta_garch Vector of GARCH beta parameters
#' @param copula Copula distribution ("mvn" or "mvt")
#' @param true_shape True shape for MVT copula
#' @param transformation PIT transformation type
#' @param seed Random seed
#' @param plot Logical: create plots?
#' @return List with all diagnostic results
#' @export
run_complete_cgarch_diagnostics <- function(
    n_sim = 50,
    n_obs = 500,
    true_alpha = 0.05,
    true_beta = 0.90,
    omega = c(0.05, 0.08),
    alpha_garch = c(0.10, 0.12),
    beta_garch = c(0.85, 0.82),
    copula = "mvn",
    true_shape = 8,
    transformation = "parametric",
    seed = 42,
    plot = TRUE
) {
  
  k <- length(omega)
  
  cat("\n")
  cat("================================================================\n")
  cat("  Complete CGARCH(1,1) Estimation Diagnostics\n")
  cat("================================================================\n\n")
  
  ## 1. Run Monte Carlo study
  cat("Step 1: Running Monte Carlo simulation...\n")
  mc_result <- run_cgarch_monte_carlo(
    n_sim = n_sim,
    n_obs = n_obs,
    k = k,
    true_alpha = true_alpha,
    true_beta = true_beta,
    omega = omega,
    alpha_garch = alpha_garch,
    beta_garch = beta_garch,
    copula = copula,
    true_shape = true_shape,
    transformation = transformation,
    seed = seed,
    verbose = TRUE
  )
  
  ## 2. Generate one dataset for surface visualization
  cat("\nStep 2: Computing likelihood surface...\n")
  y <- simulate_cgarch(
    n = n_obs,
    k = k,
    omega = omega,
    alpha_garch = alpha_garch,
    beta_garch = beta_garch,
    alpha_dcc = true_alpha,
    beta_dcc = true_beta,
    copula = copula,
    shape = true_shape,
    seed = seed
  )
  
  resid_result <- compute_cgarch_residuals(
    y, omega, alpha_garch, beta_garch,
    transformation = transformation,
    copula = copula,
    shape = true_shape
  )
  
  z_matrix <- resid_result$z_matrix
  
  surface <- compute_cgarch_nll_surface(
    z_matrix = z_matrix,
    weights = rep(1, n_obs),
    Qbar = cor(z_matrix),
    alpha_range = c(max(0.001, true_alpha - 0.10), min(0.30, true_alpha + 0.15)),
    beta_range = c(max(0.50, true_beta - 0.20), min(0.99, true_beta + 0.05)),
    n_grid = 40,
    copula = copula,
    shape = true_shape
  )
  
  ## 3. Statistical tests
  cat("\nStep 3: Running statistical tests...\n")
  normality_tests <- test_estimate_normality_generic(mc_result, param_names = c("alpha", "beta"))
  
  ## 4. Create plots
  plots <- NULL
  if (plot && requireNamespace("plotly", quietly = TRUE)) {
    cat("\nStep 4: Creating visualizations...\n")
    
    plots <- list(
      contour = plot_cgarch_nll_contours(
        surface,
        true_params = c(true_alpha, true_beta),
        title = sprintf("CGARCH(1,1) %s Copula Likelihood Surface", toupper(copula))
      ),
      mc_scatter = plot_mc_estimates_generic(mc_result, model_type = "cgarch")$scatter,
      alpha_hist = plot_mc_estimates_generic(mc_result, model_type = "cgarch")$alpha_hist,
      beta_hist = plot_mc_estimates_generic(mc_result, model_type = "cgarch")$beta_hist
    )
    
    cat("  Plots created. Access via result$plots\n")
  }
  
  cat("\n================================================================\n")
  cat("  CGARCH Diagnostics Complete\n")
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
      true_shape = if (copula == "mvt") true_shape else NULL,
      copula = copula,
      transformation = transformation,
      seed = seed
    )
  )
}


#### ______________________________________________________________________ ####
#### SECTION 6: GOGARCH Monte Carlo Simulation Study                        ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


#' @title Run Monte Carlo Study for GOGARCH Estimation
#' @description Performs a Monte Carlo simulation study to assess the accuracy
#'   of GOGARCH parameter estimation. Since GOGARCH estimates component-wise
#'   GARCH parameters, this evaluates estimation accuracy for each component.
#' @param n_sim Number of simulation replications
#' @param n_obs Number of observations per replication
#' @param k Number of series/components (default 3)
#' @param true_omega Vector of true component GARCH omega parameters (alias for omega)
#' @param true_alpha Vector of true component GARCH alpha parameters (alias for alpha_garch)
#' @param true_beta Vector of true component GARCH beta parameters (alias for beta_garch)
#' @param omega Vector of component GARCH omega parameters (alternative to true_omega)
#' @param alpha_garch Vector of component GARCH alpha parameters (alternative to true_alpha)
#' @param beta_garch Vector of component GARCH beta parameters (alternative to true_beta)
#' @param distribution Component distribution ("norm" or "std")
#' @param shape Degrees of freedom for "std" distribution
#' @param ica_method ICA algorithm ("radical" or "fastica")
#' @param confidence_level Confidence level for coverage (default 0.95)
#' @param verbose Print progress
#' @param seed Base seed for reproducibility
#' @return List with:
#'   \item{estimates}{List of data frames (one per component) with alpha, beta columns}
#'   \item{persistence}{List of vectors (one per component)}
#'   \item{convergence}{Logical vector indicating optimization convergence}
#'   \item{ica_converged}{Logical vector indicating ICA convergence}
#'   \item{bias}{List of named vectors (one per component)}
#'   \item{rmse}{List of named vectors (one per component)}
#'   \item{empirical_sd}{List of named vectors (one per component)}
#'   \item{coverage}{List of named vectors (one per component)}
#'   \item{mixing_recovery}{Statistics on mixing matrix recovery}
#'   \item{summary}{Summary data frame}
#' @export
run_gogarch_monte_carlo <- function(
    n_sim = 100,
    n_obs = 500,
    k = 3,
    true_omega = NULL,
    true_alpha = NULL,
    true_beta = NULL,
    omega = NULL,
    alpha_garch = NULL,
    beta_garch = NULL,
    distribution = "norm",
    shape = 8,
    ica_method = "radical",
    confidence_level = 0.95,
    verbose = TRUE,
    seed = 12345
) {
  
  set.seed(seed)
  
  ## Handle parameter name aliases (true_* takes precedence)
  if (!is.null(true_omega)) omega <- true_omega
  if (!is.null(true_alpha)) alpha_garch <- true_alpha
  if (!is.null(true_beta)) beta_garch <- true_beta
  
  ## Default GARCH parameters (different for each component)
  if (is.null(omega)) omega <- seq(0.03, 0.08, length.out = k)
  if (is.null(alpha_garch)) alpha_garch <- seq(0.05, 0.12, length.out = k)
  if (is.null(beta_garch)) beta_garch <- seq(0.90, 0.85, length.out = k)
  
  ## Ensure vectors have length k
  if (length(omega) < k) omega <- rep(omega[1], k)
  if (length(alpha_garch) < k) alpha_garch <- rep(alpha_garch[1], k)
  if (length(beta_garch) < k) beta_garch <- rep(beta_garch[1], k)
  
  ## Storage for each component - using the structure expected by the Rmd
  estimates <- vector("list", k)
  persistence <- vector("list", k)
  std_errors <- vector("list", k)
  
  for (c in 1:k) {
    estimates[[c]] <- data.frame(
      omega = numeric(n_sim),
      alpha = numeric(n_sim),
      beta = numeric(n_sim)
    )
    persistence[[c]] <- numeric(n_sim)
    std_errors[[c]] <- data.frame(
      omega = numeric(n_sim),
      alpha = numeric(n_sim),
      beta = numeric(n_sim)
    )
  }
  
  ## Storage for mixing matrix recovery and convergence
  mixing_angles <- numeric(n_sim)
  mixing_cond_numbers <- numeric(n_sim)
  convergence <- logical(n_sim)
  ica_converged <- logical(n_sim)
  valid_se <- logical(n_sim)
  
  ## True parameters
  true_params <- list(
    omega = omega,
    alpha = alpha_garch,
    beta = beta_garch,
    persistence = alpha_garch + beta_garch
  )
  
  if (verbose) {
    cat("Running GOGARCH Monte Carlo study:\n")
    cat(sprintf("  Replications: %d\n", n_sim))
    cat(sprintf("  Observations: %d\n", n_obs))
    cat(sprintf("  Components: %d\n", k))
    cat(sprintf("  ICA method: %s\n", ica_method))
    cat(sprintf("  Distribution: %s\n", distribution))
    cat("\n  True component parameters:\n")
    for (c in 1:k) {
      cat(sprintf("    Component %d: omega=%.4f, alpha=%.4f, beta=%.4f, persistence=%.4f\n",
                  c, omega[c], alpha_garch[c], beta_garch[c], alpha_garch[c] + beta_garch[c]))
    }
    cat("\n")
  }
  
  pb <- if (verbose) txtProgressBar(min = 0, max = n_sim, style = 3) else NULL
  
  for (i in 1:n_sim) {
    if (verbose) setTxtProgressBar(pb, i)
    
    tryCatch({
      ## Simulate GOGARCH data
      sim_result <- simulate_gogarch(
        n = n_obs,
        k = k,
        omega = omega,
        alpha_garch = alpha_garch,
        beta_garch = beta_garch,
        distribution = distribution,
        shape = shape,
        seed = seed + i
      )
      
      y <- sim_result$y
      true_A <- sim_result$A
      
      ## Fit GOGARCH model
      fit_result <- fit_gogarch_model(
        y = y,
        k = k,
        ica_method = ica_method,
        distribution = distribution
      )
      
      if (!is.null(fit_result) && fit_result$convergence) {
        convergence[i] <- TRUE
        ica_converged[i] <- fit_result$ica_converged
        
        ## Extract component estimates
        for (c in 1:k) {
          if (!is.null(fit_result$garch_pars[[c]])) {
            estimates[[c]]$omega[i] <- fit_result$garch_pars[[c]]$omega %||% NA
            estimates[[c]]$alpha[i] <- fit_result$garch_pars[[c]]$alpha1 %||% 
              fit_result$garch_pars[[c]]$alpha %||% NA
            estimates[[c]]$beta[i] <- fit_result$garch_pars[[c]]$beta1 %||% 
              fit_result$garch_pars[[c]]$beta %||% NA
            
            persistence[[c]][i] <- estimates[[c]]$alpha[i] + estimates[[c]]$beta[i]
            
            ## Extract standard errors if available
            if (!is.null(fit_result$std_errors) && !is.null(fit_result$std_errors[[c]])) {
              std_errors[[c]]$omega[i] <- fit_result$std_errors[[c]]$omega %||% NA
              std_errors[[c]]$alpha[i] <- fit_result$std_errors[[c]]$alpha1 %||% 
                fit_result$std_errors[[c]]$alpha %||% NA
              std_errors[[c]]$beta[i] <- fit_result$std_errors[[c]]$beta1 %||% 
                fit_result$std_errors[[c]]$beta %||% NA
            }
          }
        }
        
        ## Mixing matrix recovery metrics
        if (!is.null(fit_result$A)) {
          ## Compute angle between estimated and true mixing matrices
          mixing_angles[i] <- compute_mixing_angle(fit_result$A, true_A)
          mixing_cond_numbers[i] <- kappa(fit_result$A)
        }
        
        ## Check if SEs are valid
        valid_se[i] <- all(sapply(1:k, function(c) {
          all(is.finite(c(std_errors[[c]]$alpha[i], std_errors[[c]]$beta[i])))
        }))
      }
      
    }, error = function(e) {
      ## Silently continue on error
    })
  }
  
  if (verbose) {
    close(pb)
    cat("\n")
  }
  
  ## Compute summary statistics for each component
  valid_idx <- convergence
  n_valid <- sum(valid_idx)
  
  bias <- vector("list", k)
  rmse <- vector("list", k)
  empirical_sd <- vector("list", k)
  coverage <- vector("list", k)
  mean_se <- vector("list", k)
  
  z_crit <- qnorm(1 - (1 - confidence_level) / 2)
  
  for (c in 1:k) {
    ## Filter valid estimates
    alpha_est <- estimates[[c]]$alpha[valid_idx]
    beta_est <- estimates[[c]]$beta[valid_idx]
    omega_est <- estimates[[c]]$omega[valid_idx]
    
    ## Bias
    bias[[c]] <- c(
      omega = mean(omega_est, na.rm = TRUE) - omega[c],
      alpha = mean(alpha_est, na.rm = TRUE) - alpha_garch[c],
      beta = mean(beta_est, na.rm = TRUE) - beta_garch[c]
    )
    
    ## RMSE
    rmse[[c]] <- c(
      omega = sqrt(mean((omega_est - omega[c])^2, na.rm = TRUE)),
      alpha = sqrt(mean((alpha_est - alpha_garch[c])^2, na.rm = TRUE)),
      beta = sqrt(mean((beta_est - beta_garch[c])^2, na.rm = TRUE))
    )
    
    ## Empirical SD
    empirical_sd[[c]] <- c(
      omega = sd(omega_est, na.rm = TRUE),
      alpha = sd(alpha_est, na.rm = TRUE),
      beta = sd(beta_est, na.rm = TRUE)
    )
    
    ## Mean SE
    valid_se_idx <- valid_idx & valid_se
    mean_se[[c]] <- c(
      omega = mean(std_errors[[c]]$omega[valid_se_idx], na.rm = TRUE),
      alpha = mean(std_errors[[c]]$alpha[valid_se_idx], na.rm = TRUE),
      beta = mean(std_errors[[c]]$beta[valid_se_idx], na.rm = TRUE)
    )
    
    ## Coverage (if SEs are available)
    if (sum(valid_se_idx) > 0) {
      alpha_se <- std_errors[[c]]$alpha[valid_se_idx]
      beta_se <- std_errors[[c]]$beta[valid_se_idx]
      alpha_est_se <- estimates[[c]]$alpha[valid_se_idx]
      beta_est_se <- estimates[[c]]$beta[valid_se_idx]
      
      alpha_covered <- (alpha_est_se - z_crit * alpha_se <= alpha_garch[c]) & 
        (alpha_est_se + z_crit * alpha_se >= alpha_garch[c])
      beta_covered <- (beta_est_se - z_crit * beta_se <= beta_garch[c]) & 
        (beta_est_se + z_crit * beta_se >= beta_garch[c])
      
      coverage[[c]] <- c(
        alpha = mean(alpha_covered, na.rm = TRUE),
        beta = mean(beta_covered, na.rm = TRUE)
      )
    } else {
      coverage[[c]] <- c(alpha = NA, beta = NA)
    }
  }
  
  ## Mixing matrix recovery summary
  mixing_recovery <- list(
    mean_angle = mean(mixing_angles[valid_idx], na.rm = TRUE),
    sd_angle = sd(mixing_angles[valid_idx], na.rm = TRUE),
    mean_cond_number = mean(mixing_cond_numbers[valid_idx], na.rm = TRUE)
  )
  
  ## Print summary
  if (verbose) {
    cat("\n=== GOGARCH Monte Carlo Results ===\n")
    cat(sprintf("Valid replications: %d / %d (%.1f%%)\n", 
                n_valid, n_sim, 100 * n_valid / n_sim))
    cat(sprintf("ICA converged: %d / %d (%.1f%%)\n",
                sum(ica_converged), n_sim, 100 * sum(ica_converged) / n_sim))
    cat(sprintf("Valid SEs: %d / %d (%.1f%%)\n",
                sum(valid_se), n_sim, 100 * sum(valid_se) / n_sim))
    cat("\nComponent-wise bias:\n")
    for (c in 1:k) {
      cat(sprintf("  Component %d: alpha=%.4f, beta=%.4f\n",
                  c, bias[[c]]["alpha"], bias[[c]]["beta"]))
    }
    cat(sprintf("\nMixing matrix recovery: mean angle=%.4f rad, mean cond=%.2f\n",
                mixing_recovery$mean_angle, mixing_recovery$mean_cond_number))
  }
  
  ## Return structure matching Rmd expectations
  list(
    estimates = estimates,
    std_errors = std_errors,
    persistence = persistence,
    convergence = convergence,
    ica_converged = ica_converged,
    valid_se = valid_se,
    true_params = true_params,
    bias = bias,
    rmse = rmse,
    empirical_sd = empirical_sd,
    mean_se = mean_se,
    coverage = coverage,
    mixing_recovery = mixing_recovery,
    mixing_angles = mixing_angles,
    mixing_cond_numbers = mixing_cond_numbers,
    settings = list(
      n_sim = n_sim,
      n_obs = n_obs,
      k = k,
      omega = omega,
      alpha_garch = alpha_garch,
      beta_garch = beta_garch,
      distribution = distribution,
      shape = if (distribution == "std") shape else NULL,
      ica_method = ica_method,
      confidence_level = confidence_level
    )
  )
}


#' @title Compute Mixing Matrix Recovery Angle
#' @description Compute angle between true and estimated mixing matrix subspaces.
#'   Uses the Amari index or similar metric for ICA recovery assessment.
#' @param A_true True mixing matrix
#' @param A_est Estimated mixing matrix
#' @return Scalar angle metric (0 = perfect recovery, larger = worse)
#' @keywords internal
compute_mixing_recovery_angle <- function(A_true, A_est) {
  ## Compute normalized cross-product
  P <- abs(t(A_true) %*% A_est)
  
  ## Normalize rows and columns
  P_row_norm <- P / pmax(apply(P, 1, max), 1e-10)
  P_col_norm <- P / pmax(apply(P, 2, max), 1e-10)
  
  ## Amari index
  k <- nrow(P)
  amari <- (sum(P_row_norm) - k + sum(P_col_norm) - k) / (2 * k * (k - 1))
  
  amari
}


#' @title Print GOGARCH Monte Carlo Results
#' @export
print_gogarch_mc_results <- function(mc_result) {
  cat("\n")
  cat("============================================\n")
  cat("  GOGARCH Monte Carlo Simulation Study\n")
  cat("============================================\n\n")
  
  s <- mc_result$settings
  cat(sprintf("Settings:\n"))
  cat(sprintf("  Replications: %d\n", s$n_sim))
  cat(sprintf("  Observations: %d\n", s$n_obs))
  cat(sprintf("  Components: %d\n", s$k))
  cat(sprintf("  ICA method: %s\n", s$ica_method))
  cat(sprintf("  Distribution: %s\n", s$distribution))
  cat("\n")
  
  cat("Results:\n")
  cat(sprintf("  Converged: %d / %d (%.1f%%)\n",
              sum(mc_result$convergence), s$n_sim,
              100 * mean(mc_result$convergence)))
  cat("\n")
  
  cat("Component GARCH Estimation Summary:\n")
  print(mc_result$summary, digits = 4)
  cat("\n")
  
  cat("Mixing Matrix Recovery:\n")
  ms <- mc_result$mixing_summary
  cat(sprintf("  Mean Amari index: %.4f (SD: %.4f)\n", ms$mean_angle, ms$sd_angle))
  cat(sprintf("  Interpretation: %.1f%% recovery error\n", 100 * ms$mean_angle))
  cat(sprintf("  Mean condition number: %.2f (SD: %.2f)\n",
              ms$mean_cond_number, ms$sd_cond_number))
  cat("\n")
  
  invisible(mc_result)
}


#### ______________________________________________________________________ ####
#### SECTION 7: GOGARCH Parameter Trajectory and Visualization              ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


#' @title Visualize GOGARCH Component Evolution
#' @description Plot the evolution of GOGARCH component parameters across EM iterations.
#' @param diagnostics An object of class \code{ms_diagnostics}
#' @param state Integer state index (default 1)
#' @param show_persistence Logical: show persistence rather than alpha/beta separately
#' @return plotly object with component trajectory plot
#' @export
visualize_gogarch_evolution <- function(
    diagnostics,
    state = 1,
    show_persistence = TRUE
) {
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for this function")
  }
  
  ## Extract trajectory
  traj <- extract_gogarch_trajectory(diagnostics, state = state)
  
  if (is.null(traj)) {
    stop("Could not extract GOGARCH trajectory for state ", state)
  }
  
  n_iter <- nrow(traj)
  
  ## Get component persistence columns
  persist_cols <- grep("persistence", names(traj), value = TRUE)
  n_components <- length(persist_cols)
  
  if (n_components == 0) {
    stop("No component persistence data found")
  }
  
  ## Create plot
  p <- plotly::plot_ly()
  
  colors <- RColorBrewer::brewer.pal(max(3, n_components), "Set1")
  
  for (c in 1:n_components) {
    col_name <- paste0("component_", c, "_persistence")
    
    if (col_name %in% names(traj)) {
      p <- p %>% plotly::add_trace(
        x = traj$iteration,
        y = traj[[col_name]],
        type = "scatter",
        mode = "lines+markers",
        name = sprintf("Component %d", c),
        line = list(color = colors[c], width = 2),
        marker = list(color = colors[c], size = 6)
      )
    }
  }
  
  ## Add stationarity threshold
  p <- p %>% plotly::add_trace(
    x = c(1, n_iter),
    y = c(1, 1),
    type = "scatter",
    mode = "lines",
    name = "Stationarity Bound",
    line = list(color = "red", width = 2, dash = "dash")
  )
  
  ## Layout
  p <- p %>% plotly::layout(
    title = sprintf("GOGARCH Component Persistence Evolution (State %d)", state),
    xaxis = list(title = "EM Iteration"),
    yaxis = list(title = "Persistence (\u03b1 + \u03b2)", range = c(0, 1.1)),
    showlegend = TRUE
  )
  
  p
}


#' @title Plot GOGARCH Mixing Matrix Heatmap
#' @description Visualize the estimated mixing matrix from GOGARCH.
#' @param A Mixing matrix (k x k)
#' @param title Plot title
#' @return plotly object with heatmap
#' @export
plot_mixing_matrix <- function(A, title = "GOGARCH Mixing Matrix") {
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for this function")
  }
  
  k <- nrow(A)
  
  p <- plotly::plot_ly(
    x = paste0("Component ", 1:k),
    y = paste0("Series ", 1:k),
    z = A,
    type = "heatmap",
    colorscale = "RdBu",
    zmid = 0
  )
  
  p <- p %>% plotly::layout(
    title = title,
    xaxis = list(title = "Independent Component"),
    yaxis = list(title = "Observed Series")
  )
  
  p
}


#' @title Complete GOGARCH Diagnostics
#' @description Run comprehensive diagnostics for GOGARCH models.
#' @param n_sim Number of MC replications
#' @param n_obs Observations per replication
#' @param k Number of series/components
#' @param omega Vector of component GARCH omega
#' @param alpha_garch Vector of component GARCH alpha
#' @param beta_garch Vector of component GARCH beta
#' @param distribution Component distribution
#' @param ica_method ICA algorithm
#' @param seed Random seed
#' @param plot Logical: create plots?
#' @return List with diagnostic results
#' @export
run_complete_gogarch_diagnostics <- function(
    n_sim = 50,
    n_obs = 500,
    k = 3,
    omega = NULL,
    alpha_garch = NULL,
    beta_garch = NULL,
    distribution = "norm",
    ica_method = "radical",
    seed = 42,
    plot = TRUE
) {
  
  ## Default parameters
  if (is.null(omega)) omega <- seq(0.03, 0.08, length.out = k)
  if (is.null(alpha_garch)) alpha_garch <- seq(0.05, 0.12, length.out = k)
  if (is.null(beta_garch)) beta_garch <- seq(0.90, 0.85, length.out = k)
  
  cat("\n")
  cat("================================================================\n")
  cat("  Complete GOGARCH Estimation Diagnostics\n")
  cat("================================================================\n\n")
  
  ## 1. Run Monte Carlo study
  cat("Step 1: Running Monte Carlo simulation...\n")
  mc_result <- run_gogarch_monte_carlo(
    n_sim = n_sim,
    n_obs = n_obs,
    k = k,
    omega = omega,
    alpha_garch = alpha_garch,
    beta_garch = beta_garch,
    distribution = distribution,
    ica_method = ica_method,
    seed = seed,
    verbose = TRUE
  )
  
  ## 2. Create plots
  plots <- NULL
  if (plot && requireNamespace("plotly", quietly = TRUE)) {
    cat("\nStep 2: Creating visualizations...\n")
    
    valid_idx <- mc_result$convergence
    
    ## Component persistence histograms
    hist_plots <- list()
    for (c in 1:k) {
      persist_c <- mc_result$component_estimates[[c]]$persistence[valid_idx]
      true_persist <- alpha_garch[c] + beta_garch[c]
      
      hist_plots[[paste0("component_", c)]] <- plotly::plot_ly(
        x = persist_c,
        type = "histogram",
        name = paste("Component", c),
        marker = list(color = "steelblue")
      ) %>%
        plotly::add_trace(
          x = c(true_persist, true_persist),
          y = c(0, max(table(cut(persist_c, 20)))),
          type = "scatter",
          mode = "lines",
          line = list(color = "red", width = 2, dash = "dash"),
          name = sprintf("True (%.3f)", true_persist)
        ) %>%
        plotly::layout(
          title = sprintf("Component %d Persistence", c),
          xaxis = list(title = "Persistence"),
          yaxis = list(title = "Count")
        )
    }
    
    ## Mixing recovery scatter
    mixing_plot <- plotly::plot_ly(
      x = 1:sum(valid_idx),
      y = mc_result$mixing_angles[valid_idx],
      type = "scatter",
      mode = "markers",
      marker = list(color = "steelblue", size = 6),
      name = "Amari Index"
    ) %>%
      plotly::layout(
        title = "Mixing Matrix Recovery (Amari Index)",
        xaxis = list(title = "Replication"),
        yaxis = list(title = "Amari Index (lower = better)")
      )
    
    plots <- c(hist_plots, list(mixing_recovery = mixing_plot))
    
    cat("  Plots created. Access via result$plots\n")
  }
  
  cat("\n================================================================\n")
  cat("  GOGARCH Diagnostics Complete\n")
  cat("================================================================\n")
  
  list(
    mc_result = mc_result,
    plots = plots,
    settings = list(
      n_sim = n_sim,
      n_obs = n_obs,
      k = k,
      omega = omega,
      alpha_garch = alpha_garch,
      beta_garch = beta_garch,
      distribution = distribution,
      ica_method = ica_method,
      seed = seed
    )
  )
}


#### ______________________________________________________________________ ####
#### SECTION 8: Generic Helper Functions for Both Models                    ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


#' @title Test Normality of Standardized Estimates (Generic)
#' @description Test whether (theta_hat - theta_0) / SE follows N(0,1).
#'   Works for both CGARCH and DCC results.
#' @param mc_result Result from run_cgarch_monte_carlo() or run_dcc_monte_carlo()
#' @param param_names Names of parameters to test (default: c("alpha", "beta"))
#' @return List with test statistics for each parameter
#' @export
test_estimate_normality_generic <- function(mc_result, param_names = c("alpha", "beta")) {
  
  valid_idx <- mc_result$convergence & mc_result$valid_se
  
  results <- list()
  
  for (param in param_names) {
    j <- which(names(mc_result$true_params) == param)
    if (length(j) == 0) next
    
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
  cat("Under correct specification: (\u03b8\u0302 - \u03b8\u2080) / SE ~ N(0, 1)\n\n")
  
  for (param in names(results)) {
    r <- results[[param]]
    cat(sprintf("%s:\n", toupper(param)))
    cat(sprintf("  N valid: %d\n", r$n))
    cat(sprintf("  Mean z-score: %.3f (should be \u2248 0)\n", r$mean_z))
    cat(sprintf("  SD z-score: %.3f (should be \u2248 1)\n", r$sd_z))
    cat(sprintf("  Skewness: %.3f (should be \u2248 0)\n", r$skewness))
    cat(sprintf("  Excess kurtosis: %.3f (should be \u2248 0)\n", r$kurtosis))
    if (!is.na(r$shapiro_p)) {
      cat(sprintf("  Shapiro-Wilk p-value: %.4f\n", r$shapiro_p))
    }
    cat("\n")
  }
  
  invisible(results)
}


#' @title Plot Monte Carlo Estimates (Generic)
#' @description Create histograms and scatter plots for MC estimates.
#'   Works for both CGARCH and DCC results.
#' @param mc_result Result from Monte Carlo study
#' @param model_type Model type identifier for titles
#' @return List with plotly objects
#' @export
plot_mc_estimates_generic <- function(mc_result, model_type = "CGARCH") {
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for this function")
  }
  
  valid_idx <- mc_result$convergence & !is.na(mc_result$estimates[, 1])
  alpha_est <- mc_result$estimates[valid_idx, "alpha"]
  beta_est <- mc_result$estimates[valid_idx, "beta"]
  
  true_alpha <- mc_result$true_params["alpha"]
  true_beta <- mc_result$true_params["beta"]
  
  ## Scatter plot
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
      title = sprintf("%s Parameter Estimates", model_type),
      xaxis = list(title = "\u03b1"),
      yaxis = list(title = "\u03b2")
    )
  
  ## Alpha histogram
  p2 <- plotly::plot_ly(
    x = alpha_est,
    type = "histogram",
    name = "\u03b1 estimates",
    marker = list(color = "steelblue")
  ) %>%
    plotly::add_trace(
      x = c(true_alpha, true_alpha),
      y = c(0, max(table(cut(alpha_est, 30)))),
      type = "scatter",
      mode = "lines",
      line = list(color = "red", width = 2, dash = "dash"),
      name = "True \u03b1"
    ) %>%
    plotly::layout(
      title = sprintf("%s Alpha Distribution", model_type),
      xaxis = list(title = "\u03b1"),
      yaxis = list(title = "Count")
    )
  
  ## Beta histogram
  p3 <- plotly::plot_ly(
    x = beta_est,
    type = "histogram",
    name = "\u03b2 estimates",
    marker = list(color = "steelblue")
  ) %>%
    plotly::add_trace(
      x = c(true_beta, true_beta),
      y = c(0, max(table(cut(beta_est, 30)))),
      type = "scatter",
      mode = "lines",
      line = list(color = "red", width = 2, dash = "dash"),
      name = "True \u03b2"
    ) %>%
    plotly::layout(
      title = sprintf("%s Beta Distribution", model_type),
      xaxis = list(title = "\u03b2"),
      yaxis = list(title = "Count")
    )
  
  list(scatter = p1, alpha_hist = p2, beta_hist = p3)
}


## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## STUB FUNCTIONS - These reference external functions that need implementation
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#' @title Copula NLL Wrapper (MVN)
#' @description Wrapper for copula negative log-likelihood with MVN copula.
#'   This references the function in tsbs_cgarch.R
#' @keywords internal
copula_nll_mvn <- function(params, z_matrix, weights, Qbar, copula = "mvn") {
  ## This should call the actual copula_nll_mvn from tsbs_cgarch.R
  ## Placeholder implementation for standalone testing
  alpha <- params[1]
  beta <- params[2]
  
  if (alpha < 0 || beta < 0 || alpha + beta >= 1) return(1e10)
  
  n <- nrow(z_matrix)
  k <- ncol(z_matrix)
  
  ## Simple DCC recursion
  Q <- Qbar
  ll_vec <- numeric(n)
  intercept <- (1 - alpha - beta) * Qbar
  
  for (t in 1:n) {
    ## Correlation matrix
    Q_diag_inv_sqrt <- diag(1/sqrt(diag(Q)))
    R <- Q_diag_inv_sqrt %*% Q %*% Q_diag_inv_sqrt
    diag(R) <- 1
    
    ## MVN log-likelihood
    det_R <- det(R)
    if (det_R <= 0) {
      ll_vec[t] <- -1e10
    } else {
      R_inv <- solve(R)
      z_t <- z_matrix[t, ]
      mahal <- as.numeric(t(z_t) %*% R_inv %*% z_t)
      ll_vec[t] <- -0.5 * (log(det_R) + mahal - as.numeric(t(z_t) %*% z_t))
    }
    
    ## Update Q
    if (t < n) {
      z_t <- z_matrix[t, ]
      Q <- intercept + alpha * (z_t %*% t(z_t)) + beta * Q
    }
  }
  
  -sum(weights * ll_vec)
}


#' @title Copula NLL for MVT
#' @description Negative log-likelihood for Student-t copula
#' @keywords internal
copula_nll_mvt <- function(params, z_matrix, weights, Qbar) {
  alpha <- params[1]
  beta <- params[2]
  shape <- params[3]
  
  if (alpha < 0 || beta < 0 || alpha + beta >= 1 || shape <= 2) return(1e10)
  
  n <- nrow(z_matrix)
  k <- ncol(z_matrix)
  
  ## DCC recursion
  Q <- Qbar
  ll_vec <- numeric(n)
  intercept <- (1 - alpha - beta) * Qbar
  
  const_term <- lgamma(0.5 * (k + shape)) - lgamma(0.5 * shape) - 
    0.5 * k * log(pi * (shape - 2))
  
  for (t in 1:n) {
    Q_diag_inv_sqrt <- diag(1/sqrt(pmax(diag(Q), 1e-10)))
    R <- Q_diag_inv_sqrt %*% Q %*% Q_diag_inv_sqrt
    diag(R) <- 1
    
    det_R <- det(R)
    if (det_R <= 0) {
      ll_vec[t] <- -1e10
    } else {
      R_inv <- solve(R)
      z_t <- z_matrix[t, ]
      mahal <- as.numeric(t(z_t) %*% R_inv %*% z_t)
      
      ll_vec[t] <- const_term - 0.5 * log(det_R) - 
        0.5 * (shape + k) * log(1 + mahal / (shape - 2))
      
      ## Subtract marginal t log-densities
      scale <- sqrt(shape / (shape - 2))
      for (j in 1:k) {
        ll_vec[t] <- ll_vec[t] - (dt(z_t[j] * scale, df = shape, log = TRUE) + log(scale))
      }
    }
    
    if (t < n) {
      z_t <- z_matrix[t, ]
      Q <- intercept + alpha * (z_t %*% t(z_t)) + beta * Q
    }
  }
  
  -sum(weights * ll_vec)
}