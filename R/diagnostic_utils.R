#' @title Diagnostic Utility Functions for MS-VARMA-GARCH
#' @description Convenience functions for extracting and analyzing diagnostic information
#' @name diagnostic_utils
NULL

## === EXTRACTION FUNCTIONS ===

#' Extract Parameter Trajectory
#'
#' Extract the evolution of a specific parameter across EM iterations for a given state.
#'
#' @param diagnostics An object of class \code{ms_diagnostics}
#' @param state Integer state index
#' @param param_name Character string specifying parameter name (e.g., "alpha_1", "omega")
#' @param series Optional integer for series-specific parameters (for multivariate models)
#'
#' @return A data.frame with columns \code{iteration} and \code{value}, or NULL if not found
#'
#' @examples
#' \dontrun{
#' # Extract DCC alpha for state 1
#' alpha_traj <- extract_param_trajectory(diag, state = 1, param_name = "alpha_1")
#' plot(alpha_traj$iteration, alpha_traj$value, type = "b")
#' 
#' # Extract omega for state 2, series 1
#' omega_traj <- extract_param_trajectory(diag, state = 2, 
#'                                         param_name = "omega", series = 1)
#' }
#'
#' @export
extract_param_trajectory <- function(diagnostics, state, param_name, series = NULL) {
  if (!inherits(diagnostics, "ms_diagnostics")) {
    stop("diagnostics must be an object of class 'ms_diagnostics'")
  }
  
  state_key <- paste0("state_", state)
  state_data <- diagnostics$parameter_evolution[[state_key]]
  
  if (is.null(state_data) || length(state_data) == 0) {
    warning("No parameter evolution data found for state ", state)
    return(NULL)
  }
  
  trajectory <- sapply(state_data, function(iter_data) {
    params <- iter_data$parameters
    
    # Handle series-specific parameters (e.g., garch_pars[[1]]$omega)
    if (!is.null(series) && !is.null(params$garch_pars)) {
      if (series <= length(params$garch_pars)) {
        if (param_name %in% names(params$garch_pars[[series]])) {
          return(params$garch_pars[[series]][[param_name]])
        }
      }
      return(NA)
    }
    
    # Handle top-level parameters
    if (param_name %in% names(params)) {
      return(params[[param_name]])
    }
    
    # Try flattened search
    params_flat <- unlist(params)
    param_match <- grep(param_name, names(params_flat), value = TRUE, fixed = TRUE)
    if (length(param_match) > 0) {
      return(params_flat[[param_match[1]]])
    }
    
    return(NA)
  })
  
  iterations <- sapply(state_data, function(x) x$iteration)
  
  data.frame(iteration = iterations, value = trajectory)
}


#' Extract Log-Likelihood Trajectory
#'
#' Extract the complete log-likelihood evolution across EM iterations.
#'
#' @param diagnostics An object of class \code{ms_diagnostics}
#' @param type Character: "before", "after", or "change"
#'
#' @return Numeric vector of log-likelihood values
#'
#' @examples
#' \dontrun{
#' ll_after <- extract_ll_trajectory(diag, type = "after")
#' ll_changes <- extract_ll_trajectory(diag, type = "change")
#' plot(ll_after, type = "b", main = "LL Evolution")
#' }
#'
#' @export
extract_ll_trajectory <- function(diagnostics, type = c("after", "before", "change")) {
  if (!inherits(diagnostics, "ms_diagnostics")) {
    stop("diagnostics must be an object of class 'ms_diagnostics'")
  }
  
  type <- match.arg(type)
  
  field_name <- switch(type,
                       "before" = "log_lik_before_mstep",
                       "after" = "log_lik_after_mstep",
                       "change" = "ll_change"
  )
  
  sapply(diagnostics$em_iterations, function(x) x[[field_name]])
}


#' Extract Sigma Evolution Statistics
#'
#' Extract volatility evolution summary for a specific state and series.
#'
#' @param diagnostics An object of class \code{ms_diagnostics}
#' @param state Integer state index
#' @param series Integer series index
#'
#' @return A data.frame with columns: iteration, mean_sigma, sd_sigma, min_sigma, max_sigma
#'
#' @examples
#' \dontrun{
#' sigma_stats <- extract_sigma_stats(diag, state = 1, series = 1)
#' plot(sigma_stats$iteration, sigma_stats$mean_sigma, type = "b")
#' }
#'
#' @export
extract_sigma_stats <- function(diagnostics, state, series) {
  if (!inherits(diagnostics, "ms_diagnostics")) {
    stop("diagnostics must be an object of class 'ms_diagnostics'")
  }
  
  key <- paste0("state_", state, "_series_", series)
  sigma_data <- diagnostics$sigma_evolution[[key]]
  
  if (is.null(sigma_data) || length(sigma_data) == 0) {
    warning("No sigma evolution data found for state ", state, " series ", series)
    return(NULL)
  }
  
  data.frame(
    iteration = sapply(sigma_data, function(x) x$iteration),
    mean_sigma = sapply(sigma_data, function(x) x$mean_sigma),
    sd_sigma = sapply(sigma_data, function(x) x$sd_sigma),
    min_sigma = sapply(sigma_data, function(x) x$min_sigma),
    max_sigma = sapply(sigma_data, function(x) x$max_sigma)
  )
}


## === DIAGNOSTIC CHECK FUNCTIONS ===

#' Check EM Monotonicity
#'
#' Verify that the EM algorithm exhibits (near) monotonic log-likelihood improvement.
#'
#' @param diagnostics An object of class \code{ms_diagnostics}
#' @param tolerance Numeric threshold for acceptable decreases (default: 1e-4)
#'
#' @return List with \code{passed} (logical), \code{n_violations}, \code{violation_iters}
#'
#' @examples
#' \dontrun{
#' monotonicity_check <- check_em_monotonicity(diag)
#' if (!monotonicity_check$passed) {
#'   cat("Violations at iterations:", monotonicity_check$violation_iters, "\n")
#' }
#' }
#'
#' @export
check_em_monotonicity <- function(diagnostics, tolerance = 1e-4) {
  if (!inherits(diagnostics, "ms_diagnostics")) {
    stop("diagnostics must be an object of class 'ms_diagnostics'")
  }
  
  ll_changes <- extract_ll_trajectory(diagnostics, type = "change")
  
  violations <- which(ll_changes < -tolerance)
  n_violations <- length(violations)
  
  list(
    passed = n_violations <= 2,  # Allow up to 2 rare numerical issues
    n_violations = n_violations,
    violation_iters = violations,
    max_violation = if (n_violations > 0) min(ll_changes[violations]) else 0
  )
}


#' Check Parameter Stationarity
#'
#' Verify that GARCH and DCC parameters satisfy stationarity constraints.
#'
#' @param diagnostics An object of class \code{ms_diagnostics}
#' @param state Integer state index (if NULL, checks all states)
#'
#' @return List with \code{passed} (logical) and detailed violation information
#'
#' @examples
#' \dontrun{
#' stationarity_check <- check_param_stationarity(diag, state = 1)
#' if (!stationarity_check$passed) {
#'   print(stationarity_check$violations)
#' }
#' }
#'
#' @export
check_param_stationarity <- function(diagnostics, state = NULL) {
  if (!inherits(diagnostics, "ms_diagnostics")) {
    stop("diagnostics must be an object of class 'ms_diagnostics'")
  }
  
  states_to_check <- if (is.null(state)) {
    # Extract all state indices
    state_keys <- names(diagnostics$parameter_evolution)
    as.integer(gsub("state_", "", state_keys))
  } else {
    state
  }
  
  violations <- list()
  
  for (s in states_to_check) {
    state_key <- paste0("state_", s)
    state_data <- diagnostics$parameter_evolution[[state_key]]
    
    if (is.null(state_data) || length(state_data) == 0) next
    
    # Check final iteration parameters
    final_params <- state_data[[length(state_data)]]$parameters
    
    # Check GARCH parameters
    if (!is.null(final_params$garch_pars)) {
      for (i in seq_along(final_params$garch_pars)) {
        garch <- final_params$garch_pars[[i]]
        
        if (!is.null(garch$omega) && garch$omega <= 0) {
          violations[[length(violations) + 1]] <- list(
            state = s, series = i, parameter = "omega",
            value = garch$omega, constraint = "> 0"
          )
        }
        
        if (!is.null(garch$alpha1) && garch$alpha1 < 0) {
          violations[[length(violations) + 1]] <- list(
            state = s, series = i, parameter = "alpha1",
            value = garch$alpha1, constraint = ">= 0"
          )
        }
        
        if (!is.null(garch$beta1) && garch$beta1 < 0) {
          violations[[length(violations) + 1]] <- list(
            state = s, series = i, parameter = "beta1",
            value = garch$beta1, constraint = ">= 0"
          )
        }
        
        # Stationarity
        if (!is.null(garch$alpha1) && !is.null(garch$beta1)) {
          garch_sum <- garch$alpha1 + garch$beta1
          if (garch_sum >= 1) {
            violations[[length(violations) + 1]] <- list(
              state = s, series = i, parameter = "alpha1 + beta1",
              value = garch_sum, constraint = "< 1"
            )
          }
        }
      }
    }
    
    # Check DCC parameters (if dynamic)
    if (!is.null(final_params$correlation_type) && 
        final_params$correlation_type == "dynamic") {
      
      if (!is.null(final_params$alpha_1) && final_params$alpha_1 < 0) {
        violations[[length(violations) + 1]] <- list(
          state = s, parameter = "alpha_1",
          value = final_params$alpha_1, constraint = ">= 0"
        )
      }
      
      if (!is.null(final_params$beta_1) && final_params$beta_1 < 0) {
        violations[[length(violations) + 1]] <- list(
          state = s, parameter = "beta_1",
          value = final_params$beta_1, constraint = ">= 0"
        )
      }
      
      # DCC stationarity
      if (!is.null(final_params$alpha_1) && !is.null(final_params$beta_1)) {
        dcc_sum <- final_params$alpha_1 + final_params$beta_1
        if (dcc_sum >= 1) {
          violations[[length(violations) + 1]] <- list(
            state = s, parameter = "alpha_1 + beta_1",
            value = dcc_sum, constraint = "< 1"
          )
        }
      }
    }
  }
  
  list(
    passed = length(violations) == 0,
    n_violations = length(violations),
    violations = violations
  )
}


#' Check Convergence Achievement
#'
#' Verify that the model achieved convergence based on tolerance criteria.
#'
#' @param diagnostics An object of class \code{ms_diagnostics}
#' @param tolerance Numeric convergence tolerance
#'
#' @return List with \code{converged} (logical), \code{final_change}, \code{n_iterations}
#'
#' @export
check_convergence <- function(diagnostics, tolerance = 1e-4) {
  if (!inherits(diagnostics, "ms_diagnostics")) {
    stop("diagnostics must be an object of class 'ms_diagnostics'")
  }
  
  if (length(diagnostics$em_iterations) == 0) {
    return(list(converged = FALSE, final_change = NA, n_iterations = 0))
  }
  
  ll_changes <- extract_ll_trajectory(diagnostics, type = "change")
  final_change <- tail(ll_changes, 1)
  n_iters <- length(ll_changes)
  
  # Check if explicitly marked as converged
  explicitly_converged <- any(sapply(diagnostics$em_iterations, 
                                     function(x) x$converged %||% FALSE))
  
  # Or if change is below tolerance
  converged_by_tol <- abs(final_change) < tolerance
  
  list(
    converged = explicitly_converged || converged_by_tol,
    final_change = final_change,
    n_iterations = n_iters,
    explicitly_marked = explicitly_converged
  )
}


#' Check for Sigma Update Failures
#'
#' Identify cases where volatility (sigma) failed to update during estimation.
#'
#' @param diagnostics An object of class \code{ms_diagnostics}
#'
#' @return List with \code{has_failures} (logical) and details of failed updates
#'
#' @export
check_sigma_updates <- function(diagnostics) {
  if (!inherits(diagnostics, "ms_diagnostics")) {
    stop("diagnostics must be an object of class 'ms_diagnostics'")
  }
  
  failures <- list()
  
  for (key in names(diagnostics$sigma_evolution)) {
    series_data <- diagnostics$sigma_evolution[[key]]
    
    failed_iters <- which(sapply(series_data, function(x) !x$changed))
    
    if (length(failed_iters) > 0) {
      failures[[key]] <- list(
        location = key,
        n_failures = length(failed_iters),
        iterations = failed_iters
      )
    }
  }
  
  list(
    has_failures = length(failures) > 0,
    n_locations_affected = length(failures),
    details = failures
  )
}


## === ANALYSIS FUNCTIONS ===

#' Analyze Iteration in Detail
#'
#' Provide comprehensive analysis of a specific EM iteration.
#'
#' @param diagnostics An object of class \code{ms_diagnostics}
#' @param iteration Integer iteration number to analyze
#'
#' @return List with detailed iteration information
#'
#' @export
analyze_iteration <- function(diagnostics, iteration) {
  if (!inherits(diagnostics, "ms_diagnostics")) {
    stop("diagnostics must be an object of class 'ms_diagnostics'")
  }
  
  if (iteration > length(diagnostics$em_iterations)) {
    stop("Iteration ", iteration, " not found. Only ", 
         length(diagnostics$em_iterations), " iterations recorded.")
  }
  
  iter_data <- diagnostics$em_iterations[[iteration]]
  
  # Find warnings for this iteration
  iter_warnings <- Filter(function(w) w$iteration == iteration, 
                          diagnostics$warnings)
  
  # Find boundary events for this iteration
  iter_boundaries <- Filter(function(b) b$iteration == iteration,
                            diagnostics$boundary_events)
  
  # Get parameter snapshots for each state
  state_params <- list()
  for (state_key in names(diagnostics$parameter_evolution)) {
    state_data <- diagnostics$parameter_evolution[[state_key]]
    iter_match <- which(sapply(state_data, function(x) x$iteration == iteration))
    if (length(iter_match) > 0) {
      state_params[[state_key]] <- state_data[[iter_match]]$parameters
    }
  }
  
  list(
    iteration = iteration,
    log_likelihood = list(
      before_mstep = iter_data$log_lik_before_mstep,
      after_mstep = iter_data$log_lik_after_mstep,
      change = iter_data$ll_change,
      decreased = iter_data$ll_decreased
    ),
    timing = list(
      duration_seconds = iter_data$duration_seconds,
      timestamp = iter_data$timestamp
    ),
    warnings = iter_warnings,
    boundary_events = iter_boundaries,
    state_parameters = state_params,
    converged = iter_data$converged %||% FALSE
  )
}


#' Identify Problematic States
#'
#' Identify states with estimation problems based on diagnostic information.
#'
#' @param diagnostics An object of class \code{ms_diagnostics}
#'
#' @return List with problem states and issue descriptions
#'
#' @export
identify_problematic_states <- function(diagnostics) {
  if (!inherits(diagnostics, "ms_diagnostics")) {
    stop("diagnostics must be an object of class 'ms_diagnostics'")
  }
  
  problems <- list()
  
  # Check for states with constant correlation when it might be unexpected
  for (state_key in names(diagnostics$parameter_evolution)) {
    state_num <- as.integer(gsub("state_", "", state_key))
    state_data <- diagnostics$parameter_evolution[[state_key]]
    
    if (length(state_data) == 0) next
    
    final_params <- state_data[[length(state_data)]]$parameters
    
    # Constant correlation
    if (!is.null(final_params$correlation_type) && 
        final_params$correlation_type == "constant") {
      
      # Check if there were boundary events
      boundary_events <- Filter(
        function(b) b$state == state_num && grepl("alpha", b$parameter),
        diagnostics$boundary_events
      )
      
      if (length(boundary_events) > 0) {
        problems[[state_key]] <- c(
          problems[[state_key]], 
          "DCC parameters at boundary - switched to constant correlation"
        )
      }
    }
    
    # Check for parameter instability
    if (length(state_data) >= 5) {
      last_5 <- tail(state_data, 5)
      
      # Check alpha_1 variability
      if (!is.null(final_params$alpha_1)) {
        alpha_vals <- sapply(last_5, function(x) x$parameters$alpha_1 %||% NA)
        if (any(!is.na(alpha_vals))) {
          alpha_range <- max(alpha_vals, na.rm = TRUE) - min(alpha_vals, na.rm = TRUE)
          if (alpha_range > 0.05) {
            problems[[state_key]] <- c(
              problems[[state_key]],
              paste0("alpha_1 unstable (range: ", round(alpha_range, 4), ")")
            )
          }
        }
      }
    }
  }
  
  # Check for states with many warnings
  warning_by_state <- table(sapply(diagnostics$warnings, 
                                   function(w) w$details$state %||% NA))
  warning_by_state <- warning_by_state[!is.na(names(warning_by_state))]
  
  for (s in names(warning_by_state)) {
    if (warning_by_state[s] > 5) {
      state_key <- paste0("state_", s)
      problems[[state_key]] <- c(
        problems[[state_key]],
        paste0("Excessive warnings (", warning_by_state[s], ")")
      )
    }
  }
  
  list(
    has_problems = length(problems) > 0,
    n_states_affected = length(problems),
    problems = problems
  )
}


#' Compare States
#'
#' Compare parameter estimates across states to assess regime differentiation.
#'
#' @param diagnostics An object of class \code{ms_diagnostics}
#' @param states Integer vector of states to compare (if NULL, compares all)
#'
#' @return List with comparison metrics
#'
#' @export
compare_states <- function(diagnostics, states = NULL) {
  if (!inherits(diagnostics, "ms_diagnostics")) {
    stop("diagnostics must be an object of class 'ms_diagnostics'")
  }
  
  state_keys <- names(diagnostics$parameter_evolution)
  
  if (!is.null(states)) {
    state_keys <- paste0("state_", states)
    state_keys <- state_keys[state_keys %in% names(diagnostics$parameter_evolution)]
  }
  
  if (length(state_keys) < 2) {
    stop("Need at least 2 states to compare")
  }
  
  # Extract final parameters for each state
  final_params <- lapply(state_keys, function(key) {
    state_data <- diagnostics$parameter_evolution[[key]]
    state_data[[length(state_data)]]$parameters
  })
  names(final_params) <- state_keys
  
  comparisons <- list()
  
  # Compare GARCH volatility (omega for first series)
  omega_vals <- sapply(final_params, function(p) {
    if (!is.null(p$garch_pars) && length(p$garch_pars) > 0) {
      p$garch_pars[[1]]$omega %||% NA
    } else {
      NA
    }
  })
  
  if (any(!is.na(omega_vals))) {
    comparisons$omega_range <- max(omega_vals, na.rm = TRUE) - 
      min(omega_vals, na.rm = TRUE)
    comparisons$omega_relative_diff <- comparisons$omega_range / 
      mean(omega_vals, na.rm = TRUE)
  }
  
  # Compare DCC dynamics
  alpha_vals <- sapply(final_params, function(p) p$alpha_1 %||% NA)
  
  if (any(!is.na(alpha_vals))) {
    comparisons$alpha_range <- max(alpha_vals, na.rm = TRUE) - 
      min(alpha_vals, na.rm = TRUE)
  }
  
  # Check correlation types
  corr_types <- sapply(final_params, function(p) {
    p$correlation_type %||% "unknown"
  })
  comparisons$correlation_types <- table(corr_types)
  
  # Overall assessment
  comparisons$states_similar <- (
    (!is.null(comparisons$omega_relative_diff) && 
       comparisons$omega_relative_diff < 0.2) ||
      (!is.null(comparisons$alpha_range) && 
         comparisons$alpha_range < 0.05)
  )
  
  list(
    n_states_compared = length(state_keys),
    comparisons = comparisons,
    final_parameters = final_params
  )
}


## === EXPORT FUNCTIONS ===

#' Export Diagnostics to CSV Files
#'
#' Export diagnostic data to multiple CSV files for external analysis.
#'
#' @param diagnostics An object of class \code{ms_diagnostics}
#' @param prefix Character string prefix for output files
#' @param path Directory path for output files (default: current directory)
#'
#' @return Character vector of created file paths
#'
#' @export
export_diagnostics_csv <- function(diagnostics, prefix = "diagnostic", path = ".") {
  if (!inherits(diagnostics, "ms_diagnostics")) {
    stop("diagnostics must be an object of class 'ms_diagnostics'")
  }
  
  created_files <- character(0)
  
  # EM iterations
  if (length(diagnostics$em_iterations) > 0) {
    em_df <- do.call(rbind, lapply(diagnostics$em_iterations, function(x) {
      data.frame(
        iteration = x$iteration,
        ll_before = x$log_lik_before_mstep,
        ll_after = x$log_lik_after_mstep,
        ll_change = x$ll_change,
        ll_decreased = x$ll_decreased,
        duration_sec = x$duration_seconds,
        converged = x$converged %||% FALSE,
        stringsAsFactors = FALSE
      )
    }))
    
    filepath <- file.path(path, paste0(prefix, "_em_iterations.csv"))
    write.csv(em_df, filepath, row.names = FALSE)
    created_files <- c(created_files, filepath)
  }
  
  # Boundary events
  if (length(diagnostics$boundary_events) > 0) {
    boundary_df <- do.call(rbind, lapply(diagnostics$boundary_events, function(x) {
      data.frame(
        iteration = x$iteration,
        state = x$state,
        parameter = x$parameter,
        value = x$value,
        boundary_type = x$boundary_type,
        action_taken = x$action_taken,
        stringsAsFactors = FALSE
      )
    }))
    
    filepath <- file.path(path, paste0(prefix, "_boundary_events.csv"))
    write.csv(boundary_df, filepath, row.names = FALSE)
    created_files <- c(created_files, filepath)
  }
  
  # Warnings
  if (length(diagnostics$warnings) > 0) {
    warning_df <- do.call(rbind, lapply(diagnostics$warnings, function(x) {
      data.frame(
        iteration = x$iteration,
        type = x$type,
        message = x$message,
        stringsAsFactors = FALSE
      )
    }))
    
    filepath <- file.path(path, paste0(prefix, "_warnings.csv"))
    write.csv(warning_df, filepath, row.names = FALSE)
    created_files <- c(created_files, filepath)
  }
  
  invisible(created_files)
}


## === UTILITY FUNCTIONS ===

#' Get Diagnostic Summary Statistics
#'
#' Compute summary statistics from diagnostic object.
#'
#' @param diagnostics An object of class \code{ms_diagnostics}
#'
#' @return List with summary statistics
#'
#' @export
diagnostic_summary_stats <- function(diagnostics) {
  if (!inherits(diagnostics, "ms_diagnostics")) {
    stop("diagnostics must be an object of class 'ms_diagnostics'")
  }
  
  ll_changes <- extract_ll_trajectory(diagnostics, type = "change")
  
  list(
    n_iterations = length(diagnostics$em_iterations),
    n_states = length(diagnostics$parameter_evolution),
    n_warnings = length(diagnostics$warnings),
    n_boundary_events = length(diagnostics$boundary_events),
    ll_summary = if (length(ll_changes) > 0) {
      list(
        total_improvement = sum(ll_changes),
        mean_change = mean(ll_changes),
        min_change = min(ll_changes),
        max_change = max(ll_changes),
        n_decreases = sum(ll_changes < -1e-6)
      )
    } else {
      NULL
    },
    total_time_seconds = sum(sapply(diagnostics$em_iterations, 
                                    function(x) x$duration_seconds)),
    converged = check_convergence(diagnostics)$converged
  )
}


## === SIMULATION FUNCTIONS FOR TESTING ===

#' Simulate DCC-GARCH Data with Higher-Order Support
#'
#' Simulate realistic multivariate time series data from a DCC(q,p)-GARCH process.
#' Supports arbitrary DCC orders through vector-valued dcc_alpha and dcc_beta.
#' This function is primarily intended for testing and examples.
#'
#' @param n Integer number of observations to simulate
#' @param k Integer number of series (default: 2)
#' @param omega Numeric vector of length k: GARCH intercepts (default: c(0.05, 0.08))
#' @param alpha_garch Numeric vector of length k: ARCH effects (default: c(0.10, 0.12))
#' @param beta_garch Numeric vector of length k: GARCH effects (default: c(0.85, 0.82))
#' @param dcc_alpha Numeric scalar or vector: DCC alpha parameters. 
#'   Length determines q (number of alpha lags). (default: 0.04)
#' @param dcc_beta Numeric scalar or vector: DCC beta parameters.
#'   Length determines p (number of beta lags). (default: 0.93)
#' @param Qbar Matrix (k x k): Unconditional correlation matrix. If NULL, uses
#'   moderate correlation (0.5 off-diagonal)
#' @param seed Integer random seed for reproducibility
#'
#' @return A matrix of dimension (n x k) with simulated returns
#'
#' @details
#' The function simulates data from the DCC(q,p)-GARCH model:
#' \deqn{y_{i,t} = \sqrt{h_{i,t}} z_{i,t}}
#' \deqn{h_{i,t} = \omega_i + \alpha_i y_{i,t-1}^2 + \beta_i h_{i,t-1}}
#' \deqn{Q_t = \bar{Q}(1-\sum\alpha-\sum\beta) + \sum_{j=1}^{q} \alpha_j z_{t-j}z_{t-j}' + \sum_{j=1}^{p} \beta_j Q_{t-j}}
#' \deqn{R_t = \text{diag}(Q_t)^{-1/2} Q_t \text{diag}(Q_t)^{-1/2}}
#'
#' where \eqn{z_t \sim N(0, R_t)}.
#'
#' For backward compatibility, scalar dcc_alpha and dcc_beta produce DCC(1,1).
#' To simulate DCC(q,p), provide dcc_alpha as vector of length q and dcc_beta 
#' as vector of length p.
#'
#' @examples
#' \dontrun
#' # Simulate 500 observations with default DCC(1,1) parameters
#' y <- simulate_dcc_garch(n = 500, seed = 42)
#' 
#' # DCC(1,2): one alpha lag, two beta lags
#' y_dcc12 <- simulate_dcc_garch(
#'   n = 300,
#'   dcc_alpha = 0.05,
#'   dcc_beta = c(0.50, 0.40),
#'   seed = 123
#' )
#' 
#' # DCC(2,2): two lags each
#' y_dcc22 <- simulate_dcc_garch(
#'   n = 300,
#'   dcc_alpha = c(0.03, 0.02),
#'   dcc_beta = c(0.50, 0.40),
#'   seed = 456
#' )
#' 
#' # Constant correlation (set dcc_alpha = 0)
#' y_const <- simulate_dcc_garch(n = 200, dcc_alpha = 0, dcc_beta = 0, seed = 789)
#' }
#'
#' @export
simulate_dcc_garch <- function(n, 
                               k = 2,
                               omega = c(0.05, 0.08),
                               alpha_garch = c(0.10, 0.12),
                               beta_garch = c(0.85, 0.82),
                               dcc_alpha = 0.04,
                               dcc_beta = 0.93,
                               Qbar = NULL,
                               seed = NULL) {
  
  ## Input validation
  if (n <= 0 || n != round(n)) {
    stop("n must be a positive integer")
  }
  if (k <= 0 || k != round(k)) {
    stop("k must be a positive integer")
  }
  if (length(omega) != k || length(alpha_garch) != k || length(beta_garch) != k) {
    stop("omega, alpha_garch, and beta_garch must have length k")
  }
  if (any(omega <= 0)) {
    stop("All omega values must be positive")
  }
  if (any(alpha_garch < 0) || any(beta_garch < 0)) {
    stop("alpha_garch and beta_garch must be non-negative")
  }
  if (any((alpha_garch + beta_garch) >= 1)) {
    stop("GARCH process must be stationary: alpha + beta < 1")
  }
  
  ## Convert scalars to vectors for uniform handling (supports higher-order DCC)
  dcc_alpha <- as.numeric(dcc_alpha)
  dcc_beta <- as.numeric(dcc_beta)
  
  q_order <- length(dcc_alpha)  ## Number of alpha lags (ARCH-like)
  p_order <- length(dcc_beta)   ## Number of beta lags (GARCH-like)
  
  if (any(dcc_alpha < 0) || any(dcc_beta < 0)) {
    stop("All dcc_alpha and dcc_beta values must be non-negative")
  }
  
  ## Check DCC stationarity: sum(alpha) + sum(beta) < 1
  dcc_persistence <- sum(dcc_alpha) + sum(dcc_beta)
  has_dcc_dynamics <- dcc_persistence > 0
  
  if (has_dcc_dynamics && dcc_persistence >= 1) {
    stop(sprintf(
      "DCC process must be stationary: sum(alpha) + sum(beta) < 1, got %.4f",
      dcc_persistence
    ))
  }
  
  ## Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  ## Default Qbar (unconditional correlation)
  if (is.null(Qbar)) {
    Qbar <- matrix(0.5, k, k)
    diag(Qbar) <- 1
  }
  
  ## Validate Qbar
  if (!is.matrix(Qbar) || nrow(Qbar) != k || ncol(Qbar) != k) {
    stop("Qbar must be a k x k matrix")
  }
  if (!all(diag(Qbar) == 1)) {
    stop("Qbar must have 1s on the diagonal")
  }
  eig_vals <- eigen(Qbar, symmetric = TRUE, only.values = TRUE)$values
  if (any(eig_vals < 0)) {
    stop("Qbar must be positive definite")
  }
  
  ## Determine maximum lag needed
  max_lag <- max(p_order, q_order)
  
  ## Initialize storage
  y_sim <- matrix(0, n, k)
  h <- matrix(0, n, k)
  z_history <- vector("list", max_lag)  ## Store lagged z vectors
  Q_history <- vector("list", max_lag)  ## Store lagged Q matrices
  
  ## Initialize all history with Qbar and zero z
  for (j in seq_len(max_lag)) {
    z_history[[j]] <- rep(0, k)
    Q_history[[j]] <- Qbar
  }
  
  ## Current Q and R
  Q <- Qbar
  R <- Qbar
  
  ## Initialize conditional variances at unconditional level
  for (i in 1:k) {
    h[1, i] <- omega[i] / (1 - alpha_garch[i] - beta_garch[i])
  }
  
  ## Intercept weight for DCC recursion
  intercept_weight <- 1 - dcc_persistence
  
  ## Simulate
  for (t in 1:n) {
    ## Draw correlated innovations from current R
    z <- as.vector(mvtnorm::rmvnorm(1, sigma = R))
    
    ## Generate returns
    for (i in 1:k) {
      y_sim[t, i] <- sqrt(h[t, i]) * z[i]
      
      ## Update variance for next period
      if (t < n) {
        h[t+1, i] <- omega[i] + alpha_garch[i] * y_sim[t, i]^2 + 
          beta_garch[i] * h[t, i]
      }
    }
    
    ## Update DCC dynamics for next period
    if (t < n && has_dcc_dynamics) {
      
      ## Shift history: move everything back one slot
      ## Most recent goes to position 1
      if (max_lag > 1) {
        for (j in max_lag:2) {
          z_history[[j]] <- z_history[[j-1]]
          Q_history[[j]] <- Q_history[[j-1]]
        }
      }
      z_history[[1]] <- z
      Q_history[[1]] <- Q
      
      ## Compute new Q using higher-order recursion:
      ## Q_t = Qbar * (1 - sum(alpha) - sum(beta)) 
      ##     + sum_{j=1}^{q} alpha_j * z_{t-j} z_{t-j}'
      ##     + sum_{j=1}^{p} beta_j * Q_{t-j}
      
      Q_new <- Qbar * intercept_weight
      
      ## Add alpha terms (lagged outer products of z)
      for (j in seq_len(q_order)) {
        z_lag <- z_history[[j]]
        Q_new <- Q_new + dcc_alpha[j] * (z_lag %*% t(z_lag))
      }
      
      ## Add beta terms (lagged Q matrices)
      for (j in seq_len(p_order)) {
        Q_new <- Q_new + dcc_beta[j] * Q_history[[j]]
      }
      
      Q <- Q_new
      
      ## Standardize to correlation matrix
      Q_diag <- diag(Q)
      if (any(Q_diag <= 0)) {
        ## Numerical issue - regularize
        Q <- Q + diag(1e-6, k)
        Q_diag <- diag(Q)
      }
      Q_diag_inv_sqrt <- diag(1/sqrt(Q_diag))
      R <- Q_diag_inv_sqrt %*% Q %*% Q_diag_inv_sqrt
      
      ## Ensure R stays valid correlation matrix
      diag(R) <- 1
    }
  }
  
  colnames(y_sim) <- paste0("series_", 1:k)
  
  return(y_sim)
}

#' #' Simulate DCC-GARCH Data
#' #'
#' #' Simulate realistic multivariate time series data from a DCC-GARCH process.
#' #' This function is primarily intended for testing and examples.
#' #'
#' #' @param n Integer number of observations to simulate
#' #' @param k Integer number of series (default: 2)
#' #' @param omega Numeric vector of length k: GARCH intercepts (default: c(0.05, 0.08))
#' #' @param alpha_garch Numeric vector of length k: ARCH effects (default: c(0.10, 0.12))
#' #' @param beta_garch Numeric vector of length k: GARCH effects (default: c(0.85, 0.82))
#' #' @param dcc_alpha Numeric: DCC alpha parameter (default: 0.04)
#' #' @param dcc_beta Numeric: DCC beta parameter (default: 0.93)
#' #' @param Qbar Matrix (k x k): Unconditional correlation matrix. If NULL, uses
#' #'   moderate correlation (0.5 off-diagonal)
#' #' @param seed Integer random seed for reproducibility
#' #'
#' #' @return A matrix of dimension (n x k) with simulated returns
#' #'
#' #' @details
#' #' The function simulates data from the DCC-GARCH model:
#' #' \deqn{y_{i,t} = \sqrt{h_{i,t}} z_{i,t}}
#' #' \deqn{h_{i,t} = \omega_i + \alpha_i y_{i,t-1}^2 + \beta_i h_{i,t-1}}
#' #' \deqn{Q_t = \bar{Q}(1-\alpha-\beta) + \alpha z_{t-1}z_{t-1}' + \beta Q_{t-1}}
#' #' \deqn{R_t = \text{diag}(Q_t)^{-1/2} Q_t \text{diag}(Q_t)^{-1/2}}
#' #'
#' #' where \eqn{z_t \sim N(0, R_t)}.
#' #'
#' #' @examples
#' #' \dontrun{
#' #' # Simulate 500 observations with default parameters
#' #' y <- simulate_dcc_garch(n = 500, seed = 42)
#' #' 
#' #' # Custom parameters with higher volatility
#' #' y_high_vol <- simulate_dcc_garch(
#' #'   n = 300,
#' #'   omega = c(0.15, 0.20),
#' #'   alpha_garch = c(0.15, 0.18),
#' #'   beta_garch = c(0.70, 0.65),
#' #'   seed = 123
#' #' )
#' #' 
#' #' # Constant correlation (set dcc_alpha = 0)
#' #' y_const <- simulate_dcc_garch(n = 200, dcc_alpha = 0, dcc_beta = 0, seed = 456)
#' #' }
#' #'
#' #' @export
#' simulate_dcc_garch <- function(n, 
#'                                k = 2,
#'                                omega = c(0.05, 0.08),
#'                                alpha_garch = c(0.10, 0.12),
#'                                beta_garch = c(0.85, 0.82),
#'                                dcc_alpha = 0.04,
#'                                dcc_beta = 0.93,
#'                                Qbar = NULL,
#'                                seed = NULL) {
#'   
#'   # Input validation
#'   if (n <= 0 || n != round(n)) {
#'     stop("n must be a positive integer")
#'   }
#'   if (k <= 0 || k != round(k)) {
#'     stop("k must be a positive integer")
#'   }
#'   if (length(omega) != k || length(alpha_garch) != k || length(beta_garch) != k) {
#'     stop("omega, alpha_garch, and beta_garch must have length k")
#'   }
#'   if (any(omega <= 0)) {
#'     stop("All omega values must be positive")
#'   }
#'   if (any(alpha_garch < 0) || any(beta_garch < 0)) {
#'     stop("alpha_garch and beta_garch must be non-negative")
#'   }
#'   if (any((alpha_garch + beta_garch) >= 1)) {
#'     stop("GARCH process must be stationary: alpha + beta < 1")
#'   }
#'   if (dcc_alpha < 0 || dcc_beta < 0) {
#'     stop("dcc_alpha and dcc_beta must be non-negative")
#'   }
#'   if ((dcc_alpha + dcc_beta) >= 1 && (dcc_alpha > 0 || dcc_beta > 0)) {
#'     stop("DCC process must be stationary: dcc_alpha + dcc_beta < 1")
#'   }
#'   
#'   # Set seed if provided
#'   if (!is.null(seed)) {
#'     set.seed(seed)
#'   }
#'   
#'   # Default Qbar (unconditional correlation)
#'   if (is.null(Qbar)) {
#'     Qbar <- matrix(0.5, k, k)
#'     diag(Qbar) <- 1
#'   }
#'   
#'   # Validate Qbar
#'   if (!is.matrix(Qbar) || nrow(Qbar) != k || ncol(Qbar) != k) {
#'     stop("Qbar must be a k x k matrix")
#'   }
#'   if (!all(diag(Qbar) == 1)) {
#'     stop("Qbar must have 1s on the diagonal")
#'   }
#'   eig_vals <- eigen(Qbar, symmetric = TRUE, only.values = TRUE)$values
#'   if (any(eig_vals < 0)) {
#'     stop("Qbar must be positive definite")
#'   }
#'   
#'   # Initialize
#'   y_sim <- matrix(0, n, k)
#'   h <- matrix(0, n, k)
#'   Q <- Qbar
#'   R <- Qbar
#'   
#'   # Initialize conditional variances
#'   for (i in 1:k) {
#'     h[1, i] <- omega[i] / (1 - alpha_garch[i] - beta_garch[i])
#'   }
#'   
#'   # Simulate
#'   for (t in 1:n) {
#'     # Draw correlated innovations
#'     z <- as.vector(mvtnorm::rmvnorm(1, sigma = R))
#'     
#'     # Generate returns
#'     for (i in 1:k) {
#'       y_sim[t, i] <- sqrt(h[t, i]) * z[i]
#'       
#'       # Update variance for next period
#'       if (t < n) {
#'         h[t+1, i] <- omega[i] + alpha_garch[i] * y_sim[t, i]^2 + 
#'           beta_garch[i] * h[t, i]
#'       }
#'     }
#'     
#'     # Update DCC dynamics for next period
#'     if (t < n && (dcc_alpha > 0 || dcc_beta > 0)) {
#'       z_mat <- matrix(z, ncol = 1)
#'       Q <- Qbar * (1 - dcc_alpha - dcc_beta) + 
#'         dcc_alpha * (z_mat %*% t(z_mat)) + 
#'         dcc_beta * Q
#'       
#'       # Standardize to correlation
#'       Q_diag_inv_sqrt <- diag(1/sqrt(diag(Q)))
#'       R <- Q_diag_inv_sqrt %*% Q %*% Q_diag_inv_sqrt
#'     }
#'   }
#'   
#'   colnames(y_sim) <- paste0("series_", 1:k)
#'   
#'   return(y_sim)
#' }


#' Generate Model Specification for DCC-GARCH
#'
#' Generate a properly structured specification list for MS-VARMA-GARCH estimation
#' with DCC dynamics. This function is primarily intended for testing and examples.
#'
#' @param M Integer number of states
#' @param k Integer number of series (default: 2)
#' @param var_order Integer VAR order (default: 1)
#' @param garch_order Integer vector of length 2: GARCH(p,q) order (default: c(1,1))
#' @param distribution Character: distribution for DCC ("mvn" or "mvt")
#' @param seed Integer random seed for generating starting values
#' @param simple Logical: if TRUE, use identical starting values across states
#'   (faster but may not reflect regime differences). Default: FALSE
#'
#' @return A list of length M containing properly formatted specifications
#'
#' @details
#' The function generates starting parameter values with mild regime differentiation
#' unless \code{simple = TRUE}. States are ordered from low to high volatility.
#'
#' For state j:
#' \itemize{
#'   \item VAR parameters: Small values near 0.1
#'   \item GARCH omega: Increasing across states (0.05 to 0.15)
#'   \item GARCH alpha: Increasing across states (0.08 to 0.15)
#'   \item GARCH beta: Decreasing across states (0.85 to 0.75)
#'   \item DCC alpha: Increasing across states (0.03 to 0.10)
#'   \item DCC beta: Decreasing across states (0.94 to 0.85)
#' }
#'
#' @examples
#' \dontrun{
#' # Generate 2-state specification with defaults
#' spec <- generate_dcc_spec(M = 2)
#' 
#' # 3-state specification with Student-t distribution
#' spec_mvt <- generate_dcc_spec(M = 3, distribution = "mvt")
#' 
#' # Simple specification (identical starting values)
#' spec_simple <- generate_dcc_spec(M = 2, simple = TRUE)
#' 
#' # Use with simulated data
#' y <- simulate_dcc_garch(n = 300, seed = 42)
#' fit <- fit_ms_varma_garch(
#'   y = y, 
#'   M = 2, 
#'   spec = spec,
#'   model_type = "multivariate",
#'   control = list(max_iter = 20, tol = 1e-4)
#' )
#' }
#'
#' @export
generate_dcc_spec <- function(M, 
                              k = 2,
                              var_order = 1,
                              garch_order = c(1, 1),
                              distribution = c("mvn", "mvt"),
                              seed = NULL,
                              simple = FALSE) {
  
  # Input validation
  if (M < 2 || M != round(M)) {
    stop("M must be an integer >= 2")
  }
  if (k < 2 || k != round(k)) {
    stop("k must be an integer >= 2")
  }
  if (var_order < 0 || var_order != round(var_order)) {
    stop("var_order must be a non-negative integer")
  }
  if (length(garch_order) != 2 || any(garch_order < 0)) {
    stop("garch_order must be a length-2 vector of non-negative integers")
  }
  
  distribution <- match.arg(distribution)
  
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Base univariate GARCH specification (for tsmarch)
  spec_uni_garch <- list(
    model = "garch",
    garch_order = garch_order,
    distribution = "norm"
  )
  
  # DCC specification arguments (common across states)
  dcc_spec_args <- list(
    dcc_order = c(1, 1),
    dynamics = "dcc",
    distribution = distribution,
    garch_model = list(
      univariate = replicate(k, spec_uni_garch, simplify = FALSE)
    )
  )
  
  # Number of VAR parameters: intercept + k lags for each of k series
  n_var_pars <- k * (1 + k * var_order)
  
  # Generate state-specific specifications
  spec_list <- vector("list", M)
  
  for (j in 1:M) {
    if (simple) {
      # Simple mode: identical starting values
      omega_start <- rep(0.1, k)
      alpha_start <- rep(0.1, k)
      beta_start <- rep(0.8, k)
      dcc_alpha_start <- 0.05
      dcc_beta_start <- 0.90
    } else {
      # Differentiated starting values across states
      # State 1: Low volatility, low correlation dynamics
      # State M: High volatility, high correlation dynamics
      
      state_factor <- (j - 1) / max(M - 1, 1)
      
      # GARCH parameters: volatility increases with state
      omega_start <- 0.05 + state_factor * 0.10 + 
        rnorm(k, 0, 0.01) * seq_len(k) / k
      omega_start <- pmax(omega_start, 0.03)  # Ensure positive
      
      alpha_start <- 0.08 + state_factor * 0.07 + 
        rnorm(k, 0, 0.01) * seq_len(k) / k
      alpha_start <- pmax(alpha_start, 0.05)
      
      beta_start <- 0.85 - state_factor * 0.10 + 
        rnorm(k, 0, 0.01) * seq_len(k) / k
      beta_start <- pmin(beta_start, 0.88)
      
      # Ensure stationarity
      for (i in 1:k) {
        while (alpha_start[i] + beta_start[i] >= 0.98) {
          alpha_start[i] <- alpha_start[i] * 0.95
        }
      }
      
      # DCC parameters: dynamics increase with state
      dcc_alpha_start <- 0.03 + state_factor * 0.07
      dcc_beta_start <- 0.94 - state_factor * 0.09
      
      # Ensure DCC stationarity
      while (dcc_alpha_start + dcc_beta_start >= 0.98) {
        dcc_alpha_start <- dcc_alpha_start * 0.95
      }
    }
    
    # Build GARCH parameters for each series
    garch_pars_list <- vector("list", k)
    for (i in 1:k) {
      garch_pars_list[[i]] <- list(
        omega = omega_start[i],
        alpha1 = alpha_start[i],
        beta1 = beta_start[i]
      )
    }
    
    # Distribution parameters
    dist_pars <- if (distribution == "mvt") {
      list(shape = 8.0)
    } else {
      NULL
    }
    
    # Complete specification for this state
    spec_list[[j]] <- list(
      var_order = var_order,
      garch_spec_fun = "dcc_modelspec",
      garch_spec_args = dcc_spec_args,
      distribution = distribution,
      start_pars = list(
        var_pars = rep(0.1, n_var_pars),
        garch_pars = garch_pars_list,
        dcc_pars = list(alpha_1 = dcc_alpha_start, beta_1 = dcc_beta_start),
        dist_pars = dist_pars
      )
    )
  }
  
  return(spec_list)
}


#' @title Generate Single-State DCC Specification
#' @description Creates a specification for a single-state (no regime switching) 
#'   DCC-GARCH model. Useful for testing parameter recovery.
#' @keywords internal
generate_single_state_dcc_spec <- function(k = 2,
                                           var_order = 1,
                                           garch_order = c(1, 1),
                                           distribution = "mvn",
                                           omega = NULL,
                                           alpha_garch = NULL,
                                           beta_garch = NULL,
                                           dcc_alpha = NULL,
                                           dcc_beta = NULL,
                                           seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Default parameter values
  omega <- omega %||% rep(0.1, k)
  alpha_garch <- alpha_garch %||% rep(0.1, k)
  beta_garch <- beta_garch %||% rep(0.8, k)
  dcc_alpha <- dcc_alpha %||% 0.05
  dcc_beta <- dcc_beta %||% 0.90
  
  # Univariate GARCH spec
  spec_uni_garch <- list(
    model = "garch",
    garch_order = garch_order,
    distribution = "norm"
  )
  
  # DCC arguments
  dcc_spec_args <- list(
    dcc_order = c(1, 1),
    dynamics = "dcc",
    distribution = distribution,
    garch_model = list(
      univariate = replicate(k, spec_uni_garch, simplify = FALSE)
    )
  )
  
  # Build GARCH parameter list
  garch_pars_list <- lapply(1:k, function(i) {
    list(
      omega = omega[i],
      alpha1 = alpha_garch[i],
      beta1 = beta_garch[i]
    )
  })
  
  # Distribution parameters
  dist_pars <- if (distribution == "mvt") {
    list(shape = 8.0)
  } else {
    NULL
  }
  
  # Complete specification (wrapped in list for consistency with M-state case)
  spec_list <- list(
    list(
      var_order = var_order,
      garch_spec_fun = "dcc_modelspec",
      garch_spec_args = dcc_spec_args,
      distribution = distribution,
      start_pars = list(
        var_pars = rep(0.1, k * (1 + k * var_order)),
        garch_pars = garch_pars_list,
        dcc_pars = list(alpha_1 = dcc_alpha, beta_1 = dcc_beta),
        dist_pars = dist_pars
      )
    )
  )
  
  return(spec_list)
}