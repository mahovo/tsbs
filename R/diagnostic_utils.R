#' @title Diagnostic Utility Functions for MS-VARMA-GARCH
#' @description Convenience functions for extracting and analyzing diagnostic information
#' @name diagnostic_utils
NULL

#' Extract Parameter Trajectory
#'
#' Extract the evolution of parameters across EM iterations for a given state.
#' Supports DCC, CGARCH, and GOGARCH model types with automatic detection.
#'
#' @param diagnostics An object of class \code{ms_diagnostics}
#' @param state Integer state index
#' @param param_name Character string specifying parameter name. For single parameter
#'   extraction use names like "alpha_1", "beta_1", "gamma_1", "shape", "omega".
#'   Use "all" to extract all relevant correlation parameters based on model type.
#' @param series Optional integer for series-specific GARCH parameters (omega, alpha1, beta1)
#' @param component Optional integer for GOGARCH component-specific parameters
#' @param model_type Character: model type override. One of "auto" (default), "dcc", 
#'   "cgarch", or "gogarch". When "auto", attempts to detect from diagnostics.
#'
#' @return For single parameter: a data.frame with columns \code{iteration} and \code{value}.
#'   For "all": a data.frame with iteration and all relevant parameters as columns.
#'   Returns NULL if parameter not found.
#'
#' @details
#' The function extracts different parameters depending on model type:
#' \itemize{
#'   \item \strong{DCC}: alpha_1, beta_1 (correlation dynamics), shape (if MVT)
#'   \item \strong{CGARCH}: alpha_1, beta_1, gamma_1 (if ADCC), shape (if MVT copula)
#'   \item \strong{GOGARCH}: component GARCH parameters (omega, alpha1, beta1 per component)
#' }
#'
#' @examples
#' \dontrun{
#' # Extract single DCC parameter
#' alpha_traj <- extract_param_trajectory(diag, state = 1, param_name = "alpha_1")
#' plot(alpha_traj$iteration, alpha_traj$value, type = "b")
#' 
#' # Extract all correlation parameters (auto-detects model type)
#' all_traj <- extract_param_trajectory(diag, state = 1, param_name = "all")
#' 
#' # Extract GOGARCH component parameter
#' comp1_alpha <- extract_param_trajectory(diag, state = 1, param_name = "alpha1", 
#'                                          component = 1, model_type = "gogarch")
#' 
#' # Extract series-specific GARCH omega
#' omega_traj <- extract_param_trajectory(diag, state = 2, param_name = "omega", series = 1)
#' }
#'
#' @export
extract_param_trajectory <- function(diagnostics, state, param_name, 
                                     series = NULL, component = NULL,
                                     model_type = c("auto", "dcc", "cgarch", "gogarch")) {
  
  if (!inherits(diagnostics, "ms_diagnostics")) {
    stop("diagnostics must be an object of class 'ms_diagnostics'")
  }
  
  model_type <- match.arg(model_type)
  
  state_key <- paste0("state_", state)
  state_data <- diagnostics$parameter_evolution[[state_key]]
  
  if (is.null(state_data) || length(state_data) == 0) {
    warning("No parameter evolution data found for state ", state)
    return(NULL)
  }
  
  ## Auto-detect model type if needed
  if (model_type == "auto") {
    model_type <- .detect_model_type_from_params(state_data[[1]]$parameters)
  }
  
  ## Extract iterations
  iterations <- sapply(state_data, function(x) x$iteration)
  
  ## Handle "all" - extract all relevant correlation/model parameters
  if (param_name == "all") {
    return(.extract_all_params(state_data, iterations, model_type))
  }
  
  ## Single parameter extraction
  trajectory <- sapply(state_data, function(iter_data) {
    params <- iter_data$parameters
    .extract_single_param(params, param_name, series, component, model_type)
  })
  
  if (all(is.na(trajectory))) {
    warning("Parameter '", param_name, "' not found in state ", state)
    return(NULL)
  }
  
  data.frame(iteration = iterations, value = trajectory)
}


#' @keywords internal
.detect_model_type_from_params <- function(params) {
  ## GOGARCH: has ica_info or mixing_matrix
  if (!is.null(params$ica_info) || !is.null(params$mixing_matrix) ||
      isTRUE(params$model_type == "gogarch")) {
    return("gogarch")
  }
  
  ## CGARCH: has copula or transformation info
  if (!is.null(params$copula) || !is.null(params$transformation) ||
      !is.null(params$pit_method) || isTRUE(params$model_type == "cgarch")) {
    return("cgarch")
  }
  
  ## Default to DCC
  "dcc"
}


#' @keywords internal
.extract_single_param <- function(params, param_name, series, component, model_type) {
  
  ## Handle GOGARCH component-specific parameters
  if (!is.null(component) && model_type == "gogarch") {
    if (!is.null(params$garch_pars) && component <= length(params$garch_pars)) {
      garch <- params$garch_pars[[component]]
      if (param_name %in% names(garch)) {
        return(garch[[param_name]])
      }
    }
    return(NA)
  }
  
  ## Handle series-specific GARCH parameters (for DCC/CGARCH univariate GARCH)
  if (!is.null(series) && !is.null(params$garch_pars)) {
    if (series <= length(params$garch_pars)) {
      garch <- params$garch_pars[[series]]
      if (param_name %in% names(garch)) {
        return(garch[[param_name]])
      }
    }
    return(NA)
  }
  
  ## Top-level parameters (DCC/CGARCH correlation parameters)
  if (param_name %in% names(params)) {
    return(params[[param_name]])
  }
  
  ## Check in dcc_pars sublist
  if (!is.null(params$dcc_pars) && param_name %in% names(params$dcc_pars)) {
    return(params$dcc_pars[[param_name]])
  }
  
  ## Check in dist_pars sublist (for shape parameter)
  if (!is.null(params$dist_pars) && param_name %in% names(params$dist_pars)) {
    return(params$dist_pars[[param_name]])
  }
  
  ## CGARCH-specific: gamma in adcc_pars
  if (!is.null(params$adcc_pars) && param_name %in% names(params$adcc_pars)) {
    return(params$adcc_pars[[param_name]])
  }
  
  ## GOGARCH-specific: check ica_info
  if (!is.null(params$ica_info) && param_name %in% names(params$ica_info)) {
    return(params$ica_info[[param_name]])
  }
  
  ## Flattened search as last resort
  params_flat <- unlist(params)
  param_match <- grep(paste0("\\b", param_name, "$"), names(params_flat), value = TRUE)
  if (length(param_match) > 0) {
    return(params_flat[[param_match[1]]])
  }
  
  NA
}


#' @keywords internal
.extract_all_params <- function(state_data, iterations, model_type) {
  
  n_iter <- length(state_data)
  
  if (model_type == "gogarch") {
    ## GOGARCH: extract component GARCH parameters
    first_params <- state_data[[1]]$parameters
    n_components <- length(first_params$garch_pars)
    
    if (n_components == 0) {
      warning("No GARCH parameters found for GOGARCH model")
      return(NULL)
    }
    
    result <- data.frame(iteration = iterations)
    
    for (c in 1:n_components) {
      result[[paste0("component_", c, "_omega")]] <- sapply(state_data, function(x) {
        x$parameters$garch_pars[[c]]$omega %||% NA
      })
      result[[paste0("component_", c, "_alpha")]] <- sapply(state_data, function(x) {
        x$parameters$garch_pars[[c]]$alpha1 %||% NA
      })
      result[[paste0("component_", c, "_beta")]] <- sapply(state_data, function(x) {
        x$parameters$garch_pars[[c]]$beta1 %||% NA
      })
      result[[paste0("component_", c, "_persistence")]] <- 
        result[[paste0("component_", c, "_alpha")]] + 
        result[[paste0("component_", c, "_beta")]]
    }
    
    ## Add ICA info if available
    result$ica_method <- sapply(state_data, function(x) {
      x$parameters$ica_info$method %||% NA
    })
    
    return(result)
  }
  
  ## DCC / CGARCH: extract correlation dynamics parameters
  result <- data.frame(
    iteration = iterations,
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
  
  ## CGARCH-specific: gamma (ADCC leverage)
  if (model_type == "cgarch") {
    gamma_vals <- sapply(state_data, function(x) {
      p <- x$parameters
      p$gamma_1 %||% p$dcc_pars$gamma_1 %||% p$adcc_pars$gamma_1 %||% NA
    })
    if (!all(is.na(gamma_vals))) {
      result$gamma <- gamma_vals
      ## Adjust persistence for ADCC
      result$persistence <- result$alpha + result$beta + 0.5 * result$gamma
    }
  }
  
  ## Shape parameter (MVT distribution for DCC, MVT copula for CGARCH)
  shape_vals <- sapply(state_data, function(x) {
    p <- x$parameters
    p$shape %||% p$dist_pars$shape %||% NA
  })
  if (!all(is.na(shape_vals))) {
    result$shape <- shape_vals
  }
  
  ## CGARCH-specific: copula and transformation info
  if (model_type == "cgarch") {
    copula_vals <- sapply(state_data, function(x) {
      x$parameters$copula %||% NA
    })
    if (!all(is.na(copula_vals))) {
      result$copula <- copula_vals
    }
    
    transform_vals <- sapply(state_data, function(x) {
      x$parameters$transformation %||% x$parameters$pit_method %||% NA
    })
    if (!all(is.na(transform_vals))) {
      result$transformation <- transform_vals
    }
  }
  
  ## Remove columns that are all NA
  all_na_cols <- sapply(result, function(col) all(is.na(col)))
  result <- result[, !all_na_cols, drop = FALSE]
  
  result
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


## === PROBLEM DETECTION FUNCTIONS ===

#' Identify Problematic States
#'
#' Identify states with estimation problems based on diagnostic information.
#' Supports DCC, CGARCH, and GOGARCH models with model-specific checks.
#'
#' @param diagnostics An object of class \code{ms_diagnostics}
#' @param state Integer state index. If NULL (default), checks all states.
#' @param model_type Character: model type override. One of "auto" (default), "dcc", 
#'   "cgarch", or "gogarch". When "auto", attempts to detect from diagnostics.
#'
#' @return List with:
#'   \item{has_problems}{Logical indicating if any problems were found}
#'   \item{n_states_affected}{Number of states with problems}
#'   \item{problems}{Named list with problem descriptions per state}
#'
#' @details
#' The function performs different checks depending on model type:
#' \describe{
#'   \item{DCC}{
#'     \itemize{
#'       \item High persistence (alpha + beta > 0.98)
#'       \item Constant correlation fallback
#'       \item Parameter instability in final iterations
#'       \item Boundary events
#'     }
#'   }
#'   \item{CGARCH}{
#'     \itemize{
#'       \item All DCC checks
#'       \item Copula shape parameter issues (MVT: df < 3 or > 100)
#'       \item ADCC gamma constraints
#'       \item PIT transformation warnings
#'     }
#'   }
#'   \item{GOGARCH}{
#'     \itemize{
#'       \item ICA convergence failure
#'       \item Mixing matrix ill-conditioning (condition number > 1000)
#'       \item Unmixing matrix near-singularity
#'       \item Component correlation (should be < 0.2)
#'       \item Component GARCH high persistence
#'     }
#'   }
#' }
#'
#' @examples
#' \dontrun{
#' # Check all states with auto-detection
#' problems <- identify_problematic_states(diag)
#' if (problems$has_problems) {
#'   print(problems$problems)
#' }
#' 
#' # Check specific state for CGARCH model
#' problems <- identify_problematic_states(diag, state = 1, model_type = "cgarch")
#' }
#'
#' @export
identify_problematic_states <- function(diagnostics, state = NULL,
                                        model_type = c("auto", "dcc", "cgarch", "gogarch")) {
  
  if (!inherits(diagnostics, "ms_diagnostics")) {
    stop("diagnostics must be an object of class 'ms_diagnostics'")
  }
  
  model_type <- match.arg(model_type)
  
  ## Determine which states to check
  if (is.null(state)) {
    state_keys <- names(diagnostics$parameter_evolution)
    states_to_check <- as.integer(gsub("state_", "", state_keys))
  } else {
    states_to_check <- state
  }
  
  ## Auto-detect model type from first state if needed
  if (model_type == "auto" && length(diagnostics$parameter_evolution) > 0) {
    first_state <- diagnostics$parameter_evolution[[1]]
    if (length(first_state) > 0) {
      model_type <- .detect_model_type_from_params(first_state[[1]]$parameters)
    } else {
      model_type <- "dcc"  # Default
    }
  }
  
  problems <- list()
  
  for (s in states_to_check) {
    state_key <- paste0("state_", s)
    state_data <- diagnostics$parameter_evolution[[state_key]]
    
    if (is.null(state_data) || length(state_data) == 0) next
    
    final_params <- state_data[[length(state_data)]]$parameters
    state_problems <- character(0)
    
    ## Model-specific checks
    if (model_type == "gogarch") {
      state_problems <- c(state_problems, 
                          .check_gogarch_problems(final_params, state_data))
    } else if (model_type == "cgarch") {
      state_problems <- c(state_problems, 
                          .check_cgarch_problems(final_params, state_data))
    } else {
      ## DCC checks (also included in CGARCH)
      state_problems <- c(state_problems, 
                          .check_dcc_problems(final_params, state_data))
    }
    
    ## Common checks: boundary events
    boundary_events <- Filter(
      function(b) b$state == s,
      diagnostics$boundary_events
    )
    if (length(boundary_events) > 0) {
      n_alpha_boundary <- sum(sapply(boundary_events, function(b) grepl("alpha", b$parameter)))
      n_beta_boundary <- sum(sapply(boundary_events, function(b) grepl("beta", b$parameter)))
      
      if (n_alpha_boundary > 0 || n_beta_boundary > 0) {
        state_problems <- c(state_problems, 
                            sprintf("Boundary events: %d alpha, %d beta", 
                                    n_alpha_boundary, n_beta_boundary))
      }
    }
    
    ## Common checks: excessive warnings
    state_warnings <- Filter(
      function(w) (w$details$state %||% NA) == s,
      diagnostics$warnings
    )
    if (length(state_warnings) > 5) {
      state_problems <- c(state_problems,
                          sprintf("Excessive warnings (%d)", length(state_warnings)))
    }
    
    if (length(state_problems) > 0) {
      problems[[state_key]] <- state_problems
    }
  }
  
  list(
    has_problems = length(problems) > 0,
    n_states_affected = length(problems),
    problems = problems,
    model_type = model_type
  )
}


#' @keywords internal
.check_dcc_problems <- function(final_params, state_data) {
  problems <- character(0)
  
  ## Extract DCC parameters
  alpha <- final_params$alpha_1 %||% final_params$dcc_pars$alpha_1
  beta <- final_params$beta_1 %||% final_params$dcc_pars$beta_1
  
  ## Constant correlation fallback
  if (!is.null(final_params$correlation_type) && 
      final_params$correlation_type == "constant") {
    problems <- c(problems, "DCC dynamics collapsed to constant correlation")
  }
  
  ## High persistence
  if (!is.null(alpha) && !is.null(beta)) {
    persistence <- alpha + beta
    if (persistence > 0.98) {
      problems <- c(problems, 
                    sprintf("Very high DCC persistence (%.4f)", persistence))
    } else if (persistence > 0.95) {
      problems <- c(problems, 
                    sprintf("High DCC persistence (%.4f)", persistence))
    }
  }
  
  ## Parameter instability in final iterations
  if (length(state_data) >= 5) {
    last_5 <- tail(state_data, 5)
    
    if (!is.null(alpha)) {
      alpha_vals <- sapply(last_5, function(x) {
        x$parameters$alpha_1 %||% x$parameters$dcc_pars$alpha_1 %||% NA
      })
      if (any(!is.na(alpha_vals))) {
        alpha_range <- max(alpha_vals, na.rm = TRUE) - min(alpha_vals, na.rm = TRUE)
        if (alpha_range > 0.05) {
          problems <- c(problems,
                        sprintf("alpha unstable in final iterations (range: %.4f)", alpha_range))
        }
      }
    }
    
    if (!is.null(beta)) {
      beta_vals <- sapply(last_5, function(x) {
        x$parameters$beta_1 %||% x$parameters$dcc_pars$beta_1 %||% NA
      })
      if (any(!is.na(beta_vals))) {
        beta_range <- max(beta_vals, na.rm = TRUE) - min(beta_vals, na.rm = TRUE)
        if (beta_range > 0.05) {
          problems <- c(problems,
                        sprintf("beta unstable in final iterations (range: %.4f)", beta_range))
        }
      }
    }
  }
  
  ## Shape parameter issues (MVT distribution)
  shape <- final_params$shape %||% final_params$dist_pars$shape
  if (!is.null(shape)) {
    if (shape < 3) {
      problems <- c(problems, 
                    sprintf("Shape parameter (%.2f) < 3: infinite variance", shape))
    }
    if (shape > 100) {
      problems <- c(problems, 
                    "Shape parameter very high: effectively Gaussian")
    }
  }
  
  problems
}


#' @keywords internal
.check_cgarch_problems <- function(final_params, state_data) {
  ## Start with DCC checks (CGARCH includes DCC dynamics)
  problems <- .check_dcc_problems(final_params, state_data)
  
  ## ADCC gamma issues
  gamma <- final_params$gamma_1 %||% 
    final_params$dcc_pars$gamma_1 %||% 
    final_params$adcc_pars$gamma_1
  
  if (!is.null(gamma)) {
    if (gamma > 0.3) {
      problems <- c(problems, 
                    sprintf("ADCC gamma high (%.3f): strong leverage effect", gamma))
    }
    if (gamma < 0) {
      problems <- c(problems,
                    sprintf("ADCC gamma negative (%.3f): constraint violation", gamma))
    }
    
    ## Adjust persistence check for ADCC
    alpha <- final_params$alpha_1 %||% final_params$dcc_pars$alpha_1 %||% 0
    beta <- final_params$beta_1 %||% final_params$dcc_pars$beta_1 %||% 0
    adcc_persistence <- alpha + beta + 0.5 * gamma
    
    if (adcc_persistence > 0.98) {
      problems <- c(problems, 
                    sprintf("Very high ADCC persistence (%.4f)", adcc_persistence))
    }
  }
  
  ## Copula shape parameter (stricter than distribution shape)
  shape <- final_params$shape %||% final_params$dist_pars$shape
  copula <- final_params$copula %||% "mvn"
  
  if (copula == "mvt" && !is.null(shape)) {
    if (shape < 4) {
      problems <- c(problems, 
                    sprintf("MVT copula df (%.1f) < 4: very heavy tails", shape))
    }
  }
  
  ## PIT transformation warnings
  if (!is.null(final_params$pit_warning) && final_params$pit_warning) {
    problems <- c(problems, "PIT transformation generated warnings")
  }
  
  problems
}


#' @keywords internal
.check_gogarch_problems <- function(final_params, state_data) {
  problems <- character(0)
  
  ## ICA-specific checks
  if (!is.null(final_params$ica_info)) {
    
    ## ICA convergence
    if (isTRUE(final_params$ica_info$convergence_failed)) {
      problems <- c(problems, "ICA decomposition failed to converge")
    }
    
    ## Mixing matrix conditioning
    A <- final_params$ica_info$A %||% final_params$mixing_matrix
    if (!is.null(A) && is.matrix(A)) {
      cond_num <- tryCatch({
        svd_A <- svd(A)
        max(svd_A$d) / max(min(svd_A$d), 1e-10)
      }, error = function(e) Inf)
      
      if (cond_num > 1000) {
        problems <- c(problems, 
                      sprintf("Mixing matrix ill-conditioned (cond=%.0f)", cond_num))
      } else if (cond_num > 100) {
        problems <- c(problems, 
                      sprintf("Mixing matrix moderately ill-conditioned (cond=%.0f)", cond_num))
      }
    }
    
    ## Unmixing matrix near-singularity
    W <- final_params$ica_info$W
    if (!is.null(W) && is.matrix(W)) {
      det_W <- tryCatch(abs(det(W)), error = function(e) NA)
      if (is.na(det_W) || det_W < 1e-10) {
        problems <- c(problems, "Unmixing matrix near-singular")
      }
    }
    
    ## Component correlation (should be near zero for ICA)
    S <- final_params$ica_info$S
    if (!is.null(S) && is.matrix(S) && ncol(S) > 1) {
      cor_S <- tryCatch(cor(S), error = function(e) NULL)
      if (!is.null(cor_S)) {
        max_cor <- max(abs(cor_S[upper.tri(cor_S)]))
        if (max_cor > 0.2) {
          problems <- c(problems, 
                        sprintf("ICA components correlated (max |r| = %.3f)", max_cor))
        }
      }
    }
  }
  
  ## Component GARCH checks
  if (!is.null(final_params$garch_pars)) {
    n_components <- length(final_params$garch_pars)
    high_persist_comps <- integer(0)
    
    for (i in seq_len(n_components)) {
      garch <- final_params$garch_pars[[i]]
      alpha <- garch$alpha1 %||% 0
      beta <- garch$beta1 %||% 0
      persistence <- alpha + beta
      
      if (persistence > 0.99) {
        high_persist_comps <- c(high_persist_comps, i)
      }
    }
    
    if (length(high_persist_comps) > 0) {
      problems <- c(problems, 
                    sprintf("Components %s have very high GARCH persistence (>0.99)",
                            paste(high_persist_comps, collapse = ", ")))
    }
  }
  
  ## Check for parameter instability across iterations
  if (length(state_data) >= 5 && !is.null(final_params$garch_pars)) {
    last_5 <- tail(state_data, 5)
    n_components <- length(final_params$garch_pars)
    
    for (c in seq_len(n_components)) {
      alpha_vals <- sapply(last_5, function(x) {
        x$parameters$garch_pars[[c]]$alpha1 %||% NA
      })
      if (any(!is.na(alpha_vals))) {
        alpha_range <- max(alpha_vals, na.rm = TRUE) - min(alpha_vals, na.rm = TRUE)
        if (alpha_range > 0.05) {
          problems <- c(problems,
                        sprintf("Component %d: alpha unstable (range: %.4f)", c, alpha_range))
        }
      }
    }
  }
  
  problems
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


#' Check for All-States-Constant Event
#'
#' Determines if and when all states switched to constant correlation.
#'
#' @param diagnostics An object of class \code{ms_diagnostics}
#'
#' @return List with:
#'   \item{occurred}{Logical indicating if all states became constant}
#'   \item{iteration}{Iteration number when it occurred (NA if not)}
#'   \item{message}{Description of the event}
#'
#' @export
check_all_states_constant <- function(diagnostics) {
  if (!inherits(diagnostics, "ms_diagnostics")) {
    stop("diagnostics must be an object of class 'ms_diagnostics'")
  }
  
  ## Look for the "all_states_constant" warning in diagnostics
  warnings <- diagnostics$warnings
  
  if (length(warnings) == 0) {
    return(list(
      occurred = FALSE,
      iteration = NA_integer_,
      message = "All states did not become constant during estimation"
    ))
  }
  
  ## Find the all_states_constant warning
  for (w in warnings) {
    if (!is.null(w$type) && w$type == "all_states_constant") {
      iter <- w$details$iteration %||% w$iteration %||% NA_integer_
      return(list(
        occurred = TRUE,
        iteration = iter,
        message = sprintf("All states became constant at iteration %d", iter)
      ))
    }
  }
  
  return(list(
    occurred = FALSE,
    iteration = NA_integer_,
    message = "All states did not become constant during estimation"
  ))
}


#### ______________________________________________________________________ ####
#### DCC-Specific Utility Functions                                         ####

## === SIMULATION FUNCTIONS FOR TESTING ===

#' Simulate DCC-GARCH Data with Higher-Order Support
#'
#' Simulate realistic multivariate time series data from a DCC(q,p)-GARCH process.
#' Supports arbitrary DCC orders through vector-valued alpha_dcc and beta_dcc.
#' This function is primarily intended for testing and examples.
#'
#' @param n Integer number of observations to simulate
#' @param k Integer number of series (default: 2)
#' @param omega Numeric vector of length k: GARCH intercepts (default: c(0.05, 0.08))
#' @param alpha_garch Numeric vector of length k: ARCH effects (default: c(0.10, 0.12))
#' @param beta_garch Numeric vector of length k: GARCH effects (default: c(0.85, 0.82))
#' @param alpha_dcc Numeric scalar or vector: DCC alpha parameters. 
#'   Length determines q (number of alpha lags). (default: 0.04)
#' @param beta_dcc Numeric scalar or vector: DCC beta parameters.
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
#' For backward compatibility, scalar alpha_dcc and beta_dcc produce DCC(1,1).
#' To simulate DCC(q,p), provide alpha_dcc as vector of length q and beta_dcc 
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
#'   alpha_dcc = 0.05,
#'   beta_dcc = c(0.50, 0.40),
#'   seed = 123
#' )
#' 
#' # DCC(2,2): two lags each
#' y_dcc22 <- simulate_dcc_garch(
#'   n = 300,
#'   alpha_dcc = c(0.03, 0.02),
#'   beta_dcc = c(0.50, 0.40),
#'   seed = 456
#' )
#' 
#' # Constant correlation (set alpha_dcc = 0)
#' y_const <- simulate_dcc_garch(n = 200, alpha_dcc = 0, beta_dcc = 0, seed = 789)
#'
#' @export
simulate_dcc_garch <- function(
  n,
  k = 2,
  omega = c(0.05, 0.08),
  alpha_garch = c(0.10, 0.12),
  beta_garch = c(0.85, 0.82),
  alpha_dcc = 0.04,
  beta_dcc = 0.93,
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
  alpha_dcc <- as.numeric(alpha_dcc)
  beta_dcc <- as.numeric(beta_dcc)
  
  q_order <- length(alpha_dcc)  ## Number of alpha lags (ARCH-like)
  p_order <- length(beta_dcc)   ## Number of beta lags (GARCH-like)
  
  if (any(alpha_dcc < 0) || any(beta_dcc < 0)) {
    stop("All alpha_dcc and beta_dcc values must be non-negative")
  }
  
  ## Check DCC stationarity: sum(alpha) + sum(beta) < 1
  dcc_persistence <- sum(alpha_dcc) + sum(beta_dcc)
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
        Q_new <- Q_new + alpha_dcc[j] * (z_lag %*% t(z_lag))
      }
      
      ## Add beta terms (lagged Q matrices)
      for (j in seq_len(p_order)) {
        Q_new <- Q_new + beta_dcc[j] * Q_history[[j]]
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


#' Simulate DCC-GARCH data for testing
#' @param n Number of observations
#' @param k Number of series
#' @param alpha DCC alpha parameter
#' @param beta DCC beta parameter
#' @param seed Random seed
#' @return List with y (data), true_params, Qbar
simulate_dcc_garch_test_data <<- function(
    n = 200, 
    k = 2,
    alpha = 0.05, 
    beta = 0.90,
    seed = 123
  ) {
  set.seed(seed)
  
  ## GARCH parameters for each series
  omega <- rep(0.05, k)
  alpha_garch <- rep(0.10, k)
  beta_garch <- rep(0.85, k)
  
  ## Unconditional correlation
  Rbar <- diag(k)
  if (k >= 2) {
    for (i in 1:(k-1)) {
      for (j in (i+1):k) {
        Rbar[i, j] <- Rbar[j, i] <- 0.5
      }
    }
  }
  
  ## Initialize
  y <- matrix(0, n, k)
  h <- matrix(omega / (1 - alpha_garch - beta_garch), n, k)
  Q <- array(0, dim = c(k, k, n))
  R <- array(0, dim = c(k, k, n))
  z <- matrix(0, n, k)
  
  Q[,,1] <- Rbar
  R[,,1] <- Rbar
  
  ## Simulate
  for (t in 1:n) {
    ## Univariate GARCH
    if (t > 1) {
      for (i in 1:k) {
        h[t, i] <- omega[i] + alpha_garch[i] * y[t-1, i]^2 + beta_garch[i] * h[t-1, i]
      }
    }
    
    ## DCC
    if (t > 1) {
      z_lag <- z[t-1, , drop = FALSE]
      Q[,,t] <- (1 - alpha - beta) * Rbar + alpha * (t(z_lag) %*% z_lag) + beta * Q[,,t-1]
      Q_diag_inv_sqrt <- diag(1/sqrt(diag(Q[,,t])), k)
      R[,,t] <- Q_diag_inv_sqrt %*% Q[,,t] %*% Q_diag_inv_sqrt
    }
    
    ## Generate correlated innovations
    L <- tryCatch(chol(R[,,t]), error = function(e) chol(R[,,t] + diag(1e-6, k)))
    eps <- rnorm(k)
    z[t, ] <- as.vector(t(L) %*% eps)
    
    ## Returns
    y[t, ] <- sqrt(h[t, ]) * z[t, ]
  }
  
  colnames(y) <- paste0("series_", 1:k)
  
  list(
    y = y,
    true_params = list(
      garch = list(omega = omega, alpha = alpha_garch, beta = beta_garch),
      dcc = list(alpha = alpha, beta = beta)
    ),
    Qbar = Rbar,
    std_resid = z
  )
}


#' Create a standard DCC spec for testing
create_test_dcc_spec <<- function(k = 2, distribution = "mvn") {
  univariate_specs <- lapply(1:k, function(i) {
    list(model = "garch", garch_order = c(1, 1), distribution = "norm")
  })
  
  list(
    garch_spec_fun = "dcc_modelspec",
    distribution = distribution,
    garch_spec_args = list(
      dcc_order = c(1, 1),
      dynamics = "dcc",
      garch_model = list(univariate = univariate_specs)
    )
  )
}


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
      alpha_dcc_start <- 0.05
      beta_dcc_start <- 0.90
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
      alpha_dcc_start <- 0.03 + state_factor * 0.07
      beta_dcc_start <- 0.94 - state_factor * 0.09
      
      # Ensure DCC stationarity
      while (alpha_dcc_start + beta_dcc_start >= 0.98) {
        alpha_dcc_start <- alpha_dcc_start * 0.95
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
        dcc_pars = list(alpha_1 = alpha_dcc_start, beta_1 = beta_dcc_start),
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
generate_single_state_dcc_spec <- function(
    k = 2,
    var_order = 1,
    garch_order = c(1, 1),
    distribution = "mvn",
    omega = NULL,
    alpha_garch = NULL,
    beta_garch = NULL,
    alpha_dcc = NULL,
    beta_dcc = NULL,
    seed = NULL
  ) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Default parameter values
  omega <- omega %||% rep(0.1, k)
  alpha_garch <- alpha_garch %||% rep(0.1, k)
  beta_garch <- beta_garch %||% rep(0.8, k)
  alpha_dcc <- alpha_dcc %||% 0.05
  beta_dcc <- beta_dcc %||% 0.90
  
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
        dcc_pars = list(alpha_1 = alpha_dcc, beta_1 = beta_dcc),
        dist_pars = dist_pars
      )
    )
  )
  
  return(spec_list)
}



#### ______________________________________________________________________ ####
#### CGARCH-Specific Utility Functions                                      ####

#' Extract CGARCH Parameter Trajectory
#'
#' Extract the evolution of CGARCH-specific parameters across EM iterations.
#' Includes DCC/ADCC dynamics parameters and copula parameters.
#'
#' @param diagnostics An object of class \code{ms_diagnostics}
#' @param state Integer state index (default 1)
#' @return Data frame with columns: iteration, alpha, beta, gamma (if ADCC), 
#'   shape (if MVT copula). Returns NULL if not a CGARCH model.
#'
#' @examples
#' \dontrun{
#' traj <- extract_cgarch_trajectory(diag, state = 1)
#' plot(traj$iteration, traj$alpha, type = "b")
#' }
#'
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
  
  ## Extract trajectories
  n_iter <- length(state_data)
  result <- data.frame(
    iteration = integer(n_iter),
    alpha = numeric(n_iter),
    beta = numeric(n_iter),
    gamma = numeric(n_iter),
    shape = numeric(n_iter)
  )
  
  for (i in seq_len(n_iter)) {
    params <- state_data[[i]]$parameters
    result$iteration[i] <- state_data[[i]]$iteration
    
    ## DCC parameters
    result$alpha[i] <- params$alpha_1 %||% params$dcc_pars$alpha_1 %||% NA
    result$beta[i] <- params$beta_1 %||% params$dcc_pars$beta_1 %||% NA
    
    ## ADCC gamma
    result$gamma[i] <- params$gamma_1 %||% params$dcc_pars$gamma_1 %||% 
      params$adcc_gamma %||% NA
    
    ## Shape parameter (MVT copula)
    result$shape[i] <- params$shape %||% params$dist_pars$shape %||% NA
  }
  
  ## Remove columns that are all NA
  if (all(is.na(result$gamma))) result$gamma <- NULL
  if (all(is.na(result$shape))) result$shape <- NULL
  
  result
}


#' Check CGARCH-Specific Problems
#'
#' Identify CGARCH-specific estimation issues including copula parameter problems,
#' PIT transformation warnings, and ADCC dynamics issues.
#'
#' @param diagnostics An object of class \code{ms_diagnostics}
#' @param state Integer state index (if NULL, checks all states)
#'
#' @return List with \code{has_problems} (logical) and detailed problem information
#'
#' @export
check_cgarch_problems <- function(diagnostics, state = NULL) {
  if (!inherits(diagnostics, "ms_diagnostics")) {
    stop("diagnostics must be an object of class 'ms_diagnostics'")
  }
  
  states_to_check <- if (is.null(state)) {
    state_keys <- names(diagnostics$parameter_evolution)
    as.integer(gsub("state_", "", state_keys))
  } else {
    state
  }
  
  problems <- list()
  
  for (s in states_to_check) {
    state_key <- paste0("state_", s)
    state_data <- diagnostics$parameter_evolution[[state_key]]
    
    if (is.null(state_data) || length(state_data) == 0) next
    
    final_params <- state_data[[length(state_data)]]$parameters
    state_problems <- character(0)
    
    ## Shape parameter issues (MVT copula)
    shape <- final_params$shape %||% final_params$dist_pars$shape
    if (!is.null(shape)) {
      if (shape < 3) {
        state_problems <- c(state_problems, 
                            sprintf("Copula shape parameter (%.2f) < 3 - infinite variance", shape))
      }
      if (shape > 100) {
        state_problems <- c(state_problems, 
                            "Copula shape very high - effectively Gaussian")
      }
    }
    
    ## ADCC gamma issues
    gamma <- final_params$gamma_1 %||% final_params$dcc_pars$gamma_1
    if (!is.null(gamma) && gamma > 0.3) {
      state_problems <- c(state_problems, 
                          sprintf("ADCC gamma high (%.3f) - strong leverage effect", gamma))
    }
    
    ## DCC persistence
    alpha <- final_params$alpha_1 %||% final_params$dcc_pars$alpha_1 %||% 0
    beta <- final_params$beta_1 %||% final_params$dcc_pars$beta_1 %||% 0
    persistence <- alpha + beta + 0.5 * (gamma %||% 0)  # ADCC adjustment
    
    if (persistence > 0.98) {
      state_problems <- c(state_problems, 
                          sprintf("CGARCH persistence very high (%.4f)", persistence))
    }
    
    ## Constant correlation fallback
    if (!is.null(final_params$correlation_type) && 
        final_params$correlation_type == "constant") {
      state_problems <- c(state_problems, 
                          "Dynamics collapsed to constant correlation")
    }
    
    if (length(state_problems) > 0) {
      problems[[state_key]] <- state_problems
    }
  }
  
  list(
    has_problems = length(problems) > 0,
    n_states_affected = length(problems),
    problems = problems
  )
}


#' Generate CGARCH Specification
#'
#' Generate a properly structured specification list for MS-VARMA-GARCH estimation
#' with CGARCH (Copula-GARCH) dynamics. Supports DCC, ADCC dynamics and MVT copula.
#'
#' @param M Integer number of states
#' @param k Integer number of series (default: 2)
#' @param var_order Integer VAR order (default: 1)
#' @param garch_order Integer vector of length 2: GARCH(p,q) order (default: c(1,1))
#' @param dynamics Character: correlation dynamics ("dcc", "adcc", or "constant")
#' @param copula Character: copula type ("mvn" or "mvt")
#' @param transformation Character: PIT method ("parametric", "empirical", or "spd")
#' @param seed Integer random seed for generating starting values
#' @param simple Logical: if TRUE, use identical starting values across states
#'
#' @return A list of length M containing properly formatted specifications
#'
#' @examples
#' \dontrun{
#' # Generate 2-state CGARCH specification with Student-t copula
#' spec <- generate_cgarch_spec(M = 2, k = 3, copula = "mvt")
#' 
#' # ADCC dynamics with SPD transformation
#' spec_adcc <- generate_cgarch_spec(M = 2, dynamics = "adcc", transformation = "spd")
#' }
#'
#' @export
generate_cgarch_spec <- function(M,
                                 k = 2,
                                 var_order = 1,
                                 garch_order = c(1, 1),
                                 dynamics = c("dcc", "adcc", "constant"),
                                 copula = c("mvn", "mvt"),
                                 transformation = c("parametric", "empirical", "spd"),
                                 seed = NULL,
                                 simple = FALSE) {
  
  dynamics <- match.arg(dynamics)
  copula <- match.arg(copula)
  transformation <- match.arg(transformation)
  
  if (!is.null(seed)) set.seed(seed)
  
  ## Base univariate GARCH spec
  spec_uni_garch <- list(
    model = "garch",
    garch_order = garch_order,
    distribution = "norm"
  )
  
  ## CGARCH-specific arguments
  cgarch_spec_args <- list(
    dcc_order = c(1, 1),
    dynamics = dynamics,
    copula = copula,
    transformation = transformation,
    garch_model = list(
      univariate = replicate(k, spec_uni_garch, simplify = FALSE)
    )
  )
  
  n_var_pars <- k * (1 + k * var_order)
  
  spec_list <- vector("list", M)
  
  for (j in 1:M) {
    if (simple) {
      omega_start <- rep(0.1, k)
      alpha_start <- rep(0.1, k)
      beta_start <- rep(0.8, k)
      alpha_dcc_start <- 0.05
      beta_dcc_start <- 0.90
      gamma_start <- if (dynamics == "adcc") 0.02 else NULL
    } else {
      state_factor <- (j - 1) / max(M - 1, 1)
      
      omega_start <- 0.05 + state_factor * 0.10 + rnorm(k, 0, 0.01)
      omega_start <- pmax(omega_start, 0.03)
      
      alpha_start <- 0.08 + state_factor * 0.07 + rnorm(k, 0, 0.01)
      alpha_start <- pmax(alpha_start, 0.05)
      
      beta_start <- 0.85 - state_factor * 0.10 + rnorm(k, 0, 0.01)
      beta_start <- pmin(beta_start, 0.88)
      
      for (i in 1:k) {
        while (alpha_start[i] + beta_start[i] >= 0.98) {
          alpha_start[i] <- alpha_start[i] * 0.95
        }
      }
      
      alpha_dcc_start <- 0.03 + state_factor * 0.07
      beta_dcc_start <- 0.94 - state_factor * 0.09
      gamma_start <- if (dynamics == "adcc") 0.02 + state_factor * 0.03 else NULL
    }
    
    garch_pars_list <- lapply(1:k, function(i) {
      list(
        omega = omega_start[i],
        alpha1 = alpha_start[i],
        beta1 = beta_start[i]
      )
    })
    
    ## Distribution parameters
    dist_pars <- if (copula == "mvt") {
      list(shape = 8.0)
    } else {
      NULL
    }
    
    ## DCC/ADCC parameters
    dcc_pars <- list(alpha_1 = alpha_dcc_start, beta_1 = beta_dcc_start)
    if (dynamics == "adcc" && !is.null(gamma_start)) {
      dcc_pars$gamma_1 <- gamma_start
    }
    
    spec_list[[j]] <- list(
      var_order = var_order,
      garch_spec_fun = "cgarch_modelspec",
      garch_spec_args = cgarch_spec_args,
      distribution = copula,
      start_pars = list(
        var_pars = rep(0.1, n_var_pars),
        garch_pars = garch_pars_list,
        dcc_pars = dcc_pars,
        dist_pars = dist_pars
      )
    )
  }
  
  spec_list
}


#' Simulate CGARCH Data with ADCC Support
#'
#' Simulate data from a Copula-GARCH process with DCC or ADCC dynamics.
#'
#' @param n Integer number of observations
#' @param k Integer number of series (default: 2)
#' @param omega Numeric vector of GARCH omega parameters (length k)
#' @param alpha_garch Numeric vector of GARCH alpha parameters (length k)
#' @param beta_garch Numeric vector of GARCH beta parameters (length k)
#' @param alpha_dcc Numeric: DCC alpha parameter (default: 0.04)
#' @param beta_dcc Numeric: DCC beta parameter (default: 0.93)
#' @param gamma_dcc Numeric: ADCC gamma parameter for leverage (default: NULL for standard DCC)
#' @param copula Character: copula type ("mvn" or "mvt")
#' @param shape Numeric: degrees of freedom for MVT copula (default: 8)
#' @param Qbar Matrix: unconditional correlation matrix. If NULL, uses moderate correlation.
#' @param seed Integer random seed
#' @return Matrix of simulated returns (n x k)
#' 
#' @details
#' When gamma_dcc is provided (non-NULL and non-zero), the function simulates from
#' an ADCC (Asymmetric DCC) process where negative shocks have a larger impact on
#' correlation dynamics:
#' 
#' Q_t = Omega + alpha * (z_{t-1} z'_{t-1}) + gamma * (n_{t-1} n'_{t-1}) + beta * Q_{t-1}
#' 
#' where n_t = z_t * I(z_t < 0) captures negative shocks only.
#' 
#' @export
simulate_cgarch <- function(n,
                            k = 2,
                            omega = NULL,
                            alpha_garch = NULL,
                            beta_garch = NULL,
                            alpha_dcc = 0.04,
                            beta_dcc = 0.93,
                            gamma_dcc = NULL,
                            copula = "mvn",
                            shape = 8,
                            Qbar = NULL,
                            seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  ## Default parameters
  if (is.null(omega)) omega <- rep(0.05, k)
  if (is.null(alpha_garch)) alpha_garch <- rep(0.10, k)
  if (is.null(beta_garch)) beta_garch <- rep(0.85, k)
  
  ## Validate inputs
  if (length(omega) != k || length(alpha_garch) != k || length(beta_garch) != k) {
    stop("omega, alpha_garch, and beta_garch must have length k")
  }
  
  ## Check GARCH stationarity
  if (any(alpha_garch + beta_garch >= 1)) {
    stop("GARCH process must be stationary: alpha + beta < 1")
  }
  
  ## Check DCC/ADCC stationarity
  use_adcc <- !is.null(gamma_dcc) && gamma_dcc != 0
  if (use_adcc) {
    adcc_persistence <- alpha_dcc + beta_dcc + 0.5 * gamma_dcc
    if (adcc_persistence >= 1) {
      stop("ADCC process must be stationary: alpha + beta + 0.5*gamma < 1")
    }
  } else {
    if (alpha_dcc + beta_dcc >= 1) {
      stop("DCC process must be stationary: alpha + beta < 1")
    }
  }
  
  ## Default Qbar: moderate correlation
  if (is.null(Qbar)) {
    Qbar <- matrix(0.5, k, k)
    diag(Qbar) <- 1
  }
  
  ## For ADCC, compute Nbar (unconditional covariance of negative shocks)
  ## We approximate this based on the expected value for standard normal
  ## E[n_i * n_j] where n = z * I(z < 0)
  ## For standard normal: E[z * I(z<0)] = -1/sqrt(2*pi)
  ## E[z^2 * I(z<0)] = 0.5 (half the variance)
  Nbar <- NULL
  if (use_adcc) {
    ## Approximate Nbar for standard normal margins
    Nbar <- Qbar * 0.5  ## Simplified approximation
    diag(Nbar) <- 0.5   ## E[z^2 | z < 0] * P(z < 0) for each marginal
  }
  
  ## Initialize storage
  y <- matrix(0, n, k)
  h <- matrix(0, n, k)  ## Conditional variances
  z <- matrix(0, n, k)  ## Standardized residuals (copula space)
  
  ## Initialize Q and R
  Q <- Qbar
  R <- Qbar
  
  ## Compute Omega for ADCC
  if (use_adcc) {
    Omega <- (1 - alpha_dcc - beta_dcc) * Qbar - gamma_dcc * Nbar
  }
  
  ## Initialize conditional variances at unconditional level
  for (i in 1:k) {
    h[1, i] <- omega[i] / (1 - alpha_garch[i] - beta_garch[i])
  }
  
  ## Main simulation loop
  for (t in 1:n) {
    ## Generate copula innovations
    if (copula == "mvn") {
      ## Multivariate normal
      L <- tryCatch(t(chol(R)), error = function(e) {
        ## Regularize if needed
        R_reg <- R + diag(1e-6, k)
        t(chol(R_reg))
      })
      eps <- rnorm(k)
      z[t, ] <- L %*% eps
    } else {
      ## Multivariate t
      L <- tryCatch(t(chol(R)), error = function(e) {
        R_reg <- R + diag(1e-6, k)
        t(chol(R_reg))
      })
      ## t-distributed: z = sqrt((nu-2)/chi^2_nu) * N(0, R)
      chi2 <- rchisq(1, df = shape)
      eps <- rnorm(k)
      z[t, ] <- sqrt((shape - 2) / chi2) * (L %*% eps)
    }
    
    ## Generate returns
    for (i in 1:k) {
      y[t, i] <- sqrt(h[t, i]) * z[t, i]
    }
    
    ## Update conditional variances for next period
    if (t < n) {
      for (i in 1:k) {
        h[t + 1, i] <- omega[i] + alpha_garch[i] * y[t, i]^2 + beta_garch[i] * h[t, i]
      }
    }
    
    ## Update DCC/ADCC dynamics for next period
    if (t < n) {
      z_lag <- z[t, , drop = FALSE]
      
      if (use_adcc) {
        ## ADCC update
        n_lag <- z_lag * (z_lag < 0)  ## Negative shock indicator
        
        Q_new <- Omega + 
          alpha_dcc * (t(z_lag) %*% z_lag) + 
          gamma_dcc * (t(n_lag) %*% n_lag) + 
          beta_dcc * Q
      } else {
        ## Standard DCC update
        Q_new <- (1 - alpha_dcc - beta_dcc) * Qbar + 
          alpha_dcc * (t(z_lag) %*% z_lag) + 
          beta_dcc * Q
      }
      
      Q <- Q_new
      
      ## Standardize to correlation matrix
      Q_diag <- diag(Q)
      if (any(Q_diag <= 0)) {
        ## Regularize
        Q <- Q + diag(1e-6, k)
        Q_diag <- diag(Q)
      }
      Q_diag_inv_sqrt <- diag(1 / sqrt(Q_diag), k)
      R <- Q_diag_inv_sqrt %*% Q %*% Q_diag_inv_sqrt
      
      ## Ensure valid correlation matrix
      diag(R) <- 1
    }
  }
  
  colnames(y) <- paste0("series_", 1:k)
  return(y)
}


#' @title Fit ADCC Copula Model
#' @description Fits an ADCC (Asymmetric DCC) copula model to standardized residuals.
#'   Wraps estimate_adcc_copula if available, otherwise provides standalone implementation.
#' @param z_matrix Matrix of copula residuals (T x k)
#' @param weights Observation weights
#' @param Qbar Unconditional covariance matrix
#' @param copula_dist Copula distribution ("mvn" or "mvt")
#' @param shape Fixed shape parameter (if NULL and copula_dist="mvt", estimated)
#' @return List with fitted parameters and convergence info
#' @keywords internal
fit_adcc_copula <- function(
    z_matrix,
    weights,
    Qbar,
    copula_dist = "mvn",
    shape = NULL
) {
  
  ## Try to use existing estimate_adcc_copula if available
  if (exists("estimate_adcc_copula", mode = "function")) {
    result <- tryCatch({
      estimate_adcc_copula(
        z_matrix = z_matrix,
        weights = weights,
        Qbar = Qbar,
        copula_dist = copula_dist
      )
    }, error = function(e) NULL)
    
    if (!is.null(result)) {
      ## Convert to expected format
      return(list(
        convergence = TRUE,
        alpha = result$alpha,
        beta = result$beta,
        gamma = result$gamma,
        shape = result$shape,
        nll = result$nll,
        std_errors = result$std_errors
      ))
    }
  }
  
  ## Fallback: standalone implementation
  T_obs <- nrow(z_matrix)
  k <- ncol(z_matrix)
  
  ## Compute Nbar
  neg_resid <- z_matrix * (z_matrix < 0)
  Nbar <- crossprod(neg_resid) / T_obs
  
  ## Starting values
  start_alpha <- 0.05
  start_gamma <- 0.02
  start_beta <- 0.90
  
  ## Objective function
  if (copula_dist == "mvn" || !is.null(shape)) {
    ## Fixed shape or Gaussian
    obj_fn <- function(params) {
      adcc_copula_nll(
        params = c(params, if (copula_dist == "mvt") shape else NULL),
        z_matrix = z_matrix,
        weights = weights,
        Qbar = Qbar,
        Nbar = Nbar,
        copula_dist = copula_dist
      )
    }
    start_params <- c(start_alpha, start_gamma, start_beta)
    param_names <- c("alpha", "gamma", "beta")
  } else {
    ## Estimate shape
    obj_fn <- function(params) {
      adcc_copula_nll(
        params = params,
        z_matrix = z_matrix,
        weights = weights,
        Qbar = Qbar,
        Nbar = Nbar,
        copula_dist = "mvt"
      )
    }
    start_params <- c(start_alpha, start_gamma, start_beta, 8)
    param_names <- c("alpha", "gamma", "beta", "shape")
  }
  
  ## Optimization with constraints
  lower_bounds <- c(0.001, 0.001, 0.001)
  upper_bounds <- c(0.3, 0.3, 0.99)
  
  if (copula_dist == "mvt" && is.null(shape)) {
    lower_bounds <- c(lower_bounds, 2.1)
    upper_bounds <- c(upper_bounds, 100)
  }
  
  fit <- tryCatch({
    optim(
      par = start_params,
      fn = obj_fn,
      method = "L-BFGS-B",
      lower = lower_bounds,
      upper = upper_bounds,
      hessian = TRUE
    )
  }, error = function(e) NULL)
  
  if (is.null(fit) || fit$convergence != 0) {
    return(list(convergence = FALSE))
  }
  
  ## Extract results
  result <- list(
    convergence = TRUE,
    alpha = fit$par[1],
    gamma = fit$par[2],
    beta = fit$par[3],
    nll = fit$value
  )
  
  if (copula_dist == "mvt") {
    result$shape <- if (is.null(shape)) fit$par[4] else shape
  }
  
  ## Compute standard errors from Hessian
  if (!is.null(fit$hessian)) {
    se <- tryCatch({
      sqrt(diag(solve(fit$hessian)))
    }, error = function(e) NULL)
    
    if (!is.null(se) && all(is.finite(se))) {
      names(se) <- param_names
      result$std_errors <- se
    }
  }
  
  result
}


#' @title Fit Standard CGARCH Copula Model
#' @description Fits a standard DCC copula model (without ADCC asymmetry).
#' @param z_matrix Matrix of copula residuals (T x k)
#' @param weights Observation weights
#' @param Qbar Unconditional covariance matrix
#' @param copula_dist Copula distribution ("mvn" or "mvt")
#' @return List with fitted parameters and convergence info
#' @keywords internal
fit_cgarch_copula <- function(
    z_matrix,
    weights,
    Qbar,
    copula_dist = "mvn"
) {
  
  T_obs <- nrow(z_matrix)
  k <- ncol(z_matrix)
  
  ## Check if copula_nll exists (from tsbs_cgarch.R)
  if (!exists("copula_nll", mode = "function")) {
    stop("copula_nll function not found. Ensure tsbs package is loaded.")
  }
  
  ## Starting values
  start_alpha <- 0.05
  start_beta <- 0.90
  
  ## Determine number of parameters
  if (copula_dist == "mvt") {
    start_params <- c(start_alpha, start_beta, 8)
    param_names <- c("alpha", "beta", "shape")
    lower_bounds <- c(0.001, 0.001, 2.1)
    upper_bounds <- c(0.3, 0.99, 100)
  } else {
    start_params <- c(start_alpha, start_beta)
    param_names <- c("alpha", "beta")
    lower_bounds <- c(0.001, 0.001)
    upper_bounds <- c(0.3, 0.99)
  }
  
  ## Objective function
  obj_fn <- function(params) {
    copula_nll(
      params = params,
      z_matrix = z_matrix,
      weights = weights,
      Qbar = Qbar,
      copula_dist = copula_dist,
      use_reparam = FALSE
    )
  }
  
  fit <- tryCatch({
    optim(
      par = start_params,
      fn = obj_fn,
      method = "L-BFGS-B",
      lower = lower_bounds,
      upper = upper_bounds,
      hessian = FALSE  ## We'll use cgarch_standard_errors for SE
    )
  }, error = function(e) NULL)
  
  if (is.null(fit) || fit$convergence != 0) {
    return(list(convergence = FALSE))
  }
  
  ## Extract results
  result <- list(
    convergence = TRUE,
    alpha = fit$par[1],
    beta = fit$par[2],
    nll = fit$value
  )
  
  if (copula_dist == "mvt") {
    result$shape <- fit$par[3]
  }
  
  result
}


#### ______________________________________________________________________ ####
#### GOGARCH-Specific Utility Functions                                     ####

#' Extract GOGARCH Parameter Trajectory
#'
#' Extract the evolution of GOGARCH-specific parameters across EM iterations.
#' Includes ICA decomposition info and component GARCH parameters.
#'
#' @param diagnostics An object of class \code{ms_diagnostics}
#' @param state Integer state index (default 1)
#' @param component Integer ICA component index (if NULL, returns summary)
#' @return Data frame with iteration and component parameter evolution.
#'   Returns NULL if not a GOGARCH model.
#'
#' @examples
#' \dontrun{
#' traj <- extract_gogarch_trajectory(diag, state = 1)
#' plot(traj$iteration, traj$component_1_persistence, type = "b")
#' }
#'
#' @export
extract_gogarch_trajectory <- function(diagnostics, state = 1, component = NULL) {
  
  if (!inherits(diagnostics, "ms_diagnostics")) {
    stop("diagnostics must be an object of class 'ms_diagnostics'")
  }
  
  state_key <- paste0("state_", state)
  state_data <- diagnostics$parameter_evolution[[state_key]]
  
  if (is.null(state_data) || length(state_data) == 0) {
    warning("No parameter evolution data found for state ", state)
    return(NULL)
  }
  
  ## Determine number of components from first iteration
  first_params <- state_data[[1]]$parameters
  n_components <- length(first_params$garch_pars)
  
  if (n_components == 0) {
    warning("No GARCH parameters found - may not be a GOGARCH model")
    return(NULL)
  }
  
  ## Extract trajectories
  n_iter <- length(state_data)
  
  if (!is.null(component)) {
    ## Single component trajectory
    if (component > n_components) {
      stop("Component ", component, " not found. Model has ", n_components, " components.")
    }
    
    result <- data.frame(
      iteration = sapply(state_data, function(x) x$iteration),
      omega = sapply(state_data, function(x) {
        x$parameters$garch_pars[[component]]$omega %||% NA
      }),
      alpha = sapply(state_data, function(x) {
        x$parameters$garch_pars[[component]]$alpha1 %||% NA
      }),
      beta = sapply(state_data, function(x) {
        x$parameters$garch_pars[[component]]$beta1 %||% NA
      })
    )
    result$persistence <- result$alpha + result$beta
    
  } else {
    ## Summary across all components
    result <- data.frame(
      iteration = sapply(state_data, function(x) x$iteration)
    )
    
    for (c in 1:n_components) {
      col_name <- paste0("component_", c, "_persistence")
      result[[col_name]] <- sapply(state_data, function(x) {
        gp <- x$parameters$garch_pars[[c]]
        (gp$alpha1 %||% 0) + (gp$beta1 %||% 0)
      })
    }
    
    ## Add ICA info if available
    result$ica_method <- sapply(state_data, function(x) {
      x$parameters$ica_info$method %||% NA
    })
  }
  
  result
}


#' Check GOGARCH-Specific Problems
#'
#' Identify GOGARCH-specific estimation issues including ICA decomposition problems,
#' component GARCH instability, and mixing matrix ill-conditioning.
#'
#' @param diagnostics An object of class \code{ms_diagnostics}
#' @param state Integer state index (if NULL, checks all states)
#'
#' @return List with \code{has_problems} (logical) and detailed problem information
#'
#' @export
check_gogarch_problems <- function(diagnostics, state = NULL) {
  if (!inherits(diagnostics, "ms_diagnostics")) {
    stop("diagnostics must be an object of class 'ms_diagnostics'")
  }
  
  states_to_check <- if (is.null(state)) {
    state_keys <- names(diagnostics$parameter_evolution)
    as.integer(gsub("state_", "", state_keys))
  } else {
    state
  }
  
  problems <- list()
  
  for (s in states_to_check) {
    state_key <- paste0("state_", s)
    state_data <- diagnostics$parameter_evolution[[state_key]]
    
    if (is.null(state_data) || length(state_data) == 0) next
    
    final_params <- state_data[[length(state_data)]]$parameters
    state_problems <- character(0)
    
    if (!is.null(final_params$ica_info)) {
      ## ICA convergence issues
      if (isTRUE(final_params$ica_info$convergence_failed)) {
        state_problems <- c(state_problems, "ICA decomposition failed to converge")
      }
      
      ## Check for near-singular mixing matrix
      A <- final_params$ica_info$A %||% final_params$mixing_matrix
      if (!is.null(A)) {
        cond_num <- tryCatch({
          svd_A <- svd(A)
          max(svd_A$d) / min(svd_A$d)
        }, error = function(e) Inf)
        
        if (cond_num > 1000) {
          state_problems <- c(state_problems, 
                              sprintf("ICA mixing matrix ill-conditioned (cond=%.0f)", cond_num))
        }
      }
      
      ## Check for near-singular unmixing matrix
      W <- final_params$ica_info$W
      if (!is.null(W)) {
        det_W <- tryCatch(det(W), error = function(e) NA)
        if (is.na(det_W) || abs(det_W) < 1e-10) {
          state_problems <- c(state_problems, "ICA unmixing matrix near-singular")
        }
      }
      
      ## Component correlation (should be near zero)
      S <- final_params$ica_info$S
      if (!is.null(S) && is.matrix(S) && ncol(S) > 1) {
        cor_S <- cor(S)
        max_cor <- max(abs(cor_S[upper.tri(cor_S)]))
        
        if (max_cor > 0.2) {
          state_problems <- c(state_problems, 
                              sprintf("ICA components correlated (max |r| = %.3f)", max_cor))
        }
      }
    }
    
    ## Component GARCH issues
    if (!is.null(final_params$garch_pars)) {
      high_persist_comps <- 0
      for (i in seq_along(final_params$garch_pars)) {
        garch <- final_params$garch_pars[[i]]
        persist <- (garch$alpha1 %||% 0) + (garch$beta1 %||% 0)
        if (persist > 0.98) {
          high_persist_comps <- high_persist_comps + 1
        }
      }
      
      if (high_persist_comps > 0) {
        state_problems <- c(state_problems, 
                            sprintf("%d ICA component(s) have very high GARCH persistence", 
                                    high_persist_comps))
      }
    }
    
    if (length(state_problems) > 0) {
      problems[[state_key]] <- state_problems
    }
  }
  
  list(
    has_problems = length(problems) > 0,
    n_states_affected = length(problems),
    problems = problems
  )
}


#' Generate GOGARCH Specification
#'
#' Generate a properly structured specification list for MS-VARMA-GARCH estimation
#' with GOGARCH (Generalized Orthogonal GARCH) dynamics.
#'
#' @param M Integer number of states
#' @param k Integer number of series (default: 3)
#' @param var_order Integer VAR order (default: 1)
#' @param garch_order Integer vector of length 2: GARCH(p,q) order (default: c(1,1))
#' @param ica_method Character: ICA algorithm ("radical" or "fastica")
#' @param n_components Integer: number of ICA components (default: k)
#' @param distribution Character: component distribution ("norm", "nig", or "gh")
#' @param seed Integer random seed
#' @param simple Logical: if TRUE, use identical starting values across states
#'
#' @return A list of length M containing properly formatted specifications
#'
#' @examples
#' \dontrun{
#' # Generate 2-state GOGARCH specification
#' spec <- generate_gogarch_spec(M = 2, k = 3)
#' 
#' # With FastICA and NIG distribution
#' spec_nig <- generate_gogarch_spec(M = 2, ica_method = "fastica", distribution = "nig")
#' }
#'
#' @export
generate_gogarch_spec <- function(M,
                                  k = 3,
                                  var_order = 1,
                                  garch_order = c(1, 1),
                                  ica_method = c("radical", "fastica"),
                                  n_components = k,
                                  distribution = c("norm", "nig", "gh"),
                                  seed = NULL,
                                  simple = FALSE) {
  
  ica_method <- match.arg(ica_method)
  distribution <- match.arg(distribution)
  
  if (!is.null(seed)) set.seed(seed)
  
  ## GOGARCH-specific arguments
  gogarch_spec_args <- list(
    model = "garch",
    order = garch_order,
    ica = ica_method,
    components = n_components
  )
  
  n_var_pars <- k * (1 + k * var_order)
  
  spec_list <- vector("list", M)
  
  for (j in 1:M) {
    if (simple) {
      omega_start <- rep(0.05, n_components)
      alpha_start <- rep(0.08, n_components)
      beta_start <- rep(0.88, n_components)
    } else {
      state_factor <- (j - 1) / max(M - 1, 1)
      
      omega_start <- 0.03 + state_factor * 0.05 + rnorm(n_components, 0, 0.005)
      omega_start <- pmax(omega_start, 0.02)
      
      alpha_start <- 0.06 + state_factor * 0.08 + rnorm(n_components, 0, 0.01)
      alpha_start <- pmax(alpha_start, 0.04)
      
      beta_start <- 0.90 - state_factor * 0.08 + rnorm(n_components, 0, 0.01)
      beta_start <- pmin(beta_start, 0.92)
      
      for (i in 1:n_components) {
        while (alpha_start[i] + beta_start[i] >= 0.98) {
          alpha_start[i] <- alpha_start[i] * 0.95
        }
      }
    }
    
    garch_pars_list <- lapply(1:n_components, function(i) {
      list(
        omega = omega_start[i],
        alpha1 = alpha_start[i],
        beta1 = beta_start[i]
      )
    })
    
    ## Distribution parameters
    dist_pars <- if (distribution %in% c("nig", "gh")) {
      list(shape = 1, skew = 0)
    } else {
      NULL
    }
    
    spec_list[[j]] <- list(
      var_order = var_order,
      garch_spec_fun = "gogarch_modelspec",
      garch_spec_args = gogarch_spec_args,
      distribution = distribution,
      start_pars = list(
        var_pars = rep(0.1, n_var_pars),
        garch_pars = garch_pars_list,
        dist_pars = dist_pars
      )
    )
  }
  
  spec_list
}


#' Simulate GOGARCH Data
#'
#' Simulate data from a Generalized Orthogonal GARCH process.
#'
#' @param n Integer number of observations
#' @param k Integer number of series/components (default: 3)
#' @param A Matrix: mixing matrix (k x k). If NULL, a random orthogonal matrix is used.
#' @param omega Numeric vector of component GARCH omega parameters
#' @param alpha_garch Numeric vector of component GARCH alpha parameters
#' @param beta_garch Numeric vector of component GARCH beta parameters
#' @param distribution Character: component distribution ("norm" or "std")
#' @param shape Numeric: degrees of freedom for "std" distribution (default: 8)
#' @param seed Integer random seed
#'
#' @return List with:
#'   \item{y}{Matrix of simulated returns (n x k)}
#'   \item{S}{Matrix of independent components (n x k)}
#'   \item{A}{Mixing matrix used}
#'   \item{W}{Unmixing matrix}
#'   \item{true_params}{List of true parameters}
#'
#' @examples
#' \dontrun{
#' # Simulate 500 observations from 3-series GOGARCH
#' result <- simulate_gogarch(n = 500, k = 3, seed = 42)
#' y <- result$y
#' }
#'
#' @export
simulate_gogarch <- function(n,
                             k = 3,
                             A = NULL,
                             omega = NULL,
                             alpha_garch = NULL,
                             beta_garch = NULL,
                             distribution = c("norm", "std"),
                             shape = 8,
                             seed = NULL) {
  
  distribution <- match.arg(distribution)
  
  if (!is.null(seed)) set.seed(seed)
  
  ## Default mixing matrix (random orthogonal)
  if (is.null(A)) {
    M <- matrix(rnorm(k * k), k, k)
    qr_M <- qr(M)
    A <- qr.Q(qr_M)
  }
  
  ## Default GARCH parameters (different for each component)
  if (is.null(omega)) {
    omega <- seq(0.03, 0.08, length.out = k)
  }
  if (is.null(alpha_garch)) {
    alpha_garch <- seq(0.05, 0.12, length.out = k)
  }
  if (is.null(beta_garch)) {
    beta_garch <- seq(0.90, 0.85, length.out = k)
  }
  
  if (length(omega) < k) omega <- rep(omega[1], k)
  if (length(alpha_garch) < k) alpha_garch <- rep(alpha_garch[1], k)
  if (length(beta_garch) < k) beta_garch <- rep(beta_garch[1], k)
  
  ## Simulate independent components with GARCH dynamics
  S <- matrix(0, n, k)
  h <- matrix(0, n, k)
  
  for (i in 1:k) {
    ## Initialize variance
    h[1, i] <- omega[i] / (1 - alpha_garch[i] - beta_garch[i])
    
    for (t in 1:n) {
      ## Draw innovation
      if (distribution == "norm") {
        eta <- rnorm(1)
      } else {
        eta <- rt(1, df = shape) / sqrt(shape / (shape - 2))
      }
      
      S[t, i] <- sqrt(h[t, i]) * eta
      
      if (t < n) {
        h[t+1, i] <- omega[i] + alpha_garch[i] * S[t, i]^2 + beta_garch[i] * h[t, i]
      }
    }
  }
  
  ## Mix to get observed returns
  y <- S %*% t(A)
  
  colnames(y) <- paste0("series_", 1:k)
  colnames(S) <- paste0("component_", 1:k)
  
  list(
    y = y,
    S = S,
    A = A,
    W = solve(A),
    true_params = list(
      omega = omega,
      alpha_garch = alpha_garch,
      beta_garch = beta_garch,
      distribution = distribution,
      shape = if (distribution == "std") shape else NULL
    )
  )
}


#' @title Fit GOGARCH Model
#' @description Fits a GOGARCH model to multivariate data.
#' @param y Matrix of returns (T x k)
#' @param k Number of components
#' @param ica_method ICA algorithm ("radical" or "fastica")
#' @param distribution Component distribution
#' @return List with fitted parameters and convergence info
#' @keywords internal
fit_gogarch_model <- function(
    y,
    k,
    ica_method = "radical",
    distribution = "norm"
) {
  
  ## Check if gogarch_modelspec and related functions exist
  ## Try to use tsmarch if available
  if (requireNamespace("tsmarch", quietly = TRUE) && 
      requireNamespace("tsgarch", quietly = TRUE)) {
    
    fit_result <- tryCatch({
      ## Create GOGARCH specification
      gogarch_spec <- tsmarch::gogarch_modelspec(
        y = y,
        model = "garch",
        order = c(1, 1),
        ica = ica_method
      )
      
      ## Estimate
      gogarch_fit <- tsmarch::estimate(gogarch_spec)
      
      ## Extract parameters
      coefs <- coef(gogarch_fit)
      
      ## Parse coefficient structure for GOGARCH
      garch_pars <- vector("list", k)
      for (c in 1:k) {
        omega_name <- paste0("omega.", c)
        alpha_name <- paste0("alpha1.", c)
        beta_name <- paste0("beta1.", c)
        
        garch_pars[[c]] <- list(
          omega = coefs[omega_name] %||% coefs[paste0("omega_", c)] %||% NA,
          alpha1 = coefs[alpha_name] %||% coefs[paste0("alpha1_", c)] %||% NA,
          beta1 = coefs[beta_name] %||% coefs[paste0("beta1_", c)] %||% NA
        )
      }
      
      ## Get mixing matrix
      A <- tryCatch(gogarch_fit$A, error = function(e) NULL)
      
      list(
        convergence = TRUE,
        ica_converged = TRUE,
        garch_pars = garch_pars,
        A = A,
        std_errors = NULL  ## Would need to extract from model
      )
    }, error = function(e) NULL)
    
    if (!is.null(fit_result)) return(fit_result)
  }
  
  ## Fallback: manual ICA + univariate GARCH fitting
  fit_manual_gogarch(y, k, ica_method, distribution)
}


#' @title Manual GOGARCH Fitting
#' @description Fallback GOGARCH fitting using manual ICA decomposition.
#' @keywords internal
fit_manual_gogarch <- function(y, k, ica_method = "radical", distribution = "norm") {
  
  n <- nrow(y)
  
  ## Step 1: ICA decomposition
  ica_result <- tryCatch({
    if (ica_method == "fastica" && requireNamespace("fastICA", quietly = TRUE)) {
      ica_fit <- fastICA::fastICA(y, n.comp = k)
      list(
        S = ica_fit$S,
        A = ica_fit$A,
        W = ica_fit$W,
        converged = TRUE
      )
    } else {
      ## Use built-in ICA if available, otherwise simple PCA-based approach
      if (exists("perform_ica", mode = "function")) {
        perform_ica(y, n_components = k, method = ica_method)
      } else {
        ## Simple PCA-based decomposition as fallback
        pca <- prcomp(y, center = TRUE, scale. = FALSE)
        S <- pca$x[, 1:k]
        A <- pca$rotation[, 1:k]
        list(
          S = S,
          A = A,
          W = t(A),
          converged = TRUE
        )
      }
    }
  }, error = function(e) {
    list(converged = FALSE)
  })
  
  if (!ica_result$converged) {
    return(list(convergence = FALSE, ica_converged = FALSE))
  }
  
  S <- ica_result$S
  A <- ica_result$A
  
  ## Step 2: Fit univariate GARCH to each component
  garch_pars <- vector("list", k)
  std_errors <- vector("list", k)
  all_converged <- TRUE
  
  for (c in 1:k) {
    s_c <- S[, c]
    
    ## Fit GARCH(1,1) to component
    garch_result <- fit_univariate_garch(s_c, distribution)
    
    if (!is.null(garch_result) && garch_result$convergence) {
      garch_pars[[c]] <- garch_result$params
      std_errors[[c]] <- garch_result$std_errors
    } else {
      all_converged <- FALSE
      garch_pars[[c]] <- list(omega = NA, alpha1 = NA, beta1 = NA)
      std_errors[[c]] <- list(omega = NA, alpha1 = NA, beta1 = NA)
    }
  }
  
  list(
    convergence = all_converged,
    ica_converged = TRUE,
    garch_pars = garch_pars,
    std_errors = std_errors,
    A = A,
    W = ica_result$W,
    S = S
  )
}


#' @title Fit Univariate GARCH(1,1)
#' @description Fits a simple GARCH(1,1) model to a univariate series.
#' @keywords internal
fit_univariate_garch <- function(y, distribution = "norm") {
  
  n <- length(y)
  
  ## Try tsgarch if available
  if (requireNamespace("tsgarch", quietly = TRUE)) {
    fit_result <- tryCatch({
      spec <- tsgarch::garch_modelspec(
        y = y,
        model = "garch",
        garch_order = c(1, 1),
        distribution = distribution
      )
      fit <- tsgarch::estimate(spec)
      coefs <- coef(fit)
      
      list(
        convergence = TRUE,
        params = list(
          omega = coefs["omega"],
          alpha1 = coefs["alpha1"],
          beta1 = coefs["beta1"]
        ),
        std_errors = NULL  ## Would need to extract
      )
    }, error = function(e) NULL)
    
    if (!is.null(fit_result)) return(fit_result)
  }
  
  ## Fallback: simple GARCH estimation using optim
  ## Objective: GARCH(1,1) log-likelihood
  garch_nll <- function(params) {
    omega <- params[1]
    alpha <- params[2]
    beta <- params[3]
    
    ## Constraints
    if (omega <= 0 || alpha < 0 || beta < 0 || alpha + beta >= 1) {
      return(1e10)
    }
    
    ## Initialize variance
    h <- numeric(n)
    h[1] <- var(y)
    
    for (t in 2:n) {
      h[t] <- omega + alpha * y[t-1]^2 + beta * h[t-1]
      if (h[t] <= 0) return(1e10)
    }
    
    ## Negative log-likelihood (normal)
    nll <- 0.5 * sum(log(h) + y^2 / h)
    
    if (!is.finite(nll)) return(1e10)
    nll
  }
  
  ## Starting values
  start_omega <- var(y) * 0.05
  start_alpha <- 0.08
  start_beta <- 0.88
  
  fit <- tryCatch({
    optim(
      par = c(start_omega, start_alpha, start_beta),
      fn = garch_nll,
      method = "L-BFGS-B",
      lower = c(1e-8, 1e-8, 1e-8),
      upper = c(1, 0.5, 0.999),
      hessian = TRUE
    )
  }, error = function(e) NULL)
  
  if (is.null(fit) || fit$convergence != 0) {
    return(list(convergence = FALSE))
  }
  
  ## Extract standard errors
  se <- tryCatch({
    sqrt(diag(solve(fit$hessian)))
  }, error = function(e) c(NA, NA, NA))
  
  list(
    convergence = TRUE,
    params = list(
      omega = fit$par[1],
      alpha1 = fit$par[2],
      beta1 = fit$par[3]
    ),
    std_errors = list(
      omega = se[1],
      alpha1 = se[2],
      beta1 = se[3]
    )
  )
}


#' @title Compute Mixing Matrix Angle
#' @description Computes the angle between estimated and true mixing matrices.
#' @param A_est Estimated mixing matrix
#' @param A_true True mixing matrix
#' @return Scalar angle in radians
#' @keywords internal
compute_mixing_angle <- function(A_est, A_true) {
  
  ## Handle sign and permutation ambiguity in ICA
  k <- ncol(A_est)
  
  ## Normalize columns
  A_est_norm <- apply(A_est, 2, function(x) x / sqrt(sum(x^2)))
  A_true_norm <- apply(A_true, 2, function(x) x / sqrt(sum(x^2)))
  
  ## Find best matching permutation (greedy)
  best_angle <- Inf
  
  ## Simple approach: compute Frobenius norm after optimal sign adjustment
  ## For each column of A_est, find best match in A_true
  total_angle <- 0
  used_cols <- logical(k)
  
  for (i in 1:k) {
    best_match <- 0
    best_cos <- -1
    
    for (j in 1:k) {
      if (used_cols[j]) next
      
      ## Cosine similarity (absolute value for sign ambiguity)
      cos_sim <- abs(sum(A_est_norm[, i] * A_true_norm[, j]))
      
      if (cos_sim > best_cos) {
        best_cos <- cos_sim
        best_match <- j
      }
    }
    
    used_cols[best_match] <- TRUE
    total_angle <- total_angle + acos(min(best_cos, 1))
  }
  
  total_angle / k  ## Average angle
}