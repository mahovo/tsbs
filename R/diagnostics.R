#' @title Diagnostic Data Collection System
#' @description Collects comprehensive diagnostic information during EM iterations
#' @keywords internal
create_diagnostic_collector <- function() {
  diagnostics <- list(
    em_iterations = list(),
    parameter_evolution = list(),
    sigma_evolution = list(),
    convergence_info = list(),
    warnings = list(),
    boundary_events = list()
  )
  
  class(diagnostics) <- "ms_diagnostics"
  return(diagnostics)
}


#' @title Add EM Iteration Diagnostics
#' @keywords internal
add_em_iteration_diagnostic <- function(diagnostics, iteration, log_lik_before, 
                                        log_lik_after, ll_change, duration_sec,
                                        parameters, convergence_flag = FALSE) {
  
  iter_data <- list(
    iteration = iteration,
    log_lik_before_mstep = log_lik_before,
    log_lik_after_mstep = log_lik_after,
    ll_change = ll_change,
    ll_decreased = (ll_change < -1e-6),
    duration_seconds = duration_sec,
    timestamp = Sys.time(),
    converged = convergence_flag
  )
  
  diagnostics$em_iterations[[length(diagnostics$em_iterations) + 1]] <- iter_data
  
  return(diagnostics)
}


#' @title Add Parameter Evolution Data
#' @keywords internal
add_parameter_evolution <- function(diagnostics, iteration, state, parameters) {
  
  param_data <- list(
    iteration = iteration,
    state = state,
    parameters = parameters
  )
  
  key <- paste0("state_", state)
  if (is.null(diagnostics$parameter_evolution[[key]])) {
    diagnostics$parameter_evolution[[key]] <- list()
  }
  
  diagnostics$parameter_evolution[[key]][[length(diagnostics$parameter_evolution[[key]]) + 1]] <- param_data
  
  return(diagnostics)
}


#' @title Add Sigma Evolution Data
#' @keywords internal
add_sigma_evolution <- function(diagnostics, iteration, state, series, sigma_summary) {
  
  sigma_data <- list(
    iteration = iteration,
    state = state,
    series = series,
    mean_sigma = sigma_summary$mean,
    sd_sigma = sigma_summary$sd,
    min_sigma = sigma_summary$min,
    max_sigma = sigma_summary$max,
    first_5 = sigma_summary$first_5,
    last_5 = sigma_summary$last_5
  )
  
  key <- paste0("state_", state, "_series_", series)
  if (is.null(diagnostics$sigma_evolution[[key]])) {
    diagnostics$sigma_evolution[[key]] <- list()
  }
  
  diagnostics$sigma_evolution[[key]][[length(diagnostics$sigma_evolution[[key]]) + 1]] <- sigma_data
  
  return(diagnostics)
}


#' @title Add Boundary Event
#' @keywords internal
add_boundary_event <- function(diagnostics, iteration, state, parameter_name, 
                               value, boundary_type, action_taken) {
  
  boundary_event <- list(
    iteration = iteration,
    state = state,
    parameter = parameter_name,
    value = value,
    boundary_type = boundary_type,  # "lower" or "upper"
    action_taken = action_taken,     # e.g., "constant_correlation_fallback"
    timestamp = Sys.time()
  )
  
  diagnostics$boundary_events[[length(diagnostics$boundary_events) + 1]] <- boundary_event
  
  return(diagnostics)
}


#' @title Add Warning
#' @keywords internal
add_diagnostic_warning <- function(diagnostics, iteration, warning_type, message, details = NULL) {
  
  warning_entry <- list(
    iteration = iteration,
    type = warning_type,
    message = message,
    details = details,
    timestamp = Sys.time()
  )
  
  diagnostics$warnings[[length(diagnostics$warnings) + 1]] <- warning_entry
  
  return(diagnostics)
}


#' @title Summarize Diagnostics
#' @export
summary.ms_diagnostics <- function(object, ...) {
  
  cat("=== MS-VARMA-GARCH Diagnostic Summary ===\n\n")
  
  # EM Iterations
  cat("EM ITERATIONS:\n")
  cat("  Total iterations:", length(object$em_iterations), "\n")
  
  if (length(object$em_iterations) > 0) {
    ll_changes <- sapply(object$em_iterations, function(x) x$ll_change)
    ll_final <- object$em_iterations[[length(object$em_iterations)]]$log_lik_after_mstep
    ll_initial <- object$em_iterations[[1]]$log_lik_before_mstep
    
    cat("  Initial LL:", ll_initial, "\n")
    cat("  Final LL:", ll_final, "\n")
    cat("  Total LL improvement:", ll_final - ll_initial, "\n")
    cat("  LL decreased in", sum(ll_changes < -1e-6), "iterations\n")
    cat("  Mean LL change per iteration:", mean(ll_changes), "\n")
    cat("  Min LL change:", min(ll_changes), "\n")
    cat("  Max LL change:", max(ll_changes), "\n")
    
    total_time <- sum(sapply(object$em_iterations, function(x) x$duration_seconds))
    cat("  Total computation time:", round(total_time, 2), "seconds\n")
  }
  
  # Boundary Events
  cat("\nBOUNDARY EVENTS:\n")
  cat("  Total boundary events:", length(object$boundary_events), "\n")
  
  if (length(object$boundary_events) > 0) {
    for (event in object$boundary_events) {
      cat("    Iteration", event$iteration, ": State", event$state, 
          "-", event$parameter, "=", round(event$value, 4),
          "at", event$boundary_type, "boundary ->", event$action_taken, "\n")
    }
  }
  
  # Warnings
  cat("\nWARNINGS:\n")
  cat("  Total warnings:", length(object$warnings), "\n")
  
  if (length(object$warnings) > 0) {
    warning_types <- table(sapply(object$warnings, function(x) x$type))
    for (wtype in names(warning_types)) {
      cat("   ", wtype, ":", warning_types[wtype], "\n")
    }
  }
  
  # Parameter Evolution
  cat("\nPARAMETER EVOLUTION:\n")
  cat("  States tracked:", length(object$parameter_evolution), "\n")
  
  # Sigma Evolution
  cat("\nSIGMA EVOLUTION:\n")
  cat("  Series tracked:", length(object$sigma_evolution), "\n")
  
  invisible(object)
}


#' @title Plot Diagnostics for MS-VARMA-GARCH Models
#' @description Creates diagnostic plots for EM algorithm convergence and parameter evolution.
#'
#' @param x An object of class \code{ms_diagnostics} returned by the fitting procedure.
#' @param type Character string specifying which plots to produce. One of:
#'   \describe{
#'     \item{\code{"all"}}{Produce all diagnostic plots (default).}
#'     \item{\code{"ll_evolution"}}{Plot log-likelihood evolution across EM iterations.}
#'     \item{\code{"parameters"}}{Plot parameter evolution across EM iterations.}
#'     \item{\code{"sigma"}}{Plot conditional volatility (sigma) evolution.}
#'   }
#' @param parameters Optional character vector of parameter names to include in the 
#'   parameter evolution plot. If \code{NULL} (default), all parameters are plotted.
#'   Supports regex patterns when a single string containing regex metacharacters
#'   (\code{^}, \code{$}, \code{.}, \code{*}, \code{[}) is provided.
#' @param normalize Logical. If \code{TRUE}, normalize parameter values to [0, 1] 
#'   within each parameter for easier comparison across different scales. 
#'   Default is \code{FALSE}.
#' @param quiet Logical. If \code{TRUE}, suppress warnings from numeric coercion
#'   during parameter extraction. Default is \code{FALSE}.
#' @param ... Additional arguments (currently unused).
#'
#' @details
#' The function produces up to three types of diagnostic visualizations:
#' 
#' \strong{Log-Likelihood Evolution} (\code{type = "ll_evolution"}):
#' Two plots showing (1) the log-likelihood value after each M-step, and 
#' (2) the change in log-likelihood per iteration. Useful for assessing 
#' convergence and detecting any non-monotonic behavior.
#'
#' \strong{Parameter Evolution} (\code{type = "parameters"}):
#' Faceted plot showing how each parameter evolves across EM iterations,
#' with separate colors for each regime state. The \code{parameters} argument
#' can filter to specific parameters of interest.
#'
#' \strong{Sigma Evolution} (\code{type = "sigma"}):
#' Faceted plot showing the mean conditional volatility (with ± 1 SD ribbon)
#' for each series across iterations, colored by state.
#'
#' @return Invisibly returns the ggplot object(s). Called primarily for side effects
#'   (printing plots).
#'
#' @examples
#' \dontrun{
#' # After fitting a model with diagnostics enabled
#' fit <- fit_ms_varma_garch(data, n_states = 2, collect_diagnostics = TRUE)
#' diagnostics <- attr(fit, "diagnostics")
#'
#' # Plot all diagnostics
#' plot(diagnostics)
#'
#' # Plot only log-likelihood evolution
#' plot(diagnostics, type = "ll_evolution")
#'
#' # Plot specific parameters
#' plot(diagnostics, type = "parameters", parameters = c("alpha_1", "beta_1"))
#'
#' # Plot parameters matching a regex pattern
#' plot(diagnostics, type = "parameters", parameters = "^alpha")
#'
#' # Normalized parameter plot for cross-parameter comparison
#' plot(diagnostics, type = "parameters", normalize = TRUE)
#' }
#'
#' @seealso \code{\link{summary.ms_diagnostics}} for text summaries of diagnostics.
#' @export
plot.ms_diagnostics <- function(x, type = c("all", "ll_evolution", "parameters", "sigma"), 
                                parameters = NULL, normalize = FALSE, quiet = FALSE, ...) {
  
  type <- match.arg(type)
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting diagnostics")
  }
  
  if (type %in% c("all", "ll_evolution")) {
    plot_ll_evolution(x)
  }
  
  if (type %in% c("all", "parameters")) {
    plot_parameter_evolution(x, parameters = parameters, normalize = normalize, quiet = quiet)
  }
  
  if (type %in% c("all", "sigma")) {
    plot_sigma_evolution(x)
  }
}


#' @keywords internal
plot_ll_evolution <- function(diagnostics) {
  
  if (length(diagnostics$em_iterations) == 0) {
    message("No EM iterations to plot")
    return(invisible(NULL))
  }
  
  df <- data.frame(
    iteration = sapply(diagnostics$em_iterations, function(x) x$iteration),
    ll_before = sapply(diagnostics$em_iterations, function(x) x$log_lik_before_mstep),
    ll_after = sapply(diagnostics$em_iterations, function(x) x$log_lik_after_mstep),
    ll_change = sapply(diagnostics$em_iterations, function(x) x$ll_change)
  )
  
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = iteration)) +
    ggplot2::geom_line(ggplot2::aes(y = ll_after, color = "After M-step"), size = 1) +
    ggplot2::geom_point(ggplot2::aes(y = ll_after), size = 2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    ggplot2::labs(title = "Log-Likelihood Evolution",
                  x = "EM Iteration",
                  y = "Log-Likelihood",
                  color = NULL) +
    ggplot2::theme_minimal()
  
  p2 <- ggplot2::ggplot(df, ggplot2::aes(x = iteration, y = ll_change)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ggplot2::labs(title = "Log-Likelihood Change per Iteration",
                  x = "EM Iteration",
                  y = "LL Change") +
    ggplot2::theme_minimal()
  
  print(p1)
  print(p2)
  
  invisible(list(p1, p2))
}


#' @keywords internal
plot_parameter_evolution <- function(
    diagnostics,
    parameters = NULL,
    normalize = FALSE,
    quiet = FALSE
  ) {
  
  if (length(diagnostics$parameter_evolution) == 0) {
    message("No parameter evolution to plot")
    return(invisible(NULL))
  }
  
  # Extract parameter data
  param_list <- list()
  for (state_key in names(diagnostics$parameter_evolution)) {
    state_data <- diagnostics$parameter_evolution[[state_key]]
    
    for (iter_data in state_data) {
      params <- unlist(iter_data$parameters)
      
      for (pname in names(params)) {
        val <- if (quiet) {
          suppressWarnings(as.numeric(params[[pname]]))
        } else {
          as.numeric(params[[pname]])
        }
        if (!is.na(val)) {
          param_list[[length(param_list) + 1]] <- data.frame(
            iteration = iter_data$iteration,
            state = iter_data$state,
            parameter = pname,
            value = val,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  if (length(param_list) == 0) {
    message("No parameter data extracted")
    return(invisible(NULL))
  }
  
  df <- do.call(rbind, param_list)
  
  # Filter parameters if requested
  if (!is.null(parameters)) {
    # Support regex patterns
    if (length(parameters) == 1 && grepl("^\\^|\\$|\\.|\\*|\\[", parameters)) {
      df <- df[grepl(parameters, df$parameter), ]
    } else {
      df <- df[df$parameter %in% parameters, ]
    }
    
    if (nrow(df) == 0) {
      message("No parameters matched the filter")
      return(invisible(NULL))
    }
  }
  
  # Normalize to [0, 1] within each parameter
  if (normalize) {
    df <- do.call(rbind, lapply(split(df, df$parameter), function(x) {
      rng <- range(x$value, na.rm = TRUE)
      if (is.finite(rng[1]) && is.finite(rng[2]) && diff(rng) > 0) {
        x$value <- (x$value - rng[1]) / diff(rng)
      } else {
        x$value <- 0.5
      }
      x
    }))
  }
  
  # Build plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = iteration, y = value, 
                                        color = as.factor(state))) +
    ggplot2::geom_line(linewidth = 0.3) +
    ggplot2::geom_point(size = 0.6) +
    ggplot2::facet_wrap(~ parameter) +
    ggplot2::labs(title = "Parameter Evolution Across EM Iterations",
                  x = "EM Iteration",
                  y = if (normalize) "Normalized Value" else "Parameter Value",
                  color = "State") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
  
  if (normalize) {
    p <- p + 
      ggplot2::scale_y_continuous(breaks = c(0, 0.5, 1), 
                                  labels = c("0.0", "0.5", "1.0"),
                                  limits = c(0, 1))
  }
  
  print(p)
  invisible(p)
}


#' @keywords internal
plot_sigma_evolution <- function(diagnostics) {
  
  if (length(diagnostics$sigma_evolution) == 0) {
    message("No sigma evolution to plot")
    return(invisible(NULL))
  }
  
  # Extract sigma data
  sigma_list <- list()
  for (key in names(diagnostics$sigma_evolution)) {
    series_data <- diagnostics$sigma_evolution[[key]]
    
    for (iter_data in series_data) {
      sigma_list[[length(sigma_list) + 1]] <- data.frame(
        iteration = iter_data$iteration,
        state = iter_data$state,
        series = iter_data$series,
        mean_sigma = iter_data$mean_sigma,
        sd_sigma = iter_data$sd_sigma,
        min_sigma = iter_data$min_sigma,
        max_sigma = iter_data$max_sigma,
        stringsAsFactors = FALSE
      )
    }
  }
  
  df <- do.call(rbind, sigma_list)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = iteration, y = mean_sigma,
                                        color = as.factor(state),
                                        group = interaction(state, series))) +
    ggplot2::geom_line(size = 0.5) +
    ggplot2::geom_point(size = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = mean_sigma - sd_sigma,
                                      ymax = mean_sigma + sd_sigma,
                                      fill = as.factor(state)),
                         alpha = 0.2) +
    ggplot2::facet_wrap(~ series, scales = "free_y") +
    ggplot2::labs(title = "Sigma Evolution Across EM Iterations",
                  x = "EM Iteration",
                  y = "Mean Sigma (± SD)",
                  color = "State",
                  fill = "State") +
    ggplot2::theme_minimal()
  
  print(p)
  
  invisible(p)
}


#' @title Export Diagnostics to File
#' @export
save_diagnostics <- function(diagnostics, filepath = "ms_diagnostics.rds") {
  saveRDS(diagnostics, file = filepath)
  message("Diagnostics saved to: ", filepath)
}


#' @title Load Diagnostics from File
#' @export
load_diagnostics <- function(filepath = "ms_diagnostics.rds") {
  diagnostics <- readRDS(filepath)
  class(diagnostics) <- "ms_diagnostics"
  return(diagnostics)
}
