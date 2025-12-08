#' Fit a Flexible N-State Markov-Switching Vector-ARMA-GARCH Model
#'
#' This function provides a user-friendly interface to fit a general Markov-Switching
#' model. It handles data validation, pre-processing, calls the C++ estimation
#' engine, and formats the results into a well-structured object.
#'
#' @param y A numeric matrix or data frame where rows are observations and
#'   columns are the time series variables.
#' @param M An integer specifying the number of states in the Markov chain.
#' @param d An integer specifying the order of differencing to be applied to
#'   the series before fitting. This order is constant across all states.
#'   Defaults to 0 (no differencing), which is appropriate for returns data.
#' @param spec A list of model specifications, one for each of the M states.
#'   See Details for the required structure.
#' @param model_type A character string, either "univariate" or "multivariate".
#'   Defaults to "univariate".
#' @param control A list of control parameters for the EM algorithm, including:
#' \itemize{
#'   \item{max_iter}{An integer specifying the maximum number of EM iterations. Defaults to 100.}
#'   \item{tol}{A numeric value specifying the convergence tolerance for the
#'     log-likelihood. Defaults to 1e-6.}
#' }
#' @param parallel Logical.
#' @param num_cores Number of cores.
#' @param collect_diagnostics Logical. Collect diagnostics or not.
#' @param verbose Logical. If TRUE, print detailed diagnostic information during 
#'   estimation. Default is FALSE.
#' @param verbose_file Character string specifying path to file for verbose output.
#'   If NULL (default), verbose output goes to console. If specified, all verbose
#'   output is written to this file instead. Only used if verbose = TRUE.
#'
#' @details
#' The `spec` argument is a list of length `M`, where each element `spec[[j]]`
#' defines the model for state `j`. This state-specific list must contain:
#' \itemize{
#'   \item{For the mean model: `arma_order` (univariate) or `var_order` (multivariate).}
#'   \item{For the variance model: `garch_model` (e.g., "garch"), `garch_order`,
#'     `distribution`, and for multivariate models, `garch_spec_fun` (e.g., "dcc_modelspec")
#'     and `garch_spec_args`.}
#'   \item{Starting parameters: `start_pars`, a list containing `arma_pars` (or `var_pars`)
#'     and `garch_pars`.}
#' }
#' 
#' For mathematical expression of the model see [ms_varma_garch_bs()].
#'
#' @return A list containing the full results of the estimation, including model
#'   fits for each state, the transition matrix, smoothed probabilities, and
#'   information criteria.
#'
#' @export
fit_ms_varma_garch <- function(
    y, 
    M, 
    d = 0, 
    spec,
    model_type = c("univariate", "multivariate"),
    control = list(),
    parallel = FALSE,
    num_cores = 1L,
    collect_diagnostics = FALSE,
    verbose = FALSE,
    verbose_file = NULL
  ) {
  
  ## --- 1. Argument and Data Validation ---
  model_type <- match.arg(model_type)
  
  if (!is.numeric(M) || M < 2 || M != round(M)) stop("'M' must be an integer >= 2.")
  if (!is.numeric(d) || d < 0 || d != round(d)) stop("'d' must be a non-negative integer.")
  if (!is.list(spec) || length(spec) != M) stop("'spec' must be a list of length M.")
  
  if (!is.numeric(y)) {
    stop("Input 'y' must be numeric.")
  }
  
  if (any(!is.finite(y))) {
    stop("Input matrix 'y' contains non-finite values (NA, NaN, Inf).")
  }
  
  if (!is.matrix(y)) {
    if(!is.data.frame(y)) {
      stop("Input 'y' must be a numeric matrix or data frame.")
    }
    y <- as.matrix(y)
  }
  
  ## --- Setup verbose output redirection ---
  if (verbose && !is.null(verbose_file)) {
    ## Open connection to file
    verbose_con <- file(verbose_file, open = "wt")
    
    ## Redirect cat() output to file
    ## Save original connection to restore later
    sink(verbose_con, type = "output")
    
    ## Ensure we close and restore on exit
    on.exit({
      sink(type = "output")  ## Restore console output
      close(verbose_con)
    }, add = TRUE)
    
    cat("=== MS-VARMA-GARCH Verbose Output ===\n")
    cat("Started:", format(Sys.time()), "\n")
    cat("Data dimensions:", nrow(y), "x", ncol(y), "\n")
    cat("Number of states:", M, "\n")
    cat("Model type:", model_type, "\n\n")
  }

  
  ## --- 2. Set Control Parameters ---
  ## Default control parameters
  control_defaults <- list(
    max_iter = 50,
    tol = 1e-4,
    dcc_boundary_threshold = 0.02,
    dcc_boundary_criterion = "bic",  # "threshold", "aic", or "bic"
    dcc_allow_refitting = TRUE
  )
  
  control <- modifyList(control_defaults, control)
  
  ## Validate control parameters
  if (!control$dcc_boundary_criterion %in% c("threshold", "aic", "bic")) {
    stop("control$dcc_boundary_criterion must be 'threshold', 'aic', or 'bic'")
  }
  
  if (control$dcc_boundary_threshold <= 0 || control$dcc_boundary_threshold >= 1) {
    stop("control$dcc_boundary_threshold must be between 0 and 1")
  }
  
  ## --- 3. Pre-processing: Handle Differencing ---
  if (d > 0) {
    if (nrow(y) <= d) {
      stop("The number of observations must be greater than the differencing order 'd'.")
    }
    y_diff <- diff(y, differences = d)
  } else {
    y_diff <- y
  }
  
  ## --- 4. Initialize Diagnostics Collector (if requested) ---
  diagnostics <- if (collect_diagnostics) {
    create_diagnostic_collector()
  } else {
    NULL
  }
  
  ## --- 5. Call the C++ Backend ---
  #if (verbose) message("Fitting the MS-ARMA-GARCH model via C++ EM algorithm...")
  if (verbose) {
    cat("\n=== MS-VARMA-GARCH Model Fitting ===\n")
    cat("Model type:", model_type, "\n")
    cat("Number of states:", M, "\n")
    cat("Sample size:", nrow(y_diff), "\n")
    cat("Control parameters:\n")
    cat("  max_iter:", control$max_iter, "\n")
    cat("  tol:", control$tol, "\n")
    cat("  dcc_boundary_threshold:", control$dcc_boundary_threshold, "\n")
    cat("  dcc_boundary_criterion:", control$dcc_boundary_criterion, "\n")
    cat("  dcc_allow_refitting:", control$dcc_allow_refitting, "\n\n")
    cat("Fitting the MS-ARMA-GARCH model via C++ EM algorithm...\n")
  }

  ## === INITIAL FIT: Dynamic correlation for all states ===
  fit_result <- fit_ms_varma_garch_cpp(
    y = y_diff,
    M = M,
    spec = spec,
    model_type = model_type,
    control = control,
    diagnostics = diagnostics,
    verbose = verbose
  )
  
  ## === SHORT-CIRCUIT: Check if any states hit boundary ===
  if (control$dcc_allow_refitting && model_type == "multivariate") {
    
    boundary_states <- detect_boundary_states(
      fit_result, 
      threshold = control$dcc_boundary_threshold,
      criterion = control$dcc_boundary_criterion
    )
    
    if (length(boundary_states) > 0) {
      if (verbose) {
        cat("\n=== BOUNDARY DETECTION ===\n")
        cat("States", paste(boundary_states, collapse = ", "), 
            "hit DCC parameter boundary.\n")
        cat("Refitting with constant correlation for these states...\n\n")
      }
      
      ## Create modified spec with constant correlation for boundary states
      spec_refit <- spec
      for (j in boundary_states) {
        spec_refit[[j]]$force_constant_correlation <- TRUE
      }
      
      ## Refit with mixed specification
      fit_refit <- fit_ms_varma_garch_cpp(
        y = y_diff,
        M = M,
        spec = spec_refit,
        model_type = model_type,
        control = control,
        diagnostics = if (collect_diagnostics) create_diagnostic_collector() else NULL,
        verbose = verbose
      )
      
      ## Compare models using appropriate criterion
      if (control$dcc_boundary_criterion %in% c("aic", "bic")) {
        # Count DCC parameters in each model
        k_original <- count_dcc_parameters(fit_result$model_fits)
        k_refit <- count_dcc_parameters(fit_refit$model_fits)
        
        n <- nrow(y_diff)
        
        if (control$dcc_boundary_criterion == "aic") {
          ic_original <- -2 * fit_result$log_likelihood + 2 * k_original
          ic_refit <- -2 * fit_refit$log_likelihood + 2 * k_refit
        } else {  # bic
          ic_original <- -2 * fit_result$log_likelihood + log(n) * k_original
          ic_refit <- -2 * fit_refit$log_likelihood + log(n) * k_refit
        }
        
        if (verbose) {
          cat("\n=== MODEL COMPARISON ===\n")
          cat("Original (all dynamic):\n")
          cat("  Log-likelihood:", fit_result$log_likelihood, "\n")
          cat("  DCC parameters:", k_original, "\n")
          cat(" ", toupper(control$dcc_boundary_criterion), ":", ic_original, "\n")
          cat("Refit (mixed constant/dynamic):\n")
          cat("  Log-likelihood:", fit_refit$log_likelihood, "\n")
          cat("  DCC parameters:", k_refit, "\n")
          cat(" ", toupper(control$dcc_boundary_criterion), ":", ic_refit, "\n")
          cat("Winner:", ifelse(ic_refit < ic_original, "REFIT", "ORIGINAL"), "\n\n")
        }
        
        if (ic_refit < ic_original) {
          fit_result <- fit_refit
        }
      } else {
        # For threshold criterion, always use refit
        fit_result <- fit_refit
      }
    }
  }
  
  if (verbose) message("Model fitting complete.")
  
  # ## --- 6. Extract diagnostics from C++ results (if collected) ---
  # if (collect_diagnostics) {
  #   diagnostics <- cpp_results$diagnostics
  # }
  # 
  # ## --- 7. Post-processing and Formatting Results ---
  # ## Align smoothed probabilities with the original time series
  # smoothed_probs_aligned <- matrix(NA_real_, nrow = T_orig, ncol = M)
  # colnames(smoothed_probs_aligned) <- paste0("State", 1:M)
  # 
  # ## The C++ output is aligned with the effective (differenced) data.
  # ## We need to pad it to match the original data's length.
  # padding <- T_orig - T_eff
  # if (padding > 0) {
  #   smoothed_probs_aligned[(padding + 1):T_orig, ] <- cpp_results$smoothed_probabilities[1:T_eff, ]
  # } else {
  #   smoothed_probs_aligned <- cpp_results$smoothed_probabilities
  # }
  # 
  # ## Calculate AIC/BIC
  # ## Count number of estimated parameters
  # num_mean_pars <- sum(unlist(lapply(cpp_results$model_fits, function(fit) length(fit$arma_pars %||% fit$var_pars))))
  # num_garch_pars <- sum(unlist(lapply(cpp_results$model_fits, function(fit) length(fit$garch_pars))))
  # num_trans_pars <- M * (M - 1) ## M*(M-1) free parameters in transition matrix
  # num_params <- num_mean_pars + num_garch_pars + num_trans_pars
  # 
  # aic <- -2 * cpp_results$log_likelihood + 2 * num_params
  # bic <- -2 * cpp_results$log_likelihood + log(T_eff) * num_params
  # 
  # ## --- Final message to verbose output ---
  # if (verbose && !is.null(verbose_file)) {
  #   cat("\n=== Fitting Complete ===\n")
  #   cat("Finished:", format(Sys.time()), "\n")
  #   cat("Final log-likelihood:", cpp_results$log_likelihood, "\n")
  # }
  # 
  # ## Assemble the final, user-friendly object
  # result <- list(
  #   model_fits = cpp_results$model_fits,
  #   P = cpp_results$P,
  #   log_likelihood = cpp_results$log_likelihood,
  #   smoothed_probabilities = smoothed_probs_aligned,
  #   aic = aic,
  #   bic = bic,
  #   d = d,
  #   y = y_orig,
  #   call = match.call(),
  #   convergence = cpp_results$convergence,
  #   warnings = cpp_results$warnings,
  #   diagnostics = diagnostics  ## ADD diagnostics to result
  # )
  # 
  # class(result) <- "msm.fit"
  # return(result)
  
  ## Prepare output
  result <- list(
    model_fits = fit_result$model_fits,
    P = fit_result$P,
    log_likelihood = fit_result$log_likelihood,
    smoothed_probabilities = fit_result$smoothed_probabilities,
    aic = fit_result$aic,
    bic = fit_result$bic,
    d = d,
    y = y,
    call = match.call(),
    convergence = fit_result$convergence,
    warnings = fit_result$warnings
  )
  
  if (collect_diagnostics) {
    result$diagnostics <- fit_result$diagnostics
  }
  
  #class(result) <- "msm.fit"
  #class(result) <- c("ms_varma_garch_fit", "list")
  
  return(result)
}


#' @title Detect States at DCC Parameter Boundary
#' @keywords internal
detect_boundary_states <- function(fit_result, threshold, criterion) {
  
  boundary_states <- integer(0)
  
  for (j in seq_along(fit_result$model_fits)) {
    state_fit <- fit_result$model_fits[[j]]
    
    # Skip if already constant
    if (!is.null(state_fit$correlation_type) && 
        state_fit$correlation_type == "constant") {
      next
    }
    
    # Check alpha parameters
    alpha_names <- grep("^alpha_[0-9]+$", names(state_fit), value = TRUE)
    
    if (length(alpha_names) > 0) {
      alpha_values <- sapply(alpha_names, function(nm) state_fit[[nm]])
      
      if (criterion == "threshold") {
        # Simple threshold check
        if (any(alpha_values < threshold)) {
          boundary_states <- c(boundary_states, j)
        }
      } else {
        # For AIC/BIC, always recheck boundary states
        # (comparison happens later)
        if (any(alpha_values < threshold)) {
          boundary_states <- c(boundary_states, j)
        }
      }
    }
  }
  
  return(boundary_states)
}


#' @title Count DCC Parameters in Model
#' @keywords internal
count_dcc_parameters <- function(model_fits) {
  
  k <- 0
  
  for (j in seq_along(model_fits)) {
    state_fit <- model_fits[[j]]
    
    # Count dynamic correlation parameters
    if (is.null(state_fit$correlation_type) || 
        state_fit$correlation_type == "dynamic") {
      # Count alpha and beta parameters
      alpha_names <- grep("^alpha_[0-9]+$", names(state_fit), value = TRUE)
      beta_names <- grep("^beta_[0-9]+$", names(state_fit), value = TRUE)
      k <- k + length(alpha_names) + length(beta_names)
    }
    # Constant correlation adds 0 parameters
  }
  
  return(k)
}


# Helper function for parameter counting
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}
