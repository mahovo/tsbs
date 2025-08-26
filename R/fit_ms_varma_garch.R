#' Fit a Flexible N-State Markov-Switching ARMA-GARCH Model
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
#'   \item{max_iter}{An integer specifying the maximum number of EM iterations. Defaults to 100.}
#'   \item{tol}{A numeric value specifying the convergence tolerance for the
#'     log-likelihood. Defaults to 1e-6.}
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
#' @return An object of class `msm.fit` containing the full results of the
#'   estimation, including model fits for each state, the transition matrix,
#'   smoothed probabilities, and information criteria.
#'
fit_ms_varma_garch <- function(y, M, d = 0, spec,
                               model_type = c("univariate", "multivariate"),
                               control = list(),
                               parallel = FALSE,
                               num_cores = 1L) {
  
  ## --- 1. Argument and Data Validation (Unchanged) ---
  model_type <- match.arg(model_type)
  if (!is.numeric(M) || M < 2 || M != round(M)) stop("'M' must be an integer >= 2.")
  if (!is.numeric(d) || d < 0 || d != round(d)) stop("'d' must be a non-negative integer.")
  if (!is.list(spec) || length(spec) != M) stop("'spec' must be a list of length M.")
  
  if (!is.matrix(y) && !is.data.frame(y)) {
    stop("Input 'y' must be a numeric matrix or data frame.")
  }
  y_mat <- as.matrix(y)
  if (!is.numeric(y_mat)) {
    stop("Input 'y' must be numeric.")
  }
  if (any(!is.finite(y_mat))) {
    stop("Input matrix 'y' contains non-finite values (NA, NaN, Inf).")
  }
  
  ## --- Setup Parallel Backend ---
  if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("The 'future.apply' package is required for parallel execution.")
    }
    ## Set up the parallel plan. multisession is generally safer.
    future::plan(future::multisession, workers = num_cores)
    ## Ensure the plan is reset when the function exits
    on.exit(future::plan(future::sequential), add = TRUE)
  }
  
  ## --- 2. Set Control Parameters ---
  ctrl <- list(max_iter = 100, tol = 1e-6)
  ctrl[names(control)] <- control 
  
  ## --- 3. Pre-processing: Handle Differencing ---
  y_orig <- as.matrix(y)
  T_orig <- nrow(y_orig)
  if (d > 0) {
    y_effective <- as.matrix(diff(y_orig, differences = d))
  } else {
    y_effective <- y_orig
  }
  T_eff <- nrow(y_effective)
  
  ## --- 4. Call the C++ Backend ---
  message("Fitting the MS-ARMA-GARCH model via C++ EM algorithm...")
  cpp_results <- fit_ms_varma_garch_cpp(
    y = y_effective,
    M = M,
    spec = spec,
    model_type = model_type,
    control = ctrl
  )
  message("Model fitting complete.")
  
  ## --- 5. Post-processing and Formatting Results ---
  ## Align smoothed probabilities with the original time series
  smoothed_probs_aligned <- matrix(NA_real_, nrow = T_orig, ncol = M)
  colnames(smoothed_probs_aligned) <- paste0("State", 1:M)
  
  ## The C++ output is aligned with the effective (differenced) data.
  ## We need to pad it to match the original data's length.
  padding <- T_orig - T_eff
  if (padding > 0) {
    smoothed_probs_aligned[(padding + 1):T_orig, ] <- cpp_results$smoothed_probabilities[1:T_eff, ]
  } else {
    smoothed_probs_aligned <- cpp_results$smoothed_probabilities
  }
  
  ## Calculate AIC/BIC
  ## Count number of estimated parameters
  num_mean_pars <- sum(unlist(lapply(cpp_results$model_fits, function(fit) length(fit$arma_pars %||% fit$var_pars))))
  num_garch_pars <- sum(unlist(lapply(cpp_results$model_fits, function(fit) length(fit$garch_pars))))
  num_trans_pars <- M * (M - 1) ## M*(M-1) free parameters in transition matrix
  num_params <- num_mean_pars + num_garch_pars + num_trans_pars
  
  aic <- -2 * cpp_results$log_likelihood + 2 * num_params
  bic <- -2 * cpp_results$log_likelihood + log(T_eff) * num_params
  
  ## Assemble the final, user-friendly object
  result <- list(
    model_fits = cpp_results$model_fits,
    P = cpp_results$P,
    log_likelihood = cpp_results$log_likelihood,
    smoothed_probabilities = smoothed_probs_aligned,
    aic = aic,
    bic = bic,
    d = d,
    y = y_orig,
    call = match.call(),
    convergence = cpp_results$convergence,
    warnings = cpp_results$warnings
  )
  
  class(result) <- "msm.fit"
  return(result)
}


#' Fit a Flexible N-State Markov-Switching ARMA-GARCH Model
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
#'   \item{max_iter}{An integer specifying the maximum number of EM iterations. Defaults to 100.}
#'   \item{tol}{A numeric value specifying the convergence tolerance for the
#'     log-likelihood. Defaults to 1e-6.}
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
#' @return An object of class `msm.fit` containing the full results of the
#'   estimation, including model fits for each state, the transition matrix,
#'   smoothed probabilities, and information criteria.
#'
fit_ms_varma_garch_old <- function(y, M, d = 0, spec,
                               model_type = c("univariate", "multivariate"),
                               control = list()) {
  
  # --- 1. Argument and Data Validation ---
  model_type <- match.arg(model_type)
  if (!is.numeric(M) || M < 2 || M != round(M)) stop("'M' must be an integer >= 2.")
  if (!is.numeric(d) || d < 0 || d != round(d)) stop("'d' must be a non-negative integer.")
  if (!is.list(spec) || length(spec) != M) stop("'spec' must be a list of length M.")
  
  if (!is.matrix(y) && !is.data.frame(y)) {
    stop("Input 'y' must be a numeric matrix or data frame.")
  }
  y_mat <- as.matrix(y)
  if (!is.numeric(y_mat)) {
    stop("Input 'y' must be numeric.")
  }
  if (any(!is.finite(y_mat))) {
    stop("Input matrix 'y' contains non-finite values (NA, NaN, Inf).")
  }
  
  # --- 2. Set Control Parameters ---
  ctrl <- list(max_iter = 100, tol = 1e-6)
  # Overwrite defaults with any user-provided values
  ctrl[names(control)] <- control 
  
  # --- 3. Pre-processing: Handle Differencing ---
  y_orig <- y_mat
  T_orig <- nrow(y_orig)
  
  if (d > 0) {
    if (T_orig <= d) {
      stop("The number of observations must be greater than the differencing order 'd'.")
    }
    # The core step: difference the data
    y_effective <- as.matrix(diff(y_orig, differences = d))
  } else {
    y_effective <- y_orig
  }
  
  T_eff <- nrow(y_effective)
  
  # --- 4. Call the C++ Backend ---
  # The C++ function performs the entire EM algorithm on y_effective.
  # It returns a raw list of results.
  message("Fitting the MS-ARMA-GARCH model via C++ EM algorithm...")
  cpp_results <- fit_ms_varma_garch_cpp(
    y = y_effective,
    M = M,
    spec = spec,
    model_type = model_type,
    control = ctrl
  )
  message("Model fitting complete.")
  
  # --- 5. Post-processing and Formatting Results ---
  # Align smoothed probabilities with the original time series
  smoothed_probs_aligned <- matrix(NA_real_, nrow = T_orig, ncol = M)
  colnames(smoothed_probs_aligned) <- paste0("State", 1:M)
  
  # The C++ output is aligned with the effective (differenced) data.
  # We need to pad it to match the original data's length.
  padding <- T_orig - T_eff
  if (padding > 0) {
    smoothed_probs_aligned[(padding + 1):T_orig, ] <- cpp_results$smoothed_probabilities[1:T_eff, ]
  } else {
    smoothed_probs_aligned <- cpp_results$smoothed_probabilities
  }
  
  
  # Calculate AIC/BIC
  # Count number of estimated parameters
  num_mean_pars <- sum(unlist(lapply(cpp_results$model_fits, function(fit) length(fit$arma_pars %||% fit$var_pars))))
  num_garch_pars <- sum(unlist(lapply(cpp_results$model_fits, function(fit) length(fit$garch_pars))))
  num_trans_pars <- M * (M - 1) # M*(M-1) free parameters in transition matrix
  num_params <- num_mean_pars + num_garch_pars + num_trans_pars
  
  aic <- -2 * cpp_results$log_likelihood + 2 * num_params
  bic <- -2 * cpp_results$log_likelihood + log(T_eff) * num_params
  
  # Assemble the final, user-friendly object
  result <- list(
    model_fits = cpp_results$model_fits,
    P = cpp_results$P,
    log_likelihood = cpp_results$log_likelihood,
    smoothed_probabilities = smoothed_probs_aligned,
    aic = aic,
    bic = bic,
    d = d,
    y = y_orig,
    call = match.call(),
    convergence = cpp_results$convergence,
    warnings = cpp_results$warnings
  )
  
  class(result) <- "msm.fit"
  return(result)
}

# Helper function for parameter counting
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}
