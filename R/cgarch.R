## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## Copula GARCH Implementation for MS-VARMA-GARCH Bootstrap (tsbs)
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
##
## This file implements the Copula GARCH model (cgarch_modelspec) as an 
## alternative to the DCC model for the tsbs package's Markov-Switching 
## VARMA-GARCH framework.
##
## OVERVIEW
## - - - - 
## Copula GARCH models separate marginal distributions from the dependence
## structure, providing greater flexibility than DCC for modeling multivariate
## financial time series. This implementation supports:
##
##   - Three PIT transformation methods: parametric, empirical, SPD
##   - Two copula families: Gaussian (MVN) and Student-t (MVT)
##   - Three dynamics types: constant, DCC, ADCC (asymmetric)
##   - Analytical gradients for DCC(1,1) optimization
##   - Full integration with MS-VARMA-GARCH bootstrap framework
##
## KEY DIFFERENCES FROM DCC
## 
## | Aspect            | DCC                    | Copula GARCH           |
## |- - - - - - - - - -|- - - - - - - - - - - - |- - - - - - - - - - - - |
## | Marginals         | Normal assumed         | Flexible (via PIT)     |
## | Tail dependence   | Limited                | MVT copula captures    |
## | Transformation    | None                   | PIT to uniform         |
## | Distribution      | MVN or MVT on z        | Copula on U → Z        |
##
## USAGE IN tsbs()
## - - - - - - - -
## To use CGARCH instead of DCC, specify in your model:
##
##   spec <- list(
##     list(
##       var_order = 1,
##       garch_spec_fun = "cgarch_modelspec",    # Key setting
##       distribution = "mvn",                   # or "mvt"
##       garch_spec_args = list(
##         dcc_order = c(1, 1),
##         dynamics = "dcc",                     # or "adcc", "constant"
##         transformation = "parametric",        # or "empirical", "spd"
##         copula = "mvn",                       # or "mvt"
##         garch_model = list(univariate = list(...))
##       ),
##       start_pars = list(...)
##     )
##   )
##
## FILE STRUCTURE
## - - - - - - - 
## PART 1:  Copula GARCH Weighted Estimation (main entry point)
##          - estimate_garch_weighted_cgarch()
##
## PART 2:  Copula Parameters Weighted Estimation (Stage 2)
##          - estimate_copula_parameters_weighted()
##
## PART 3:  Helper Functions for Copula Computations
##          - compute_pit_transform()
##          - fit_spd_transform()
##          - compute_spd_manual()
##          - compute_copula_residuals()
##
## PART 4:  Copula Log-Likelihood Functions
##          - copula_nll()
##          - compute_constant_copula_likelihood()
##
## PART 4b: Analytical Gradient Implementation
##          - copula_gradient()
##          - copula_gradient_original()
##          - grad_nll_wrt_R_copula_mvn/mvt()
##          - grad_R_to_Q_copula()
##          - copula_gradient_shape()
##
## PART 4c: ADCC (Asymmetric DCC) Support
##          - adcc_recursion()
##          - adcc_stationarity()
##          - adcc_copula_nll()
##          - adcc_copula_gradient()
##          - estimate_adcc_copula()
##
## DEPENDENCIES
## - - - - - - 
## Requires: tsgarch, tsmarch, tsmethods, xts, data.table
## Optional: tsdistributions (for full SPD support)
##
## The following functions from dcc_gradient.R are used:
##   - dcc11_to_unconstrained(), dcc11_from_unconstrained()
##   - dcc11_reparam_jacobian()
##   - dcc_recursion()
##   - dcc11_recursion_with_grad()
##   - is_dcc11(), get_dcc_order()
##
## SEE ALSO
## - - - - 
## - tsbs(): Main bootstrap function (set garch_spec_fun = "cgarch_modelspec")
## - dcc_gradient.R: DCC gradient and recursion functions
## - ms-varma-garch_helper_functions.R: Integration with MS framework
## - vignette("cgarch_vs_dcc"): Comparison of DCC and CGARCH models
## - vignette("Diagnostics"): Diagnostic system documentation
##
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


#### ______________________________________________________________________ ####
#### PART 1: Copula GARCH Weighted Estimation                               ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


#' @title Copula GARCH Weighted Estimation (Two-Stage)
#' @description Implements weighted MLE for Copula GARCH models.
#'   Stage 1: Estimate univariate GARCH parameters for each series
#'   Stage 2: Estimate copula dependence parameters (DCC dynamics + shape)
#'
#' @details
#' The Copula GARCH model differs from DCC in how it handles the dependence
#' structure:
#' 
#' 1. **Marginal Transformation**: Standardized residuals are transformed to 
#'    uniform \[0,1\] margins via probability integral transform (PIT). The 
#'    transformation can be:
#'    - "parametric": Uses the estimated univariate distribution's CDF
#'    - "empirical": Uses the empirical CDF
#'    - "spd": Uses semi-parametric distribution (see \code{\link{fit_spd_transform}})
#'
#' 2. **Copula Specification**: The uniform margins are then transformed
#'    according to the copula distribution:
#'    - "mvn": Multivariate Normal copula (Gaussian copula)
#'    - "mvt": Multivariate Student-t copula
#'
#' 3. **Correlation Dynamics**: Same as DCC - can be "constant", "dcc", or "adcc"
#'    For ADCC, see \code{\link{adcc_recursion}}.
#'
#' This function is called by the M-step of the EM algorithm in
#' \code{\link{fit_ms_varma_garch}} when \code{garch_spec_fun = "cgarch_modelspec"}.
#'
#' @param residuals Matrix of residuals (T x k)
#' @param weights Vector of observation weights (length T)
#' @param spec Model specification list containing:
#'   - garch_spec_fun: "cgarch_modelspec"
#'   - garch_spec_args: list with dcc_order, dynamics, transformation, copula, 
#'     garch_model (univariate specs)
#'   - start_pars: starting parameter values
#'   - distribution: copula distribution ("mvn" or "mvt")
#' @param diagnostics Optional diagnostics collector object
#' @param iteration Current EM iteration number
#' @param state Current regime state index
#' @param verbose Logical; print diagnostic information
#' @param dcc_threshold Threshold for DCC degeneracy detection (default 0.02)
#' @param dcc_criterion Selection criterion for constant vs dynamic ("bic", "aic", "threshold")
#' @param force_constant Logical; if TRUE, skip dynamic estimation
#'
#' @return List with:
#'   - coefficients: list with garch_pars, dcc_pars, dist_pars, correlation_type
#'   - warnings: list of warning messages
#'   - diagnostics: updated diagnostics object
#'
#' @seealso 
#' \itemize{
#'   \item \code{\link{tsbs}}: Main bootstrap function (use \code{garch_spec_fun = "cgarch_modelspec"})
#'   \item \code{\link{estimate_copula_parameters_weighted}}: Stage 2 estimation
#'   \item \code{\link{compute_pit_transform}}: PIT transformation methods
#'   \item \code{\link{copula_nll}}: Copula log-likelihood
#'   \item \code{\link{adcc_recursion}}: ADCC dynamics
#'   \item \code{\link{estimate_garch_weighted_dcc}}: Alternative DCC estimator
#' }
#'
#' @keywords internal
estimate_garch_weighted_cgarch <- function(
    residuals, 
    weights, 
    spec, 
    diagnostics = NULL, 
    iteration = NULL, 
    state = NULL,
    verbose = FALSE,
    dcc_threshold = 0.02,
    dcc_criterion = "bic",
    force_constant = FALSE
) {
  
  k <- ncol(residuals)
  T_obs <- nrow(residuals)
  
  ## Adjust weights to match residuals length
  if (length(weights) > T_obs) {
    n_to_remove <- length(weights) - T_obs
    w_target <- weights[(n_to_remove + 1):length(weights)]
  } else if (length(weights) < T_obs) {
    stop("Weights vector is shorter than residuals - this should not happen")
  } else {
    w_target <- weights
  }
  
  ## === STAGE 1: Estimate Univariate GARCH for Each Series ===
  garch_pars_list <- list()
  warnings_stage1 <- list()
  
  for (i in 1:k) {
    series_residuals <- residuals[, i]
    series_spec <- spec$garch_spec_args$garch_model$univariate[[i]]
    
    uni_spec <- list(
      garch_model = series_spec$model,
      garch_order = series_spec$garch_order,
      distribution = series_spec$distribution,
      start_pars = list(
        garch_pars = spec$start_pars$garch_pars[[i]],
        dist_pars = NULL
      )
    )
    
    uni_result <- estimate_garch_weighted_univariate(
      residuals = series_residuals, 
      weights = w_target, 
      spec = uni_spec,
      verbose = verbose
    )
    garch_pars_list[[i]] <- uni_result$coefficients
    warnings_stage1 <- c(warnings_stage1, uni_result$warnings)
  }
  
  ## === CHECK FOR OMEGA BOUNDARY CONDITIONS ===
  omega_boundary_threshold <- 1e-8
  
  for (i in 1:k) {
    omega_est <- garch_pars_list[[i]]$omega
    
    if (!is.null(omega_est) && omega_est < omega_boundary_threshold) {
      warning(sprintf(
        "\nState %s, Series %d: GARCH omega near boundary (%.2e < %.2e). Volatility dynamics may be degenerate.\n",
        if (!is.null(state)) state else "?", i, omega_est, omega_boundary_threshold
      ))
      
      if (!is.null(diagnostics) && !is.null(iteration) && !is.null(state)) {
        diagnostics <- add_boundary_event(
          diagnostics,
          iteration = iteration,
          state = state,
          parameter_name = sprintf("omega_series_%d", i),
          value = omega_est,
          boundary_type = "lower",
          action_taken = "warning_issued"
        )
      }
    }
  }
  
  ## === STAGE 2: Copula Parameters ===
  
  dcc_start_pars <- spec$start_pars$dcc_pars
  dist_start_pars <- spec$start_pars$dist_pars
  transformation <- spec$garch_spec_args$transformation %||% "parametric"
  copula_dist <- spec$garch_spec_args$copula %||% spec$distribution %||% "mvn"
  
  ## === SHORT-CIRCUIT: If forced to constant, skip dynamic estimation ===
  if (force_constant) {
    if (verbose) {
      cat(sprintf("\n=== State %d: Forced to CONSTANT correlation (Copula) ===\n", state))
    }
    
    return(list(
      coefficients = list(
        garch_pars = garch_pars_list,
        dcc_pars = list(),
        dist_pars = dist_start_pars,
        correlation_type = "constant",
        degeneracy_reason = "forced_constant"
      ),
      warnings = warnings_stage1,
      diagnostics = diagnostics
    ))
  }
  
  ## === HANDLE MISSING DCC PARAMETERS (Constant Copula) ===
  if (is.null(dcc_start_pars) || length(dcc_start_pars) == 0) {
    return(list(
      coefficients = list(
        garch_pars = garch_pars_list,
        dcc_pars = list(),
        dist_pars = dist_start_pars,
        correlation_type = "constant"
      ),
      warnings = warnings_stage1,
      diagnostics = diagnostics
    ))
  }
  
  ## === ESTIMATE DYNAMIC COPULA MODEL ===
  copula_result <- estimate_copula_parameters_weighted(
    residuals = residuals,
    weights = w_target,
    garch_pars = garch_pars_list,
    dcc_start_pars = dcc_start_pars,
    dist_start_pars = dist_start_pars,
    spec = spec,
    transformation = transformation,
    copula_dist = copula_dist,
    diagnostics = diagnostics,
    iteration = iteration,
    state = state,
    verbose = verbose
  )
  
  ## Update diagnostics
  if (!is.null(copula_result$diagnostics)) {
    diagnostics <- copula_result$diagnostics
  }
  
  ## GET THE WEIGHTED LL
  ll_dynamic <- copula_result$weighted_ll
  
  if (is.null(ll_dynamic)) {
    stop("Internal error: weighted_ll not returned from Copula estimation")
  }
  
  ## === DECIDE: Dynamic vs Constant ===
  
  alpha_params <- copula_result$dcc_pars[grepl("alpha", names(copula_result$dcc_pars))]
  beta_params <- copula_result$dcc_pars[grepl("beta", names(copula_result$dcc_pars))]
  
  ## Check for boundary conditions
  dcc_boundary_threshold <- 1e-4
  
  near_boundary <- !is.null(alpha_params) && 
    length(alpha_params) > 0 && 
    any(unlist(alpha_params) < dcc_threshold)
  
  if (!near_boundary) {
    ## Parameters clearly away from boundary - use dynamic
    diagnostics <- warn_dcc_boundary_params(
      alpha_params, beta_params, dcc_boundary_threshold,
      state, iteration, diagnostics
    )
    
    return(list(
      coefficients = list(
        garch_pars = garch_pars_list,
        dcc_pars = copula_result$dcc_pars,
        dist_pars = copula_result$dist_pars,
        correlation_type = "dynamic"
      ),
      warnings = c(warnings_stage1, copula_result$warnings),
      diagnostics = diagnostics
    ))
  }
  
  ## === NEAR BOUNDARY: Apply selection criterion ===
  
  if (dcc_criterion == "threshold") {
    degeneracy_reason <- sprintf("threshold (alpha=%.6f < %.3f)", 
                                 min(unlist(alpha_params)), dcc_threshold)
    
    if (verbose) {
      cat(sprintf("\n=== COPULA DEGENERACY (State %d, Iter %d) ===\n", state, iteration))
      cat("Criterion: threshold\n")
      cat("Alpha:", paste(round(unlist(alpha_params), 6), collapse = ", "), "\n")
      cat("Decision: CONSTANT\n\n")
    }
    
    if (!is.null(diagnostics) && !is.null(iteration) && !is.null(state)) {
      for (pname in names(alpha_params)) {
        diagnostics <- add_boundary_event(
          diagnostics,
          iteration = iteration,
          state = state,
          parameter_name = pname,
          value = alpha_params[[pname]],
          boundary_type = "lower",
          action_taken = paste0("constant_correlation: ", degeneracy_reason)
        )
      }
    }
    
    return(list(
      coefficients = list(
        garch_pars = garch_pars_list,
        dcc_pars = list(),
        dist_pars = copula_result$dist_pars,
        correlation_type = "constant",
        degeneracy_reason = degeneracy_reason
      ),
      warnings = c(warnings_stage1, copula_result$warnings),
      diagnostics = diagnostics
    ))
  }
  
  ## === AIC/BIC CRITERION: Compare models ===
  
  if (is.null(ll_dynamic) || !is.finite(ll_dynamic)) {
    warning("\nCould not retrieve dynamic likelihood for comparison. Falling back to threshold criterion.\n")
    degeneracy_reason <- "comparison_failed_fallback_to_threshold"
    use_constant <- TRUE
  } else {
    const_result <- compute_constant_copula_likelihood(
      residuals = residuals,
      weights = w_target,
      garch_pars = garch_pars_list,
      dist_pars = dist_start_pars,
      spec = spec,
      transformation = transformation,
      copula_dist = copula_dist
    )
    
    ll_constant <- const_result$weighted_ll
    
    k_dyn <- length(alpha_params) + length(beta_params)
    k_const <- 0
    
    n_eff <- (sum(w_target))^2 / sum(w_target^2)
    
    if (dcc_criterion == "aic") {
      ic_dynamic <- -2 * ll_dynamic + 2 * k_dyn
      ic_constant <- -2 * ll_constant + 2 * k_const
    } else {  # bic
      ic_dynamic <- -2 * ll_dynamic + log(n_eff) * k_dyn
      ic_constant <- -2 * ll_constant + log(n_eff) * k_const
    }
    
    use_constant <- (ic_constant < ic_dynamic)
    
    if (verbose) {
      cat(sprintf("\n=== Copula Model Selection (State %d, Iter %d) ===\n", state, iteration))
      cat(sprintf("Criterion: %s (n_eff=%.1f)\n", toupper(dcc_criterion), n_eff))
      cat(sprintf("Dynamic:  LL=%.4f, k=%d, %s=%.4f\n", 
                  ll_dynamic, k_dyn, toupper(dcc_criterion), ic_dynamic))
      cat(sprintf("Constant: LL=%.4f, k=%d, %s=%.4f\n", 
                  ll_constant, k_const, toupper(dcc_criterion), ic_constant))
      cat(sprintf("Decision: %s\n\n", ifelse(use_constant, "CONSTANT", "DYNAMIC")))
    }
    
    degeneracy_reason <- sprintf("%s selection (IC_const=%.2f < IC_dyn=%.2f)",
                                 toupper(dcc_criterion), ic_constant, ic_dynamic)
  }
  
  if (use_constant) {
    if (!is.null(diagnostics) && !is.null(iteration) && !is.null(state)) {
      for (pname in names(alpha_params)) {
        diagnostics <- add_boundary_event(
          diagnostics,
          iteration = iteration,
          state = state,
          parameter_name = pname,
          value = alpha_params[[pname]],
          boundary_type = "lower",
          action_taken = paste0("constant_correlation: ", degeneracy_reason)
        )
      }
    }
    
    return(list(
      coefficients = list(
        garch_pars = garch_pars_list,
        dcc_pars = list(),
        dist_pars = const_result$dist_pars,
        correlation_type = "constant",
        degeneracy_reason = degeneracy_reason
      ),
      warnings = c(warnings_stage1, const_result$warnings),
      diagnostics = diagnostics
    ))
  } else {
    diagnostics <- warn_dcc_boundary_params(
      alpha_params, beta_params, dcc_boundary_threshold,
      state, iteration, diagnostics
    )
    
    return(list(
      coefficients = list(
        garch_pars = garch_pars_list,
        dcc_pars = copula_result$dcc_pars,
        dist_pars = copula_result$dist_pars,
        correlation_type = "dynamic"
      ),
      warnings = c(warnings_stage1, copula_result$warnings),
      diagnostics = diagnostics
    ))
  }
}


#### ______________________________________________________________________ ####
#### PART 2: Copula Parameters Weighted Estimation                          ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


#' @title Estimate Copula Parameters with Weighted Likelihood
#' @description 
#' Estimates copula correlation parameters using weighted maximum likelihood.
#' This function handles the copula-specific transformation from standardized
#' residuals to uniform margins before applying the DCC-style correlation dynamics.
#'
#' @param residuals T x k matrix of residuals
#' @param weights T-vector of observation weights
#' @param garch_pars List of GARCH parameters per series
#' @param dcc_start_pars Named list of starting DCC parameters
#' @param dist_start_pars Named list of distribution parameters (e.g., shape)
#' @param spec Model specification list
#' @param transformation Type of PIT transformation ("parametric", "empirical", "spd")
#' @param copula_dist Copula distribution ("mvn" or "mvt")
#' @param diagnostics Diagnostics object (optional)
#' @param iteration Current EM iteration (optional)
#' @param state Current regime state (optional)
#' @param verbose Logical; print diagnostic information
#'
#' @return List with dcc_pars, dist_pars, weighted_ll, warnings, diagnostics
#' @keywords internal
estimate_copula_parameters_weighted <- function(
    residuals, 
    weights, 
    garch_pars,
    dcc_start_pars, 
    dist_start_pars, 
    spec,
    transformation = "parametric",
    copula_dist = "mvn",
    dynamics = "dcc",
    diagnostics = NULL,
    iteration = NULL,
    state = NULL,
    verbose = FALSE
) {
  
  k <- ncol(residuals)
  T_obs <- nrow(residuals)
  
  ## Get dynamics from spec if not explicitly provided
  if (missing(dynamics) && !is.null(spec$garch_spec_args$dynamics)) {
    dynamics <- spec$garch_spec_args$dynamics
  }
  
  if (length(weights) != T_obs) {
    stop(sprintf("Dimension mismatch: residuals has %d rows but weights has %d elements", 
                 T_obs, length(weights)))
  }
  
  ## 1. Get standardized residuals using Stage 1 GARCH parameters ==============
  std_residuals <- matrix(0, nrow = T_obs, ncol = k)
  sigma_list <- list()
  uni_fit_list <- list()
  
  for (i in 1:k) {
    series_residuals <- residuals[, i]
    series_spec <- spec$garch_spec_args$garch_model$univariate[[i]]
    
    uni_spec_obj <- tsgarch::garch_modelspec(
      y = xts::xts(series_residuals, order.by = Sys.Date() - (T_obs:1)),
      model = series_spec$model,
      garch_order = series_spec$garch_order,
      distribution = series_spec$distribution
    )
    
    for (par_name in names(garch_pars[[i]])) {
      if (par_name %in% uni_spec_obj$parmatrix$parameter) {
        row_idx <- which(uni_spec_obj$parmatrix$parameter == par_name)
        uni_spec_obj$parmatrix[row_idx, "value"] <- garch_pars[[i]][[par_name]]
      }
    }
    
    uni_fit <- tsmethods::tsfilter(uni_spec_obj)
    sigma_vec <- as.numeric(uni_fit$sigma)
    std_residuals[, i] <- series_residuals / sigma_vec
    sigma_list[[i]] <- sigma_vec
    uni_fit_list[[i]] <- uni_fit
  }
  
  ## 2. Transform to Uniform Margins (Copula-specific) =========================
  ## This is the key difference from DCC - we apply PIT transformation
  
  u_matrix <- compute_pit_transform(
    std_residuals = std_residuals,
    uni_fit_list = uni_fit_list,
    transformation = transformation,
    copula_dist = copula_dist,
    dist_pars = dist_start_pars
  )
  
  ## Bound u away from 0 and 1 for numerical stability
  u_matrix[u_matrix < 3.330669e-16] <- 2.220446e-16
  u_matrix[u_matrix > 0.99999] <- 0.99999
  
  ## 3. Transform u to copula residuals Z =====================================
  z_matrix <- compute_copula_residuals(
    u_matrix = u_matrix,
    copula_dist = copula_dist,
    dist_pars = dist_start_pars
  )
  
  ## 4. Handle empty parameter case ============================================
  all_stage2_pars <- c(dcc_start_pars, dist_start_pars)
  
  if (length(all_stage2_pars) == 0) {
    return(list(
      dcc_pars = list(), 
      dist_pars = list(), 
      weighted_ll = 0,
      warnings = list(),
      diagnostics = diagnostics
    ))
  }
  
  ## 5. Compute Qbar (weighted unconditional covariance of Z) ==================
  Qbar <- tryCatch({
    stats::cov.wt(z_matrix, wt = weights, method = "ML")$cov
  }, error = function(e) {
    if (verbose) cat("ERROR in cov.wt:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(Qbar)) {
    return(list(
      dcc_pars = dcc_start_pars,
      dist_pars = dist_start_pars,
      weighted_ll = -1e10,
      warnings = list("Weighted covariance calculation failed"),
      diagnostics = diagnostics
    ))
  }
  
  ## Ensure positive definite
  eig <- eigen(Qbar, symmetric = TRUE)
  if (any(eig$values < 1e-8)) {
    Qbar <- Qbar + diag(1e-6, k)
  }
  
  ## 6. Determine optimization strategy based on DCC order =====================
  ## 6. Determine optimization strategy based on dynamics and DCC order ========
  is_adcc <- (dynamics == "adcc")
  use_analytical_gradient <- is_dcc11(dcc_start_pars) && !is_adcc
  
  if (verbose) {
    if (is_adcc) {
      cat("Copula GARCH ADCC detected: Using numerical gradient\n")
    } else if (use_analytical_gradient) {
      cat("Copula GARCH DCC(1,1) detected: Using analytical gradient with reparameterization\n")
    } else {
      order <- get_dcc_order(dcc_start_pars)
      cat(sprintf("Copula GARCH DCC(%d,%d) detected: Using finite differences\n",
                  order["p"], order["q"]))
    }
  }
  
  warnings_list <- list()
  
  ## 7. Optimization ===========================================================
  
  if (is_adcc) {
    ## ADCC optimization
    adcc_result <- estimate_adcc_copula(
      z_matrix = z_matrix,
      weights = weights,
      Qbar = Qbar,
      copula_dist = copula_dist,
      start_pars = NULL  # Use defaults
    )
    
    ## Convert to output format
    dcc_pars_final <- list(
      alpha_1 = adcc_result$alpha,
      gamma_1 = adcc_result$gamma,
      beta_1 = adcc_result$beta
    )
    
    if (copula_dist == "mvt") {
      dist_pars_final <- list(shape = adcc_result$shape)
    } else {
      dist_pars_final <- list()
    }
    
    final_nll <- adcc_result$nll
    
    if (adcc_result$stationarity >= 0.999) {
      warnings_list <- c(warnings_list, 
                         list("ADCC stationarity constraint near boundary"))
    }
    
    return(list(
      dcc_pars = dcc_pars_final,
      dist_pars = dist_pars_final,
      weighted_ll = -final_nll,
      warnings = warnings_list,
      diagnostics = diagnostics,
      Nbar = adcc_result$Nbar
    ))
  } else if (use_analytical_gradient) {
    ## DCC(1,1) with reparameterization
    alpha_start <- as.numeric(dcc_start_pars$alpha_1)
    beta_start <- as.numeric(dcc_start_pars$beta_1)
    
    unconstrained_start <- dcc11_to_unconstrained(alpha_start, beta_start)
    psi_start <- unconstrained_start["psi"]
    phi_start <- unconstrained_start["phi"]
    
    if (copula_dist == "mvt" && !is.null(dist_start_pars$shape)) {
      opt_params <- c(psi = psi_start, phi = phi_start, shape = dist_start_pars$shape)
    } else {
      opt_params <- c(psi = psi_start, phi = phi_start)
    }
    
    ## Define objective function for copula
    objective_fn <- function(params) {
      copula_nll(
        params = params,
        z_matrix = z_matrix,
        weights = weights,
        Qbar = Qbar,
        copula_dist = copula_dist,
        use_reparam = TRUE
      )
    }
    
    ## Define gradient function
    gradient_fn <- function(params) {
      copula_gradient(
        params = params,
        z_matrix = z_matrix,
        weights = weights,
        Qbar = Qbar,
        copula_dist = copula_dist,
        use_reparam = TRUE
      )
    }
    
    ## Set bounds
    if (copula_dist == "mvt") {
      lower_bounds <- c(-Inf, -Inf, 2.01)
      upper_bounds <- c(Inf, Inf, Inf)
    } else {
      lower_bounds <- c(-Inf, -Inf)
      upper_bounds <- c(Inf, Inf)
    }
    
    ## Run optimization
    opt_result <- withCallingHandlers({
      stats::optim(
        par = opt_params,
        fn = objective_fn,
        gr = gradient_fn,
        method = "L-BFGS-B",
        lower = lower_bounds,
        upper = upper_bounds,
        control = list(maxit = 500)
      )
    }, warning = function(w) {
      warnings_list <<- c(warnings_list, list(w))
      invokeRestart("muffleWarning")
    })
    
    ## Transform results back
    final_psi <- as.numeric(opt_result$par[1])
    final_phi <- as.numeric(opt_result$par[2])
    final_orig <- dcc11_from_unconstrained(final_psi, final_phi)
    
    dcc_pars_final <- list(
      alpha_1 = as.numeric(final_orig["alpha"]),
      beta_1 = as.numeric(final_orig["beta"])
    )
    
    if (copula_dist == "mvt") {
      dist_pars_final <- list(shape = as.numeric(opt_result$par[3]))
    } else {
      dist_pars_final <- dist_start_pars
    }
    
    weighted_ll <- -opt_result$value
    
  } else {
    ## Higher-order DCC with finite differences
    
    param_vector <- unlist(all_stage2_pars)
    
    weighted_copula_loglik <- function(params, z_mat, w, copula_dist, Qbar, dcc_start_pars, 
                                       dist_start_pars, verbose = FALSE) {
      
      param_list <- as.list(params)
      names(param_list) <- names(c(dcc_start_pars, dist_start_pars))
      
      dcc_param_names <- names(dcc_start_pars)
      dist_param_names <- names(dist_start_pars)
      
      dcc_params_current <- param_list[dcc_param_names]
      dist_params_current <- param_list[dist_param_names]
      
      pers <- compute_dcc_persistence(dcc_params_current)
      
      if (pers$persistence >= 1 || any(pers$alphas < 0) || any(pers$betas < 0)) {
        return(1e10)
      }
      
      recursion <- dcc_recursion(z_mat, Qbar, pers$alphas, pers$betas)
      
      if (!recursion$success) {
        return(1e10)
      }
      
      R <- recursion$R
      T_eff <- nrow(z_mat)
      ll_vec <- numeric(T_eff)
      
      if (copula_dist == "mvn") {
        for (t in 1:T_eff) {
          R_t <- R[,,t]
          if (any(!is.finite(R_t))) {
            ll_vec[t] <- -1e10
            next
          }
          det_R <- det(R_t)
          if (det_R <= 0) {
            ll_vec[t] <- -1e10
            next
          }
          R_inv <- tryCatch(solve(R_t), error = function(e) NULL)
          if (is.null(R_inv)) {
            ll_vec[t] <- -1e10
            next
          }
          z_t <- z_mat[t, ]
          mahal <- as.numeric(t(z_t) %*% R_inv %*% z_t)
          ll_vec[t] <- -0.5 * (log(det_R) + mahal - as.numeric(t(z_t) %*% z_t))
        }
      } else {
        ## MVT copula
        shape <- dist_params_current$shape
        if (is.null(shape) || shape <= 2) {
          return(1e10)
        }
        k <- ncol(z_mat)
        const_term <- lgamma(0.5 * (k + shape)) - lgamma(0.5 * shape) - 
          0.5 * k * log(pi * (shape - 2))
        
        for (t in 1:T_eff) {
          R_t <- R[,,t]
          if (any(!is.finite(R_t))) {
            ll_vec[t] <- -1e10
            next
          }
          det_R <- det(R_t)
          if (det_R <= 0) {
            ll_vec[t] <- -1e10
            next
          }
          R_inv <- tryCatch(solve(R_t), error = function(e) NULL)
          if (is.null(R_inv)) {
            ll_vec[t] <- -1e10
            next
          }
          z_t <- z_mat[t, ]
          mahal <- as.numeric(t(z_t) %*% R_inv %*% z_t)
          
          ## MVT log-density for copula
          ll_vec[t] <- const_term - 0.5 * log(det_R) - 
            0.5 * (shape + k) * log(1 + mahal / (shape - 2))
          
          ## Subtract marginal t densities (already in z_t)
          for (j in 1:k) {
            ll_vec[t] <- ll_vec[t] - dt_scaled_log(z_t[j], shape)
          }
        }
      }
      
      ll_vec[!is.finite(ll_vec)] <- -1e10
      nll <- -sum(w * ll_vec)
      
      return(nll)
    }
    
    ## Get bounds
    lower_bounds <- sapply(names(all_stage2_pars), function(pn) {
      if (pn == "shape") return(2.01)
      if (grepl("alpha|beta|gamma", pn)) return(1e-8)
      return(-Inf)
    })
    
    upper_bounds <- sapply(names(all_stage2_pars), function(pn) {
      if (grepl("alpha|beta|gamma", pn)) return(0.999)
      if (pn == "shape") return(100)
      return(Inf)
    })
    
    opt_result <- withCallingHandlers({
      stats::optim(
        par = param_vector,
        fn = weighted_copula_loglik,
        method = "L-BFGS-B",
        lower = lower_bounds,
        upper = upper_bounds,
        control = list(ndeps = rep(1e-12, length(param_vector))),
        z_mat = z_matrix,
        w = weights,
        copula_dist = copula_dist,
        Qbar = Qbar,
        dcc_start_pars = dcc_start_pars,
        dist_start_pars = dist_start_pars,
        verbose = verbose
      )
    }, warning = function(w) {
      warnings_list <<- c(warnings_list, list(w))
      invokeRestart("muffleWarning")
    })
    
    final_params <- as.list(opt_result$par)
    names(final_params) <- names(all_stage2_pars)
    
    dcc_pars_final <- final_params[names(dcc_start_pars)]
    dist_pars_final <- final_params[names(dist_start_pars)]
    if (length(dist_pars_final) == 0) dist_pars_final <- dist_start_pars
    
    weighted_ll <- -opt_result$value
  }
  
  if (verbose) {
    cat("\n=== COPULA M-STEP DIAGNOSTIC ===\n")
    cat("Estimated DCC params:", paste(names(dcc_pars_final), "=", 
                                       round(unlist(dcc_pars_final), 4), collapse = ", "), "\n")
    if (!is.null(dist_pars_final$shape)) {
      cat("Estimated shape:", round(dist_pars_final$shape, 2), "\n")
    }
    cat("Weighted log-likelihood:", round(weighted_ll, 4), "\n")
  }
  
  return(list(
    dcc_pars = dcc_pars_final,
    dist_pars = dist_pars_final,
    weighted_ll = weighted_ll,
    warnings = warnings_list,
    diagnostics = diagnostics
  ))
}


#### ______________________________________________________________________ ####
#### PART 3: Helper Functions for Copula Computations                       ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


#' @title Compute Probability Integral Transform (PIT)
#' @description Transforms standardized residuals to uniform \[0,1\] margins
#' @param std_residuals Matrix of standardized residuals (T x k)
#' @param uni_fit_list List of fitted univariate GARCH models
#' @param transformation Type: "parametric", "empirical", or "spd"
#' @param copula_dist Copula distribution ("mvn" or "mvt")
#' @param dist_pars Distribution parameters
#' @return Matrix of uniform margins (T x k)
#' @keywords internal
compute_pit_transform <- function(
    std_residuals,
    uni_fit_list,
    transformation = "parametric",
    copula_dist = "mvn",
    dist_pars = NULL
) {
  T_obs <- nrow(std_residuals)
  k <- ncol(std_residuals)
  u_matrix <- matrix(0, nrow = T_obs, ncol = k)
  
  if (transformation == "parametric") {
    ## Use the fitted distribution's CDF
    for (i in 1:k) {
      uni_fit <- uni_fit_list[[i]]
      dist_name <- uni_fit$spec$distribution
      
      ## Get distribution parameters from the fitted model
      skew <- uni_fit$parmatrix[parameter == "skew"]$value
      shape <- uni_fit$parmatrix[parameter == "shape"]$value
      lambda <- uni_fit$parmatrix[parameter == "lambda"]$value
      
      if (is.null(skew) || length(skew) == 0) skew <- 1
      if (is.null(shape) || length(shape) == 0) shape <- 4
      if (is.null(lambda) || length(lambda) == 0) lambda <- 0
      
      u_matrix[, i] <- tsdistributions::pdist(
        distribution = dist_name,
        q = std_residuals[, i],
        mu = 0,
        sigma = 1,
        skew = skew,
        shape = shape,
        lambda = lambda
      )
    }
  } else if (transformation == "empirical") {
    ## Use empirical CDF
    for (i in 1:k) {
      ecdf_fn <- ecdf(std_residuals[, i])
      u_matrix[, i] <- ecdf_fn(std_residuals[, i])
    }
    ## Adjust to avoid exact 0 and 1
    u_matrix <- u_matrix * (T_obs / (T_obs + 1))
  } else if (transformation == "spd") {
    ## Semi-Parametric Distribution (SPD)
    ## Combines parametric tails with empirical/kernel-smoothed center
    ## Requires tsdistributions package
    
    if (!requireNamespace("tsdistributions", quietly = TRUE)) {
      warning("tsdistributions package not available for SPD transformation, using empirical.\n")
      for (i in 1:k) {
        ecdf_fn <- ecdf(std_residuals[, i])
        u_matrix[, i] <- ecdf_fn(std_residuals[, i])
      }
      u_matrix <- u_matrix * (T_obs / (T_obs + 1))
    } else {
      ## Use SPD transformation for each series
      for (i in 1:k) {
        spd_result <- fit_spd_transform(
          z = std_residuals[, i],
          uni_fit = uni_fit_list[[i]],
          dist_pars = dist_pars
        )
        u_matrix[, i] <- spd_result$u
      }
    }
  }
  
  return(u_matrix)
}


#' @title Fit SPD Transformation for a Single Series
#' @description Fits a Semi-Parametric Distribution and computes PIT values.
#'   SPD combines:
#'   - Parametric tails from the fitted univariate distribution
#'   - Kernel-smoothed empirical distribution in the center
#' @param z Standardized residuals for one series
#' @param uni_fit Univariate GARCH fit object (optional, for parametric tails)
#' @param dist_pars Distribution parameters (optional)
#' @param lower_threshold Lower quantile for parametric tail (default 0.1)
#' @param upper_threshold Upper quantile for parametric tail (default 0.9)
#' @return List with u (uniform values) and spd_model (fitted SPD object)
#' @keywords internal
fit_spd_transform <- function(
    z,
    uni_fit = NULL,
    dist_pars = NULL,
    lower_threshold = 0.1,
    upper_threshold = 0.9
) {
  n <- length(z)
  
  ## Try to use tsdistributions::spd_modelspec if available
  if (requireNamespace("tsdistributions", quietly = TRUE)) {
    tryCatch({
      ## Fit SPD model
      ## SPD uses kernel density estimation in the interior
      ## and GPD (Generalized Pareto) for the tails
      spd_spec <- tsdistributions::spd_modelspec(
        y = z,
        lower = lower_threshold,
        upper = upper_threshold
      )
      
      spd_fit <- tsdistributions::estimate(spd_spec)
      
      ## Get PIT values using the fitted SPD
      u <- tsdistributions::pspd(z, object = spd_fit)
      
      ## Bound away from 0 and 1
      u <- pmax(pmin(u, 1 - 1e-10), 1e-10)
      
      return(list(u = u, spd_model = spd_fit))
    }, error = function(e) {
      ## Fall back to manual SPD implementation
      warning(sprintf("SPD fitting failed: %s. Using manual SPD.\n", e$message))
    })
  }
  
  ## Manual SPD implementation as fallback
  ## This implements a simplified version of SPD:
  ## - Uses kernel density estimation (KDE) for the interior
  ## - Uses empirical quantiles for the tails
  
  u <- compute_spd_manual(z, lower_threshold, upper_threshold)
  
  return(list(u = u, spd_model = NULL))
}


#' @title Manual SPD Computation
#' @description Computes SPD transformation manually without tsdistributions.
#'   Uses kernel density estimation for interior and empirical tails.
#' @param z Standardized residuals
#' @param lower_threshold Lower quantile threshold
#' @param upper_threshold Upper quantile threshold
#' @return Vector of uniform values
#' @keywords internal
compute_spd_manual <- function(z, lower_threshold = 0.1, upper_threshold = 0.9) {
  n <- length(z)
  
  ## Get threshold values
  q_lower <- quantile(z, lower_threshold)
  q_upper <- quantile(z, upper_threshold)
  
  ## Identify regions
  idx_lower <- z < q_lower
  idx_upper <- z > q_upper
  idx_middle <- !idx_lower & !idx_upper
  
  u <- numeric(n)
  
  ## Lower tail: Use empirical CDF scaled to [0, lower_threshold]
  if (any(idx_lower)) {
    z_lower <- z[idx_lower]
    ## Rank within lower tail
    ranks_lower <- rank(z_lower, ties.method = "average")
    n_lower <- sum(idx_lower)
    ## Scale to [0, lower_threshold]
    u[idx_lower] <- (ranks_lower / (n_lower + 1)) * lower_threshold
  }
  
  ## Upper tail: Use empirical CDF scaled to [upper_threshold, 1]
  if (any(idx_upper)) {
    z_upper <- z[idx_upper]
    ## Rank within upper tail
    ranks_upper <- rank(z_upper, ties.method = "average")
    n_upper <- sum(idx_upper)
    ## Scale to [upper_threshold, 1]
    u[idx_upper] <- upper_threshold + (ranks_upper / (n_upper + 1)) * (1 - upper_threshold)
  }
  
  ## Middle: Use kernel density estimation
  if (any(idx_middle)) {
    z_middle <- z[idx_middle]
    
    ## Fit kernel density to middle region
    kde <- tryCatch({
      density(z_middle, bw = "SJ", n = 512)  # Sheather-Jones bandwidth
    }, error = function(e) {
      density(z_middle, n = 512)  # Fall back to default bandwidth
    })
    
    ## Create CDF from KDE by numerical integration
    kde_cdf <- approxfun(
      x = kde$x,
      y = cumsum(kde$y) / sum(kde$y),
      rule = 2  # Return boundary values for out-of-range
    )
    
    ## Get raw CDF values
    cdf_raw <- kde_cdf(z_middle)
    
    ## Scale to [lower_threshold, upper_threshold]
    u[idx_middle] <- lower_threshold + cdf_raw * (upper_threshold - lower_threshold)
  }
  
  ## Ensure values are in (0, 1)
  u <- pmax(pmin(u, 1 - 1e-10), 1e-10)
  
  return(u)
}


#' @title Compute Copula Residuals from Uniform Margins
#' @description Transforms uniform margins to copula residuals (Z) based on
#'   the copula distribution
#' @param u_matrix Matrix of uniform margins (T x k)
#' @param copula_dist Copula distribution ("mvn" or "mvt")
#' @param dist_pars Distribution parameters (shape for mvt)
#' @return Matrix of copula residuals (T x k)
#' @keywords internal
compute_copula_residuals <- function(
    u_matrix,
    copula_dist = "mvn",
    dist_pars = NULL
) {
  T_obs <- nrow(u_matrix)
  k <- ncol(u_matrix)
  z_matrix <- matrix(0, nrow = T_obs, ncol = k)
  
  if (copula_dist == "mvn") {
    ## Gaussian copula: Z = Phi^{-1}(U)
    for (i in 1:k) {
      z_matrix[, i] <- qnorm(u_matrix[, i])
    }
  } else if (copula_dist == "mvt") {
    ## Student-t copula: Z = t_{nu}^{-1}(U)
    shape <- dist_pars$shape %||% 4
    scale <- sqrt(shape / (shape - 2))
    
    for (i in 1:k) {
      z_matrix[, i] <- qt(u_matrix[, i], df = shape) / scale
    }
  }
  
  ## Handle any infinite values
  z_matrix[!is.finite(z_matrix)] <- 0
  
  return(z_matrix)
}


#' @title Copula Negative Log-Likelihood for DCC(1,1)
#' @description Computes the weighted copula negative log-likelihood
#' @param params Parameter vector (psi, phi) or (psi, phi, shape)
#' @param z_matrix Matrix of copula residuals (T x k)
#' @param weights Observation weights
#' @param Qbar Unconditional covariance matrix
#' @param copula_dist Copula distribution ("mvn" or "mvt")
#' @param use_reparam Logical; if TRUE, params are in reparameterized space
#' @return Scalar negative log-likelihood
#' @keywords internal
copula_nll <- function(
    params,
    z_matrix,
    weights,
    Qbar,
    copula_dist = "mvn",
    use_reparam = TRUE
) {
  ## Extract parameters
  if (use_reparam) {
    psi <- params[1]
    phi <- params[2]
    ab <- dcc11_from_unconstrained(psi, phi)
    alpha <- ab["alpha"]
    beta <- ab["beta"]
  } else {
    alpha <- params[1]
    beta <- params[2]
  }
  
  if (copula_dist == "mvt" && length(params) >= 3) {
    shape <- params[3]
    if (shape <= 2) return(1e10)
  } else {
    shape <- NULL
  }
  
  ## Check stationarity
  persistence <- alpha + beta
  if (persistence >= 1 || alpha < 0 || beta < 0) {
    return(1e10)
  }
  
  ## DCC recursion
  recursion <- dcc_recursion(
    std_resid = z_matrix,
    Qbar = Qbar,
    alphas = alpha,
    betas = beta,
    verbose = FALSE
  )
  
  if (!recursion$success) {
    return(1e10)
  }
  
  R <- recursion$R
  T_obs <- nrow(z_matrix)
  k <- ncol(z_matrix)
  ll_vec <- numeric(T_obs)
  
  if (copula_dist == "mvn") {
    ## Gaussian copula log-likelihood
    for (t in 1:T_obs) {
      R_t <- R[,,t]
      if (any(!is.finite(R_t))) {
        ll_vec[t] <- -1e10
        next
      }
      
      det_R <- det(R_t)
      if (det_R <= 0) {
        ll_vec[t] <- -1e10
        next
      }
      
      R_inv <- tryCatch(solve(R_t), error = function(e) NULL)
      if (is.null(R_inv)) {
        ll_vec[t] <- -1e10
        next
      }
      
      z_t <- z_matrix[t, ]
      mahal <- as.numeric(t(z_t) %*% R_inv %*% z_t)
      
      ## Copula LL = -0.5 * (log|R| + z'R^{-1}z - z'z)
      ll_vec[t] <- -0.5 * (log(det_R) + mahal - as.numeric(t(z_t) %*% z_t))
    }
  } else {
    ## Student-t copula
    if (is.null(shape)) shape <- 4
    const_term <- lgamma(0.5 * (k + shape)) - lgamma(0.5 * shape) - 
      0.5 * k * log(pi * (shape - 2))
    
    for (t in 1:T_obs) {
      R_t <- R[,,t]
      if (any(!is.finite(R_t))) {
        ll_vec[t] <- -1e10
        next
      }
      
      det_R <- det(R_t)
      if (det_R <= 0) {
        ll_vec[t] <- -1e10
        next
      }
      
      R_inv <- tryCatch(solve(R_t), error = function(e) NULL)
      if (is.null(R_inv)) {
        ll_vec[t] <- -1e10
        next
      }
      
      z_t <- z_matrix[t, ]
      mahal <- as.numeric(t(z_t) %*% R_inv %*% z_t)
      
      ll_vec[t] <- const_term - 0.5 * log(det_R) - 
        0.5 * (shape + k) * log(1 + mahal / (shape - 2))
      
      ## Subtract marginal t log-densities
      for (j in 1:k) {
        ll_vec[t] <- ll_vec[t] - dt_scaled_log(z_t[j], shape)
      }
    }
  }
  
  ll_vec[!is.finite(ll_vec)] <- -1e10
  nll <- -sum(weights * ll_vec)
  
  return(nll)
}


#### ______________________________________________________________________ ####
#### PART 4c: ADCC (Asymmetric DCC) Support                                 ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


#' @title ADCC(1,1) Recursion
#' @description Computes Q and R matrices for ADCC(1,1) model.
#'   ADCC adds asymmetric response to negative shocks via gamma parameter.
#'   
#'   The ADCC recursion is:
#'   Q_t = Ω + α*(z_{t-1}z'_{t-1}) + γ*(n_{t-1}n'_{t-1}) + β*Q_{t-1}
#'   
#'   where n_t = z_t * I(z_t < 0) (element-wise negative indicator)
#'   and Ω = (1 - α - β)*Qbar - γ*Nbar
#'   
#' @param std_resid Matrix of standardized residuals (T x k)
#' @param Qbar Unconditional covariance matrix (k x k)
#' @param alpha ADCC alpha parameter (response to shocks)
#' @param gamma ADCC gamma parameter (asymmetric response to negative shocks)
#' @param beta ADCC beta parameter (persistence)
#' @param Nbar Unconditional covariance of negative shocks (k x k), computed if NULL
#' @return List with success, Q, R, Nbar, and error info if failed
#' @keywords internal
adcc_recursion <- function(std_resid, Qbar, alpha, gamma, beta, Nbar = NULL) {
  T_obs <- nrow(std_resid)
  k <- ncol(std_resid)
  
  ## Ensure scalar parameters
  alpha <- as.numeric(alpha)[1]
  gamma <- as.numeric(gamma)[1]
  beta <- as.numeric(beta)[1]
  
  ## Compute negative indicator shocks: n_t = z_t * I(z_t < 0)
  neg_resid <- std_resid * (std_resid < 0)
  
  ## Compute Nbar (unconditional covariance of negative shocks) if not provided
  if (is.null(Nbar)) {
    Nbar <- crossprod(neg_resid) / T_obs
  }
  
  ## Compute Omega: Ω = (1 - α - β)*Qbar - γ*Nbar
  persistence <- alpha + beta
  Omega <- (1 - persistence) * Qbar - gamma * Nbar
  
  ## Initialize arrays
  Q <- array(0, dim = c(k, k, T_obs))
  R <- array(0, dim = c(k, k, T_obs))
  
  ## Initialize first observation with Qbar
  Q[,,1] <- Qbar
  Qbar_diag <- diag(Qbar)
  if (any(Qbar_diag <= 0)) {
    return(list(success = FALSE, error_type = "non_positive_Qbar_diagonal",
                error_time = 0, Q = Q, R = R, Nbar = Nbar))
  }
  Qbar_diag_inv_sqrt <- diag(1/sqrt(Qbar_diag), k)
  R[,,1] <- Qbar_diag_inv_sqrt %*% Qbar %*% Qbar_diag_inv_sqrt
  
  ## Main recursion for t >= 2
  for (t in 2:T_obs) {
    z_lag <- std_resid[t - 1, , drop = FALSE]
    n_lag <- neg_resid[t - 1, , drop = FALSE]
    
    Q_t <- Omega + 
      alpha * (t(z_lag) %*% z_lag) + 
      gamma * (t(n_lag) %*% n_lag) + 
      beta * Q[,,t - 1]
    
    if (any(!is.finite(Q_t))) {
      return(list(success = FALSE, error_type = "non_finite_Q",
                  error_time = t, Q = Q, R = R, Nbar = Nbar))
    }
    
    Q_diag <- diag(Q_t)
    if (any(Q_diag <= 0)) {
      return(list(success = FALSE, error_type = "non_positive_Q_diagonal",
                  error_time = t, Q = Q, R = R, Nbar = Nbar))
    }
    
    Q[,,t] <- Q_t
    Q_diag_inv_sqrt <- diag(1/sqrt(Q_diag), k)
    R[,,t] <- Q_diag_inv_sqrt %*% Q_t %*% Q_diag_inv_sqrt
  }
  
  return(list(success = TRUE, Q = Q, R = R, Nbar = Nbar))
}


#' @title ADCC Stationarity Constraint
#' @description Computes the ADCC stationarity constraint value.
#'   For ADCC: α + β + δ*γ < 1, where δ depends on the data.
#'   A simplified constraint is: α + β + 0.5*γ < 1
#' @param alpha ADCC alpha
#' @param gamma ADCC gamma  
#' @param beta ADCC beta
#' @param delta Asymmetry scaling factor (default 0.5)
#' @return Stationarity measure (should be < 1 for stationarity)
#' @keywords internal
adcc_stationarity <- function(alpha, gamma, beta, delta = 0.5) {
  alpha + beta + delta * gamma
}


#' @title ADCC Copula Negative Log-Likelihood
#' @description Computes the weighted copula NLL for ADCC model.
#' @param params Parameter vector: (alpha, gamma, beta) or with shape for MVT
#' @param z_matrix Matrix of copula residuals (T x k)
#' @param weights Observation weights
#' @param Qbar Unconditional covariance matrix
#' @param Nbar Unconditional covariance of negative shocks (computed if NULL)
#' @param copula_dist Copula distribution ("mvn" or "mvt")
#' @return Scalar negative log-likelihood
#' @keywords internal
adcc_copula_nll <- function(
    params,
    z_matrix,
    weights,
    Qbar,
    Nbar = NULL,
    copula_dist = "mvn"
) {
  ## Extract parameters
  alpha <- params[1]
  gamma <- params[2]
  beta <- params[3]
  
  if (copula_dist == "mvt" && length(params) >= 4) {
    shape <- params[4]
    if (shape <= 2) return(1e10)
  } else {
    shape <- NULL
  }
  
  ## Check stationarity (simplified constraint)
  if (adcc_stationarity(alpha, gamma, beta) >= 1 || 
      alpha < 0 || gamma < 0 || beta < 0) {
    return(1e10)
  }
  
  ## ADCC recursion
  recursion <- adcc_recursion(
    std_resid = z_matrix,
    Qbar = Qbar,
    alpha = alpha,
    gamma = gamma,
    beta = beta,
    Nbar = Nbar
  )
  
  if (!recursion$success) {
    return(1e10)
  }
  
  R <- recursion$R
  T_obs <- nrow(z_matrix)
  k <- ncol(z_matrix)
  ll_vec <- numeric(T_obs)
  
  if (copula_dist == "mvn") {
    ## Gaussian copula log-likelihood
    for (t in 1:T_obs) {
      R_t <- R[,,t]
      if (any(!is.finite(R_t))) {
        ll_vec[t] <- -1e10
        next
      }
      
      det_R <- det(R_t)
      if (det_R <= 0) {
        ll_vec[t] <- -1e10
        next
      }
      
      R_inv <- tryCatch(solve(R_t), error = function(e) NULL)
      if (is.null(R_inv)) {
        ll_vec[t] <- -1e10
        next
      }
      
      z_t <- z_matrix[t, ]
      mahal <- as.numeric(t(z_t) %*% R_inv %*% z_t)
      
      ## Copula LL = -0.5 * (log|R| + z'R^{-1}z - z'z)
      ll_vec[t] <- -0.5 * (log(det_R) + mahal - as.numeric(t(z_t) %*% z_t))
    }
  } else {
    ## Student-t copula
    if (is.null(shape)) shape <- 4
    const_term <- lgamma(0.5 * (k + shape)) - lgamma(0.5 * shape) - 
      0.5 * k * log(pi * (shape - 2))
    
    for (t in 1:T_obs) {
      R_t <- R[,,t]
      if (any(!is.finite(R_t))) {
        ll_vec[t] <- -1e10
        next
      }
      
      det_R <- det(R_t)
      if (det_R <= 0) {
        ll_vec[t] <- -1e10
        next
      }
      
      R_inv <- tryCatch(solve(R_t), error = function(e) NULL)
      if (is.null(R_inv)) {
        ll_vec[t] <- -1e10
        next
      }
      
      z_t <- z_matrix[t, ]
      mahal <- as.numeric(t(z_t) %*% R_inv %*% z_t)
      
      ll_vec[t] <- const_term - 0.5 * log(det_R) - 
        0.5 * (shape + k) * log(1 + mahal / (shape - 2))
      
      ## Subtract marginal t log-densities
      for (j in 1:k) {
        ll_vec[t] <- ll_vec[t] - dt_scaled_log(z_t[j], shape)
      }
    }
  }
  
  ll_vec[!is.finite(ll_vec)] <- -1e10
  nll <- -sum(weights * ll_vec)
  
  return(nll)
}


#' @title ADCC Copula Gradient (Numerical)
#' @description Computes the gradient of ADCC copula NLL using numerical differentiation.
#' @param params Parameter vector: (alpha, gamma, beta) or with shape for MVT
#' @param z_matrix Matrix of copula residuals
#' @param weights Observation weights
#' @param Qbar Unconditional covariance matrix
#' @param Nbar Unconditional covariance of negative shocks
#' @param copula_dist Copula distribution
#' @return Gradient vector
#' @keywords internal
adcc_copula_gradient <- function(
    params,
    z_matrix,
    weights,
    Qbar,
    Nbar = NULL,
    copula_dist = "mvn"
) {
  ## Use numerical gradient for ADCC (analytical is more complex due to gamma)
  eps <- 1e-6
  n_params <- length(params)
  grad <- numeric(n_params)
  
  f0 <- adcc_copula_nll(
    params = params, 
    z_matrix = z_matrix, 
    weights = weights, 
    Qbar = Qbar, 
    Nbar = Nbar, 
    copula_dist = copula_dist
  )
  
  for (i in 1:n_params) {
    params_plus <- params
    params_plus[i] <- params_plus[i] + eps
    f_plus <- adcc_copula_nll(
      params = params, 
      z_matrix = z_matrix, 
      weights = weights, 
      Qbar = Qbar, 
      Nbar = Nbar, 
      copula_dist = copula_dist
    )
    grad[i] <- (f_plus - f0) / eps
  }
  
  return(grad)
}


#' @title Estimate ADCC Copula Parameters
#' @description Estimates ADCC parameters (alpha, gamma, beta) for copula model.
#' @param z_matrix Matrix of copula residuals (T x k)
#' @param weights Observation weights
#' @param Qbar Unconditional covariance matrix
#' @param copula_dist Copula distribution ("mvn" or "mvt")
#' @param start_pars Optional starting values
#' @return List with alpha, gamma, beta, shape (if MVT), nll, convergence
#' @keywords internal
estimate_adcc_copula <- function(
    z_matrix,
    weights,
    Qbar,
    copula_dist = "mvn",
    start_pars = NULL
) {
  T_obs <- nrow(z_matrix)
  k <- ncol(z_matrix)
  
  ## Compute Nbar
  neg_resid <- z_matrix * (z_matrix < 0)
  Nbar <- crossprod(neg_resid) / T_obs
  
  ## Set up starting values
  if (is.null(start_pars)) {
    start_pars <- c(alpha = 0.05, gamma = 0.02, beta = 0.90)
    if (copula_dist == "mvt") {
      start_pars <- c(start_pars, shape = 8)
    }
  }
  
  ## Objective function
  obj_fn <- function(params) {
    adcc_copula_nll(
      params = params, 
      z_matrix = z_matrix, 
      weights = weights, 
      Qbar = Qbar, 
      Nbar = Nbar, 
      copula_dist = copula_dist
    )
  }
  
  ## Gradient function
  grad_fn <- function(params) {
    adcc_copula_gradient(
      params = params, 
      z_matrix = z_matrix, 
      weights = weights, 
      Qbar = Qbar, 
      Nbar = Nbar, 
      copula_dist = copula_dist
    )
  }
  
  ## Set bounds
  n_dcc <- 3  # alpha, gamma, beta
  lower <- c(rep(1e-8, n_dcc))
  upper <- c(rep(0.999, n_dcc))
  
  if (copula_dist == "mvt") {
    lower <- c(lower, 2.01)
    upper <- c(upper, 50)
  }
  
  ## Optimize
  result <- tryCatch({
    optim(
      par = start_pars,
      fn = obj_fn,
      gr = grad_fn,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(maxit = 500)
    )
  }, error = function(e) {
    list(par = start_pars, value = Inf, convergence = 1, message = e$message)
  })
  
  ## Extract results
  out <- list(
    alpha = result$par[1],
    gamma = result$par[2],
    beta = result$par[3],
    nll = result$value,
    convergence = result$convergence
  )
  
  if (copula_dist == "mvt") {
    out$shape <- result$par[4]
  }
  
  out$Nbar <- Nbar
  out$stationarity <- adcc_stationarity(out$alpha, out$gamma, out$beta)
  
  return(out)
}


#' @title Copula Gradient for DCC(1,1) - Analytical
#' @description Computes the analytical gradient of copula NLL w.r.t. (psi, phi)
#'   or (alpha, beta) depending on use_reparam flag.
#' @param params Parameter vector: (psi, phi) if use_reparam=TRUE, (alpha, beta) otherwise.
#'   For MVT copula, params\[3\] is the shape parameter.
#' @param z_matrix T x k matrix of copula residuals
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional covariance matrix
#' @param copula_dist "mvn" or "mvt"
#' @param use_reparam Logical; if TRUE, use reparameterized (psi, phi) space
#' @return Gradient vector (same length as params)
#' @keywords internal
copula_gradient <- function(
    params,
    z_matrix,
    weights,
    Qbar,
    copula_dist = "mvn",
    use_reparam = TRUE
) {
  ## Extract parameters
  if (use_reparam) {
    psi <- params[1]
    phi <- params[2]
    ab <- dcc11_from_unconstrained(psi, phi)
    alpha <- ab["alpha"]
    beta <- ab["beta"]
  } else {
    alpha <- params[1]
    beta <- params[2]
  }
  
  ## Shape parameter for MVT
  if (copula_dist == "mvt" && length(params) >= 3) {
    shape <- params[3]
    if (shape <= 2) {
      ## Return large gradient pointing toward valid region
      grad <- rep(0, length(params))
      grad[3] <- -1e10  # Push shape up
      return(grad)
    }
  } else {
    shape <- NULL
  }
  
  ## Check stationarity - if violated, return gradient pointing to valid region
  persistence <- alpha + beta
  if (persistence >= 1 || alpha < 0 || beta < 0) {
    if (use_reparam) {
      ## Reparameterization should prevent this, but just in case
      return(rep(0, length(params)))
    } else {
      grad <- c(0, 0)
      if (persistence >= 1) {
        grad <- c(-1e10, -1e10)  # Push both down
      }
      if (alpha < 0) grad[1] <- 1e10  # Push alpha up
      if (beta < 0) grad[2] <- 1e10   # Push beta up
      if (copula_dist == "mvt") grad <- c(grad, 0)
      return(grad)
    }
  }
  
  ## Compute gradient in original (alpha, beta) space
  grad_orig <- copula_gradient_original(
    alpha = alpha,
    beta = beta,
    z = z_matrix,
    weights = weights,
    Qbar = Qbar,
    copula_dist = copula_dist,
    shape = shape
  )
  
  ## Check for NA gradients (recursion failure)
  if (any(is.na(grad_orig))) {
    return(rep(0, length(params)))
  }
  
  ## Transform gradient if using reparameterization
  if (use_reparam) {
    ## Compute Jacobian d(alpha, beta)/d(psi, phi)
    J <- dcc11_reparam_jacobian(psi, phi)
    
    ## Chain rule: grad_reparam = J^T * grad_orig
    grad_reparam <- as.vector(t(J) %*% grad_orig[1:2])
    names(grad_reparam) <- c("psi", "phi")
    
    ## Add shape gradient if MVT
    if (copula_dist == "mvt" && length(params) >= 3) {
      grad_shape <- copula_gradient_shape(
        shape = shape,
        z = z_matrix,
        weights = weights,
        alpha = alpha,
        beta = beta,
        Qbar = Qbar
      )
      grad_reparam <- c(grad_reparam, shape = grad_shape)
    }
    
    return(grad_reparam)
  } else {
    ## Return gradient in original space
    if (copula_dist == "mvt" && length(params) >= 3) {
      grad_shape <- copula_gradient_shape(
        shape = shape,
        z = z_matrix,
        weights = weights,
        alpha = alpha,
        beta = beta,
        Qbar = Qbar
      )
      grad_orig <- c(grad_orig, shape = grad_shape)
    }
    return(grad_orig)
  }
}


#' @title Copula Gradient in Original (alpha, beta) Space
#' @description Computes analytical gradient of copula NLL w.r.t. (alpha, beta)
#' @param alpha DCC alpha parameter
#' @param beta DCC beta parameter
#' @param z T x k matrix of copula residuals
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional covariance
#' @param copula_dist "mvn" or "mvt"
#' @param shape Degrees of freedom (MVT only)
#' @return Named vector c(alpha, beta)
#' @keywords internal
copula_gradient_original <- function(alpha, beta, z, weights, Qbar,
                                     copula_dist = "mvn", shape = NULL) {
  ## Ensure plain numeric scalars
  alpha <- as.numeric(alpha)[1]
  beta <- as.numeric(beta)[1]
  
  T_obs <- nrow(z)
  k <- ncol(z)
  
  ## Forward pass with gradient storage
  fwd <- dcc11_recursion_with_grad(z, alpha, beta, Qbar)
  
  ## Check for recursion failure
  if (!isTRUE(fwd$success)) {
    return(c(alpha = NA_real_, beta = NA_real_))
  }
  
  ## Accumulate gradients
  grad_alpha <- 0
  grad_beta <- 0
  
  for (t in 1:T_obs) {
    w_t <- weights[t]
    z_t <- z[t, ]
    R_t <- fwd$R[,,t]
    R_inv_t <- fwd$R_inv[,,t]
    Q_t <- fwd$Q[,,t]
    
    ## Skip if any matrix is invalid
    if (any(!is.finite(R_t)) || any(!is.finite(R_inv_t))) {
      next
    }
    
    ## Gradient of NLL w.r.t. R_t
    if (copula_dist == "mvn") {
      grad_R <- grad_nll_wrt_R_copula_mvn(z_t, R_t, R_inv_t)
    } else {
      grad_R <- grad_nll_wrt_R_copula_mvt(z_t, R_t, R_inv_t, shape)
    }
    
    ## Backpropagate to Q_t
    grad_Q <- grad_R_to_Q_copula(grad_R, Q_t, R_t)
    
    ## Chain rule using Frobenius inner product
    grad_alpha <- grad_alpha + w_t * sum(grad_Q * fwd$dQ_dalpha[,,t])
    grad_beta <- grad_beta + w_t * sum(grad_Q * fwd$dQ_dbeta[,,t])
  }
  
  c(alpha = grad_alpha, beta = grad_beta)
}


#' @title Gradient of Copula NLL w.r.t. R_t (MVN)
#' @description Compute gradient of Gaussian copula NLL contribution w.r.t. R_t.
#'   For Gaussian copula: LL_t = -0.5 * (log|R| + z'R^{-1}z - z'z)
#'   Gradient: d(-LL_t)/dR = 0.5 * (R^{-1} - R^{-1}zz'R^{-1})
#' @param z_t k-vector of copula residuals at time t
#' @param R_t k x k correlation matrix
#' @param R_inv_t k x k inverse correlation matrix
#' @return k x k gradient matrix
#' @keywords internal
grad_nll_wrt_R_copula_mvn <- function(z_t, R_t, R_inv_t) {
  u <- as.vector(R_inv_t %*% z_t)
  0.5 * (R_inv_t - outer(u, u))
}


#' @title Gradient of Copula NLL w.r.t. R_t (MVT)
#' @description Compute gradient of Student-t copula NLL contribution w.r.t. R_t.
#'   For Student-t copula, the gradient has an additional scaling factor.
#' @param z_t k-vector of copula residuals at time t
#' @param R_t k x k correlation matrix
#' @param R_inv_t k x k inverse correlation matrix
#' @param shape Degrees of freedom
#' @return k x k gradient matrix
#' @keywords internal
grad_nll_wrt_R_copula_mvt <- function(z_t, R_t, R_inv_t, shape) {
  k <- length(z_t)
  u <- as.vector(R_inv_t %*% z_t)
  q_t <- sum(z_t * u)  # Mahalanobis distance
  kappa_t <- 1 + q_t / (shape - 2)
  weight <- (shape + k) / (2 * (shape - 2) * kappa_t)
  0.5 * R_inv_t - weight * outer(u, u)
}


#' @title Gradient of NLL w.r.t. Q_t (Copula)
#' @description Backpropagate gradient from R_t to Q_t through the normalization.
#'   R_t = D^{-1} Q_t D^{-1} where D = diag(sqrt(diag(Q_t)))
#' @param grad_R k x k gradient w.r.t. R_t
#' @param Q_t k x k Q matrix
#' @param R_t k x k correlation matrix
#' @return k x k gradient w.r.t. Q_t
#' @keywords internal
grad_R_to_Q_copula <- function(grad_R, Q_t, R_t) {
  k <- nrow(Q_t)
  
  d <- sqrt(diag(Q_t))
  d[d <= 0] <- 1e-6
  D_inv <- diag(1 / d, k)
  
  ## Main term: D^{-1} grad_R D^{-1}
  grad_Q <- D_inv %*% grad_R %*% D_inv
  
  ## Diagonal correction for the chain rule through D^{-1}
  for (i in 1:k) {
    correction <- sum(grad_R[i,] * R_t[i,]) / Q_t[i, i]
    grad_Q[i, i] <- grad_Q[i, i] - correction
  }
  
  ## Ensure symmetry
  (grad_Q + t(grad_Q)) / 2
}


#' @title Gradient of Copula NLL w.r.t. Shape (MVT)
#' @description Computes gradient of MVT copula NLL w.r.t. degrees of freedom.
#'   Uses numerical differentiation for simplicity.
#' @param shape Current degrees of freedom
#' @param z T x k matrix of copula residuals
#' @param weights Observation weights
#' @param alpha DCC alpha
#' @param beta DCC beta
#' @param Qbar Unconditional covariance
#' @return Scalar gradient
#' @keywords internal
copula_gradient_shape <- function(shape, z, weights, alpha, beta, Qbar) {
  ## Use numerical differentiation for shape parameter
  ## (analytical gradient is complex due to digamma functions)
  eps <- 1e-6
  
  nll_base <- copula_nll_fixed_dcc(
    shape = shape,
    z_matrix = z,
    weights = weights,
    alpha = alpha,
    beta = beta,
    Qbar = Qbar
  )
  
  nll_plus <- copula_nll_fixed_dcc(
    shape = shape + eps,
    z_matrix = z,
    weights = weights,
    alpha = alpha,
    beta = beta,
    Qbar = Qbar
  )
  
  (nll_plus - nll_base) / eps
}


#' @title Copula NLL with Fixed DCC Parameters
#' @description Computes copula NLL for given shape with fixed (alpha, beta).
#'   Helper for shape gradient computation.
#' @param shape Degrees of freedom
#' @param z_matrix Copula residuals
#' @param weights Observation weights
#' @param alpha DCC alpha
#' @param beta DCC beta
#' @param Qbar Unconditional covariance
#' @return Scalar NLL
#' @keywords internal
copula_nll_fixed_dcc <- function(shape, z_matrix, weights, alpha, beta, Qbar) {
  if (shape <= 2) return(1e10)
  
  ## DCC recursion
  recursion <- dcc_recursion(
    std_resid = z_matrix,
    Qbar = Qbar,
    alphas = alpha,
    betas = beta,
    verbose = FALSE
  )
  
  if (!recursion$success) return(1e10)
  
  R <- recursion$R
  T_obs <- nrow(z_matrix)
  k <- ncol(z_matrix)
  
  ## Student-t copula log-likelihood
  const_term <- lgamma(0.5 * (k + shape)) - lgamma(0.5 * shape) - 
    0.5 * k * log(pi * (shape - 2))
  
  ll_vec <- numeric(T_obs)
  
  for (t in 1:T_obs) {
    R_t <- R[,,t]
    if (any(!is.finite(R_t))) {
      ll_vec[t] <- -1e10
      next
    }
    
    det_R <- det(R_t)
    if (det_R <= 0) {
      ll_vec[t] <- -1e10
      next
    }
    
    R_inv <- tryCatch(solve(R_t), error = function(e) NULL)
    if (is.null(R_inv)) {
      ll_vec[t] <- -1e10
      next
    }
    
    z_t <- z_matrix[t, ]
    mahal <- as.numeric(t(z_t) %*% R_inv %*% z_t)
    
    ll_vec[t] <- const_term - 0.5 * log(det_R) - 
      0.5 * (shape + k) * log(1 + mahal / (shape - 2))
    
    ## Subtract marginal t log-densities
    for (j in 1:k) {
      ll_vec[t] <- ll_vec[t] - dt_scaled_log(z_t[j], shape)
    }
  }
  
  ll_vec[!is.finite(ll_vec)] <- -1e10
  -sum(weights * ll_vec)
}


#' @title Log-density of Scaled Student-t
#' @description Computes log-density of standardized Student-t (mean 0, var 1)
#' @param x Value
#' @param shape Degrees of freedom
#' @return Log-density
#' @keywords internal
dt_scaled_log <- function(x, shape) {
  scale <- sqrt(shape / (shape - 2))
  dt(x * scale, df = shape, log = TRUE) + log(scale)
}


#### ______________________________________________________________________ ####
#### PART 4: Constant Correlation Copula Likelihood                         ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


#' @title Compute Constant Copula Correlation Likelihood
#' @description Computes weighted log-likelihood for constant correlation copula
#' @param residuals Matrix of residuals
#' @param weights Observation weights
#' @param garch_pars List of GARCH parameters
#' @param dist_pars Distribution parameters
#' @param spec Model specification
#' @param transformation PIT transformation type
#' @param copula_dist Copula distribution
#' @return List with weighted_ll, dist_pars, warnings
#' @keywords internal
compute_constant_copula_likelihood <- function(
    residuals,
    weights,
    garch_pars,
    dist_pars,
    spec,
    transformation = "parametric",
    copula_dist = "mvn"
) {
  k <- ncol(residuals)
  T_obs <- nrow(residuals)
  
  ## Get standardized residuals
  std_residuals <- matrix(0, nrow = T_obs, ncol = k)
  uni_fit_list <- list()
  
  for (i in 1:k) {
    series_residuals <- residuals[, i]
    series_spec <- spec$garch_spec_args$garch_model$univariate[[i]]
    
    uni_spec_obj <- tsgarch::garch_modelspec(
      y = xts::xts(series_residuals, order.by = Sys.Date() - (T_obs:1)),
      model = series_spec$model,
      garch_order = series_spec$garch_order,
      distribution = series_spec$distribution
    )
    
    for (par_name in names(garch_pars[[i]])) {
      if (par_name %in% uni_spec_obj$parmatrix$parameter) {
        row_idx <- which(uni_spec_obj$parmatrix$parameter == par_name)
        uni_spec_obj$parmatrix[row_idx, "value"] <- garch_pars[[i]][[par_name]]
      }
    }
    
    uni_fit <- tsmethods::tsfilter(uni_spec_obj)
    sigma_vec <- as.numeric(uni_fit$sigma)
    std_residuals[, i] <- series_residuals / sigma_vec
    uni_fit_list[[i]] <- uni_fit
  }
  
  ## PIT transform
  u_matrix <- compute_pit_transform(
    std_residuals = std_residuals,
    uni_fit_list = uni_fit_list,
    transformation = transformation,
    copula_dist = copula_dist,
    dist_pars = dist_pars
  )
  
  u_matrix[u_matrix < 3.330669e-16] <- 2.220446e-16
  u_matrix[u_matrix > 0.99999] <- 0.99999
  
  ## Transform to copula residuals
  z_matrix <- compute_copula_residuals(
    u_matrix = u_matrix,
    copula_dist = copula_dist,
    dist_pars = dist_pars
  )
  
  ## Compute constant correlation
  R_const <- stats::cov.wt(z_matrix, wt = weights, cor = TRUE)$cor
  
  ## Ensure positive definite
  eig <- eigen(R_const, symmetric = TRUE)
  if (any(eig$values < 1e-8)) {
    R_const <- R_const + diag(1e-6, k)
    R_const <- cov2cor(R_const)
  }
  
  ## Compute log-likelihood
  ll_vec <- numeric(T_obs)
  
  det_R <- det(R_const)
  R_inv <- solve(R_const)
  
  if (copula_dist == "mvn") {
    for (t in 1:T_obs) {
      z_t <- z_matrix[t, ]
      mahal <- as.numeric(t(z_t) %*% R_inv %*% z_t)
      ll_vec[t] <- -0.5 * (log(det_R) + mahal - as.numeric(t(z_t) %*% z_t))
    }
  } else {
    shape <- dist_pars$shape %||% 4
    const_term <- lgamma(0.5 * (k + shape)) - lgamma(0.5 * shape) - 
      0.5 * k * log(pi * (shape - 2))
    
    for (t in 1:T_obs) {
      z_t <- z_matrix[t, ]
      mahal <- as.numeric(t(z_t) %*% R_inv %*% z_t)
      
      ll_vec[t] <- const_term - 0.5 * log(det_R) - 
        0.5 * (shape + k) * log(1 + mahal / (shape - 2))
      
      for (j in 1:k) {
        ll_vec[t] <- ll_vec[t] - dt_scaled_log(z_t[j], shape)
      }
    }
  }
  
  ll_vec[!is.finite(ll_vec)] <- -1e10
  weighted_ll <- sum(weights * ll_vec)
  
  return(list(
    weighted_ll = weighted_ll,
    dist_pars = dist_pars,
    warnings = list()
  ))
}


## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## END OF COPULA GARCH IMPLEMENTATION
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =