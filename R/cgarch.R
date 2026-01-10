## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## Copula GARCH Implementation for MS-VARMA-GARCH Bootstrap (tsbs)
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
##
## This file implements the Copula GARCH model (cgarch_modelspec) as an 
## alternative to the DCC model for the tsbs package's Markov-Switching 
## VARMA-GARCH framework.
##
## Key differences from DCC:
## - Copula GARCH uses probability integral transform (PIT) on standardized
##   residuals to create uniform margins, then applies copula for dependence
## - Supports "parametric", "empirical", and "spd" transformations
## - Available copula families: "mvn" (Gaussian) and "mvt" (Student-t)
## - Log-likelihood = sum(univariate GARCH LL) + copula LL
##   (whereas DCC combines them differently)
##
## Structure:
## - - - - - -
## PART 1: Copula GARCH Weighted Estimation (Stage 1 + Stage 2)
## PART 2: Copula Parameters Weighted Estimation 
## PART 3: Helper Functions for Copula Computations
## PART 4: Constant Correlation Copula Likelihood
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
#'    uniform [0,1] margins via probability integral transform (PIT). The 
#'    transformation can be:
#'    - "parametric": Uses the estimated univariate distribution's CDF
#'    - "empirical": Uses the empirical CDF
#'    - "spd": Uses semi-parametric distribution
#'
#' 2. **Copula Specification**: The uniform margins are then transformed
#'    according to the copula distribution:
#'    - "mvn": Multivariate Normal copula (Gaussian copula)
#'    - "mvt": Multivariate Student-t copula
#'
#' 3. **Correlation Dynamics**: Same as DCC - can be "constant", "dcc", or "adcc"
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
    diagnostics = NULL,
    iteration = NULL,
    state = NULL,
    verbose = FALSE
) {
  
  k <- ncol(residuals)
  T_obs <- nrow(residuals)
  
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
  use_analytical_gradient <- is_dcc11(dcc_start_pars)
  
  if (verbose) {
    if (use_analytical_gradient) {
      cat("Copula GARCH DCC(1,1) detected: Using analytical gradient with reparameterization\n")
    } else {
      order <- get_dcc_order(dcc_start_pars)
      cat(sprintf("Copula GARCH DCC(%d,%d) detected: Using finite differences\n",
                  order["p"], order["q"]))
    }
  }
  
  warnings_list <- list()
  
  ## 7. Optimization ===========================================================
  
  if (use_analytical_gradient) {
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
#' @description Transforms standardized residuals to uniform [0,1] margins
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
    ## Semi-parametric distribution - would need tsdistributions::spd
    ## For now, fall back to empirical
    warning("SPD transformation not fully implemented, using empirical.\n")
    for (i in 1:k) {
      ecdf_fn <- ecdf(std_residuals[, i])
      u_matrix[, i] <- ecdf_fn(std_residuals[, i])
    }
    u_matrix <- u_matrix * (T_obs / (T_obs + 1))
  }
  
  return(u_matrix)
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


#' @title Copula Gradient for DCC(1,1)
#' @description Computes the gradient of copula NLL using numerical differentiation
#' @param params Parameter vector
#' @param z_matrix Copula residuals
#' @param weights Observation weights
#' @param Qbar Unconditional covariance
#' @param copula_dist Copula distribution
#' @param use_reparam Logical; reparameterization flag
#' @return Gradient vector
#' @keywords internal
copula_gradient <- function(
    params,
    z_matrix,
    weights,
    Qbar,
    copula_dist = "mvn",
    use_reparam = TRUE
) {
  ## Use numerical gradient (can be replaced with analytical later)
  eps <- 1e-6
  n_params <- length(params)
  grad <- numeric(n_params)
  
  f0 <- copula_nll(params, z_matrix, weights, Qbar, copula_dist, use_reparam)
  
  for (i in 1:n_params) {
    params_plus <- params
    params_plus[i] <- params_plus[i] + eps
    f_plus <- copula_nll(params_plus, z_matrix, weights, Qbar, copula_dist, use_reparam)
    grad[i] <- (f_plus - f0) / eps
  }
  
  return(grad)
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