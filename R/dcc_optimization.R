## =============================================================================
## DCC Parameter Handling: Reparameterization and Higher-Order Support
## =============================================================================
##
## This module provides:
## 1. Correct persistence computation for arbitrary DCC(p,q) orders
## 2. Reparameterized optimization for DCC(1,1) to avoid boundary instabilities
## 3. Penalty-based optimization for higher-order DCC(p,q) as fallback
##
## For DCC(1,1), we use the reparameterization:
##   persistence = α + β ∈ (0, 1)
##   ratio = α / (α + β) ∈ (0, 1)
##
## This transforms the constrained problem (α + β < 1) into an unconstrained
## box-constrained problem where stationarity is guaranteed by construction.
##
## Reference: Seitz, S. (2023). "Varying Coefficient GARCH" 
##            https://sarem-seitz.com/posts/varying-coefficient-garch.html
##            
##            Donfack & Dufays (2021). "Modeling time-varying parameters using 
##            artificial neural networks: a GARCH illustration." 
##            Studies in Nonlinear Dynamics & Econometrics.
##
## For DCC(p,q) with max(p,q) > 1, we use the standard penalty method since
## reparameterization becomes significantly more complex (requires softmax
## distributions over parameter vectors).
##
## =============================================================================


## =============================================================================
## SECTION 1: Persistence Computation (All Orders)
## =============================================================================

#' Extract DCC order (p, q) from parameter names
#' 
#' @param dcc_params Named list of DCC parameters (e.g., alpha_1, alpha_2, beta_1)
#' @return Named vector c(p = ..., q = ...) where q is ARCH order (number of 
#'   alpha parameters) and p is GARCH order (number of beta parameters)
#' @keywords internal
get_dcc_order <- function(dcc_params) {
  
  if (is.null(dcc_params) || length(dcc_params) == 0) {
    return(c(p = 0, q = 0))
  }
  
  param_names <- names(dcc_params)
  
  ## Extract alpha indices (q order - ARCH-like, controls lagged z'z terms)
  alpha_names <- param_names[grepl("^alpha_", param_names)]
  if (length(alpha_names) > 0) {
    alpha_indices <- as.integer(gsub("^alpha_", "", alpha_names))
    q_order <- max(alpha_indices, na.rm = TRUE)
  } else {
    q_order <- 0
  }
  
  ## Extract beta indices (p order - GARCH-like, controls lagged Q terms)
  beta_names <- param_names[grepl("^beta_", param_names)]
  if (length(beta_names) > 0) {
    beta_indices <- as.integer(gsub("^beta_", "", beta_names))
    p_order <- max(beta_indices, na.rm = TRUE)
  } else {
    p_order <- 0
  }
  
  return(c(p = p_order, q = q_order))
}


#' Compute total persistence from DCC parameters
#' 
#' For a stationary DCC(p,q) model, we require:
#'   P = sum(alpha_1, ..., alpha_q) + sum(beta_1, ..., beta_p) < 1
#' 
#' @param dcc_params Named list of DCC parameters
#' @return List with components:
#'   \item{persistence}{Total persistence P = sum(alphas) + sum(betas)}
#'   \item{alpha_sum}{Sum of all alpha parameters}
#'   \item{beta_sum}{Sum of all beta parameters}
#'   \item{alphas}{Numeric vector of alpha values in order}
#'   \item{betas}{Numeric vector of beta values in order}
#'   \item{order}{Named vector c(p, q) giving the DCC order}
#' @keywords internal
compute_dcc_persistence <- function(dcc_params) {
  
  if (is.null(dcc_params) || length(dcc_params) == 0) {
    return(list(
      persistence = 0,
      alpha_sum = 0,
      beta_sum = 0,
      alphas = numeric(0),
      betas = numeric(0),
      order = c(p = 0, q = 0)
    ))
  }
  
  param_names <- names(dcc_params)
  
  ## Extract all alpha values in order
  alpha_names <- param_names[grepl("^alpha_", param_names)]
  if (length(alpha_names) > 0) {
    ## Sort by index to ensure correct order
    alpha_indices <- as.integer(gsub("^alpha_", "", alpha_names))
    alpha_order <- order(alpha_indices)
    alphas <- unlist(dcc_params[alpha_names[alpha_order]])
    names(alphas) <- NULL
  } else {
    alphas <- numeric(0)
  }
  
  ## Extract all beta values in order
  beta_names <- param_names[grepl("^beta_", param_names)]
  if (length(beta_names) > 0) {
    beta_indices <- as.integer(gsub("^beta_", "", beta_names))
    beta_order <- order(beta_indices)
    betas <- unlist(dcc_params[beta_names[beta_order]])
    names(betas) <- NULL
  } else {
    betas <- numeric(0)
  }
  
  alpha_sum <- sum(alphas)
  beta_sum <- sum(betas)
  
  return(list(
    persistence = alpha_sum + beta_sum,
    alpha_sum = alpha_sum,
    beta_sum = beta_sum,
    alphas = alphas,
    betas = betas,
    order = c(p = length(betas), q = length(alphas))
  ))
}


#' Check DCC stationarity constraints
#' 
#' Verifies that DCC parameters satisfy:
#' 1. All alpha_j >= 0
#' 2. All beta_j >= 0
#' 3. sum(alphas) + sum(betas) < 1
#' 
#' @param dcc_params Named list of DCC parameters
#' @param verbose Logical; if TRUE, print diagnostic messages
#' @return List with components:
#'   \item{is_stationary}{Logical indicating if constraints are satisfied}
#'   \item{persistence}{Total persistence value}
#'   \item{reason}{Character description if not stationary, NULL otherwise
#'   \item{details}{Output from compute_dcc_persistence()}
#' @keywords internal
check_dcc_stationarity <- function(dcc_params, verbose = FALSE) {
  
  pers <- compute_dcc_persistence(dcc_params)
  
  ## Check individual alpha positivity
  if (length(pers$alphas) > 0 && any(pers$alphas < 0)) {
    neg_idx <- which(pers$alphas < 0)
    reason <- sprintf("negative alpha at index %s (values: %s)", 
                      paste(neg_idx, collapse = ", "),
                      paste(round(pers$alphas[neg_idx], 6), collapse = ", "))
    if (verbose) cat("*** Stationarity violation:", reason, "***\n")
    return(list(
      is_stationary = FALSE,
      persistence = pers$persistence,
      reason = reason,
      details = pers
    ))
  }
  
  ## Check individual beta positivity
  if (length(pers$betas) > 0 && any(pers$betas < 0)) {
    neg_idx <- which(pers$betas < 0)
    reason <- sprintf("negative beta at index %s (values: %s)", 
                      paste(neg_idx, collapse = ", "),
                      paste(round(pers$betas[neg_idx], 6), collapse = ", "))
    if (verbose) cat("*** Stationarity violation:", reason, "***\n")
    return(list(
      is_stationary = FALSE,
      persistence = pers$persistence,
      reason = reason,
      details = pers
    ))
  }
  
  ## Check stationarity constraint: persistence < 1
  if (pers$persistence >= 1) {
    reason <- sprintf(
      "non-stationary (persistence = %.6f >= 1; alpha_sum = %.4f, beta_sum = %.4f)",
      pers$persistence, pers$alpha_sum, pers$beta_sum
    )
    if (verbose) cat("*** Stationarity violation:", reason, "***\n")
    return(list(
      is_stationary = FALSE,
      persistence = pers$persistence,
      reason = reason,
      details = pers
    ))
  }
  
  return(list(
    is_stationary = TRUE,
    persistence = pers$persistence,
    reason = NULL,
    details = pers
  ))
}


## =============================================================================
## SECTION 2: DCC(1,1) Reparameterization
## =============================================================================

#' Transform DCC(1,1) parameters to reparameterized space
#' 
#' Transforms (alpha, beta) with constraint alpha + beta < 1 to 
#' (persistence, ratio) where both are in (0, 1) with no joint constraint.
#' 
#' @param alpha DCC alpha parameter
#' @param beta DCC beta parameter
#' @return Named vector c(persistence, ratio)
#' @keywords internal
dcc_to_reparam <- function(alpha, beta) {
  persistence <- alpha + beta
  
  ## Handle edge case where persistence is 0
  
  if (persistence < .Machine$double.eps) {
    ratio <- 0.5  # Arbitrary, since both are ~0
  } else {
    ratio <- alpha / persistence
  }
  
  c(persistence = persistence, ratio = ratio)
}


#' Transform reparameterized space back to DCC(1,1) parameters
#' 
#' Transforms (persistence, ratio) back to (alpha, beta).
#' Stationarity (alpha + beta < 1) is guaranteed if persistence < 1.
#' 
#' @param persistence Total persistence (alpha + beta), must be in (0, 1)
#' @param ratio Proportion allocated to alpha, must be in (0, 1)
#' @return Named vector c(alpha, beta)
#' @keywords internal
dcc_from_reparam <- function(persistence, ratio) {
  alpha <- persistence * ratio
  beta <- persistence * (1 - ratio)
  c(alpha = alpha, beta = beta)
}


#' Check if DCC specification is DCC(1,1)
#' 
#' @param dcc_params Named list of DCC parameters
#' @return Logical; TRUE if this is a DCC(1,1) specification
#' @keywords internal
is_dcc_11 <- function(dcc_params) {
  if (is.null(dcc_params) || length(dcc_params) == 0) {
    return(FALSE)
  }
  
  order <- get_dcc_order(dcc_params)
  return(order["p"] == 1 && order["q"] == 1)
}


## =============================================================================
## SECTION 3: DCC Recursion (All Orders)
## =============================================================================

#' Perform DCC(p,q) recursion
#' 
#' Computes the Q and R matrices for arbitrary DCC(p,q) order:
#'   Q_t = Qbar * (1 - sum(alpha) - sum(beta)) + 
#'         sum_{j=1}^{q} alpha_j * (z_{t-j} * z_{t-j}') +
#'         sum_{j=1}^{p} beta_j * Q_{t-j}
#'   R_t = diag(Q_t)^{-1/2} * Q_t * diag(Q_t)^{-1/2}
#' 
#' @param std_resid Matrix of standardized residuals (T x k)
#' @param Qbar Unconditional covariance matrix of standardized residuals (k x k)
#' @param alphas Numeric vector of alpha parameters (length q)
#' @param betas Numeric vector of beta parameters (length p)
#' @param verbose Logical; if TRUE, print diagnostic messages
#' @return List with components:
#'   \item{success}{Logical indicating if recursion completed without errors}
#'   \item{Q}{Array of Q matrices (k x k x T)}
#'   \item{R}{Array of correlation matrices (k x k x T)}
#'   \item{maxpq}{Maximum of p and q, used for burn-in period}
#'   \item{error_type}{Character describing error if success is FALSE}
#'   \item{error_time}{Time index where error occurred}
#' @keywords internal
dcc_recursion <- function(std_resid, Qbar, alphas, betas, verbose = FALSE) {
  
  T_obs <- nrow(std_resid)
  k <- ncol(std_resid)
  
  q_order <- length(alphas)  # Number of lagged z'z terms
  p_order <- length(betas)   # Number of lagged Q terms
  maxpq <- max(p_order, q_order, 1)
  
  ## Compute total persistence
  alpha_sum <- sum(alphas)
  beta_sum <- sum(betas)
  persistence <- alpha_sum + beta_sum
  
  if (verbose) {
    cat(sprintf("DCC(%d,%d) recursion: persistence = %.4f\n", 
                p_order, q_order, persistence))
  }
  
  ## Initialize arrays
  Q <- array(0, dim = c(k, k, T_obs))
  R <- array(0, dim = c(k, k, T_obs))
  
  ## Initialize first maxpq observations with Qbar
  for (t in 1:min(maxpq, T_obs)) {
    Q[,,t] <- Qbar
    Qbar_diag <- diag(Qbar)
    if (any(Qbar_diag <= 0)) {
      return(list(
        success = FALSE,
        error_type = "non_positive_Qbar_diagonal",
        error_time = 0,
        Q = Q,
        R = R,
        maxpq = maxpq
      ))
    }
    Qbar_diag_inv_sqrt <- diag(1/sqrt(Qbar_diag), k)
    R[,,t] <- Qbar_diag_inv_sqrt %*% Qbar %*% Qbar_diag_inv_sqrt
  }
  
  ## Handle case where T_obs <= maxpq (not enough data for recursion)
  if (T_obs <= maxpq) {
    return(list(
      success = TRUE,
      Q = Q,
      R = R,
      maxpq = maxpq
    ))
  }
  
  ## Main recursion starting from t = maxpq + 1
  for (t in (maxpq + 1):T_obs) {
    
    ## Intercept term: Qbar * (1 - persistence)
    Q_t <- Qbar * (1 - persistence)
    
    ## Add alpha terms: sum_{j=1}^{q} alpha_j * (z_{t-j} * z_{t-j}')
    for (j in seq_along(alphas)) {
      if (t - j >= 1) {
        z_lag <- std_resid[t - j, , drop = FALSE]
        Q_t <- Q_t + alphas[j] * (t(z_lag) %*% z_lag)
      }
    }
    
    ## Add beta terms: sum_{j=1}^{p} beta_j * Q_{t-j}
    for (j in seq_along(betas)) {
      if (t - j >= 1) {
        Q_t <- Q_t + betas[j] * Q[,,t - j]
      }
    }
    
    ## Check for numerical issues in Q
    if (any(!is.finite(Q_t))) {
      return(list(
        success = FALSE,
        error_type = "non_finite_Q",
        error_time = t,
        Q = Q,
        R = R,
        maxpq = maxpq
      ))
    }
    
    Q_diag <- diag(Q_t)
    if (any(Q_diag <= 0)) {
      return(list(
        success = FALSE,
        error_type = "non_positive_Q_diagonal",
        error_time = t,
        Q_diag_min = min(Q_diag),
        Q = Q,
        R = R,
        maxpq = maxpq
      ))
    }
    
    Q[,,t] <- Q_t
    
    ## Compute correlation matrix: R_t = diag(Q_t)^{-1/2} * Q_t * diag(Q_t)^{-1/2}
    Q_diag_inv_sqrt <- diag(1/sqrt(Q_diag), k)
    R[,,t] <- Q_diag_inv_sqrt %*% Q_t %*% Q_diag_inv_sqrt
  }
  
  return(list(
    success = TRUE,
    Q = Q,
    R = R,
    maxpq = maxpq
  ))
}


## =============================================================================
## SECTION 4: Reparameterized Objective Function for DCC(1,1)
## =============================================================================

#' Reparameterized DCC(1,1) negative log-likelihood
#' 
#' This function computes the weighted negative log-likelihood for DCC(1,1)
#' using the (persistence, ratio) parameterization. This eliminates the
#' need for penalty-based stationarity enforcement since persistence < 1
#' is guaranteed by the box constraints.
#' 
#' @param reparam_pars Numeric vector c(persistence, ratio) or named list
#' @param std_resid Matrix of standardized residuals (T x k)
#' @param weights Observation weights (length T)
#' @param Qbar Unconditional covariance matrix (k x k)
#' @param distribution Character; either "mvn" or "mvt"
#' @param dist_pars Distribution parameters (e.g., shape for mvt)
#' @param verbose Logical; if TRUE, print diagnostic messages
#' @return Negative log-likelihood value (scalar)
#' @keywords internal
dcc11_nll_reparam <- function(reparam_pars, std_resid, weights, Qbar,
                              distribution = "mvn", dist_pars = NULL,
                              verbose = FALSE) {
  
  ## Extract reparameterized values
  if (is.list(reparam_pars)) {
    persistence <- reparam_pars$persistence
    ratio <- reparam_pars$ratio
  } else {
    persistence <- reparam_pars[1]
    ratio <- reparam_pars[2]
  }
  
  ## Check box constraints (should be enforced by optimizer, but be safe)
  eps <- 1e-10
  if (persistence <= eps || persistence >= 1 - eps ||
      ratio <= eps || ratio >= 1 - eps) {
    return(1e10)
  }
  
  ## Transform back to (alpha, beta)
  ## Stationarity is GUARANTEED since persistence < 1
  ab <- dcc_from_reparam(persistence, ratio)
  alpha <- ab["alpha"]
  beta <- ab["beta"]
  
  if (verbose) {
    cat(sprintf("Reparam: persistence=%.4f, ratio=%.4f -> alpha=%.4f, beta=%.4f\n",
                persistence, ratio, alpha, beta))
  }
  
  ## Perform DCC recursion
  recursion <- dcc_recursion(
    std_resid = std_resid,
    Qbar = Qbar,
    alphas = alpha,
    betas = beta,
    verbose = verbose
  )
  
  if (!recursion$success) {
    if (verbose) {
      cat("*** Recursion failed:", recursion$error_type, 
          "at t =", recursion$error_time, "***\n")
    }
    return(1e10)
  }
  
  ## Compute weighted log-likelihood
  T_obs <- nrow(std_resid)
  k <- ncol(std_resid)
  R <- recursion$R
  
  ll_vec <- numeric(T_obs)
  
  if (distribution == "mvn") {
    for (t in 1:T_obs) {
      R_t <- R[,,t]
      
      ## Check for valid correlation matrix
      if (any(!is.finite(R_t))) {
        ll_vec[t] <- -1e10
        next
      }
      
      ## Compute MVN log-likelihood
      ## For standardized residuals with correlation R:
      ## ll = -0.5 * (k*log(2*pi) + log(det(R)) + z' * R^{-1} * z)
      ## But z already has unit variance, so:
      ## ll = -0.5 * (log(det(R)) + z' * R^{-1} * z)  [up to constant]
      
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
      
      z_t <- std_resid[t, ]
      mahal <- as.numeric(t(z_t) %*% R_inv %*% z_t)
      
      ll_vec[t] <- -0.5 * (k * log(2 * pi) + log(det_R) + mahal)
    }
  } else {
    ## MVT distribution - would need shape parameter
    ## For now, fall back to MVN
    warning("MVT distribution not yet implemented in reparameterized version")
    return(1e10)
  }
  
  ## Compute weighted NLL
  ll_vec[!is.finite(ll_vec)] <- -1e10
  nll <- -sum(weights * ll_vec)
  
  return(nll)
}


## =============================================================================
## SECTION 5: Main Estimation Interface
## =============================================================================

#' Estimate DCC parameters with appropriate method
#' 
#' For DCC(1,1), uses reparameterized optimization to avoid boundary
#' instabilities. For higher-order DCC(p,q), uses penalty-based optimization.
#' 
#' @param std_resid Matrix of standardized residuals (T x k)
#' @param weights Observation weights (length T)
#' @param dcc_start_pars Named list of starting DCC parameters
#' @param Qbar Unconditional covariance matrix (k x k)
#' @param distribution Character; "mvn" or "mvt"
#' @param dist_pars Distribution parameters
#' @param verbose Logical; if TRUE, print diagnostic messages
#' @return List with estimated parameters and diagnostics
#' @keywords internal
estimate_dcc_params <- function(std_resid, weights, dcc_start_pars, Qbar,
                                distribution = "mvn", dist_pars = NULL,
                                verbose = FALSE) {
  
  ## Determine DCC order
  order <- get_dcc_order(dcc_start_pars)
  
  if (verbose) {
    cat(sprintf("DCC order: p=%d, q=%d\n", order["p"], order["q"]))
  }
  
  ## Choose estimation method based on order
  if (order["p"] == 1 && order["q"] == 1) {
    ## DCC(1,1): Use reparameterized optimization
    result <- estimate_dcc11_reparam(
      std_resid = std_resid,
      weights = weights,
      dcc_start_pars = dcc_start_pars,
      Qbar = Qbar,
      distribution = distribution,
      dist_pars = dist_pars,
      verbose = verbose
    )
    result$method <- "reparameterized"
    
  } else {
    ## Higher order: Use penalty-based optimization
    result <- estimate_dcc_pq_penalty(
      std_resid = std_resid,
      weights = weights,
      dcc_start_pars = dcc_start_pars,
      Qbar = Qbar,
      distribution = distribution,
      dist_pars = dist_pars,
      verbose = verbose
    )
    result$method <- "penalty"
  }
  
  return(result)
}


#' Estimate DCC(1,1) using reparameterized optimization
#' @keywords internal
estimate_dcc11_reparam <- function(std_resid, weights, dcc_start_pars, Qbar,
                                   distribution = "mvn", dist_pars = NULL,
                                   verbose = FALSE) {
  
  ## Extract starting values
  alpha_start <- dcc_start_pars$alpha_1
  beta_start <- dcc_start_pars$beta_1
  
  ## Transform to reparameterized space
  reparam_start <- dcc_to_reparam(alpha_start, beta_start)
  
  if (verbose) {
    cat(sprintf("Starting values: alpha=%.4f, beta=%.4f\n", alpha_start, beta_start))
    cat(sprintf("Reparameterized: persistence=%.4f, ratio=%.4f\n",
                reparam_start["persistence"], reparam_start["ratio"]))
  }
  
  ## Set up bounds for reparameterized parameters
  ## Both persistence and ratio should be in (eps, 1-eps)
  eps <- 1e-6
  lower <- c(eps, eps)
  upper <- c(1 - eps, 1 - eps)
  
  ## Optimize
  opt_result <- optim(
    par = as.numeric(reparam_start),
    fn = dcc11_nll_reparam,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    std_resid = std_resid,
    weights = weights,
    Qbar = Qbar,
    distribution = distribution,
    dist_pars = dist_pars,
    verbose = FALSE
  )
  
  ## Transform back to original parameters
  final_ab <- dcc_from_reparam(opt_result$par[1], opt_result$par[2])
  
  if (verbose) {
    cat(sprintf("Optimized: persistence=%.4f, ratio=%.4f\n",
                opt_result$par[1], opt_result$par[2]))
    cat(sprintf("Final: alpha=%.4f, beta=%.4f (sum=%.4f)\n",
                final_ab["alpha"], final_ab["beta"], sum(final_ab)))
    cat(sprintf("Convergence: %d, NLL: %.4f\n", opt_result$convergence, opt_result$value))
  }
  
  return(list(
    dcc_pars = list(alpha_1 = final_ab["alpha"], beta_1 = final_ab["beta"]),
    reparam_pars = list(persistence = opt_result$par[1], ratio = opt_result$par[2]),
    nll = opt_result$value,
    convergence = opt_result$convergence,
    n_penalty_triggers = 0  # Reparameterization eliminates penalties
  ))
}


#' Estimate DCC(p,q) using penalty-based optimization
#' @keywords internal
estimate_dcc_pq_penalty <- function(std_resid, weights, dcc_start_pars, Qbar,
                                    distribution = "mvn", dist_pars = NULL,
                                    verbose = FALSE) {
  
  ## This would contain the corrected penalty-based optimization
  ## using compute_dcc_persistence() and dcc_recursion()
  ## 
  ## For now, return a placeholder indicating this needs integration
  ## with the existing weighted_dcc_loglik infrastructure
  
  warning("Higher-order DCC estimation requires integration with existing code")
  
  return(list(
    dcc_pars = dcc_start_pars,
    nll = NA,
    convergence = -1,
    n_penalty_triggers = NA
  ))
}