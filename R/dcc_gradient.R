## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## DCC(1,1) Analytical Gradient Implementation
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
##
## This file implements analytical gradients for DCC(1,1) models.
## The gradients support:
##   - Multivariate Normal (MVN) distribution
##   - Multivariate Student-t (MVT) distribution
##   - Weighted log-likelihood (for EM algorithm M-step)
##   - Reparameterized space (psi/phi) for unconstrained optimization
##
## For higher-order DCC(p,q) with p > 1 or q > 1, the optimizer falls back to
## finite difference approximation with ndeps = 1e-12.
##
## Reparameterization scheme (for DCC(1,1)):
##   psi = logit(alpha + beta)  -- persistence in logit space
##   phi = log(alpha / beta)    -- ratio in log space
##
## This transforms the constrained problem (alpha > 0, beta > 0, alpha + beta < 1)
## into a fully unconstrained problem where stationarity is guaranteed by
## construction for any (psi, phi) in R^2.
##
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


## SECTION 1: DCC Order Detection and Persistence Computation ==================

# #' @title Extract DCC Order from Parameter Names
# #' @description Determines (p, q) order from parameter naming convention.
# #' @param dcc_params Named list of DCC parameters (e.g., list(alpha_1 = 0.05, beta_1 = 0.9))
# #' @return Named vector c(p = ..., q = ...) where p is GARCH order (beta count)
# #'   and q is ARCH order (alpha count)
# #' @keywords internal
# get_dcc_order <- function(dcc_params) {
#   if (is.null(dcc_params) || length(dcc_params) == 0) {
#     return(c(p = 0, q = 0))
#   }
#   
#   param_names <- names(dcc_params)
#   
#   alpha_names <- param_names[grepl("^alpha_[0-9]+$", param_names)]
#   if (length(alpha_names) > 0) {
#     alpha_indices <- as.integer(gsub("^alpha_", "", alpha_names))
#     q_order <- max(alpha_indices, na.rm = TRUE)
#   } else {
#     q_order <- 0
#   }
#   
#   beta_names <- param_names[grepl("^beta_[0-9]+$", param_names)]
#   if (length(beta_names) > 0) {
#     beta_indices <- as.integer(gsub("^beta_", "", beta_names))
#     p_order <- max(beta_indices, na.rm = TRUE)
#   } else {
#     p_order <- 0
#   }
#   
#   return(c(p = p_order, q = q_order))
# }


#' @title Check if DCC Order is (1,1)
#' @description Determines if the DCC model is of order (1,1), which allows
#'   for analytical gradient computation and reparameterized optimization.
#' @param dcc_pars Named list of DCC parameters
#' @return Logical: TRUE if DCC(1,1), FALSE otherwise
#' @keywords internal
is_dcc11 <- function(dcc_pars) {
  if (is.null(dcc_params) || length(dcc_params) == 0) {
    return(FALSE)
  }
  
  order <- get_dcc_order(dcc_pars)
  return(order["p"] == 1 && order["q"] == 1)
}


#' @title Compute DCC Persistence and Extract Parameters
#' @description Extracts all alpha/beta parameters and computes total persistence.
#'   For a stationary DCC(p,q) model, we require:
#'   P = sum(alpha_1, ..., alpha_q) + sum(beta_1, ..., beta_p) < 1
#' @param dcc_params Named list of DCC parameters
#' @return List with components:
#'   \item{persistence}{Total persistence = sum(alphas) + sum(betas)}
#'   \item{alpha_sum}{Sum of all alpha parameters}
#'   \item{beta_sum}{Sum of all beta parameters}
#'   \item{alphas}{Numeric vector of alpha values in order}
#'   \item{betas}{Numeric vector of beta values in order}
#'   \item{order}{Named vector c(p, q)}
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
  
  ## Extract alphas in order
  alpha_names <- param_names[grepl("^alpha_[0-9]+$", param_names)]
  if (length(alpha_names) > 0) {
    alpha_indices <- as.integer(gsub("^alpha_", "", alpha_names))
    alpha_order <- order(alpha_indices)
    alphas <- unlist(dcc_params[alpha_names[alpha_order]])
    names(alphas) <- NULL
  } else {
    alphas <- numeric(0)
  }
  
  ## Extract betas in order
  beta_names <- param_names[grepl("^beta_[0-9]+$", param_names)]
  if (length(beta_names) > 0) {
    beta_indices <- as.integer(gsub("^beta_", "", beta_names))
    beta_order <- order(beta_indices)
    betas <- unlist(dcc_params[beta_names[beta_order]])
    names(betas) <- NULL
  } else {
    betas <- numeric(0)
  }
  
  return(list(
    persistence = sum(alphas) + sum(betas),
    alpha_sum = sum(alphas),
    beta_sum = sum(betas),
    alphas = alphas,
    betas = betas,
    order = c(p = length(betas), q = length(alphas))
  ))
}


# #' @title Check DCC Stationarity Constraints
# #' @description Verifies that DCC parameters satisfy stationarity conditions.
# #' @param dcc_params Named list of DCC parameters
# #' @param verbose Logical; if TRUE, print diagnostic messages
# #' @return List with components:
# #'   \item{is_stationary}{Logical}
# #'   \item{persistence}{Total persistence value}
# #'   \item{reason}{Character description if not stationary, NULL otherwise}
# #'   \item{details}{Output from compute_dcc_persistence()}
# #' @keywords internal
# check_dcc_stationarity <- function(dcc_params, verbose = FALSE) {
#   pers <- compute_dcc_persistence(dcc_params)
#   
#   ## Check alpha positivity
#   if (length(pers$alphas) > 0 && any(pers$alphas < 0)) {
#     neg_idx <- which(pers$alphas < 0)
#     reason <- sprintf("negative alpha at index %s (values: %s)",
#                       paste(neg_idx, collapse = ", "),
#                       paste(round(pers$alphas[neg_idx], 6), collapse = ", "))
#     if (verbose) cat("*** Stationarity violation:", reason, "***\n")
#     return(list(is_stationary = FALSE, persistence = pers$persistence,
#                 reason = reason, details = pers))
#   }
#   
#   ## Check beta positivity
#   if (length(pers$betas) > 0 && any(pers$betas < 0)) {
#     neg_idx <- which(pers$betas < 0)
#     reason <- sprintf("negative beta at index %s (values: %s)",
#                       paste(neg_idx, collapse = ", "),
#                       paste(round(pers$betas[neg_idx], 6), collapse = ", "))
#     if (verbose) cat("*** Stationarity violation:", reason, "***\n")
#     return(list(is_stationary = FALSE, persistence = pers$persistence,
#                 reason = reason, details = pers))
#   }
#   
#   ## Check persistence < 1
#   if (pers$persistence >= 1) {
#     reason <- sprintf("non-stationary (persistence = %.6f >= 1)",
#                       pers$persistence)
#     if (verbose) cat("*** Stationarity violation:", reason, "***\n")
#     return(list(is_stationary = FALSE, persistence = pers$persistence,
#                 reason = reason, details = pers))
#   }
#   
#   return(list(is_stationary = TRUE, persistence = pers$persistence,
#               reason = NULL, details = pers))
# }


## SECTION 2: DCC(1,1) Reparameterization ======================================

#' @title Reparameterize DCC(1,1) to Unconstrained Space
#' @description Transform constrained (alpha, beta) to unconstrained (psi, phi).
#'   psi = logit(alpha + beta) controls persistence
#'   phi = log(alpha / beta) controls the ratio
#' @param alpha DCC alpha parameter (0 < alpha)
#' @param beta DCC beta parameter (0 < beta)
#' @return Named vector c(psi, phi)
#' @keywords internal
dcc11_to_unconstrained <- function(alpha, beta) {
  ## Strip names and ensure numeric
  alpha <- as.numeric(alpha)
  beta <- as.numeric(beta)
  
  ## Guard against boundary values
  eps <- 1e-10
  alpha <- max(eps, alpha)
  beta <- max(eps, beta)
  
  persistence <- alpha + beta
  persistence <- max(eps, min(1 - eps, persistence))
  
  ## psi = logit(persistence) = log(p / (1-p))
  psi <- log(persistence / (1 - persistence))
  
  ## phi = log(alpha / beta)
  phi <- log(alpha / beta)
  
  return(c(psi = psi, phi = phi))
}


#' @title Inverse Reparameterization for DCC(1,1)
#' @description Transform unconstrained (psi, phi) back to constrained (alpha, beta).
#'   Stationarity is guaranteed for any finite (psi, phi).
#' @param psi Unconstrained persistence parameter (can be any real number)
#' @param phi Unconstrained ratio parameter (can be any real number)
#' @return Named vector c(alpha, beta)
#' @keywords internal
dcc11_from_unconstrained <- function(psi, phi) {
  ## Strip names and ensure numeric scalars
  psi <- as.numeric(psi)[1]
  phi <- as.numeric(phi)[1]
  
  ## persistence = sigmoid(psi) = 1 / (1 + exp(-psi))
  ## Use stable computation to avoid overflow
  if (psi > 0) {
    persistence <- 1 / (1 + exp(-psi))
  } else {
    exp_psi <- exp(psi)
    persistence <- exp_psi / (1 + exp_psi)
  }
  
  ## ratio = exp(phi) = alpha / beta
  ## alpha = persistence * ratio / (1 + ratio)
  ## beta = persistence / (1 + ratio)
  ## Use stable computation
  if (phi > 0) {
    exp_phi <- exp(phi)
    alpha <- persistence * exp_phi / (1 + exp_phi)
    beta <- persistence / (1 + exp_phi)
  } else {
    exp_neg_phi <- exp(-phi)
    alpha <- persistence / (1 + exp_neg_phi)
    beta <- persistence * exp_neg_phi / (1 + exp_neg_phi)
  }
  
  return(c(alpha = alpha, beta = beta))
}


#' @title Jacobian of Reparameterization
#' @description Compute d(alpha, beta) / d(psi, phi) for chain rule.
#' @param psi Unconstrained persistence parameter
#' @param phi Unconstrained ratio parameter
#' @return 2x2 matrix: rows = (alpha, beta), cols = (psi, phi)
#' @keywords internal
dcc11_reparam_jacobian <- function(psi, phi) {
  ## Strip names and ensure numeric scalars
  psi <- as.numeric(psi)[1]
  phi <- as.numeric(phi)[1]
  
  ## Get current values using stable computation
  if (psi > 0) {
    persistence <- 1 / (1 + exp(-psi))
  } else {
    exp_psi <- exp(psi)
    persistence <- exp_psi / (1 + exp_psi)
  }
  
  if (phi > 0) {
    exp_phi <- exp(phi)
    alpha <- persistence * exp_phi / (1 + exp_phi)
    beta <- persistence / (1 + exp_phi)
  } else {
    exp_neg_phi <- exp(-phi)
    alpha <- persistence / (1 + exp_neg_phi)
    beta <- persistence * exp_neg_phi / (1 + exp_neg_phi)
  }
  
  ## d(persistence)/d(psi) = persistence * (1 - persistence)
  d_pers_d_psi <- persistence * (1 - persistence)
  
  ## d(alpha)/d(psi) = (alpha / persistence) * d_pers_d_psi
  d_alpha_d_psi <- (alpha / persistence) * d_pers_d_psi
  
  ## d(beta)/d(psi) = (beta / persistence) * d_pers_d_psi
  d_beta_d_psi <- (beta / persistence) * d_pers_d_psi
  
  ## d(alpha)/d(phi) = alpha * beta / persistence
  d_alpha_d_phi <- alpha * beta / persistence
  
  ## d(beta)/d(phi) = -alpha * beta / persistence
  d_beta_d_phi <- -alpha * beta / persistence
  
  J <- matrix(c(d_alpha_d_psi, d_beta_d_psi,
                d_alpha_d_phi, d_beta_d_phi),
              nrow = 2, ncol = 2, byrow = FALSE)
  rownames(J) <- c("alpha", "beta")
  colnames(J) <- c("psi", "phi")
  
  return(J)
}


## SECTION 3: DCC Recursion ====================================================

#' @title Make Matrix Positive Definite (Safe Version)
#' @description Ensures a matrix is positive definite by adjusting eigenvalues.
#' @param R Symmetric matrix
#' @param min_eig Minimum eigenvalue (default 1e-8)
#' @return Positive definite matrix, or NULL if input contains NA/Inf
#' @keywords internal
make_pd_safe <- function(R, min_eig = 1e-8) {
  ## Check for NA/Inf values first
  if (any(!is.finite(R))) {
    return(NULL)
  }
  
  ## Ensure symmetry
  R <- (R + t(R)) / 2
  
  ## Check eigenvalues
  eig <- eigen(R, symmetric = TRUE)
  
  if (all(eig$values > min_eig)) {
    return(R)
  }
  
  ## Adjust small/negative eigenvalues
  eig$values[eig$values < min_eig] <- min_eig
  
  ## Reconstruct
  R_new <- eig$vectors %*% diag(eig$values, nrow = length(eig$values)) %*% t(eig$vectors)
  
  ## Ensure unit diagonal (for correlation matrix)
  d <- sqrt(diag(R_new))
  if (all(d > 0)) {
    R_new <- diag(1/d) %*% R_new %*% diag(1/d)
  }
  
  return(R_new)
}


#' @title DCC(1,1) Q-Matrix Recursion with Gradient Storage
#' @description Compute Q_t matrices and store gradient information for backpropagation.
#' @param z T x k matrix of standardized residuals
#' @param alpha DCC alpha parameter
#' @param beta DCC beta parameter
#' @param Qbar k x k unconditional covariance matrix
#' @return List with Q, R, R_inv, log_det_R, dQ_dalpha, dQ_dbeta, success
#' @keywords internal
dcc11_recursion_with_grad <- function(z, alpha, beta, Qbar) {
  ## Ensure alpha and beta are plain numeric scalars
  alpha <- as.numeric(alpha)[1]
  beta <- as.numeric(beta)[1]
  
  T_obs <- nrow(z)
  k <- ncol(z)
  
  ## Pre-allocate arrays
  Q <- array(0, dim = c(k, k, T_obs))
  R <- array(0, dim = c(k, k, T_obs))
  R_inv <- array(0, dim = c(k, k, T_obs))
  log_det_R <- numeric(T_obs)
  
  ## Gradient storage: dQ_t/dalpha and dQ_t/dbeta
  dQ_dalpha <- array(0, dim = c(k, k, T_obs))
  dQ_dbeta <- array(0, dim = c(k, k, T_obs))
  
  ## Initialize Q_1 = Qbar
  Q[,,1] <- Qbar
  
  ## Normalize to get R_1
  d1 <- sqrt(diag(Q[,,1]))
  d1[d1 <= 0] <- 1e-6
  D1_inv <- diag(1 / d1, k)
  R[,,1] <- D1_inv %*% Q[,,1] %*% D1_inv
  
  ## Ensure PD and compute inverse
  R_1_pd <- make_pd_safe(R[,,1])
  if (is.null(R_1_pd)) {
    R_1_pd <- diag(1, k)  ## Fallback to identity
  }
  
  chol_R1 <- tryCatch(chol(R_1_pd), error = function(e) NULL)
  if (is.null(chol_R1)) {
    R_1_pd <- R_1_pd + diag(1e-6, k)
    chol_R1 <- chol(R_1_pd)
  }
  
  R[,,1] <- R_1_pd
  R_inv[,,1] <- chol2inv(chol_R1)
  log_det_R[1] <- 2 * sum(log(diag(chol_R1)))
  
  ## dQ_1/dalpha = dQ_1/dbeta = 0 (Q_1 = Qbar is constant)
  ## Already initialized to 0
  
  ## Recursion for t = 2, ..., T
  for (t in 2:T_obs) {
    z_lag <- z[t-1, , drop = FALSE]
    zz <- t(z_lag) %*% z_lag  ## k x k outer product
    
    ## Q_t = (1 - alpha - beta) * Qbar + alpha * z_{t-1} z_{t-1}' + beta * Q_{t-1}
    Q[,,t] <- (1 - alpha - beta) * Qbar + alpha * zz + beta * Q[,,t-1]
    
    ## Check for numerical issues
    if (any(!is.finite(Q[,,t]))) {
      ## Return what we have with a failure flag
      return(list(
        Q = Q, R = R, R_inv = R_inv, log_det_R = log_det_R,
        dQ_dalpha = dQ_dalpha, dQ_dbeta = dQ_dbeta,
        success = FALSE, error_time = t
      ))
    }
    
    ## Gradient recursion
    ## dQ_t/dalpha = -Qbar + z_{t-1} z_{t-1}' + beta * dQ_{t-1}/dalpha
    dQ_dalpha[,,t] <- -Qbar + zz + beta * dQ_dalpha[,,t-1]
    
    ## dQ_t/dbeta = -Qbar + Q_{t-1} + beta * dQ_{t-1}/dbeta
    dQ_dbeta[,,t] <- -Qbar + Q[,,t-1] + beta * dQ_dbeta[,,t-1]
    
    ## Normalize Q_t to get R_t
    d_t <- sqrt(diag(Q[,,t]))
    d_t[!is.finite(d_t) | d_t <= 0] <- 1e-6
    
    D_t_inv <- diag(1 / d_t, k)
    R[,,t] <- D_t_inv %*% Q[,,t] %*% D_t_inv
    
    ## Make PD and compute inverse
    R_t_pd <- make_pd_safe(R[,,t])
    if (is.null(R_t_pd)) {
      R_t_pd <- diag(1, k)  ## Fallback to identity
    }
    
    chol_Rt <- tryCatch(chol(R_t_pd), error = function(e) NULL)
    if (is.null(chol_Rt)) {
      R_t_pd <- R_t_pd + diag(1e-6, k)
      chol_Rt <- chol(R_t_pd)
    }
    
    R[,,t] <- R_t_pd
    R_inv[,,t] <- chol2inv(chol_Rt)
    log_det_R[t] <- 2 * sum(log(diag(chol_Rt)))
  }
  
  return(list(
    Q = Q,
    R = R,
    R_inv = R_inv,
    log_det_R = log_det_R,
    dQ_dalpha = dQ_dalpha,
    dQ_dbeta = dQ_dbeta,
    success = TRUE
  ))
}


# See dcc_recursion() in dcc.R
#
# #' @title General DCC(p,q) Recursion (No Gradient)
# #' @description Compute Q and R matrices for arbitrary DCC order.
# #' @param std_resid Matrix of standardized residuals (T x k)
# #' @param Qbar Unconditional covariance matrix (k x k)
# #' @param alphas Numeric vector of alpha parameters (length q)
# #' @param betas Numeric vector of beta parameters (length p)
# #' @return List with success, Q, R, maxpq, and error info if failed
# #' @keywords internal
# dcc_recursion <- function(std_resid, Qbar, alphas, betas) {
#   T_obs <- nrow(std_resid)
#   k <- ncol(std_resid)
#   
#   q_order <- length(alphas)
#   p_order <- length(betas)
#   maxpq <- max(p_order, q_order, 1)
#   
#   persistence <- sum(alphas) + sum(betas)
#   
#   ## Initialize arrays
#   Q <- array(0, dim = c(k, k, T_obs))
#   R <- array(0, dim = c(k, k, T_obs))
#   
#   ## Initialize first maxpq observations with Qbar
#   for (t in 1:min(maxpq, T_obs)) {
#     Q[,,t] <- Qbar
#     Qbar_diag <- diag(Qbar)
#     if (any(Qbar_diag <= 0)) {
#       return(list(success = FALSE, error_type = "non_positive_Qbar_diagonal",
#                   error_time = 0, Q = Q, R = R, maxpq = maxpq))
#     }
#     Qbar_diag_inv_sqrt <- diag(1/sqrt(Qbar_diag), k)
#     R[,,t] <- Qbar_diag_inv_sqrt %*% Qbar %*% Qbar_diag_inv_sqrt
#   }
#   
#   if (T_obs <= maxpq) {
#     return(list(success = TRUE, Q = Q, R = R, maxpq = maxpq))
#   }
#   
#   ## Main recursion
#   for (t in (maxpq + 1):T_obs) {
#     Q_t <- Qbar * (1 - persistence)
#     
#     for (j in seq_along(alphas)) {
#       if (t - j >= 1) {
#         z_lag <- std_resid[t - j, , drop = FALSE]
#         Q_t <- Q_t + alphas[j] * (t(z_lag) %*% z_lag)
#       }
#     }
#     
#     for (j in seq_along(betas)) {
#       if (t - j >= 1) {
#         Q_t <- Q_t + betas[j] * Q[,,t - j]
#       }
#     }
#     
#     if (any(!is.finite(Q_t))) {
#       return(list(success = FALSE, error_type = "non_finite_Q",
#                   error_time = t, Q = Q, R = R, maxpq = maxpq))
#     }
#     
#     Q_diag <- diag(Q_t)
#     if (any(Q_diag <= 0)) {
#       return(list(success = FALSE, error_type = "non_positive_Q_diagonal",
#                   error_time = t, Q = Q, R = R, maxpq = maxpq))
#     }
#     
#     Q[,,t] <- Q_t
#     Q_diag_inv_sqrt <- diag(1/sqrt(Q_diag), k)
#     R[,,t] <- Q_diag_inv_sqrt %*% Q_t %*% Q_diag_inv_sqrt
#   }
#   
#   return(list(success = TRUE, Q = Q, R = R, maxpq = maxpq))
# }


## SECTION 4: Gradient Computation =============================================

#' @title Gradient of NLL w.r.t. R_t (MVN)
#' @description Compute gradient of negative log-likelihood contribution w.r.t. R_t.
#' @param z_t k-vector of standardized residuals at time t
#' @param R_t k x k correlation matrix
#' @param R_inv_t k x k inverse correlation matrix
#' @return k x k gradient matrix
#' @keywords internal
grad_nll_wrt_R_mvn <- function(z_t, R_t, R_inv_t) {
  u <- as.vector(R_inv_t %*% z_t)
  0.5 * (R_inv_t - outer(u, u))
}


#' @title Gradient of NLL w.r.t. R_t (MVT)
#' @description Compute gradient of negative log-likelihood contribution w.r.t. R_t.
#' @param z_t k-vector of standardized residuals at time t
#' @param R_t k x k correlation matrix
#' @param R_inv_t k x k inverse correlation matrix
#' @param shape Degrees of freedom
#' @return k x k gradient matrix
#' @keywords internal
grad_nll_wrt_R_mvt <- function(z_t, R_t, R_inv_t, shape) {
  k <- length(z_t)
  u <- as.vector(R_inv_t %*% z_t)
  q_t <- sum(z_t * u)
  kappa_t <- 1 + q_t / (shape - 2)
  weight <- (shape + k) / (2 * (shape - 2) * kappa_t)
  0.5 * R_inv_t - weight * outer(u, u)
}


#' @title Gradient of NLL w.r.t. Q_t (via R_t normalization)
#' @description Backpropagate gradient from R_t to Q_t through the normalization.
#' @param grad_R k x k gradient w.r.t. R_t
#' @param Q_t k x k Q matrix
#' @param R_t k x k correlation matrix
#' @return k x k gradient w.r.t. Q_t
#' @keywords internal
grad_R_to_Q <- function(grad_R, Q_t, R_t) {
  k <- nrow(Q_t)
  
  d <- sqrt(diag(Q_t))
  d[d <= 0] <- 1e-6
  D_inv <- diag(1 / d, k)
  
  ## Main term: D^{-1} grad_R D^{-1}
  grad_Q <- D_inv %*% grad_R %*% D_inv
  
  ## Diagonal correction
  for (i in 1:k) {
    correction <- sum(grad_R[i,] * R_t[i,]) / Q_t[i, i]
    grad_Q[i, i] <- grad_Q[i, i] - correction
  }
  
  ## Ensure symmetry
  (grad_Q + t(grad_Q)) / 2
}


#' @title DCC(1,1) Weighted NLL Gradient (Original Parameterization)
#' @description Compute gradient of weighted NLL w.r.t. (alpha, beta).
#' @param alpha DCC alpha parameter
#' @param beta DCC beta parameter
#' @param z T x k matrix of standardized residuals
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional covariance
#' @param distribution "mvn" or "mvt"
#' @param shape Degrees of freedom (MVT only)
#' @return Named vector c(alpha, beta)
#' @keywords internal
dcc11_gradient_original <- function(alpha, beta, z, weights, Qbar,
                                    distribution = "mvn", shape = NULL) {
  ## Ensure plain numeric scalars
  alpha <- as.numeric(alpha)[1]
  beta <- as.numeric(beta)[1]
  
  T_obs <- nrow(z)
  
  ## Forward pass
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
    
    ## Gradient of NLL w.r.t. R_t
    if (distribution == "mvn") {
      grad_R <- grad_nll_wrt_R_mvn(z_t, R_t, R_inv_t)
    } else {
      grad_R <- grad_nll_wrt_R_mvt(z_t, R_t, R_inv_t, shape)
    }
    
    ## Backpropagate to Q_t
    grad_Q <- grad_R_to_Q(grad_R, Q_t, R_t)
    
    ## Chain rule using Frobenius inner product
    grad_alpha <- grad_alpha + w_t * sum(grad_Q * fwd$dQ_dalpha[,,t])
    grad_beta <- grad_beta + w_t * sum(grad_Q * fwd$dQ_dbeta[,,t])
  }
  
  c(alpha = grad_alpha, beta = grad_beta)
}


#' @title DCC(1,1) Weighted NLL Gradient (Reparameterized)
#' @description Compute gradient of weighted NLL w.r.t. (psi, phi).
#' @param psi Unconstrained persistence parameter
#' @param phi Unconstrained ratio parameter
#' @param z T x k matrix of standardized residuals
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional covariance
#' @param distribution "mvn" or "mvt"
#' @param shape Degrees of freedom (MVT only)
#' @return Named vector c(psi, phi)
#' @keywords internal
dcc11_gradient_reparam <- function(psi, phi, z, weights, Qbar,
                                   distribution = "mvn", shape = NULL) {
  ## Ensure plain numeric scalars
  psi <- as.numeric(psi)[1]
  phi <- as.numeric(phi)[1]
  
  ## Convert to original parameters
  orig <- dcc11_from_unconstrained(psi, phi)
  
  ## Get gradient in original space
  grad_orig <- dcc11_gradient_original(orig["alpha"], orig["beta"], z, weights, Qbar,
                                       distribution, shape)
  
  ## Check for failure
  if (any(is.na(grad_orig))) {
    return(c(psi = NA_real_, phi = NA_real_))
  }
  
  ## Get Jacobian and apply chain rule
  J <- dcc11_reparam_jacobian(psi, phi)
  grad_reparam <- as.vector(t(J) %*% grad_orig)
  names(grad_reparam) <- c("psi", "phi")
  
  grad_reparam
}


#' @title Gradient of NLL w.r.t. Shape Parameter (MVT)
#' @description Compute gradient of weighted NLL w.r.t. degrees of freedom.
#'   
#'   The MVT log-likelihood (per observation) is:
#'   ll_t = lgamma((nu+k)/2) - lgamma(nu/2) - (k/2)*log(pi*(nu-2)) 
#'          - 0.5*log|R_t| - ((nu+k)/2)*log(kappa_t)
#'   
#'   where kappa_t = 1 + q_t/(nu-2) and q_t = z_t' R_inv z_t.
#'   
#' @param shape Degrees of freedom (nu)
#' @param z T x k matrix of standardized residuals
#' @param weights T-vector of observation weights
#' @param R_inv Array of inverse correlation matrices
#' @return Scalar gradient (unnamed)
#' @keywords internal
dcc_gradient_shape <- function(shape, z, weights, R_inv) {
  T_obs <- nrow(z)
  k <- ncol(z)
  nu <- shape
  
  grad_shape <- 0
  
  for (t in 1:T_obs) {
    z_t <- z[t, ]
    R_inv_t <- R_inv[,,t]
    
    q_t <- as.numeric(t(z_t) %*% R_inv_t %*% z_t)
    kappa_t <- 1 + q_t / (nu - 2)
    
    ## Gradient of LL w.r.t. nu:
    ## d(ll_t)/d(nu) = 0.5*psi((nu+k)/2) - 0.5*psi(nu/2) - k/(2*(nu-2))
    ##                 - 0.5*log(kappa_t) - ((nu+k)/2) * d(log(kappa_t))/d(nu)
    ##
    ## where d(log(kappa_t))/d(nu) = -q_t / ((nu-2)^2 * kappa_t)
    
    d_lgamma_term <- 0.5 * (digamma((nu + k) / 2) - digamma(nu / 2))
    d_const_term <- -k / (2 * (nu - 2))
    d_log_kappa <- -q_t / ((nu - 2)^2 * kappa_t)
    d_kappa_term <- -0.5 * log(kappa_t) - ((nu + k) / 2) * d_log_kappa
    
    ## d(ll_t)/d(nu)
    dll_dnu <- d_lgamma_term + d_const_term + d_kappa_term
    
    ## Gradient of NLL is negative of gradient of LL
    grad_shape <- grad_shape - weights[t] * dll_dnu
  }
  
  as.numeric(grad_shape)
}


## SECTION 5: Main Interface Functions =========================================

#' @title Full DCC(1,1) Gradient Function for Optimizer
#' @description Compute gradient of weighted NLL w.r.t. all DCC parameters.
#'   This is the main interface for use with optim().
#' @param params Named vector of parameters to optimize
#' @param std_resid T x k matrix of standardized residuals
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional covariance
#' @param distribution "mvn" or "mvt"
#' @param use_reparam Logical: use reparameterized (psi, phi) space?
#' @return Vector of gradients (same length and order as params)
#' @export
dcc11_gradient <- function(params, std_resid, weights, Qbar,
                           distribution = "mvn", use_reparam = TRUE) {
  
  if (use_reparam) {
    psi <- as.numeric(params[1])
    phi <- as.numeric(params[2])
    
    grad_dcc <- dcc11_gradient_reparam(psi, phi, std_resid, weights, Qbar,
                                       distribution,
                                       shape = if (distribution == "mvt") as.numeric(params[3]) else NULL)
    
    if (distribution == "mvt") {
      shape <- as.numeric(params[3])
      orig <- dcc11_from_unconstrained(psi, phi)
      fwd <- dcc11_recursion_with_grad(std_resid, orig["alpha"], orig["beta"], Qbar)
      grad_shape <- dcc_gradient_shape(shape, std_resid, weights, fwd$R_inv)
      result <- c(grad_dcc, grad_shape)
      names(result) <- c("psi", "phi", "shape")
      return(result)
    } else {
      return(grad_dcc)
    }
    
  } else {
    alpha <- as.numeric(params[1])
    beta <- as.numeric(params[2])
    
    grad_dcc <- dcc11_gradient_original(alpha, beta, std_resid, weights, Qbar,
                                        distribution,
                                        shape = if (distribution == "mvt") as.numeric(params[3]) else NULL)
    
    if (distribution == "mvt") {
      shape <- as.numeric(params[3])
      fwd <- dcc11_recursion_with_grad(std_resid, alpha, beta, Qbar)
      grad_shape <- dcc_gradient_shape(shape, std_resid, weights, fwd$R_inv)
      result <- c(grad_dcc, grad_shape)
      names(result) <- c("alpha", "beta", "shape")
      return(result)
    } else {
      return(grad_dcc)
    }
  }
}


#' @title DCC(1,1) Weighted Negative Log-Likelihood
#' @description Compute weighted NLL for DCC(1,1) model.
#' @param params Named vector: c(psi, phi, [shape]) or c(alpha, beta, [shape])
#' @param std_resid T x k matrix of standardized residuals
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional covariance
#' @param distribution "mvn" or "mvt"
#' @param use_reparam Logical: use reparameterized space?
#' @return Scalar negative log-likelihood
#' @export
dcc11_nll <- function(params, std_resid, weights, Qbar,
                      distribution = "mvn", use_reparam = TRUE) {
  
  T_obs <- nrow(std_resid)
  k <- ncol(std_resid)
  
  ## Get alpha, beta - ensure plain numeric scalars
  if (use_reparam) {
    psi <- as.numeric(params[1])
    phi <- as.numeric(params[2])
    orig <- dcc11_from_unconstrained(psi, phi)
    alpha <- as.numeric(orig["alpha"])
    beta <- as.numeric(orig["beta"])
  } else {
    alpha <- as.numeric(params[1])
    beta <- as.numeric(params[2])
  }
  
  shape <- if (distribution == "mvt") as.numeric(params[length(params)]) else NULL
  
  ## Check constraints (should not trigger with reparameterization)
  if (!is.finite(alpha) || !is.finite(beta) || 
      alpha + beta >= 1 || alpha < 0 || beta < 0) {
    return(1e10)
  }
  
  if (distribution == "mvt" && (!is.null(shape) && shape <= 2)) {
    return(1e10)
  }
  
  ## Run DCC recursion
  fwd <- dcc11_recursion_with_grad(std_resid, alpha, beta, Qbar)
  
  ## Check for recursion failure
  if (!isTRUE(fwd$success)) {
    return(1e10)
  }
  
  ## Compute weighted NLL
  nll <- 0
  
  for (t in 1:T_obs) {
    z_t <- std_resid[t, ]
    R_inv_t <- fwd$R_inv[,,t]
    log_det_R_t <- fwd$log_det_R[t]
    
    q_t <- as.numeric(t(z_t) %*% R_inv_t %*% z_t)
    
    if (distribution == "mvn") {
      ll_t <- -0.5 * (log_det_R_t + q_t - sum(z_t^2))
    } else {
      const_term <- lgamma((shape + k) / 2) - lgamma(shape / 2) -
        (k / 2) * log(pi * (shape - 2))
      kappa_t <- 1 + q_t / (shape - 2)
      ll_t <- const_term - 0.5 * log_det_R_t - ((shape + k) / 2) * log(kappa_t)
    }
    
    nll <- nll - weights[t] * ll_t
  }
  
  nll
}


#' @title Numerical Gradient (for Verification)
#' @description Compute numerical gradient using central differences.
#' @param fn Objective function
#' @param params Parameter vector
#' @param eps Step size (default 1e-6)
#' @param ... Additional arguments to fn
#' @return Vector of numerical gradients
#' @keywords internal
numerical_gradient <- function(fn, params, eps = 1e-6, ...) {
  n <- length(params)
  grad <- numeric(n)
  
  for (i in 1:n) {
    params_plus <- params
    params_minus <- params
    params_plus[i] <- params[i] + eps
    params_minus[i] <- params[i] - eps
    grad[i] <- (fn(params_plus, ...) - fn(params_minus, ...)) / (2 * eps)
  }
  
  names(grad) <- names(params)
  grad
}