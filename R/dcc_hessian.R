## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## DCC(1,1) Hessian and Standard Error Computation
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
##
## This file implements Hessian computation and standard errors for DCC(1,1).
##
## Methods supported:
##   1. Numerical Hessian via finite differences (primary)
##   2. Analytical Hessian via DCC recursion derivatives
##   3. Observed Information Matrix (negative Hessian at MLE)
##   4. Sandwich estimator for robust standard errors
##
## IMPORTANT: Inference Recommendations
## ─────────────────────────────────────
## Validation studies show that the DCC likelihood surface is highly anisotropic:
##   - Curvature in alpha direction is ~8x steeper than beta direction
##   - Hessian-based SEs for beta are severely underestimated (by ~5-6x)
##   - This is a fundamental property of DCC with high persistence, not a bug
##
## Recommended inference methods:
##   - Alpha: Hessian-based SE is acceptable (ratio to true SD ~0.8)
##   - Beta:  Use Bootstrap SE or Profile Likelihood CI
##   - Persistence (α+β): May be better identified than individual parameters
##
## See dcc_inference_guide.Rmd vignette for detailed findings and recommendations.
## See dcc_inference.R for bootstrap and profile likelihood implementations.
##
## Robust handling:
##   - Boundary estimates (alpha → 0 or beta → 0): Returns NA with warning
##   - Non-PD Hessian: Uses pseudo-inverse with warning
##   - Constant correlation: Returns NULL (no DCC parameters to estimate)
##
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


## SECTION 0: Robust SE Computation Wrapper ====================================

#' @title Compute DCC Standard Errors with Robust Handling
#' @description High-level wrapper that computes standard errors for DCC parameters
#'   with appropriate handling of edge cases (boundary estimates, constant correlation,
#'   non-positive-definite Hessians).
#'   
#' @param dcc_result Result from estimate_dcc_parameters_weighted() or similar
#' @param std_resid T x k matrix of standardized residuals
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional correlation matrix
#' @param distribution "mvn" or "mvt"
#' @param boundary_threshold Threshold for detecting boundary estimates (default 1e-4)
#' @param method "hessian" or "sandwich"
#' @return List with:
#'   \item{se}{Named vector of standard errors (NA for invalid/boundary)}
#'   \item{vcov}{Variance-covariance matrix}
#'   \item{valid}{Logical: are the SEs valid/reliable?}
#'   \item{reason}{Character: reason if SEs are invalid}
#'   \item{correlation_type}{"dynamic" or "constant"}
#' @export
compute_dcc_standard_errors_robust <- function(
    dcc_result,
    std_resid,
    weights,
    Qbar,
    distribution = "mvn",
    boundary_threshold = 1e-4,
    method = "hessian"
) {
  
  ## Initialize result structure
  result <- list(
    se = NULL,
    vcov = NULL,
    valid = FALSE,
    reason = NULL,
    correlation_type = NULL,
    info = NULL,
    info_eigenvalues = NULL
  )
  

  ## Case 1: Constant correlation - no DCC parameters to estimate --------------
  
  correlation_type <- dcc_result$correlation_type %||% 
    dcc_result$coefficients$correlation_type %||%
    "dynamic"
  
  result$correlation_type <- correlation_type
  
  if (correlation_type == "constant") {
    result$reason <- "constant_correlation"
    result$valid <- TRUE  ## Valid in the sense that there's nothing to estimate
    return(result)
  }
  

  ## Case 2: Extract DCC parameters --------------------------------------------
  
  dcc_pars <- dcc_result$dcc_pars %||% dcc_result$coefficients$dcc_pars
  
  if (is.null(dcc_pars) || length(dcc_pars) == 0) {
    result$reason <- "no_dcc_parameters"
    return(result)
  }
  
  ## Extract alpha and beta
  alpha <- dcc_pars$alpha_1 %||% dcc_pars[["alpha_1"]]
  beta <- dcc_pars$beta_1 %||% dcc_pars[["beta_1"]]
  
  if (is.null(alpha) || is.null(beta)) {
    result$reason <- "missing_alpha_or_beta"
    return(result)
  }
  

  ## Case 3: Check for boundary estimates --------------------------------------
  
  at_lower_boundary <- alpha < boundary_threshold || beta < boundary_threshold
  at_upper_boundary <- (alpha + beta) > (1 - boundary_threshold)
  
  if (at_lower_boundary) {
    result$reason <- sprintf("boundary_lower (alpha=%.2e, beta=%.2e)", alpha, beta)
    result$se <- c(alpha = NA_real_, beta = NA_real_)
    if (distribution == "mvt") {
      result$se <- c(result$se, shape = NA_real_)
    }
    return(result)
  }
  
  if (at_upper_boundary) {
    result$reason <- sprintf("boundary_upper (persistence=%.6f)", alpha + beta)
    result$se <- c(alpha = NA_real_, beta = NA_real_)
    if (distribution == "mvt") {
      result$se <- c(result$se, shape = NA_real_)
    }
    return(result)
  }
  

  ## Case 4: Compute standard errors -------------------------------------------
  
  ## Build parameter vector
  params <- c(alpha = alpha, beta = beta)
  if (distribution == "mvt") {
    shape <- dcc_result$dist_pars$shape %||% 
      dcc_result$coefficients$dist_pars$shape %||% 
      8.0
    params <- c(params, shape = shape)
  }
  
  ## Compute SEs
  se_result <- tryCatch({
    if (method == "sandwich") {
      dcc11_sandwich_se(
        params = params,
        std_resid = std_resid,
        weights = weights,
        Qbar = Qbar,
        distribution = distribution,
        use_reparam = FALSE
      )
    } else {
      dcc11_standard_errors(
        params = params,
        std_resid = std_resid,
        weights = weights,
        Qbar = Qbar,
        distribution = distribution,
        use_reparam = FALSE
      )
    }
  }, error = function(e) {
    list(se = NULL, vcov = NULL, info = NULL, error = e$message)
  })
  
  ## Check for computation error
  if (!is.null(se_result$error)) {
    result$reason <- paste("computation_error:", se_result$error)
    result$se <- c(alpha = NA_real_, beta = NA_real_)
    if (distribution == "mvt") {
      result$se <- c(result$se, shape = NA_real_)
    }
    return(result)
  }
  

  ## Case 5: Validate results --------------------------------------------------
  
  result$se <- se_result$se
  result$vcov <- se_result$vcov
  result$info <- se_result$info
  
  ## Check eigenvalues of information matrix
  if (!is.null(se_result$info)) {
    eig <- eigen(se_result$info, symmetric = TRUE, only.values = TRUE)$values
    result$info_eigenvalues <- eig
    
    if (any(eig <= 0)) {
      result$reason <- sprintf("non_pd_information (min_eig=%.2e)", min(eig))
      result$valid <- FALSE
      return(result)
    }
  }
  
  ## Check for NA or non-positive SEs
  if (any(is.na(se_result$se))) {
    result$reason <- "na_standard_errors"
    result$valid <- FALSE
    return(result)
  }
  
  if (any(se_result$se <= 0)) {
    result$reason <- "non_positive_standard_errors"
    result$valid <- FALSE
    return(result)
  }
  
  ## All checks passed
  result$valid <- TRUE
  result$reason <- "ok"
  
  return(result)
}


## Helper for NULL coalescing
`%||%` <- function(a, b) if (is.null(a)) b else a


## SECTION 1: Numerical Hessian Computation ====================================

#' @title Compute Numerical Hessian via Finite Differences
#' @description Compute the Hessian matrix using central finite differences.
#'   Uses Richardson extrapolation for improved accuracy when available.
#' @param fn Objective function (returns scalar)
#' @param params Parameter vector at which to evaluate Hessian
#' @param eps Step size for finite differences (default 1e-5)
#' @param ... Additional arguments passed to fn
#' @return Square matrix of second derivatives
#' @keywords internal
numerical_hessian <- function(fn, params, eps = 1e-5, ...) {
  n <- length(params)
  H <- matrix(0, n, n)
  
  f0 <- fn(params, ...)
  
  for (i in 1:n) {
    for (j in i:n) {
      if (i == j) {
        ## Diagonal: use standard second derivative formula
        ## f''(x) ≈ [f(x+h) - 2f(x) + f(x-h)] / h^2
        params_plus <- params_minus <- params
        params_plus[i] <- params[i] + eps
        params_minus[i] <- params[i] - eps
        
        f_plus <- fn(params_plus, ...)
        f_minus <- fn(params_minus, ...)
        
        H[i, i] <- (f_plus - 2 * f0 + f_minus) / (eps^2)
      } else {
        ## Off-diagonal: use mixed partial derivative formula
        ## f_xy ≈ [f(x+h,y+h) - f(x+h,y-h) - f(x-h,y+h) + f(x-h,y-h)] / (4h^2)
        params_pp <- params_pm <- params_mp <- params_mm <- params
        
        params_pp[i] <- params[i] + eps
        params_pp[j] <- params[j] + eps
        
        params_pm[i] <- params[i] + eps
        params_pm[j] <- params[j] - eps
        
        params_mp[i] <- params[i] - eps
        params_mp[j] <- params[j] + eps
        
        params_mm[i] <- params[i] - eps
        params_mm[j] <- params[j] - eps
        
        f_pp <- fn(params_pp, ...)
        f_pm <- fn(params_pm, ...)
        f_mp <- fn(params_mp, ...)
        f_mm <- fn(params_mm, ...)
        
        H[i, j] <- H[j, i] <- (f_pp - f_pm - f_mp + f_mm) / (4 * eps^2)
      }
    }
  }
  
  ## Ensure perfect symmetry
  H <- (H + t(H)) / 2
  
  ## Add names if params are named
  if (!is.null(names(params))) {
    rownames(H) <- colnames(H) <- names(params)
  }
  
  return(H)
}


#' @title Compute Numerical Hessian with Richardson Extrapolation
#' @description More accurate Hessian using Richardson extrapolation.
#'   Combines estimates at different step sizes for higher-order accuracy.
#' @param fn Objective function
#' @param params Parameter vector
#' @param eps Base step size (default 1e-4)
#' @param r Reduction factor for Richardson extrapolation (default 2)
#' @param ... Additional arguments to fn
#' @return Hessian matrix with improved accuracy
#' @keywords internal
numerical_hessian_richardson <- function(fn, params, eps = 1e-4, r = 2, ...) {
  ## Compute Hessian at two step sizes
  H1 <- numerical_hessian(fn, params, eps = eps, ...)
  H2 <- numerical_hessian(fn, params, eps = eps / r, ...)
  
  ## Richardson extrapolation: H = (r^2 * H2 - H1) / (r^2 - 1)
  ## This eliminates the leading error term
  H_extrap <- (r^2 * H2 - H1) / (r^2 - 1)
  
  return(H_extrap)
}


## SECTION 2: Observed Information Matrix ======================================

#' @title Compute Observed Information Matrix for DCC(1,1)
#' @description Compute the observed Fisher information matrix, which is the
#'   negative Hessian of the log-likelihood evaluated at the MLE.
#' @param params MLE parameter estimates (alpha, beta) or (psi, phi)
#' @param std_resid T x k matrix of standardized residuals
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional covariance matrix
#' @param distribution "mvn" or "mvt"
#' @param use_reparam Logical: parameters in (psi, phi) space?
#' @param eps Step size for numerical differentiation (default 1e-5)
#' @return Observed information matrix (positive semi-definite if at MLE)
#' @export
dcc11_observed_information <- function(params, std_resid, weights, Qbar,
                                       distribution = "mvn",
                                       use_reparam = FALSE,
                                       eps = 1e-5) {
  
  ## Compute Hessian of NLL (using Richardson for accuracy)
  H <- numerical_hessian_richardson(
    fn = dcc11_nll,
    params = params,
    eps = eps,
    std_resid = std_resid,
    weights = weights,
    Qbar = Qbar,
    distribution = distribution,
    use_reparam = use_reparam
  )
  
  ## Information = negative Hessian of LL = Hessian of NLL
  ## (since we minimize NLL, Hessian of NLL at minimum should be positive definite)
  I_obs <- H
  
  return(I_obs)
}


#' @title Compute Hessian of DCC(1,1) NLL
#' @description Compute the Hessian matrix of the negative log-likelihood.
#'   This is a convenience wrapper that returns both the Hessian and derived
#'   quantities (standard errors, eigenvalues).
#' @param params MLE parameter estimates c(alpha, beta) or c(alpha, beta, shape)
#' @param std_resid T x k matrix of standardized residuals
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional correlation matrix
#' @param distribution "mvn" or "mvt"
#' @param use_reparam Logical: parameters in reparameterized space?
#' @param method "numerical" (default) or "analytical"
#' @param eps Step size for numerical differentiation
#' @return List with:
#'   \item{hessian}{Hessian matrix of NLL}
#'   \item{info}{Observed information matrix (= Hessian for NLL)}
#'   \item{vcov}{Variance-covariance matrix (inverse of info)}
#'   \item{se}{Standard errors}
#'   \item{eigenvalues}{Eigenvalues of Hessian}
#'   \item{condition_number}{Condition number of Hessian}
#' @export
dcc11_hessian <- function(params, std_resid, weights, Qbar,
                          distribution = "mvn",
                          use_reparam = FALSE,
                          method = "numerical",
                          eps = 1e-5) {
  
  n_params <- length(params)
  
  ## Compute Hessian
  if (method == "analytical") {
    H <- dcc11_analytical_hessian(
      params = params,
      std_resid = std_resid,
      weights = weights,
      Qbar = Qbar,
      distribution = distribution
    )
  } else {
    H <- numerical_hessian_richardson(
      fn = dcc11_nll,
      params = params,
      eps = eps,
      std_resid = std_resid,
      weights = weights,
      Qbar = Qbar,
      distribution = distribution,
      use_reparam = use_reparam
    )
  }
  
  ## Set parameter names
  if (n_params == 2) {
    if (use_reparam) {
      param_names <- c("psi", "phi")
    } else {
      param_names <- c("alpha", "beta")
    }
  } else if (n_params == 3) {
    if (use_reparam) {
      param_names <- c("psi", "phi", "shape")
    } else {
      param_names <- c("alpha", "beta", "shape")
    }
  } else {
    param_names <- paste0("param", 1:n_params)
  }
  
  rownames(H) <- colnames(H) <- param_names
  
  ## Eigendecomposition
  eig <- eigen(H, symmetric = TRUE)
  
  ## Inverse for variance-covariance (with regularization if needed)
  vcov <- tryCatch({
    if (all(eig$values > 1e-10)) {
      solve(H)
    } else {
      ## Pseudo-inverse for near-singular case
      warning("Hessian is near-singular, using pseudo-inverse")
      eig_reg <- pmax(eig$values, 1e-10)
      eig$vectors %*% diag(1 / eig_reg) %*% t(eig$vectors)
    }
  }, error = function(e) {
    matrix(NA, n_params, n_params)
  })
  
  rownames(vcov) <- colnames(vcov) <- param_names
  
  ## Standard errors
  se <- sqrt(diag(vcov))
  names(se) <- param_names
  
  list(
    hessian = H,
    info = H,  ## For NLL, Hessian = observed information
    vcov = vcov,
    se = se,
    eigenvalues = eig$values,
    eigenvectors = eig$vectors,
    condition_number = max(eig$values) / min(pmax(eig$values, 1e-16)),
    method = method,
    params = params,
    param_names = param_names
  )
}


#' @title Analytical Hessian of DCC(1,1) NLL
#' @description Compute the analytical Hessian matrix of the DCC(1,1) negative
#'   log-likelihood. This requires second derivatives through the DCC recursion.
#' @param params Parameter vector c(alpha, beta)
#' @param std_resid T x k matrix of standardized residuals
#' @param weights T-vector of observation weights  
#' @param Qbar k x k unconditional correlation matrix
#' @param distribution "mvn" or "mvt"
#' @return k x k Hessian matrix
#' @keywords internal
dcc11_analytical_hessian <- function(params, std_resid, weights, Qbar,
                                     distribution = "mvn") {
  
  alpha <- params[1]
  beta <- params[2]
  
  n <- nrow(std_resid)
  k <- ncol(std_resid)
  
  ## Initialize storage for second derivatives
  ## We need d²L/dα², d²L/dβ², d²L/dαdβ
  
  ## Initialize Q and its derivatives
  Q <- Qbar
  dQ_dalpha <- matrix(0, k, k)
  dQ_dbeta <- matrix(0, k, k)
  
  ## Second derivatives of Q
  d2Q_dalpha2 <- matrix(0, k, k)
  d2Q_dbeta2 <- matrix(0, k, k)
  d2Q_dalphadbeta <- matrix(0, k, k)
  
  ## Hessian components
  H11 <- 0  ## d²NLL/dα²
  H22 <- 0  ## d²NLL/dβ²
  H12 <- 0  ## d²NLL/dαdβ
  
  for (t in 1:n) {
    z_t <- std_resid[t, , drop = FALSE]
    w_t <- weights[t]
    
    ## Normalize Q to get R
    q_diag <- diag(Q)
    q_diag_safe <- pmax(q_diag, 1e-10)
    d_inv <- 1 / sqrt(q_diag_safe)
    D_inv <- diag(d_inv, k)
    R <- D_inv %*% Q %*% D_inv
    
    ## Ensure R is valid correlation
    diag(R) <- 1
    R <- (R + t(R)) / 2
    
    ## Compute R inverse
    R_inv <- tryCatch(solve(R), error = function(e) {
      eig <- eigen(R, symmetric = TRUE)
      eig$vectors %*% diag(1 / pmax(eig$values, 1e-8)) %*% t(eig$vectors)
    })
    
    ## First derivatives of D_inv w.r.t. Q diagonal
    ## dD_inv/dq_ii = -0.5 * q_ii^(-3/2)
    
    ## First derivatives of R w.r.t. alpha and beta (through Q)
    dR_dalpha <- compute_dR_dQ(Q, dQ_dalpha, k)
    dR_dbeta <- compute_dR_dQ(Q, dQ_dbeta, k)
    
    ## Second derivatives of R
    d2R_dalpha2 <- compute_d2R_dQ2(Q, dQ_dalpha, dQ_dalpha, d2Q_dalpha2, k)
    d2R_dbeta2 <- compute_d2R_dQ2(Q, dQ_dbeta, dQ_dbeta, d2Q_dbeta2, k)
    d2R_dalphadbeta <- compute_d2R_dQ2(Q, dQ_dalpha, dQ_dbeta, d2Q_dalphadbeta, k)
    
    ## Likelihood contribution for MVN: 
    ## L_t = -0.5 * (log|R| + z'R^{-1}z)
    ## 
    ## First derivatives (we already compute these in gradient):
    ## dL/dθ = -0.5 * tr(R^{-1} dR/dθ) + 0.5 * z'R^{-1}(dR/dθ)R^{-1}z
    ##
    ## Second derivatives:
    ## d²L/dθ₁dθ₂ = 0.5 * tr(R^{-1}(dR/dθ₁)R^{-1}(dR/dθ₂))
    ##             - 0.5 * tr(R^{-1}(d²R/dθ₁dθ₂))
    ##             - z'R^{-1}(dR/dθ₁)R^{-1}(dR/dθ₂)R^{-1}z
    ##             + 0.5 * z'R^{-1}(d²R/dθ₁dθ₂)R^{-1}z
    
    z_vec <- as.vector(z_t)
    Rz <- R_inv %*% z_vec
    
    ## d²NLL/dα²
    RdRa <- R_inv %*% dR_dalpha
    RdRaRz <- RdRa %*% Rz
    
    term1_aa <- 0.5 * sum(diag(RdRa %*% RdRa))
    term2_aa <- -0.5 * sum(diag(R_inv %*% d2R_dalpha2))
    term3_aa <- -as.numeric(t(Rz) %*% dR_dalpha %*% RdRaRz)
    term4_aa <- 0.5 * as.numeric(t(Rz) %*% d2R_dalpha2 %*% Rz)
    
    H11 <- H11 + w_t * (term1_aa + term2_aa + term3_aa + term4_aa)
    
    ## d²NLL/dβ²
    RdRb <- R_inv %*% dR_dbeta
    RdRbRz <- RdRb %*% Rz
    
    term1_bb <- 0.5 * sum(diag(RdRb %*% RdRb))
    term2_bb <- -0.5 * sum(diag(R_inv %*% d2R_dbeta2))
    term3_bb <- -as.numeric(t(Rz) %*% dR_dbeta %*% RdRbRz)
    term4_bb <- 0.5 * as.numeric(t(Rz) %*% d2R_dbeta2 %*% Rz)
    
    H22 <- H22 + w_t * (term1_bb + term2_bb + term3_bb + term4_bb)
    
    ## d²NLL/dαdβ
    term1_ab <- 0.5 * sum(diag(RdRa %*% RdRb))
    term2_ab <- -0.5 * sum(diag(R_inv %*% d2R_dalphadbeta))
    term3_ab <- -as.numeric(t(Rz) %*% dR_dalpha %*% RdRbRz)
    term4_ab <- 0.5 * as.numeric(t(Rz) %*% d2R_dalphadbeta %*% Rz)
    
    H12 <- H12 + w_t * (term1_ab + term2_ab + term3_ab + term4_ab)
    
    ## Update Q and derivatives for next period
    zz <- t(z_t) %*% z_t
    
    ## Q_{t+1} = (1-α-β)Qbar + α*z_t*z_t' + β*Q_t
    ## dQ/dα = -Qbar + z_t*z_t' + β*(dQ/dα)
    ## dQ/dβ = -Qbar + Q_t + β*(dQ/dβ)
    
    ## Second derivatives:
    ## d²Q/dα² = β*(d²Q/dα²)
    ## d²Q/dβ² = 2*(dQ/dβ) + β*(d²Q/dβ²)  -- wait, need to be careful
    ## Actually: d²Q/dβ² = β*(d²Q/dβ²)  since dQ_t/dβ appears in next step
    ## d²Q/dαdβ = -I + (dQ/dα) + β*(d²Q/dαdβ)  -- also needs care
    
    ## Let me be more careful. Q_{t+1} depends on Q_t:
    ## dQ_{t+1}/dα = -Qbar + zz + β*(dQ_t/dα)
    ## d²Q_{t+1}/dα² = β*(d²Q_t/dα²)
    
    ## dQ_{t+1}/dβ = -Qbar + Q_t + β*(dQ_t/dβ)
    ## d²Q_{t+1}/dβ² = 2*(dQ_t/dβ) + β*(d²Q_t/dβ²)
    ## No wait: d/dβ[dQ_{t+1}/dβ] = d/dβ[-Qbar + Q_t + β*(dQ_t/dβ)]
    ##                            = dQ_t/dβ + dQ_t/dβ + β*(d²Q_t/dβ²)
    ##                            = 2*(dQ_t/dβ) + β*(d²Q_t/dβ²)
    
    ## d²Q_{t+1}/dαdβ = d/dα[-Qbar + Q_t + β*(dQ_t/dβ)]
    ##                = dQ_t/dα + β*(d²Q_t/dαdβ)
    ## No: d/dβ[dQ_{t+1}/dα] = d/dβ[-Qbar + zz + β*(dQ_t/dα)]
    ##                       = dQ_t/dα + β*(d²Q_t/dαdβ)
    
    d2Q_dalpha2_new <- beta * d2Q_dalpha2
    d2Q_dbeta2_new <- 2 * dQ_dbeta + beta * d2Q_dbeta2
    d2Q_dalphadbeta_new <- dQ_dalpha + beta * d2Q_dalphadbeta
    
    dQ_dalpha_new <- -Qbar + zz + beta * dQ_dalpha
    dQ_dbeta_new <- -Qbar + Q + beta * dQ_dbeta
    Q_new <- (1 - alpha - beta) * Qbar + alpha * zz + beta * Q
    
    Q <- Q_new
    dQ_dalpha <- dQ_dalpha_new
    dQ_dbeta <- dQ_dbeta_new
    d2Q_dalpha2 <- d2Q_dalpha2_new
    d2Q_dbeta2 <- d2Q_dbeta2_new
    d2Q_dalphadbeta <- d2Q_dalphadbeta_new
  }
  
  ## Assemble Hessian (for NLL, so signs are already correct)
  H <- matrix(c(H11, H12, H12, H22), 2, 2)
  
  return(H)
}


#' @title Compute dR/dQ (derivative of correlation w.r.t. Q matrix)
#' @keywords internal
compute_dR_dQ <- function(Q, dQ, k) {
  ## R = D^{-1} Q D^{-1} where D = diag(sqrt(diag(Q)))
  ## dR/dQ_ij requires chain rule through D
  
  q_diag <- diag(Q)
  q_diag_safe <- pmax(q_diag, 1e-10)
  d_inv <- 1 / sqrt(q_diag_safe)
  
  dq_diag <- diag(dQ)
  
  ## dD^{-1}/dq_ii = -0.5 * q_ii^{-3/2} * dq_ii
  dd_inv <- -0.5 * d_inv^3 * dq_diag
  
  ## dR = dD^{-1} Q D^{-1} + D^{-1} dQ D^{-1} + D^{-1} Q dD^{-1}
  D_inv <- diag(d_inv, k)
  dD_inv <- diag(dd_inv, k)
  
  dR <- dD_inv %*% Q %*% D_inv + D_inv %*% dQ %*% D_inv + D_inv %*% Q %*% dD_inv
  
  ## Force diagonal to be zero (correlation diagonal is always 1)
  diag(dR) <- 0
  
  ## Symmetrize
  dR <- (dR + t(dR)) / 2
  
  return(dR)
}


#' @title Compute d²R/dQ² (second derivative of correlation w.r.t. Q)
#' @keywords internal  
compute_d2R_dQ2 <- function(Q, dQ1, dQ2, d2Q, k) {
  ## Second derivative through the normalization
  ## This is quite involved...
  
  q_diag <- diag(Q)
  q_diag_safe <- pmax(q_diag, 1e-10)
  d_inv <- 1 / sqrt(q_diag_safe)
  
  dq1_diag <- diag(dQ1)
  dq2_diag <- diag(dQ2)
  d2q_diag <- diag(d2Q)
  
  ## First derivatives of d_inv
  dd1_inv <- -0.5 * d_inv^3 * dq1_diag
  dd2_inv <- -0.5 * d_inv^3 * dq2_diag
  
  ## Second derivative of d_inv
  ## d²(q^{-1/2})/dq² = (3/4) * q^{-5/2}
  ## d²d_inv/d(param1)d(param2) = (3/4) * d_inv^5 * dq1 * dq2 + (-1/2) * d_inv^3 * d2q
  d2d_inv <- 0.75 * d_inv^5 * dq1_diag * dq2_diag - 0.5 * d_inv^3 * d2q_diag
  
  D_inv <- diag(d_inv, k)
  dD1_inv <- diag(dd1_inv, k)
  dD2_inv <- diag(dd2_inv, k)
  d2D_inv <- diag(d2d_inv, k)
  
  ## d²R = d²D^{-1} Q D^{-1} + dD1^{-1} dQ2 D^{-1} + dD1^{-1} Q dD2^{-1}
  ##     + dD2^{-1} dQ1 D^{-1} + D^{-1} d2Q D^{-1} + D^{-1} dQ1 dD2^{-1}
  ##     + dD2^{-1} Q dD1^{-1} + D^{-1} dQ2 dD1^{-1} + D^{-1} Q d²D^{-1}
  
  d2R <- d2D_inv %*% Q %*% D_inv +
    dD1_inv %*% dQ2 %*% D_inv +
    dD1_inv %*% Q %*% dD2_inv +
    dD2_inv %*% dQ1 %*% D_inv +
    D_inv %*% d2Q %*% D_inv +
    D_inv %*% dQ1 %*% dD2_inv +
    dD2_inv %*% Q %*% dD1_inv +
    D_inv %*% dQ2 %*% dD1_inv +
    D_inv %*% Q %*% d2D_inv
  
  ## Force diagonal to zero
  diag(d2R) <- 0
  
  ## Symmetrize
  d2R <- (d2R + t(d2R)) / 2
  
  return(d2R)
}


## SECTION 3: Standard Error Computation =======================================

#' @title Compute Standard Errors for DCC(1,1) Parameters
#' @description Compute asymptotic standard errors from the observed information
#'   matrix. Returns NA for parameters where the information matrix is singular.
#' @param params MLE parameter estimates
#' @param std_resid T x k matrix of standardized residuals
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional covariance matrix
#' @param distribution "mvn" or "mvt"
#' @param use_reparam Logical: parameters in (psi, phi) space?
#' @param method "hessian" (default), "sandwich", or "both"
#' @return List with:
#'   \item{se}{Standard errors}
#'   \item{vcov}{Variance-covariance matrix}
#'   \item{info}{Observed information matrix}
#'   \item{info_inv}{Inverse information (variance-covariance)}
#'   \item{method}{Method used}
#' @export
dcc11_standard_errors <- function(params, std_resid, weights, Qbar,
                                  distribution = "mvn",
                                  use_reparam = FALSE,
                                  method = "hessian") {
  
  n_params <- length(params)
  param_names <- names(params)
  if (is.null(param_names)) {
    if (use_reparam) {
      param_names <- if (distribution == "mvt") c("psi", "phi", "shape") else c("psi", "phi")
    } else {
      param_names <- if (distribution == "mvt") c("alpha", "beta", "shape") else c("alpha", "beta")
    }
  }
  
  result <- list(
    se = rep(NA_real_, n_params),
    vcov = matrix(NA_real_, n_params, n_params),
    info = NULL,
    info_inv = NULL,
    method = method,
    params = params,
    param_names = param_names
  )
  names(result$se) <- param_names
  rownames(result$vcov) <- colnames(result$vcov) <- param_names
  
  ## Compute observed information
  I_obs <- tryCatch({
    dcc11_observed_information(
      params = params,
      std_resid = std_resid,
      weights = weights,
      Qbar = Qbar,
      distribution = distribution,
      use_reparam = use_reparam
    )
  }, error = function(e) {
    warning("Failed to compute observed information: ", e$message)
    return(NULL)
  })
  
  if (is.null(I_obs)) {
    return(result)
  }
  
  result$info <- I_obs
  rownames(result$info) <- colnames(result$info) <- param_names
  
  ## Check if information matrix is positive definite
  eig <- eigen(I_obs, symmetric = TRUE, only.values = TRUE)$values
  
  if (any(eig <= 0)) {
    warning("Information matrix is not positive definite. ",
            "Eigenvalues: ", paste(round(eig, 6), collapse = ", "),
            ". Standard errors may be unreliable.")
  }
  
  ## Invert to get variance-covariance matrix
  I_inv <- tryCatch({
    solve(I_obs)
  }, error = function(e) {
    warning("Information matrix is singular, using pseudo-inverse")
    ## Use Moore-Penrose pseudo-inverse
    MASS::ginv(I_obs)
  })
  
  result$info_inv <- I_inv
  result$vcov <- I_inv
  rownames(result$vcov) <- colnames(result$vcov) <- param_names
  
  ## Standard errors = sqrt of diagonal elements
  var_diag <- diag(I_inv)
  
  ## Handle any negative variances (shouldn't happen at true MLE)
  if (any(var_diag < 0)) {
    warning("Negative variance estimates detected. ",
            "This may indicate the optimizer did not find a true minimum.")
    var_diag[var_diag < 0] <- NA_real_
  }
  
  result$se <- sqrt(var_diag)
  names(result$se) <- param_names
  
  return(result)
}


#' @title Compute Sandwich (Robust) Standard Errors
#' @description Compute heteroskedasticity-robust standard errors using the
#'   sandwich estimator: Var(theta) = I^{-1} J I^{-1} where J is the outer
#'   product of gradients.
#' @param params MLE parameter estimates
#' @param std_resid T x k matrix of standardized residuals
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional covariance matrix
#' @param distribution "mvn" or "mvt"
#' @param use_reparam Logical: parameters in (psi, phi) space?
#' @return List with se, vcov, bread (I_inv), meat (J), sandwich
#' @export
dcc11_sandwich_se <- function(params, std_resid, weights, Qbar,
                              distribution = "mvn",
                              use_reparam = FALSE) {
  
  T_obs <- nrow(std_resid)
  n_params <- length(params)
  
  ## Get standard (Hessian-based) results first
  hessian_result <- dcc11_standard_errors(
    params = params,
    std_resid = std_resid,
    weights = weights,
    Qbar = Qbar,
    distribution = distribution,
    use_reparam = use_reparam,
    method = "hessian"
  )
  
  if (is.null(hessian_result$info_inv)) {
    warning("Cannot compute sandwich SE: information matrix is singular")
    return(hessian_result)
  }
  
  ## Compute the "meat" matrix: J = sum of outer products of score contributions
  ## We need per-observation gradients
  
  ## For efficiency, compute numerical gradient of each observation's contribution
  J <- matrix(0, n_params, n_params)
  eps <- 1e-6
  
  for (t in 1:T_obs) {
    ## Compute gradient for observation t
    ## This requires computing NLL for just observation t
    
    grad_t <- numerical_gradient_single_obs(
      params = params,
      std_resid = std_resid,
      weights = weights,
      Qbar = Qbar,
      t = t,
      distribution = distribution,
      use_reparam = use_reparam,
      eps = eps
    )
    
    ## Outer product contribution
    J <- J + weights[t]^2 * outer(grad_t, grad_t)
  }
  
  ## Sandwich: I^{-1} J I^{-1}
  bread <- hessian_result$info_inv
  sandwich <- bread %*% J %*% bread
  
  ## Standard errors
  var_diag <- diag(sandwich)
  var_diag[var_diag < 0] <- NA_real_
  se_sandwich <- sqrt(var_diag)
  
  result <- list(
    se = se_sandwich,
    vcov = sandwich,
    bread = bread,
    meat = J,
    info = hessian_result$info,
    info_inv = bread,
    method = "sandwich",
    params = params,
    param_names = hessian_result$param_names
  )
  
  names(result$se) <- hessian_result$param_names
  rownames(result$vcov) <- colnames(result$vcov) <- hessian_result$param_names
  
  return(result)
}


#' @title Numerical Gradient for Single Observation
#' @description Helper function to compute the gradient contribution from a
#'   single observation, needed for sandwich estimator.
#' @keywords internal
numerical_gradient_single_obs <- function(params, std_resid, weights, Qbar,
                                          t, distribution, use_reparam, eps = 1e-6) {
  n_params <- length(params)
  grad <- numeric(n_params)
  
  ## We need to compute d(nll_t)/d(params)
  ## This requires running the full recursion but only using observation t's likelihood
  
  for (i in 1:n_params) {
    params_plus <- params_minus <- params
    params_plus[i] <- params[i] + eps
    params_minus[i] <- params[i] - eps
    
    ## Compute contribution at t for params_plus and params_minus
    nll_plus <- dcc11_nll_single_obs(params_plus, std_resid, weights, Qbar, t,
                                     distribution, use_reparam)
    nll_minus <- dcc11_nll_single_obs(params_minus, std_resid, weights, Qbar, t,
                                      distribution, use_reparam)
    
    grad[i] <- (nll_plus - nll_minus) / (2 * eps)
  }
  
  return(grad)
}


#' @title NLL Contribution from Single Observation
#' @description Compute the NLL contribution from observation t.
#' @keywords internal
dcc11_nll_single_obs <- function(params, std_resid, weights, Qbar, t,
                                 distribution, use_reparam) {
  
  k <- ncol(std_resid)
  
  ## Get alpha, beta
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
  
  ## Check constraints
  if (!is.finite(alpha) || !is.finite(beta) || 
      alpha + beta >= 1 || alpha < 0 || beta < 0) {
    return(1e10)
  }
  
  ## Run full recursion to get R_t at time t
  fwd <- dcc11_recursion_with_grad(std_resid, alpha, beta, Qbar)
  
  if (!isTRUE(fwd$success)) {
    return(1e10)
  }
  
  ## Get observation t's likelihood contribution
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
  
  ## Return weighted NLL contribution
  -weights[t] * ll_t
}


## SECTION 4: Transformation Between Parameterizations =========================

#' @title Transform Standard Errors Between Parameterizations
#' @description Transform standard errors from (psi, phi) to (alpha, beta) or
#'   vice versa using the delta method.
#' @param se_result Result from dcc11_standard_errors()
#' @param to_reparam Logical: transform TO (psi, phi)? If FALSE, transform to (alpha, beta)
#' @return Updated se_result with transformed standard errors
#' @export
dcc11_transform_se <- function(se_result, to_reparam = FALSE) {
  
  params <- se_result$params
  vcov_orig <- se_result$vcov
  
  ## Get the Jacobian of the transformation
  if (to_reparam) {
    ## (alpha, beta) -> (psi, phi)
    ## Need d(psi, phi) / d(alpha, beta)
    alpha <- params[1]
    beta <- params[2]
    
    ## psi = logit(alpha + beta)
    ## phi = log(alpha / beta)
    persistence <- alpha + beta
    
    ## d(psi)/d(alpha) = 1 / (persistence * (1 - persistence))
    ## d(psi)/d(beta) = 1 / (persistence * (1 - persistence))
    d_psi_d_alpha <- 1 / (persistence * (1 - persistence))
    d_psi_d_beta <- 1 / (persistence * (1 - persistence))
    
    ## d(phi)/d(alpha) = 1 / alpha
    ## d(phi)/d(beta) = -1 / beta
    d_phi_d_alpha <- 1 / alpha
    d_phi_d_beta <- -1 / beta
    
    J <- matrix(c(d_psi_d_alpha, d_phi_d_alpha,
                  d_psi_d_beta, d_phi_d_beta),
                nrow = 2, ncol = 2, byrow = FALSE)
    
    new_param_names <- c("psi", "phi")
    
  } else {
    ## (psi, phi) -> (alpha, beta)
    ## Need d(alpha, beta) / d(psi, phi)
    ## This is the Jacobian we already have in dcc11_reparam_jacobian()
    psi <- params[1]
    phi <- params[2]
    
    J <- dcc11_reparam_jacobian(psi, phi)
    new_param_names <- c("alpha", "beta")
  }
  
  ## Handle shape parameter if present
  if (length(params) > 2) {
    ## Extend Jacobian to include shape (unchanged)
    n <- length(params)
    J_full <- diag(n)
    J_full[1:2, 1:2] <- J
    J <- J_full
    new_param_names <- c(new_param_names, "shape")
  }
  
  ## Delta method: Var(g(theta)) ≈ J Var(theta) J'
  vcov_new <- J %*% vcov_orig %*% t(J)
  
  ## Update result
  se_result$vcov <- vcov_new
  rownames(se_result$vcov) <- colnames(se_result$vcov) <- new_param_names
  
  se_result$se <- sqrt(pmax(0, diag(vcov_new)))
  names(se_result$se) <- new_param_names
  
  se_result$param_names <- new_param_names
  se_result$jacobian <- J
  
  return(se_result)
}


## SECTION 5: Confidence Intervals =============================================

#' @title Compute Confidence Intervals for DCC(1,1) Parameters
#' @description Compute asymptotic confidence intervals based on standard errors.
#' @param se_result Result from dcc11_standard_errors() or dcc11_sandwich_se()
#' @param level Confidence level (default 0.95)
#' @return Matrix with columns: estimate, se, lower, upper
#' @export
dcc11_confint <- function(se_result, level = 0.95) {
  
  alpha <- 1 - level
  z <- qnorm(1 - alpha / 2)
  
  params <- se_result$params
  se <- se_result$se
  
  lower <- params - z * se
  upper <- params + z * se
  
  ci <- cbind(
    estimate = params,
    se = se,
    lower = lower,
    upper = upper
  )
  
  rownames(ci) <- se_result$param_names
  
  attr(ci, "level") <- level
  attr(ci, "method") <- se_result$method
  
  return(ci)
}


#' @title Compute Profile Likelihood Confidence Intervals
#' @description Compute confidence intervals by inverting the likelihood ratio
#'   test. More accurate than Wald intervals, especially for bounded parameters.
#' @param params MLE estimates
#' @param std_resid Standardized residuals
#' @param weights Observation weights
#' @param Qbar Unconditional covariance
#' @param distribution "mvn" or "mvt"
#' @param use_reparam Logical: use reparameterized space?
#' @param param_idx Index of parameter for which to compute CI
#' @param level Confidence level (default 0.95)
#' @param n_grid Number of grid points for profile (default 50)
#' @return Named vector c(lower, upper)
#' @export
dcc11_profile_ci <- function(params, std_resid, weights, Qbar,
                             distribution = "mvn", use_reparam = FALSE,
                             param_idx = 1, level = 0.95, n_grid = 50) {
  
  ## Critical value for chi-squared(1)
  chi2_crit <- qchisq(level, df = 1)
  
  ## NLL at MLE
  nll_mle <- dcc11_nll(params, std_resid, weights, Qbar, distribution, use_reparam)
  
  ## Threshold for profile
  threshold <- nll_mle + chi2_crit / 2
  
  ## Search for lower bound
  param_mle <- params[param_idx]
  se_approx <- dcc11_standard_errors(params, std_resid, weights, Qbar,
                                     distribution, use_reparam)$se[param_idx]
  
  ## Grid search around MLE
  if (is.na(se_approx)) se_approx <- abs(param_mle) * 0.1 + 0.01
  
  search_range <- 4 * se_approx
  lower_grid <- seq(param_mle - search_range, param_mle, length.out = n_grid)
  upper_grid <- seq(param_mle, param_mle + search_range, length.out = n_grid)
  
  ## Profile likelihood on lower side
  profile_lower <- function(p) {
    params_test <- params
    params_test[param_idx] <- p
    dcc11_nll(params_test, std_resid, weights, Qbar, distribution, use_reparam)
  }
  
  nll_lower <- sapply(lower_grid, profile_lower)
  below_threshold <- which(nll_lower <= threshold)
  
  if (length(below_threshold) > 0) {
    lower_ci <- lower_grid[min(below_threshold)]
  } else {
    lower_ci <- lower_grid[1]
  }
  
  ## Profile likelihood on upper side
  nll_upper <- sapply(upper_grid, profile_lower)
  below_threshold <- which(nll_upper <= threshold)
  
  if (length(below_threshold) > 0) {
    upper_ci <- upper_grid[max(below_threshold)]
  } else {
    upper_ci <- upper_grid[n_grid]
  }
  
  return(c(lower = lower_ci, upper = upper_ci))
}


## SECTION 6: Summary and Printing Functions ===================================

#' @title Summarize DCC(1,1) Estimation Results
#' @description Create a comprehensive summary of DCC estimation with
#'   standard errors, confidence intervals, and diagnostic information.
#' @param params MLE estimates
#' @param std_resid Standardized residuals
#' @param weights Observation weights
#' @param Qbar Unconditional covariance
#' @param distribution "mvn" or "mvt"
#' @param use_reparam Logical: use reparameterized space?
#' @param level Confidence level (default 0.95)
#' @param se_method "hessian" or "sandwich" (default "hessian")
#' @return List with comprehensive estimation summary
#' @export
dcc11_estimation_summary <- function(params, std_resid, weights, Qbar,
                                     distribution = "mvn",
                                     use_reparam = FALSE,
                                     level = 0.95,
                                     se_method = "hessian") {
  
  T_obs <- nrow(std_resid)
  k <- ncol(std_resid)
  
  ## Compute NLL at MLE
  nll <- dcc11_nll(params, std_resid, weights, Qbar, distribution, use_reparam)
  ll <- -nll
  
  ## Compute standard errors
  if (se_method == "sandwich") {
    se_result <- dcc11_sandwich_se(params, std_resid, weights, Qbar,
                                   distribution, use_reparam)
  } else {
    se_result <- dcc11_standard_errors(params, std_resid, weights, Qbar,
                                       distribution, use_reparam)
  }
  
  ## Confidence intervals
  ci <- dcc11_confint(se_result, level = level)
  
  ## Get original parameters if reparameterized
  if (use_reparam) {
    orig <- dcc11_from_unconstrained(params[1], params[2])
    alpha <- orig["alpha"]
    beta <- orig["beta"]
    persistence <- alpha + beta
    
    ## Transform SE to original space
    se_orig <- dcc11_transform_se(se_result, to_reparam = FALSE)
    ci_orig <- dcc11_confint(se_orig, level = level)
  } else {
    alpha <- params[1]
    beta <- params[2]
    persistence <- alpha + beta
    se_orig <- se_result
    ci_orig <- ci
  }
  
  ## Information criteria (approximate)
  n_params <- length(params)
  aic <- 2 * nll + 2 * n_params
  bic <- 2 * nll + log(T_obs) * n_params
  
  ## Assemble result
  result <- list(
    ## Basic info
    n_obs = T_obs,
    n_series = k,
    n_params = n_params,
    distribution = distribution,
    
    ## Likelihood
    log_likelihood = ll,
    nll = nll,
    aic = aic,
    bic = bic,
    
    ## Parameters (original space)
    alpha = as.numeric(alpha),
    beta = as.numeric(beta),
    persistence = as.numeric(persistence),
    shape = if (distribution == "mvt") params[length(params)] else NULL,
    
    ## Standard errors
    se_method = se_method,
    se = se_orig$se,
    vcov = se_orig$vcov,
    
    ## Confidence intervals
    ci = ci_orig,
    ci_level = level,
    
    ## Information matrix diagnostics
    info_eigenvalues = eigen(se_result$info, symmetric = TRUE, only.values = TRUE)$values,
    info_condition = if (!is.null(se_result$info)) kappa(se_result$info) else NA
  )
  
  class(result) <- "dcc11_summary"
  return(result)
}


#' @title Print DCC(1,1) Estimation Summary
#' @export
print.dcc11_summary <- function(x, digits = 4, ...) {
  cat("\n")
  cat("===== DCC(1,1) Estimation Summary =====\n")
  cat("\n")
  
  cat("Data:\n")
  cat(sprintf("  Observations: %d\n", x$n_obs))
  cat(sprintf("  Series: %d\n", x$n_series))
  cat(sprintf("  Distribution: %s\n", toupper(x$distribution)))
  cat("\n")
  
  cat("Likelihood:\n")
  cat(sprintf("  Log-likelihood: %.4f\n", x$log_likelihood))
  cat(sprintf("  AIC: %.4f\n", x$aic))
  cat(sprintf("  BIC: %.4f\n", x$bic))
  cat("\n")
  
  cat(sprintf("Parameter Estimates (SE method: %s):\n", x$se_method))
  cat(sprintf("  alpha:       %.%df  (SE: %.%df)\n", digits, x$alpha, digits, x$se["alpha"]))
  cat(sprintf("  beta:        %.%df  (SE: %.%df)\n", digits, x$beta, digits, x$se["beta"]))
  cat(sprintf("  persistence: %.%df\n", digits, x$persistence))
  
  if (!is.null(x$shape)) {
    cat(sprintf("  shape:       %.%df  (SE: %.%df)\n", digits, x$shape, digits, x$se["shape"]))
  }
  cat("\n")
  
  cat(sprintf("%.0f%% Confidence Intervals:\n", x$ci_level * 100))
  print(round(x$ci, digits))
  cat("\n")
  
  cat("Information Matrix Diagnostics:\n")
  cat(sprintf("  Condition number: %.2e\n", x$info_condition))
  cat(sprintf("  Min eigenvalue: %.2e\n", min(x$info_eigenvalues)))
  cat("\n")
  
  invisible(x)
}