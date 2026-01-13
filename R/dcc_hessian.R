## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## DCC(1,1) Hessian and Standard Error Computation
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
##
## This file implements Hessian computation and standard errors for DCC(1,1).
##
## ARCHITECTURE
## - - - - - - 
## The functions are organized in layers:
##
##   Layer 1 (Low-level): Core numerical routines
##     - numerical_hessian()
##     - numerical_hessian_richardson()
##     - dcc11_analytical_hessian()
##     - .compute_dR_dQ(), .compute_d2R_dQ2()
##
##   Layer 2 (Mid-level): SE computation engines
##     - dcc11_hessian()          -- Hessian + vcov + SE + diagnostics
##     - dcc11_sandwich_se()      -- Sandwich/robust SE
##
##   Layer 3 (High-level): User-facing functions
##     - dcc11_standard_errors()  -- Simple interface, dispatches by method
##     - compute_dcc_standard_errors_robust() -- Wrapper with edge case handling
##     - dcc11_estimation_summary() -- Full summary with diagnostics
##
## SE METHODS SUPPORTED
## - - - - - - - - - - 
##   1. "hessian"  - Inverse observed information (fast, may underestimate beta SE)
##   2. "sandwich" - Robust/heteroskedasticity-consistent (slower, more robust)
##
## IMPORTANT: Inference Recommendations
## - - - - - - - - - - - - - - - - - - 
## Validation studies show that the DCC likelihood surface is highly anisotropic:
##   - Curvature in alpha direction is ~8x steeper than beta direction
##   - Hessian-based SEs for beta are severely underestimated (by ~5-6x)
##   - This is a fundamental property of DCC with high persistence, not a bug
##
## Recommended inference methods:
##   - Alpha: Hessian-based SE is acceptable (ratio to true SD ~0.8)
##   - Beta:  Use Bootstrap SE or Profile Likelihood CI (see dcc_inference.R)
##   - Persistence (alpha+beta): May be better identified than individual parameters
##
## See dcc_inference_guide.Rmd vignette for detailed findings and recommendations.
## See dcc_inference.R for bootstrap and profile likelihood implementations.
##
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


## SECTION 1: Low-Level Numerical Hessian Computation ==========================

#' @title Compute Numerical Hessian via Finite Differences
#'
#' @param fn Objective function (returns scalar)
#' @param params Parameter vector at which to evaluate Hessian
#' @param eps Step size for finite differences (default 1e-5)
#' @param ... Additional arguments passed to fn
#' 
#' @return Square matrix of second derivatives
#' 
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
        # Off-diagonal: use mixed partial derivative formula
        ## f_xy ≈ [f(x+h,y+h) - f(x+h,y-h) - f(x-h,y+h) + f(x-h,y-h)] / (4h^2)
        params_pp <- params_pm <- params_mp <- params_mm <- params
        params_pp[i] <- params[i] + eps; params_pp[j] <- params[j] + eps
        params_pm[i] <- params[i] + eps; params_pm[j] <- params[j] - eps
        params_mp[i] <- params[i] - eps; params_mp[j] <- params[j] + eps
        params_mm[i] <- params[i] - eps; params_mm[j] <- params[j] - eps
        f_pp <- fn(params_pp, ...); f_pm <- fn(params_pm, ...)
        f_mp <- fn(params_mp, ...); f_mm <- fn(params_mm, ...)
        H[i, j] <- H[j, i] <- (f_pp - f_pm - f_mp + f_mm) / (4 * eps^2)
      }
    }
  }
  
  ## Ensure perfect symmetry
  H <- (H + t(H)) / 2
  
  ## Add names if params are named
  if (!is.null(names(params))) rownames(H) <- colnames(H) <- names(params)
  return(H)
}


#' @title Compute Numerical Hessian with Richardson Extrapolation
#' 
#' @description More accurate Hessian using Richardson extrapolation.
#'   Combines estimates at different step sizes for higher-order accuracy.
#'   
#' @param fn Objective function
#' @param params Parameter vector
#' @param eps Base step size (default 1e-4)
#' @param r Reduction factor (default 2)
#' @param ... Additional arguments to fn
#' 
#' @return Hessian matrix with improved accuracy
#' 
#' @keywords internal
numerical_hessian_richardson <- function(fn, params, eps = 1e-4, r = 2, ...) {
  ## Compute Hessian at two step sizes
  H1 <- numerical_hessian(fn, params, eps = eps, ...)
  H2 <- numerical_hessian(fn, params, eps = eps / r, ...)
  
  ## Richardson extrapolation: H = (r^2 * H2 - H1) / (r^2 - 1)
  ## This eliminates the leading error term
  (r^2 * H2 - H1) / (r^2 - 1)
}


## SECTION 1b: Observed Information Matrix (convenience wrapper) ===============

#' Compute Observed Information Matrix for DCC(1,1)
#'
#' @description Compute the observed Fisher information matrix, which is the
#'   Hessian of the negative log-likelihood evaluated at the MLE.
#'   This is a convenience wrapper around numerical_hessian_richardson().
#'
#' @param params MLE parameter estimates (alpha, beta) or (psi, phi)
#' @param std_resid T x k matrix of standardized residuals
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional covariance matrix
#' @param distribution "mvn" or "mvt"
#' @param use_reparam Logical: parameters in (psi, phi) space?
#' @param eps Step size for numerical differentiation (default 1e-5)
#'
#' @return Observed information matrix (positive semi-definite if at MLE)
#'
#' @export
dcc11_observed_information <- function(
    params, 
    std_resid, 
    weights, 
    Qbar,
    distribution = "mvn", 
    use_reparam = FALSE,
    eps = 1e-5
  ) {
  ## Compute Hessian of NLL (using Richardson for accuracy)
  numerical_hessian_richardson(
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


## SECTION 2: Analytical Hessian (DCC-specific) ================================
##
## This section implements the analytical Hessian of the DCC(1,1) log-likelihood.
## The DCC model has a two-stage structure:
##   1. Q_t dynamics (pseudo-correlation process)
##   2. R_t = D_t^{-1} Q_t D_t^{-1} normalization to correlation matrix
##
## We need derivatives of both stages, composed via chain rule.
##

#' @title Analytical Hessian of DCC(1,1) Negative Log-Likelihood
#' 
#' @description Computes the exact Hessian matrix of the DCC(1,1) negative
#'   log-likelihood using analytical derivatives. This is more accurate than
#'   numerical differentiation but only implemented for the 2-parameter (alpha,
#'   beta) case.
#'   
#' @param params Parameter vector c(alpha, beta)
#' @param std_resid T x k matrix of standardized residuals from first-stage
#'   GARCH
#' @param weights T-vector of observation weights (typically all 1s)
#' @param Qbar k x k unconditional correlation matrix (sample correlation of
#'   std_resid)
#' @param distribution "mvn" (multivariate normal) - mvt not yet supported
#'   analytically
#'   
#' @return 2 x 2 Hessian matrix of the NLL with respect to (alpha, beta)
#'
#' @details The DCC(1,1) model specifies: \deqn{Q_t = (1 - \alpha - \beta)
#' \bar{Q} + \alpha z_{t-1} z_{t-1}' + \beta Q_{t-1}} \deqn{R_t =
#' diag(Q_t)^{-1/2} Q_t diag(Q_t)^{-1/2}}
#'
#' The correlation log-likelihood contribution at time t is: \deqn{\ell_t =
#' -\frac{1}{2}[\log|R_t| + z_t' R_t^{-1} z_t - z_t' z_t]}
#'
#' This function computes \eqn{\partial^2 (-\sum_t \ell_t) / \partial \theta
#' \partial \theta'} where \eqn{\theta = (\alpha, \beta)'}.
#'
#' @keywords internal
dcc11_analytical_hessian <- function(
    params, 
    std_resid, 
    weights, 
    Qbar,
    distribution = "mvn"
  ) {
  alpha <- params[1]
  beta <- params[2]
  n <- nrow(std_resid)
  k <- ncol(std_resid)
  
  
  ## Initialize state variables for the forward recursion ======================
  
  ## Q_t and its derivatives are propagated forward through time.
  ## At t=1, we use Q_1 = Qbar with zero derivatives (no dependence on params
  ## for the initial condition).
  Q <- Qbar
  dQ_dalpha <- matrix(0, k, k)
  dQ_dbeta <- matrix(0, k, k)
  d2Q_dalpha2 <- matrix(0, k, k)
  d2Q_dbeta2 <- matrix(0, k, k)
  d2Q_dalphadbeta <- matrix(0, k, k)
  
  ## Accumulators for Hessian elements (of the LOG-LIKELIHOOD, not NLL)
  H11 <- 0  # d²ℓ/dα²
  H22 <- 0  # d²ℓ/dβ²
  H12 <- 0  # d²ℓ/dαdβ
  
  
  ## Main forward recursion ====================================================

  for (t in 1:n) {
    z_t <- std_resid[t, , drop = FALSE]
    w_t <- weights[t]
    
    ## *** Step 1: Compute R_t from Q_t via normalization ***
    ## Normalize Q to get R
    ## R_t = D_t^{-1} Q_t D_t^{-1}, where D_t = diag(sqrt(diag(Q_t)))
    
    q_diag <- pmax(diag(Q), 1e-10)  ## Diagonal elements of Q_t
    d_inv <- 1 / sqrt(q_diag)       ## D_t^{-1} diagonal elements
    D_inv <- diag(d_inv, k)
    R <- D_inv %*% Q %*% D_inv
    
    ## Ensure R is a valid correlation matrix
    diag(R) <- 1
    R <- (R + t(R)) / 2
    
    ## Invert R (with fallback to pseudo-inverse for near-singular cases)
    R_inv <- tryCatch(
      solve(R),
      error = function(e) {
        eig <- eigen(R, symmetric = TRUE)
        eig$vectors %*% diag(1 / pmax(eig$values, 1e-8)) %*% t(eig$vectors)
      }
    )
    
    
    ## *** Step 2: Compute derivatives of R_t w.r.t. parameters ***
    
    ## These use the chain rule: dR/dθ = (dR/dQ)(dQ/dθ)
    ## The helper functions handle the Q -> R transformation derivatives.
    
    ## First derivatives of R w.r.t. alpha and beta (through Q)
    dR_dalpha <- .compute_dR_dQ(Q, dQ_dalpha, k)
    dR_dbeta <- .compute_dR_dQ(Q, dQ_dbeta, k)
    
    ## Second derivatives of R
    d2R_dalpha2 <- .compute_d2R_dQ2(Q, dQ_dalpha, dQ_dalpha, d2Q_dalpha2, k)
    d2R_dbeta2 <- .compute_d2R_dQ2(Q, dQ_dbeta, dQ_dbeta, d2Q_dbeta2, k)
    d2R_dalphadbeta <- .compute_d2R_dQ2(Q, dQ_dalpha, dQ_dbeta, d2Q_dalphadbeta, k)
    
    
    ## *** Step 3: Accumulate Hessian contributions from time t ***
    
    ## The log-likelihood contribution at t is:
    ##   ℓ_t = -0.5 * [log|R_t| + z_t' R_t^{-1} z_t - z_t' z_t]
    ##
    ## The Hessian of ℓ_t w.r.t. θ_i, θ_j involves:
    ##   d²ℓ_t/dθ_i dθ_j = 0.5 * tr(R^{-1} dR/dθ_i R^{-1} dR/dθ_j)
    ##                   - 0.5 * tr(R^{-1} d²R/dθ_i dθ_j)
    ##                   - z' (dR/dθ_i) R^{-1} (dR/dθ_j) R^{-1} z
    ##                   + 0.5 * z' (d²R/dθ_i dθ_j) R^{-1} z
    ##                   + [additional terms from R^{-1} derivatives]
    ##
    ## We precompute common quantities for efficiency.
    
    z_vec <- as.vector(z_t)
    Rz <- R_inv %*% z_vec  ## R^{-1} z
    
    ## --- Alpha-alpha (H11) ---
    ## d²NLL/dα²
    RdRa <- R_inv %*% dR_dalpha
    RdRaRz <- RdRa %*% Rz
    
    H11 <- H11 + w_t * (
      0.5 * sum(diag(RdRa %*% RdRa)) -
        0.5 * sum(diag(R_inv %*% d2R_dalpha2)) -
        as.numeric(t(Rz) %*% dR_dalpha %*% RdRaRz) +
        0.5 * as.numeric(t(Rz) %*% d2R_dalpha2 %*% Rz)
    )
    
    ## --- Beta-beta (H22) ---
    ## d²NLL/dβ²
    RdRb <- R_inv %*% dR_dbeta
    RdRbRz <- RdRb %*% Rz
    
    H22 <- H22 + w_t * (
      0.5 * sum(diag(RdRb %*% RdRb)) -
        0.5 * sum(diag(R_inv %*% d2R_dbeta2)) -
        as.numeric(t(Rz) %*% dR_dbeta %*% RdRbRz) +
        0.5 * as.numeric(t(Rz) %*% d2R_dbeta2 %*% Rz)
    )
    
    ## --- Alpha-beta (H12) ---
    ## d²NLL/dαdβ
    term1_ab <- 0.5 * sum(diag(RdRa %*% RdRb))
    term2_ab <- -0.5 * sum(diag(R_inv %*% d2R_dalphadbeta))
    term3_ab <- -as.numeric(t(Rz) %*% dR_dalpha %*% RdRbRz)
    term4_ab <- 0.5 * as.numeric(t(Rz) %*% d2R_dalphadbeta %*% Rz)
    
    H12 <- H12 + w_t * (term1_ab + term2_ab + term3_ab + term4_ab)
    
    
    ## *** Step 4: Update Q_t and its derivatives for next iteration ***
    
    ## The DCC recursion is:
    ##   Q_{t+1} = (1 - α - β) Q̄ + α z_t z_t' + β Q_t
    ##
    ## Differentiating w.r.t. α:
    ##   dQ_{t+1}/dα = -Q̄ + z_t z_t' + β (dQ_t/dα)
    ##
    ## Differentiating w.r.t. β:
    ##   dQ_{t+1}/dβ = -Q̄ + Q_t + β (dQ_t/dβ)
    ##
    ## Second derivatives (note the order of operations matters!):
    ##   d²Q_{t+1}/dα² = β (d²Q_t/dα²)
    ##
    ##   d²Q_{t+1}/dβ² = d/dβ[-Q̄ + Q_t + β(dQ_t/dβ)]
    ##                 = dQ_t/dβ + (dQ_t/dβ) + β(d²Q_t/dβ²)
    ##                 = 2(dQ_t/dβ) + β(d²Q_t/dβ²)
    ##
    ##   d²Q_{t+1}/dαdβ = d/dβ[-Q̄ + z_t z_t' + β(dQ_t/dα)]
    ##                  = dQ_t/dα + β(d²Q_t/dαdβ)
    ##
    ## IMPORTANT: We must update second derivatives BEFORE first derivatives,
    ## and first derivatives BEFORE Q, because each update uses the current
    ## (not yet updated) values.
    
    zz <- t(z_t) %*% z_t
    
    ## Second derivatives first (depend on current dQ values)
    d2Q_dalpha2 <- beta * d2Q_dalpha2
    d2Q_dbeta2 <- 2 * dQ_dbeta + beta * d2Q_dbeta2
    d2Q_dalphadbeta <- dQ_dalpha + beta * d2Q_dalphadbeta
    
    ## Then first derivatives (depend on current Q)
    dQ_dalpha <- -Qbar + zz + beta * dQ_dalpha
    dQ_dbeta <- -Qbar + Q + beta * dQ_dbeta
    
    ## Finally update Q itself
    Q <- (1 - alpha - beta) * Qbar + alpha * zz + beta * Q
  }
  
  
  ## Construct and return the Hessian matrix
  
  ## The loop above accumulated the Hessian of the LOG-LIKELIHOOD (ℓ).
  ## Since we want the Hessian of the NEGATIVE log-likelihood (NLL = -ℓ),
  ## we negate the result. This ensures the Hessian is positive definite
  ## at the MLE (since we're minimizing NLL).
  
  H <- -matrix(c(H11, H12, H12, H22), nrow = 2, ncol = 2)
  rownames(H) <- colnames(H) <- c("alpha", "beta")
  
  return(H)
}



## Helper: Derivative of R w.r.t. Q (first-order)
##
## R = D^{-1} Q D^{-1} where D = diag(sqrt(diag(Q)))
##
## By the product rule:
##   dR = d(D^{-1}) Q D^{-1} + D^{-1} dQ D^{-1} + D^{-1} Q d(D^{-1})
##
## where d(D^{-1})/dQ_ii = -0.5 * Q_ii^{-3/2} * dQ_ii
##
#' @keywords internal
.compute_dR_dQ <- function(Q, dQ, k) {
  q_diag <- pmax(diag(Q), 1e-10)
  d_inv <- 1 / sqrt(q_diag)
  
  ## Derivative of D^{-1} diagonal: d(q^{-1/2})/dq = -0.5 * q^{-3/2}
  dd_inv <- -0.5 * d_inv^3 * diag(dQ)
  
  D_inv <- diag(d_inv, k)
  dD_inv <- diag(dd_inv, k)
  
  ## Product rule expansion
  dR <- dD_inv %*% Q %*% D_inv +
    D_inv %*% dQ %*% D_inv +
    D_inv %*% Q %*% dD_inv
  
  ## R is a correlation matrix: diagonal must be zero derivative
  ## (since R_ii = 1 always), and ensure symmetry
  diag(dR) <- 0
  dR <- (dR + t(dR)) / 2
  
  return(dR)
}



## Helper: Second derivative of R w.r.t. Q (for mixed partials)
##
## This computes d²R/(dθ_i dθ_j) given dQ/dθ_i, dQ/dθ_j, and d²Q/(dθ_i dθ_j).
##
## The calculation involves second derivatives of D^{-1}:
##   d²(q^{-1/2})/(dq)² = 0.75 * q^{-5/2}
##
## and the full product rule expansion for R = D^{-1} Q D^{-1}.
##
#' @keywords internal
.compute_d2R_dQ2 <- function(Q, dQ1, dQ2, d2Q, k) {
  q_diag <- pmax(diag(Q), 1e-10)
  d_inv <- 1 / sqrt(q_diag)
  
  ## Extract diagonal elements of the dQ matrices
  dq1 <- diag(dQ1)
  dq2 <- diag(dQ2)
  d2q <- diag(d2Q)
  
  ## First derivatives of D^{-1} diagonal
  dd1 <- -0.5 * d_inv^3 * dq1
  dd2 <- -0.5 * d_inv^3 * dq2
  
  ## Second derivative of D^{-1} diagonal:
  ## d²(q^{-1/2})/dθ_i dθ_j = 0.75 * q^{-5/2} * (dq/dθ_i)(dq/dθ_j)
  ##                        - 0.5 * q^{-3/2} * (d²q/dθ_i dθ_j)
  d2d <- 0.75 * d_inv^5 * dq1 * dq2 - 0.5 * d_inv^3 * d2q
  
  ## Build diagonal matrices
  D <- diag(d_inv, k)
  dD1 <- diag(dd1, k)
  dD2 <- diag(dd2, k)
  d2D <- diag(d2d, k)
  
  ## Full product rule expansion for d²(D^{-1} Q D^{-1})
  ## There are 9 terms from differentiating each of the 3 factors twice
  d2R <- d2D %*% Q %*% D +      # d²D * Q * D
    dD1 %*% dQ2 %*% D +    # dD_1 * dQ_2 * D
    dD1 %*% Q %*% dD2 +    # dD_1 * Q * dD_2
    dD2 %*% dQ1 %*% D +    # dD_2 * dQ_1 * D
    D %*% d2Q %*% D +      # D * d²Q * D
    D %*% dQ1 %*% dD2 +    # D * dQ_1 * dD_2
    dD2 %*% Q %*% dD1 +    # dD_2 * Q * dD_1
    D %*% dQ2 %*% dD1 +    # D * dQ_2 * dD_1
    D %*% Q %*% d2D        # D * Q * d²D
  
  ## Enforce correlation matrix structure
  diag(d2R) <- 0
  d2R <- (d2R + t(d2R)) / 2
  
  return(d2R)
}


## SECTION 3: Mid-Level SE Computation Engines =================================

#' Compute DCC(1,1) Hessian with Full Diagnostics
#'
#' @description Core function for Hessian-based inference. Returns Hessian,
#'   variance-covariance matrix, standard errors, and diagnostic information
#'   (eigenvalues, condition number) useful for detecting the "flat beta problem".
#'
#' @param params MLE parameter estimates c(alpha, beta) or c(alpha, beta, shape)
#' @param std_resid T x k matrix of standardized residuals
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional correlation matrix
#' @param distribution "mvn" or "mvt"
#' @param use_reparam Logical: parameters in reparameterized (psi, phi) space?
#' @param hessian_method "numerical" (default) or "analytical"
#' @param eps Step size for numerical differentiation
#'
#' @return List with: 
#'   \item{hessian}{Hessian matrix of NLL}
#'   \item{info}{Observed information matrix (= Hessian for NLL)}
#'   \item{vcov}{Variance-covariance matrix (inverse of info)}
#'   \item{se}{Standard errors}
#'   \item{eigenvalues}{Eigenvalues of Hessian}
#'   \item{eigenvectors}{Eigenvectors of Hessian}
#'   \item{condition_number}{Condition number of Hessian}
#'   \item{params}{Parameters}
#'   \item{param_names}{Para,eter names}
#'   \item{method}{"hessian"}
#'
#' @export
dcc11_hessian <- function(
    params, 
    std_resid, 
    weights, 
    Qbar,
    distribution = "mvn", 
    use_reparam = FALSE,
    hessian_method = "numerical", 
    eps = 1e-5
  ) {
  n_params <- length(params)
  param_names <- if (n_params == 2) {
    if (use_reparam) c("psi", "phi") else c("alpha", "beta")
  } else if (n_params == 3) {
    if (use_reparam) c("psi", "phi", "shape") else c("alpha", "beta", "shape")
  } else paste0("param", 1:n_params)
  
  H <- if (hessian_method == "analytical" && n_params == 2) {
    dcc11_analytical_hessian(params, std_resid, weights, Qbar, distribution)
  } else {
    numerical_hessian_richardson(
      dcc11_nll, 
      params, 
      eps, 
      std_resid = std_resid,
      weights = weights, 
      Qbar = Qbar,
      distribution = distribution,
      use_reparam = use_reparam
    )
  }
  rownames(H) <- colnames(H) <- param_names
  
  eig <- eigen(H, symmetric = TRUE)
  vcov <- tryCatch({
    if (all(eig$values > 1e-10)) solve(H)
    else {
      warning("Hessian is near-singular, using pseudo-inverse")
      eig$vectors %*% diag(1 / pmax(eig$values, 1e-10)) %*% t(eig$vectors)
    }
  }, error = function(e) matrix(NA_real_, n_params, n_params))
  rownames(vcov) <- colnames(vcov) <- param_names
  
  se <- sqrt(pmax(0, diag(vcov)))
  names(se) <- param_names
  
  list(
    hessian = H, 
    info = H, 
    vcov = vcov, 
    se = se,
    eigenvalues = eig$values, 
    eigenvectors = eig$vectors,
    condition_number = max(eig$values) / max(min(eig$values), 1e-16),
    params = params, 
    param_names = param_names,
    method = "hessian"
  )
}


#' @title Compute Sandwich (Robust) Standard Errors
#'
#' @description Heteroskedasticity-robust SEs using sandwich estimator:
#'   Var(theta) = I^{-1} J I^{-1}
#'
#' @param params MLE parameter estimates
#' @param std_resid T x k matrix of standardized residuals
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional covariance matrix
#' @param distribution "mvn" or "mvt"
#' @param use_reparam Logical: parameters in (psi, phi) space?
#' 
#' @details
#' The sandwich (robust) variance estimator is:
#' \deqn{Var(\hat{\theta}) = H^{-1} J H^{-1}}
#' where H is the Hessian and J is the outer product of score contributions.
#' This is consistent under heteroskedasticity when the standard Hessian-based
#' estimator may not be.
#'
#' @return List with: 
#'   \item{se}{Standard errors (square root of diagonal of vcov)}
#'   \item{vcov}{Sandwich variance-covariance matrix: (\eqn{H^{-1} J H^{-1}}}
#'   \item{bread}{\eqn{H^{−1}}, inverse Hessian of the negative log-likelihood}
#'   \item{meat}{J, outer product of score vectors: \eqn{sum_t (grad_t)(grad_t)'}}
#'   \item{params}{Parameter values at which SEs were computed}
#'   \item{param_names}{Names of parameters}
#'   \item{method}{"sandwich"}
#'
#' @export
dcc11_sandwich_se <- function(
    params, 
    std_resid, 
    weights, 
    Qbar,
    distribution = "mvn", 
    use_reparam = FALSE
  ) {
  T_obs <- nrow(std_resid); n_params <- length(params)
  
  hess_result <- dcc11_hessian(
    params, 
    std_resid, 
    weights, 
    Qbar, 
    distribution, 
    use_reparam
  )
  
  if (is.null(hess_result$vcov) || any(is.na(hess_result$vcov))) {
    warning("Cannot compute sandwich SE: information matrix is singular")
    hess_result$method <- "sandwich"
    return(hess_result)
  }
  
  J <- matrix(0, n_params, n_params)
  for (t in 1:T_obs) {
    grad_t <- .numerical_gradient_single_obs(
      params, 
      std_resid, 
      weights, 
      Qbar, 
      t,
      distribution, 
      use_reparam, 
      eps = 1e-6
    )
    J <- J + weights[t]^2 * outer(grad_t, grad_t)
  }
  
  bread <- hess_result$vcov
  sandwich <- bread %*% J %*% bread
  var_diag <- diag(sandwich); var_diag[var_diag < 0] <- NA_real_
  se_sandwich <- sqrt(var_diag)
  names(se_sandwich) <- hess_result$param_names
  rownames(sandwich) <- colnames(sandwich) <- hess_result$param_names
  
  list(
    se = se_sandwich, 
    vcov = sandwich, 
    bread = bread, 
    meat = J,
    info = hess_result$info, 
    params = params,
    param_names = hess_result$param_names,
    method = "sandwich"
  )
}


#' @keywords internal
.numerical_gradient_single_obs <- function(
    params, 
    std_resid, 
    weights, 
    Qbar,
    t, 
    distribution, 
    use_reparam, 
    eps = 1e-6
  ) {
  n_params <- length(params); grad <- numeric(n_params)
  for (i in 1:n_params) {
    p_plus <- p_minus <- params
    p_plus[i] <- params[i] + eps; p_minus[i] <- params[i] - eps
    grad[i] <- (.dcc11_nll_single_obs(
      p_plus, 
      std_resid, 
      weights, 
      Qbar, 
      t, 
      distribution, 
      use_reparam
    ) - .dcc11_nll_single_obs(
      p_minus, 
      std_resid, 
      weights, 
      Qbar, 
      t, 
      distribution, 
      use_reparam
    )
  ) / (2 * eps)
  }
  grad
}


#' @keywords internal
.dcc11_nll_single_obs <- function(
    params, 
    std_resid, 
    weights, 
    Qbar, 
    t,
    distribution, 
    use_reparam
  ) {
  k <- ncol(std_resid)
  if (use_reparam) {
    orig <- dcc11_from_unconstrained(params[1], params[2])
    alpha <- as.numeric(orig["alpha"]); beta <- as.numeric(orig["beta"])
  } else {
    alpha <- as.numeric(params[1]); beta <- as.numeric(params[2])
  }
  shape <- if (distribution == "mvt") as.numeric(params[length(params)]) else NULL
  
  if (!is.finite(alpha) || !is.finite(beta) || alpha + beta >= 1 || alpha < 0 || beta < 0)
    return(1e10)
  
  fwd <- dcc11_recursion_with_grad(std_resid, alpha, beta, Qbar)
  if (!isTRUE(fwd$success)) return(1e10)
  
  z_t <- std_resid[t, ]; R_inv_t <- fwd$R_inv[, , t]; log_det_R_t <- fwd$log_det_R[t]
  q_t <- as.numeric(t(z_t) %*% R_inv_t %*% z_t)
  
  ll_t <- if (distribution == "mvn") {
    -0.5 * (log_det_R_t + q_t - sum(z_t^2))
  } else {
    lgamma((shape + k) / 2) - lgamma(shape / 2) - (k / 2) * log(pi * (shape - 2)) -
      0.5 * log_det_R_t - ((shape + k) / 2) * log(1 + q_t / (shape - 2))
  }
  -weights[t] * ll_t
}


## SECTION 4: High-Level User-Facing Functions =================================

#' @title Compute Standard Errors for DCC(1,1) Parameters
#'
#' @description Main user-facing function for SE computation. Dispatches to
#'   appropriate method based on the \code{method} argument.
#'
#' @param params MLE parameter estimates
#' @param std_resid T x k matrix of standardized residuals
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional covariance matrix
#' @param distribution "mvn" or "mvt"
#' @param use_reparam Logical: parameters in (psi, phi) space?
#' @param method SE method: "hessian" (default) or "sandwich"
#'
#' @return List with: 
#'   If \code{method="hessian"}:
#'   \item{hessian}{Hessian matrix of NLL}
#'   \item{info}{Observed information matrix (= Hessian for NLL)}
#'   \item{vcov}{Variance-covariance matrix (inverse of info)}
#'   \item{se}{Standard errors}
#'   \item{eigenvalues}{Eigenvalues of Hessian}
#'   \item{eigenvectors}{Eigenvectors of Hessian}
#'   \item{condition_number}{Condition number of Hessian}
#'   \item{params}{Parameters}
#'   \item{param_names}{Para,eter names}
#'   \item{method}{"hessian"}
#'
#'   If \code{method="sandwich"}:
#'   \item{se}{Standard errors}
#'   \item{vcov}{Variance-covariance matrix (inverse of info)}
#'   \item{bread}{Bread}
#'   \item{meatinfo}{Meat}
#'   \item{vcov}{Variance-covariance matrix (inverse of info)}
#'   \item{params}{Parameters}
#'   \item{param_names}{Parameter names}
#'   \item{method}{"sandwich"}

#'
#' @details
#' For reliable beta inference in high-persistence DCC models, consider using
#' bootstrap SEs or profile likelihood CIs from \code{dcc_inference.R}.
#'
#' @export
dcc11_standard_errors <- function(
    params, 
    std_resid, 
    weights, 
    Qbar,
    distribution = "mvn", 
    use_reparam = FALSE,
    method = c("hessian", "sandwich")
  ) {
  method <- match.arg(method)
  if (method == "sandwich") {
    dcc11_sandwich_se(
      params, 
      std_resid, 
      weights, 
      Qbar, 
      distribution, 
      use_reparam
    )
  } else {
    dcc11_hessian(
      params, 
      std_resid, 
      weights, 
      Qbar, 
      distribution, 
      use_reparam
    )
  }
}


#' Compute DCC Standard Errors with Robust Edge Case Handling
#'
#' @description High-level wrapper with handling for boundary estimates,
#'   constant correlation, and non-positive-definite Hessians.
#'
#' @param dcc_result Result from estimate_dcc_parameters_weighted() or similar
#' @param std_resid T x k matrix of standardized residuals
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional correlation matrix
#' @param distribution "mvn" or "mvt"
#' @param boundary_threshold Threshold for boundary detection (default 1e-4)
#' @param method "hessian" or "sandwich"
#'
#' @return List with: 
#'   \item{se}{se}
#'   \item{vcov}{vcov}
#'   \item{valid}{valid}
#'   \item{reason}{reason}
#'   \item{correlation_type}{correlation_type}
#'   \item{info}{info}
#'   \item{info_eigenvalues}{info_eigenvalues}
#'
#' @export
compute_dcc_standard_errors_robust <- function(
    dcc_result, 
    std_resid, 
    weights, 
    Qbar,
    distribution = "mvn",
    boundary_threshold = 1e-4,
    method = c("hessian", "sandwich")
  ) {
  method <- match.arg(method)
  result <- list(se = NULL, vcov = NULL, valid = FALSE, reason = NULL,
                 correlation_type = NULL, info = NULL, info_eigenvalues = NULL)
  
  correlation_type <- dcc_result$correlation_type %||%
    dcc_result$coefficients$correlation_type %||% "dynamic"
  result$correlation_type <- correlation_type
  
  if (correlation_type == "constant") {
    result$reason <- "constant_correlation"; result$valid <- TRUE
    return(result)
  }
  
  dcc_pars <- dcc_result$dcc_pars %||% dcc_result$coefficients$dcc_pars
  if (is.null(dcc_pars) || length(dcc_pars) == 0) {
    result$reason <- "no_dcc_parameters"; return(result)
  }
  
  alpha <- dcc_pars$alpha_1 %||% dcc_pars[["alpha_1"]]
  beta <- dcc_pars$beta_1 %||% dcc_pars[["beta_1"]]
  if (is.null(alpha) || is.null(beta)) {
    result$reason <- "missing_alpha_or_beta"; return(result)
  }
  
  if (alpha < boundary_threshold || beta < boundary_threshold) {
    result$reason <- sprintf("boundary_lower (alpha=%.2e, beta=%.2e)", alpha, beta)
    result$se <- c(alpha = NA_real_, beta = NA_real_)
    if (distribution == "mvt") result$se <- c(result$se, shape = NA_real_)
    return(result)
  }
  if ((alpha + beta) > (1 - boundary_threshold)) {
    result$reason <- sprintf("boundary_upper (persistence=%.6f)", alpha + beta)
    result$se <- c(alpha = NA_real_, beta = NA_real_)
    if (distribution == "mvt") result$se <- c(result$se, shape = NA_real_)
    return(result)
  }
  
  params <- c(alpha = alpha, beta = beta)
  if (distribution == "mvt") {
    shape <- dcc_result$dist_pars$shape %||% dcc_result$coefficients$dist_pars$shape %||% 8.0
    params <- c(params, shape = shape)
  }
  
  se_result <- tryCatch(
    dcc11_standard_errors(params, std_resid, weights, Qbar, distribution, FALSE, method),
    error = function(e) list(se = NULL, vcov = NULL, info = NULL, error = e$message)
  )
  
  if (!is.null(se_result$error)) {
    result$reason <- paste("computation_error:", se_result$error)
    result$se <- c(alpha = NA_real_, beta = NA_real_)
    if (distribution == "mvt") result$se <- c(result$se, shape = NA_real_)
    return(result)
  }
  
  result$se <- se_result$se
  result$vcov <- se_result$vcov
  result$info <- se_result$info
  
  if (!is.null(se_result$info)) {
    eig <- eigen(se_result$info, symmetric = TRUE, only.values = TRUE)$values
    result$info_eigenvalues <- eig
    if (any(eig <= 0)) {
      result$reason <- sprintf("non_pd_information (min_eig=%.2e)", min(eig))
      return(result)
    }
  }
  
  if (any(is.na(se_result$se))) { result$reason <- "na_standard_errors"; return(result) }
  if (any(se_result$se <= 0)) { result$reason <- "non_positive_standard_errors"; return(result) }
  
  result$valid <- TRUE; result$reason <- "ok"
  result
}


## SECTION 5: SE Transformation and Confidence Intervals =======================

#' @title Transform Standard Errors Between Parameterizations
#'
#' @param se_result Result from dcc11_standard_errors() or dcc11_hessian()
#' @param to_reparam Logical: transform TO (psi, phi)?
#'
#' @return Updated se_result with transformed SEs
#'
#' @export
dcc11_transform_se <- function(se_result, to_reparam = FALSE) {
  params <- se_result$params; vcov_orig <- se_result$vcov
  
  if (to_reparam) {
    alpha <- params[1]; beta <- params[2]; pers <- alpha + beta
    d_psi <- 1 / (pers * (1 - pers))
    J <- matrix(c(d_psi, 1/alpha, d_psi, -1/beta), 2, 2, byrow = FALSE)
    new_names <- c("psi", "phi")
  } else {
    J <- dcc11_reparam_jacobian(params[1], params[2])
    new_names <- c("alpha", "beta")
  }
  
  if (length(params) > 2) {
    J_full <- diag(length(params)); J_full[1:2, 1:2] <- J; J <- J_full
    new_names <- c(new_names, "shape")
  }
  
  vcov_new <- J %*% vcov_orig %*% t(J)
  se_result$vcov <- vcov_new
  rownames(se_result$vcov) <- colnames(se_result$vcov) <- new_names
  se_result$se <- sqrt(pmax(0, diag(vcov_new)))
  names(se_result$se) <- new_names
  se_result$param_names <- new_names
  se_result$jacobian <- J
  se_result
}


#' @title Compute Confidence Intervals for DCC(1,1) Parameters
#'
#' @param se_result Result from dcc11_standard_errors() or dcc11_hessian()
#' @param level Confidence level (default 0.95)
#'
#' @return Matrix with columns: estimate, se, lower, upper
#'
#' @export
dcc11_confint <- function(se_result, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)
  params <- se_result$params; se <- se_result$se
  ci <- cbind(estimate = params, se = se, lower = params - z * se, upper = params + z * se)
  rownames(ci) <- se_result$param_names
  attr(ci, "level") <- level; attr(ci, "method") <- se_result$method
  ci
}


#' @title Compute Profile Likelihood Confidence Intervals
#'
#' @param params MLE estimates
#' @param std_resid Standardized residuals
#' @param weights Observation weights
#' @param Qbar Unconditional covariance
#' @param distribution "mvn" or "mvt"
#' @param use_reparam Logical
#' @param param_idx Index of parameter for CI
#' @param level Confidence level (default 0.95)
#' @param n_grid Grid points (default 50)
#'
#' @return Named vector c(lower, upper)
#'
#' @export
dcc11_profile_ci <- function(
    params, 
    std_resid, 
    weights, 
    Qbar,
    distribution = "mvn", 
    use_reparam = FALSE,
    param_idx = 1, 
    level = 0.95, 
    n_grid = 50
  ) {
  chi2_crit <- qchisq(level, df = 1)
  nll_mle <- dcc11_nll(params, std_resid, weights, Qbar, distribution, use_reparam)
  threshold <- nll_mle + chi2_crit / 2
  
  param_mle <- params[param_idx]
  se_approx <- tryCatch(
    dcc11_standard_errors(
      params, 
      std_resid, 
      weights, 
      Qbar, 
      distribution,
      use_reparam
    )$se[param_idx],
    error = function(e) NA_real_
  )
  if (is.na(se_approx)) se_approx <- abs(param_mle) * 0.1 + 0.01
  
  search_range <- 4 * se_approx
  lower_grid <- seq(param_mle - search_range, param_mle, length.out = n_grid)
  upper_grid <- seq(param_mle, param_mle + search_range, length.out = n_grid)
  
  profile_nll <- function(p) {
    p_test <- params; p_test[param_idx] <- p
    dcc11_nll(p_test, std_resid, weights, Qbar, distribution, use_reparam)
  }
  
  nll_lower <- sapply(lower_grid, profile_nll)
  below <- which(nll_lower <= threshold)
  lower_ci <- if (length(below) > 0) lower_grid[min(below)] else lower_grid[1]
  
  nll_upper <- sapply(upper_grid, profile_nll)
  below <- which(nll_upper <= threshold)
  upper_ci <- if (length(below) > 0) upper_grid[max(below)] else upper_grid[n_grid]
  
  c(lower = lower_ci, upper = upper_ci)
}


## SECTION 6: Summary and Printing =============================================

#' @title Summarize DCC(1,1) Estimation Results
#'
#' @param params MLE estimates
#' @param std_resid Standardized residuals
#' @param weights Observation weights
#' @param Qbar Unconditional covariance
#' @param distribution "mvn" or "mvt"
#' @param use_reparam Logical
#' @param level Confidence level (default 0.95)
#' @param method SE method: "hessian" or "sandwich"
#'
#' @return List of class "dcc11_summary"
#'
#' @export
dcc11_estimation_summary <- function(
    params, 
    std_resid, 
    weights, 
    Qbar,
    distribution = "mvn", 
    use_reparam = FALSE,
    level = 0.95, 
    method = c("hessian", "sandwich")
  ) {
  method <- match.arg(method)
  T_obs <- nrow(std_resid); k <- ncol(std_resid)
  
  nll <- dcc11_nll(params, std_resid, weights, Qbar, distribution, use_reparam)
  
  se_result <- if (method == "sandwich") {
    dcc11_sandwich_se(params, std_resid, weights, Qbar, distribution, use_reparam)
  } else {
    dcc11_hessian(params, std_resid, weights, Qbar, distribution, use_reparam)
  }
  
  ci <- dcc11_confint(se_result, level)
  
  if (use_reparam) {
    orig <- dcc11_from_unconstrained(params[1], params[2])
    alpha <- orig["alpha"]; beta <- orig["beta"]
    se_orig <- dcc11_transform_se(se_result, to_reparam = FALSE)
    ci_orig <- dcc11_confint(se_orig, level)
  } else {
    alpha <- params[1]; beta <- params[2]
    se_orig <- se_result; ci_orig <- ci
  }
  
  n_params <- length(params)
  eig_values <- se_result$eigenvalues %||%
    (if (!is.null(se_result$info)) eigen(se_result$info, symmetric = TRUE, only.values = TRUE)$values else NA_real_)
  
  result <- list(
    n_obs = T_obs, n_series = k, n_params = n_params, distribution = distribution,
    log_likelihood = -nll, nll = nll,
    aic = 2 * nll + 2 * n_params, bic = 2 * nll + log(T_obs) * n_params,
    alpha = as.numeric(alpha), beta = as.numeric(beta),
    persistence = as.numeric(alpha + beta),
    shape = if (distribution == "mvt") params[length(params)] else NULL,
    method = method, se = se_orig$se, vcov = se_orig$vcov,
    ci = ci_orig, ci_level = level,
    info_eigenvalues = eig_values,
    info_condition = if (!is.null(se_result$info)) kappa(se_result$info) else NA
  )
  class(result) <- "dcc11_summary"
  result
}


#' @export
print.dcc11_summary <- function(x, digits = 4, ...) {
  cat("\n===== DCC(1,1) Estimation Summary =====\n\n")
  cat("Data:\n")
  cat(sprintf("  Observations: %d\n  Series: %d\n  Distribution: %s\n\n",
              x$n_obs, x$n_series, toupper(x$distribution)))
  cat("Likelihood:\n")
  cat(sprintf("  Log-likelihood: %.*f\n  AIC: %.*f\n  BIC: %.*f\n\n",
              digits, x$log_likelihood, digits, x$aic, digits, x$bic))
  cat(sprintf("Parameter Estimates (SE method: %s):\n", x$method))
  cat(sprintf("  alpha:       %.*f  (SE: %.*f)\n", digits, x$alpha, digits, x$se["alpha"]))
  cat(sprintf("  beta:        %.*f  (SE: %.*f)\n", digits, x$beta, digits, x$se["beta"]))
  cat(sprintf("  persistence: %.*f\n", digits, x$persistence))
  if (!is.null(x$shape))
    cat(sprintf("  shape:       %.*f  (SE: %.*f)\n", digits, x$shape, digits, x$se["shape"]))
  cat(sprintf("\n%.0f%% Confidence Intervals:\n", x$ci_level * 100))
  print(round(x$ci, digits))
  cat("\nInformation Matrix Diagnostics:\n")
  cat(sprintf("  Condition number: %.2e\n", x$info_condition))
  if (length(x$info_eigenvalues) > 1 && all(is.finite(x$info_eigenvalues))) {
    ratio <- max(x$info_eigenvalues) / max(min(x$info_eigenvalues), 1e-16)
    cat(sprintf("  Eigenvalue ratio: %.1f\n", ratio))
    if (ratio > 10) cat("  WARNING: High eigenvalue ratio - consider bootstrap SE for beta.\n")
  }
  cat("\n")
  invisible(x)
}