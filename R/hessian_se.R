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


#### ______________________________________________________________________ ####
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
        ## f''(x) â‰ˆ [f(x+h) - 2f(x) + f(x-h)] / h^2
        params_plus <- params_minus <- params
        params_plus[i] <- params[i] + eps
        params_minus[i] <- params[i] - eps
        f_plus <- fn(params_plus, ...)
        f_minus <- fn(params_minus, ...)
        H[i, i] <- (f_plus - 2 * f0 + f_minus) / (eps^2)
      } else {
        # Off-diagonal: use mixed partial derivative formula
        ## f_xy â‰ˆ [f(x+h,y+h) - f(x+h,y-h) - f(x-h,y+h) + f(x-h,y-h)] / (4h^2)
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


#### ______________________________________________________________________ ####
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
#' @param distribution "mvn" (multivariate normal). For "mvt", use 
#'   numerical_hessian_richardson() instead as analytical MVT Hessian
#'   is not yet implemented.
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
  H11 <- 0  # dÂ²â„“/dÎ±Â²
  H22 <- 0  # dÂ²â„“/dÎ²Â²
  H12 <- 0  # dÂ²â„“/dÎ±dÎ²
  
  
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
    
    ## These use the chain rule: dR/dÎ¸ = (dR/dQ)(dQ/dÎ¸)
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
    ##   â„“_t = -0.5 * [log|R_t| + z_t' R_t^{-1} z_t - z_t' z_t]
    ##
    ## The Hessian of â„“_t w.r.t. Î¸_i, Î¸_j involves:
    ##   dÂ²â„“_t/dÎ¸_i dÎ¸_j = 0.5 * tr(R^{-1} dR/dÎ¸_i R^{-1} dR/dÎ¸_j)
    ##                   - 0.5 * tr(R^{-1} dÂ²R/dÎ¸_i dÎ¸_j)
    ##                   - z' (dR/dÎ¸_i) R^{-1} (dR/dÎ¸_j) R^{-1} z
    ##                   + 0.5 * z' (dÂ²R/dÎ¸_i dÎ¸_j) R^{-1} z
    ##                   + [additional terms from R^{-1} derivatives]
    ##
    ## We precompute common quantities for efficiency.
    
    z_vec <- as.vector(z_t)
    Rz <- R_inv %*% z_vec  ## R^{-1} z
    
    ## --- Alpha-alpha (H11) ---
    ## dÂ²NLL/dÎ±Â²
    RdRa <- R_inv %*% dR_dalpha
    RdRaRz <- RdRa %*% Rz
    
    H11 <- H11 + w_t * (
      0.5 * sum(diag(RdRa %*% RdRa)) -
        0.5 * sum(diag(R_inv %*% d2R_dalpha2)) -
        as.numeric(t(Rz) %*% dR_dalpha %*% RdRaRz) +
        0.5 * as.numeric(t(Rz) %*% d2R_dalpha2 %*% Rz)
    )
    
    ## --- Beta-beta (H22) ---
    ## dÂ²NLL/dÎ²Â²
    RdRb <- R_inv %*% dR_dbeta
    RdRbRz <- RdRb %*% Rz
    
    H22 <- H22 + w_t * (
      0.5 * sum(diag(RdRb %*% RdRb)) -
        0.5 * sum(diag(R_inv %*% d2R_dbeta2)) -
        as.numeric(t(Rz) %*% dR_dbeta %*% RdRbRz) +
        0.5 * as.numeric(t(Rz) %*% d2R_dbeta2 %*% Rz)
    )
    
    ## --- Alpha-beta (H12) ---
    ## dÂ²NLL/dÎ±dÎ²
    term1_ab <- 0.5 * sum(diag(RdRa %*% RdRb))
    term2_ab <- -0.5 * sum(diag(R_inv %*% d2R_dalphadbeta))
    term3_ab <- -as.numeric(t(Rz) %*% dR_dalpha %*% RdRbRz)
    term4_ab <- 0.5 * as.numeric(t(Rz) %*% d2R_dalphadbeta %*% Rz)
    
    H12 <- H12 + w_t * (term1_ab + term2_ab + term3_ab + term4_ab)
    
    
    ## *** Step 4: Update Q_t and its derivatives for next iteration ***
    
    ## The DCC recursion is:
    ##   Q_{t+1} = (1 - Î± - Î²) QÌ„ + Î± z_t z_t' + Î² Q_t
    ##
    ## Differentiating w.r.t. Î±:
    ##   dQ_{t+1}/dÎ± = -QÌ„ + z_t z_t' + Î² (dQ_t/dÎ±)
    ##
    ## Differentiating w.r.t. Î²:
    ##   dQ_{t+1}/dÎ² = -QÌ„ + Q_t + Î² (dQ_t/dÎ²)
    ##
    ## Second derivatives (note the order of operations matters!):
    ##   dÂ²Q_{t+1}/dÎ±Â² = Î² (dÂ²Q_t/dÎ±Â²)
    ##
    ##   dÂ²Q_{t+1}/dÎ²Â² = d/dÎ²[-QÌ„ + Q_t + Î²(dQ_t/dÎ²)]
    ##                 = dQ_t/dÎ² + (dQ_t/dÎ²) + Î²(dÂ²Q_t/dÎ²Â²)
    ##                 = 2(dQ_t/dÎ²) + Î²(dÂ²Q_t/dÎ²Â²)
    ##
    ##   dÂ²Q_{t+1}/dÎ±dÎ² = d/dÎ²[-QÌ„ + z_t z_t' + Î²(dQ_t/dÎ±)]
    ##                  = dQ_t/dÎ± + Î²(dÂ²Q_t/dÎ±dÎ²)
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
  
  ## The loop above accumulated the Hessian of the LOG-LIKELIHOOD (â„“).
  ## Since we want the Hessian of the NEGATIVE log-likelihood (NLL = -â„“),
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
## This computes dÂ²R/(dÎ¸_i dÎ¸_j) given dQ/dÎ¸_i, dQ/dÎ¸_j, and dÂ²Q/(dÎ¸_i dÎ¸_j).
##
## The calculation involves second derivatives of D^{-1}:
##   dÂ²(q^{-1/2})/(dq)Â² = 0.75 * q^{-5/2}
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
  ## dÂ²(q^{-1/2})/dÎ¸_i dÎ¸_j = 0.75 * q^{-5/2} * (dq/dÎ¸_i)(dq/dÎ¸_j)
  ##                        - 0.5 * q^{-3/2} * (dÂ²q/dÎ¸_i dÎ¸_j)
  d2d <- 0.75 * d_inv^5 * dq1 * dq2 - 0.5 * d_inv^3 * d2q
  
  ## Build diagonal matrices
  D <- diag(d_inv, k)
  dD1 <- diag(dd1, k)
  dD2 <- diag(dd2, k)
  d2D <- diag(d2d, k)
  
  ## Full product rule expansion for dÂ²(D^{-1} Q D^{-1})
  ## There are 9 terms from differentiating each of the 3 factors twice
  d2R <- d2D %*% Q %*% D +      # dÂ²D * Q * D
    dD1 %*% dQ2 %*% D +    # dD_1 * dQ_2 * D
    dD1 %*% Q %*% dD2 +    # dD_1 * Q * dD_2
    dD2 %*% dQ1 %*% D +    # dD_2 * dQ_1 * D
    D %*% d2Q %*% D +      # D * dÂ²Q * D
    D %*% dQ1 %*% dD2 +    # D * dQ_1 * dD_2
    dD2 %*% Q %*% dD1 +    # dD_2 * Q * dD_1
    D %*% dQ2 %*% dD1 +    # D * dQ_2 * dD_1
    D %*% Q %*% d2D        # D * Q * dÂ²D
  
  ## Enforce correlation matrix structure
  diag(d2R) <- 0
  d2R <- (d2R + t(d2R)) / 2
  
  return(d2R)
}


#### ______________________________________________________________________ ####
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
#'   \item{bread}{\eqn{H^{âˆ’1}}, inverse Hessian of the negative log-likelihood}
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


#### ______________________________________________________________________ ####
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


#### ______________________________________________________________________ ####
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


#### ______________________________________________________________________ ####
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


#### ______________________________________________________________________ ####
## SECTION 7: CGARCH Standard Error Computation ================================
##
## CGARCH uses copula-based dependence, where the log-likelihood has a different
## structure than DCC. The copula likelihood subtracts marginal t densities.
##

#' @title Compute Standard Errors for CGARCH Models
#'
#' @description Computes standard errors for Copula GARCH DCC parameters using
#'   numerical differentiation of the copula log-likelihood. Supports both
#'   Gaussian (MVN) and Student-t (MVT) copulas.
#'
#' @param params Parameter vector:
#'   - For MVN copula: c(alpha, beta)
#'   - For MVT copula: c(alpha, beta, shape)
#' @param z_matrix T x k matrix of copula residuals (transformed standardized
#'   residuals). These are the result of the PIT transformation followed by
#'   inverse normal/t quantile transform.
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional covariance matrix of z_matrix
#' @param copula_dist Copula distribution: "mvn" or "mvt"
#' @param use_reparam Logical: parameters in (psi, phi) space?
#' @param method SE method: "hessian" (default) or "sandwich"
#'
#' @return List with:
#'   \item{se}{Standard errors}
#'   \item{vcov}{Variance-covariance matrix}
#'   \item{params}{Parameter values}
#'   \item{param_names}{Parameter names}
#'   \item{method}{Method used}
#'   \item{hessian}{Hessian matrix (if method = "hessian")}
#'   \item{eigenvalues}{Hessian eigenvalues (if method = "hessian")}
#'
#' @details
#' The copula log-likelihood differs from DCC in that it subtracts the marginal
#' log-densities. For MVN copula:
#' \deqn{\ell_t = -0.5[\log|R_t| + z_t'R_t^{-1}z_t - z_t'z_t]}
#'
#' For MVT copula:

#' \deqn{\ell_t = c(\nu) - 0.5\log|R_t| - \frac{\nu+k}{2}\log(1 + q_t/(\nu-2)) - \sum_j \log f_t(z_{jt})}
#' where \eqn{f_t} is the marginal Student-t density.
#'
#' @seealso \code{\link{estimate_garch_weighted_cgarch}}, \code{\link{copula_nll}}
#'
#' @export
cgarch_standard_errors <- function(
    params,
    z_matrix,
    weights,
    Qbar,
    copula_dist = "mvn",
    use_reparam = FALSE,
    method = c("hessian", "sandwich")
) {
  method <- match.arg(method)
  
  n_params <- length(params)
  param_names <- if (n_params == 2) {
    if (use_reparam) c("psi", "phi") else c("alpha", "beta")
  } else if (n_params == 3) {
    if (use_reparam) c("psi", "phi", "shape") else c("alpha", "beta", "shape")
  } else {
    paste0("param", 1:n_params)
  }
  
  if (method == "sandwich") {
    return(cgarch_sandwich_se(params, z_matrix, weights, Qbar, copula_dist, use_reparam))
  }
  
  ## Compute numerical Hessian of the copula NLL
  H <- numerical_hessian_richardson(
    fn = cgarch_nll_for_hessian,
    params = params,
    eps = 1e-5,
    z_matrix = z_matrix,
    weights = weights,
    Qbar = Qbar,
    copula_dist = copula_dist,
    use_reparam = use_reparam
  )
  rownames(H) <- colnames(H) <- param_names
  
  ## Compute eigendecomposition and inverse
  eig <- eigen(H, symmetric = TRUE)
  vcov <- tryCatch({
    if (all(eig$values > 1e-10)) {
      solve(H)
    } else {
      warning("CGARCH Hessian is near-singular, using pseudo-inverse")
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


#' @title CGARCH Negative Log-Likelihood for Hessian Computation
#' @description Wrapper around copula NLL for use with numerical_hessian.
#' @keywords internal
cgarch_nll_for_hessian <- function(
    params,
    z_matrix,
    weights,
    Qbar,
    copula_dist = "mvn",
    use_reparam = FALSE
) {
  T_obs <- nrow(z_matrix)
  k <- ncol(z_matrix)
  
  ## Extract alpha, beta
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
  
  shape <- if (copula_dist == "mvt" && length(params) >= 3) {
    as.numeric(params[3])
  } else {
    NULL
  }
  
  ## Check constraints
  if (!is.finite(alpha) || !is.finite(beta) ||
      alpha + beta >= 1 || alpha < 0 || beta < 0) {
    return(1e10)
  }
  
  if (copula_dist == "mvt" && (!is.null(shape) && shape <= 2)) {
    return(1e10)
  }
  
  ## Run DCC recursion (CGARCH uses same Q dynamics as DCC)
  fwd <- dcc11_recursion_with_grad(z_matrix, alpha, beta, Qbar)
  
  if (!isTRUE(fwd$success)) {
    return(1e10)
  }
  
  ## Compute copula log-likelihood
  nll <- 0
  
  for (t in 1:T_obs) {
    z_t <- z_matrix[t, ]
    R_inv_t <- fwd$R_inv[,,t]
    log_det_R_t <- fwd$log_det_R[t]
    
    q_t <- as.numeric(t(z_t) %*% R_inv_t %*% z_t)
    
    if (copula_dist == "mvn") {
      ## Gaussian copula: -0.5 * (log|R| + z'R^-1 z - z'z)
      ll_t <- -0.5 * (log_det_R_t + q_t - sum(z_t^2))
    } else {
      ## Student-t copula
      const_term <- lgamma((shape + k) / 2) - lgamma(shape / 2) -
        (k / 2) * log(pi * (shape - 2))
      
      ll_t <- const_term - 0.5 * log_det_R_t -
        ((shape + k) / 2) * log(1 + q_t / (shape - 2))
      
      ## Subtract marginal t log-densities
      scale <- sqrt(shape / (shape - 2))
      for (j in 1:k) {
        ll_t <- ll_t - (dt(z_t[j] * scale, df = shape, log = TRUE) + log(scale))
      }
    }
    
    nll <- nll - weights[t] * ll_t
  }
  
  nll
}


#' @title CGARCH Sandwich (Robust) Standard Errors
#' @description Computes sandwich estimator for CGARCH parameters.
#' @param params Parameter vector
#' @param z_matrix Copula residuals
#' @param weights Observation weights
#' @param Qbar Unconditional covariance
#' @param copula_dist "mvn" or "mvt"
#' @param use_reparam Use reparameterization?
#' @return List with SE results
#' @keywords internal
cgarch_sandwich_se <- function(
    params,
    z_matrix,
    weights,
    Qbar,
    copula_dist = "mvn",
    use_reparam = FALSE
) {
  T_obs <- nrow(z_matrix)
  n_params <- length(params)
  
  ## First get Hessian-based result
  hess_result <- cgarch_standard_errors(
    params, z_matrix, weights, Qbar, copula_dist, use_reparam, method = "hessian"
  )
  
  if (is.null(hess_result$vcov) || any(is.na(hess_result$vcov))) {
    warning("Cannot compute CGARCH sandwich SE: information matrix is singular")
    hess_result$method <- "sandwich"
    return(hess_result)
  }
  
  ## Compute meat matrix (outer product of score contributions)
  J <- matrix(0, n_params, n_params)
  
  for (t in 1:T_obs) {
    grad_t <- .cgarch_gradient_single_obs(
      params, z_matrix, weights, Qbar, t, copula_dist, use_reparam, eps = 1e-6
    )
    J <- J + weights[t]^2 * outer(grad_t, grad_t)
  }
  
  ## Sandwich: bread %*% meat %*% bread
  bread <- hess_result$vcov
  sandwich <- bread %*% J %*% bread
  var_diag <- diag(sandwich)
  var_diag[var_diag < 0] <- NA_real_
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
.cgarch_gradient_single_obs <- function(
    params,
    z_matrix,
    weights,
    Qbar,
    t,
    copula_dist,
    use_reparam,
    eps = 1e-6
) {
  n_params <- length(params)
  grad <- numeric(n_params)
  
  for (i in 1:n_params) {
    p_plus <- p_minus <- params
    p_plus[i] <- params[i] + eps
    p_minus[i] <- params[i] - eps
    
    grad[i] <- (
      .cgarch_nll_single_obs(p_plus, z_matrix, weights, Qbar, t, copula_dist, use_reparam) -
        .cgarch_nll_single_obs(p_minus, z_matrix, weights, Qbar, t, copula_dist, use_reparam)
    ) / (2 * eps)
  }
  
  grad
}


#' @keywords internal
.cgarch_nll_single_obs <- function(
    params,
    z_matrix,
    weights,
    Qbar,
    t,
    copula_dist,
    use_reparam
) {
  k <- ncol(z_matrix)
  
  ## Extract parameters
  if (use_reparam) {
    orig <- dcc11_from_unconstrained(params[1], params[2])
    alpha <- as.numeric(orig["alpha"])
    beta <- as.numeric(orig["beta"])
  } else {
    alpha <- as.numeric(params[1])
    beta <- as.numeric(params[2])
  }
  
  shape <- if (copula_dist == "mvt" && length(params) >= 3) {
    as.numeric(params[3])
  } else {
    NULL
  }
  
  if (!is.finite(alpha) || !is.finite(beta) ||
      alpha + beta >= 1 || alpha < 0 || beta < 0) {
    return(1e10)
  }
  
  ## Run recursion
  fwd <- dcc11_recursion_with_grad(z_matrix, alpha, beta, Qbar)
  if (!isTRUE(fwd$success)) return(1e10)
  
  z_t <- z_matrix[t, ]
  R_inv_t <- fwd$R_inv[,,t]
  log_det_R_t <- fwd$log_det_R[t]
  q_t <- as.numeric(t(z_t) %*% R_inv_t %*% z_t)
  
  if (copula_dist == "mvn") {
    ll_t <- -0.5 * (log_det_R_t + q_t - sum(z_t^2))
  } else {
    const_term <- lgamma((shape + k) / 2) - lgamma(shape / 2) -
      (k / 2) * log(pi * (shape - 2))
    ll_t <- const_term - 0.5 * log_det_R_t -
      ((shape + k) / 2) * log(1 + q_t / (shape - 2))
    
    scale <- sqrt(shape / (shape - 2))
    for (j in 1:k) {
      ll_t <- ll_t - (dt(z_t[j] * scale, df = shape, log = TRUE) + log(scale))
    }
  }
  
  -weights[t] * ll_t
}


#' @title Robust CGARCH Standard Errors with Edge Case Handling
#'
#' @description High-level wrapper for CGARCH SE computation with handling for
#'   boundary estimates and degenerate cases.
#'
#' @param cgarch_result Result from estimate_garch_weighted_cgarch() or similar
#' @param z_matrix Copula residuals (T x k)
#' @param weights Observation weights
#' @param Qbar Unconditional covariance
#' @param copula_dist "mvn" or "mvt"
#' @param boundary_threshold Threshold for boundary detection
#' @param method "hessian" or "sandwich"
#'
#' @return List with SE results and validity flags
#'
#' @export
compute_cgarch_standard_errors_robust <- function(
    cgarch_result,
    z_matrix,
    weights,
    Qbar,
    copula_dist = "mvn",
    boundary_threshold = 1e-4,
    method = c("hessian", "sandwich")
) {
  method <- match.arg(method)
  result <- list(se = NULL, vcov = NULL, valid = FALSE, reason = NULL,
                 correlation_type = NULL, info = NULL, info_eigenvalues = NULL)
  
  correlation_type <- cgarch_result$correlation_type %||%
    cgarch_result$coefficients$correlation_type %||% "dynamic"
  result$correlation_type <- correlation_type
  
  if (correlation_type == "constant") {
    result$reason <- "constant_correlation"
    result$valid <- TRUE
    return(result)
  }
  
  dcc_pars <- cgarch_result$dcc_pars %||% cgarch_result$coefficients$dcc_pars
  if (is.null(dcc_pars) || length(dcc_pars) == 0) {
    result$reason <- "no_dcc_parameters"
    return(result)
  }
  
  alpha <- dcc_pars$alpha_1 %||% dcc_pars[["alpha_1"]]
  beta <- dcc_pars$beta_1 %||% dcc_pars[["beta_1"]]
  if (is.null(alpha) || is.null(beta)) {
    result$reason <- "missing_alpha_or_beta"
    return(result)
  }
  
  ## Check boundaries
  if (alpha < boundary_threshold || beta < boundary_threshold) {
    result$reason <- sprintf("boundary_lower (alpha=%.2e, beta=%.2e)", alpha, beta)
    result$se <- c(alpha = NA_real_, beta = NA_real_)
    if (copula_dist == "mvt") result$se <- c(result$se, shape = NA_real_)
    return(result)
  }
  if ((alpha + beta) > (1 - boundary_threshold)) {
    result$reason <- sprintf("boundary_upper (persistence=%.6f)", alpha + beta)
    result$se <- c(alpha = NA_real_, beta = NA_real_)
    if (copula_dist == "mvt") result$se <- c(result$se, shape = NA_real_)
    return(result)
  }
  
  ## Build parameter vector
  params <- c(alpha = alpha, beta = beta)
  if (copula_dist == "mvt") {
    shape <- cgarch_result$dist_pars$shape %||%
      cgarch_result$coefficients$dist_pars$shape %||% 8.0
    params <- c(params, shape = shape)
  }
  
  ## Compute SE
  se_result <- tryCatch(
    cgarch_standard_errors(params, z_matrix, weights, Qbar, copula_dist, FALSE, method),
    error = function(e) list(se = NULL, vcov = NULL, info = NULL, error = e$message)
  )
  
  if (!is.null(se_result$error)) {
    result$reason <- paste("computation_error:", se_result$error)
    result$se <- c(alpha = NA_real_, beta = NA_real_)
    if (copula_dist == "mvt") result$se <- c(result$se, shape = NA_real_)
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
  
  if (any(is.na(se_result$se))) {
    result$reason <- "na_standard_errors"
    return(result)
  }
  if (any(se_result$se <= 0)) {
    result$reason <- "non_positive_standard_errors"
    return(result)
  }
  
  result$valid <- TRUE
  result$reason <- "ok"
  result
}


#### ______________________________________________________________________ ####
## SECTION 8: GOGARCH Standard Error Computation ===============================
##
## GOGARCH uses ICA decomposition followed by univariate GARCH on components.
## Standard errors are computed for the component GARCH parameters.
## The ICA mixing matrix is typically treated as fixed after estimation.
##

#' @title Compute Standard Errors for GOGARCH Models
#'
#' @description Computes standard errors for GOGARCH component GARCH parameters.
#'   GOGARCH differs from DCC/CGARCH in that it uses ICA decomposition followed
#'   by univariate GARCH on independent components. SE computation focuses on
#'   the component-level GARCH parameters.
#'
#' @param garch_pars List of GARCH parameters for each component:
#'   Each element is a list with omega, alpha1, beta1, etc.
#' @param ica_info ICA decomposition results (A, W, K matrices, S components)
#' @param residuals Original residuals matrix (T x k)
#' @param weights Observation weights (length T)
#' @param distribution Component distribution: "norm", "std", "nig", "gh"
#' @param method SE method: "hessian" (default) or "sandwich"
#'
#' @return List with:
#'   \item{component_se}{List of SE for each component}
#'   \item{vcov_blocks}{Block-diagonal vcov matrix (component-wise)}
#'   \item{valid}{Logical: all SEs computed successfully}
#'   \item{n_components}{Number of components}
#'   \item{method}{Method used}
#'
#' @details
#' GOGARCH models the observation vector as: Y = A * S, where S contains
#' independent components each following univariate GARCH. The log-likelihood
#' decomposes as:
#' \deqn{LL = \sum_i LL_i(S_i; \theta_i) + \log|det(K)|}
#'
#' Standard errors are computed independently for each component's GARCH
#' parameters using the component-wise Hessian. This is justified by the
#' independence assumption of ICA.
#'
#' Note: SEs for the ICA mixing matrix A are not provided as A is typically
#' treated as a fixed transformation after ICA estimation.
#'
#' @seealso \code{\link{estimate_garch_weighted_gogarch}}, 
#'   \code{\link{compute_gogarch_loglik_ms}}
#'
#' @export
gogarch_standard_errors <- function(
    garch_pars,
    ica_info,
    residuals,
    weights,
    distribution = "norm",
    method = c("hessian", "sandwich")
) {
  method <- match.arg(method)
  
  n_components <- length(garch_pars)
  T_obs <- nrow(residuals)
  
  ## Transform residuals to independent components
  S <- residuals %*% t(ica_info$W)
  
  ## Compute SE for each component independently
  component_se_list <- list()
  vcov_list <- list()
  all_valid <- TRUE
  
  for (i in 1:n_components) {
    pars_i <- garch_pars[[i]]
    component_i <- S[, i]
    
    se_i <- tryCatch({
      gogarch_component_se(
        pars = pars_i,
        component = component_i,
        weights = weights,
        distribution = distribution,
        method = method
      )
    }, error = function(e) {
      list(se = NULL, vcov = NULL, valid = FALSE, error = e$message)
    })
    
    component_se_list[[i]] <- se_i
    
    if (!is.null(se_i$vcov)) {
      vcov_list[[i]] <- se_i$vcov
    }
    
    if (!isTRUE(se_i$valid)) {
      all_valid <- FALSE
    }
  }
  
  list(
    component_se = component_se_list,
    vcov_blocks = vcov_list,
    valid = all_valid,
    n_components = n_components,
    method = method
  )
}


#' @title Compute SE for Single GOGARCH Component
#' @description Computes SE for univariate GARCH parameters on one ICA component.
#' @param pars GARCH parameters (omega, alpha1, beta1, possibly shape, skew)
#' @param component Vector of ICA component values
#' @param weights Observation weights
#' @param distribution Component distribution
#' @param method "hessian" or "sandwich"
#' @return List with se, vcov, valid flag
#' @keywords internal
gogarch_component_se <- function(
    pars,
    component,
    weights,
    distribution = "norm",
    method = "hessian"
) {
  ## Convert list pars to named vector
  if (is.list(pars)) {
    param_vec <- unlist(pars)
  } else {
    param_vec <- pars
  }
  
  n_params <- length(param_vec)
  param_names <- names(param_vec)
  if (is.null(param_names)) {
    param_names <- paste0("param", 1:n_params)
  }
  
  ## Define component NLL function
  component_nll <- function(params, resid, w, dist) {
    if (length(params) < 3) {
      return(1e10)
    }
    
    omega <- params[1]
    
    ## Find alpha and beta parameters (may be indexed)
    alpha_idx <- grep("alpha", names(params))
    beta_idx <- grep("beta", names(params))
    
    if (length(alpha_idx) == 0) {
      alpha <- params[2]
      alpha_idx <- 2
    } else {
      alpha <- params[alpha_idx]
    }
    
    if (length(beta_idx) == 0) {
      beta <- params[3]
      beta_idx <- 3
    } else {
      beta <- params[beta_idx]
    }
    
    ## Check constraints
    if (omega <= 0 || any(alpha < 0) || any(beta < 0)) {
      return(1e10)
    }
    if (sum(alpha) + sum(beta) >= 1) {
      return(1e10)
    }
    
    ## Compute GARCH variance path
    n <- length(resid)
    sigma2 <- rep(var(resid), n)
    p <- length(alpha)
    q <- length(beta)
    maxpq <- max(p, q, 1)
    
    for (t in (maxpq + 1):n) {
      sigma2[t] <- omega
      for (j in 1:p) {
        sigma2[t] <- sigma2[t] + alpha[j] * resid[t - j]^2
      }
      for (j in 1:q) {
        sigma2[t] <- sigma2[t] + beta[j] * sigma2[t - j]
      }
      if (sigma2[t] <= 0) sigma2[t] <- 1e-10
    }
    
    sig <- sqrt(sigma2)
    
    ## Compute log-likelihood
    if (dist == "norm") {
      ll <- dnorm(resid, mean = 0, sd = sig, log = TRUE)
    } else if (dist == "std") {
      shape <- if (length(params) > 3) params["shape"] else 4
      if (is.na(shape) || shape <= 2) shape <- 4
      ll <- tryCatch(
        tsdistributions::dstd(resid, mu = 0, sigma = sig, shape = shape, log = TRUE),
        error = function(e) dnorm(resid, mean = 0, sd = sig, log = TRUE)
      )
    } else {
      ll <- dnorm(resid, mean = 0, sd = sig, log = TRUE)
    }
    
    ll[!is.finite(ll)] <- -1e10
    -sum(w * ll, na.rm = TRUE)
  }
  
  ## Compute numerical Hessian
  H <- tryCatch({
    numerical_hessian_richardson(
      fn = component_nll,
      params = param_vec,
      eps = 1e-5,
      resid = component,
      w = weights,
      dist = distribution
    )
  }, error = function(e) {
    matrix(NA_real_, n_params, n_params)
  })
  
  if (any(is.na(H))) {
    return(list(
      se = rep(NA_real_, n_params),
      vcov = H,
      valid = FALSE,
      reason = "hessian_computation_failed"
    ))
  }
  
  rownames(H) <- colnames(H) <- param_names
  
  ## Invert Hessian
  eig <- eigen(H, symmetric = TRUE)
  vcov <- tryCatch({
    if (all(eig$values > 1e-10)) {
      solve(H)
    } else {
      eig$vectors %*% diag(1 / pmax(eig$values, 1e-10)) %*% t(eig$vectors)
    }
  }, error = function(e) matrix(NA_real_, n_params, n_params))
  
  rownames(vcov) <- colnames(vcov) <- param_names
  
  se <- sqrt(pmax(0, diag(vcov)))
  names(se) <- param_names
  
  ## Validate
  valid <- !any(is.na(se)) && all(se > 0) && all(eig$values > 0)
  
  list(
    se = se,
    vcov = vcov,
    hessian = H,
    eigenvalues = eig$values,
    valid = valid,
    reason = if (valid) "ok" else "numerical_issues",
    method = method
  )
}


#' @title Robust GOGARCH Standard Errors with Edge Case Handling
#'
#' @description High-level wrapper for GOGARCH SE computation with validation
#'   and error handling.
#'
#' @param gogarch_result Result from estimate_garch_weighted_gogarch()
#' @param residuals Original residuals (T x k)
#' @param weights Observation weights
#' @param distribution Component distribution
#' @param method "hessian" or "sandwich"
#'
#' @return List with SE results and validity information
#'
#' @export
compute_gogarch_standard_errors_robust <- function(
    gogarch_result,
    residuals,
    weights,
    distribution = "norm",
    method = c("hessian", "sandwich")
) {
  method <- match.arg(method)
  
  result <- list(
    component_se = NULL,
    valid = FALSE,
    reason = NULL,
    n_components = NULL
  )
  
  ## Extract coefficients
  coeffs <- gogarch_result$coefficients %||% gogarch_result
  
  garch_pars <- coeffs$garch_pars
  if (is.null(garch_pars) || length(garch_pars) == 0) {
    result$reason <- "no_garch_parameters"
    return(result)
  }
  
  ica_info <- coeffs$ica_info
  if (is.null(ica_info)) {
    result$reason <- "no_ica_info"
    return(result)
  }
  
  ## Compute SE
  se_result <- tryCatch(
    gogarch_standard_errors(
      garch_pars = garch_pars,
      ica_info = ica_info,
      residuals = residuals,
      weights = weights,
      distribution = distribution,
      method = method
    ),
    error = function(e) {
      list(component_se = NULL, valid = FALSE, error = e$message)
    }
  )
  
  if (!is.null(se_result$error)) {
    result$reason <- paste("computation_error:", se_result$error)
    return(result)
  }
  
  result$component_se <- se_result$component_se
  result$vcov_blocks <- se_result$vcov_blocks
  result$n_components <- se_result$n_components
  result$valid <- se_result$valid
  result$reason <- if (se_result$valid) "ok" else "some_components_failed"
  result$method <- method
  
  result
}


#' @title GOGARCH Estimation Summary
#'
#' @description Provides a summary of GOGARCH estimation including parameter
#'   estimates and standard errors for all components.
#'
#' @param gogarch_result Result from estimate_garch_weighted_gogarch()
#' @param residuals Original residuals
#' @param weights Observation weights
#' @param distribution Component distribution
#' @param level Confidence level
#' @param method SE method
#'
#' @return Object of class "gogarch_summary"
#'
#' @export
gogarch_estimation_summary <- function(
    gogarch_result,
    residuals,
    weights,
    distribution = "norm",
    level = 0.95,
    method = c("hessian", "sandwich")
) {
  method <- match.arg(method)
  
  coeffs <- gogarch_result$coefficients %||% gogarch_result
  garch_pars <- coeffs$garch_pars
  ica_info <- coeffs$ica_info
  n_components <- length(garch_pars)
  
  ## Compute SE
  se_result <- compute_gogarch_standard_errors_robust(
    gogarch_result, residuals, weights, distribution, method
  )
  
  ## Build summary for each component
  z <- qnorm(1 - (1 - level) / 2)
  component_summaries <- list()
  
  for (i in 1:n_components) {
    pars_i <- garch_pars[[i]]
    se_i <- se_result$component_se[[i]]
    
    if (!is.null(se_i) && isTRUE(se_i$valid)) {
      param_names <- names(se_i$se)
      estimates <- unlist(pars_i)[param_names]
      ci <- cbind(
        estimate = estimates,
        se = se_i$se,
        lower = estimates - z * se_i$se,
        upper = estimates + z * se_i$se
      )
      rownames(ci) <- param_names
    } else {
      estimates <- unlist(pars_i)
      ci <- cbind(
        estimate = estimates,
        se = NA,
        lower = NA,
        upper = NA
      )
    }
    
    component_summaries[[i]] <- list(
      params = pars_i,
      ci = ci,
      valid = isTRUE(se_i$valid)
    )
  }
  
  result <- list(
    n_obs = nrow(residuals),
    n_series = ncol(residuals),
    n_components = n_components,
    distribution = distribution,
    ica_method = ica_info$method %||% "unknown",
    component_summaries = component_summaries,
    se_method = method,
    ci_level = level,
    overall_valid = se_result$valid
  )
  
  class(result) <- "gogarch_summary"
  result
}


#' @export
print.gogarch_summary <- function(x, digits = 4, ...) {
  cat("\n===== GOGARCH Estimation Summary =====\n\n")
  cat("Data:\n")
  cat(sprintf("  Observations: %d\n  Series: %d\n  Components: %d\n",
              x$n_obs, x$n_series, x$n_components))
  cat(sprintf("  Distribution: %s\n  ICA method: %s\n\n",
              x$distribution, x$ica_method))
  
  cat(sprintf("Component Parameter Estimates (SE method: %s):\n", x$se_method))
  cat(sprintf("%.0f%% Confidence Intervals:\n\n", x$ci_level * 100))
  
  for (i in seq_along(x$component_summaries)) {
    cat(sprintf("--- Component %d ---\n", i))
    comp <- x$component_summaries[[i]]
    
    if (comp$valid) {
      print(round(comp$ci, digits))
    } else {
      cat("  SE computation failed\n")
      cat("  Estimates:", paste(names(comp$params), "=",
                                round(unlist(comp$params), digits), collapse = ", "), "\n")
    }
    cat("\n")
  }
  
  if (!x$overall_valid) {
    cat("WARNING: SE computation failed for some components.\n\n")
  }
  
  invisible(x)
}


#### ______________________________________________________________________ ####
## SECTION 9: ADCC Standard Errors and Inference ===============================

## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## This section implements Hessian-based standard errors for ADCC (Asymmetric DCC)
## models. ADCC adds a gamma parameter capturing asymmetric response to negative
## shocks, making it a 3-parameter model: (alpha, gamma, beta).
##
## The ADCC(1,1) model is:
##   Q_t = (1 - α - β)Q̄ - γN̄ + α(z_{t-1}z'_{t-1}) + γ(n_{t-1}n'_{t-1}) + βQ_{t-1}
##
## where:
##   - n_t = z_t * I(z_t < 0) captures negative shocks
##   - N̄ = E[n_t n'_t] is the average outer product of negative shocks
##   - Stationarity requires: α + β + δγ < 1 (where δ ≈ 0.5 typically)
##
## Dependencies:
##   - tsbs_cgarch.R (adcc_copula_nll, adcc_recursion, estimate_adcc_copula)
##   - hessian_se.R (numerical_hessian_richardson)
##
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


## SECTION 9a: ADCC Hessian-Based Standard Errors ==============================

#' @title Compute Standard Errors for ADCC Parameters
#' @description Computes Hessian-based (or sandwich) standard errors for the
#'   3-parameter ADCC model: (alpha, gamma, beta), optionally with shape for MVT.
#'
#' @param params Parameter vector: c(alpha, gamma, beta) or c(alpha, gamma, beta, shape)
#' @param z_matrix T x k matrix of copula residuals (PIT-transformed)
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional correlation matrix
#' @param Nbar k x k average outer product of negative residuals (optional)
#' @param copula_dist "mvn" or "mvt"
#' @param method SE method: "hessian" (default) or "sandwich"
#'
#' @return List with:
#'   \item{se}{Standard errors for each parameter}
#'   \item{vcov}{Variance-covariance matrix}
#'   \item{hessian}{Hessian matrix of NLL}
#'   \item{eigenvalues}{Eigenvalues of Hessian (for diagnostics)}
#'   \item{condition_number}{Condition number of Hessian}
#'   \item{param_names}{Parameter names}
#'   \item{method}{"hessian" or "sandwich"}
#'   \item{valid}{Logical: are SEs valid?}
#'
#' @details
#' For ADCC, the gamma parameter introduces asymmetric correlation dynamics.
#' The Hessian-based SEs may underestimate uncertainty for high-persistence
#' models (similar to standard DCC). Bootstrap SEs are recommended for
#' publication-quality inference.
#'
#' @export
adcc_standard_errors <- function(
    params,
    z_matrix,
    weights,
    Qbar,
    Nbar = NULL,
    copula_dist = "mvn",
    method = c("hessian", "sandwich")
) {
  method <- match.arg(method)
  
  n_params <- length(params)
  T_obs <- nrow(z_matrix)
  k <- ncol(z_matrix)
  
  ## Determine parameter names
  param_names <- if (n_params == 3) {
    c("alpha", "gamma", "beta")
  } else if (n_params == 4) {
    c("alpha", "gamma", "beta", "shape")
  } else {
    paste0("param", 1:n_params)
  }
  
  ## Compute Nbar if not provided
  if (is.null(Nbar)) {
    neg_resid <- z_matrix * (z_matrix < 0)
    Nbar <- crossprod(neg_resid) / T_obs
  }
  
  if (method == "sandwich") {
    return(adcc_sandwich_se(params, z_matrix, weights, Qbar, Nbar, copula_dist))
  }
  
  ## Compute numerical Hessian of ADCC NLL
  H <- numerical_hessian_richardson(
    fn = adcc_copula_nll,
    params = params,
    eps = 1e-5,
    z_matrix = z_matrix,
    weights = weights,
    Qbar = Qbar,
    Nbar = Nbar,
    copula_dist = copula_dist
  )
  rownames(H) <- colnames(H) <- param_names
  
  ## Compute eigendecomposition and inverse
  eig <- eigen(H, symmetric = TRUE)
  
  vcov <- tryCatch({
    if (all(eig$values > 1e-10)) {
      solve(H)
    } else {
      warning("ADCC Hessian is near-singular, using pseudo-inverse")
      eig$vectors %*% diag(1 / pmax(eig$values, 1e-10)) %*% t(eig$vectors)
    }
  }, error = function(e) matrix(NA_real_, n_params, n_params))
  rownames(vcov) <- colnames(vcov) <- param_names
  
  se <- sqrt(pmax(0, diag(vcov)))
  names(se) <- param_names
  
  ## Validate results
  valid <- !any(is.na(se)) && all(se > 0) && all(eig$values > 0)
  
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
    method = "hessian",
    valid = valid,
    reason = if (valid) "ok" else "numerical_issues"
  )
}


#' @title ADCC Sandwich (Robust) Standard Errors
#' @description Computes sandwich estimator for ADCC parameters.
#' @param params Parameter vector
#' @param z_matrix Copula residuals
#' @param weights Observation weights
#' @param Qbar Unconditional covariance
#' @param Nbar Average outer product of negative residuals
#' @param copula_dist "mvn" or "mvt"
#' @return List with SE results
#' @keywords internal
adcc_sandwich_se <- function(
    params,
    z_matrix,
    weights,
    Qbar,
    Nbar,
    copula_dist = "mvn"
) {
  T_obs <- nrow(z_matrix)
  n_params <- length(params)
  
  ## First get Hessian-based result
  hess_result <- adcc_standard_errors(
    params, z_matrix, weights, Qbar, Nbar, copula_dist, method = "hessian"
  )
  
  if (is.null(hess_result$vcov) || any(is.na(hess_result$vcov))) {
    warning("Cannot compute ADCC sandwich SE: information matrix is singular")
    hess_result$method <- "sandwich"
    return(hess_result)
  }
  
  ## Compute meat matrix (outer product of score contributions)
  J <- matrix(0, n_params, n_params)
  
  for (t in 1:T_obs) {
    grad_t <- .adcc_gradient_single_obs(
      params, z_matrix, weights, Qbar, Nbar, t, copula_dist, eps = 1e-6
    )
    J <- J + weights[t]^2 * outer(grad_t, grad_t)
  }
  
  ## Sandwich: bread %*% meat %*% bread
  bread <- hess_result$vcov
  sandwich <- bread %*% J %*% bread
  var_diag <- diag(sandwich)
  var_diag[var_diag < 0] <- NA_real_
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
    method = "sandwich",
    valid = !any(is.na(se_sandwich)) && all(se_sandwich > 0),
    reason = if (!any(is.na(se_sandwich)) && all(se_sandwich > 0)) "ok" else "numerical_issues"
  )
}


#' @keywords internal
.adcc_gradient_single_obs <- function(
    params,
    z_matrix,
    weights,
    Qbar,
    Nbar,
    t,
    copula_dist,
    eps = 1e-6
) {
  n_params <- length(params)
  grad <- numeric(n_params)
  
  for (i in 1:n_params) {
    p_plus <- p_minus <- params
    p_plus[i] <- params[i] + eps
    p_minus[i] <- params[i] - eps
    
    grad[i] <- (
      .adcc_nll_single_obs(p_plus, z_matrix, weights, Qbar, Nbar, t, copula_dist) -
        .adcc_nll_single_obs(p_minus, z_matrix, weights, Qbar, Nbar, t, copula_dist)
    ) / (2 * eps)
  }
  
  grad
}


#' @keywords internal
.adcc_nll_single_obs <- function(
    params,
    z_matrix,
    weights,
    Qbar,
    Nbar,
    t,
    copula_dist
) {
  k <- ncol(z_matrix)
  
  alpha <- as.numeric(params[1])
  gamma <- as.numeric(params[2])
  beta <- as.numeric(params[3])
  shape <- if (copula_dist == "mvt" && length(params) >= 4) as.numeric(params[4]) else NULL
  
  ## Check constraints using ADCC stationarity
  stationarity <- alpha + beta + 0.5 * gamma  # Approximate stationarity
  if (!is.finite(alpha) || !is.finite(gamma) || !is.finite(beta) ||
      stationarity >= 1 || alpha < 0 || gamma < 0 || beta < 0) {
    return(1e10)
  }
  
  if (copula_dist == "mvt" && (!is.null(shape) && shape <= 2)) {
    return(1e10)
  }
  
  ## Run ADCC recursion
  recursion <- adcc_recursion(
    std_resid = z_matrix,
    Qbar = Qbar,
    alpha = alpha,
    gamma = gamma,
    beta = beta,
    Nbar = Nbar
  )
  
  if (!isTRUE(recursion$success)) {
    return(1e10)
  }
  
  ## Extract R_t and its inverse for observation t
  R_t <- recursion$R[, , t]
  R_inv_t <- tryCatch(solve(R_t), error = function(e) {
    eig <- eigen(R_t, symmetric = TRUE)
    eig$vectors %*% diag(1 / pmax(eig$values, 1e-8)) %*% t(eig$vectors)
  })
  log_det_R_t <- log(max(det(R_t), 1e-100))
  
  z_t <- z_matrix[t, ]
  q_t <- as.numeric(t(z_t) %*% R_inv_t %*% z_t)
  
  if (copula_dist == "mvn") {
    ll_t <- -0.5 * (log_det_R_t + q_t - sum(z_t^2))
  } else {
    const_term <- lgamma((shape + k) / 2) - lgamma(shape / 2) -
      (k / 2) * log(pi * (shape - 2))
    
    ll_t <- const_term - 0.5 * log_det_R_t -
      ((shape + k) / 2) * log(1 + q_t / (shape - 2))
    
    ## Subtract marginal t log-densities
    scale <- sqrt(shape / (shape - 2))
    for (j in 1:k) {
      ll_t <- ll_t - (dt(z_t[j] * scale, df = shape, log = TRUE) + log(scale))
    }
  }
  
  -weights[t] * ll_t
}


## SECTION 9b: Robust ADCC SE with Edge Case Handling ==========================

#' @title Compute ADCC Standard Errors with Robust Edge Case Handling
#' @description High-level wrapper for ADCC SE computation with validation
#'   and edge case handling (boundary estimates, near-singular Hessians).
#'
#' @param adcc_result Result from estimate_adcc_copula() or similar
#' @param z_matrix T x k matrix of copula residuals
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional correlation matrix
#' @param copula_dist "mvn" or "mvt"
#' @param boundary_threshold Threshold for boundary detection (default 1e-4)
#' @param method "hessian" or "sandwich"
#'
#' @return List with SE results and validity information
#'
#' @export
compute_adcc_standard_errors_robust <- function(
    adcc_result,
    z_matrix,
    weights,
    Qbar,
    copula_dist = "mvn",
    boundary_threshold = 1e-4,
    method = c("hessian", "sandwich")
) {
  method <- match.arg(method)
  
  result <- list(
    se = NULL, vcov = NULL, valid = FALSE, reason = NULL,
    info = NULL, info_eigenvalues = NULL
  )
  
  ## Extract parameters from adcc_result
  if (is.list(adcc_result)) {
    alpha <- adcc_result$alpha %||% adcc_result$dcc_pars$alpha_1
    gamma <- adcc_result$gamma %||% adcc_result$dcc_pars$gamma_1
    beta <- adcc_result$beta %||% adcc_result$dcc_pars$beta_1
    shape <- adcc_result$shape %||% adcc_result$dist_pars$shape
    Nbar <- adcc_result$Nbar
  } else {
    result$reason <- "invalid_adcc_result"
    return(result)
  }
  
  if (is.null(alpha) || is.null(gamma) || is.null(beta)) {
    result$reason <- "missing_adcc_parameters"
    result$se <- c(alpha = NA_real_, gamma = NA_real_, beta = NA_real_)
    return(result)
  }
  
  ## Check boundary conditions
  if (alpha < boundary_threshold || gamma < boundary_threshold || beta < boundary_threshold) {
    result$reason <- sprintf("boundary_lower (alpha=%.2e, gamma=%.2e, beta=%.2e)", 
                             alpha, gamma, beta)
    result$se <- c(alpha = NA_real_, gamma = NA_real_, beta = NA_real_)
    if (copula_dist == "mvt") result$se <- c(result$se, shape = NA_real_)
    return(result)
  }
  
  stationarity <- alpha + beta + 0.5 * gamma
  if (stationarity > (1 - boundary_threshold)) {
    result$reason <- sprintf("boundary_upper (stationarity=%.6f)", stationarity)
    result$se <- c(alpha = NA_real_, gamma = NA_real_, beta = NA_real_)
    if (copula_dist == "mvt") result$se <- c(result$se, shape = NA_real_)
    return(result)
  }
  
  ## Build parameter vector
  params <- c(alpha = alpha, gamma = gamma, beta = beta)
  if (copula_dist == "mvt" && !is.null(shape)) {
    params <- c(params, shape = shape)
  }
  
  ## Compute Nbar if not available
  if (is.null(Nbar)) {
    T_obs <- nrow(z_matrix)
    neg_resid <- z_matrix * (z_matrix < 0)
    Nbar <- crossprod(neg_resid) / T_obs
  }
  
  ## Compute SEs
  se_result <- tryCatch(
    adcc_standard_errors(params, z_matrix, weights, Qbar, Nbar, copula_dist, method),
    error = function(e) list(se = NULL, vcov = NULL, info = NULL, error = e$message)
  )
  
  if (!is.null(se_result$error)) {
    result$reason <- paste("computation_error:", se_result$error)
    result$se <- c(alpha = NA_real_, gamma = NA_real_, beta = NA_real_)
    if (copula_dist == "mvt") result$se <- c(result$se, shape = NA_real_)
    return(result)
  }
  
  result$se <- se_result$se
  result$vcov <- se_result$vcov
  result$info <- se_result$info
  
  if (!is.null(se_result$eigenvalues)) {
    result$info_eigenvalues <- se_result$eigenvalues
    if (any(se_result$eigenvalues <= 0)) {
      result$reason <- sprintf("non_pd_information (min_eig=%.2e)", min(se_result$eigenvalues))
      return(result)
    }
  }
  
  if (any(is.na(se_result$se))) {
    result$reason <- "na_standard_errors"
    return(result)
  }
  
  if (any(se_result$se <= 0)) {
    result$reason <- "non_positive_standard_errors"
    return(result)
  }
  
  result$valid <- TRUE
  result$reason <- "ok"
  result$method <- method
  result
}


## SECTION 9c: ADCC Bootstrap Standard Errors ==================================

#' @title Bootstrap Standard Errors for ADCC Parameters
#' @description Compute bootstrap standard errors for ADCC correlation
#'   parameters (alpha, gamma, beta) using residual or parametric resampling.
#'
#' @param z_matrix T x k matrix of copula residuals (PIT-transformed)
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional correlation matrix
#' @param mle_params MLE estimates c(alpha, gamma, beta) or with shape for MVT
#' @param Nbar k x k average outer product of negative residuals (optional)
#' @param n_boot Number of bootstrap replications (default 200)
#' @param method Bootstrap method: "residual" or "parametric"
#' @param copula_dist "mvn" or "mvt"
#' @param verbose Print progress
#' @param seed Random seed for reproducibility
#'
#' @return List with bootstrap results including SEs for alpha, gamma, beta
#'
#' @export
adcc_bootstrap_se <- function(
    z_matrix,
    weights,
    Qbar,
    mle_params,
    Nbar = NULL,
    n_boot = 200,
    method = c("residual", "parametric"),
    copula_dist = "mvn",
    verbose = TRUE,
    seed = NULL
) {
  method <- match.arg(method)
  
  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(z_matrix)
  k <- ncol(z_matrix)
  
  alpha_mle <- mle_params[1]
  gamma_mle <- mle_params[2]
  beta_mle <- mle_params[3]
  shape_mle <- if (length(mle_params) > 3) mle_params[4] else NULL
  
  n_params <- if (copula_dist == "mvt" && !is.null(shape_mle)) 4 else 3
  param_names <- if (n_params == 4) c("alpha", "gamma", "beta", "shape") else c("alpha", "gamma", "beta")
  
  ## Compute Nbar if not provided
  if (is.null(Nbar)) {
    neg_resid <- z_matrix * (z_matrix < 0)
    Nbar <- crossprod(neg_resid) / n
  }
  
  boot_estimates <- matrix(NA, n_boot, n_params)
  colnames(boot_estimates) <- param_names
  
  if (verbose) {
    cat(sprintf("Running ADCC %s bootstrap with %d replications...\n", method, n_boot))
    pb <- txtProgressBar(min = 0, max = n_boot, style = 3)
  }
  
  for (b in 1:n_boot) {
    
    if (method == "residual") {
      ## Residual bootstrap: resample rows
      boot_idx <- sample(1:n, n, replace = TRUE)
      boot_z <- z_matrix[boot_idx, , drop = FALSE]
      boot_weights <- weights[boot_idx]
      boot_Qbar <- cor(boot_z)
      
      ## Recompute Nbar for bootstrap sample
      neg_resid_boot <- boot_z * (boot_z < 0)
      boot_Nbar <- crossprod(neg_resid_boot) / n
      
    } else {
      ## Parametric bootstrap: simulate from fitted ADCC model
      boot_z <- simulate_adcc_residuals(
        n = n, k = k,
        alpha = alpha_mle, gamma = gamma_mle, beta = beta_mle,
        Qbar = Qbar, Nbar = Nbar,
        copula_dist = copula_dist, shape = shape_mle
      )
      boot_weights <- weights
      boot_Qbar <- Qbar
      boot_Nbar <- Nbar
    }
    
    ## Re-estimate ADCC parameters
    start_pars <- c(alpha = alpha_mle, gamma = gamma_mle, beta = beta_mle)
    if (copula_dist == "mvt" && !is.null(shape_mle)) {
      start_pars <- c(start_pars, shape = shape_mle)
    }
    
    opt_result <- tryCatch({
      estimate_adcc_copula(
        z_matrix = boot_z,
        weights = boot_weights,
        Qbar = boot_Qbar,
        copula_dist = copula_dist,
        start_pars = start_pars
      )
    }, error = function(e) NULL)
    
    if (!is.null(opt_result) && opt_result$convergence == 0) {
      boot_estimates[b, 1] <- opt_result$alpha
      boot_estimates[b, 2] <- opt_result$gamma
      boot_estimates[b, 3] <- opt_result$beta
      if (n_params == 4) boot_estimates[b, 4] <- opt_result$shape
    }
    
    if (verbose) setTxtProgressBar(pb, b)
  }
  
  if (verbose) {
    close(pb)
    cat("\n")
  }
  
  ## Process results
  boot_estimates_valid <- boot_estimates[complete.cases(boot_estimates), , drop = FALSE]
  n_valid <- nrow(boot_estimates_valid)
  
  if (n_valid < 10) {
    warning(sprintf("Only %d valid bootstrap replications. Results may be unreliable.", n_valid))
  }
  
  boot_se <- apply(boot_estimates_valid, 2, sd)
  boot_mean <- colMeans(boot_estimates_valid)
  boot_bias <- boot_mean - mle_params[1:n_params]
  
  names(boot_se) <- param_names
  names(boot_mean) <- param_names
  names(boot_bias) <- param_names
  
  ## Percentile CIs
  ci_percentile <- apply(boot_estimates_valid, 2, quantile, probs = c(0.025, 0.975))
  colnames(ci_percentile) <- param_names
  
  if (verbose) {
    cat(sprintf("Valid bootstrap replications: %d/%d\n", n_valid, n_boot))
    cat(sprintf("Bootstrap SE: %s\n",
                paste(param_names, "=", round(boot_se, 4), collapse = ", ")))
  }
  
  list(
    se = boot_se,
    boot_estimates = boot_estimates,
    boot_estimates_valid = boot_estimates_valid,
    n_valid = n_valid,
    mean = boot_mean,
    bias = boot_bias,
    ci_percentile = ci_percentile,
    mle_params = mle_params,
    param_names = param_names,
    valid = n_valid >= 10,
    reason = if (n_valid >= 10) "ok" else "insufficient_valid_replicates"
  )
}


#' @title Simulate ADCC Residuals
#' @description Simulate residuals from an ADCC model for parametric bootstrap.
#' @keywords internal
simulate_adcc_residuals <- function(
    n,
    k,
    alpha,
    gamma,
    beta,
    Qbar,
    Nbar,
    copula_dist = "mvn",
    shape = NULL
) {
  ## Initialize Q
  Q <- Qbar
  z <- matrix(0, n, k)
  
  persistence <- alpha + beta + 0.5 * gamma
  Omega <- (1 - persistence) * Qbar - gamma * Nbar
  
  for (t in 1:n) {
    ## Normalize Q to get R
    d <- sqrt(pmax(diag(Q), 1e-10))
    R <- Q / outer(d, d)
    diag(R) <- 1
    R <- (R + t(R)) / 2
    
    ## Ensure R is positive definite
    eig <- eigen(R, symmetric = TRUE)
    if (any(eig$values < 1e-8)) {
      eig$values <- pmax(eig$values, 1e-8)
      R <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
      R <- cov2cor(R)
    }
    
    ## Generate observation
    if (copula_dist == "mvn" || is.null(shape)) {
      z[t, ] <- MASS::mvrnorm(1, mu = rep(0, k), Sigma = R)
    } else {
      ## MVT
      z[t, ] <- mvtnorm::rmvt(1, sigma = R * (shape - 2) / shape, df = shape)
    }
    
    ## Update Q for next period
    if (t < n) {
      z_lag <- z[t, ]
      n_lag <- z_lag * (z_lag < 0)
      Q <- Omega + alpha * outer(z_lag, z_lag) + gamma * outer(n_lag, n_lag) + beta * Q
    }
  }
  
  z
}


## SECTION 9d: ADCC Comprehensive Inference ====================================

#' @title Comprehensive ADCC Inference
#' @description Compute Hessian-based and bootstrap standard errors for ADCC
#'   parameters, with comparison and diagnostics.
#'
#' @param z_matrix T x k matrix of copula residuals
#' @param weights T-vector of observation weights
#' @param adcc_result Result from estimate_adcc_copula()
#' @param Qbar k x k unconditional correlation matrix
#' @param copula_dist "mvn" or "mvt"
#' @param n_boot Number of bootstrap replications
#' @param boot_method "residual" or "parametric"
#' @param conf_level Confidence level
#' @param verbose Print progress
#' @param seed Random seed
#'
#' @return List with comprehensive inference results
#'
#' @export
adcc_comprehensive_inference <- function(
    z_matrix,
    weights,
    adcc_result,
    Qbar = NULL,
    copula_dist = "mvn",
    n_boot = 200,
    boot_method = "residual",
    conf_level = 0.95,
    verbose = TRUE,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(z_matrix)
  k <- ncol(z_matrix)
  
  if (is.null(Qbar)) Qbar <- cor(z_matrix)
  
  ## Extract MLE parameters
  alpha_mle <- adcc_result$alpha
  gamma_mle <- adcc_result$gamma
  beta_mle <- adcc_result$beta
  shape_mle <- adcc_result$shape
  Nbar <- adcc_result$Nbar
  
  mle_params <- c(alpha = alpha_mle, gamma = gamma_mle, beta = beta_mle)
  if (copula_dist == "mvt" && !is.null(shape_mle)) {
    mle_params <- c(mle_params, shape = shape_mle)
  }
  
  if (verbose) {
    cat("=== ADCC Comprehensive Inference ===\n")
    cat(sprintf("Observations: %d, Series: %d\n", n, k))
    cat(sprintf("MLE: %s\n\n", paste(names(mle_params), "=", round(mle_params, 4), collapse = ", ")))
  }
  
  ## Step 1: Hessian-based SEs
  if (verbose) cat("Computing Hessian-based SEs...\n")
  
  hessian_result <- compute_adcc_standard_errors_robust(
    adcc_result = adcc_result,
    z_matrix = z_matrix,
    weights = weights,
    Qbar = Qbar,
    copula_dist = copula_dist,
    method = "hessian"
  )
  
  ## Step 2: Bootstrap SEs
  if (verbose) cat("\nComputing bootstrap SEs...\n")
  
  boot_result <- adcc_bootstrap_se(
    z_matrix = z_matrix,
    weights = weights,
    Qbar = Qbar,
    mle_params = mle_params,
    Nbar = Nbar,
    n_boot = n_boot,
    method = boot_method,
    copula_dist = copula_dist,
    verbose = verbose,
    seed = seed
  )
  
  ## Build comparison table
  z_crit <- qnorm(1 - (1 - conf_level) / 2)
  param_names <- names(mle_params)
  
  comparison <- data.frame(
    Parameter = param_names,
    Estimate = mle_params,
    Hessian_SE = hessian_result$se,
    Boot_SE = boot_result$se,
    Hess_CI_lower = mle_params - z_crit * hessian_result$se,
    Hess_CI_upper = mle_params + z_crit * hessian_result$se,
    Boot_CI_lower = boot_result$ci_percentile[1, ],
    Boot_CI_upper = boot_result$ci_percentile[2, ],
    SE_Ratio = hessian_result$se / boot_result$se,
    row.names = NULL
  )
  
  if (verbose) {
    cat("\n=== Inference Comparison ===\n")
    cat(sprintf("Confidence level: %.0f%%\n\n", conf_level * 100))
    print(comparison, digits = 4)
    
    ## Check for SE underestimation
    ratio <- comparison$SE_Ratio
    if (any(ratio < 0.7, na.rm = TRUE)) {
      cat("\nWARNING: Hessian SEs appear underestimated for some parameters (ratio < 0.7).\n")
      cat("         Bootstrap SEs recommended for inference.\n")
    }
  }
  
  list(
    comparison = comparison,
    mle = mle_params,
    hessian = hessian_result,
    bootstrap = boot_result,
    conf_level = conf_level,
    model_type = "adcc",
    settings = list(
      n = n,
      k = k,
      n_boot = n_boot,
      boot_method = boot_method,
      copula_dist = copula_dist
    )
  )
}


