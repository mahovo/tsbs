## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## MS-VARMA-GARCH Inference: Bootstrap SEs and Profile Likelihood CIs
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
##
## This file provides alternative inference methods for multivariate GARCH
## correlation parameters that do not rely solely on Hessian-based standard
## errors. Supports three correlation model types:
##
##   - DCC (Dynamic Conditional Correlation)
##   - CGARCH (Copula GARCH)
##   - GOGARCH (Generalized Orthogonal GARCH)
##
## Methods provided:
##
##   1. Bootstrap Standard Errors
##      - Parametric bootstrap (resample from fitted model)
##      - Residual bootstrap (resample standardized residuals)
##      - BCa intervals for bias correction
##
##   2. Profile Likelihood Confidence Intervals
##      - Does not assume quadratic likelihood
##      - Based on likelihood ratio test inversion
##      - Handles asymmetric uncertainty correctly
##
##   3. Comprehensive Inference
##      - Computes all three methods (Hessian, Bootstrap, Profile)
##      - Provides comparison and recommendations
##
## WHY THESE METHODS ARE NEEDED
## ────────────────────────────
## Validation studies demonstrate that:
##
##   1. The DCC/CGARCH likelihood is ~8x flatter in beta than alpha direction
##   2. Hessian-based SE for beta is often only ~16% of bootstrap SE
##   3. This is a fundamental property of DCC-type models
##   4. GOGARCH has different issues due to ICA estimation uncertainty
##
## RECOMMENDATIONS
## ───────────────
##   - Alpha: Hessian-based SE is often acceptable
##   - Beta:  ALWAYS use bootstrap or profile likelihood
##   - For publication: Report bootstrap SE and/or profile CI for beta
##   - GOGARCH: Bootstrap recommended due to ICA uncertainty
##
## Dependencies:
##   - dcc_gradient.R (DCC NLL computation)
##   - hessian_se.R (Hessian-based SE computation)
##   - tsbs_cgarch.R (CGARCH functions)
##   - tsbs_gogarch.R (GOGARCH functions)
##
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


#### ______________________________________________________________________ ####
#### SECTION 1: UNIFIED BOOTSTRAP FRAMEWORK                                 ####


#' @title Bootstrap Standard Errors for Correlation Parameters
#' @description Unified bootstrap SE computation for DCC, CGARCH, and GOGARCH
#'   correlation models.
#' @param model_type Type of correlation model: "dcc", "cgarch", or "gogarch"
#' @param residuals T x k matrix of standardized residuals (or copula residuals)
#' @param weights T-vector of observation weights
#' @param mle_result List containing MLE estimates (model-specific structure)
#' @param n_boot Number of bootstrap replications (default 200)
#' @param method Bootstrap method: "residual" or "parametric"
#' @param distribution Distribution: "mvn", "mvt", "norm", "std"
#' @param verbose Print progress
#' @param seed Random seed for reproducibility
#' @param ... Additional model-specific arguments
#' @return List with bootstrap results
#' @export
correlation_bootstrap_se <- function(
    model_type = c("dcc", "cgarch", "gogarch"),
    residuals,
    weights,
    mle_result,
    n_boot = 200,
    method = c("residual", "parametric"),
    distribution = "mvn",
    verbose = TRUE,
    seed = NULL,
    ...
) {
  
  model_type <- match.arg(model_type)
  method <- match.arg(method)
  
  if (!is.null(seed)) set.seed(seed)
  
  ## Dispatch to model-specific bootstrap
  result <- switch(model_type,
                   "dcc" = dcc_bootstrap_se_internal(
                     residuals = residuals,
                     weights = weights,
                     mle_result = mle_result,
                     n_boot = n_boot,
                     method = method,
                     distribution = distribution,
                     verbose = verbose,
                     ...
                   ),
                   "cgarch" = cgarch_bootstrap_se_internal(
                     residuals = residuals,
                     weights = weights,
                     mle_result = mle_result,
                     n_boot = n_boot,
                     method = method,
                     distribution = distribution,
                     verbose = verbose,
                     ...
                   ),
                   "gogarch" = gogarch_bootstrap_se_internal(
                     residuals = residuals,
                     weights = weights,
                     mle_result = mle_result,
                     n_boot = n_boot,
                     method = method,
                     distribution = distribution,
                     verbose = verbose,
                     ...
                   )
  )
  
  result$model_type <- model_type
  result$boot_method <- method
  result$n_boot <- n_boot
  
  result
}


#### ______________________________________________________________________ ####
#### SECTION 2: DCC BOOTSTRAP (preserved from original)                     ####

#' @title Bootstrap Standard Errors for DCC(1,1) Parameters
#' @description Compute bootstrap standard errors using either parametric or
#'   residual resampling.
#' @param std_resid T x k matrix of standardized residuals
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional correlation matrix
#' @param mle_params MLE estimates c(alpha, beta) or c(alpha, beta, shape)
#' @param n_boot Number of bootstrap replications (default 200)
#' @param method Bootstrap method: "parametric" or "residual"
#' @param distribution "mvn" or "mvt"
#' @param shape Shape parameter for MVT (required if distribution="mvt")
#' @param verbose Print progress
#' @param seed Random seed for reproducibility
#' @return List with bootstrap results
#' @export
dcc11_bootstrap_se <- function(
    std_resid,
    weights,
    Qbar,
    mle_params,
    n_boot = 200,
    method = c("residual", "parametric"),
    distribution = "mvn",
    shape = NULL,
    verbose = TRUE,
    seed = NULL
) {
  
  method <- match.arg(method)
  
  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(std_resid)
  k <- ncol(std_resid)
  
  alpha_mle <- mle_params[1]
  beta_mle <- mle_params[2]
  shape_mle <- if (length(mle_params) > 2) mle_params[3] else shape
  
  n_params <- if (distribution == "mvt") 3 else 2
  param_names <- if (distribution == "mvt") c("alpha", "beta", "shape") else c("alpha", "beta")
  
  ## Storage for bootstrap estimates
  boot_estimates <- matrix(NA, n_boot, n_params)
  colnames(boot_estimates) <- param_names
  
  if (verbose) {
    cat(sprintf("Running %s bootstrap with %d replications...\n", method, n_boot))
    pb <- txtProgressBar(min = 0, max = n_boot, style = 3)
  }
  
  for (b in 1:n_boot) {
    
    if (method == "residual") {
      ## Residual bootstrap: resample rows of standardized residuals
      boot_idx <- sample(1:n, n, replace = TRUE)
      boot_resid <- std_resid[boot_idx, , drop = FALSE]
      boot_weights <- weights[boot_idx]
      boot_Qbar <- cor(boot_resid)
      
    } else {
      ## Parametric bootstrap: simulate from fitted DCC model
      boot_resid <- simulate_dcc_residuals(
        n = n,
        k = k,
        alpha = alpha_mle,
        beta = beta_mle,
        Qbar = Qbar,
        distribution = distribution,
        shape = shape_mle
      )
      boot_weights <- weights
      boot_Qbar <- Qbar
    }
    
    ## Re-estimate DCC parameters
    start_params <- if (distribution == "mvt") {
      c(alpha_mle, beta_mle, shape_mle)
    } else {
      c(alpha_mle, beta_mle)
    }
    
    lower <- if (distribution == "mvt") c(1e-6, 1e-6, 2.5) else c(1e-6, 1e-6)
    upper <- if (distribution == "mvt") c(0.5, 0.999, 50) else c(0.5, 0.999)
    
    opt_result <- tryCatch({
      optim(
        par = start_params,
        fn = dcc11_nll,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper,
        std_resid = boot_resid,
        weights = boot_weights,
        Qbar = boot_Qbar,
        distribution = distribution,
        use_reparam = FALSE
      )
    }, error = function(e) NULL)
    
    if (!is.null(opt_result) && opt_result$convergence == 0) {
      boot_estimates[b, ] <- opt_result$par
    }
    
    if (verbose) setTxtProgressBar(pb, b)
  }
  
  if (verbose) {
    close(pb)
    cat("\n")
  }
  
  ## Process bootstrap results
  process_bootstrap_results(boot_estimates, mle_params, param_names, n_boot, verbose)
}


#' @title Internal DCC bootstrap for unified interface
#' @keywords internal
dcc_bootstrap_se_internal <- function(
    residuals,
    weights,
    mle_result,
    n_boot,
    method,
    distribution,
    verbose,
    ...
) {
  
  args <- list(...)
  Qbar <- args$Qbar %||% cor(residuals)
  
  ## Extract MLE params from result structure
  if (is.list(mle_result) && !is.null(mle_result$params)) {
    mle_params <- mle_result$params
  } else if (is.numeric(mle_result)) {
    mle_params <- mle_result
  } else {
    stop("mle_result must contain 'params' or be a numeric vector")
  }
  
  shape <- if (distribution == "mvt" && length(mle_params) > 2) mle_params[3] else args$shape
  
  dcc11_bootstrap_se(
    std_resid = residuals,
    weights = weights,
    Qbar = Qbar,
    mle_params = mle_params,
    n_boot = n_boot,
    method = method,
    distribution = distribution,
    shape = shape,
    verbose = verbose
  )
}


#### ______________________________________________________________________ ####
#### SECTION 3: CGARCH BOOTSTRAP                                            ####

#' @title Bootstrap Standard Errors for CGARCH Parameters
#' @description Compute bootstrap standard errors for Copula GARCH correlation
#'   parameters using residual resampling.
#' @param z_matrix T x k matrix of copula residuals (PIT-transformed)
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional correlation matrix
#' @param mle_params MLE estimates c(alpha, beta) or c(alpha, beta, shape)
#' @param n_boot Number of bootstrap replications (default 200)
#' @param method Bootstrap method: "residual" or "parametric"
#' @param copula_dist Copula distribution: "mvn" or "mvt"
#' @param verbose Print progress
#' @param seed Random seed for reproducibility
#' @return List with bootstrap results
#' @export
cgarch_bootstrap_se <- function(
    z_matrix,
    weights,
    Qbar,
    mle_params,
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
  beta_mle <- mle_params[2]
  shape_mle <- if (length(mle_params) > 2) mle_params[3] else NULL
  
  n_params <- if (copula_dist == "mvt") 3 else 2
  param_names <- if (copula_dist == "mvt") c("alpha", "beta", "shape") else c("alpha", "beta")
  
  boot_estimates <- matrix(NA, n_boot, n_params)
  colnames(boot_estimates) <- param_names
  
  if (verbose) {
    cat(sprintf("Running CGARCH %s bootstrap with %d replications...\n", method, n_boot))
    pb <- txtProgressBar(min = 0, max = n_boot, style = 3)
  }
  
  for (b in 1:n_boot) {
    
    if (method == "residual") {
      ## Residual bootstrap: resample rows
      boot_idx <- sample(1:n, n, replace = TRUE)
      boot_z <- z_matrix[boot_idx, , drop = FALSE]
      boot_weights <- weights[boot_idx]
      boot_Qbar <- cor(boot_z)
      
    } else {
      ## Parametric bootstrap: simulate from fitted copula
      boot_z <- simulate_copula_residuals(
        n = n,
        k = k,
        alpha = alpha_mle,
        beta = beta_mle,
        Qbar = Qbar,
        copula_dist = copula_dist,
        shape = shape_mle
      )
      boot_weights <- weights
      boot_Qbar <- Qbar
    }
    
    ## Re-estimate CGARCH parameters
    start_params <- if (copula_dist == "mvt") {
      c(alpha_mle, beta_mle, shape_mle %||% 8)
    } else {
      c(alpha_mle, beta_mle)
    }
    
    lower <- if (copula_dist == "mvt") c(1e-6, 1e-6, 2.5) else c(1e-6, 1e-6)
    upper <- if (copula_dist == "mvt") c(0.4, 0.95, 30) else c(0.4, 0.95)
    
    opt_result <- tryCatch({
      optim(
        par = start_params,
        fn = cgarch_nll_for_hessian,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper,
        z_matrix = boot_z,
        weights = boot_weights,
        Qbar = boot_Qbar,
        copula_dist = copula_dist,
        use_reparam = FALSE
      )
    }, error = function(e) NULL)
    
    if (!is.null(opt_result) && opt_result$convergence == 0) {
      boot_estimates[b, ] <- opt_result$par
    }
    
    if (verbose) setTxtProgressBar(pb, b)
  }
  
  if (verbose) {
    close(pb)
    cat("\n")
  }
  
  process_bootstrap_results(boot_estimates, mle_params, param_names, n_boot, verbose)
}


#' @title Internal CGARCH bootstrap for unified interface
#' @keywords internal
cgarch_bootstrap_se_internal <- function(
    residuals,
    weights,
    mle_result,
    n_boot,
    method,
    distribution,
    verbose,
    ...
) {
  
  args <- list(...)
  Qbar <- args$Qbar %||% cor(residuals)
  
  ## Extract MLE params
  if (is.list(mle_result)) {
    alpha <- mle_result$dcc_pars$alpha_1 %||% mle_result$params[1]
    beta <- mle_result$dcc_pars$beta_1 %||% mle_result$params[2]
    shape <- mle_result$dist_pars$shape
    mle_params <- if (!is.null(shape)) c(alpha, beta, shape) else c(alpha, beta)
  } else {
    mle_params <- mle_result
  }
  
  copula_dist <- if (distribution %in% c("mvt", "std")) "mvt" else "mvn"
  
  cgarch_bootstrap_se(
    z_matrix = residuals,
    weights = weights,
    Qbar = Qbar,
    mle_params = mle_params,
    n_boot = n_boot,
    method = method,
    copula_dist = copula_dist,
    verbose = verbose
  )
}


#' @title Simulate Copula Residuals
#' @description Simulate residuals from a DCC-style copula for parametric bootstrap.
#' @keywords internal
simulate_copula_residuals <- function(
    n,
    k,
    alpha,
    beta,
    Qbar,
    copula_dist = "mvn",
    shape = NULL
) {
  ## Same structure as DCC simulation
  simulate_dcc_residuals(
    n = n,
    k = k,
    alpha = alpha,
    beta = beta,
    Qbar = Qbar,
    distribution = if (copula_dist == "mvt") "mvt" else "mvn",
    shape = shape
  )
}


#### ______________________________________________________________________ ####
#### SECTION 4: GOGARCH BOOTSTRAP                                           ####

#' @title Bootstrap Standard Errors for GOGARCH Parameters
#' @description Compute bootstrap standard errors for GOGARCH component GARCH
#'   parameters using residual resampling.
#' @param residuals T x k matrix of observed residuals
#' @param weights T-vector of observation weights
#' @param garch_pars List of GARCH parameters for each component
#' @param ica_info ICA decomposition information
#' @param n_boot Number of bootstrap replications (default 200)
#' @param method Bootstrap method: "residual" (recommended) or "parametric"
#' @param distribution GARCH distribution: "norm" or "std"
#' @param verbose Print progress
#' @param seed Random seed for reproducibility
#' @return List with bootstrap results per component
#' @export
gogarch_bootstrap_se <- function(
    residuals,
    weights,
    garch_pars,
    ica_info,
    n_boot = 200,
    method = c("residual", "parametric"),
    distribution = "norm",
    verbose = TRUE,
    seed = NULL
) {
  
  method <- match.arg(method)
  
  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(residuals)
  k <- ncol(residuals)
  n_components <- length(garch_pars)
  
  if (verbose) {
    cat(sprintf("Running GOGARCH %s bootstrap with %d replications...\n", method, n_boot))
    cat(sprintf("Components: %d\n", n_components))
  }
  
  ## Get independent components
  W <- ica_info$W
  S <- residuals %*% t(W)  ## T x k independent components
  
  ## Bootstrap each component separately (they're independent)
  component_results <- vector("list", n_components)
  
  for (i in 1:n_components) {
    if (verbose) cat(sprintf("\nComponent %d:\n", i))
    
    pars_i <- garch_pars[[i]]
    S_i <- S[, i]
    
    ## Get parameter vector and names
    par_vec <- unlist(pars_i)
    param_names <- names(par_vec)
    n_params <- length(par_vec)
    
    boot_estimates <- matrix(NA, n_boot, n_params)
    colnames(boot_estimates) <- param_names
    
    if (verbose) {
      pb <- txtProgressBar(min = 0, max = n_boot, style = 3)
    }
    
    for (b in 1:n_boot) {
      
      if (method == "residual") {
        ## Resample the component series
        boot_idx <- sample(1:n, n, replace = TRUE)
        boot_S <- S_i[boot_idx]
        boot_weights <- weights[boot_idx]
        
      } else {
        ## Parametric: simulate from fitted GARCH
        boot_S <- simulate_garch_component(
          n = n,
          pars = pars_i,
          distribution = distribution
        )
        boot_weights <- weights
      }
      
      ## Re-estimate GARCH parameters for this component
      opt_result <- tryCatch({
        estimate_component_garch(
          component = boot_S,
          weights = boot_weights,
          start_pars = pars_i,
          distribution = distribution
        )
      }, error = function(e) NULL)
      
      if (!is.null(opt_result) && opt_result$convergence == 0) {
        boot_estimates[b, ] <- opt_result$par
      }
      
      if (verbose) setTxtProgressBar(pb, b)
    }
    
    if (verbose) {
      close(pb)
      cat("\n")
    }
    
    component_results[[i]] <- process_bootstrap_results(
      boot_estimates, par_vec, param_names, n_boot, verbose = FALSE
    )
    component_results[[i]]$component <- i
    
    if (verbose) {
      cat(sprintf("  Valid replications: %d/%d\n", 
                  component_results[[i]]$n_valid, n_boot))
      cat(sprintf("  Bootstrap SEs: %s\n",
                  paste(param_names, "=", round(component_results[[i]]$se, 4), 
                        collapse = ", ")))
    }
  }
  
  list(
    component_results = component_results,
    n_components = n_components,
    method = method,
    n_boot = n_boot,
    distribution = distribution
  )
}


#' @title Internal GOGARCH bootstrap for unified interface
#' @keywords internal
gogarch_bootstrap_se_internal <- function(
    residuals,
    weights,
    mle_result,
    n_boot,
    method,
    distribution,
    verbose,
    ...
) {
  
  ## Extract from mle_result
  if (is.list(mle_result) && !is.null(mle_result$coefficients)) {
    garch_pars <- mle_result$coefficients$garch_pars
    ica_info <- mle_result$coefficients$ica_info
  } else {
    garch_pars <- mle_result$garch_pars
    ica_info <- mle_result$ica_info
  }
  
  dist <- if (distribution %in% c("std", "mvt")) "std" else "norm"
  
  gogarch_bootstrap_se(
    residuals = residuals,
    weights = weights,
    garch_pars = garch_pars,
    ica_info = ica_info,
    n_boot = n_boot,
    method = method,
    distribution = dist,
    verbose = verbose
  )
}


#' @title Simulate GARCH Component
#' @description Simulate a univariate GARCH(1,1) series for parametric bootstrap.
#' @keywords internal
simulate_garch_component <- function(n, pars, distribution = "norm") {
  
  omega <- pars$omega %||% pars[1]
  alpha <- pars$alpha1 %||% pars[2]
  beta <- pars$beta1 %||% pars[3]
  shape <- pars$shape
  
  ## Initialize
  sigma2 <- omega / (1 - alpha - beta)
  x <- numeric(n)
  
  for (t in 1:n) {
    if (t == 1) {
      sigma2_t <- sigma2
    } else {
      sigma2_t <- omega + alpha * x[t-1]^2 + beta * sigma2
    }
    sigma2 <- sigma2_t
    
    if (distribution == "std" && !is.null(shape)) {
      x[t] <- sqrt(sigma2_t) * rt(1, df = shape) * sqrt((shape - 2) / shape)
    } else {
      x[t] <- sqrt(sigma2_t) * rnorm(1)
    }
  }
  
  x
}


#' @title Estimate Component GARCH
#' @description Estimate GARCH(1,1) parameters for a single component.
#' @keywords internal
estimate_component_garch <- function(
    component,
    weights,
    start_pars,
    distribution = "norm"
) {
  
  par_vec <- unlist(start_pars)
  n_params <- length(par_vec)
  
  ## Simple GARCH(1,1) NLL
  garch_nll <- function(pars) {
    omega <- pars[1]
    alpha <- pars[2]
    beta <- pars[3]
    shape <- if (n_params > 3) pars[4] else NULL
    
    if (omega <= 0 || alpha < 0 || beta < 0 || alpha + beta >= 1) {
      return(1e10)
    }
    
    n <- length(component)
    sigma2 <- omega / (1 - alpha - beta)
    nll <- 0
    
    for (t in 1:n) {
      if (t > 1) {
        sigma2 <- omega + alpha * component[t-1]^2 + beta * sigma2
      }
      sigma2 <- max(sigma2, 1e-10)
      
      if (distribution == "std" && !is.null(shape)) {
        ## Student-t log-likelihood
        nll <- nll + weights[t] * (
          0.5 * log(sigma2) + 
            (shape + 1) / 2 * log(1 + component[t]^2 / (sigma2 * (shape - 2)))
        )
      } else {
        ## Normal log-likelihood
        nll <- nll + weights[t] * (0.5 * log(sigma2) + 0.5 * component[t]^2 / sigma2)
      }
    }
    
    nll
  }
  
  lower <- if (n_params > 3) c(1e-8, 1e-8, 1e-8, 2.5) else c(1e-8, 1e-8, 1e-8)
  upper <- if (n_params > 3) c(1, 0.5, 0.999, 50) else c(1, 0.5, 0.999)
  
  optim(
    par = par_vec,
    fn = garch_nll,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper
  )
}


#### ______________________________________________________________________ ####
#### SECTION 5: SHARED UTILITIES                                            ####

#' @title Simulate DCC Standardized Residuals
#' @description Simulate standardized residuals from a DCC(1,1) model.
#' @keywords internal
simulate_dcc_residuals <- function(
    n,
    k,
    alpha,
    beta,
    Qbar,
    distribution = "mvn",
    shape = NULL
) {
  
  z <- matrix(0, n, k)
  Q <- Qbar
  
  for (t in 1:n) {
    ## Normalize Q to correlation R
    d <- sqrt(diag(Q))
    d[d < 1e-8] <- 1e-8
    D_inv <- diag(1 / d, k)
    R <- D_inv %*% Q %*% D_inv
    
    ## Ensure PD
    R <- (R + t(R)) / 2
    eig <- eigen(R, symmetric = TRUE)
    if (any(eig$values < 1e-8)) {
      eig$values[eig$values < 1e-8] <- 1e-8
      R <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
      d <- sqrt(diag(R))
      D_inv <- diag(1 / d, k)
      R <- D_inv %*% R %*% D_inv
    }
    
    ## Draw from multivariate distribution
    if (distribution == "mvn") {
      L <- chol(R)
      eps <- rnorm(k)
      z[t, ] <- as.vector(t(L) %*% eps)
    } else {
      ## MVT
      if (is.null(shape)) shape <- 8
      L <- chol(R)
      eps <- rnorm(k)
      chi_sq <- rchisq(1, df = shape)
      z[t, ] <- as.vector(t(L) %*% eps) * sqrt(shape / chi_sq)
    }
    
    ## Update Q for next period
    if (t < n) {
      z_t <- z[t, , drop = FALSE]
      Q <- (1 - alpha - beta) * Qbar + alpha * (t(z_t) %*% z_t) + beta * Q
    }
  }
  
  z
}


#' @title Process Bootstrap Results
#' @description Compute bootstrap statistics from replications.
#' @keywords internal
process_bootstrap_results <- function(
    boot_estimates,
    mle_params,
    param_names,
    n_boot,
    verbose = TRUE
) {
  
  ## Remove failed replications
  valid_idx <- complete.cases(boot_estimates)
  n_valid <- sum(valid_idx)
  boot_estimates_valid <- boot_estimates[valid_idx, , drop = FALSE]
  
  if (n_valid < n_boot * 0.5) {
    warning(sprintf("Only %d/%d bootstrap replications succeeded", n_valid, n_boot))
  }
  
  if (n_valid < 3) {
    return(list(
      se = rep(NA, length(mle_params)),
      boot_estimates = boot_estimates,
      n_valid = n_valid,
      mean = rep(NA, length(mle_params)),
      bias = rep(NA, length(mle_params)),
      ci_percentile = NULL,
      ci_basic = NULL,
      ci_bca = NULL,
      valid = FALSE,
      reason = "insufficient_valid_replications"
    ))
  }
  
  ## Compute bootstrap statistics
  boot_mean <- colMeans(boot_estimates_valid)
  boot_se <- apply(boot_estimates_valid, 2, sd)
  boot_bias <- boot_mean - mle_params[1:length(boot_mean)]
  
  names(boot_se) <- param_names
  names(boot_mean) <- param_names
  names(boot_bias) <- param_names
  
  ## Percentile confidence intervals (95%)
  ci_percentile <- apply(boot_estimates_valid, 2, quantile, probs = c(0.025, 0.975))
  colnames(ci_percentile) <- param_names
  
  ## Basic bootstrap CI
  ci_basic <- matrix(NA, 2, length(mle_params))
  rownames(ci_basic) <- c("2.5%", "97.5%")
  colnames(ci_basic) <- param_names
  ci_basic[1, ] <- 2 * mle_params[1:ncol(ci_basic)] - ci_percentile[2, ]
  ci_basic[2, ] <- 2 * mle_params[1:ncol(ci_basic)] - ci_percentile[1, ]
  
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
    ci_basic = ci_basic,
    ci_bca = NULL,  ## BCa computed separately if needed
    valid = TRUE,
    reason = "ok"
  )
}


#### ______________________________________________________________________ ####
#### SECTION 6: PROFILE LIKELIHOOD (DCC and CGARCH)                         ####

#' @title Profile Likelihood Confidence Intervals
#' @description Compute profile likelihood CIs for correlation parameters.
#'   Supports DCC and CGARCH models.
#' @param model_type "dcc" or "cgarch"
#' @param residuals T x k matrix of residuals
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional correlation matrix
#' @param mle_params MLE estimates
#' @param mle_nll NLL at MLE
#' @param conf_level Confidence level (default 0.95)
#' @param param Which parameter: "alpha", "beta", or "both"
#' @param n_points Number of profile points
#' @param distribution Distribution type
#' @param verbose Print progress
#' @return List with profile likelihood results
#' @export
correlation_profile_likelihood_ci <- function(
    model_type = c("dcc", "cgarch"),
    residuals,
    weights,
    Qbar,
    mle_params,
    mle_nll = NULL,
    conf_level = 0.95,
    param = c("both", "alpha", "beta"),
    n_points = 50,
    distribution = "mvn",
    verbose = TRUE
) {
  
  model_type <- match.arg(model_type)
  param <- match.arg(param)
  
  ## Select NLL function based on model type
  nll_fn <- switch(model_type,
                   "dcc" = function(params, ...) {
                     dcc11_nll(params = params, std_resid = residuals, weights = weights,
                               Qbar = Qbar, distribution = distribution, use_reparam = FALSE)
                   },
                   "cgarch" = function(params, ...) {
                     copula_dist <- if (distribution %in% c("mvt", "std")) "mvt" else "mvn"
                     cgarch_nll_for_hessian(params = params, z_matrix = residuals, 
                                            weights = weights, Qbar = Qbar,
                                            copula_dist = copula_dist, use_reparam = FALSE)
                   }
  )
  
  alpha_mle <- mle_params[1]
  beta_mle <- mle_params[2]
  
  ## Compute MLE NLL if not provided
  if (is.null(mle_nll)) {
    mle_nll <- nll_fn(mle_params)
  }
  
  ## Critical value for LR test
  chi_sq_crit <- qchisq(conf_level, df = 1)
  nll_threshold <- mle_nll + chi_sq_crit / 2
  
  results <- list()
  
  ## Profile alpha
  if (param %in% c("both", "alpha")) {
    if (verbose) cat("Profiling alpha...\n")
    
    results$alpha <- profile_one_param_generic(
      param_idx = 1,
      param_name = "alpha",
      mle_params = mle_params,
      mle_nll = mle_nll,
      nll_threshold = nll_threshold,
      nll_fn = nll_fn,
      n_points = n_points,
      model_type = model_type
    )
  }
  
  ## Profile beta
  if (param %in% c("both", "beta")) {
    if (verbose) cat("Profiling beta...\n")
    
    results$beta <- profile_one_param_generic(
      param_idx = 2,
      param_name = "beta",
      mle_params = mle_params,
      mle_nll = mle_nll,
      nll_threshold = nll_threshold,
      nll_fn = nll_fn,
      n_points = n_points,
      model_type = model_type
    )
  }
  
  results$mle_params <- mle_params
  results$mle_nll <- mle_nll
  results$conf_level <- conf_level
  results$nll_threshold <- nll_threshold
  results$model_type <- model_type
  
  if (verbose) {
    cat("\nProfile Likelihood CIs:\n")
    if (!is.null(results$alpha)) {
      cat(sprintf("  alpha: [%.4f, %.4f]\n", 
                  results$alpha$ci[1], results$alpha$ci[2]))
    }
    if (!is.null(results$beta)) {
      cat(sprintf("  beta:  [%.4f, %.4f]\n", 
                  results$beta$ci[1], results$beta$ci[2]))
    }
  }
  
  results
}


#' @title Profile One Parameter (Generic)
#' @description Internal function to compute profile likelihood for one parameter.
#' @keywords internal
profile_one_param_generic <- function(
    param_idx,
    param_name,
    mle_params,
    mle_nll,
    nll_threshold,
    nll_fn,
    n_points,
    model_type
) {
  
  mle_value <- mle_params[param_idx]
  other_idx <- if (param_idx == 1) 2 else 1
  
  ## Determine search range
  if (param_name == "alpha") {
    search_min <- 1e-6
    search_max <- min(0.4, 0.98 - mle_params[2])
  } else {
    search_min <- max(0.01, 1e-6)
    search_max <- min(0.98, 0.98 - mle_params[1])
  }
  
  ## Create parameter grid
  grid_lower <- seq(search_min, mle_value, length.out = n_points / 2)
  grid_upper <- seq(mle_value, search_max, length.out = n_points / 2)
  param_grid <- unique(sort(c(grid_lower, grid_upper)))
  
  ## Compute profile likelihood
  profile_nll <- numeric(length(param_grid))
  profile_other <- numeric(length(param_grid))
  
  for (i in seq_along(param_grid)) {
    fixed_value <- param_grid[i]
    
    ## Optimize over other parameter
    opt_result <- tryCatch({
      optimize(
        f = function(other_val) {
          params <- if (param_idx == 1) c(fixed_value, other_val) else c(other_val, fixed_value)
          if (sum(params[1:2]) >= 0.999) return(1e10)
          nll_fn(params)
        },
        interval = c(1e-6, 0.98 - fixed_value),
        maximum = FALSE
      )
    }, error = function(e) list(objective = NA, minimum = NA))
    
    profile_nll[i] <- opt_result$objective
    profile_other[i] <- opt_result$minimum
  }
  
  ## Find CI bounds
  valid_idx <- is.finite(profile_nll)
  if (sum(valid_idx) < 3) {
    return(list(ci = c(NA, NA), grid = param_grid, nll = profile_nll))
  }
  
  ## Interpolate to find crossing points
  ci_lower <- NA
  ci_upper <- NA
  
  ## Lower bound
  lower_idx <- which(param_grid < mle_value & valid_idx)
  if (length(lower_idx) > 0) {
    for (i in rev(lower_idx)) {
      if (profile_nll[i] > nll_threshold) {
        if (i < length(param_grid) && profile_nll[i+1] <= nll_threshold) {
          ## Linear interpolation
          ci_lower <- param_grid[i] + (param_grid[i+1] - param_grid[i]) *
            (nll_threshold - profile_nll[i]) / (profile_nll[i+1] - profile_nll[i])
        } else {
          ci_lower <- param_grid[i]
        }
        break
      }
    }
    if (is.na(ci_lower)) ci_lower <- search_min
  }
  
  ## Upper bound
  upper_idx <- which(param_grid > mle_value & valid_idx)
  if (length(upper_idx) > 0) {
    for (i in upper_idx) {
      if (profile_nll[i] > nll_threshold) {
        if (i > 1 && profile_nll[i-1] <= nll_threshold) {
          ci_upper <- param_grid[i-1] + (param_grid[i] - param_grid[i-1]) *
            (nll_threshold - profile_nll[i-1]) / (profile_nll[i] - profile_nll[i-1])
        } else {
          ci_upper <- param_grid[i]
        }
        break
      }
    }
    if (is.na(ci_upper)) ci_upper <- search_max
  }
  
  list(
    ci = c(ci_lower, ci_upper),
    mle = mle_value,
    grid = param_grid,
    nll = profile_nll,
    other_params = profile_other,
    nll_threshold = nll_threshold
  )
}


#### ______________________________________________________________________ ####
#### SECTION 7: COMPREHENSIVE INFERENCE                                     ####


#' @title Comprehensive Correlation Model Inference
#' @description Compute all three types of inference (Hessian, Bootstrap, Profile)
#'   for DCC or CGARCH correlation parameters.
#' @param model_type "dcc" or "cgarch"
#' @param residuals T x k matrix of residuals
#' @param weights T-vector of observation weights
#' @param Qbar k x k unconditional correlation matrix (optional)
#' @param mle_params MLE estimates (optional, computed if NULL)
#' @param distribution "mvn" or "mvt"
#' @param n_boot Number of bootstrap replications
#' @param boot_method Bootstrap method
#' @param n_profile_points Points for profile likelihood
#' @param conf_level Confidence level
#' @param verbose Print progress
#' @param seed Random seed
#' @return List with all inference results
#' @export
correlation_comprehensive_inference <- function(
    model_type = c("dcc", "cgarch", "adcc"),
    residuals,
    weights,
    Qbar = NULL,
    mle_params = NULL,
    distribution = "mvn",
    n_boot = 200,
    boot_method = "residual",
    n_profile_points = 50,
    conf_level = 0.95,
    verbose = TRUE,
    seed = NULL
) {
  
  model_type <- match.arg(model_type)
  
  if (model_type == "adcc") {
    return(adcc_comprehensive_inference(...))
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(residuals)
  if (is.null(Qbar)) Qbar <- cor(residuals)
  
  ## Select NLL function
  nll_fn <- switch(model_type,
                   "dcc" = function(params) {
                     dcc11_nll(params, residuals, weights, Qbar, distribution, FALSE)
                   },
                   "cgarch" = function(params) {
                     copula_dist <- if (distribution %in% c("mvt", "std")) "mvt" else "mvn"
                     cgarch_nll_for_hessian(params, residuals, weights, Qbar, copula_dist, FALSE)
                   }
  )
  
  ## Step 1: Find MLE if not provided
  if (is.null(mle_params)) {
    if (verbose) cat("Finding MLE...\n")
    
    start <- c(0.05, 0.90)
    lower <- c(1e-6, 1e-6)
    upper <- c(0.4, 0.95)
    
    if (distribution %in% c("mvt", "std")) {
      start <- c(start, 8)
      lower <- c(lower, 2.5)
      upper <- c(upper, 30)
    }
    
    opt_result <- optim(
      par = start,
      fn = nll_fn,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper
    )
    mle_params <- opt_result$par
    mle_nll <- opt_result$value
  } else {
    mle_nll <- nll_fn(mle_params)
  }
  
  param_names <- if (length(mle_params) > 2) c("alpha", "beta", "shape") else c("alpha", "beta")
  names(mle_params) <- param_names
  
  if (verbose) {
    cat(sprintf("MLE: %s\n", paste(param_names, "=", round(mle_params, 4), collapse = ", ")))
    cat("\n")
  }
  
  ## Step 2: Hessian-based SEs
  if (verbose) cat("Computing Hessian-based SEs...\n")
  
  hessian_result <- tryCatch({
    if (model_type == "dcc") {
      dcc11_standard_errors(mle_params, residuals, weights, Qbar, 
                            distribution, use_reparam = FALSE)
    } else {
      copula_dist <- if (distribution %in% c("mvt", "std")) "mvt" else "mvn"
      cgarch_standard_errors(mle_params, residuals, weights, Qbar,
                             copula_dist, use_reparam = FALSE)
    }
  }, error = function(e) list(se = rep(NA, length(mle_params)), warning = e$message))
  
  z_crit <- qnorm(1 - (1 - conf_level) / 2)
  hessian_ci <- matrix(NA, 2, length(mle_params))
  rownames(hessian_ci) <- c("lower", "upper")
  colnames(hessian_ci) <- param_names
  
  if (all(is.finite(hessian_result$se))) {
    hessian_ci[1, ] <- mle_params - z_crit * hessian_result$se
    hessian_ci[2, ] <- mle_params + z_crit * hessian_result$se
  }
  
  ## Step 3: Bootstrap SEs
  if (verbose) cat("\nComputing bootstrap SEs...\n")
  
  boot_result <- if (model_type == "dcc") {
    dcc11_bootstrap_se(residuals, weights, Qbar, mle_params, n_boot,
                       boot_method, distribution, verbose = verbose, seed = seed)
  } else {
    copula_dist <- if (distribution %in% c("mvt", "std")) "mvt" else "mvn"
    cgarch_bootstrap_se(residuals, weights, Qbar, mle_params, n_boot,
                        boot_method, copula_dist, verbose = verbose, seed = seed)
  }
  
  ## Step 4: Profile likelihood CIs (only for alpha and beta)
  if (verbose) cat("\nComputing profile likelihood CIs...\n")
  
  profile_result <- correlation_profile_likelihood_ci(
    model_type = model_type,
    residuals = residuals,
    weights = weights,
    Qbar = Qbar,
    mle_params = mle_params[1:2],  ## Only alpha and beta
    mle_nll = mle_nll,
    conf_level = conf_level,
    param = "both",
    n_points = n_profile_points,
    distribution = distribution,
    verbose = verbose
  )
  
  ## Compile summary
  summary_table <- data.frame(
    Parameter = param_names,
    MLE = mle_params,
    Hessian_SE = hessian_result$se,
    Hessian_CI_lower = hessian_ci[1, ],
    Hessian_CI_upper = hessian_ci[2, ],
    Boot_SE = boot_result$se[param_names],
    Boot_CI_lower = boot_result$ci_percentile[1, param_names],
    Boot_CI_upper = boot_result$ci_percentile[2, param_names],
    row.names = NULL
  )
  
  ## Add profile CIs for alpha and beta only
  summary_table$Profile_CI_lower <- c(
    profile_result$alpha$ci[1],
    profile_result$beta$ci[1],
    if (length(mle_params) > 2) NA else NULL
  )
  summary_table$Profile_CI_upper <- c(
    profile_result$alpha$ci[2],
    profile_result$beta$ci[2],
    if (length(mle_params) > 2) NA else NULL
  )
  
  if (verbose) {
    cat("\n")
    cat("=== Inference Summary ===\n")
    cat(sprintf("Model type: %s\n", toupper(model_type)))
    cat(sprintf("Confidence level: %.0f%%\n\n", conf_level * 100))
    print(summary_table, digits = 4)
  }
  
  list(
    mle = mle_params,
    mle_nll = mle_nll,
    hessian = list(
      se = hessian_result$se,
      ci = hessian_ci,
      eigenvalues = hessian_result$eigenvalues
    ),
    bootstrap = boot_result,
    profile = profile_result,
    summary = summary_table,
    conf_level = conf_level,
    model_type = model_type,
    settings = list(
      n = n,
      n_boot = n_boot,
      boot_method = boot_method,
      n_profile_points = n_profile_points,
      distribution = distribution
    )
  )
}


#' @title GOGARCH Comprehensive Inference
#' @description Compute bootstrap-based inference for GOGARCH model.
#'   Profile likelihood is not available for GOGARCH due to the ICA structure.
#' @param residuals T x k matrix of observed residuals
#' @param weights T-vector of observation weights
#' @param garch_pars List of GARCH parameters per component
#' @param ica_info ICA decomposition information
#' @param n_boot Number of bootstrap replications
#' @param boot_method Bootstrap method
#' @param distribution "norm" or "std"
#' @param conf_level Confidence level
#' @param verbose Print progress
#' @param seed Random seed
#' @return List with inference results
#' @export
gogarch_comprehensive_inference <- function(
    residuals,
    weights,
    garch_pars,
    ica_info,
    n_boot = 200,
    boot_method = "residual",
    distribution = "norm",
    conf_level = 0.95,
    verbose = TRUE,
    seed = NULL
) {
  
  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(residuals)
  k <- ncol(residuals)
  n_components <- length(garch_pars)
  
  if (verbose) {
    cat("=== GOGARCH Comprehensive Inference ===\n")
    cat(sprintf("Observations: %d, Series: %d, Components: %d\n", n, k, n_components))
    cat("\n")
  }
  
  ## Step 1: Hessian-based SEs (per component)
  if (verbose) cat("Computing Hessian-based SEs...\n")
  
  hessian_result <- gogarch_standard_errors(
    garch_pars = garch_pars,
    ica_info = ica_info,
    residuals = residuals,
    weights = weights,
    distribution = distribution,
    method = "hessian"
  )
  
  ## Step 2: Bootstrap SEs
  if (verbose) cat("\nComputing bootstrap SEs...\n")
  
  boot_result <- gogarch_bootstrap_se(
    residuals = residuals,
    weights = weights,
    garch_pars = garch_pars,
    ica_info = ica_info,
    n_boot = n_boot,
    method = boot_method,
    distribution = distribution,
    verbose = verbose,
    seed = seed
  )
  
  ## Compile component summaries
  z_crit <- qnorm(1 - (1 - conf_level) / 2)
  component_summaries <- vector("list", n_components)
  
  for (i in 1:n_components) {
    pars <- unlist(garch_pars[[i]])
    param_names <- names(pars)
    
    hess_se <- hessian_result$component_se[[i]]$se
    boot_se <- boot_result$component_results[[i]]$se
    
    ## CIs
    hess_ci_lower <- pars - z_crit * hess_se
    hess_ci_upper <- pars + z_crit * hess_se
    boot_ci <- boot_result$component_results[[i]]$ci_percentile
    
    component_summaries[[i]] <- data.frame(
      Parameter = param_names,
      Estimate = pars,
      Hessian_SE = hess_se,
      Hessian_CI_lower = hess_ci_lower,
      Hessian_CI_upper = hess_ci_upper,
      Boot_SE = boot_se,
      Boot_CI_lower = boot_ci[1, ],
      Boot_CI_upper = boot_ci[2, ],
      row.names = NULL
    )
  }
  
  if (verbose) {
    cat("\n=== Inference Summary ===\n")
    cat(sprintf("Confidence level: %.0f%%\n\n", conf_level * 100))
    
    for (i in 1:n_components) {
      cat(sprintf("--- Component %d ---\n", i))
      print(component_summaries[[i]], digits = 4)
      cat("\n")
    }
  }
  
  list(
    component_summaries = component_summaries,
    hessian = hessian_result,
    bootstrap = boot_result,
    conf_level = conf_level,
    n_components = n_components,
    settings = list(
      n = n,
      n_boot = n_boot,
      boot_method = boot_method,
      distribution = distribution
    )
  )
}


#' @title Print Inference Comparison
#' @description Print a formatted comparison of inference methods.
#' @param inf_result Result from correlation_comprehensive_inference or 
#'   gogarch_comprehensive_inference
#' @export
print_inference_comparison <- function(inf_result) {
  
  if (!is.null(inf_result$n_components)) {
    ## GOGARCH result
    print_gogarch_inference(inf_result)
  } else {
    ## DCC/CGARCH result
    print_correlation_inference(inf_result)
  }
  
  invisible(inf_result)
}


#' @keywords internal
print_correlation_inference <- function(inf_result) {
  
  model_name <- toupper(inf_result$model_type %||% "DCC")
  
  cat("\n")
  cat("================================================================\n")
  cat(sprintf("  %s Inference Comparison\n", model_name))
  cat("================================================================\n\n")
  
  cat(sprintf("Sample size: %d\n", inf_result$settings$n))
  cat(sprintf("MLE: %s\n", 
              paste(names(inf_result$mle), "=", round(inf_result$mle, 4), collapse = ", ")))
  cat(sprintf("Confidence level: %.0f%%\n\n", inf_result$conf_level * 100))
  
  ## Standard Errors comparison
  cat("Standard Errors:\n")
  param_names <- names(inf_result$mle)
  header <- sprintf("  %-10s", "Method")
  for (p in param_names) header <- paste0(header, sprintf(" %12s", p))
  cat(header, "\n")
  
  hess_line <- sprintf("  %-10s", "Hessian")
  boot_line <- sprintf("  %-10s", "Bootstrap")
  ratio_line <- sprintf("  %-10s", "Ratio")
  
  for (i in seq_along(param_names)) {
    hess_line <- paste0(hess_line, sprintf(" %12.4f", inf_result$hessian$se[i]))
    boot_line <- paste0(boot_line, sprintf(" %12.4f", inf_result$bootstrap$se[i]))
    ratio <- inf_result$hessian$se[i] / inf_result$bootstrap$se[i]
    ratio_line <- paste0(ratio_line, sprintf(" %12.2f", ratio))
  }
  cat(hess_line, "\n")
  cat(boot_line, "\n")
  cat(ratio_line, "\n\n")
  
  ## Confidence Intervals
  cat("Confidence Intervals:\n")
  cat(sprintf("  %-10s", "Method"))
  for (p in param_names[1:2]) cat(sprintf(" %20s", p))
  cat("\n")
  
  cat(sprintf("  %-10s [%8.4f, %8.4f] [%8.4f, %8.4f]\n", "Hessian",
              inf_result$hessian$ci[1, 1], inf_result$hessian$ci[2, 1],
              inf_result$hessian$ci[1, 2], inf_result$hessian$ci[2, 2]))
  cat(sprintf("  %-10s [%8.4f, %8.4f] [%8.4f, %8.4f]\n", "Bootstrap",
              inf_result$bootstrap$ci_percentile[1, 1], 
              inf_result$bootstrap$ci_percentile[2, 1],
              inf_result$bootstrap$ci_percentile[1, 2], 
              inf_result$bootstrap$ci_percentile[2, 2]))
  
  if (!is.null(inf_result$profile$alpha$ci)) {
    cat(sprintf("  %-10s [%8.4f, %8.4f] [%8.4f, %8.4f]\n", "Profile",
                inf_result$profile$alpha$ci[1], inf_result$profile$alpha$ci[2],
                inf_result$profile$beta$ci[1], inf_result$profile$beta$ci[2]))
  }
  cat("\n")
  
  ## Warning about Hessian reliability
  hess_boot_ratio <- inf_result$hessian$se[1:2] / inf_result$bootstrap$se[1:2]
  if (any(hess_boot_ratio < 0.7, na.rm = TRUE)) {
    cat("NOTE: Hessian SEs appear underestimated (ratio < 0.7).\n")
    cat("      Bootstrap or Profile CIs recommended for inference.\n")
  }
}


#' @keywords internal
print_gogarch_inference <- function(inf_result) {
  
  cat("\n")
  cat("================================================================\n")
  cat("  GOGARCH Inference Comparison\n")
  cat("================================================================\n\n")
  
  cat(sprintf("Sample size: %d\n", inf_result$settings$n))
  cat(sprintf("Components: %d\n", inf_result$n_components))
  cat(sprintf("Confidence level: %.0f%%\n\n", inf_result$conf_level * 100))
  
  for (i in seq_along(inf_result$component_summaries)) {
    cat(sprintf("--- Component %d ---\n", i))
    
    summ <- inf_result$component_summaries[[i]]
    
    cat("  SE Comparison:\n")
    cat(sprintf("    %-10s %12s %12s %8s\n", "Parameter", "Hessian", "Bootstrap", "Ratio"))
    
    for (j in 1:nrow(summ)) {
      ratio <- summ$Hessian_SE[j] / summ$Boot_SE[j]
      cat(sprintf("    %-10s %12.4f %12.4f %8.2f\n",
                  summ$Parameter[j], summ$Hessian_SE[j], summ$Boot_SE[j], ratio))
    }
    cat("\n")
  }
  
  cat("NOTE: For GOGARCH, bootstrap SEs are generally recommended\n")
  cat("      due to ICA estimation uncertainty.\n")
}


