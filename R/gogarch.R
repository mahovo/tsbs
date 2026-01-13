## ============================================================================
## tsbs_gogarch.R
## GOGARCH-specific functions for MS-VARMA-GARCH framework
## ============================================================================

#' @title Weighted GARCH Estimation for GOGARCH Models
#' 
#' @description
#' Implements weighted maximum likelihood estimation for GOGARCH models in the
#' context of Markov-Switching frameworks. GOGARCH differs from DCC/CGARCH in
#' that the correlation structure comes from the ICA mixing matrix, not from
#' estimated correlation parameters.
#' 
#' The estimation proceeds as follows:
#' 1. Perform ICA decomposition on the residuals to obtain independent components
#' 2. Estimate univariate GARCH models on each ICA component using weighted MLE
#' 
#' Unlike DCC/CGARCH, there are no correlation dynamics parameters to estimate.
#' The time-varying correlation comes from the time-varying component volatilities
#' transformed through the fixed ICA mixing matrix.
#'
#' @param residuals Matrix of residuals (T x k)
#' @param weights Vector of state probabilities/weights (length T)
#' @param spec Model specification list containing garch_spec_args and start_pars
#' @param diagnostics Optional diagnostics object for logging
#' @param verbose Logical; print progress information
#' 
#' @return A list containing:
#'   \item{coefficients}{List with garch_pars (one per ICA component) and ica_info}
#'   \item{warnings}{List of any warnings generated during estimation}
#'   \item{diagnostics}{Updated diagnostics object}
#'   
#' @keywords internal
estimate_garch_weighted_gogarch <- function(
    residuals, 
    weights, 
    spec, 
    diagnostics = NULL, 
    verbose = FALSE
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
  
  ## === STEP 1: Perform ICA Decomposition ===
  ## ICA extracts independent components from the residuals
  ## The mixing matrix A and unmixing matrix W are key for GOGARCH
  
  ## Convert to xts for tsmarch compatibility
  residuals_xts <- xts::xts(residuals, order.by = seq.Date(
    from = as.Date("2000-01-01"), 
    by = "day", 
    length.out = T_obs
  ))
  colnames(residuals_xts) <- paste0("series_", 1:k)
  
  ## Get ICA parameters from spec
  ica_method <- spec$garch_spec_args$ica %||% "radical"
  n_components <- spec$garch_spec_args$components %||% k
  
  ## Perform ICA decomposition
  ## Use tsmarch's radical or fastica implementation
  ic <- tryCatch({
    if (ica_method == "radical") {
      tsmarch::radical(residuals, components = n_components, demean = FALSE, trace = verbose)
    } else if (ica_method == "fastica") {
      tsmarch::fastica(residuals, components = n_components, demean = FALSE, trace = verbose)
    } else {
      stop("Unsupported ICA method: ", ica_method)
    }
  }, error = function(e) {
    warning("ICA decomposition failed: ", e$message, ". Using identity transformation.")
    ## Fallback: identity transformation (no rotation)
    list(
      S = residuals,
      A = diag(k),
      W = diag(k),
      K = diag(k)
    )
  })
  
  ## Extract independent components
  S <- ic$S  # T x n_components matrix of independent components
  
  if (verbose) {
    cat(sprintf("GOGARCH: ICA decomposition complete. %d components extracted.\n", ncol(S)))
  }
  
  ## === STEP 2: Estimate Univariate GARCH on Each ICA Component (Weighted) ===
  
  garch_pars_list <- list()
  warnings_list <- list()
  
  garch_model <- spec$garch_spec_args$model %||% "garch"
  garch_order <- spec$garch_spec_args$order %||% c(1, 1)
  distribution <- spec$distribution %||% "norm"
  
  for (i in 1:n_components) {
    component_data <- S[, i]
    
    ## Create spec for this component
    comp_spec <- list(
      garch_model = garch_model,
      garch_order = garch_order,
      distribution = distribution,
      start_pars = list(
        garch_pars = spec$start_pars$garch_pars[[i]] %||% 
          list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
        dist_pars = spec$start_pars$dist_pars
      )
    )
    
    ## Estimate using weighted univariate function
    uni_result <- estimate_garch_weighted_univariate_gogarch(
      residuals = component_data,
      weights = w_target,
      spec = comp_spec,
      verbose = verbose
    )
    
    garch_pars_list[[i]] <- uni_result$coefficients
    warnings_list <- c(warnings_list, uni_result$warnings)
  }
  
  ## === STEP 3: Check for Boundary Conditions ===
  omega_boundary_threshold <- 1e-8
  
  for (i in 1:n_components) {
    omega_est <- garch_pars_list[[i]]$omega
    
    if (!is.null(omega_est) && omega_est < omega_boundary_threshold) {
      warning(sprintf(
        "GOGARCH Component %d: omega near boundary (%.2e < %.2e). Volatility dynamics may be degenerate.",
        i, omega_est, omega_boundary_threshold
      ))
    }
  }
  
  ## === Return Results ===
  return(list(
    coefficients = list(
      garch_pars = garch_pars_list,
      ica_info = list(
        A = ic$A,        # Mixing matrix
        W = ic$W,        # Unmixing matrix  
        K = ic$K,        # Pre-whitening matrix
        S = ic$S,        # Independent components (for reference)
        method = ica_method,
        n_components = n_components
      ),
      dist_pars = spec$start_pars$dist_pars,
      correlation_type = "gogarch"  # Always GOGARCH structure
    ),
    warnings = warnings_list,
    diagnostics = diagnostics
  ))
}


#' @title Weighted Univariate GARCH Estimation for GOGARCH Components
#' 
#' @description
#' Estimates univariate GARCH parameters for a single ICA component using
#' weighted maximum likelihood. This is a simplified version of
#' estimate_garch_weighted_univariate that works with ICA component data.
#'
#' @param residuals Vector of ICA component values
#' @param weights Vector of weights (state probabilities)
#' @param spec Component specification with garch_model, garch_order, distribution
#' @param verbose Logical; print progress
#' 
#' @return List with coefficients and warnings
#' @keywords internal
estimate_garch_weighted_univariate_gogarch <- function(
    residuals,
    weights,
    spec,
    verbose = FALSE
) {
  
  ## Extract parameters
  start_pars <- c(spec$start_pars$garch_pars, spec$start_pars$dist_pars)
  
  if (length(start_pars) == 0) {
    return(list(coefficients = list(), warnings = list()))
  }
  
  ## Create univariate GARCH spec using tsgarch
  garch_spec <- tsgarch::garch_modelspec(
    y = xts::xts(residuals, order.by = seq.Date(
      from = as.Date("2000-01-01"),
      by = "day", 
      length.out = length(residuals)
    )),
    model = spec$garch_model %||% "garch",
    order = spec$garch_order %||% c(1, 1),
    constant = FALSE,  # ICA components are demeaned
    distribution = spec$distribution %||% "norm"
  )
  
  ## Extract bounds
  parmatrix <- garch_spec$parmatrix
  pars_to_estimate <- names(start_pars)
  
  ## Match parameters - handle naming differences
  param_mapping <- list(
    omega = "omega",
    alpha1 = "alpha1", alpha2 = "alpha2",
    beta1 = "beta1", beta2 = "beta2",
    skew = "skew", shape = "shape", lambda = "lambda"
  )
  
  lower_bounds <- numeric(length(pars_to_estimate))
  upper_bounds <- numeric(length(pars_to_estimate))
  
  for (j in seq_along(pars_to_estimate)) {
    pname <- pars_to_estimate[j]
    if (pname %in% parmatrix$parameter) {
      row_idx <- which(parmatrix$parameter == pname)
      lower_bounds[j] <- parmatrix$lower[row_idx]
      upper_bounds[j] <- parmatrix$upper[row_idx]
    } else {
      ## Default bounds
      lower_bounds[j] <- 1e-12
      upper_bounds[j] <- 1 - 1e-6
    }
  }
  
  ## Weighted log-likelihood objective
  weighted_loglik <- function(params, resid, w, dist) {
    omega <- params["omega"]
    alpha <- params[grepl("alpha", names(params))]
    beta <- params[grepl("beta", names(params))]
    
    ## Simple GARCH(p,q) variance recursion
    n <- length(resid)
    sigma2 <- rep(var(resid), n)  # Initialize with unconditional variance
    
    p <- length(alpha)
    q <- length(beta)
    maxpq <- max(p, q)
    
    for (t in (maxpq + 1):n) {
      sigma2[t] <- omega
      for (i in 1:p) {
        sigma2[t] <- sigma2[t] + alpha[i] * resid[t - i]^2
      }
      for (j in 1:q) {
        sigma2[t] <- sigma2[t] + beta[j] * sigma2[t - j]
      }
      ## Ensure positive variance
      if (sigma2[t] <= 0) sigma2[t] <- 1e-10
    }
    
    sig <- sqrt(sigma2)
    
    ## Compute log-likelihood based on distribution
    if (dist == "norm") {
      ll <- dnorm(resid, mean = 0, sd = sig, log = TRUE)
    } else if (dist == "std") {
      shape <- params["shape"]
      if (is.na(shape) || shape <= 2) shape <- 4
      ll <- tsdistributions::dstd(resid, mu = 0, sigma = sig, shape = shape, log = TRUE)
    } else if (dist == "sstd") {
      shape <- params["shape"]
      skew <- params["skew"]
      if (is.na(shape) || shape <= 2) shape <- 4
      if (is.na(skew)) skew <- 1
      ll <- tsdistributions::dsstd(resid, mu = 0, sigma = sig, shape = shape, skew = skew, log = TRUE)
    } else if (dist == "nig") {
      shape <- params["shape"]
      skew <- params["skew"]
      if (is.na(shape) || shape <= 0) shape <- 1
      if (is.na(skew)) skew <- 0
      ll <- tsdistributions::dnig(resid, mu = 0, sigma = sig, shape = shape, skew = skew, log = TRUE)
    } else if (dist == "gh") {
      shape <- params["shape"]
      skew <- params["skew"]
      lambda <- params["lambda"]
      if (is.na(shape) || shape <= 0) shape <- 1
      if (is.na(skew)) skew <- 0
      if (is.na(lambda)) lambda <- -0.5
      ll <- tsdistributions::dgh(resid, mu = 0, sigma = sig, shape = shape, skew = skew, 
                                  lambda = lambda, log = TRUE)
    } else {
      ## Fallback to normal
      ll <- dnorm(resid, mean = 0, sd = sig, log = TRUE)
    }
    
    ll[!is.finite(ll)] <- -1e10
    
    ## Return weighted negative log-likelihood
    return(-sum(w * ll, na.rm = TRUE))
  }
  
  ## Optimize
  warnings_list <- list()
  
  opt_result <- tryCatch({
    withCallingHandlers({
      stats::optim(
        par = unlist(start_pars),
        fn = weighted_loglik,
        lower = lower_bounds,
        upper = upper_bounds,
        method = "L-BFGS-B",
        control = list(ndeps = rep(1e-8, length(start_pars))),
        resid = residuals,
        w = weights,
        dist = spec$distribution %||% "norm"
      )
    }, warning = function(w) {
      warnings_list <<- c(warnings_list, list(w))
      invokeRestart("muffleWarning")
    })
  }, error = function(e) {
    warning("GARCH optimization failed: ", e$message)
    list(par = unlist(start_pars), convergence = 1)
  })
  
  estimated_coeffs <- as.list(opt_result$par)
  names(estimated_coeffs) <- names(start_pars)
  
  return(list(coefficients = estimated_coeffs, warnings = warnings_list))
}


#' @title Compute GOGARCH Log-Likelihood for MS Framework
#' 
#' @description
#' Computes the GOGARCH log-likelihood given estimated parameters and ICA info.
#' Used in the E-step of the EM algorithm.
#'
#' @param residuals Matrix of residuals
#' @param garch_pars List of GARCH parameters for each component
#' @param ica_info List with ICA matrices (A, W, K)
#' @param distribution Distribution name
#' @param return_vector If TRUE, return vector of per-observation log-likelihoods
#' 
#' @return Scalar log-likelihood or vector of per-observation values
#' @keywords internal
compute_gogarch_loglik_ms <- function(
    residuals,
    garch_pars,
    ica_info,
    distribution = "norm",
    return_vector = FALSE
) {
  
  T_obs <- nrow(residuals)
  n_components <- length(garch_pars)
  
  ## Transform residuals to independent components using W
  S <- residuals %*% t(ica_info$W)
  
  ## Compute component-wise log-likelihoods
  ll_matrix <- matrix(0, nrow = T_obs, ncol = n_components)
  
  for (i in 1:n_components) {
    pars <- garch_pars[[i]]
    component <- S[, i]
    
    ## Compute GARCH variance path
    omega <- pars$omega %||% 0.1
    alpha <- unlist(pars[grepl("alpha", names(pars))])
    beta <- unlist(pars[grepl("beta", names(pars))])
    
    if (length(alpha) == 0) alpha <- 0.1
    if (length(beta) == 0) beta <- 0.8
    
    sigma2 <- rep(var(component), T_obs)
    p <- length(alpha)
    q <- length(beta)
    maxpq <- max(p, q, 1)
    
    for (t in (maxpq + 1):T_obs) {
      sigma2[t] <- omega
      for (j in 1:p) {
        sigma2[t] <- sigma2[t] + alpha[j] * component[t - j]^2
      }
      for (j in 1:q) {
        sigma2[t] <- sigma2[t] + beta[j] * sigma2[t - j]
      }
      if (sigma2[t] <= 0) sigma2[t] <- 1e-10
    }
    
    sig <- sqrt(sigma2)
    
    ## Compute log-likelihood
    if (distribution == "norm") {
      ll_matrix[, i] <- dnorm(component, mean = 0, sd = sig, log = TRUE)
    } else {
      ## For other distributions, extract distribution parameters
      shape <- pars$shape %||% 4
      skew <- pars$skew %||% 1
      
      ll_matrix[, i] <- tryCatch({
        switch(distribution,
          "std" = tsdistributions::dstd(component, mu = 0, sigma = sig, shape = shape, log = TRUE),
          "nig" = tsdistributions::dnig(component, mu = 0, sigma = sig, shape = shape, skew = skew, log = TRUE),
          "gh" = tsdistributions::dgh(component, mu = 0, sigma = sig, shape = shape, skew = skew, 
                                       lambda = pars$lambda %||% -0.5, log = TRUE),
          dnorm(component, mean = 0, sd = sig, log = TRUE)
        )
      }, error = function(e) {
        dnorm(component, mean = 0, sd = sig, log = TRUE)
      })
    }
  }
  
  ## Sum across components for each observation
  ll_vector <- rowSums(ll_matrix)
  
  ## Add Jacobian adjustment: log|det(K)|
  K <- ica_info$K
  if (nrow(K) == ncol(K)) {
    jacobian_adj <- log(abs(det(K)))
  } else {
    jacobian_adj <- log(abs(det(K %*% t(K))))
  }
  
  if (return_vector) {
    ## Distribute Jacobian across observations
    return(ll_vector + jacobian_adj / T_obs)
  } else {
    return(sum(ll_vector) + jacobian_adj)
  }
}
