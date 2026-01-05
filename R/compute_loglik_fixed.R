#' Compute Log-Likelihood with Fixed Parameters
#'
#' @description
#' Evaluates the log-likelihood for tsmarch models using fixed (non-estimated)
#' parameters. This function uses the existing tsmarch estimation object and
#' replaces the estimated parameters with user-provided fixed values.
#'
#' @param object An estimated tsmarch object of class "dcc.estimate", 
#'   "cgarch.estimate", or "gogarch.estimate"
#' @param params A named list of fixed parameters. The names should match the
#'   parameter names in the model's parmatrix:
#'   \itemize{
#'     \item For DCC models: `list(alpha_1 = 0.05, beta_1 = 0.90)` for DCC(1,1),
#'       or `list(alpha_1 = 0.05, alpha_2 = 0.03, beta_1 = 0.85)` for DCC(2,1).
#'       For ADCC models, add `gamma_1`, `gamma_2`, etc.
#'     \item For Student-t DCC: Add `shape = 3` to specify degrees of freedom.
#'     \item For Copula-GARCH models: Same as DCC models above. For Student-t copula,
#'       use `shape = 3` for degrees of freedom.
#'     \item For GOGARCH models: Parameters for each independent component GARCH model,
#'       e.g., `list(omega_1 = 0.01, alpha_1 = 0.05, beta_1 = 0.90, omega_2 = 0.02, ...)`.
#'   }
#'   Parameter names use underscore notation where the number indicates the lag order
#'   or component index. Use `coef(object)` to see the exact parameter names for your model.
#' @param return_components Logical. If TRUE, returns both the total log-likelihood
#'   and its components (univariate GARCH + multivariate). Default is FALSE.
#' @param ll_vec Logical. If TRUE, returns per-observation log-likelihood vector
#'   instead of the total. Cannot be used with return_components = TRUE. Default is FALSE.
#' @param ... Additional arguments (currently unused)
#'
#' @return
#'   If ll_vec = TRUE, returns a numeric vector of per-observation log-likelihoods.
#'   The vector has length n-1 for DCC models where n is the number of observations,
#'   because the first observation serves as initialization for the DCC recursion
#'   (tsmarch returns a zero placeholder for this observation which is removed).
#'   The sum of the returned vector equals the scalar log-likelihood returned when
#'   ll_vec = FALSE, and also equals logLik(object) when using estimated parameters.
#'
#'   If ll_vec = FALSE and return_components = FALSE (default), returns a single
#'   numeric value representing the total log-likelihood.
#'
#'   If ll_vec = FALSE and return_components = TRUE, returns a list with components:
#'
#'   \itemize{
#'     \item{loglik}{Total log-likelihood}
#'     \item{garch_loglik}{Univariate GARCH component (DCC/Copula only)}
#'     \item{multivariate_loglik}{Multivariate component log-likelihood}
#'  }
#' 
#' @details
#' This function extracts the specification from an estimated tsmarch object,
#' replaces the estimated parameters with the user-provided fixed parameters,
#' and computes the log-likelihood without re-estimating.
#'
#' The function works by:
#' \enumerate{
#'   \item Extracting the model specification from the estimated object
#'   \item Updating the parmatrix with the fixed parameter values
#'   \item Calling the internal tsmarch likelihood computation functions
#' }
#'
#' @export
#' @examples
#' \donttest{
#' # This example requires tsmarch and tsgarch packages
#' if (require(tsmarch) && require(tsgarch) && require(xts)) {
#'   # Generate small sample data
#'   set.seed(123)
#'   n <- 500
#'   returns <- matrix(rnorm(n * 2), ncol = 2)
#'   returns <- xts(returns, order.by = seq.Date(Sys.Date() - n + 1, 
#'                                                 Sys.Date(), by = "day"))
#'   colnames(returns) <- c("series1", "series2")
#'   
#'   # Estimate univariate GARCH models
#'   spec1 <- garch_modelspec(returns[,1], model = "garch", order = c(1,1))
#'   spec2 <- garch_modelspec(returns[,2], model = "garch", order = c(1,1))
#'   fit1 <- estimate(spec1, keep_tmb = TRUE)
#'   fit2 <- estimate(spec2, keep_tmb = TRUE)
#'   
#'   # Combine into multivariate
#'   garch_fits <- to_multi_estimate(list(fit1, fit2))
#'   
#'   # Estimate DCC model
#'   dcc_spec <- dcc_modelspec(garch_fits, dynamics = "dcc", dcc_order = c(1,1))
#'   dcc_fit <- estimate(dcc_spec)
#'   
#'   # Get estimated parameters
#'   est_params <- coef(dcc_fit)
#'   
#'   # Compute log-likelihood at estimated parameters
#'   ll_at_est <- compute_loglik_fixed(dcc_fit, params = as.list(est_params))
#'   
#'   # Compute at alternative parameters
#'   ll_alt <- compute_loglik_fixed(dcc_fit, 
#'                                   params = list(alpha_1 = 0.05, beta_1 = 0.90))
#'   
#'   # The estimated parameters should give higher likelihood
#'   print(paste("LL at estimated:", ll_at_est))
#'   print(paste("LL at alternative:", ll_alt))
#'   print(paste("Difference:", ll_at_est - ll_alt))
#' }
#' }
#' 
#' @examples
#' \dontrun{
#' # Extended example with profile likelihood and LR tests
#' library(tsmarch)
#' library(tsgarch)
#' library(xts)
#' 
#' # Generate sample data with correlation structure
#' set.seed(100)
#' n <- 1500
#' # Create correlated innovations
#' rho <- 0.6
#' Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
#' L <- chol(Sigma)
#' z <- matrix(rnorm(n * 2), ncol = 2) %*% L
#' 
#' # Add GARCH dynamics
#' returns <- z
#' for (i in 2:n) {
#'   h <- 0.01 + 0.08 * returns[i-1,]^2 + 0.90 * returns[i-1,]^2
#'   returns[i,] <- returns[i,] * sqrt(pmax(h, 0.001))
#' }
#' 
#' returns <- xts(returns, order.by = seq(Sys.Date() - n + 1, Sys.Date(), by = "day"))
#' colnames(returns) <- c("asset1", "asset2")
#' 
#' # Estimate univariate GARCH models
#' spec1 <- garch_modelspec(returns[,1], model = "garch", order = c(1,1))
#' spec2 <- garch_modelspec(returns[,2], model = "garch", order = c(1,1))
#' fit1 <- estimate(spec1, keep_tmb = TRUE)
#' fit2 <- estimate(spec2, keep_tmb = TRUE)
#' 
#' # Combine into multivariate
#' garch_fits <- to_multi_estimate(list(fit1, fit2))
#' 
#' # Estimate DCC model
#' dcc_spec <- dcc_modelspec(garch_fits, dynamics = "dcc", dcc_order = c(1,1))
#' dcc_fit <- estimate(dcc_spec)
#' 
#' # Get estimated parameters
#' est_params <- coef(dcc_fit)
#' cat("Estimated DCC parameters:\n")
#' print(est_params)
#' 
#' # Compute log-likelihood at estimated parameters
#' ll_at_estimated <- compute_loglik_fixed(dcc_fit, params = as.list(est_params))
#' cat("\nLL at estimated parameters:", ll_at_estimated, "\n")
#' 
#' # Test with alternative parameter values
#' ll_alternative <- compute_loglik_fixed(
#'   dcc_fit,
#'   params = list(alpha_1 = 0.03, beta_1 = 0.95)
#' )
#' cat("LL at alternative parameters:", ll_alternative, "\n")
#' 
#' # Likelihood ratio test
#' lr_stat <- 2 * (ll_at_estimated - ll_alternative)
#' p_value <- pchisq(lr_stat, df = 2, lower.tail = FALSE)
#' cat("\nLikelihood Ratio Test:\n")
#' cat("  LR statistic:", round(lr_stat, 4), "\n")
#' cat("  P-value:", format.pval(p_value, digits = 4), "\n")
#' 
#' # Profile likelihood for alpha
#' # Only compute if estimated alpha is not at boundary
#' alpha_est <- est_params["alpha_1"]
#' beta_est <- est_params["beta_1"]
#' 
#' if (alpha_est > 0.01 && alpha_est < 0.2) {
#'   cat("\nComputing profile likelihood for alpha...\n")
#'   
#'   # Create grid around estimated alpha
#'   alpha_range <- seq(max(0.001, alpha_est - 0.03), 
#'                      min(0.3, alpha_est + 0.03), 
#'                      length.out = 20)
#'   
#'   profile_ll <- sapply(alpha_range, function(a) {
#'     # Skip if would violate stationarity
#'     if (a + beta_est >= 0.999) return(NA_real_)
#'     
#'     tryCatch({
#'       compute_loglik_fixed(dcc_fit, 
#'                            params = list(alpha_1 = a, beta_1 = beta_est))
#'     }, error = function(e) NA_real_)
#'   })
#'   
#'   # Plot profile likelihood
#'   valid_idx <- !is.na(profile_ll)
#'   if (sum(valid_idx) > 5) {
#'     plot(alpha_range[valid_idx], profile_ll[valid_idx], type = "b", 
#'          xlab = expression(alpha), ylab = "Log-Likelihood",
#'          main = "Profile Likelihood for Alpha Parameter",
#'          pch = 19, col = "blue")
#'     abline(v = alpha_est, col = "red", lty = 2, lwd = 2)
#'     
#'     # Add confidence interval based on chi-squared cutoff
#'     ll_cutoff <- max(profile_ll, na.rm = TRUE) - qchisq(0.95, 1)/2
#'     abline(h = ll_cutoff, col = "gray", lty = 3)
#'     
#'     legend("bottomright", 
#'            legend = c("Profile LL", "Estimated alpha", "95% CI cutoff"), 
#'            col = c("blue", "red", "gray"), 
#'            lty = c(1, 2, 3), pch = c(19, NA, NA), cex = 0.8)
#'   }
#' } else {
#'   cat("\nSkipping profile likelihood (alpha at boundary)\n")
#'   cat("For a better example, try a different seed or larger sample.\n")
#' }
#' 
#' # Get component-wise log-likelihoods
#' ll_components <- compute_loglik_fixed(
#'   dcc_fit,
#'   params = as.list(est_params),
#'   return_components = TRUE
#' )
#' cat("\nComponent-wise log-likelihoods:\n")
#' cat("  Total:", round(ll_components$loglik, 2), "\n")
#' cat("  GARCH:", round(ll_components$garch_loglik, 2), "\n")
#' cat("  DCC:  ", round(ll_components$multivariate_loglik, 2), "\n")
#' cat("  Sum:  ", round(ll_components$garch_loglik + 
#'                        ll_components$multivariate_loglik, 2), "\n")
#' }
compute_loglik_fixed <- function(
    object, 
    params, 
    return_components = FALSE, 
    ll_vec = FALSE,
    ...) {
  
  ## Validate input object
  if (!inherits(object, c("dcc.estimate", "cgarch.estimate", "gogarch.estimate"))) {
    stop("object must be an estimated tsmarch model of class 'dcc.estimate', 'cgarch.estimate', or 'gogarch.estimate'")
  }
  
  if (!is.list(params)) {
    stop("params must be a named list")
  }
  
  ## Cannot use both return_components and ll_vec
  if (return_components && ll_vec) {
    stop("Cannot use both return_components = TRUE and ll_vec = TRUE")
  }

  ## Dispatch to appropriate method
  if (inherits(object, "dcc.estimate")) {
    if (inherits(object, "dcc.dynamic")) {
      result <- .compute_dcc_dynamic_loglik(object, params, return_components, ll_vec)
    } else if (inherits(object, "dcc.constant")) {
      result <- .compute_dcc_constant_loglik(object, params, return_components, ll_vec)
    } else {
      stop("Unknown DCC model type")
    }
  } else if (inherits(object, "cgarch.estimate")) {
    if (inherits(object, "cgarch.dynamic")) {
      result <- .compute_copula_dynamic_loglik(object, params, return_components, ll_vec)
    } else if (inherits(object, "cgarch.constant")) {
      result <- .compute_copula_constant_loglik(object, params, return_components, ll_vec)
    } else {
      stop("Unknown Copula-GARCH model type")
    }
  } else if (inherits(object, "gogarch.estimate")) {
    result <- .compute_gogarch_loglik(object, params, return_components, ll_vec)
  }
  
  return(result)
}

## Internal method for DCC dynamic
.compute_dcc_dynamic_loglik <- function(
    object, 
    params, 
    return_components = FALSE,
    ll_vec = FALSE
  ) {
  
  ## Create a copy of the spec
  spec <- object$spec
  
  ## Update parameters in the parmatrix
  #spec$parmatrix <- .update_parmatrix(spec$parmatrix, params)
  ## If no params provided, use the fitted values from the object
  if (length(params) == 0) {
    spec$parmatrix <- object$parmatrix
  } else {
    spec$parmatrix <- .update_parmatrix(object$parmatrix, params)
  }
  
  ## Get parameter values for the multivariate component
  estimate <- NULL
  pars <- spec$parmatrix[estimate == 1]$value
  
  if (ll_vec) {
    ## Get vector of per-observation negative log-likelihoods
    ## NOTE: For DCC models, ll_vec includes BOTH GARCH and DCC components
    total_nll_vec <- tsmarch:::.dcc_dynamic_values(pars, spec, type = "ll_vec")

    # ## tsmarch:::.dcc_dynamic_values returns n+1 values:
    # ##   - Row 1: initialization placeholder (always 0)
    # ##   - Rows 2:(n+1): actual log-likelihoods for observations 1:n
    # ## We remove the first row plus any additional DCC burn-in
    # dccorder <- spec$dynamics$order
    # maxpq <- max(dccorder)
    # n_remove <- 1 + maxpq  ## 1 for placeholder + maxpq for DCC burn-in
    # total_nll_vec <- total_nll_vec[-(1:n_remove), , drop = TRUE]
    
    ## *** FIXED ***
    ## Remove only the initialization placeholder (first row, always 0)
    ## Note: tsmarch's type="nll" includes all subsequent observations, so we do too
    total_nll_vec <- total_nll_vec[-1, , drop = TRUE]
    
    ## Return positive per-observation log-likelihoods
    return(-total_nll_vec)
  }
  
  # Compute total negative log-likelihood using internal tsmarch function
  # Note: For DCC models, type = "nll" returns the TOTAL (GARCH + DCC) NLL
  total_nll <- tsmarch:::.dcc_dynamic_values(pars, spec, type = "nll")
  
  if (return_components) {
    # Get DCC component only
    dcc_nll <- tsmarch:::.dcc_dynamic_values(pars, spec, type = "dcc_nll")
    # Calculate GARCH component as difference
    garch_nll <- total_nll - dcc_nll
    
    return(list(
      loglik = -total_nll,        # Convert to positive LL
      garch_loglik = -garch_nll,  # Convert to positive LL
      multivariate_loglik = -dcc_nll  # Convert to positive LL
    ))
  } else {
    # Return positive log-likelihood (to match logLik() method)
    return(-total_nll)
  }
}

## Internal method for DCC constant
.compute_dcc_constant_loglik <- function(
    object, 
    params, 
    return_components = FALSE,
    ll_vec = FALSE
  ) {
  
  ## For constant correlation with mvn, there are no parameters to estimate
  if (object$spec$distribution == "mvn" && length(params) == 0) {
    if (ll_vec) {
      ## Get vector of per-observation likelihoods
      ## NOTE: For DCC models, ll_vec includes BOTH GARCH and DCC components
      total_nll_vec <- tsmarch:::.dcc_constant_values(NULL, object$spec, type = "ll_vec")
      return(-as.vector(total_nll_vec))
    }
    
    total_nll <- object$loglik
    
    if (return_components) {
      ## Get DCC component
      dcc_nll <- tsmarch:::.dcc_constant_values(NULL, object$spec, type = "dcc_nll")
      ## Calculate GARCH component as difference
      garch_nll <- total_nll - dcc_nll
      
      return(list(
        loglik = -total_nll,
        garch_loglik = -garch_nll,
        multivariate_loglik = -dcc_nll
      ))
    } else {
      return(-total_nll)
    }
  }
  
  ## Create a copy of the spec
  spec <- object$spec
  
  ## Update parameters in the parmatrix
  spec$parmatrix <- .update_parmatrix(spec$parmatrix, params)
  
  ## Get parameter values
  estimate <- NULL
  pars <- spec$parmatrix[estimate == 1]$value
  
  if (ll_vec) {
    ## Get vector of per-observation negative log-likelihoods
    ## NOTE: For DCC models, ll_vec includes BOTH GARCH and DCC components
    total_nll_vec <- tsmarch:::.dcc_constant_values(pars, spec, type = "ll_vec")
    return(-as.vector(total_nll_vec))
  }
  
  ## Compute total negative log-likelihood
  ## Note: For DCC models, type = "nll" returns the TOTAL (GARCH + DCC) NLL
  total_nll <- tsmarch:::.dcc_constant_values(pars, spec, type = "nll")
  
  if (return_components) {
    ## Get DCC component only  
    dcc_nll <- tsmarch:::.dcc_constant_values(pars, spec, type = "dcc_nll")
    ## Calculate GARCH component as difference
    garch_nll <- total_nll - dcc_nll
    
    return(list(
      loglik = -total_nll,
      garch_loglik = -garch_nll,
      multivariate_loglik = -dcc_nll
    ))
  } else {
    return(-total_nll)
  }
}

.compute_copula_dynamic_loglik <- function(
    object, 
    params, 
    return_components = FALSE,
    ll_vec = FALSE
  ) {

  ## Create a copy of the spec
  spec <- object$spec

  ## Update parameters in the parmatrix
  spec$parmatrix <- .update_parmatrix(spec$parmatrix, params)

  ## Get parameter values
  estimate <- NULL
  pars <- spec$parmatrix[estimate == 1]$value

  if (ll_vec) {
    ## Get vector of per-observation negative log-likelihoods
    ## NOTE: For copula models, ll_vec contains ONLY the copula component
    ## We need to add the GARCH component separately
    copula_nll_vec <- tsmarch:::.copula_dynamic_values(pars, spec, type = "ll_vec")
    
    ## Strip off the first max(p, q) observations (initialization period)
    dccorder <- spec$dynamics$order
    maxpq <- max(dccorder)
    if (maxpq > 0) {
      copula_nll_vec <- copula_nll_vec[-(1:maxpq), , drop = TRUE]
    } else {
      ## Convert to vector (it's a matrix with 1 column)
      copula_nll_vec <- as.vector(copula_nll_vec)
    }
    
    ## Get per-observation univariate GARCH negative log-likelihoods
    garch_nll_vec <- .get_garch_nll_vec(object$spec$univariate)
    
    ## Total per-observation negative log-likelihood
    total_nll_vec <- garch_nll_vec + copula_nll_vec
    
    ## Return positive per-observation log-likelihoods
    return(-total_nll_vec)
  }
  
  ## Compute multivariate negative log-likelihood (copula component)
  mv_nll <- tsmarch:::.copula_dynamic_values(pars, spec, type = "nll")

  ## Get univariate GARCH negative log-likelihoods
  garch_nll <- sum(sapply(object$spec$univariate, function(x) x$loglik))

  ## Total negative log-likelihood
  total_nll <- garch_nll + mv_nll

  if (return_components) {
    return(list(
      loglik = -total_nll,
      garch_loglik = -garch_nll,
      multivariate_loglik = -mv_nll
    ))
  } else {
    return(-total_nll)
  }
}

## Internal method for Copula constant
.compute_copula_constant_loglik <- function(
    object, 
    params, 
    return_components = FALSE,
    ll_vec = FALSE
  ) {
  
  ## For constant correlation with mvn, there are no parameters to estimate
  if (object$spec$copula == "mvn" && length(params) == 0) {
    if (ll_vec) {
      ## NOTE: For copula models, ll_vec contains ONLY the copula component
      copula_nll_vec <- tsmarch:::.copula_constant_values(NULL, object$spec, type = "ll_vec")
      copula_nll_vec <- as.vector(copula_nll_vec)
      
      ## Get per-observation univariate GARCH negative log-likelihoods
      garch_nll_vec <- .get_garch_nll_vec(object$spec$univariate)
      
      ## Total per-observation negative log-likelihood
      total_nll_vec <- garch_nll_vec + copula_nll_vec
      
      return(-total_nll_vec)
    }
    
    if (return_components) {
      garch_ll <- sum(sapply(object$spec$univariate, function(x) as.numeric(logLik(x))))
      mv_nll <- tsmarch:::.copula_constant_values(NULL, object$spec, type = "nll")
      mv_ll <- -mv_nll
      return(list(
        loglik = garch_ll + mv_ll,
        garch_loglik = garch_ll,
        multivariate_loglik = mv_ll
      ))
    } else {
      return(-object$loglik)
    }
  }
  
  ## Create a copy of the spec
  spec <- object$spec
  
  ## Update parameters in the parmatrix
  spec$parmatrix <- .update_parmatrix(spec$parmatrix, params)
  
  ## Get parameter values
  estimate <- NULL
  pars <- spec$parmatrix[estimate == 1]$value
  
  ## Compute multivariate log-likelihood
  mv_nll <- tsmarch:::.copula_constant_values(pars, spec, type = "nll")
  
  if (return_components) {
    garch_nll <- sum(sapply(object$spec$univariate, function(x) as.numeric(logLik(x))))
    total_nll <- garch_nll + mv_nll
    
    return(list(
      loglik = -total_nll,
      garch_loglik = -garch_nll,
      multivariate_loglik = -mv_nll
    ))
  } else {
    garch_nll <- sum(sapply(object$spec$univariate, function(x) as.numeric(logLik(x))))
    total_nll <- garch_nll + mv_nll
    return(-total_nll)
  }
}

# Internal method for GOGARCH
.compute_gogarch_loglik <- function(object, params = list(), 
                                    return_components = FALSE, 
                                    ll_vec = FALSE) {
  
  n_components <- length(object$univariate)
  K <- object$ica$K
  
  # Initialize storage for component likelihoods
  component_lls <- vector("list", n_components)
  component_ll_vecs <- vector("list", n_components)
  
  # Compute likelihood for each independent component
  for (i in 1:n_components) {
    comp <- object$univariate[[i]]
    tmb_env <- comp$tmb$env
    
    # Determine which parameters to use
    if (length(params) == 0) {
      # Use estimated parameters from last.par.best
      comp_params <- get("last.par.best", envir = tmb_env)
    } else {
      # Extract parameters for this component from params list
      param_names <- c(paste0("omega_", i), 
                       paste0("alpha_", i), 
                       paste0("beta_", i))
      
      # Get last.par.best as default
      last_par_best <- get("last.par.best", envir = tmb_env)
      comp_params <- last_par_best
      names(comp_params) <- c("omega", "alpha", "beta")
      
      # Override with provided parameters if available
      for (j in 1:3) {
        if (param_names[j] %in% names(params)) {
          comp_params[j] <- params[[param_names[j]]]
        }
      }
    }
    
    # Evaluate TMB function and get report (returns negative log-likelihood)
    nll <- comp$tmb$fn(comp_params)
    rep <- comp$tmb$report(comp_params)
    
    # Store results
    component_lls[[i]] <- -nll  # Convert to log-likelihood
    
    # ll_vector from report is the likelihood contribution, need to take log
    component_ll_vecs[[i]] <- log(rep$ll_vector)
  }
  
  # Compute Jacobian adjustment for transformation
  # log|det(K)| accounts for the change of variables from independent components
  # to the original observed variables
  if (is_square(K)) {
    jacobian_adj <- log(abs(det(K)))
  } else {
    jacobian_adj <- log(abs(det(K %*% t(K))))
  }
  
  if (ll_vec) {
    # Return observation-wise log-likelihoods
    n_obs <- length(component_ll_vecs[[1]])
    
    # Create matrix of log-likelihood vectors
    log_lik_matrix <- matrix(0, nrow = n_obs, ncol = n_components)
    for (i in 1:n_components) {
      log_lik_matrix[, i] <- component_ll_vecs[[i]]
    }
    
    # Sum log-likelihoods across components for each observation
    # Add Jacobian adjustment (constant, distributed equally across observations)
    obs_ll <- rowSums(log_lik_matrix) + jacobian_adj / n_obs
    
    if (return_components) {
      return(list(
        loglik = sum(obs_ll),  # Total log-likelihood
        lik_vector = obs_ll,
        component_logliks = unlist(component_lls),
        jacobian_adjustment = jacobian_adj
      ))
    } else {
      return(obs_ll)
    }
    
  } else {
    # Return scalar log-likelihood
    total_ll <- sum(unlist(component_lls)) + jacobian_adj
    
    if (return_components) {
      return(list(
        loglik = total_ll,
        component_logliks = unlist(component_lls),
        jacobian_adjustment = jacobian_adj
      ))
    } else {
      return(total_ll)
    }
  }
}

# Helper function to check if matrix is square
is_square <- function(K) {
  nrow(K) == ncol(K)
}



# Helper function to parse GOGARCH parameters
.parse_gogarch_params <- function(params, object) {
  # Parse parameter list into component-specific vectors
  # Expected format: list(omega_1 = ..., alpha_1 = ..., beta_1 = ..., 
  #                       omega_2 = ..., alpha_2 = ..., beta_2 = ...)
  
  n_components <- length(object$univariate)
  params_by_comp <- vector("list", n_components)
  
  # Get parameter names for each component
  for (i in 1:n_components) {
    comp_parmatrix <- object$univariate[[i]]$parmatrix
    estimate <- NULL
    param_names <- comp_parmatrix[estimate == 1]$parameter
    
    # Extract values for this component
    # Try with suffix _i first (e.g., omega_1, alpha_1)
    comp_params <- numeric(length(param_names))
    names(comp_params) <- param_names
    
    for (j in seq_along(param_names)) {
      param_base <- param_names[j]
      # Try with component suffix
      param_with_suffix <- paste0(param_base, "_", i)
      
      if (param_with_suffix %in% names(params)) {
        comp_params[j] <- params[[param_with_suffix]]
      } else if (param_base %in% names(params)) {
        # Fallback: use base name (assumes same value for all components)
        comp_params[j] <- params[[param_base]]
      } else {
        # Use estimated value if not provided
        comp_params[j] <- comp_parmatrix[estimate == 1]$value[j]
      }
    }
    
    params_by_comp[[i]] <- comp_params
  }
  
  return(params_by_comp)
}

## Helper function to get per-observation GARCH negative log-likelihoods
.get_garch_nll_vec <- function(univariate_list) {
  ## Extract per-observation NLL from each univariate GARCH model
  ## Sum across series to get total GARCH contribution per time point

  n_obs <- length(univariate_list[[1]]$spec$target$y_orig)
  n_series <- length(univariate_list)

  ## Initialize matrix: rows = observations, cols = series
  garch_nll_matrix <- matrix(0, nrow = n_obs, ncol = n_series)

  for (i in seq_along(univariate_list)) {
    ## Extract the lik_vector from the univariate model
    ## This contains per-observation NEGATIVE log-likelihoods
    if (!is.null(univariate_list[[i]]$lik_vector)) {
      garch_nll_matrix[, i] <- univariate_list[[i]]$lik_vector
    } else {
      stop("Per-observation likelihoods (lik_vector) not found in univariate model ", i,
           ". The model may need to be re-estimated.")
    }
  }

  ## Sum across series to get per-observation total GARCH NLL contribution
  return(rowSums(garch_nll_matrix))
}

## Helper function to update parmatrix with new parameter values
.update_parmatrix <- function(parmatrix, params) {
  
  parameter <- NULL
  
  ## Make a copy to avoid modifying the original
  pmatrix <- data.table::copy(parmatrix)
  
  ## Update each parameter
  for (param_name in names(params)) {
    ## Find matching parameter(s) in parmatrix
    matches <- which(pmatrix$parameter == param_name)
    
    if (length(matches) == 0) {
      warning("Parameter '", param_name, "' not found in model parmatrix. Ignoring.")
    } else if (length(matches) > 1) {
      ## Multiple parameters with same name (shouldn't happen in DCC/Copula)
      warning("Multiple parameters named '", param_name, "' found. Updating all.")
      pmatrix[matches, ]$value <- params[[param_name]]
    } else {
      pmatrix[matches, ]$value <- params[[param_name]]
    }
  }
  
  return(pmatrix)
}