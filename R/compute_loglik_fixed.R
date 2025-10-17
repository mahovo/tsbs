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
#'   
#'   If ll_vec = FALSE and return_components = FALSE (default), returns a single 
#'   numeric value representing the total log-likelihood. 
#'   
#'   If ll_vec = FALSE and return_components = TRUE, returns a list with components:
#'   \item{loglik}{Total log-likelihood}
#'   \item{garch_loglik}{Univariate GARCH component (DCC/Copula only)}
#'   \item{multivariate_loglik}{Multivariate component log-likelihood}
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
    ll_vec = FALSE,
    return_components = FALSE
  ) {
  
  ## Create a copy of the spec
  spec <- object$spec
  
  ## Update parameters in the parmatrix
  spec$parmatrix <- .update_parmatrix(spec$parmatrix, params)
  
  ## Get parameter values for the multivariate component
  estimate <- NULL
  pars <- spec$parmatrix[estimate == 1]$value
  
  if (ll_vec) {
    ## Get vector of per-observation negative log-likelihoods
    mv_nll_vec <- tsmarch:::.dcc_dynamic_values(pars, spec, type = "ll_vec")
    
    ## Get per-observation univariate GARCH negative log-likelihoods
    garch_nll_vec <- .get_garch_nll_vec(object$spec$univariate)
    
    ## Total per-observation negative log-likelihood
    total_nll_vec <- garch_nll_vec + mv_nll_vec
    
    ## Return positive per-observation log-likelihoods
    return(-total_nll_vec)
  }
  
  ## Compute multivariate log-likelihood using internal tsmarch function
  ## Note: .dcc_dynamic_values returns NEGATIVE log-likelihood (NLL)
  mv_nll <- tsmarch:::.dcc_dynamic_values(pars, spec, type = "nll")
  
  ## Get univariate GARCH negative log-likelihoods
  garch_nll <- sum(sapply(object$spec$univariate, function(x) x$loglik))
  
  ## Total negative log-likelihood (what tsmarch stores in object$loglik)
  total_nll <- garch_nll + mv_nll
  
  if (return_components) {
    return(list(
      loglik = -total_nll,          ## Convert to positive ll
      garch_loglik = -garch_nll,    ## Convert to positive ll
      multivariate_loglik = -mv_nll ## Convert to positive ll
    ))
  } else {
    ## Return positive log-likelihood (to match logLik() method)
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
      mv_nll_vec <- tsmarch:::.dcc_constant_values(NULL, object$spec, type = "ll_vec")
      garch_nll_vec <- .get_garch_nll_vec(object$spec$univariate)
      total_nll_vec <- garch_nll_vec + mv_nll_vec
      return(-total_nll_vec)
    }
    
    total_nll <- object$loglik
    
    if (return_components) {
      garch_nll <- sum(sapply(object$spec$univariate, function(x) as.numeric(logLik(x))))
      mv_nll <- total_nll - garch_nll
      return(list(
        loglik = -total_nll,
        garch_loglik = -garch_nll,
        multivariate_loglik = -mv_nll
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
    # Get per-observation negative log-likelihoods
    mv_nll_vec <- tsmarch:::.dcc_constant_values(pars, spec, type = "ll_vec")
    garch_nll_vec <- .get_garch_nll_vec(object$spec$univariate)
    total_nll_vec <- garch_nll_vec + mv_nll_vec
    return(-total_nll_vec)
  }
  
  ## Compute multivariate log-likelihood
  mv_nll <- tsmarch:::.dcc_constant_values(pars, spec, type = "nll")
  
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
    mv_nll_vec <- tsmarch:::.copula_dynamic_values(pars, spec, type = "ll_vec")
    garch_nll_vec <- .get_garch_nll_vec(object$spec$univariate)
    total_nll_vec <- garch_nll_vec + mv_nll_vec
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
      mv_nll_vec <- tsmarch:::.copula_constant_values(NULL, object$spec, type = "ll_vec")
      garch_nll_vec <- .get_garch_nll_vec(object$spec$univariate)
      total_nll_vec <- garch_nll_vec + mv_nll_vec
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

## Internal method for GOGARCH
.compute_gogarch_loglik <- function(object, params, return_components = FALSE) {
  
  ## GOGARCH is different - we need to update the univariate GARCH models
  ## This is more complex because we need to re-estimate with fixed parameters
  ## For now, we'll use the estimated models as-is and compute the likelihood
  
  stop("GOGARCH fixed parameter log-likelihood computation is not yet implemented.\n",
       "GOGARCH requires re-estimation of independent component GARCH models with fixed parameters.\n",
       "This functionality will be added in a future version.")
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