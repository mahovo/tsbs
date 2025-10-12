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
#'   parameter names in the model's parmatrix. For DCC/Copula models, typically
#'   includes "alpha", "beta", and optionally "gamma" (for ADCC) and "shape"
#'   (for Student-t). For GOGARCH, this should include the GARCH parameters
#'   for each independent component.
#' @param return_components Logical. If TRUE, returns both the total log-likelihood
#'   and its components (univariate GARCH + multivariate). Default is FALSE.
#' @param ... Additional arguments (currently unused)
#'
#' @return If return_components = FALSE (default), returns a single numeric value
#'   representing the total log-likelihood. If return_components = TRUE, returns
#'   a list with components:
#'   \item{loglik}{Total log-likelihood (note: NOT negative log-likelihood)}
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
#' \dontrun{
#' # Estimate a DCC model
#' library(tsmarch)
#' data(dmbp)
#' spec <- dcc_modelspec(dmbp, dcc_order = c(1,1))
#' fit <- estimate(spec)
#' 
#' # Compute log-likelihood at different parameter values
#' ll1 <- compute_loglik_fixed(fit, params = list(alpha_1 = 0.05, beta_1 = 0.90))
#' ll2 <- compute_loglik_fixed(fit, params = list(alpha_1 = 0.03, beta_1 = 0.95))
#' 
#' # Compare with estimated log-likelihood
#' logLik(fit)
compute_loglik_fixed <- function(object, params, return_components = FALSE, ...) {
  
  # Validate input object
  if (!inherits(object, c("dcc.estimate", "cgarch.estimate", "gogarch.estimate"))) {
    stop("object must be an estimated tsmarch model of class 'dcc.estimate', 'cgarch.estimate', or 'gogarch.estimate'")
  }
  
  if (!is.list(params)) {
    stop("params must be a named list")
  }
  
  # Dispatch to appropriate method
  if (inherits(object, "dcc.estimate")) {
    if (inherits(object, "dcc.dynamic")) {
      result <- .compute_dcc_dynamic_loglik(object, params, return_components)
    } else if (inherits(object, "dcc.constant")) {
      result <- .compute_dcc_constant_loglik(object, params, return_components)
    } else {
      stop("Unknown DCC model type")
    }
  } else if (inherits(object, "cgarch.estimate")) {
    if (inherits(object, "cgarch.dynamic")) {
      result <- .compute_copula_dynamic_loglik(object, params, return_components)
    } else if (inherits(object, "cgarch.constant")) {
      result <- .compute_copula_constant_loglik(object, params, return_components)
    } else {
      stop("Unknown Copula-GARCH model type")
    }
  } else if (inherits(object, "gogarch.estimate")) {
    result <- .compute_gogarch_loglik(object, params, return_components)
  }
  
  return(result)
}

# Internal method for DCC dynamic
.compute_dcc_dynamic_loglik <- function(object, params, return_components = FALSE) {
  
  # Create a copy of the spec
  spec <- object$spec
  
  # Update parameters in the parmatrix
  spec$parmatrix <- .update_parmatrix(spec$parmatrix, params)
  
  # Get parameter values for the multivariate component
  estimate <- NULL
  pars <- spec$parmatrix[estimate == 1]$value
  
  # Compute multivariate log-likelihood using internal tsmarch function
  # Note: .dcc_dynamic_values returns NEGATIVE log-likelihood (NLL)
  mv_nll <- tsmarch:::.dcc_dynamic_values(pars, spec, type = "nll")
  
  # Get univariate GARCH negative log-likelihoods
  garch_nll <- sum(sapply(object$spec$univariate, function(x) x$loglik))
  
  # Total negative log-likelihood (what tsmarch stores in object$loglik)
  total_nll <- garch_nll + mv_nll
  
  if (return_components) {
    return(list(
      loglik = -total_nll,          ## Convert to positive ll
      garch_loglik = -garch_nll,    ## Convert to positive ll
      multivariate_loglik = -mv_nll ## Convert to positive ll
    ))
  } else {
    # Return positive log-likelihood (to match logLik() method)
    return(-total_nll)
  }
}

# Internal method for DCC constant
.compute_dcc_constant_loglik <- function(object, params, return_components = FALSE) {
  
  # For constant correlation with mvn, there are no parameters to estimate
  if (object$spec$distribution == "mvn" && length(params) == 0) {
    if (return_components) {
      garch_nll <- sum(sapply(object$spec$univariate, function(x) as.numeric(logLik(x))))
      mv_nll <- tsmarch:::.dcc_constant_values(NULL, object$spec, type = "nll")
      total_nll <- garch_nll + mv_nll
      return(list(
        loglik = -total_nll,
        garch_loglik = -garch_nll,
        multivariate_loglik = -mv_nll
      ))
    } else {
      return(-object$loglik)
    }
  }
  
  # Create a copy of the spec
  spec <- object$spec
  
  # Update parameters in the parmatrix
  spec$parmatrix <- .update_parmatrix(spec$parmatrix, params)
  
  # Get parameter values
  estimate <- NULL
  pars <- spec$parmatrix[estimate == 1]$value
  
  # Compute multivariate log-likelihood
  mv_nll <- tsmarch:::.dcc_constant_values(pars, spec, type = "nll")
  
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

# Internal method for Copula dynamic
# .compute_copula_dynamic_loglik <- function(object, params, return_components = FALSE) {
# 
#   # Create a copy of the spec
#   spec <- object$spec
# 
#   # Update parameters in the parmatrix
#   spec$parmatrix <- .update_parmatrix(spec$parmatrix, params)
# 
#   # Get parameter values
#   estimate <- NULL
#   pars <- spec$parmatrix[estimate == 1]$value
# 
#   # Compute multivariate log-likelihood (copula component)
#   # .copula_dynamic_values returns NEGATIVE log-likelihood
#   mv_nll <- tsmarch:::.copula_dynamic_values(pars, spec, type = "nll")
# 
#   if (return_components) {
#     garch_nll <- sum(sapply(object$spec$univariate, function(x) as.numeric(logLik(x))))
#     total_nll <- garch_nll + mv_nll
# 
#     return(list(
#       loglik = -total_nll,
#       garch_loglik = -garch_nll,
#       multivariate_loglik = -mv_nll
#     ))
#   } else {
#     garch_nll <- sum(sapply(object$spec$univariate, function(x) as.numeric(logLik(x))))
#     total_nll <- garch_nll + mv_nll
#     return(-total_nll)
#   }
# }
.compute_copula_dynamic_loglik <- function(object, params, return_components = FALSE) {

  # Create a copy of the spec
  spec <- object$spec

  # Update parameters in the parmatrix
  spec$parmatrix <- .update_parmatrix(spec$parmatrix, params)

  # Get parameter values
  estimate <- NULL
  pars <- spec$parmatrix[estimate == 1]$value

  # Compute multivariate negative log-likelihood (copula component)
  mv_nll <- tsmarch:::.copula_dynamic_values(pars, spec, type = "nll")

  # Get univariate GARCH negative log-likelihoods
  garch_nll <- sum(sapply(object$spec$univariate, function(x) x$loglik))

  # Total negative log-likelihood
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

# Internal method for Copula constant
.compute_copula_constant_loglik <- function(object, params, return_components = FALSE) {
  
  # For constant correlation with mvn, there are no parameters to estimate
  if (object$spec$copula == "mvn" && length(params) == 0) {
    if (return_components) {
      garch_nll <- sum(sapply(object$spec$univariate, function(x) as.numeric(logLik(x))))
      mv_nll <- tsmarch:::.copula_constant_values(NULL, object$spec, type = "nll")
      total_nll <- garch_nll + mv_nll
      return(list(
        loglik = -total_nll,
        garch_loglik = -garch_nll,
        multivariate_loglik = -mv_nll
      ))
    } else {
      return(-object$loglik)
    }
  }
  
  # Create a copy of the spec
  spec <- object$spec
  
  # Update parameters in the parmatrix
  spec$parmatrix <- .update_parmatrix(spec$parmatrix, params)
  
  # Get parameter values
  estimate <- NULL
  pars <- spec$parmatrix[estimate == 1]$value
  
  # Compute multivariate log-likelihood
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
.compute_gogarch_loglik <- function(object, params, return_components = FALSE) {
  
  # GOGARCH is different - we need to update the univariate GARCH models
  # This is more complex because we need to re-estimate with fixed parameters
  # For now, we'll use the estimated models as-is and compute the likelihood
  
  stop("GOGARCH fixed parameter log-likelihood computation is not yet implemented.\n",
       "GOGARCH requires re-estimation of independent component GARCH models with fixed parameters.\n",
       "This functionality will be added in a future version.")
}

# Helper function to update parmatrix with new parameter values
.update_parmatrix <- function(parmatrix, params) {
  
  parameter <- NULL
  
  # Make a copy to avoid modifying the original
  pmatrix <- data.table::copy(parmatrix)
  
  # Update each parameter
  for (param_name in names(params)) {
    # Find matching parameter(s) in parmatrix
    matches <- which(pmatrix$parameter == param_name)
    
    if (length(matches) == 0) {
      warning("Parameter '", param_name, "' not found in model parmatrix. Ignoring.")
    } else if (length(matches) > 1) {
      # Multiple parameters with same name (shouldn't happen in DCC/Copula)
      warning("Multiple parameters named '", param_name, "' found. Updating all.")
      pmatrix[matches, ]$value <- params[[param_name]]
    } else {
      pmatrix[matches, ]$value <- params[[param_name]]
    }
  }
  
  return(pmatrix)
}