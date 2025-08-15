## R helper functions to be called from the C++ EM orchestrator.
## These functions provide the R-side statistical computations needed for the M-step,
## correctly interfacing with the tsgarch/tsmarch packages and handling the
## ARMA/VARMA estimation without non-CRAN dependencies.


#' Calculate the Log-Likelihood Vector
#'
#' This function calculates the vector of log-likelihood contributions for a given
#' set of mean and variance parameters. It's called during the E-step.
#'
#' @param y The time series data (a numeric vector or matrix).
#' @param current_pars A list containing the current 'arma_pars' or 'var_pars', and 'garch_pars'.
#' @param spec The model specification list for the current state.
#' @param model_type A string, either "univariate" or "multivariate".
#' @return A numeric vector of log-likelihoods, one for each time point.
calculate_loglik_vector_r <- function(y, current_pars, spec, model_type = "univariate") {
  ## 1. Get Residuals from the Conditional Mean Model
  if (model_type == "univariate") {
    arma_order <- spec$arma_order
    arma_pars <- current_pars$arma_pars
    
    ## Centering y is not needed when we set include.mean = FALSE in stats::arima()
    # n_arma <- sum(arma_order)
    # has_intercept <- length(arma_pars) > n_arma
    # if (has_intercept) {
    #   mean_val <- arma_pars[n_arma + 1]
    #   arma_filter_pars <- arma_pars[1:n_arma]
    # } else {
    #   mean_val <- 0
    #   arma_filter_pars <- arma_pars
    # }
    # y_centered <- y - mean_val

    # model_residuals <- stats::arima(y_centered, order = c(arma_order[1], 0, arma_order[2]),
    #                                 fixed = arma_filter_pars, include.mean = FALSE)$residuals
    
    model_residuals <- stats::arima(y, order = c(arma_order[1], 0, arma_order[2]),
                                    fixed = arma_pars, include.mean = FALSE)$residuals
  } else {
    ## Multivariate case: VAR(p) residuals
    var_order <- spec$var_order ## e.g., 1 for VAR(1)
    k <- ncol(y)
    T_obs <- nrow(y)
    X_lagged <- matrix(1, nrow = T_obs - var_order, ncol = 1 + k * var_order)
    for (i in 1:var_order) {
      X_lagged[, (2 + (i - 1) * k):(1 + i * k)] <- y[(var_order - i + 1):(T_obs - i), ]
    }
    y_target <- y[(var_order + 1):T_obs, ]
    beta_mat <- matrix(current_pars$var_pars, nrow = 1 + k * var_order, ncol = k)
    model_residuals <- y_target - X_lagged %*% beta_mat
  }

  ## 2. Get GARCH log-likelihood from the residuals
  if (model_type == "univariate") {
    if (!requireNamespace("tsgarch", quietly = TRUE)) stop("The 'tsgarch' package is required.")
    # Convert to xts for tsgarch (dates are arbitrary)
    residuals_xts <- xts::xts(model_residuals, order.by = Sys.Date() - (length(model_residuals):1))
    
    garch_spec_obj <- tsgarch::garch_modelspec(y = residuals_xts,
                                               model = spec$garch_model,
                                               garch_order = spec$garch_order,
                                               #fixed_pars = current_pars$garch_pars,
                                               distribution = spec$distribution)
    fixed_garch_pars <- current_pars$garch_pars
    for (par_name in names(fixed_garch_pars)) {
      garch_spec_obj$parmatrix[parameter == par_name, value := fixed_garch_pars[[par_name]]]
      garch_spec_obj$parmatrix[parameter == par_name, estimate := 0]
    }

    # garch_model_fit <- estimate(garch_spec_obj)
    # report <- garch_model_fit$TMB_OBJECT$report()
    
    ## tsfilter() is to apply a model with already-specified parameters to a 
    ## dataset to "filter" it, producing the residuals and conditional 
    ## volatility (sigma) for each time point.
    ## In the E-step, the calculate_loglik_vector_r function receives the 
    ## current parameters for a state. For this step, these parameters are 
    ## considered fixed. We just need to see how well they explain the data.
    ## 
    ## (When we pass a specification to estimate() where all parameters are 
    ## fixed (by setting estimate := 0 in the parmatrix), it recognizes that 
    ## there's nothing to optimize and internally dispatches to tsfilter() 
    ## anyway.)
    garch_model_fit <- tsfilter(garch_spec_obj)
    res <- as.numeric(garch_model_fit$residuals)
    sig <- as.numeric(garch_model_fit$sigma)
    ll_vector <- dnorm(res, mean = 0, sd = sig, log = TRUE)
  } else {
    if (!requireNamespace("tsmarch", quietly = TRUE)) stop("The 'tsmarch' package is required.")
    # Convert to xts for tsgarch (dates are arbitrary)
    residuals_xts <- xts::xts(model_residuals, order.by = Sys.Date() - (nrow(model_residuals):1))
    
    ## Dynamically call the specified tsmarch model function (e.g., dcc_modelspec)
    spec_fun_name <- spec$garch_spec_fun ## e.g., "dcc_modelspec"
    spec_fun <- get(spec_fun_name, asNamespace("tsmarch"))
    
    ## Prepare arguments for the spec function
    spec_args <- c(list(y = residuals_xts), spec$garch_spec_args)

    ## Create the spec object
    garch_spec_obj <- do.call(spec_fun, spec_args)

    # garch_model_fit <- tsmarch::estimate(garch_spec_obj)
    # report <- garch_model_fit$TMB_OBJECT$report()
    fixed_garch_pars <- current_pars$garch_pars
    for (par_name in names(fixed_garch_pars)) {
      garch_spec_obj$parmatrix[parameter == par_name, value := fixed_garch_pars[[par_name]]]
      garch_spec_obj$parmatrix[parameter == par_name, estimate := 0]
    }
    garch_model_fit <- tsfilter(garch_spec_obj)
    # Placeholder for multivariate log-likelihood calculation
    ll_vector <- rep(mean(dnorm(garch_model_fit$residuals, log = TRUE), na.rm = TRUE), nrow(garch_model_fit$residuals))
  }

  ## Pad the ll_vector to match original data length
  # ll_vector <- report$ll_vector
  # padding <- nrow(y) - length(ll_vector)
  # return(c(rep(0, padding), ll_vector))
  ll_vector[!is.finite(ll_vector)] <- 0
  if (length(ll_vector) < NROW(y)) {
    ## For vectors, NROW() returns length, while nrow() returns NULL
    padding <- NROW(y) - length(ll_vector)
    ll_vector <- c(rep(0, padding), ll_vector)
  }
  return(ll_vector)
}


#' Estimate Conditional Mean Parameters with Weighted Likelihood
#'
#' This function estimates ARMA(p,q) or VAR(p) parameters using a weighted
#' objective function. It's called during the M-step.
#'
#' @param y The time series data.
#' @param weights A numeric vector of weights from the E-step.
#' @param spec The model specification for the current state.
#' @param model_type A string, either "univariate" or "multivariate".
#' @return A list containing the estimated 'coefficients' and the 'residuals'.
estimate_arma_weighted_r <- function(y, weights, spec, model_type = "univariate") {
  if (model_type == "univariate") {
    arma_order <- spec$arma_order
    #start_pars <- spec$start_pars$arma_pars
    start_pars <- unlist(spec$start_pars$arma_pars)
    
    # Objective function for optim, using Kalman Filter for likelihood
    # weighted_arma_loglik <- function(params, y_data, arma_order, w) {
    #   mod <- try(stats::KalmanRun(y_data, model = list(AR = params[1:arma_order[1]], MA = params[(arma_order[1]+1):sum(arma_order)])), silent = TRUE)
    #   if (inherits(mod, "try-error")) return(1e10)
    #   ll_vec <- dnorm(mod$resid, mean = 0, sd = sqrt(mod$var), log = TRUE)
    #   return(-sum(w * ll_vec, na.rm = TRUE))
    # }
    
    weighted_arma_loglik <- function(params, y_data, arma_order, w) {
      p <- arma_order[1]
      q <- arma_order[2]
      ar_params <- if (p > 0) params[1:p] else numeric(0)
      ma_params <- if (q > 0) params[(p+1):(p+q)] else numeric(0)
      
      mod <- try(stats::KalmanRun(y_data, model = list(AR = ar_params, MA = ma_params)), silent = TRUE)
      if (inherits(mod, "try-error")) return(1e10)
      
      ll_vec <- dnorm(mod$resid, mean = 0, sd = sqrt(mod$var), log = TRUE)
      return(-sum(w * ll_vec, na.rm = TRUE))
    }
    
    opt_result <- try(stats::optim(
      par = start_pars, 
      fn = weighted_arma_loglik, 
      y_data = y, 
      arma_order = arma_order, 
      w = weights, 
      method = "BFGS"
    ), silent = TRUE)
    
    #estimated_coeffs <- opt_result$par
    if (inherits(opt_result, "try-error")) {
      estimated_coeffs <- start_pars
    } else {
      estimated_coeffs <- opt_result$par
    }
    
    # Get final residuals
    final_residuals <- stats::arima(
      y, 
      order = c(arma_order[1], 0, arma_order[2]), 
      fixed = estimated_coeffs,
      include.mean = FALSE
    )$residuals
    final_residuals <- as.matrix(final_residuals)
  } else {
    ## Multivariate: Weighted Least Squares for VAR(p)
    var_order <- spec$var_order
    k <- ncol(y)
    T_obs <- nrow(y)
    padding <- var_order
    X_lagged <- matrix(1, nrow = T_obs - padding, ncol = 1 + k * var_order)
    for (i in 1:var_order) {
      X_lagged[, (2 + (i-1)*k):(1 + i*k)] <- y[(padding-i+1):(T_obs-i), ]
    }
    y_target <- y[(padding+1):T_obs, ]
    w_target <- weights[(padding+1):T_obs]
    
    ## Apply weights
    X_w <- X_lagged * sqrt(w_target)
    y_w <- y_target * sqrt(w_target)
    
    ## Solve WLS
    #beta_mat <- solve(crossprod(X_w), crossprod(X_w, y_w))
    beta_mat <- try(solve(crossprod(X_w), crossprod(X_w, y_w)), silent = TRUE)
    if (inherits(beta_mat, "try-error")) {
      # Fallback if matrix is singular
      beta_mat <- matrix(0, nrow = 1 + k*var_order, ncol = k)
    }
    estimated_coeffs <- as.vector(beta_mat)
    final_residuals <- y_target - X_lagged %*% beta_mat
  }
  
  return(list(coefficients = estimated_coeffs, residuals = final_residuals))
}


#' Estimate GARCH Parameters with Weighted Likelihood
#'
#' This function estimates GARCH parameters using weighted maximum likelihood on a
#' given series of residuals. It's called during the M-step.
#'
#' @param residuals The residuals from the updated mean model.
#' @param weights A numeric vector of weights from the E-step.
#' @param spec The model specification for the current state.
#' @param model_type A string, either "univariate" or "multivariate".
#' @return A list containing the estimated 'coefficients'.
estimate_garch_weighted_r <- function(residuals, weights, spec, model_type = "univariate") {
  start_pars <- spec$start_pars$garch_pars
  padding <- NROW(residuals) - length(weights)
  w_target <- weights ## Assume weights are already aligned with residuals
  if(padding > 0) w_target <- weights[(padding + 1):length(weights)]
  
  ## Get bounds from a temporary spec object ---
  ## Create a dummy xts object to initialize the spec
  temp_residuals_xts <- xts::xts(residuals, order.by = Sys.Date() - (length(residuals):1))
  #temp_spec_obj <- tsgarch::garch_modelspec(y = temp_residuals_xts, model = spec$garch_model, garch_order = spec$garch_order)
  if (model_type == "univariate") {
    temp_spec_obj <- tsgarch::garch_modelspec(y = temp_residuals_xts, model = spec$garch_model, garch_order = spec$garch_order)
  } else {
    spec_fun_name <- spec$garch_spec_fun
    spec_fun <- get(spec_fun_name, asNamespace("tsmarch"))
    spec_args <- c(list(y = temp_residuals_xts), spec$garch_spec_args)
    temp_spec_obj <- do.call(spec_fun, spec_args)
  }
  
  ## Extract the parameter matrix
  parmatrix <- temp_spec_obj$parmatrix
  
  ## Filter for the parameters we are actually estimating
  pars_to_estimate <- names(start_pars)
  bounds_matrix <- parmatrix[parameter %in% pars_to_estimate]
  
  ## Ensure the order is correct
  bounds_matrix <- bounds_matrix[match(pars_to_estimate, parameter),]
  
  lower_bounds <- bounds_matrix$lower
  upper_bounds <- bounds_matrix$upper
  
  weighted_garch_loglik <- function(params, residuals_data, w, spec, model_type) {
    param_list <- as.list(params)
    names(param_list) <- names(spec$start_pars$garch_pars)
    
    ## Convert to xts for the garch spec functions
    residuals_xts <- xts::xts(residuals_data, order.by = Sys.Date() - (length(residuals_data):1))
    
    if (model_type == "univariate") {
      if (!requireNamespace("tsgarch", quietly = TRUE)) stop("The 'tsgarch' package is required.")
      garch_spec_obj <- tsgarch::garch_modelspec(y = residuals_xts, model = spec$garch_model, garch_order = spec$garch_order, fixed_pars = param_list, distribution = spec$distribution)
      for (par_name in names(param_list)) {
        garch_spec_obj$parmatrix[parameter == par_name, value := param_list[[par_name]]]
      }
      fit <- try(estimate(garch_spec_obj), silent = TRUE)
    } else {
      if (!requireNamespace("tsmarch", quietly = TRUE)) stop("The 'tsmarch' package is required.")
      spec_fun_name <- spec$garch_spec_fun
      spec_fun <- get(spec_fun_name, asNamespace("tsmarch"))
      #spec_args <- cc(list(y = residuals_xts, fixed_pars = param_list), spec$garch_spec_args)
      spec_args <- cc(list(y = residuals_xts), spec$garch_spec_args)
      garch_spec_obj <- do.call(spec_fun, spec_args)
      for (par_name in names(param_list)) {
        garch_spec_obj$parmatrix[parameter == par_name, value := param_list[[par_name]]]
      }
      fit <- try(estimate(garch_spec_obj), silent = TRUE)
    }
    
    if (inherits(fit, "try-error")) return(1e10)
    
    #ll_vector <- fit$TMB_OBJECT$report()$ll_vector
    if ("TMB_OBJECT" %in% names(fit)) {
      ll_vector <- fit$TMB_OBJECT$report()$ll_vector
    } else {
      res <- as.numeric(fit$residuals)
      sig <- as.numeric(fit$sigma)
      ll_vector <- dnorm(res, mean = 0, sd = sig, log = TRUE)
    }
    
    return(-sum(w * ll_vector, na.rm = TRUE))
  }
  
  ## Run the optimization using a constrained method with the correct bounds
  opt_result <- stats::optim(par = unlist(start_pars), 
                             fn = weighted_garch_loglik,
                             lower = lower_bounds,
                             upper = upper_bounds,
                             method = "L-BFGS-B",
                             residuals_data = residuals, 
                             w = w_target, 
                             spec = spec, 
                             model_type = model_type)
  
  estimated_coeffs <- as.list(opt_result$par)
  names(estimated_coeffs) <- names(start_pars)
  
  return(list(coefficients = estimated_coeffs))
}
