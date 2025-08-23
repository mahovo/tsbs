## ===================================================================
## Generalized R Helper Functions for MS-ARMA-GARCH Fitting
## ===================================================================
## These functions are designed to be called from the C++ EM orchestrator.
## They contain the full logic for handling general ARMA(p,q) and VAR(p)
## models, and interface with the tsgarch/tsmarch packages.


#' @title Create a GARCH Specification Object (Convenience Function)
#' @description This is a new convenience function that handles the complex,
#' model-specific logic for creating both univariate and multivariate GARCH
#' specification objects. It isolates the complexity of the tsmarch workflow.
create_garch_spec_object_r <- function(residuals, spec, model_type) {
  residuals_xts <- xts::xts(residuals, order.by = Sys.Date() - (NROW(residuals):1))
  
  if (model_type == "univariate") {
    garch_spec_obj <- tsgarch::garch_modelspec(y = residuals_xts, 
                                               model = spec$garch_model,
                                               garch_order = spec$garch_order, 
                                               distribution = spec$distribution)
  } else {
    spec_fun_name <- spec$garch_spec_fun
    spec_fun <- get(spec_fun_name, asNamespace("tsmarch"))
    spec_args <- spec$garch_spec_args
    
    if (spec_fun_name %in% c("dcc_modelspec", "cgarch_modelspec")) {
      univariate_models <- lapply(1:ncol(residuals_xts), function(i) {
        uni_spec <- spec_args$garch_model$univariate[[i]]
        suppressWarnings({
          estimate(tsgarch::garch_modelspec(y = residuals_xts[,i], 
                                            model = uni_spec$model, 
                                            garch_order = uni_spec$garch_order,
                                            distribution = "norm"),
                   keep_tmb = TRUE)
        })
      })
      garch_model_est <- to_multi_estimate(univariate_models)
      names(garch_model_est) <- paste0("series_", 1:ncol(residuals_xts))
      
      final_args <- c(list(object = garch_model_est), 
                      spec_args[names(spec_args) != "garch_model"])
      garch_spec_obj <- do.call(spec_fun, final_args)
      
    } else {
      final_args <- c(list(y = residuals_xts), spec_args)
      garch_spec_obj <- do.call(spec_fun, final_args)
    }
  }
  return(garch_spec_obj)
}


#' @title Calculate the Log-Likelihood Vector (R Helper)
calculate_loglik_vector_r <- function(y, current_pars, spec, model_type = "univariate") {
  
  ## 1. Get Residuals from the Conditional Mean Model
  if (model_type == "univariate") {
    arma_pars <- current_pars$arma_pars
    model_residuals <- stats::arima(y, order = c(spec$arma_order[1], 0, spec$arma_order[2]),
                                    fixed = arma_pars, include.mean = FALSE)$residuals
  } else {
    var_order <- spec$var_order
    k <- ncol(y)
    T_obs <- nrow(y)
    X_lagged <- matrix(1, nrow = T_obs - var_order, ncol = 1 + k * var_order)
    for (i in 1:var_order) {
      X_lagged[, (2 + (i-1)*k):(1 + i*k)] <- y[(var_order-i+1):(T_obs-i), ]
    }
    y_target <- y[(var_order+1):T_obs, ]
    beta_mat <- matrix(current_pars$var_pars, nrow = 1 + k * var_order, ncol = k)
    model_residuals <- y_target - X_lagged %*% beta_mat
  }
  
  ## 2. Get GARCH log-likelihood from the residuals
  garch_spec_obj <- create_garch_spec_object_r(model_residuals, spec, model_type)
  
  fixed_garch_pars <- current_pars$garch_pars
  for (par_name in names(fixed_garch_pars)) {
    garch_spec_obj$parmatrix[parameter == par_name, value := fixed_garch_pars[[par_name]]]
    garch_spec_obj$parmatrix[parameter == par_name, estimate := 0]
  }
  
  if (model_type == "univariate") {
    garch_model_fit <- tsfilter(garch_spec_obj)
    res <- as.numeric(garch_model_fit$residuals)
    sig <- as.numeric(garch_model_fit$sigma)
    ll_vector <- dnorm(res, mean = 0, sd = sig, log = TRUE)
  } else {
    garch_model_fit <- estimate(garch_spec_obj)
    k <- ncol(model_residuals)
    T_res <- nrow(model_residuals)
    ll_vector <- numeric(T_res)
    
    H_vectorized <- garch_model_fit$H
    for (t in 1:T_res) {
      cov_mat <- matrix(0, k, k)
      cov_mat[upper.tri(cov_mat, diag = TRUE)] <- H_vectorized[t, ]
      cov_mat <- cov_mat + t(cov_mat)
      diag(cov_mat) <- diag(cov_mat) / 2
      
      ll_vector[t] <- mvtnorm::dmvnorm(model_residuals[t,], mean = rep(0, k), 
                                       sigma = cov_mat, log = TRUE)
    }
  }
  
  ## Sanitize and pad the vector before returning to C++
  ll_vector[!is.finite(ll_vector)] <- 0
  if (length(ll_vector) < NROW(y)) {
    padding <- NROW(y) - length(ll_vector)
    ll_vector <- c(rep(0, padding), ll_vector)
  }
  return(ll_vector)
}


#' @title Estimate Conditional Mean Parameters (R Helper)
estimate_arma_weighted_r <- function(y, weights, spec, model_type = "univariate") {
  if (model_type == "univariate") {
    arma_order <- spec$arma_order
    start_pars <- unlist(spec$start_pars$arma_pars)
    weighted_arma_loglik <- function(params, y_data, arma_order, w) {
      p <- arma_order[1]; q <- arma_order[2]
      ar_params <- if (p > 0) params[1:p] else numeric(0)
      ma_params <- if (q > 0) params[(p+1):(p+q)] else numeric(0)
      mod <- try(stats::KalmanRun(y_data, model = list(AR = ar_params, MA = ma_params)), silent = TRUE)
      if (inherits(mod, "try-error")) return(1e10)
      ll_vec <- dnorm(mod$resid, mean = 0, sd = sqrt(mod$var), log = TRUE)
      return(-sum(w * ll_vec, na.rm = TRUE))
    }
    opt_result <- try(stats::optim(par = start_pars, fn = weighted_arma_loglik, y_data = y,
                                   arma_order = arma_order, w = weights, method = "BFGS"), silent = TRUE)
    estimated_coeffs <- if (inherits(opt_result, "try-error")) start_pars else opt_result$par
    final_residuals <- stats::arima(y, order = c(arma_order[1], 0, arma_order[2]),
                                    fixed = estimated_coeffs, include.mean = FALSE)$residuals
    final_residuals <- as.matrix(final_residuals)
  } else {
    var_order <- spec$var_order
    k <- ncol(y); T_obs <- nrow(y); padding <- var_order
    X_lagged <- matrix(1, nrow = T_obs - padding, ncol = 1 + k * var_order)
    for (i in 1:var_order) {
      X_lagged[, (2 + (i-1)*k):(1 + i*k)] <- y[(padding-i+1):(T_obs-i), ]
    }
    y_target <- y[(padding+1):T_obs, ]; w_target <- weights[(padding+1):T_obs]
    X_w <- X_lagged * sqrt(w_target); y_w <- y_target * sqrt(w_target)
    beta_mat <- try(solve(crossprod(X_w), crossprod(X_w, y_w)), silent = TRUE)
    if (inherits(beta_mat, "try-error")) {
      beta_mat <- matrix(0, nrow = 1 + k*var_order, ncol = k)
    }
    estimated_coeffs <- as.vector(beta_mat)
    final_residuals <- y_target - X_lagged %*% beta_mat
  }
  return(list(coefficients = estimated_coeffs, residuals = final_residuals))
}


#' @title Estimate GARCH Parameters (R Helper)
estimate_garch_weighted_r <- function(residuals, weights, spec, model_type = "univariate") {
  start_pars <- spec$start_pars$garch_pars
  padding <- NROW(residuals) - length(weights)
  w_target <- if(padding > 0) weights[(padding+1):length(weights)] else weights
  
  temp_spec_obj <- create_garch_spec_object_r(residuals, spec, model_type)
  
  parmatrix <- temp_spec_obj$parmatrix
  pars_to_estimate <- names(start_pars)
  bounds_matrix <- parmatrix[parameter %in% pars_to_estimate]
  bounds_matrix <- bounds_matrix[match(pars_to_estimate, parameter),]
  
  lower_bounds <- bounds_matrix$lower
  upper_bounds <- bounds_matrix$upper
  
  weighted_garch_loglik <- function(params, residuals_data, w, spec, model_type) {
    param_list <- as.list(params)
    names(param_list) <- names(spec$start_pars$garch_pars)
    
    garch_spec_obj <- create_garch_spec_object_r(residuals_data, spec, model_type)
    
    for (par_name in names(param_list)) {
      garch_spec_obj$parmatrix[parameter == par_name, value := param_list[[par_name]]]
    }
    
    fit <- try(estimate(garch_spec_obj), silent = TRUE)
    if (inherits(fit, "try-error")) return(1e10)
    
    if ("TMB_OBJECT" %in% names(fit)) {
      ll_vector <- fit$TMB_OBJECT$report()$ll_vector
    } else {
      res <- as.numeric(fit$residuals)
      sig <- as.numeric(fit$sigma)
      ll_vector <- dnorm(res, mean = 0, sd = sig, log = TRUE)
    }
    return(-sum(w * ll_vector, na.rm = TRUE))
  }
  
  ## ---- Capture and summarize warnings from optim ----
  warnings_list <- list()
  opt_result <- withCallingHandlers({
    stats::optim(par = unlist(start_pars), 
                 fn = weighted_garch_loglik,
                 lower = lower_bounds,
                 upper = upper_bounds,
                 method = "L-BFGS-B",
                 residuals_data = residuals, 
                 w = w_target, 
                 spec = spec, 
                 model_type = model_type)
  }, warning = function(w) {
    warnings_list <<- c(warnings_list, list(w))
    invokeRestart("muffleWarning")
  })
  
  estimated_coeffs <- as.list(opt_result$par)
  names(estimated_coeffs) <- names(start_pars)
  
  return(list(coefficients = estimated_coeffs, warnings = warnings_list))
}
