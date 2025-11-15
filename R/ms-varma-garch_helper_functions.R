## === === === === === === === === === === === === === === === === === 
## Generalized R Helper Functions for MS-ARMA-GARCH Fitting
## === === === === === === === === === === === === === === === === === 
## These functions are designed to be called from the C++ EM orchestrator
## fit_ms_varma_garch_cpp(). They contain the full logic for handling general 
## ARMA(p,q) and VAR(p) models, and interface with the tsgarch/tsmarch packages.


#' @title Generate Correct Parameter Names for tsmarch
#' @description Translates a nested parameter list into the flat named list
#'              that tsmarch's parmatrix expects (e.g., "omega[1]").
#' @param pars A nested parameter list.
#' @return A flat named list.
#' @keywords internal
generate_tsmarch_parnames <- function(pars) {
  ## Handle the nested GARCH parameters, adding the "[i]" suffix.
  garch_pars_flat <- list()
  if (!is.null(pars$garch_pars)) {
    for (i in 1:length(pars$garch_pars)) {
      series_pars <- pars$garch_pars[[i]]
      if (!is.null(series_pars) && length(series_pars) > 0) {
        names(series_pars) <- paste0(names(series_pars), "[", i, "]")
        garch_pars_flat <- c(garch_pars_flat, series_pars)
      }
    }
  }
  
  ## Handle all other parameters, which are assumed to be in a
  ## flat structure already and do not need suffixing.
  other_pars <- pars[!names(pars) %in% c("var_pars", "garch_pars")]
  
  return(c(garch_pars_flat, other_pars))
}


#' @title Create a GARCH Specification Object (Convenience Function)
#' @description This is a convenience function that handles the complex,
#'   model-specific logic for creating both univariate and multivariate GARCH
#'   specification objects.
#' @param residuals Numeric
#' @param spec A spec list
#' @param model_type Character string
#' @param current_pars A list of parameters
#' @return A GARCH specification object
create_garch_spec_object_r <- function(
    residuals, 
    spec, 
    model_type, 
    current_pars
) {
  residuals_xts <- xts::xts(residuals, order.by = Sys.Date() - (NROW(residuals):1))
  
  if (model_type == "univariate") {
    ## --- Univariate ---
    garch_spec_obj <- tsgarch::garch_modelspec(
      y = residuals_xts,
      model = spec$garch_model,
      garch_order = spec$garch_order,
      distribution = spec$distribution
    )
     
    ## Set parameters from previous M-step
    all_fixed_pars <- c(current_pars$garch_pars, current_pars$dist_pars)
    for (par_name in names(all_fixed_pars)) {
      if (par_name %in% garch_spec_obj$parmatrix$parameter) {
        #garch_spec_obj$parmatrix[parameter == par_name, value := all_fixed_pars[[par_name]]]
        row_idx <- which(garch_spec_obj$parmatrix$parameter == par_name)
        garch_spec_obj$parmatrix[row_idx, "value"] <- all_fixed_pars[[par_name]]
      }
    }
    
    return(garch_spec_obj)
    
  } else {
    ## --- Multivariate ---
    spec_fun_name <- spec$garch_spec_fun
    spec_fun <- get(spec_fun_name, asNamespace("tsmarch"))
    spec_args <- spec$garch_spec_args
    
    if (spec_fun_name %in% c("dcc_modelspec", "cgarch_modelspec")) {
      ## STEP 1: Create specs and filter with fixed parameters (NO estimation)
      univariate_models <- lapply(1:ncol(residuals_xts), function(i) {
        uni_spec_details <- spec_args$garch_model$univariate[[i]]
        uni_spec_obj <- tsgarch::garch_modelspec(
          y = residuals_xts[,i],
          model = uni_spec_details$model,
          garch_order = uni_spec_details$garch_order,
          distribution = uni_spec_details$distribution
        )
        
        ## Estimate to get TMB object (with original parameters)
        estimated_obj <- suppressWarnings(estimate(uni_spec_obj, keep_tmb = TRUE))
        
        ## Now UPDATE the parameter values (but keep estimate = 1)
        series_garch_pars <- current_pars$garch_pars[[i]]
        
        for (par_name in names(series_garch_pars)) {
          if (par_name %in% estimated_obj$parmatrix$parameter) {
            value <- series_garch_pars[[par_name]]
            
            ## Ensure value respects bounds
            row_match <- estimated_obj$parmatrix$parameter == par_name
            lower <- estimated_obj$parmatrix$lower[row_match][1]
            upper <- estimated_obj$parmatrix$upper[row_match][1]
            value <- max(value, lower + 1e-8)
            value <- min(value, upper - 1e-8)
            
            ## Update ONLY the value, keep estimate flag unchanged
            estimated_obj$parmatrix$value[row_match] <- value
          }
        }
        
        ## CRITICAL: Recompute sigma using TMB with updated parameters
        ## Get parameters in the order TMB expects (estimate == 1)
        estimate_col <- NULL
        pars_for_tmb <- estimated_obj$parmatrix[estimate_col == 1]$value
        
        if (!is.null(estimated_obj$tmb) && length(pars_for_tmb) == length(estimated_obj$tmb$par)) {
          ## Evaluate TMB at new parameters
          tmb_report <- estimated_obj$tmb$report(pars_for_tmb)
          
          ## Update sigma (strip initialization period if needed)
          maxpq <- max(estimated_obj$spec$model$order)
          if (maxpq > 0) {
            estimated_obj$sigma <- xts::xts(
              tmb_report$sigma[-(1:maxpq)],
              order.by = index(estimated_obj$sigma)
            )
          } else {
            estimated_obj$sigma <- xts::xts(
              tmb_report$sigma,
              order.by = index(estimated_obj$sigma)
            )
          }
        }
        
        return(estimated_obj)
      })
      
      ## STEP 2: Combine into multi_estimate object
      multi_estimate_object <- tsgarch::to_multi_estimate(univariate_models)
      names(multi_estimate_object) <- paste0("series_", 1:ncol(residuals_xts))
      
      ## STEP 3: Create DCC spec from the multi_estimate object
      final_args <- c(
        list(object = multi_estimate_object), 
        spec_args[names(spec_args) != "garch_model"]
      )
      final_args$distribution <- spec$distribution
      
      ## Must explicitly specify dynamics = "dcc" for dynamic correlation
      ## Without this, tsmarch defaults to constant correlation!
      if (is.null(final_args$dynamics)) {
        final_args$dynamics <- "dcc"
      }
      
      garch_spec_obj <- do.call(spec_fun, final_args)

      
      ## STEP 4: Set DCC-level parameters (now that spec is created)
      
      ## 4a. Set DCC parameters (alpha_1, beta_1, gamma_1, etc.)
      dcc_param_names <- grep("^(alpha|beta|gamma)_[0-9]+$", names(current_pars), value = TRUE)
      
      for (par_name in dcc_param_names) {
        if (par_name %in% garch_spec_obj$parmatrix$parameter) {
          value <- current_pars[[par_name]]
          
          ## Ensure value respects bounds
          lower <- garch_spec_obj$parmatrix[parameter == par_name]$lower
          upper <- garch_spec_obj$parmatrix[parameter == par_name]$upper
          value <- max(value, lower + 1e-8)
          value <- min(value, upper - 1e-8)
          
          #garch_spec_obj$parmatrix[parameter == par_name, value := value]
          row_idx <- which(garch_spec_obj$parmatrix$parameter == par_name)
          garch_spec_obj$parmatrix[row_idx, "value"] <- value
        }
      }
      
      ## 4b. Set distribution parameters (e.g., shape for MVT)
      if (!is.null(current_pars$dist_pars)) {
        for (par_name in names(current_pars$dist_pars)) {
          if (par_name %in% garch_spec_obj$parmatrix$parameter) {
            value <- current_pars$dist_pars[[par_name]]
            
            lower <- garch_spec_obj$parmatrix[parameter == par_name]$lower
            upper <- garch_spec_obj$parmatrix[parameter == par_name]$upper
            value <- max(value, lower + 1e-8)
            value <- min(value, upper - 1e-8)
            
            #garch_spec_obj$parmatrix[parameter == par_name, value := value]
            row_idx <- which(garch_spec_obj$parmatrix$parameter == par_name)
            garch_spec_obj$parmatrix[row_idx, "value"] <- value
          }
        }
      }
      
      return(garch_spec_obj)
      
    } else {
      ## GOGARCH or other model types
      final_args <- c(list(y = residuals_xts), spec_args)
      garch_spec_obj <- do.call(spec_fun, final_args)
      return(garch_spec_obj)
    }
  }
}


#' @title Calculate the Log-Likelihood Vector (R Helper)
#'
#' @param y 
#' @param current_pars 
#' @param spec 
#' @param model_type 
#' 
#' @import data.table
#' @importFrom tsmethods tsfilter
calculate_loglik_vector_r <- function(
    y, 
    current_pars, 
    spec, 
    model_type = "univariate"
) {
  
  ## 1. Get Residuals from the Conditional Mean Model
  if (model_type == "univariate") {
    arma_pars <- current_pars$arma_pars
    model_residuals <- stats::arima(
      y, 
      order = c(spec$arma_order[1], 0, spec$arma_order[2]),
      fixed = arma_pars, 
      include.mean = FALSE
    )$residuals
  } else {
    var_order <- spec$var_order
    k <- ncol(y)
    T_obs <- nrow(y)
    X_lagged <- matrix(1, nrow = T_obs - var_order, ncol = 1 + k * var_order)
    for (i in 1:var_order) {
      X_lagged[, (2 + (i - 1) * k):(1 + i*k)] <- y[(var_order - i + 1):(T_obs - i), ]
    }
    y_target <- y[(var_order + 1):T_obs, ]
    beta_mat <- matrix(current_pars$var_pars, nrow = 1 + k * var_order, ncol = k)
    model_residuals <- y_target - X_lagged %*% beta_mat
  }
  
  ## 2. Create the GARCH spec object with current parameters
  garch_spec_obj <- create_garch_spec_object_r(
    model_residuals, 
    spec, 
    model_type, 
    current_pars
  )
  
  if (model_type == "univariate") {
    ## --- UNIVARIATE ---
    
    ## Use tsfilter() to get the sigma path
    garch_model_fit <- tsmethods::tsfilter(garch_spec_obj)
    sig <- as.numeric(garch_model_fit$sigma)
    
    ## Call the correct density function
    dist_fun <- switch(spec$distribution,
                       "norm" = stats::dnorm,
                       "snorm" = tsdistributions::dsnorm,
                       "std"  = tsdistributions::dstd,
                       "sstd" = tsdistributions::dsstd,
                       "ged"  = tsdistributions::dged,
                       "sged"  = tsdistributions::dsged,
                       "ghyp"  = tsdistributions::dghyp,
                       "ghst"  = tsdistributions::dghst,
                       "jsu"  = tsdistributions::djsu,
                       stop(paste("Unsupported univariate distribution:", spec$distribution))
    )
    
    if (spec$distribution == "norm") {
      dist_args <- list(x = model_residuals, mean = 0, sd = sig, log = TRUE)
    } else {
      dist_args <- c(
        list(x = model_residuals, mu = 0, sigma = sig, log = TRUE),
        current_pars$dist_pars
      )
    }
    
    ll_vector <- do.call(dist_fun, dist_args)
    
  } else {
    ## --- MULTIVARIATE ---
    ## CRITICAL: We must evaluate likelihood at FIXED parameters,
    ## NOT re-estimate them. This is for the EM algorithm's E-step.
    
    spec_fun_name <- spec$garch_spec_fun
    
    if (spec_fun_name %in% c("dcc_modelspec", "cgarch_modelspec")) {
      ## ==== DCC or Copula-GARCH Models ====
      
      ## Extract DCC-layer parameters from current EM iteration
      other_pars <- current_pars[!names(current_pars) %in% c("var_pars", "garch_pars")]

      ## Inject these into the spec's parmatrix for the internal functions to use
      for (par_name in names(other_pars)) {
        if (par_name %in% garch_spec_obj$parmatrix$parameter) {
          #garch_spec_obj$parmatrix[parameter == par_name, value := other_pars[[par_name]]]
          row_idx <- which(garch_spec_obj$parmatrix$parameter == par_name)
          garch_spec_obj$parmatrix[row_idx, "value"] <- other_pars[[par_name]]
        }
      }
      
      ## Extract parameter vector for internal functions
      estimate <- NULL
      pars <- garch_spec_obj$parmatrix[estimate == 1]$value
      
      ## Determine if dynamic or constant
      is_dynamic <- !is.null(garch_spec_obj$dynamics) && 
        garch_spec_obj$dynamics$model %in% c("dcc", "adcc")
      
      ## Call the appropriate internal tsmarch function
      ## These functions EVALUATE likelihood at given parameters (no optimization!)
      if (spec_fun_name == "dcc_modelspec") {
        if (is_dynamic) {
          ## DCC Dynamic: Call internal function directly
          total_nll_vec <- tsmarch:::.dcc_dynamic_values(
            pars, 
            garch_spec_obj, 
            type = "ll_vec"
          )

          
          ## Strip off initialization period
          dcc_order <- garch_spec_obj$dynamics$order
          maxpq <- max(dcc_order)
          if (maxpq > 0) {
            total_nll_vec <- total_nll_vec[-(1:maxpq), , drop = TRUE]
          } else {
            total_nll_vec <- as.vector(total_nll_vec)
          }
          
          ll_vector <- -total_nll_vec
          
        } else {
          ## DCC Constant
          if (garch_spec_obj$distribution == "mvn" && length(pars) == 0) {
            ## MVN constant correlation has no parameters
            total_nll_vec <- tsmarch:::.dcc_constant_values(
              NULL, 
              garch_spec_obj, 
              type = "ll_vec"
            )
          } else {
            ## MVT constant correlation has shape parameter
            total_nll_vec <- tsmarch:::.dcc_constant_values(
              pars, 
              garch_spec_obj, 
              type = "ll_vec"
            )
          }
          
          ll_vector <- -as.vector(total_nll_vec)
        }
        
      } else if (spec_fun_name == "cgarch_modelspec") {
        if (is_dynamic) {
          ## Copula Dynamic
          copula_nll_vec <- tsmarch:::.copula_dynamic_values(
            pars, 
            garch_spec_obj, 
            type = "ll_vec"
          )
          
          ## Strip off initialization period
          dcc_order <- garch_spec_obj$dynamics$order
          maxpq <- max(dcc_order)
          if (maxpq > 0) {
            copula_nll_vec <- copula_nll_vec[-(1:maxpq), , drop = TRUE]
          } else {
            copula_nll_vec <- as.vector(copula_nll_vec)
          }
          
          ## Get univariate GARCH component
          garch_nll_vec <- .get_garch_nll_vec_from_univariate(garch_spec_obj$univariate)
          
          ## Total = GARCH + Copula
          total_nll_vec <- garch_nll_vec + copula_nll_vec
          ll_vector <- -total_nll_vec
          
        } else {
          ## Copula Constant
          if (garch_spec_obj$copula == "mvn" && length(pars) == 0) {
            copula_nll_vec <- tsmarch:::.copula_constant_values(
              NULL, 
              garch_spec_obj, 
              type = "ll_vec"
            )
          } else {
            copula_nll_vec <- tsmarch:::.copula_constant_values(
              pars, 
              garch_spec_obj, 
              type = "ll_vec"
            )
          }
          
          copula_nll_vec <- as.vector(copula_nll_vec)
          garch_nll_vec <- .get_garch_nll_vec_from_univariate(garch_spec_obj$univariate)
          
          total_nll_vec <- garch_nll_vec + copula_nll_vec
          ll_vector <- -total_nll_vec
        }
      }
      
    } else if (spec_fun_name == "gogarch_modelspec") {
      ## ==== GOGARCH Model ====
      ## GOGARCH is fundamentally different - uses ICA decomposition
      ## We need to estimate once to get the structure, but the parameters
      ## are already fixed by create_garch_spec_object_r()
      
      ## Estimate GOGARCH (this finds ICA decomposition with fixed GARCH params)
      garch_model_fit <- suppressWarnings(estimate(garch_spec_obj, trace = FALSE))
      
      ## Now extract likelihood at these fixed parameters
      ## For GOGARCH, we can use compute_loglik_fixed since the structure is 
      ## different
      ll_vector <- compute_loglik_fixed(
        object = garch_model_fit,
        params = list(),  # Parameters already in the object
        ll_vec = TRUE
      )
      
    } else {
      stop("Unsupported multivariate model type: ", spec_fun_name)
    }
  }
  
  ## Sanitize and pad the vector before returning to C++
  ll_vector[!is.finite(ll_vector)] <- -1e10
  if (length(ll_vector) < NROW(y)) {
    padding <- NROW(y) - length(ll_vector)
    ll_vector <- c(rep(0, padding), ll_vector)
  }
  return(ll_vector)
}


#' @title Estimate Conditional Mean Parameters (R Helper)
#' 
#' @param y 
#' @param weights 
#' @param spec 
#' @param model_type 
estimate_arma_weighted_r <- function(y, weights, spec, model_type = "univariate") {
  if (model_type == "univariate") {
    arma_order <- spec$arma_order
    start_pars <- unlist(spec$start_pars$arma_pars)
    weighted_arma_loglik <- function(params, y_data, arma_order, w) {
      p <- arma_order[1]; q <- arma_order[2]
      ar_params <- if (p > 0) params[1:p] else numeric(0)
      ma_params <- if (q > 0) params[(p+1):(p+q)] else numeric(0)
      mod <- try(
        stats::KalmanRun(y_data, model = list(AR = ar_params, MA = ma_params)), 
        silent = TRUE
      )
      if (inherits(mod, "try-error")) return(1e10)
      ll_vec <- dnorm(mod$resid, mean = 0, sd = sqrt(mod$var), log = TRUE)
      return(-sum(w * ll_vec, na.rm = TRUE))
    }
    opt_result <- try(
      stats::optim(
        par = start_pars, 
        fn = weighted_arma_loglik, 
        y_data = y,
        arma_order = arma_order, 
        w = weights, 
        method = "BFGS"
      ), 
      silent = TRUE
    )
    estimated_coeffs <- if (inherits(opt_result, "try-error")) start_pars else opt_result$par
    final_residuals <- stats::arima(
      y, 
      order = c(arma_order[1], 0, arma_order[2]),
      fixed = estimated_coeffs, 
      include.mean = FALSE
    )$residuals
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


#' @title Estimate GARCH Parameters for Multivariate Models (R Helper)
#' @description Two-stage weighted estimation for multivariate GARCH models:
#'   Stage 1: Univariate GARCH parameters for each series
#'   Stage 2: Multivariate dependence parameters (DCC/Copula) or decomposition (GOGARCH)
#' @param residuals Matrix of residuals (T x k) for multivariate case
#' @param weights Vector of weights from E-step
#' @param spec Model specification
#' @param model_type Either "univariate" or "multivariate"
#' @return List with coefficients and warnings
estimate_garch_weighted_r <- function(residuals, weights, spec, model_type = "univariate") {
  
  if (model_type == "univariate") {
    ## === UNIVARIATE CASE ===
    return(estimate_garch_weighted_univariate(residuals, weights, spec))
  } else {
    ## === MULTIVARIATE CASE ===
    return(estimate_garch_weighted_multivariate(residuals, weights, spec))
  }
}


#' @title Univariate GARCH + Distribution Parameter Estimation
#' @keywords internal
estimate_garch_weighted_univariate <- function(residuals, weights, spec) {
  ## Combine GARCH and Distribution starting parameters into one vector
  start_pars <- c(spec$start_pars$garch_pars, spec$start_pars$dist_pars)
  
  ## Handle case where there are no GARCH parameters to estimate
  if (length(start_pars) == 0) {
    return(list(coefficients = list(), warnings = list()))
  }
  
  padding <- NROW(residuals) - length(weights)
  w_target <- if(padding > 0) weights[(padding+1):length(weights)] else weights
  
  temp_spec_obj <- create_garch_spec_object_r(residuals, spec, "univariate", list())
  
  ## Extract bounds for ALL parameters to be estimated (GARCH + dist)
  parmatrix <- temp_spec_obj$parmatrix
  pars_to_estimate <- names(start_pars)
  bounds_matrix <- parmatrix[parameter %in% pars_to_estimate]
  bounds_matrix <- bounds_matrix[match(pars_to_estimate, parameter),]
  lower_bounds <- bounds_matrix$lower
  upper_bounds <- bounds_matrix$upper
  
  ## Generalized objective function
  weighted_garch_loglik <- function(params, residuals_data, w, spec) {
    param_list <- as.list(params)
    names(param_list) <- names(start_pars)
    
    garch_spec_obj <- create_garch_spec_object_r(residuals_data, spec, "univariate", list())
    
    ## Set all parameter values (GARCH + dist) for this iteration
    for (par_name in names(param_list)) {
      if (par_name %in% garch_spec_obj$parmatrix$parameter) {
        #garch_spec_obj$parmatrix[parameter == par_name, value := param_list[[par_name]]]
        row_idx <- which(garch_spec_obj$parmatrix$parameter == par_name)
        garch_spec_obj$parmatrix[row_idx, "value"] <- param_list[[par_name]]
      }
    }
    
    ## Get the sigma path using the current parameter set
    fit <- try(tsmethods::tsfilter(garch_spec_obj), silent = TRUE)
    if (inherits(fit, "try-error")) return(1e10)
    sig <- as.numeric(fit$sigma)
    
    ## Calculate log-likelihood using the specified distribution
    dist_fun <- switch(spec$distribution,
                       "norm"  = stats::dnorm, "snorm" = tsdistributions::dsnorm,
                       "std"   = tsdistributions::dstd,  "sstd"  = tsdistributions::dsstd,
                       "ged"   = tsdistributions::dged,  "sged"  = tsdistributions::dsged,
                       "ghyp"  = tsdistributions::dghyp, "ghst"  = tsdistributions::dghst,
                       "jsu"   = tsdistributions::djsu,
                       stop(paste("Unsupported univariate distribution:", spec$distribution)))
    
    dist_param_names <- names(spec$start_pars$dist_pars)
    dist_pars_current <- param_list[dist_param_names]
    
    if (spec$distribution == "norm") {
      dist_args <- list(x = residuals_data, mean = 0, sd = sig, log = TRUE)
    } else {
      dist_args <- c(list(x = residuals_data, mu = 0, sigma = sig, log = TRUE), dist_pars_current)
    }
    
    ll_vector <- do.call(dist_fun, dist_args)
    ll_vector[!is.finite(ll_vector)] <- -1e10
    
    ## Return the weighted negative log-likelihood
    return(-sum(w * ll_vector, na.rm = TRUE))
  }
  
  ## Capture warnings
  warnings_list <- list()
  opt_result <- withCallingHandlers({
    stats::optim(par = unlist(start_pars),
                 fn = weighted_garch_loglik,
                 lower = lower_bounds,
                 upper = upper_bounds,
                 method = "L-BFGS-B",
                 residuals_data = residuals,
                 w = w_target,
                 spec = spec)
  }, warning = function(w) {
    warnings_list <<- c(warnings_list, list(w))
    invokeRestart("muffleWarning")
  })
  
  estimated_coeffs <- as.list(opt_result$par)
  names(estimated_coeffs) <- names(start_pars)
  
  return(list(coefficients = estimated_coeffs, warnings = warnings_list))
}


#' @title Multivariate GARCH Parameter Estimation (Two-Stage)
#' @description Implements weighted MLE for multivariate GARCH models
#'   Stage 1: Estimate univariate GARCH parameters for each series
#'   Stage 2: Estimate dependence parameters (DCC/Copula) or GOGARCH rotation
#' @keywords internal
estimate_garch_weighted_multivariate <- function(residuals, weights, spec) {
  
  ## Determine model type
  model_type <- spec$garch_spec_fun
  
  if (model_type == "gogarch_modelspec") {
    return(estimate_garch_weighted_gogarch(residuals, weights, spec))
  } else if (model_type %in% c("dcc_modelspec", "cgarch_modelspec")) {
    return(estimate_garch_weighted_dcc(residuals, weights, spec))
  } else if (model_type == "copula_modelspec") {
    return(estimate_garch_weighted_copula(residuals, weights, spec))
  } else {
    stop(paste("Unsupported multivariate model type:", model_type))
  }
}


#' @title DCC Weighted Estimation
#' @keywords internal
estimate_garch_weighted_dcc <- function(residuals, weights, spec) {
  
  k <- ncol(residuals)
  T_obs <- nrow(residuals)
  
  ## Adjust weights to match residuals length
  ## If residuals are shorter (due to VAR padding), truncate weights from the beginning
  if (length(weights) > T_obs) {
    ## Remove the first (length(weights) - T_obs) elements
    n_to_remove <- length(weights) - T_obs
    w_target <- weights[(n_to_remove + 1):length(weights)]
  } else if (length(weights) < T_obs) {
    stop("Weights vector is shorter than residuals - this should not happen")
  } else {
    w_target <- weights
  }
  
  ## === STAGE 1: Estimate Univariate GARCH for Each Series ===
  garch_pars_list <- list()
  warnings_stage1 <- list()
  
  for (i in 1:k) {
    series_residuals <- residuals[, i]
    series_spec <- spec$garch_spec_args$garch_model$univariate[[i]]
    
    ## Create a univariate spec structure
    uni_spec <- list(
      garch_model = series_spec$model,
      garch_order = series_spec$garch_order,
      distribution = series_spec$distribution,
      start_pars = list(
        garch_pars = spec$start_pars$garch_pars[[i]],
        dist_pars = NULL  ## Marginals in DCC are always normal
      )
    )
    
    ## Estimate using univariate function
    uni_result <- estimate_garch_weighted_univariate(series_residuals, w_target, uni_spec)
    garch_pars_list[[i]] <- uni_result$coefficients
    warnings_stage1 <- c(warnings_stage1, uni_result$warnings)
  }
  
  ## === STAGE 2: Estimate DCC Parameters ===
  
  ## Check if DCC parameters need estimation
  dcc_start_pars <- spec$start_pars$dcc_pars
  dist_start_pars <- spec$start_pars$dist_pars
  
  if (is.null(dcc_start_pars) || length(dcc_start_pars) == 0) {
    ## No DCC parameters to estimate (e.g., constant correlation)
    dcc_pars <- list()
  } else {
    ## Estimate DCC dynamics
    dcc_result <- estimate_dcc_parameters_weighted(
      residuals = residuals,
      weights = w_target,  ## Pass adjusted weights
      garch_pars = garch_pars_list,
      dcc_start_pars = dcc_start_pars,
      dist_start_pars = dist_start_pars,
      spec = spec
    )
    
    dcc_pars <- dcc_result$dcc_pars
    dist_pars <- dcc_result$dist_pars
    warnings_stage1 <- c(warnings_stage1, dcc_result$warnings)
  }
  
  ## === Combine Results ===
  return(list(
    coefficients = list(
      garch_pars = garch_pars_list,
      dcc_pars = dcc_pars,
      dist_pars = if(exists("dist_pars")) dist_pars else dist_start_pars
    ),
    warnings = warnings_stage1
  ))
}


#' @title Estimate DCC Correlation Dynamics (Weighted)
#' @keywords internal
estimate_dcc_parameters_weighted <- function(
    residuals, 
    weights, 
    garch_pars,
    dcc_start_pars, 
    dist_start_pars, 
    spec
) {
  
  k <- ncol(residuals)
  T_obs <- nrow(residuals)
  
  ## Verify dimensions match
  if (length(weights) != T_obs) {
    stop(sprintf("Dimension mismatch in estimate_dcc_parameters_weighted: residuals has %d rows but weights has %d elements", 
                 T_obs, length(weights)))
  }
  
  ## 1. Get standardized residuals using Stage 1 GARCH parameters
  std_residuals <- matrix(0, nrow = T_obs, ncol = k)
  for (i in 1:k) {
    series_residuals <- residuals[, i]
    series_spec <- spec$garch_spec_args$garch_model$univariate[[i]]
    
    ## Create spec with estimated GARCH parameters
    uni_spec_obj <- tsgarch::garch_modelspec(
      y = xts::xts(series_residuals, order.by = Sys.Date() - (T_obs:1)),
      model = series_spec$model,
      garch_order = series_spec$garch_order,
      distribution = series_spec$distribution
    )
    
    ## Set estimated parameters
    for (par_name in names(garch_pars[[i]])) {
      if (par_name %in% uni_spec_obj$parmatrix$parameter) {
        row_idx <- which(uni_spec_obj$parmatrix$parameter == par_name)
        uni_spec_obj$parmatrix[row_idx, "value"] <- garch_pars[[i]][[par_name]]
      }
    }
    
    ## Filter to get conditional volatilities
    uni_fit <- tsmethods::tsfilter(uni_spec_obj)
    std_residuals[, i] <- series_residuals / as.numeric(uni_fit$sigma)
  }
  
  ## 2. Combine DCC and distribution parameters for joint optimization
  all_stage2_pars <- c(dcc_start_pars, dist_start_pars)
  
  if (length(all_stage2_pars) == 0) {
    return(list(dcc_pars = list(), dist_pars = list(), warnings = list()))
  }
  
  ## 3. Define objective function for DCC parameters
  weighted_dcc_loglik <- function(params, std_resid, w, spec, k, debug = FALSE) {
    
    param_list <- as.list(params)
    names(param_list) <- names(all_stage2_pars)
    
    ## Separate DCC and distribution parameters
    dcc_param_names <- names(dcc_start_pars)
    dist_param_names <- names(dist_start_pars)
    
    dcc_params_current <- param_list[dcc_param_names]
    dist_params_current <- param_list[dist_param_names]
    
    ## Extract DCC order
    dcc_order <- spec$garch_spec_args$dcc_order
    
    ## Initialize matrices
    T_eff <- nrow(std_resid)
    
    ## Verify dimensions one more time inside objective function
    if (length(w) != T_eff) {
      stop(sprintf("Dimension mismatch: std_resid has %d rows but weights has %d elements", 
                   T_eff, length(w)))
    }
    
    ## ===== DIAGNOSTIC: Check weight distribution =====
    if (debug) {
      cat("\n=== OBJECTIVE FUNCTION DIAGNOSTIC ===\n")
      cat("Weight summary: min =", min(w), ", max =", max(w), ", mean =", mean(w), "\n")
      cat("Effective sample size:", sum(w)^2 / sum(w^2), "\n")
      cat("Params being evaluated: alpha =", dcc_params_current$alpha_1, 
          ", beta =", dcc_params_current$beta_1, "\n")
    }
    
    ## Weighted covariance using smoothed probabilities
    Qbar_result <- tryCatch({
      stats::cov.wt(std_resid, wt = w, method = "ML")$cov
    }, error = function(e) {
      if (debug) cat("ERROR in cov.wt:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(Qbar_result)) {
      if (debug) cat("Returning penalty: cov.wt failed\n")
      return(1e10)
    }
    
    Qbar <- Qbar_result
    
    ## Ensure positive definite
    eig <- eigen(Qbar, symmetric = TRUE)
    if (debug) {
      cat("Qbar eigenvalues:", eig$values, "\n")
    }
    
    if (any(eig$values < 1e-8)) {
      if (debug) cat("Qbar not PD, regularizing\n")
      ## Regularize: add small diagonal
      Qbar <- Qbar + diag(1e-6, k)
    }
    
    ## DCC recursion
    Q <- array(0, dim = c(k, k, T_eff))
    R <- array(0, dim = c(k, k, T_eff))
    
    ## Initialize Q and R
    Q[,,1] <- Qbar
    
    ## Properly initialize R[,,1] by standardizing Qbar
    Qbar_diag_inv_sqrt <- diag(1/sqrt(diag(Qbar)), k)
    R[,,1] <- Qbar_diag_inv_sqrt %*% Qbar %*% Qbar_diag_inv_sqrt
    
    alpha <- if(!is.null(dcc_params_current$alpha_1)) dcc_params_current$alpha_1 else 0
    beta <- if(!is.null(dcc_params_current$beta_1)) dcc_params_current$beta_1 else 0
    
    ## Check stationarity
    if ((alpha + beta) >= 1 || alpha < 0 || beta < 0) {
      return(1e10)
    }
    
    for (t in 2:T_eff) {
      z_lag <- std_resid[t-1, , drop = FALSE]
      Q[,,t] <- Qbar * (1 - alpha - beta) + 
        alpha * (t(z_lag) %*% z_lag) + 
        beta * Q[,,t-1]
      
      ## Check for invalid Q
      if (any(!is.finite(Q[,,t]))) {
        if (debug) cat("Non-finite Q at t =", t, "\n")
        return(1e10)
      }
      
      ## Standardize to get correlation
      Q_diag <- diag(Q[,,t])
      if (any(Q_diag <= 0)) {
        if (debug) cat("Non-positive diagonal in Q at t =", t, "\n")
        return(1e10)
      }
      
      Q_diag_inv_sqrt <- diag(1/sqrt(Q_diag), k)
      R[,,t] <- Q_diag_inv_sqrt %*% Q[,,t] %*% Q_diag_inv_sqrt
    }
    
    ## Calculate weighted log-likelihood
    ll_vec <- numeric(T_eff)
    n_bad <- 0
    
    if (spec$distribution == "mvn") {
      for (t in 1:T_eff) {
        R_t <- R[,,t]
        ## Ensure positive definite
        if (any(is.na(R_t)) || any(!is.finite(R_t))) {
          ll_vec[t] <- -1e10
          n_bad <- n_bad + 1
          if (debug && n_bad <= 3) cat("Bad R_t at t =", t, "(NA or non-finite)\n")
        } else {
          eig <- try(eigen(R_t, symmetric = TRUE, only.values = TRUE)$values, silent = TRUE)
          if (inherits(eig, "try-error") || any(eig <= 0)) {
            ll_vec[t] <- -1e10
            n_bad <- n_bad + 1
            if (debug && n_bad <= 3) {
              cat("Bad R_t at t =", t, "(not PD), eigenvalues:", 
                  if(inherits(eig, "try-error")) "ERROR" else paste(eig, collapse=", "), "\n")
            }
          } else {
            ll_vec[t] <- mvtnorm::dmvnorm(std_resid[t,], mean = rep(0, k), 
                                          sigma = R_t, log = TRUE)
          }
        }
      }
    } else if (spec$distribution == "mvt") {
      shape <- dist_params_current$shape
      if (shape <= 2) return(1e10)  ## Need finite variance
      
      for (t in 1:T_eff) {
        R_t <- R[,,t]
        if (any(is.na(R_t)) || any(!is.finite(R_t))) {
          ll_vec[t] <- -1e10
          n_bad <- n_bad + 1
        } else {
          eig <- try(eigen(R_t, symmetric = TRUE, only.values = TRUE)$values, silent = TRUE)
          if (inherits(eig, "try-error") || any(eig <= 0)) {
            ll_vec[t] <- -1e10
            n_bad <- n_bad + 1
          } else {
            ll_vec[t] <- mvtnorm::dmvt(std_resid[t,], delta = rep(0, k), 
                                       sigma = R_t, df = shape, log = TRUE)
          }
        }
      }
    }
    
    ll_vec[!is.finite(ll_vec)] <- -1e10
    
    nll <- -sum(w * ll_vec, na.rm = TRUE)
    
    if (debug) {
      cat("Negative log-likelihood:", nll, "\n")
      cat("Number of valid obs:", sum(ll_vec > -1e10), "/", T_eff, "\n")
    }
    
    return(nll)
  }
  
  ## 4. Get bounds
  lower_bounds <- numeric(length(all_stage2_pars))
  upper_bounds <- numeric(length(all_stage2_pars))
  
  for (i in seq_along(all_stage2_pars)) {
    par_name <- names(all_stage2_pars)[i]
    if (grepl("alpha", par_name)) {
      lower_bounds[i] <- 0.01  ## Enforce minimum alpha to avoid degeneracy
      upper_bounds[i] <- 0.99
    } else if (grepl("beta", par_name)) {
      lower_bounds[i] <- 1e-6
      upper_bounds[i] <- 0.99
    } else if (par_name == "shape") {
      lower_bounds[i] <- 2.1
      upper_bounds[i] <- 100
    }
  }
  
  ## 5. Optimize
  warnings_list <- list()
  
  ## First, do a test evaluation at starting parameters with debug=TRUE
  cat("\n=== TESTING OBJECTIVE AT START PARAMS ===\n")
  test_nll <- weighted_dcc_loglik(
    params = unlist(all_stage2_pars),
    std_resid = std_residuals,
    w = weights,
    spec = spec,
    k = k,
    debug = TRUE
  )
  cat("Start NLL:", test_nll, "\n")
  
  opt_result <- withCallingHandlers({
    stats::optim(
      par = unlist(all_stage2_pars),
      fn = weighted_dcc_loglik,
      lower = lower_bounds,
      upper = upper_bounds,
      method = "L-BFGS-B",
      std_resid = std_residuals,  ## Pass std_residuals
      w = weights,                ## Pass weights (now guaranteed to match)
      spec = spec,
      k = k,
      debug = FALSE  ## Turn off debug for actual optimization
    )
  }, warning = function(w) {
    warnings_list <<- c(warnings_list, list(w))
    invokeRestart("muffleWarning")
  })
  
  ## ======================= DIAGNOSTIC begin ==========================
  
  cat("\n=== DCC M-STEP DIAGNOSTIC ===\n")
  cat("Starting DCC params:\n")
  print(unlist(dcc_start_pars))
  cat("\nOptimized DCC params:\n")
  print(opt_result$par)
  cat("\nOptimization convergence:", opt_result$convergence, "\n")
  cat("Final objective value:", opt_result$value, "\n")
  
  ## Check if at boundary
  dcc_param_names <- names(dcc_start_pars)
  alpha_params <- opt_result$par[grepl("alpha", names(opt_result$par))]
  at_boundary <- any(alpha_params < 0.02)
  
  cat("At boundary:", at_boundary, "\n")
  
  if (at_boundary) {
    cat("\n*** WARNING: DCC alpha parameter at/near boundary ***\n")
    cat("Interpretation: This state may have CONSTANT (not dynamic) correlation.\n")
    cat("Consider:\n")
    cat("  1. This is a valid finding - the state genuinely lacks correlation dynamics\n")
    cat("  2. Using a constant correlation model for this state\n")
    cat("  3. Reducing the number of states if multiple states collapse\n")
  }
  
  ## ======================= DIAGNOSTIC end ==========================  
  
  ## 6. Extract results
  estimated_pars <- as.list(opt_result$par)
  names(estimated_pars) <- names(all_stage2_pars)
  
  dcc_param_names <- names(dcc_start_pars)
  dist_param_names <- names(dist_start_pars)
  
  dcc_pars_final <- estimated_pars[dcc_param_names]
  dist_pars_final <- estimated_pars[dist_param_names]
  
  return(list(
    dcc_pars = dcc_pars_final,
    dist_pars = dist_pars_final,
    warnings = warnings_list
  ))
}


#' @title GOGARCH Weighted Estimation (Placeholder)
#' @keywords internal
estimate_garch_weighted_gogarch <- function(residuals, weights, spec) {
  ## TODO: Implement GOGARCH estimation
  ## For now, return structure with empty parameters
  stop("GOGARCH estimation not yet implemented. Use DCC or Copula models.")
}


#' @title Copula Weighted Estimation (Placeholder)
#' @keywords internal
estimate_garch_weighted_copula <- function(residuals, weights, spec) {
  ## TODO: Implement Copula estimation
  ## For now, return structure with empty parameters
  stop("Copula estimation not yet implemented. Use DCC models for now.")
}


#' @title Perform the M-Step in Parallel (R Helper)
#' @description This function is called once per EM iteration from C++. It uses
#' the 'future' framework to estimate the parameters for all M states in parallel.
#' Handles both univariate and multivariate models with proper parameter structuring.
#' @param y The time series data.
#' @param weights The (T x M) matrix of smoothed probabilities from the E-step.
#' @param spec The full list of model specifications.
#' @param model_type "univariate" or "multivariate".
#' @return A list of length M containing the updated model fits for each state.
perform_m_step_parallel_r <- function(y, weights, spec, model_type) {
  
  ## Required packages for parallel workers
  required_packages <- c("data.table", "xts", "tsgarch", "tsmarch", 
                         "tsdistributions", "mvtnorm")
  
  ## Iterate over all states in parallel
  updated_fits <- future.apply::future_lapply(1:length(spec), function(j) {
    
    ## Explicitly load packages on each parallel worker
    library(data.table)
    library(xts)
    library(tsgarch)
    library(tsmarch)
    library(tsdistributions)
    library(mvtnorm)
    
    ## Extract state-specific data
    state_weights <- weights[, j]
    state_spec <- spec[[j]]
    
    ## === M-Step Stage 1: Update Mean Parameters ===
    new_mean_fit <- estimate_arma_weighted_r(
      y = y,
      weights = state_weights,
      spec = state_spec,
      model_type = model_type
    )
    
    ## === M-Step Stage 2: Update Variance Parameters ===
    new_variance_fit <- estimate_garch_weighted_r(
      residuals = new_mean_fit$residuals,
      weights = state_weights,
      spec = state_spec,
      model_type = model_type
    )
    
    ## === Structure the output based on model type ===
    if (model_type == "univariate") {
      ## For univariate: separate GARCH and distribution parameters
      all_params <- new_variance_fit$coefficients
      dist_param_names <- names(state_spec$start_pars$dist_pars)
      
      if (length(dist_param_names) > 0) {
        estimated_dist_pars <- all_params[dist_param_names]
        estimated_garch_pars <- all_params[!names(all_params) %in% dist_param_names]
      } else {
        estimated_dist_pars <- list()
        estimated_garch_pars <- all_params
      }
      
      return(list(
        arma_pars = new_mean_fit$coefficients,
        garch_pars = estimated_garch_pars,
        dist_pars = estimated_dist_pars
      ))
      
    } else {
      ## For multivariate: structure depends on model type
      variance_coeffs <- new_variance_fit$coefficients
      
      ## The coefficients already come structured from estimate_garch_weighted_multivariate
      ## For DCC: list(garch_pars = list(...), dcc_pars = list(...), dist_pars = list(...))
      ## For GOGARCH: list(garch_pars = list(...), rotation_pars = list(...))
      ## For Copula: list(garch_pars = list(...), copula_pars = list(...), dist_pars = list(...))
      
      ## Flatten DCC/Copula/GOGARCH-specific parameters into the top level
      ## This matches the structure expected by calculate_loglik_vector_r and C++ code
      
      result <- list(var_pars = new_mean_fit$coefficients)
      
      ## Add univariate GARCH parameters (always present)
      if (!is.null(variance_coeffs$garch_pars)) {
        result$garch_pars <- variance_coeffs$garch_pars
      }
      
      ## Add model-specific parameters (DCC, Copula, GOGARCH)
      if (!is.null(variance_coeffs$dcc_pars)) {
        ## Flatten DCC parameters to top level
        for (par_name in names(variance_coeffs$dcc_pars)) {
          result[[par_name]] <- variance_coeffs$dcc_pars[[par_name]]
        }
      }
      
      if (!is.null(variance_coeffs$copula_pars)) {
        ## Flatten Copula parameters to top level
        for (par_name in names(variance_coeffs$copula_pars)) {
          result[[par_name]] <- variance_coeffs$copula_pars[[par_name]]
        }
      }
      
      if (!is.null(variance_coeffs$rotation_pars)) {
        ## Flatten GOGARCH rotation parameters
        for (par_name in names(variance_coeffs$rotation_pars)) {
          result[[par_name]] <- variance_coeffs$rotation_pars[[par_name]]
        }
      }
      
      ## Add distribution parameters (if present)
      if (!is.null(variance_coeffs$dist_pars)) {
        result$dist_pars <- variance_coeffs$dist_pars
      }
      
      return(result)
    }
    
  }, future.seed = TRUE, future.packages = required_packages)
  
  return(updated_fits)
}


## Helper function for parameter counting
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}


## Helper function to extract GARCH NLL from univariate fitted objects
.get_garch_nll_vec_from_univariate <- function(univariate_list) {
  ## For Copula models, the univariate component is stored as fitted objects
  ## Extract per-observation negative log-likelihoods
  
  if (is.null(univariate_list) || length(univariate_list) == 0) {
    stop("No univariate GARCH models found")
  }
  
  n_obs <- length(univariate_list[[1]]$spec$target$y_orig)
  n_series <- length(univariate_list)
  
  garch_nll_matrix <- matrix(0, nrow = n_obs, ncol = n_series)
  
  for (i in seq_along(univariate_list)) {
    if (!is.null(univariate_list[[i]]$lik_vector)) {
      garch_nll_matrix[, i] <- univariate_list[[i]]$lik_vector
    } else {
      stop("lik_vector not found in univariate model ", i)
    }
  }
  
  return(rowSums(garch_nll_matrix))
}