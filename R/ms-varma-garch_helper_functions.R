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
    
    ## For univariate, we need to set the pars from the previous M-step before 
    ## filtering
    all_fixed_pars <- c(current_pars$garch_pars, current_pars$dist_pars)
    for (par_name in names(all_fixed_pars)) {
      if (par_name %in% garch_spec_obj$parmatrix$parameter) {
        garch_spec_obj$parmatrix[parameter == par_name, value := all_fixed_pars[[par_name]]]
      }
    }
  } else {
    ## --- Multivariate ---
    spec_fun_name <- spec$garch_spec_fun
    spec_fun <- get(spec_fun_name, asNamespace("tsmarch"))
    spec_args <- spec$garch_spec_args
    
    if (spec_fun_name %in% c("dcc_modelspec", "cgarch_modelspec")) {
      ## STEP 1: Create and parameterize each univariate GARCH model.
      univariate_models <- lapply(1:ncol(residuals_xts), function(i) {
        uni_spec_details <- spec_args$garch_model$univariate[[i]]
        uni_spec_obj <- tsgarch::garch_modelspec(
          y = residuals_xts[,i],
          model = uni_spec_details$model,
          garch_order = uni_spec_details$garch_order,
          distribution = uni_spec_details$distribution
        )
        
        ## Inject the current parameters for this series
        series_garch_pars <- current_pars$garch_pars[[i]]
        for (par_name in names(series_garch_pars)) {
          uni_spec_obj$parmatrix[parameter == par_name, value := series_garch_pars[[par_name]]]
        }
        
        ## We "estimate" with fixed parameters from the previous M-step to get 
        ## a valid fitted object
        suppressWarnings(estimate(uni_spec_obj, keep_tmb = TRUE))
      })
      
      ## STEP 2: Combine the univariate objects
      multi_estimate_object <- tsgarch::to_multi_estimate(univariate_models)
      names(multi_estimate_object) <- paste0("series_", 1:ncol(residuals_xts))
      
      ## STEP 3: Build the final DCC spec
      final_args <- c(
        list(object = multi_estimate_object), 
        spec_args[names(spec_args) != "garch_model"]
      )
      final_args$distribution <- spec$distribution
      
      garch_spec_obj <- do.call(spec_fun, final_args)
      
    } else {
      final_args <- c(list(y = residuals_xts), spec_args)
      garch_spec_obj <- do.call(spec_fun, final_args)
    }
  }
  return(garch_spec_obj)
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
          garch_spec_obj$parmatrix[parameter == par_name, value := other_pars[[par_name]]]
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
            pars, garch_spec_obj, 
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


#' @title Estimate GARCH Parameters (R Helper)
#' This function is the core of the M-step; it takes the weights (smoothed 
#' probabilities) from the E-step and finds the GARCH and distribution 
#' parameters that maximize the weighted log-likelihood.
#' @param residuals 
#' @param weights 
#' @param spec 
#' @param model_type 
#' @title Estimate GARCH and Distribution Parameters (R Helper) - GENERALIZED
#' @description Estimate both GARCH and distribution parameters (e.g., shape, 
#'   skew) simultaneously for the univariate case.
estimate_garch_weighted_r <- function(residuals, weights, spec, model_type = "univariate") {
  ## Combine GARCH and Distribution starting parameters into one vector
  start_pars <- c(spec$start_pars$garch_pars, spec$start_pars$dist_pars)
  
  ## Handle case where there are no GARCH parameters to estimate
  if (length(start_pars) == 0) {
    return(list(coefficients = list(), warnings = list()))
  }
  
  padding <- NROW(residuals) - length(weights)
  w_target <- if(padding > 0) weights[(padding+1):length(weights)] else weights
  
  ## Create a temporary current_pars for the initial spec object
  ## This uses the start_pars as initial values
  temp_current_pars <- list(
    garch_pars = spec$start_pars$garch_pars,
    dist_pars = spec$start_pars$dist_pars
  )
  
  temp_spec_obj <- create_garch_spec_object_r(
    residuals, 
    spec, 
    model_type, 
    temp_current_pars
  )
  
  ## Extract bounds for ALL parameters to be estimated (GARCH + dist)
  parmatrix <- temp_spec_obj$parmatrix
  pars_to_estimate <- names(start_pars)
  bounds_matrix <- parmatrix[parameter %in% pars_to_estimate]
  ## Ensure the bounds are in the same order as the parameters
  bounds_matrix <- bounds_matrix[match(pars_to_estimate, parameter),]
  lower_bounds <- bounds_matrix$lower
  upper_bounds <- bounds_matrix$upper
  
  ## Generalized objective function
  weighted_garch_loglik <- function(params, residuals_data, w, spec) {
    param_list <- as.list(params)
    names(param_list) <- names(start_pars)
    
    ## Separate GARCH and dist parameters for this iteration
    dist_param_names <- names(spec$start_pars$dist_pars)
    garch_param_names <- names(spec$start_pars$garch_pars)
    
    iter_current_pars <- list(
      garch_pars = param_list[garch_param_names],
      dist_pars = param_list[dist_param_names]
    )
    
    garch_spec_obj <- create_garch_spec_object_r(
      residuals_data, 
      spec, 
      "univariate",
      iter_current_pars  ## <-- Pass current iteration parameters
    )
    
    ## Set all parameter values (GARCH + dist) for this iteration of the optimizer
    for (par_name in names(param_list)) {
      if (par_name %in% garch_spec_obj$parmatrix$parameter) {
        garch_spec_obj$parmatrix[parameter == par_name, value := param_list[[par_name]]]
      }
    }
    
    ## Get the sigma path using the current parameter set
    fit <- try(tsmethods::tsfilter(garch_spec_obj), silent = TRUE)
    if (inherits(fit, "try-error")) return(1e10) # Penalize invalid parameter sets
    sig <- as.numeric(fit$sigma)
    
    ## Calculate log-likelihood using the specified distribution and its parameters
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
      dist_args <- c(
        list(
          x = residuals_data, 
          mu = 0, 
          sigma = sig, 
          log = TRUE
        ), 
        dist_pars_current
      )
    }
    
    ll_vector <- do.call(dist_fun, dist_args)
    ll_vector[!is.finite(ll_vector)] <- -1e10
    
    ## Return the weighted negative log-likelihood
    return(-sum(w * ll_vector, na.rm = TRUE))
  }
  
  ## Capture and summarize warnings from optim
  warnings_list <- list()
  opt_result <- withCallingHandlers({
    stats::optim(
      par = unlist(start_pars),
      fn = weighted_garch_loglik,
      lower = lower_bounds,
      upper = upper_bounds,
      method = "L-BFGS-B",
      # Pass additional arguments to the objective function
      residuals_data = residuals,
      w = w_target,
      spec = spec
    )
  }, warning = function(w) {
    warnings_list <<- c(warnings_list, list(w))
    invokeRestart("muffleWarning")
  })
  
  estimated_coeffs <- as.list(opt_result$par)
  names(estimated_coeffs) <- names(start_pars)
  
  return(list(coefficients = estimated_coeffs, warnings = warnings_list))
}


#' @title Perform the M-Step in Parallel (R Helper)
#' @description This function is called once per EM iteration from C++. It uses
#' the 'future' framework to estimate the parameters for all M states in parallel.
#' @param y The time series data.
#' @param weights The (T x M) matrix of smoothed probabilities from the E-step.
#' @param spec The full list of model specifications.
#' @param model_type "univariate" or "multivariate".
#' @return A list of length M containing the updated model fits for each state.
perform_m_step_parallel_r <- function(y, weights, spec, model_type) {
  
  ## --- Tell the parallel workers which packages to load ---
  ## The data.table syntax (e.g., parmatrix[parameter == ...]) and the
  ## tsgarch/tsmarch functions must be available on each worker.
  required_packages <- c("data.table", "xts", "tsgarch", "tsmarch", "tsdistributions")
  
  ## future_lapply will iterate from 1 to M (the number of states) in parallel.
  ## Each worker gets the index 'j' for the state it's responsible for.
  updated_fits <- future.apply::future_lapply(1:length(spec), function(j) {
    
    ## --- Explicitly load packages on each parallel worker ---
    ## This is a more robust approach than relying on future.packages, as it
    ## ensures the packages are fully attached, making special syntax like
    ## data.table's `[...]` available.
    library(data.table)
    library(xts)
    library(tsgarch)
    library(tsmarch)
    library(tsdistributions)
    
    ## Extract the data for this specific state
    state_weights <- weights[, j]
    state_spec <- spec[[j]]
    
    ## M-Step 1: Update Mean Parameters
    new_arma_fit <- estimate_arma_weighted_r(
      y = y,
      weights = state_weights,
      spec = state_spec,
      model_type = model_type
    )
    
    ## M-Step 2: Update Variance Parameters
    new_garch_fit <- estimate_garch_weighted_r(
      residuals = new_arma_fit$residuals,
      weights = state_weights,
      spec = state_spec,
      model_type = model_type
    )
    
    ## The returned coefficients are a flat list of garch + dist pars.
    ## We must separate them back into their structured groups.
    all_params <- new_garch_fit$coefficients
    dist_param_names <- names(state_spec$start_pars$dist_pars)
    
    if (length(dist_param_names) > 0) {
      ## Case: Distribution has shape/skew parameters
      estimated_dist_pars <- all_params[dist_param_names]
      estimated_garch_pars <- all_params[!names(all_params) %in% dist_param_names]
    } else {
      ## Case: Normal distribution (no dist_pars)
      estimated_dist_pars <- list() ## Return an empty list for type consistency
      estimated_garch_pars <- all_params
    }

    ## Return the fully structured list for the next EM iteration
    if (model_type == "univariate") {
      return(list(
        arma_pars = new_arma_fit$coefficients,
        garch_pars = estimated_garch_pars,
        dist_pars = estimated_dist_pars
      ))
    } else {
      ## Placeholder for the future multivariate refactoring
      return(list(
        var_pars = new_arma_fit$coefficients,
        garch_pars = estimated_garch_pars, ## Will need more complex logic later
        dist_pars = estimated_dist_pars   ## Will need more complex logic later
      ))
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