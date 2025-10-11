## === === === === === === === === === === === === === === === === === 
## Generalized R Helper Functions for MS-ARMA-GARCH Fitting
## === === === === === === === === === === === === === === === === === 
## These functions are designed to be called from the C++ EM orchestrator
## fit_ms_varma_garch_cpp(). They contain the full logic for handling general 
## ARMA(p,q) and VAR(p) models, and interface with the tsgarch/tsmarch packages.


#' @title Generate Correct Parameter Names for tsmarch
#' @description Translates a nested parameter list into the flat named list
#'              that tsmarch's parmatrix expects (e.g., "omega[1]").
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
#' @param residuals 
#' @param spec 
#' @param model_type 
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
        uni_spec_obj <- tsgarch::garch_modelspec(y = residuals_xts[,i],
                                                 model = uni_spec_details$model,
                                                 garch_order = uni_spec_details$garch_order,
                                                 distribution = uni_spec_details$distribution)
        
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
      final_args <- c(list(object = multi_estimate_object), spec_args[names(spec_args) != "garch_model"])
      final_args$distribution <- spec$distribution
      
      garch_spec_obj <- do.call(spec_fun, final_args)
      
    } else {
      final_args <- c(list(y = residuals_xts), spec_args)
      garch_spec_obj <- do.call(spec_fun, final_args)
    }
  }
  return(garch_spec_obj)
}

# create_garch_spec_object_r <- function(residuals, spec, model_type) {
#   residuals_xts <- xts::xts(residuals, order.by = Sys.Date() - (NROW(residuals):1))
#   
#   if (model_type == "univariate") {
#     ## --- Univariate ---
#     garch_spec_obj <- tsgarch::garch_modelspec(
#       y = residuals_xts,
#       model = spec$garch_model,
#       garch_order = spec$garch_order,
#       distribution = spec$distribution
#     )
#   } else {
#     ## --- Multivariate ---
#     spec_fun_name <- spec$garch_spec_fun
#     spec_fun <- get(spec_fun_name, asNamespace("tsmarch"))
#     spec_args <- spec$garch_spec_args
#     
#     if (spec_fun_name %in% c("dcc_modelspec", "cgarch_modelspec")) {
#       ## Fit the first-stage univariate models using their specified 
#       ## distributions
#       # univariate_models <- lapply(1:ncol(residuals_xts), function(i) {
#       #   uni_spec <- spec_args$garch_model$univariate[[i]]
#       #   suppressWarnings({
#       #     estimate(tsgarch::garch_modelspec(y = residuals_xts[,i], 
#       #                                       model = uni_spec$model, 
#       #                                       garch_order = uni_spec$garch_order,
#       #                                       distribution = "norm"),
#       #              keep_tmb = TRUE)
#       #   })
#       # })
#       # garch_model_est <- to_multi_estimate(univariate_models)
#       # names(garch_model_est) <- paste0("series_", 1:ncol(residuals_xts))
#       # 
#       # unnamed_arg <- garch_model_est
#       # named_args <- spec_args[names(spec_args) != "garch_model"]
#       # 
#       # ## Workaround for tsmarch bug where default distribution is invalid
#       # if (is.null(named_args$distribution)) {
#       #   ## If user did not specify a distribution, provide a valid default ("mvn")
#       #   ## to prevent tsmarch from using its own invalid default ("mvnorm").
#       #   named_args$distribution <- "mvn"
#       # }
#       # 
#       # final_args <- c(list(unnamed_arg), named_args)
#       # garch_spec_obj <- do.call(spec_fun, final_args)
#       # 
#       # # final_args <- c(list(object = garch_model_est), 
#       # #                 spec_args[names(spec_args) != "garch_model"])
#       # # garch_spec_obj <- do.call(spec_fun, final_args)
#       
#       ## STEP 1: Estimate each univariate GARCH model individually.
#       univariate_models <- lapply(1:ncol(residuals_xts), function(i) {
#         uni_spec_details <- spec_args$garch_model$univariate[[i]]
#         suppressWarnings(
#           estimate(
#             tsgarch::garch_modelspec(
#               y = residuals_xts[,i],
#               model = uni_spec_details$model,
#               garch_order = uni_spec_details$garch_order,
#               distribution = uni_spec_details$distribution
#             ),
#             keep_tmb = TRUE)
#         )
#       })
#       
#       ## STEP 2: Use 'to_multi_estimate' to create the special object class.
#       multi_estimate_object <- tsgarch::to_multi_estimate(univariate_models)
#       names(multi_estimate_object) <- paste0("series_", 1:ncol(residuals_xts))
#       
#       ## STEP 3: Call the spec function on the object.
#       final_args <- c(
#         list(object = multi_estimate_object),
#         spec_args[names(spec_args) != "garch_model"]
#       )
#       final_args$distribution <- spec$distribution
#       
#       garch_spec_obj <- do.call(spec_fun, final_args)
#       
#     } else {
#       final_args <- c(list(y = residuals_xts), spec_args)
#       garch_spec_obj <- do.call(spec_fun, final_args)
#     }
#   }
#   return(garch_spec_obj)
# }

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
  
  ## 2. Create the fully-parameterized Spec Object to get GARCH log-likelihood 
  ## from the residuals
  garch_spec_obj <- create_garch_spec_object_r(
    model_residuals, 
    spec, 
    model_type, 
    current_pars
  )

  if (model_type == "univariate") {
    ## --- UNIVARIATE ---
    
    ## Step 3a: Use tsfilter() to get the sigma path
    garch_model_fit <- tsmethods::tsfilter(garch_spec_obj)
    sig <- as.numeric(garch_model_fit$sigma)
    
    ## Step 3b: Call the correct density function based on the distribution
    dist_fun <- switch(spec$distribution,
                       "norm" = stats::dnorm, ## Normal
                       "snorm" = tsdistributions::dsnorm, ## Skew normal
                       "std"  = tsdistributions::dstd, ## Student t
                       "sstd" = tsdistributions::dsstd, ## Skew Student
                       "ged"  = tsdistributions::dged, ## Generalized error
                       "sged"  = tsdistributions::dsged, ## Skew generalized error
                       "ghyp"  = tsdistributions::dghyp, ## Generalized hyperbolic
                       "ghst"  = tsdistributions::dghst, ## Generalized hyperbolic skew Student
                       "jsu"  = tsdistributions::djsu, ## Johnson reparameterized SU
                       stop(paste("Unsupported univariate distribution:", spec$distribution))
    )
    
    ## Step 3c: Build the argument list specifically for the chosen distribution
    if (spec$distribution == "norm") {
      # Use the correct residuals from Step 1
      dist_args <- list(x = model_residuals, mean = 0, sd = sig, log = TRUE)
    } else {
      dist_args <- c(
        ## Use the correct residuals from Step 1
        list(x = model_residuals, mu = 0, sigma = sig, log = TRUE),
        current_pars$dist_pars
      )
    }
    
    ## The function will only use the arguments it needs (e.g., dnorm ignores 
    ## 'sigma' and 'shape')
    ll_vector <- do.call(dist_fun, dist_args)
    
  } else {
    ## --- MULTIVARIATE ---
    ## We need to get the estimated params from the previous M-step. To this end
    ## we force tsmarch to provide the loglikelihood given the already estimated
    ## parameters rather than estimating the parameters. This is done by
    ## setting (garch_spec_obj$parmatrix$estimate := 0), which will bypass
    ## the solver when we run estimate() and return a list containing the 
    ## loglikelihood and holding the params fixed.
    
    ## Get all DCC-layer and distribution parameters that need to be fixed.
    other_pars <- current_pars[!names(current_pars) %in% c("var_pars", "garch_pars")]
    
    ## Set the parameter values and mark them as NOT for estimation.
    for (par_name in names(other_pars)) {
      if (par_name %in% garch_spec_obj$parmatrix$parameter) {
        garch_spec_obj$parmatrix[parameter == par_name, value := other_pars[[par_name]]]
        garch_spec_obj$parmatrix[parameter == par_name, estimate := 0]
      }
    }
    
    ## TRICK: Re-enable estimation for ONE parameter, but with fixed bounds.
    ## This forces estimate() to run the full TMB path without actually changing 
    ## the value.
    if ("shape" %in% names(other_pars)) {
      shape_val <- other_pars$shape
      garch_spec_obj$parmatrix[parameter == "shape", estimate := 1]
      garch_spec_obj$parmatrix[parameter == "shape", lower := shape_val]
      garch_spec_obj$parmatrix[parameter == "shape", upper := shape_val]
    } else if ("alpha_1" %in% names(other_pars)) {
      ## Fallback if there is no shape parameter (e.g., mvn distribution)
      alpha_val <- other_pars$alpha_1
      garch_spec_obj$parmatrix[parameter == "alpha_1", estimate := 1]
      garch_spec_obj$parmatrix[parameter == "alpha_1", lower := alpha_val]
      garch_spec_obj$parmatrix[parameter == "alpha_1", upper := alpha_val]
    }
    
    ## 3d. Now, estimate() will build the TMB object and "converge" instantly.
    garch_model_fit <- suppressWarnings(estimate(garch_spec_obj, keep_tmb = TRUE))
    
    if (!is.null(garch_model_fit$TMB_OBJECT)) {
      ll_vector <- garch_model_fit$TMB_OBJECT$report()$ll_vector
    } else {
      stop("Could not extract log-likelihood vector: TMB object not found.")
    }
  }
  
  ## Sanitize and pad the vector before returning to C++
  ll_vector[!is.finite(ll_vector)] <- -1e10 # 0
  if (length(ll_vector) < NROW(y)) {
    padding <- NROW(y) - length(ll_vector)
    ll_vector <- c(rep(0, padding), ll_vector)
  }
  return(ll_vector)
}

# calculate_loglik_vector_r <- function(
#     y, 
#     current_pars, 
#     spec, 
#     model_type = "univariate"
#   ) {
#   
#   ## 1. Get Residuals from the Conditional Mean Model
#   if (model_type == "univariate") {
#     arma_pars <- current_pars$arma_pars
#     model_residuals <- stats::arima(
#       y, 
#       order = c(spec$arma_order[1], 0, spec$arma_order[2]),
#       fixed = arma_pars, 
#       include.mean = FALSE
#     )$residuals
#   } else {
#     var_order <- spec$var_order
#     k <- ncol(y)
#     T_obs <- nrow(y)
#     X_lagged <- matrix(1, nrow = T_obs - var_order, ncol = 1 + k * var_order)
#     for (i in 1:var_order) {
#       X_lagged[, (2 + (i - 1) * k):(1 + i*k)] <- y[(var_order - i + 1):(T_obs - i), ]
#     }
#     y_target <- y[(var_order + 1):T_obs, ]
#     beta_mat <- matrix(current_pars$var_pars, nrow = 1 + k * var_order, ncol = k)
#     model_residuals <- y_target - X_lagged %*% beta_mat
#   }
#   
#   ## 2. Get GARCH log-likelihood from the residuals
#   garch_spec_obj <- create_garch_spec_object_r(model_residuals, spec, model_type)
#   
#   # fixed_garch_pars <- current_pars$garch_pars
#   # for (par_name in names(fixed_garch_pars)) {
#   #   garch_spec_obj$parmatrix[parameter == par_name, value := fixed_garch_pars[[par_name]]]
#   #   garch_spec_obj$parmatrix[parameter == par_name, estimate := 0]
#   # }
#   
#   ## --- Combine GARCH and Distribution parameters for fixing ---
#   all_fixed_pars <- c(current_pars$garch_pars, current_pars$dist_pars)
#   for (par_name in names(all_fixed_pars)) {
#     if (par_name %in% garch_spec_obj$parmatrix$parameter) {
#       garch_spec_obj$parmatrix[parameter == par_name, value := all_fixed_pars[[par_name]]]
#       garch_spec_obj$parmatrix[parameter == par_name, estimate := 0]
#     }
#   }
#   
#   if (model_type == "univariate") {
#     ## --- UNIVARIATE ---
#     ## garch_model_fit <- tsmethods::tsfilter(garch_spec_obj)
#     ## res <- as.numeric(garch_model_fit$residuals)
#     ## sig <- as.numeric(garch_model_fit$sigma)
#     ## ll_vector <- dnorm(res, mean = 0, sd = sig, log = TRUE)
#     
#     ## Step 3a: Use tsfilter() to get the sigma path
#     garch_model_fit <- tsmethods::tsfilter(garch_spec_obj)
#     sig <- as.numeric(garch_model_fit$sigma)
#     
#     ## Step 3b: Call the correct density function based on the distribution
#     dist_fun <- switch(spec$distribution,
#                        "norm" = stats::dnorm, ## Normal
#                        "snorm" = tsdistributions::dsnorm, ## Skew normal
#                        "std"  = tsdistributions::dstd, ## Student t
#                        "sstd" = tsdistributions::dsstd, ## Skew Student
#                        "ged"  = tsdistributions::dged, ## Generalized error
#                        "sged"  = tsdistributions::dsged, ## Skew generalized error
#                        "ghyp"  = tsdistributions::dghyp, ## Generalized hyperbolic
#                        "ghst"  = tsdistributions::dghst, ## Generalized hyperbolic skew Student
#                        "jsu"  = tsdistributions::djsu, ## Johnson reparameterized SU
#                        stop(paste("Unsupported univariate distribution:", spec$distribution))
#     )
#     
#     ## Step 3c: Build the argument list specifically for the chosen distribution
#     if (spec$distribution == "norm") {
#       # Use the correct residuals from Step 1
#       dist_args <- list(x = model_residuals, mean = 0, sd = sig, log = TRUE)
#     } else {
#       dist_args <- c(
#         ## Use the correct residuals from Step 1
#         list(x = model_residuals, mu = 0, sigma = sig, log = TRUE),
#         current_pars$dist_pars
#       )
#     }
#     
#     ## The function will only use the arguments it needs (e.g., dnorm ignores 
#     ## 'sigma' and 'shape')
#     ll_vector <- do.call(dist_fun, dist_args)
#     
#   } else {
#     ## --- MULTIVARIATE ---
#     ## Combine all parameters into a single flat vector for fixing
#     ## Use the helper to generate the correctly named parameter list
#     all_fixed_pars <- generate_tsmarch_parnames(current_pars)
#     
#     for (par_name in names(all_fixed_pars)) {
#       if (par_name %in% garch_spec_obj$parmatrix$parameter) {
#         garch_spec_obj$parmatrix[parameter == par_name, value := all_fixed_pars[[par_name]]]
#         garch_spec_obj$parmatrix[parameter == par_name, estimate := 0]
#       }
#     }
#     
#     ## "Estimate" with all parameters fixed to get the TMB object
#     ## Suppress non-critical warnings during filtering
#     garch_model_fit <- suppressWarnings(estimate(garch_spec_obj, keep_tmb = TRUE))
#     
#     ## Extract the per-observation log-likelihood vector from the TMB report
#     if (!is.null(garch_model_fit$TMB_OBJECT)) {
#       ll_vector <- garch_model_fit$TMB_OBJECT$report()$ll_vector
#     } else {
#       stop("Could not extract log-likelihood vector: TMB object not found in multivariate fit.")
#     }
#     
#     # k <- ncol(model_residuals)
#     # T_res <- nrow(model_residuals)
#     # ll_vector <- numeric(T_res)
#     
#     # H_vectorized <- garch_model_fit$H
#     # for (t in 1:T_res) {
#     #   cov_mat <- matrix(0, k, k)
#     #   cov_mat[upper.tri(cov_mat, diag = TRUE)] <- H_vectorized[t, ]
#     #   cov_mat <- cov_mat + t(cov_mat)
#     #   diag(cov_mat) <- diag(cov_mat) / 2
#     #   
#     #   ll_vector[t] <- mvtnorm::dmvnorm(model_residuals[t,], mean = rep(0, k), 
#     #                                    sigma = cov_mat, log = TRUE)
#     # }
#   }
#   
#   ## Sanitize and pad the vector before returning to C++
#   ll_vector[!is.finite(ll_vector)] <- -1e10 # 0
#   if (length(ll_vector) < NROW(y)) {
#     padding <- NROW(y) - length(ll_vector)
#     ll_vector <- c(rep(0, padding), ll_vector)
#   }
#   return(ll_vector)
# }


# calculate_loglik_vector_r <- function(y, current_pars, spec, model_type = "univariate") {
#   
#   ## 1. Get Residuals from the Conditional Mean Model
#   if (model_type == "univariate") {
#     arma_pars <- current_pars$arma_pars
#     model_residuals <- stats::arima(y, order = c(spec$arma_order[1], 0, spec$arma_order[2]),
#                                     fixed = arma_pars, include.mean = FALSE)$residuals
#   } else {
#     var_order <- spec$var_order
#     k <- ncol(y)
#     T_obs <- nrow(y)
#     X_lagged <- matrix(1, nrow = T_obs - var_order, ncol = 1 + k * var_order)
#     for (i in 1:var_order) {
#       X_lagged[, (2 + (i - 1) * k):(1 + i * k)] <- y[(var_order - i + 1):(T_obs - i), ]
#     }
#     y_target <- y[(var_order + 1):T_obs, ]
#     beta_mat <- matrix(current_pars$var_pars, nrow = 1 + k * var_order, ncol = k)
#     model_residuals <- y_target - X_lagged %*% beta_mat
#   }
#   
#   ## 2. Get GARCH log-likelihood from the residuals
#   garch_spec_obj <- create_garch_spec_object_r(model_residuals, spec, model_type)
#   
#   # fixed_garch_pars <- current_pars$garch_pars
#   # for (par_name in names(fixed_garch_pars)) {
#   #   garch_spec_obj$parmatrix[parameter == par_name, value := fixed_garch_pars[[par_name]]]
#   #   garch_spec_obj$parmatrix[parameter == par_name, estimate := 0]
#   # }
#   
#   ## --- Combine GARCH and Distribution parameters for fixing ---
#   all_fixed_pars <- c(current_pars$garch_pars, current_pars$dist_pars)
#   for (par_name in names(all_fixed_pars)) {
#     if (par_name %in% garch_spec_obj$parmatrix$parameter) {
#       garch_spec_obj$parmatrix[parameter == par_name, value := all_fixed_pars[[par_name]]]
#       garch_spec_obj$parmatrix[parameter == par_name, estimate := 0]
#     }
#   }
#   
#   if (model_type == "univariate") {
#     ## garch_model_fit <- tsmethods::tsfilter(garch_spec_obj)
#     ## res <- as.numeric(garch_model_fit$residuals)
#     ## sig <- as.numeric(garch_model_fit$sigma)
#     ## ll_vector <- dnorm(res, mean = 0, sd = sig, log = TRUE)
#     
#     # Step 3a: Use tsfilter() to get the sigma path
#     garch_model_fit <- tsmethods::tsfilter(garch_spec_obj)
#     sig <- as.numeric(garch_model_fit$sigma)
#     
#     ## Step 3b: Call the correct density function based on the distribution
#     dist_fun <- switch(spec$distribution,
#      "norm" = stats::dnorm, ## Normal
#      "snorm" = tsdistributions::dsnorm, ## Skew normal
#      "std"  = tsdistributions::dstd, ## Student t
#      "sstd" = tsdistributions::dsstd, ## Skew Student
#      "ged"  = tsdistributions::dged, ## Generalized error
#      "sged"  = tsdistributions::dsged, ## Skew generalized error
#      "ghyp"  = tsdistributions::dghyp, ## Generalized hyperbolic
#      "ghst"  = tsdistributions::dghst, ## Generalized hyperbolic skew Student
#      "jsu"  = tsdistributions::djsu, ## Johnson reparameterized SU
#      stop(paste("Unsupported univariate distribution:", spec$distribution))
#     )
#     
#     ## Step 3c: Build the argument list specifically for the chosen distribution
#     if (spec$distribution == "norm") {
#       # Use the correct residuals from Step 1
#       dist_args <- list(x = model_residuals, mean = 0, sd = sig, log = TRUE)
#     } else {
#       dist_args <- c(
#         ## Use the correct residuals from Step 1
#         list(x = model_residuals, mu = 0, sigma = sig, log = TRUE),
#         current_pars$dist_pars
#       )
#     }
#     
#     ## The function will only use the arguments it needs (e.g., dnorm ignores 
#     ## 'sigma' and 'shape')
#     ll_vector <- do.call(dist_fun, dist_args)
#     
#   } else {
#     ## Suppress non-critical warnings during filtering
#     garch_model_fit <- suppressWarnings(estimate(garch_spec_obj))
#     k <- ncol(model_residuals)
#     T_res <- nrow(model_residuals)
#     ll_vector <- numeric(T_res)
#     
#     H_vectorized <- garch_model_fit$H
#     for (t in 1:T_res) {
#       cov_mat <- matrix(0, k, k)
#       cov_mat[upper.tri(cov_mat, diag = TRUE)] <- H_vectorized[t, ]
#       cov_mat <- cov_mat + t(cov_mat)
#       diag(cov_mat) <- diag(cov_mat) / 2
#       
#       ll_vector[t] <- mvtnorm::dmvnorm(model_residuals[t,], mean = rep(0, k), 
#                                        sigma = cov_mat, log = TRUE)
#     }
#   }
#   
#   ## Sanitize and pad the vector before returning to C++
#   ll_vector[!is.finite(ll_vector)] <- -1e10 # 0
#   if (length(ll_vector) < NROW(y)) {
#     padding <- NROW(y) - length(ll_vector)
#     ll_vector <- c(rep(0, padding), ll_vector)
#   }
#   return(ll_vector)
# }



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
  
  temp_spec_obj <- create_garch_spec_object_r(residuals, spec, model_type)
  
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
    
    garch_spec_obj <- create_garch_spec_object_r(residuals_data, spec, "univariate")
    
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
      dist_args <- c(list(x = residuals_data, mu = 0, sigma = sig, log = TRUE), dist_pars_current)
    }
    
    ll_vector <- do.call(dist_fun, dist_args)
    ll_vector[!is.finite(ll_vector)] <- -1e10
    
    ## Return the weighted negative log-likelihood
    return(-sum(w * ll_vector, na.rm = TRUE))
  }
  
  ## Capture and summarize warnings from optim
  warnings_list <- list()
  opt_result <- withCallingHandlers({
    stats::optim(par = unlist(start_pars),
                 fn = weighted_garch_loglik,
                 lower = lower_bounds,
                 upper = upper_bounds,
                 method = "L-BFGS-B",
                 # Pass additional arguments to the objective function
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

# estimate_garch_weighted_r <- function(residuals, weights, spec, model_type = "univariate") {
#   start_pars <- spec$start_pars$garch_pars
#   
#   ## Handle case where there are no GARCH parameters to estimate
#   if (length(start_pars) == 0) {
#     return(list(coefficients = list(), warnings = list()))
#   }
#   
#   padding <- NROW(residuals) - length(weights)
#   w_target <- if(padding > 0) weights[(padding+1):length(weights)] else weights
#   
#   temp_spec_obj <- create_garch_spec_object_r(residuals, spec, model_type)
#   
#   parmatrix <- temp_spec_obj$parmatrix
#   pars_to_estimate <- names(start_pars)
#   
#   bounds_matrix <- parmatrix[parameter %in% pars_to_estimate]
#   bounds_matrix <- bounds_matrix[match(pars_to_estimate, parameter),]
#   
#   lower_bounds <- bounds_matrix$lower
#   upper_bounds <- bounds_matrix$upper
#   
#   weighted_garch_loglik <- function(params, residuals_data, w, spec, model_type) {
#     param_list <- as.list(params)
#     names(param_list) <- names(spec$start_pars$garch_pars)
#     
#     garch_spec_obj <- create_garch_spec_object_r(residuals_data, spec, model_type)
#     
#     for (par_name in names(param_list)) {
#       garch_spec_obj$parmatrix[parameter == par_name, value := param_list[[par_name]]]
#     }
#     
#     fit <- try(estimate(garch_spec_obj), silent = TRUE)
#     if (inherits(fit, "try-error")) return(1e10)
#     
#     if ("TMB_OBJECT" %in% names(fit)) {
#       ll_vector <- fit$TMB_OBJECT$report()$ll_vector
#     } else {
#       res <- as.numeric(fit$residuals)
#       sig <- as.numeric(fit$sigma)
#       ll_vector <- dnorm(res, mean = 0, sd = sig, log = TRUE)
#     }
#     return(-sum(w * ll_vector, na.rm = TRUE))
#   }
#   
#   ## --- Capture and summarize warnings from optim ---
#   warnings_list <- list()
#   opt_result <- withCallingHandlers({
#     stats::optim(par = unlist(start_pars), 
#                  fn = weighted_garch_loglik,
#                  lower = lower_bounds,
#                  upper = upper_bounds,
#                  method = "L-BFGS-B",
#                  residuals_data = residuals, 
#                  w = w_target, 
#                  spec = spec, 
#                  model_type = model_type)
#   }, warning = function(w) {
#     warnings_list <<- c(warnings_list, list(w))
#     invokeRestart("muffleWarning")
#   })
#   
#   estimated_coeffs <- as.list(opt_result$par)
#   names(estimated_coeffs) <- names(start_pars)
#   
#   return(list(coefficients = estimated_coeffs, warnings = warnings_list))
# }



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
    
  #   ## Combine the results for this state into a single list
  #   if (model_type == "univariate") {
  #     return(list(
  #       arma_pars = new_arma_fit$coefficients,
  #       garch_pars = new_garch_fit$coefficients
  #     ))
  #   } else {
  #     return(list(
  #       var_pars = new_arma_fit$coefficients,
  #       garch_pars = new_garch_fit$coefficients
  #     ))
  #   }
  # }, future.seed = TRUE, future.packages = required_packages) ## <-- Pass the packages here
    
    ## The returned coefficients are a flat list of garch + dist pars.
    ## We must separate them back into their structured groups.
    all_params <- new_garch_fit$coefficients
    dist_param_names <- names(state_spec$start_pars$dist_pars)
    
    if (length(dist_param_names) > 0) {
      # Case: Distribution has shape/skew parameters
      estimated_dist_pars <- all_params[dist_param_names]
      estimated_garch_pars <- all_params[!names(all_params) %in% dist_param_names]
    } else {
      # Case: Normal distribution (no dist_pars)
      estimated_dist_pars <- list() # Return an empty list for type consistency
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



# perform_m_step_parallel_r <- function(y, weights, spec, model_type) {
#   
#   ## --- Tell the parallel workers which packages to load ---
#   ## The data.table syntax (e.g., parmatrix[parameter == ...]) and the
#   ## tsgarch/tsmarch functions must be available on each worker.
#   required_packages <- c("data.table", "xts", "tsgarch", "tsmarch", "tsdistributions")
#   
#   ## future_lapply will iterate from 1 to M (the number of states) in parallel.
#   ## Each worker gets the index 'j' for the state it's responsible for.
#   updated_fits <- future.apply::future_lapply(1:length(spec), function(j) {
#     
#     ## --- Explicitly load packages on each parallel worker ---
#     ## This is a more robust approach than relying on future.packages, as it
#     ## ensures the packages are fully attached, making special syntax like
#     ## data.table's `[...]` available.
#     library(data.table)
#     library(tsgarch)
#     library(tsmarch)
#     library(xts)
#     
#     ## Extract the data for this specific state
#     state_weights <- weights[, j]
#     state_spec <- spec[[j]]
#     
#     ## M-Step 1: Update Mean Parameters
#     new_arma_fit <- estimate_arma_weighted_r(
#       y = y,
#       weights = state_weights,
#       spec = state_spec,
#       model_type = model_type
#     )
#     
#     ## M-Step 2: Update Variance Parameters
#     new_garch_fit <- estimate_garch_weighted_r(
#       residuals = new_arma_fit$residuals,
#       weights = state_weights,
#       spec = state_spec,
#       model_type = model_type
#     )
#     
#     ## Combine the results for this state into a single list
#     if (model_type == "univariate") {
#       return(list(
#         arma_pars = new_arma_fit$coefficients,
#         garch_pars = new_garch_fit$coefficients
#       ))
#     } else {
#       return(list(
#         var_pars = new_arma_fit$coefficients,
#         garch_pars = new_garch_fit$coefficients
#       ))
#     }
#   }, future.seed = TRUE, future.packages = required_packages) ## <-- Pass the packages here
#   
#   return(updated_fits)
# }
# 
# 
# ## Helper function for parameter counting
# `%||%` <- function(a, b) {
#   if (is.null(a)) b else a
# }



