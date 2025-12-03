## === === === === === === === === === === === === === === === === === 
## Generalized R Helper Functions for MS-ARMA-GARCH Fitting
## === === === === === === === === === === === === === === === === === 
## These functions are designed to be called from the C++ EM orchestrator
## fit_ms_varma_garch_cpp(). They contain the full logic for handling general 
## ARMA(p,q) and VAR(p) models, and interface with the tsgarch/tsmarch packages.


#' @title Generate Correct Parameter Names for tsmarch
#' @description Translates a nested parameter list into the flat named list
#'              that tsmarch's parmatrix expects (e.g., "omega\[1\]").
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
    current_pars,
    diagnostics = NULL,
    iteration = NULL,
    state = NULL,
    verbose = FALSE
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
        
        ## Update the parameter values (but keep estimate = 1)
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
            
            ## Update only the value, keep estimate flag unchanged
            estimated_obj$parmatrix$value[row_match] <- value
          }
        }
        
        ## Recompute sigma using TMB with updated parameters.
        ## Get parameters in the order TMB expects (estimate == 1)
        estimate <- NULL
        pars_for_tmb <- estimated_obj$parmatrix[estimate == 1]$value
        
        if (verbose) {
          cat("\n=== TMB SIGMA RECOMPUTATION (Series", i, ") ===\n")
          cat("Parameters for TMB (should reflect NEW values):\n")
          print(pars_for_tmb)
          cat("Omega we WANT:", series_garch_pars$omega, "\n")
          cat("Omega in pars_for_tmb:", pars_for_tmb[1], "\n")
        }
        
        if (!is.null(estimated_obj$tmb) && length(pars_for_tmb) == length(estimated_obj$tmb$par)) {
          if (verbose) cat("TMB exists and lengths match. Calling tmb$report()...\n")
          
          ## Get sigma before recomputation
          sigma_before <- as.numeric(estimated_obj$sigma)
          if (verbose) {
            cat("Sigma BEFORE tmb$report (first 5):", head(sigma_before, 5), "\n")
            cat("Mean sigma BEFORE:", mean(sigma_before), "\n")
          }
          
          ## Evaluate TMB at new parameters
          tmb_report <- estimated_obj$tmb$report(pars_for_tmb)
          
          ## Update sigma (strip initialization period to match the data length)
          maxpq <- max(estimated_obj$spec$model$order)
          
          if (verbose) {
            cat("\n=== XTS CREATION DIAGNOSTIC ===\n")
            cat("Length of tmb_report$sigma:", length(tmb_report$sigma), "\n")
            cat("maxpq:", maxpq, "\n")
            cat("Length of residuals_xts[,", i, "]:", length(residuals_xts[,i]), "\n")
            cat("Class of estimated_obj$sigma:", class(estimated_obj$sigma), "\n")  ## KEEP THIS
          }
          
          ## Create proper dates for the sigma vector
          ## The sigma from TMB matches the full residuals length
          sigma_dates <- zoo::index(residuals_xts[,i])
          
          if (maxpq > 0) {
            ## Strip initialization period from both sigma and dates
            new_sigma <- tmb_report$sigma[-(1:maxpq)]
            if (verbose) {
              cat("After stripping", maxpq, "init periods:\n")
              cat("  Length of new_sigma:", length(new_sigma), "\n")
              cat("  Length of sigma_dates:", length(sigma_dates), "\n")
            }
          } else {
            new_sigma <- tmb_report$sigma
            if (verbose) {
              cat("No initialization period to strip\n")
              cat("  Length of new_sigma:", length(new_sigma), "\n")
              cat("  Length of sigma_dates:", length(sigma_dates), "\n")
            }
          }
          
          if (verbose) cat("=== END XTS DIAGNOSTIC ===\n\n")
          
          ## Create xts object with proper dates
          estimated_obj$sigma <- xts::xts(new_sigma, order.by = sigma_dates)
          
          ## Get sigma AFTER recomputation
          sigma_after <- as.numeric(estimated_obj$sigma)
          sigma_changed <- !identical(sigma_before, sigma_after)
          
          if (verbose) {
            cat("Sigma AFTER tmb$report (first 5):", head(sigma_after, 5), "\n")
            cat("Mean sigma AFTER:", mean(sigma_after), "\n")
            cat("Did sigma change?", sigma_changed, "\n")
          }
          
          ## DIAGNOSTIC: Collect sigma evolution
          if (!is.null(diagnostics) && !is.null(iteration) && !is.null(state)) {
            sigma_summary <- list(
              mean = mean(sigma_after, na.rm = TRUE),
              sd = sd(sigma_after, na.rm = TRUE),
              min = min(sigma_after, na.rm = TRUE),
              max = max(sigma_after, na.rm = TRUE),
              first_5 = head(sigma_after, 5),
              last_5 = tail(sigma_after, 5),
              changed = sigma_changed
            )
            
            diagnostics <- add_sigma_evolution(diagnostics, iteration, state, i, sigma_summary)
          }
          
          ## Also update standardized residuals
          ## Make sure we're working with the numeric values
          if (inherits(estimated_obj$residuals, "xts")) {
            resid_numeric <- as.numeric(estimated_obj$residuals)
          } else {
            resid_numeric <- estimated_obj$residuals
          }
          
          sigma_numeric <- as.numeric(estimated_obj$sigma)
          
          ## Create xts object for standardized residuals
          estimated_obj$std_residuals <- xts::xts(
            resid_numeric / sigma_numeric,
            order.by = sigma_dates
          )
        } else {
          if (verbose) {
            cat("WARNING: TMB recomputation SKIPPED!\n")
            cat("  tmb is NULL?", is.null(estimated_obj$tmb), "\n")
            if (!is.null(estimated_obj$tmb)) {
              cat("  Length mismatch: pars_for_tmb has", length(pars_for_tmb), 
                  "but tmb$par has", length(estimated_obj$tmb$par), "\n")
            }
          }
          
          ## Log warning in diagnostics
          if (!is.null(diagnostics) && !is.null(iteration) && !is.null(state)) {
            diagnostics <- add_diagnostic_warning(
              diagnostics,
              iteration,
              "tmb_skip",
              paste0("TMB recomputation skipped for state ", state, " series ", i),
              list(tmb_null = is.null(estimated_obj$tmb))
            )
          }
        }
        
        if (verbose) cat("=== END TMB RECOMPUTATION ===\n\n")
        
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
calculate_loglik_vector_r <- function(
    y, 
    current_pars, 
    spec, 
    model_type = "univariate",
    diagnostics = NULL,
    iteration = NULL,
    state = NULL,
    verbose = FALSE
) {
  
  ## 1. Get Residuals from Conditional Mean Model
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
  
  ## 2. Create GARCH spec object
  garch_spec_obj <- create_garch_spec_object_r(
    model_residuals, 
    spec, 
    model_type, 
    current_pars,
    diagnostics = diagnostics,
    iteration = iteration,
    state = state,
    verbose = verbose
  )
  
  if (model_type == "univariate") {
    ## === UNIVARIATE ===
    garch_model_fit <- tsmethods::tsfilter(garch_spec_obj)
    sig <- as.numeric(garch_model_fit$sigma)
    
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
    ## === MULTIVARIATE ===
    spec_fun_name <- spec$garch_spec_fun
    
    ## ==================== DIAGNOSTIC begin ====================
    if (verbose && model_type == "multivariate" && spec_fun_name == "dcc_modelspec") {
      cat("\n=== DIAGNOSTIC: Sigma values in univariate models ===\n")
      for (i in 1:ncol(model_residuals)) {
        uni_model <- garch_spec_obj$univariate[[paste0("series_", i)]]
        cat("Series", i, ":\n")
        cat("  omega from parmatrix:", uni_model$parmatrix[parameter == "omega"]$value, "\n")
        cat("  omega from current_pars:", current_pars$garch_pars[[i]]$omega, "\n")
        cat("  First 5 sigma values:", head(as.numeric(uni_model$sigma), 5), "\n")
        cat("  Mean sigma:", mean(as.numeric(uni_model$sigma)), "\n")
      }
      cat("=== END DIAGNOSTIC ===\n\n")
    }
    ## ==================== DIAGNOSTIC end ====================
    
    if (spec_fun_name %in% c("dcc_modelspec", "cgarch_modelspec")) {
      ## ==== DCC or Copula-GARCH ====
      
      ## Check if this state has constant or dynamic correlation
      has_dcc_params <- any(grepl("^(alpha|beta)_[0-9]+$", names(current_pars)))
      is_constant_corr <- !has_dcc_params || 
        (!is.null(current_pars$correlation_type) && current_pars$correlation_type == "constant")
      
      if (is_constant_corr) {
        ## === CONSTANT CORRELATION ===
        if(verbose) {
          cat("Using constant correlation (DCC parameters absent or at boundary)\n")
        }
        
        ## Inject distribution parameters only
        if (!is.null(current_pars$dist_pars)) {
          for (par_name in names(current_pars$dist_pars)) {
            if (par_name %in% garch_spec_obj$parmatrix$parameter) {
              row_idx <- which(garch_spec_obj$parmatrix$parameter == par_name)
              garch_spec_obj$parmatrix[row_idx, "value"] <- current_pars$dist_pars[[par_name]]
            }
          }
        }
        
        ## Extract parameter vector (for dist params only, if any)
        estimate <- NULL
        pars <- garch_spec_obj$parmatrix[estimate == 1]$value
        
        ## Call constant correlation function
        if (spec_fun_name == "dcc_modelspec") {
          if (spec$distribution == "mvn" && length(pars) == 0) {
            total_nll_vec <- tsmarch:::.dcc_constant_values(
              NULL, 
              garch_spec_obj, 
              type = "ll_vec"
            )
          } else {
            total_nll_vec <- tsmarch:::.dcc_constant_values(
              pars, 
              garch_spec_obj, 
              type = "ll_vec"
            )
          }
          ll_vector <- -as.vector(total_nll_vec)
          
        } else if (spec_fun_name == "cgarch_modelspec") {
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
        
      } else {
        ## === DYNAMIC CORRELATION ===
        
        ## Inject DCC and distribution parameters
        other_pars <- current_pars[!names(current_pars) %in% c("var_pars", "garch_pars", "correlation_type")]
        
        for (par_name in names(other_pars)) {
          if (par_name %in% garch_spec_obj$parmatrix$parameter) {
            row_idx <- which(garch_spec_obj$parmatrix$parameter == par_name)
            garch_spec_obj$parmatrix[row_idx, "value"] <- other_pars[[par_name]]
          }
        }
        
        estimate <- NULL
        pars <- garch_spec_obj$parmatrix[estimate == 1]$value
        
        is_dynamic <- !is.null(garch_spec_obj$dynamics) && 
          garch_spec_obj$dynamics$model %in% c("dcc", "adcc")
        
        if (spec_fun_name == "dcc_modelspec") {
          if (is_dynamic) {
            
            ## ==================== DIAGNOSTIC begin ====================
            if (verbose) {
              cat("\n=== BEFORE calling .dcc_dynamic_values() ===\n")
              cat("Parameters being passed:\n")
              print(pars)
              cat("\nSigma from univariate[['series_1']]:\n")
              print(head(as.numeric(garch_spec_obj$univariate[['series_1']]$sigma), 10))
              cat("\n")
            }
            ## ==================== DIAGNOSTIC end ====================
            
            total_nll_vec <- tsmarch:::.dcc_dynamic_values(
              pars, 
              garch_spec_obj, 
              type = "ll_vec"
            )
            
            ## ==================== DIAGNOSTIC begin ====================
            if (verbose) {
              cat("=== AFTER .dcc_dynamic_values() ===\n")
              cat("Returned NLL vector (first 10):\n")
              print(head(total_nll_vec, 10))
              cat("\n")
            }
            ## ==================== DIAGNOSTIC end ====================
            
            dcc_order <- garch_spec_obj$dynamics$order
            maxpq <- max(dcc_order)
            if (maxpq > 0) {
              total_nll_vec <- total_nll_vec[-(1:maxpq), , drop = TRUE]
            } else {
              total_nll_vec <- as.vector(total_nll_vec)
            }
            ll_vector <- -total_nll_vec
          } else {
            ## Shouldn't reach here, but fallback to constant
            if (spec$distribution == "mvn" && length(pars) == 0) {
              total_nll_vec <- tsmarch:::.dcc_constant_values(NULL, garch_spec_obj, type = "ll_vec")
            } else {
              total_nll_vec <- tsmarch:::.dcc_constant_values(pars, garch_spec_obj, type = "ll_vec")
            }
            ll_vector <- -as.vector(total_nll_vec)
          }
          
        } else if (spec_fun_name == "cgarch_modelspec") {
          if (is_dynamic) {
            copula_nll_vec <- tsmarch:::.copula_dynamic_values(
              pars, 
              garch_spec_obj, 
              type = "ll_vec"
            )
            
            dcc_order <- garch_spec_obj$dynamics$order
            maxpq <- max(dcc_order)
            if (maxpq > 0) {
              copula_nll_vec <- copula_nll_vec[-(1:maxpq), , drop = TRUE]
            } else {
              copula_nll_vec <- as.vector(copula_nll_vec)
            }
            
            garch_nll_vec <- .get_garch_nll_vec_from_univariate(garch_spec_obj$univariate)
            total_nll_vec <- garch_nll_vec + copula_nll_vec
            ll_vector <- -total_nll_vec
          } else {
            if (garch_spec_obj$copula == "mvn" && length(pars) == 0) {
              copula_nll_vec <- tsmarch:::.copula_constant_values(NULL, garch_spec_obj, type = "ll_vec")
            } else {
              copula_nll_vec <- tsmarch:::.copula_constant_values(pars, garch_spec_obj, type = "ll_vec")
            }
            copula_nll_vec <- as.vector(copula_nll_vec)
            garch_nll_vec <- .get_garch_nll_vec_from_univariate(garch_spec_obj$univariate)
            total_nll_vec <- garch_nll_vec + copula_nll_vec
            ll_vector <- -total_nll_vec
          }
        }
      }
      
    } else if (spec_fun_name == "gogarch_modelspec") {
      ## ==== GOGARCH ====
      garch_model_fit <- suppressWarnings(estimate(garch_spec_obj, trace = FALSE))
      ll_vector <- compute_loglik_fixed(
        object = garch_model_fit,
        params = list(),
        ll_vec = TRUE
      )
    } else {
      stop("Unsupported multivariate model type: ", spec_fun_name)
    }
  }
  
  ## Sanitize and pad
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
estimate_garch_weighted_r <- function(
    residuals, 
    weights, 
    spec, 
    model_type = "univariate",
    diagnostics = NULL, 
    iteration = NULL, 
    state = NULL,
    verbose = FALSE,
    dcc_threshold = 0.02,
    dcc_criterion = "bic",
    force_constant = FALSE
  ) {
  
  if (model_type == "univariate") {
    ## === UNIVARIATE CASE ===
    return(
      estimate_garch_weighted_univariate(
        residuals = residuals, 
        weights = weights, 
        spec = spec,
        verbose = verbose
      )
    )
  } else {
    ## === MULTIVARIATE CASE ===
    return(estimate_garch_weighted_multivariate(
      residuals = residuals, 
      weights = weights, 
      spec = spec,
      diagnostics = diagnostics, 
      iteration = iteration, 
      state = state,
      verbose = verbose,
      dcc_threshold = dcc_threshold,
      dcc_criterion = dcc_criterion,
      force_constant = force_constant
    ))
  }
}


#' @title Univariate GARCH + Distribution Parameter Estimation
#' @keywords internal
estimate_garch_weighted_univariate <- function(
    residuals, 
    weights, 
    spec,
    verbose
  ) {
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
    
    if (inherits(fit, "try-error")) {
      if (verbose) cat("*** tsfilter failed in weighted_garch_loglik ***\n")
      return(1e10)
    }
    
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
estimate_garch_weighted_multivariate <- function(
    residuals,
    weights,
    spec,
    diagnostics = NULL,
    iteration = NULL,
    state = NULL,
    verbose = FALSE,
    dcc_threshold = 0.02,
    dcc_criterion = "bic",
    force_constant = FALSE
  ) {

  ## Determine model type
  model_type <- spec$garch_spec_fun

  if (model_type %in% c("dcc_modelspec", "cgarch_modelspec")) {
    return(estimate_garch_weighted_dcc(
      residuals = residuals,
      weights = weights,
      spec = spec,
      diagnostics = diagnostics,
      iteration = iteration,
      state = state,
      verbose = verbose,
      dcc_threshold = 0.02,
      dcc_criterion = "bic",
      force_constant = FALSE
    ))
  } else if (model_type == "copula_modelspec") {
    return(estimate_garch_weighted_copula(
      residuals = residuals,
      weights = weights,
      spec = spec,
      diagnostics = diagnostics,
      verbose = verbose
    ))
  } else if (model_type == "gogarch_modelspec") {
    return(estimate_garch_weighted_gogarch(
      residuals = residuals,
      weights = weights,
      spec = spec,
      diagnostics = diagnostics,
      verbose = verbose
    ))
  } else {
    stop(paste("Unsupported multivariate model type:", model_type))
  }
}


#' #' @title DCC Weighted Estimation
#' #' @keywords internal
#' estimate_garch_weighted_dcc <- function(
#'     residuals, 
#'     weights, 
#'     spec, 
#'     diagnostics = NULL, 
#'     iteration = NULL, 
#'     state = NULL,
#'     verbose = FALSE
#'   ) {
#'   
#'   k <- ncol(residuals)
#'   T_obs <- nrow(residuals)
#'   
#'   ## Adjust weights to match residuals length
#'   if (length(weights) > T_obs) {
#'     n_to_remove <- length(weights) - T_obs
#'     w_target <- weights[(n_to_remove + 1):length(weights)]
#'   } else if (length(weights) < T_obs) {
#'     stop("Weights vector is shorter than residuals - this should not happen")
#'   } else {
#'     w_target <- weights
#'   }
#'   
#'   ## === STAGE 1: Estimate Univariate GARCH for Each Series ===
#'   garch_pars_list <- list()
#'   warnings_stage1 <- list()
#'   
#'   for (i in 1:k) {
#'     series_residuals <- residuals[, i]
#'     series_spec <- spec$garch_spec_args$garch_model$univariate[[i]]
#'     
#'     uni_spec <- list(
#'       garch_model = series_spec$model,
#'       garch_order = series_spec$garch_order,
#'       distribution = series_spec$distribution,
#'       start_pars = list(
#'         garch_pars = spec$start_pars$garch_pars[[i]],
#'         dist_pars = NULL
#'       )
#'     )
#'     
#'     uni_result <- estimate_garch_weighted_univariate(
#'       residuals = series_residuals, 
#'       weights = w_target, 
#'       spec = uni_spec,
#'       verbose = verbose
#'     )
#'     garch_pars_list[[i]] <- uni_result$coefficients
#'     warnings_stage1 <- c(warnings_stage1, uni_result$warnings)
#'   }
#'   
#'   ## === STAGE 2: Estimate DCC Parameters (with degeneracy detection) ===
#'   
#'   dcc_start_pars <- spec$start_pars$dcc_pars
#'   dist_start_pars <- spec$start_pars$dist_pars
#'   
#'   ## Check if this state was constant in the previous iteration
#'   was_constant_before <- !is.null(spec$start_pars$correlation_type) && 
#'     spec$start_pars$correlation_type == "constant"
#'   
#'   if (was_constant_before) {
#'     ## Once constant, stay constant - don't try to re-estimate DCC
#'     if (verbose) {
#'       cat(sprintf("\n=== State %d: Maintaining CONSTANT correlation (was constant before) ===\n",
#'                   state))
#'     }
#'     
#'     return(list(
#'       coefficients = list(
#'         garch_pars = garch_pars_list,
#'         dcc_pars = list(),
#'         dist_pars = dist_start_pars,
#'         correlation_type = "constant",
#'         degeneracy_reason = "maintained_from_previous_iteration"
#'       ),
#'       warnings = warnings_stage1,
#'       diagnostics = diagnostics
#'     ))
#'   }
#' 
#'   if (is.null(dcc_start_pars) || length(dcc_start_pars) == 0) {
#'     ## No DCC parameters - already constant correlation
#'     return(list(
#'       coefficients = list(
#'         garch_pars = garch_pars_list,
#'         dcc_pars = list(),
#'         dist_pars = dist_start_pars,
#'         correlation_type = "constant"
#'       ),
#'       warnings = warnings_stage1,
#'       diagnostics = diagnostics
#'     ))
#'   }
#'   
#'   ## Attempt DCC estimation
#'   dcc_result <- estimate_dcc_parameters_weighted(
#'     residuals = residuals,
#'     weights = w_target,
#'     garch_pars = garch_pars_list,
#'     dcc_start_pars = dcc_start_pars,
#'     dist_start_pars = dist_start_pars,
#'     spec = spec,
#'     diagnostics = diagnostics,
#'     iteration = iteration,
#'     state = state,
#'     verbose = verbose
#'   )
#'   
#'   ## Update diagnostics from the result
#'   if (!is.null(dcc_result$diagnostics)) {
#'     diagnostics <- dcc_result$diagnostics
#'   }
#'   
#'   ## Check for degeneracy
#'   # alpha_params <- dcc_result$dcc_pars[grepl("alpha", names(dcc_result$dcc_pars))]
#'   
#'   alpha_params <- dcc_result$dcc_pars[grepl("alpha", names(dcc_result$dcc_pars))]
#'   beta_params <- dcc_result$dcc_pars[grepl("beta", names(dcc_result$dcc_pars))]
#' 
#'   # Degeneracy if no alpha params OR alpha near zero
#'   # is_degenerate <- any(unlist(alpha_params) < 0.015)  # Near lower bound
#'   is_degenerate <- is.null(alpha_params) || 
#'     length(alpha_params) == 0 || 
#'     any(unlist(alpha_params) <= 0.0101)
#'   
#'   # if (is_degenerate) {
#'   #   
#'   #   # DEBUG: Add this
#'   #   if (!is.null(iteration)) {
#'   #     cat(sprintf("  ACTION: Switching to constant correlation\n"))
#'   #   }
#'   #   
#'   #   if(verbose) {
#'   #     cat("\n=== DCC DEGENERACY DETECTED ===\n")
#'   #     cat("Alpha parameter(s) at lower bound:", unlist(alpha_params), "\n")
#'   #     cat("Switching to CONSTANT CORRELATION model for this state.\n")
#'   #     cat("This is appropriate when correlation lacks meaningful dynamics.\n\n")
#'   #   }
#'   #   
#'   #   ## Return constant correlation specification
#'   #   return(list(
#'   #     coefficients = list(
#'   #       garch_pars = garch_pars_list,
#'   #       dcc_pars = list(),  ## Empty - signals constant correlation
#'   #       dist_pars = dcc_result$dist_pars,
#'   #       correlation_type = "constant"
#'   #     ),
#'   #     warnings = c(warnings_stage1, dcc_result$warnings),
#'   #     diagnostics = diagnostics  # DEBUG: Is this being returned?
#'   #   ))
#'   # }
#'   
#'   if (is_degenerate) {
#'     # Determine reason
#'     degeneracy_reason <- if (is.null(alpha_params) || length(alpha_params) == 0) {
#'       "no DCC parameters returned"
#'     } else {
#'       sprintf("alpha near zero (min: %.6f)", min(unlist(alpha_params)))
#'     }
#'     
#'     if(verbose || !is.null(diagnostics)) {
#'       cat(sprintf("\n=== DCC DEGENERACY DETECTED (State %d, Iteration %d) ===\n",
#'                   state, iteration))
#'       cat("Reason:", degeneracy_reason, "\n")
#'       if (!is.null(alpha_params) && length(alpha_params) > 0) {
#'         cat("Alpha:", paste(round(unlist(alpha_params), 6), collapse = ", "), "\n")
#'         cat("Beta:", paste(round(unlist(beta_params), 6), collapse = ", "), "\n")
#'       }
#'       cat("Action: Switching to CONSTANT CORRELATION model.\n")
#'       cat("Interpretation: Correlations are stable over time in this regime.\n\n")
#'     }
#'     
#'     ## RECORD BOUNDARY EVENT
#'     if (!is.null(diagnostics) && !is.null(iteration) && !is.null(state)) {
#'       if (!is.null(alpha_params) && length(alpha_params) > 0) {
#'         # Alpha params exist but are small
#'         for (pname in names(alpha_params)) {
#'           diagnostics <- add_boundary_event(
#'             diagnostics,
#'             iteration = iteration,
#'             state = state,
#'             parameter_name = pname,
#'             value = alpha_params[[pname]],
#'             boundary_type = "lower",
#'             action_taken = paste0("constant_correlation_fallback: ", degeneracy_reason)
#'           )
#'         }
#'       } else {
#'         # No alpha params returned at all
#'         diagnostics <- add_boundary_event(
#'           diagnostics,
#'           iteration = iteration,
#'           state = state,
#'           parameter_name = "alpha_1",
#'           value = NA,
#'           boundary_type = "lower",
#'           action_taken = paste0("constant_correlation_fallback: ", degeneracy_reason)
#'         )
#'       }
#'     }
#'     
#'     ## Return constant correlation specification
#'     return(list(
#'       coefficients = list(
#'         garch_pars = garch_pars_list,
#'         dcc_pars = list(),
#'         dist_pars = dcc_result$dist_pars,
#'         correlation_type = "constant",
#'         degeneracy_reason = degeneracy_reason
#'       ),
#'       warnings = c(warnings_stage1, dcc_result$warnings),
#'       diagnostics = diagnostics
#'     ))
#'   }
#'   
#'   ## DCC dynamics are meaningful
#'   return(list(
#'     coefficients = list(
#'       garch_pars = garch_pars_list,
#'       dcc_pars = dcc_result$dcc_pars,
#'       dist_pars = dcc_result$dist_pars,
#'       correlation_type = "dynamic"
#'     ),
#'     warnings = c(warnings_stage1, dcc_result$warnings),
#'     diagnostics = diagnostics
#'   ))
#' }


#' @title DCC Weighted Estimation
#' @keywords internal
estimate_garch_weighted_dcc <- function(
    residuals, 
    weights, 
    spec, 
    diagnostics = NULL, 
    iteration = NULL, 
    state = NULL,
    verbose = FALSE,
    dcc_threshold = 0.02,
    dcc_criterion = "bic",
    force_constant = FALSE
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
  
  ## === STAGE 1: Estimate Univariate GARCH for Each Series ===
  garch_pars_list <- list()
  warnings_stage1 <- list()
  
  for (i in 1:k) {
    series_residuals <- residuals[, i]
    series_spec <- spec$garch_spec_args$garch_model$univariate[[i]]
    
    uni_spec <- list(
      garch_model = series_spec$model,
      garch_order = series_spec$garch_order,
      distribution = series_spec$distribution,
      start_pars = list(
        garch_pars = spec$start_pars$garch_pars[[i]],
        dist_pars = NULL
      )
    )
    
    uni_result <- estimate_garch_weighted_univariate(
      residuals = series_residuals, 
      weights = w_target, 
      spec = uni_spec,
      verbose = verbose
    )
    garch_pars_list[[i]] <- uni_result$coefficients
    warnings_stage1 <- c(warnings_stage1, uni_result$warnings)
  }
  
  ## === STAGE 2: DCC Parameters ===
  
  dcc_start_pars <- spec$start_pars$dcc_pars
  dist_start_pars <- spec$start_pars$dist_pars
  
  ## === SHORT-CIRCUIT: If forced to constant, skip DCC estimation ===
  if (force_constant) {
    if (verbose) {
      cat(sprintf("\n=== State %d: Forced to CONSTANT correlation ===\n", state))
    }
    
    return(list(
      coefficients = list(
        garch_pars = garch_pars_list,
        dcc_pars = list(),
        dist_pars = dist_start_pars,
        correlation_type = "constant",
        degeneracy_reason = "forced_constant"
      ),
      warnings = warnings_stage1,
      diagnostics = diagnostics
    ))
  }
  
  ## === HANDLE MISSING DCC PARAMETERS ===
  if (is.null(dcc_start_pars) || length(dcc_start_pars) == 0) {
    return(list(
      coefficients = list(
        garch_pars = garch_pars_list,
        dcc_pars = list(),
        dist_pars = dist_start_pars,
        correlation_type = "constant"
      ),
      warnings = warnings_stage1,
      diagnostics = diagnostics
    ))
  }
  
  ## === ESTIMATE DYNAMIC CORRELATION MODEL ===
  dcc_result <- estimate_dcc_parameters_weighted(
    residuals = residuals,
    weights = w_target,
    garch_pars = garch_pars_list,
    dcc_start_pars = dcc_start_pars,
    dist_start_pars = dist_start_pars,
    spec = spec,
    diagnostics = diagnostics,
    iteration = iteration,
    state = state,
    verbose = verbose
  )
  
  ## Update diagnostics from DCC estimation
  if (!is.null(dcc_result$diagnostics)) {
    diagnostics <- dcc_result$diagnostics
  }
  
  ## Update diagnostics from DCC estimation
  if (!is.null(dcc_result$diagnostics)) {
    diagnostics <- dcc_result$diagnostics
  }
  
  ## GET THE WEIGHTED LL
  ll_dynamic <- dcc_result$weighted_ll
  
  ## Check if it exists
  if (is.null(ll_dynamic)) {
    stop("Internal error: weighted_ll not returned from DCC estimation")
  }
  
  ## === DECIDE: Dynamic vs Constant ===
  
  alpha_params <- dcc_result$dcc_pars[grepl("alpha", names(dcc_result$dcc_pars))]
  beta_params <- dcc_result$dcc_pars[grepl("beta", names(dcc_result$dcc_pars))]
  
  ## Check if parameters are near boundary
  near_boundary <- !is.null(alpha_params) && 
    length(alpha_params) > 0 && 
    any(unlist(alpha_params) < dcc_threshold)
  
  if (!near_boundary) {
    ## Parameters clearly away from boundary - use dynamic
    return(list(
      coefficients = list(
        garch_pars = garch_pars_list,
        dcc_pars = dcc_result$dcc_pars,
        dist_pars = dcc_result$dist_pars,
        correlation_type = "dynamic"
      ),
      warnings = c(warnings_stage1, dcc_result$warnings),
      diagnostics = diagnostics
    ))
  }
  
  ## === NEAR BOUNDARY: Apply selection criterion ===
  
  if (dcc_criterion == "threshold") {
    ## Simple threshold - switch to constant
    degeneracy_reason <- sprintf("threshold (alpha=%.6f < %.3f)", 
                                 min(unlist(alpha_params)), dcc_threshold)
    
    if (verbose) {
      cat(sprintf("\n=== DCC DEGENERACY (State %d, Iter %d) ===\n", state, iteration))
      cat("Criterion: threshold\n")
      cat("Alpha:", paste(round(unlist(alpha_params), 6), collapse = ", "), "\n")
      cat("Decision: CONSTANT\n\n")
    }
    
    ## Record boundary event
    if (!is.null(diagnostics) && !is.null(iteration) && !is.null(state)) {
      for (pname in names(alpha_params)) {
        diagnostics <- add_boundary_event(
          diagnostics,
          iteration = iteration,
          state = state,
          parameter_name = pname,
          value = alpha_params[[pname]],
          boundary_type = "lower",
          action_taken = paste0("constant_correlation: ", degeneracy_reason)
        )
      }
    }
    
    return(list(
      coefficients = list(
        garch_pars = garch_pars_list,
        dcc_pars = list(),
        dist_pars = dcc_result$dist_pars,
        correlation_type = "constant",
        degeneracy_reason = degeneracy_reason
      ),
      warnings = c(warnings_stage1, dcc_result$warnings),
      diagnostics = diagnostics
    ))
  }
  
  ## === AIC/BIC CRITERION: Compare models ===
  
  ## Compute constant correlation likelihood
  const_result <- compute_constant_correlation_likelihood(
    residuals = residuals,
    weights = w_target,
    garch_pars = garch_pars_list,
    dist_pars = dist_start_pars,
    spec = spec
  )
  
  ## Get dynamic likelihood (already computed)
  ll_dynamic <- dcc_result$weighted_ll
  
  # ll_constant <- const_result$weighted_ll
  # 
  # ## Count parameters
  # k_dyn <- length(alpha_params) + length(beta_params)
  # k_const <- 0
  # 
  # ## Compute effective sample size for BIC
  # n_eff <- sum(w_target^2) / sum(w_target)
  # 
  # ## Compute information criterion
  # if (dcc_criterion == "aic") {
  #   ic_dynamic <- -2 * ll_dynamic + 2 * k_dyn
  #   ic_constant <- -2 * ll_constant + 2 * k_const
  # } else {  # bic
  #   ic_dynamic <- -2 * ll_dynamic + log(n_eff) * k_dyn
  #   ic_constant <- -2 * ll_constant + log(n_eff) * k_const
  # }
  # 
  # use_constant <- (ic_constant < ic_dynamic)
  # 
  # if (verbose) {
  #   cat(sprintf("\n=== Model Selection (State %d, Iter %d) ===\n", state, iteration))
  #   cat(sprintf("Criterion: %s\n", toupper(dcc_criterion)))
  #   cat(sprintf("Dynamic:  LL=%.4f, k=%d, %s=%.4f\n", 
  #               ll_dynamic, k_dyn, toupper(dcc_criterion), ic_dynamic))
  #   cat(sprintf("Constant: LL=%.4f, k=%d, %s=%.4f\n", 
  #               ll_constant, k_const, toupper(dcc_criterion), ic_constant))
  #   cat(sprintf("Decision: %s\n\n", ifelse(use_constant, "CONSTANT", "DYNAMIC")))
  # }
  
  ## Safety check - if weighted_ll wasn't returned, we have a problem
  if (is.null(ll_dynamic) || !is.finite(ll_dynamic)) {
    warning("Could not retrieve dynamic likelihood for comparison. Falling back to threshold criterion.")
    
    # Use threshold fallback
    degeneracy_reason <- "comparison_failed_fallback_to_threshold"
    use_constant <- TRUE
    
  } else {
    ## Compute constant correlation likelihood
    const_result <- compute_constant_correlation_likelihood(
      residuals = residuals,
      weights = w_target,
      garch_pars = garch_pars_list,
      dist_pars = dist_start_pars,
      spec = spec
    )
    
    ll_constant <- const_result$weighted_ll
    
    ## Count parameters
    k_dyn <- length(alpha_params) + length(beta_params)
    k_const <- 0
    
    ## Compute effective sample size for BIC
    n_eff <- sum(w_target^2) / sum(w_target)
    
    ## Compute information criterion
    if (dcc_criterion == "aic") {
      ic_dynamic <- -2 * ll_dynamic + 2 * k_dyn
      ic_constant <- -2 * ll_constant + 2 * k_const
    } else {  # bic
      ic_dynamic <- -2 * ll_dynamic + log(n_eff) * k_dyn
      ic_constant <- -2 * ll_constant + log(n_eff) * k_const
    }
    
    use_constant <- (ic_constant < ic_dynamic)
    
    if (verbose) {
      cat(sprintf("\n=== Model Selection (State %d, Iter %d) ===\n", state, iteration))
      cat(sprintf("Criterion: %s (n_eff=%.1f)\n", toupper(dcc_criterion), n_eff))
      cat(sprintf("Dynamic:  LL=%.4f, k=%d, %s=%.4f\n", 
                  ll_dynamic, k_dyn, toupper(dcc_criterion), ic_dynamic))
      cat(sprintf("Constant: LL=%.4f, k=%d, %s=%.4f\n", 
                  ll_constant, k_const, toupper(dcc_criterion), ic_constant))
      cat(sprintf("Decision: %s\n\n", ifelse(use_constant, "CONSTANT", "DYNAMIC")))
    }
    
    degeneracy_reason <- sprintf("%s selection (IC_const=%.2f < IC_dyn=%.2f)",
                                 toupper(dcc_criterion), ic_constant, ic_dynamic)
  }
  
  if (use_constant) {
    if (!is.null(diagnostics) && !is.null(iteration) && !is.null(state)) {
      for (pname in names(alpha_params)) {
        diagnostics <- add_boundary_event(
          diagnostics,
          iteration = iteration,
          state = state,
          parameter_name = pname,
          value = alpha_params[[pname]],
          boundary_type = "lower",
          action_taken = paste0("constant_correlation: ", degeneracy_reason)
        )
      }
    }
    
    return(list(
      coefficients = list(
        garch_pars = garch_pars_list,
        dcc_pars = list(),
        dist_pars = const_result$dist_pars,
        correlation_type = "constant",
        degeneracy_reason = degeneracy_reason
      ),
      warnings = c(warnings_stage1, const_result$warnings),
      diagnostics = diagnostics
    ))
  } else {
    ## Keep dynamic
    return(list(
      coefficients = list(
        garch_pars = garch_pars_list,
        dcc_pars = dcc_result$dcc_pars,
        dist_pars = dcc_result$dist_pars,
        correlation_type = "dynamic"
      ),
      warnings = c(warnings_stage1, dcc_result$warnings),
      diagnostics = diagnostics
    ))
  }
}


#' @title Estimate DCC Correlation Dynamics (Weighted)
#' @keywords internal
estimate_dcc_parameters_weighted <- function(
    residuals, 
    weights, 
    garch_pars,
    dcc_start_pars, 
    dist_start_pars, 
    spec,
    diagnostics = NULL,
    iteration = NULL,
    state = NULL,
    verbose = FALSE
) {
  
  k <- ncol(residuals)
  T_obs <- nrow(residuals)
  
  if (length(weights) != T_obs) {
    stop(sprintf("Dimension mismatch: residuals has %d rows but weights has %d elements", 
                 T_obs, length(weights)))
  }
  
  ## 1. Get standardized residuals using Stage 1 GARCH parameters
  std_residuals <- matrix(0, nrow = T_obs, ncol = k)
  
  for (i in 1:k) {
    series_residuals <- residuals[, i]
    series_spec <- spec$garch_spec_args$garch_model$univariate[[i]]
    
    uni_spec_obj <- tsgarch::garch_modelspec(
      y = xts::xts(series_residuals, order.by = Sys.Date() - (T_obs:1)),
      model = series_spec$model,
      garch_order = series_spec$garch_order,
      distribution = series_spec$distribution
    )
    
    for (par_name in names(garch_pars[[i]])) {
      if (par_name %in% uni_spec_obj$parmatrix$parameter) {
        row_idx <- which(uni_spec_obj$parmatrix$parameter == par_name)
        uni_spec_obj$parmatrix[row_idx, "value"] <- garch_pars[[i]][[par_name]]
      }
    }
    
    uni_fit <- tsmethods::tsfilter(uni_spec_obj)
    sigma_vec <- as.numeric(uni_fit$sigma)
    std_residuals[, i] <- series_residuals / sigma_vec
    
    ## DIAGNOSTIC: Collect sigma evolution
    if (!is.null(diagnostics) && !is.null(iteration)) {
      sigma_summary <- list(
        mean = mean(sigma_vec, na.rm = TRUE),
        sd = sd(sigma_vec, na.rm = TRUE),
        min = min(sigma_vec, na.rm = TRUE),
        max = max(sigma_vec, na.rm = TRUE),
        first_5 = head(sigma_vec, 5),
        last_5 = tail(sigma_vec, 5)
      )
      
      diagnostics <- add_sigma_evolution(diagnostics, iteration, state, i, sigma_summary)
    }
  }
  
  ## 2. Combine DCC and distribution parameters
  all_stage2_pars <- c(dcc_start_pars, dist_start_pars)
  
  if (length(all_stage2_pars) == 0) {
    return(list(
      dcc_pars = list(), 
      dist_pars = list(), 
      warnings = list(),
      diagnostics = diagnostics
    ))
  }
  
  ## 3. Define objective function
  ## NOTE: This function modifies 'diagnostics' in parent scope using <<-
  weighted_dcc_loglik <- function(
      params, 
      std_resid, 
      w, 
      spec, 
      k, 
      verbose = FALSE
    ) {
    
    param_list <- as.list(params)
    names(param_list) <- names(all_stage2_pars)
    
    dcc_param_names <- names(dcc_start_pars)
    dist_param_names <- names(dist_start_pars)
    
    dcc_params_current <- param_list[dcc_param_names]
    dist_params_current <- param_list[dist_param_names]
    
    T_eff <- nrow(std_resid)
    
    if (length(w) != T_eff) {
      stop(sprintf("Dimension mismatch: std_resid has %d rows but weights has %d elements", 
                   T_eff, length(w)))
    }
    
    if (verbose) {
      cat("\n=== OBJECTIVE FUNCTION DIAGNOSTIC ===\n")
      cat("Weight summary: min =", min(w), ", max =", max(w), ", mean =", mean(w), "\n")
      cat("Effective sample size:", sum(w)^2 / sum(w^2), "\n")
      cat("Params: alpha =", dcc_params_current$alpha_1, ", beta =", dcc_params_current$beta_1, "\n")
    }
    
    ## Weighted covariance
    Qbar_result <- tryCatch({
      stats::cov.wt(std_resid, wt = w, method = "ML")$cov
    }, error = function(e) {
      if (verbose) cat("ERROR in cov.wt:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(Qbar_result)) {
      if (verbose) cat("Returning penalty: cov.wt failed\n")
      
      ## Update diagnostics in parent scope
      if (!is.null(diagnostics) && !is.null(iteration) && !is.null(state)) {
        diagnostics <<- add_diagnostic_warning(
          diagnostics,
          iteration,
          "dcc_penalty",
          "Weighted covariance calculation failed",
          list(reason = "cov.wt returned NULL")
        )
      }
      
      return(1e10)
    }
    
    Qbar <- Qbar_result
    
    ## Ensure positive definite
    eig <- eigen(Qbar, symmetric = TRUE)
    if (verbose) cat("Qbar eigenvalues:", eig$values, "\n")
    
    if (any(eig$values < 1e-8)) {
      if (verbose) cat("Qbar not PD, regularizing\n")
      Qbar <- Qbar + diag(1e-6, k)
    }
    
    ## DCC recursion
    Q <- array(0, dim = c(k, k, T_eff))
    R <- array(0, dim = c(k, k, T_eff))
    
    ## Initialize Q and R properly
    Q[,,1] <- Qbar
    Qbar_diag_inv_sqrt <- diag(1/sqrt(diag(Qbar)), k)
    R[,,1] <- Qbar_diag_inv_sqrt %*% Qbar %*% Qbar_diag_inv_sqrt
    
    alpha <- if(!is.null(dcc_params_current$alpha_1)) dcc_params_current$alpha_1 else 0
    beta <- if(!is.null(dcc_params_current$beta_1)) dcc_params_current$beta_1 else 0
    
    ## Check stationarity
    if ((alpha + beta) >= 1 || alpha < 0 || beta < 0) {
      reason <- if ((alpha + beta) >= 1) {
        "non-stationary (alpha + beta >= 1)"
      } else if (alpha < 0) {
        "negative alpha"
      } else {
        "negative beta"
      }
      
      if (verbose) {
        cat("*** Returning penalty 1e10: DCC parameters", reason, "***\n")
        cat("  alpha =", alpha, ", beta =", beta, ", sum =", alpha + beta, "\n")
      }
      
      if (!is.null(diagnostics) && !is.null(iteration) && !is.null(state)) {
        diagnostics <- add_diagnostic_warning(
          diagnostics,
          iteration,
          "dcc_penalty",
          paste0("DCC parameters ", reason),
          list(alpha = alpha, beta = beta, sum = alpha + beta)
        )
      }
      
      return(1e10)
    }
    
    for (t in 2:T_eff) {
      z_lag <- std_resid[t-1, , drop = FALSE]
      Q[,,t] <- Qbar * (1 - alpha - beta) + 
        alpha * (t(z_lag) %*% z_lag) + 
        beta * Q[,,t-1]
      
      if (any(!is.finite(Q[,,t]))) {
        if (verbose) {
          cat("*** Returning penalty 1e10: Non-finite Q at t =", t, "***\n")
        }
        
        ## Update diagnostics in parent scope
        if (!is.null(diagnostics) && !is.null(iteration) && !is.null(state)) {
          diagnostics <<- add_diagnostic_warning(
            diagnostics,
            iteration,
            "dcc_penalty",
            paste0("Non-finite Q matrix at t=", t),
            list(t = t)
          )
        }
        
        return(1e10)
      }
      
      Q_diag <- diag(Q[,,t])
      
      if (any(Q_diag <= 0)) {
        if (verbose) {
          cat("*** Returning penalty 1e10: Non-positive diagonal in Q at t =", t, "***\n")
        }
        
        if (!is.null(diagnostics) && !is.null(iteration) && !is.null(state)) {
          diagnostics <<- add_diagnostic_warning(
            diagnostics,
            iteration,
            "dcc_penalty",
            paste0("Non-positive Q diagonal at t=", t),
            list(t = t, Q_diag_min = min(Q_diag))
          )
        }
        
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
        
        if (any(is.na(R_t)) || any(!is.finite(R_t))) {
          ll_vec[t] <- -1e10
          n_bad <- n_bad + 1
          
          if (verbose && n_bad <= 3) {
            cat("*** Setting ll_vec[", t, "] = -1e10: Bad R_t (NA or non-finite) ***\n")
          }
        } else {
          eig <- try(eigen(R_t, symmetric = TRUE, only.values = TRUE)$values, silent = TRUE)
          
          if (inherits(eig, "try-error") || any(eig <= 0)) {
            ll_vec[t] <- -1e10
            n_bad <- n_bad + 1
            
            if (verbose && n_bad <= 3) {
              cat("Bad R_t at t =", t, "(not PD), eigenvalues:", 
                  if(inherits(eig, "try-error")) "ERROR" else paste(eig, collapse=", "), "\n")
            }
          } else {
            ll_vec[t] <- mvtnorm::dmvnorm(
              std_resid[t,], 
              mean = rep(0, k),
              sigma = R_t, 
              log = TRUE
            )
          }
        }
      }
      
      if (verbose && n_bad > 0) {
        cat("Total observations with bad R_t:", n_bad, "/", T_eff, "\n")
      }
      
      if (!is.null(diagnostics) && !is.null(iteration) && !is.null(state) && n_bad > 0) {
        diagnostics <- add_diagnostic_warning(
          diagnostics, iteration, "dcc_bad_correlation",
          paste0("Bad correlation matrices in ", n_bad, " observations"),
          list(n_bad = n_bad, T_eff = T_eff)
        )
      }
    } else if (spec$distribution == "mvt") {
      shape <- dist_params_current$shape
      
      if (shape <= 2) {
        if (verbose) cat("*** Returning penalty 1e10: shape <= 2 ***\n")
        
        if (!is.null(diagnostics) && !is.null(iteration) && !is.null(state)) {
          diagnostics <<- add_diagnostic_warning(
            diagnostics,
            iteration,
            "dcc_penalty",
            "MVT shape parameter too small (must be > 2)",
            list(shape = shape)
          )
        }
        
        return(1e10)
      }
      
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
    
    if (verbose) {
      cat("Negative log-likelihood:", nll, "\n")
      cat("Number of valid obs:", sum(ll_vec > -1e10), "/", T_eff, "\n")
    }
    
    ## Log bad correlations to diagnostics
    if (!is.null(diagnostics) && !is.null(iteration) && !is.null(state) && n_bad > 0) {
      diagnostics <<- add_diagnostic_warning(
        diagnostics,
        iteration,
        "dcc_bad_correlation",
        paste0("Bad correlation matrices in ", n_bad, " observations"),
        list(n_bad = n_bad, T_eff = T_eff)
      )
    }
    
    return(nll)
  }
  
  ## 4. Get bounds
  lower_bounds <- numeric(length(all_stage2_pars))
  upper_bounds <- numeric(length(all_stage2_pars))
  
  for (i in seq_along(all_stage2_pars)) {
    par_name <- names(all_stage2_pars)[i]
    if (grepl("alpha", par_name)) {
      lower_bounds[i] <- 0.01
      upper_bounds[i] <- 0.99
    } else if (grepl("beta", par_name)) {
      lower_bounds[i] <- 1e-6 #0.01 
      upper_bounds[i] <- 0.99
    } else if (par_name == "shape") {
      lower_bounds[i] <- 2.1
      upper_bounds[i] <- 100
    }
  }
  
  ## 5. Optimize
  warnings_list <- list()
  
  if (verbose) {
    cat("\n=== TESTING OBJECTIVE AT START PARAMS ===\n")
    test_nll <- weighted_dcc_loglik(
      params = unlist(all_stage2_pars),
      std_resid = std_residuals,
      w = weights,
      spec = spec,
      k = k,
      verbose = TRUE
    )
    cat("Start NLL:", test_nll, "\n")
  }
  
  opt_result <- withCallingHandlers({
    stats::optim(
      par = unlist(all_stage2_pars),
      fn = weighted_dcc_loglik,
      lower = lower_bounds,
      upper = upper_bounds,
      method = "L-BFGS-B",
      std_resid = std_residuals,
      w = weights,
      spec = spec,
      k = k,
      verbose = FALSE
    )
  }, warning = function(w) {
    warnings_list <<- c(warnings_list, list(w))
    invokeRestart("muffleWarning")
  })
  
  ## Diagnostics
  if (verbose) {
    cat("\n=== DCC M-STEP DIAGNOSTIC ===\n")
    cat("Starting DCC params:\n")
    print(unlist(dcc_start_pars))
    cat("\nOptimized DCC params:\n")
    print(opt_result$par)
    cat("\nOptimization convergence:", opt_result$convergence, "\n")
    cat("Final objective value:", opt_result$value, "\n")
    cat("Optimizer convergence code:", opt_result$convergence, "\n")
    cat("  0 = success, 1 = maxit reached, 51/52 = warning\n")
    cat("Starting NLL:", test_nll, "\n")
    cat("Final NLL:", opt_result$value, "\n")
    cat("Improvement:", test_nll - opt_result$value, "\n")
    cat("Function evaluations:", opt_result$counts[1], "\n")
  }
  
  alpha_params <- opt_result$par[grepl("alpha", names(opt_result$par))]
  at_boundary <- any(alpha_params < 0.02)
  
  if (verbose) cat("At boundary:", at_boundary, "\n")
  
  ## DIAGNOSTIC: Record boundary event
  # if (at_boundary && !is.null(diagnostics) && !is.null(iteration)) {
  #   for (pname in names(alpha_params)) {
  #     if (alpha_params[[pname]] < 0.02) {
  #       diagnostics <- add_boundary_event(
  #         diagnostics,
  #         iteration = iteration,
  #         state = state,
  #         parameter_name = pname,
  #         value = alpha_params[[pname]],
  #         boundary_type = "lower",
  #         action_taken = "constant_correlation_fallback"
  #       )
  #     }
  #   }
  #   
  #   if (verbose) {
  #     cat("\n*** WARNING: DCC alpha parameter at/near boundary ***\n")
  #     cat("This state appears to have CONSTANT (not dynamic) correlation.\n")
  #   }
  # }
  
  ## Compute the weighted log-likelihood for model comparison
  ## This is the objective we just optimized (but positive, since we minimized NLL)
  weighted_ll <- -opt_result$value
  
  ## Extract results
  estimated_pars <- as.list(opt_result$par)
  names(estimated_pars) <- names(all_stage2_pars)
  
  dcc_param_names <- names(dcc_start_pars)
  dist_param_names <- names(dist_start_pars)
  
  dcc_pars_final <- estimated_pars[dcc_param_names]
  dist_pars_final <- estimated_pars[dist_param_names]
  
  return(list(
    dcc_pars = dcc_pars_final,
    dist_pars = dist_pars_final,
    weighted_ll = weighted_ll,  # ADD THIS
    warnings = warnings_list,
    diagnostics = diagnostics
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
#' @description
#' *** THIS FUNCTION IS NOW DEPRECATED. Using perform_m_step_r() instead ***
#' This function is called once per EM iteration from C++. It uses the 'future'
#' framework to estimate the parameters for all M states in parallel. Handles
#' both univariate and multivariate models with proper parameter structuring.
#' @param y The time series data.
#' @param weights The (T x M) matrix of smoothed probabilities from the E-step.
#' @param spec The full list of model specifications.
#' @param model_type "univariate" or "multivariate".
#' @return A list of length M containing the updated model fits for each state.
perform_m_step_parallel_r <- function(
    y, 
    weights, 
    spec, 
    model_type,
    diagnostics = NULL, 
    iteration = NULL,
    verbose = FALSE
  ) {
  
  required_packages <- c("data.table", "xts", "tsgarch", "tsmarch", 
                         "tsdistributions", "mvtnorm")
  
  updated_fits <- future.apply::future_lapply(1:length(spec), function(j) {
    
    library(data.table)
    library(xts)
    library(tsgarch)
    library(tsmarch)
    library(tsdistributions)
    library(mvtnorm)
    
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
      model_type = model_type,
      diagnostics = diagnostics,
      iteration = iteration,
      state = j,  # Current state number
      verbose = verbose
    )
    
    ## === Structure the output ===
    if (model_type == "univariate") {
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
      ## === MULTIVARIATE ===
      variance_coeffs <- new_variance_fit$coefficients
      
      result <- list(var_pars = new_mean_fit$coefficients)
      
      ## Add univariate GARCH parameters
      if (!is.null(variance_coeffs$garch_pars)) {
        result$garch_pars <- variance_coeffs$garch_pars
      }
      
      ## Handle DCC parameters - may be empty for constant correlation
      if (!is.null(variance_coeffs$dcc_pars) && length(variance_coeffs$dcc_pars) > 0) {
        ## Dynamic correlation: flatten DCC parameters to top level
        for (par_name in names(variance_coeffs$dcc_pars)) {
          result[[par_name]] <- variance_coeffs$dcc_pars[[par_name]]
        }
        result$correlation_type <- "dynamic"
      } else {
        ## Constant correlation: no DCC parameters
        result$correlation_type <- "constant"
      }
      
      ## Handle copula/GOGARCH parameters (future)
      if (!is.null(variance_coeffs$copula_pars)) {
        for (par_name in names(variance_coeffs$copula_pars)) {
          result[[par_name]] <- variance_coeffs$copula_pars[[par_name]]
        }
      }
      
      if (!is.null(variance_coeffs$rotation_pars)) {
        for (par_name in names(variance_coeffs$rotation_pars)) {
          result[[par_name]] <- variance_coeffs$rotation_pars[[par_name]]
        }
      }
      
      ## Add distribution parameters
      if (!is.null(variance_coeffs$dist_pars)) {
        result$dist_pars <- variance_coeffs$dist_pars
      }
      
      return(result)
    }
    
  }, future.seed = TRUE, future.packages = required_packages)
  
  return(updated_fits)
}


#' @title Perform the M-Step (Sequential Execution)
#' @description This function is called once per EM iteration from C++. It
#' estimates the parameters for all M states sequentially, properly handling
#' diagnostics collection.
#' @param y The time series data.
#' @param weights The (T x M) matrix of smoothed probabilities from the E-step.
#' @param spec The full list of model specifications.
#' @param model_type "univariate" or "multivariate".
#' @param diagnostics Diagnostic collector object (ms_diagnostics class).
#' @param iteration Current EM iteration number.
#' @param verbose Logical indicating whether to print progress.
#' @return A list of length M containing the updated model fits for each state.
perform_m_step_r <- function(
    y, 
    weights, 
    spec, 
    model_type,
    diagnostics = NULL, 
    iteration = NULL,
    verbose = FALSE,
    dcc_threshold = 0.02,
    dcc_criterion = "bic"
) {
  
  M <- length(spec)
  updated_fits <- vector("list", M)
  
  ## Process each state sequentially
  for (j in 1:M) {
    
    state_weights <- weights[, j]
    state_spec <- spec[[j]]
    
    # ## === M-Step Stage 1: Update Mean Parameters ===
    # new_mean_fit <- estimate_arma_weighted_r(
    #   y = y,
    #   weights = state_weights,
    #   spec = state_spec,
    #   model_type = model_type
    # )
    # 
    # ## === M-Step Stage 2: Update Variance Parameters ===
    # new_variance_fit <- estimate_garch_weighted_r(
    #   residuals = new_mean_fit$residuals,
    #   weights = state_weights,
    #   spec = state_spec,
    #   model_type = model_type,
    #   diagnostics = diagnostics,  # Pass diagnostics directly
    #   iteration = iteration,
    #   state = j,
    #   verbose = verbose
    # )
    
    ## Check if this state is forced to constant correlation
    force_constant <- !is.null(state_spec$force_constant_correlation) && 
      state_spec$force_constant_correlation
    
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
      model_type = model_type,
      diagnostics = diagnostics,
      iteration = iteration,
      state = j,
      verbose = verbose,
      dcc_threshold = dcc_threshold,     # Pass through
      dcc_criterion = dcc_criterion,     # Pass through
      force_constant = force_constant    # NEW
    )
    
    ## Update diagnostics from this state
    if (!is.null(new_variance_fit$diagnostics)) {
      diagnostics <- new_variance_fit$diagnostics
    }
    
    ## === Structure the output ===
    if (model_type == "univariate") {
      all_params <- new_variance_fit$coefficients
      dist_param_names <- names(state_spec$start_pars$dist_pars)
      
      if (length(dist_param_names) > 0) {
        estimated_dist_pars <- all_params[dist_param_names]
        estimated_garch_pars <- all_params[!names(all_params) %in% dist_param_names]
      } else {
        estimated_dist_pars <- list()
        estimated_garch_pars <- all_params
      }
      
      updated_fits[[j]] <- list(
        arma_pars = new_mean_fit$coefficients,
        garch_pars = estimated_garch_pars,
        dist_pars = estimated_dist_pars
      )
      
    } else {
      ## === MULTIVARIATE ===
      variance_coeffs <- new_variance_fit$coefficients
      
      result <- list(var_pars = new_mean_fit$coefficients)
      
      ## Add univariate GARCH parameters
      if (!is.null(variance_coeffs$garch_pars)) {
        result$garch_pars <- variance_coeffs$garch_pars
      }
      
      ## Handle DCC parameters - may be empty for constant correlation
      if (!is.null(variance_coeffs$dcc_pars) && length(variance_coeffs$dcc_pars) > 0) {
        ## Dynamic correlation: flatten DCC parameters to top level
        for (par_name in names(variance_coeffs$dcc_pars)) {
          result[[par_name]] <- variance_coeffs$dcc_pars[[par_name]]
        }
        result$correlation_type <- "dynamic"
      } else {
        ## Constant correlation: no DCC parameters
        result$correlation_type <- "constant"
      }
      
      ## Handle copula/GOGARCH parameters (future)
      if (!is.null(variance_coeffs$copula_pars)) {
        for (par_name in names(variance_coeffs$copula_pars)) {
          result[[par_name]] <- variance_coeffs$copula_pars[[par_name]]
        }
      }
      
      if (!is.null(variance_coeffs$rotation_pars)) {
        for (par_name in names(variance_coeffs$rotation_pars)) {
          result[[par_name]] <- variance_coeffs$rotation_pars[[par_name]]
        }
      }
      
      ## Add distribution parameters
      if (!is.null(variance_coeffs$dist_pars)) {
        result$dist_pars <- variance_coeffs$dist_pars
      }
      
      ## Store degeneracy reason if present
      if (!is.null(variance_coeffs$degeneracy_reason)) {
        result$degeneracy_reason <- variance_coeffs$degeneracy_reason
      }
      
      updated_fits[[j]] <- result
    }
  }
  
  ## Return both fits and updated diagnostics
  return(list(
    fits = updated_fits,
    diagnostics = diagnostics
  ))
}


#' @title Compute Constant Correlation Weighted Likelihood
#' @keywords internal
compute_constant_correlation_likelihood <- function(
    residuals,
    weights,
    garch_pars,
    dist_pars,
    spec
) {
  
  k <- ncol(residuals)
  T_obs <- nrow(residuals)
  
  ## Get standardized residuals using GARCH parameters
  std_residuals <- matrix(0, nrow = T_obs, ncol = k)
  
  for (i in 1:k) {
    series_residuals <- residuals[, i]
    series_spec <- spec$garch_spec_args$garch_model$univariate[[i]]
    
    uni_spec_obj <- tsgarch::garch_modelspec(
      y = xts::xts(series_residuals, order.by = Sys.Date() - (T_obs:1)),
      model = series_spec$model,
      garch_order = series_spec$garch_order,
      distribution = series_spec$distribution
    )
    
    for (par_name in names(garch_pars[[i]])) {
      if (par_name %in% uni_spec_obj$parmatrix$parameter) {
        row_idx <- which(uni_spec_obj$parmatrix$parameter == par_name)
        uni_spec_obj$parmatrix[row_idx, "value"] <- garch_pars[[i]][[par_name]]
      }
    }
    
    uni_fit <- tsmethods::tsfilter(uni_spec_obj)
    sigma_vec <- as.numeric(uni_fit$sigma)
    std_residuals[, i] <- series_residuals / sigma_vec
  }
  
  ## Compute weighted unconditional correlation
  Qbar <- tryCatch({
    stats::cov.wt(std_residuals, wt = weights, method = "ML")$cov
  }, error = function(e) {
    # Fallback to unweighted if weighted fails
    cov(std_residuals)
  })
  
  ## Ensure positive definite
  eig <- eigen(Qbar, symmetric = TRUE)
  if (any(eig$values < 1e-8)) {
    Qbar <- Qbar + diag(1e-6, k)
  }
  
  ## Compute weighted log-likelihood with constant correlation
  ll_vec <- numeric(T_obs)
  
  if (spec$distribution == "mvn") {
    for (t in 1:T_obs) {
      ll_vec[t] <- mvtnorm::dmvnorm(
        std_residuals[t,], 
        mean = rep(0, k),
        sigma = Qbar, 
        log = TRUE
      )
    }
  } else if (spec$distribution == "mvt") {
    shape <- dist_pars$shape
    for (t in 1:T_obs) {
      ll_vec[t] <- mvtnorm::dmvt(
        std_residuals[t,], 
        delta = rep(0, k),
        sigma = Qbar,
        df = shape,
        log = TRUE
      )
    }
  }
  
  ll_vec[!is.finite(ll_vec)] <- -1e10
  weighted_ll <- sum(weights * ll_vec)
  
  return(list(
    weighted_ll = weighted_ll,
    dist_pars = dist_pars,
    warnings = list()
  ))
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
