## === === === === === === === === === === === === === === === === === 
## Generalized R Helper Functions for MS-ARMA-GARCH Fitting
## === === === === === === === === === === === === === === === === === 
## These functions are designed to be called from the C++ EM orchestrator
## fit_ms_varma_garch_cpp(). They contain the full logic for handling general 
## ARMA(p,q) and VAR(p) models, and interface with the tsgarch/tsmarch packages.



## tsmarch Version Compatibility ===============================================

#' Check if tsmarch supports higher-order DCC
#' 
#' @description tsmarch v1.0.0 has a bug in `.copula_parameters` where 
#'   `paste0("beta_",1:order[1])` should be `paste0("beta_",1:order[2])`.
#'   This causes DCC(p,q) with p!=q to create wrong number of parameters.
#'   Fixed in v1.0.1.
#'
#' @return Logical TRUE if tsmarch >= 1.0.1, FALSE otherwise
#' @keywords internal
tsmarch_supports_higher_order_dcc <- function() {
  tsmarch_version <- utils::packageVersion("tsmarch")
  return(tsmarch_version >= "1.0.1")
}


#' Validate DCC order against tsmarch version
#'
#' @description Checks if requested DCC order is supported by installed tsmarch.
#'   Issues informative error if higher-order DCC requested with buggy tsmarch.
#'
#' @param dcc_order Integer vector c(p, q) for DCC(p,q) order
#' @param action Character: "error" to stop, "warn" to warn and fall back to (1,1)
#' @return Validated/adjusted dcc_order
#' @keywords internal
validate_dcc_order <- function(dcc_order, action = c("error", "warn")) {
  action <- match.arg(action)
  
  
  ## Handle NULL or missing dcc_order
  if (is.null(dcc_order)) {
    return(c(1, 1))
  }
  
  ## DCC(1,1) always works
  if (length(dcc_order) >= 2 && all(dcc_order[1:2] == c(1, 1))) {
    return(dcc_order)
  }
  
  ## DCC(0,0) means constant correlation - always works
  if (length(dcc_order) >= 2 && all(dcc_order[1:2] == c(0, 0))) {
    return(dcc_order)
  }
  
  ## Check if higher-order is supported
  if (!tsmarch_supports_higher_order_dcc()) {
    tsmarch_version <- utils::packageVersion("tsmarch")
    
    msg <- sprintf(
      paste0(
        "Higher-order DCC(%d,%d) requested but tsmarch v%s has a bug with asymmetric orders.\n",
        "Options:\n",
        "  1. Install tsmarch >= 1.0.1: remotes::install_github('tsmodels/tsmarch')\n",
        "  2. Use DCC(1,1) instead"
      ),
      dcc_order[1], dcc_order[2], as.character(tsmarch_version)
    )
    
    if (action == "error") {
      stop(msg, call. = FALSE)
    } else {
      warning("\n", msg, "\nFalling back to DCC(1,1).\n", call. = FALSE)
      return(c(1, 1))
    }
  }
  
  return(dcc_order)
}


#' Validate DCC orders in model specification
#'
#' @description Validates all DCC orders in a multi-state spec list.
#'   Called during model fitting to ensure compatibility with installed tsmarch.
#'
#' @param spec List of state specifications
#' @param action Character: "error" to stop, "warn" to warn and fall back
#' @return Validated/adjusted spec (modified in place if needed)
#' @keywords internal
validate_spec_dcc_orders <- function(spec, action = c("error", "warn")) {
  action <- match.arg(action)
  
  for (j in seq_along(spec)) {
    state_spec <- spec[[j]]
    
    ## Check if this is a DCC model
    if (!is.null(state_spec$garch_spec_fun) && 
        state_spec$garch_spec_fun == "dcc_modelspec") {
      
      dcc_order <- state_spec$garch_spec_args$dcc_order
      
      if (!is.null(dcc_order)) {
        validated_order <- validate_dcc_order(dcc_order, action = action)
        
        ## Update spec if order was changed
        if (!identical(dcc_order, validated_order)) {
          spec[[j]]$garch_spec_args$dcc_order <- validated_order
        }
      }
    }
  }
  
  return(spec)
}


## DCC Boundary Warning Helper =================================================

#' @title Warn About DCC Parameters Near Boundary
#' @description Issues warnings when DCC parameters are near the lower boundary,
#'   which may indicate poorly identified correlation dynamics. Only called when
#'   returning dynamic correlation (not when falling back to constant).
#' @param alpha_params Named list of DCC alpha parameters
#' @param beta_params Named list of DCC beta parameters  
#' @param threshold Numeric threshold for boundary warning (default 1e-4)
#' @param state State index (for warning message)
#' @param iteration EM iteration (for diagnostics)
#' @param diagnostics Diagnostics collector object
#' @return Updated diagnostics object
#' @keywords internal
warn_dcc_boundary_params <- function(alpha_params, beta_params, threshold,
                                     state, iteration, diagnostics) {
  
  ## Check alpha parameters
  if (length(alpha_params) > 0) {
    for (pname in names(alpha_params)) {
      pval <- alpha_params[[pname]]
      if (!is.null(pval) && pval < threshold) {
        warning(sprintf(
          "\nState %s: DCC %s near boundary (%.2e < %.2e). Correlation dynamics may be poorly identified.\n",
          if (!is.null(state)) state else "?", pname, pval, threshold
        ))
        
        ## Log boundary event
        if (!is.null(diagnostics) && !is.null(iteration) && !is.null(state)) {
          diagnostics <- add_boundary_event(
            diagnostics,
            iteration = iteration,
            state = state,
            parameter_name = pname,
            value = pval,
            boundary_type = "lower",
            action_taken = "warning_issued_dynamic_returned"
          )
        }
      }
    }
  }
  
  ## Check beta parameters (important for higher-order DCC)
  if (length(beta_params) > 0) {
    for (pname in names(beta_params)) {
      pval <- beta_params[[pname]]
      if (!is.null(pval) && pval < threshold) {
        warning(sprintf(
          "\nState %s: DCC %s near boundary (%.2e < %.2e). Higher-order DCC parameters may be poorly identified.\n",
          if (!is.null(state)) state else "?", pname, pval, threshold
        ))
        
        ## Log boundary event
        if (!is.null(diagnostics) && !is.null(iteration) && !is.null(state)) {
          diagnostics <- add_boundary_event(
            diagnostics,
            iteration = iteration,
            state = state,
            parameter_name = pname,
            value = pval,
            boundary_type = "lower",
            action_taken = "warning_issued_dynamic_returned"
          )
        }
      }
    }
  }
  
  return(diagnostics)
}


## DCC Helper functions ========================================================

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
#' @description Convenience function that handles the complex, model-specific
#'   logic for creating both univariate and multivariate GARCH specification
#'   objects for the purpose of log-likelihood calculation during the
#'   EM-algorithm.
#' @param residuals Numeric
#' @param spec A spec list
#' @param model_type Character string
#' @param current_pars A list of parameters
#' @param diagnostics diagnostics
#' @param iteration iteration
#' @param state state
#' @param verbose verbose
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
            cat("Class of estimated_obj$sigma:", class(estimated_obj$sigma), "\n")
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
#' @param y y
#' @param current_pars current_pars 
#' @param spec spec
#' @param model_type model_type
#' @param diagnostics diagnostics
#' @param iteration iteration
#' @param state state
#' @param verbose verbose
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
  
    if (var_order > 0) {
      X_lagged <- matrix(1, nrow = T_obs - var_order, ncol = 1 + k * var_order)
      for (i in 1:var_order) {
        X_lagged[, (2 + (i - 1) * k):(1 + i*k)] <- y[(var_order - i + 1):(T_obs - i), ]
      }
      y_target <- y[(var_order + 1):T_obs, ]
      beta_mat <- matrix(current_pars$var_pars, nrow = 1 + k * var_order, ncol = k)
      model_residuals <- y_target - X_lagged %*% beta_mat
    } else {
      ## var_order = 0: No VAR dynamics, just intercept
      ## Residuals = y - intercept (or just y if no intercept)
      if (length(current_pars$var_pars) == k) {
        ## Intercept only
        intercepts <- matrix(current_pars$var_pars, nrow = 1, ncol = k)
        model_residuals <- sweep(y, 2, intercepts, "-")
      } else if (length(current_pars$var_pars) == 0 || is.null(current_pars$var_pars)) {
        ## No conditional mean parameters - use raw data
        model_residuals <- y
      } else {
        ## Fallback
        model_residuals <- y
      }
    }
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
            
            # cat("=== DEBUG: Qbar in garch_spec_obj ===\n")
            # print(garch_spec_obj$Qbar)
            # cat("=== DEBUG: class of garch_spec_obj ===\n")
            # print(class(garch_spec_obj))
            # cat("=== DEBUG: names in garch_spec_obj ===\n")
            # print(names(garch_spec_obj))
            
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
            
            # dcc_order <- garch_spec_obj$dynamics$order
            # maxpq <- max(dcc_order)
            # if (maxpq > 0) {
            #   total_nll_vec <- total_nll_vec[-(1:maxpq), , drop = TRUE]
            # } else {
            #   total_nll_vec <- as.vector(total_nll_vec)
            # }
            # ll_vector <- -total_nll_vec
            
            # dcc_order <- garch_spec_obj$dynamics$order
            # maxpq <- max(dcc_order)
            # ## Remove initialization placeholder (1 row) plus DCC burn-in (maxpq rows)
            # n_remove <- 1 + maxpq
            # total_nll_vec <- total_nll_vec[-(1:n_remove), , drop = TRUE]
            
            
            ## *** FIXED ***
            ## Remove only the initialization placeholder (first row, always 0)
            ## Note: tsmarch's type="nll" includes all subsequent observations, so we do too
            total_nll_vec <- total_nll_vec[-1, , drop = TRUE]
            
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
#' @param y y
#' @param weights weights
#' @param spec spec
#' @param model_type model_type
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
        method = "BFGS"#,
        #control = list(ndeps = rep(1e-12, length(start_pars))),
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
    k <- ncol(y)
    T_obs <- nrow(y)

    if (var_order > 0) {
      padding <- var_order
      X_lagged <- matrix(1, nrow = T_obs - padding, ncol = 1 + k * var_order)
      for (i in 1:var_order) {
        X_lagged[, (2 + (i-1)*k):(1 + i*k)] <- y[(padding-i+1):(T_obs-i), ]
      }
      y_target <- y[(padding+1):T_obs, ]
      w_target <- weights[(padding+1):T_obs]
      X_w <- X_lagged * sqrt(w_target)
      y_w <- y_target * sqrt(w_target)
      beta_mat <- try(solve(crossprod(X_w), crossprod(X_w, y_w)), silent = TRUE)
      if (inherits(beta_mat, "try-error")) {
        beta_mat <- matrix(0, nrow = 1 + k*var_order, ncol = k)
      }
      estimated_coeffs <- as.vector(beta_mat)
      final_residuals <- y_target - X_lagged %*% beta_mat
    } else {
      ## var_order = 0: No VAR dynamics
      ## Just compute weighted mean as intercept (or use zero intercept)
      y_target <- y
      w_target <- weights

      ## Weighted mean for each series as intercept
      w_sum <- sum(w_target)
      intercepts <- colSums(y_target * w_target) / w_sum

      estimated_coeffs <- intercepts  ## k intercepts
      final_residuals <- sweep(y_target, 2, intercepts, "-")
    }
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
#' @param diagnostics diagnostics
#' @param iteration iteration
#' @param state state
#' @param verbose verbose
#' @param dcc_threshold dcc_threshold
#' @param dcc_criterion dcc_criterion
#' @param force_constant force_constant
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
    stats::optim(
      par = unlist(start_pars),
      fn = weighted_garch_loglik,
      lower = lower_bounds,
      upper = upper_bounds,
      method = "L-BFGS-B",
      control = list(ndeps = rep(1e-12, length(start_pars))),
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


#' @title Multivariate GARCH Parameter Estimation (Two-Stage)
#' @description Implements weighted MLE for multivariate GARCH models
#'   Stage 1: Estimate univariate GARCH parameters for each series
#'   Stage 2: Estimate dependence parameters (DCC/Copula) or GOGARCH rotation
#' @param residuals residuals
#' @param weights weights
#' @param spec spec
#' @param diagnostics diagnostics
#' @param iteration iteration
#' @param state state
#' @param verbose verbose
#' @param dcc_threshold dcc_threshold
#' @param dcc_criterion dcc_criterion
#' @param force_constant force_constant
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
      dcc_threshold = dcc_threshold,
      dcc_criterion = dcc_criterion,
      force_constant = force_constant
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


#' @title DCC Weighted Estimation
#' @param residuals residuals
#' @param weights weights
#' @param spec spec
#' @param diagnostics diagnostics
#' @param iteration iteration
#' @param state state
#' @param verbose verbose
#' @param dcc_threshold dcc_threshold
#' @param dcc_criterion dcc_criterion
#' @param force_constant force_constant
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
  
  ## === CHECK FOR OMEGA BOUNDARY CONDITIONS ===
  ## tsgarch applies an omega boundary threshold of 1e-12, so 1e-8 is more
  ## tight, serving as a kind of "early warning".
  omega_boundary_threshold <- 1e-8

for (i in 1:k) {
  omega_est <- garch_pars_list[[i]]$omega
  
  if (!is.null(omega_est) && omega_est < omega_boundary_threshold) {
    ## Issue warning
    warning(sprintf(
      "\nState %s, Series %d: GARCH omega near boundary (%.2e < %.2e). Volatility dynamics may be degenerate.\n",
      if (!is.null(state)) state else "?", i, omega_est, omega_boundary_threshold
    ))
    
    ## Log boundary event
    if (!is.null(diagnostics) && !is.null(iteration) && !is.null(state)) {
      diagnostics <- add_boundary_event(
        diagnostics,
        iteration = iteration,
        state = state,
        parameter_name = sprintf("omega_series_%d", i),
        value = omega_est,
        boundary_type = "lower",
        action_taken = "warning_issued"
      )
    }
  }
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
  
  ## === CHECK FOR DCC PARAMETER BOUNDARY CONDITIONS ===
  ## Warn when any DCC parameter is near the lower boundary (suggests optimization issues)
  dcc_boundary_threshold <- 1e-4  ## More lenient than omega, but catches 0.0000 cases
  
  ## Check if parameters are near boundary (for constant correlation decision)
  near_boundary <- !is.null(alpha_params) && 
    length(alpha_params) > 0 && 
    any(unlist(alpha_params) < dcc_threshold)
  
  if (!near_boundary) {
    ## Parameters clearly away from boundary - use dynamic
    ## Check for boundary warnings before returning
    diagnostics <- warn_dcc_boundary_params(
      alpha_params, beta_params, dcc_boundary_threshold,
      state, iteration, diagnostics
    )
    
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
  
  ## Safety check - if weighted_ll wasn't returned, we have a problem
  if (is.null(ll_dynamic) || !is.finite(ll_dynamic)) {
    warning("\nCould not retrieve dynamic likelihood for comparison. Falling back to threshold criterion.\n")
    
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
    ### Use total sample size, sum(w_target^2)/sum(w_target)
    #n_eff <- length(w_target)
    n_eff <- (sum(w_target))^2 /sum(w_target^2)
    
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
    ## Keep dynamic - but warn if parameters are near boundary
    diagnostics <- warn_dcc_boundary_params(
      alpha_params, beta_params, dcc_boundary_threshold,
      state, iteration, diagnostics
    )
    
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


#' @title Estimate DCC Parameters with Weighted Likelihood
#' @description 
#' Estimates DCC parameters using weighted maximum likelihood, with specialized

#' handling for DCC(1,1) models (analytical gradients) and higher-order models
#' (finite difference with fine step size).
#'
#' For DCC(1,1), the optimization is performed in reparameterized space:
#'   psi = logit(alpha + beta)  -- persistence in logit space
#'   phi = log(alpha / beta)    -- ratio in log space
#'
#' This ensures stationarity by construction and provides smooth gradients.
#'
#' @param residuals T x k matrix of residuals
#' @param weights T-vector of observation weights
#' @param garch_pars List of GARCH parameters per series
#' @param dcc_start_pars Named list of starting DCC parameters
#' @param dist_start_pars Named list of distribution parameters
#' @param spec Model specification list
#' @param diagnostics Diagnostics object (optional)
#' @param iteration Current EM iteration (optional)
#' @param state Current regime state (optional)
#' @param verbose Logical; print diagnostic information
#' @return List with dcc_pars, dist_pars, weighted_ll, warnings, diagnostics
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
  
  
  ## 1. Get standardized residuals using Stage 1 GARCH parameters ==============
  
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
  
  
  ## 2. Handle empty parameter case ============================================
  
  all_stage2_pars <- c(dcc_start_pars, dist_start_pars)
  
  if (length(all_stage2_pars) == 0) {
    return(list(
      dcc_pars = list(), 
      dist_pars = list(), 
      weighted_ll = 0,
      warnings = list(),
      diagnostics = diagnostics
    ))
  }
  

  ## 3. Compute Qbar (weighted unconditional covariance) =======================
  
  Qbar <- tryCatch({
    stats::cov.wt(std_residuals, wt = weights, method = "ML")$cov
  }, error = function(e) {
    if (verbose) cat("ERROR in cov.wt:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(Qbar)) {
    return(list(
      dcc_pars = dcc_start_pars,
      dist_pars = dist_start_pars,
      weighted_ll = -1e10,
      warnings = list("Weighted covariance calculation failed"),
      diagnostics = diagnostics
    ))
  }
  
  ## Ensure positive definite
  eig <- eigen(Qbar, symmetric = TRUE)
  if (any(eig$values < 1e-8)) {
    Qbar <- Qbar + diag(1e-6, k)
  }
  
  
  ## 4. Determine optimization strategy based on DCC order =====================
  
  use_analytical_gradient <- is_dcc11(dcc_start_pars)
  
  if (verbose) {
    if (use_analytical_gradient) {
      cat("DCC(1,1) detected: Using analytical gradient with reparameterization\n")
    } else {
      order <- get_dcc_order(dcc_start_pars)
      cat(sprintf("DCC(%d,%d) detected: Using finite differences with ndeps=1e-12\n",
                  order["p"], order["q"]))
    }
  }
  
  warnings_list <- list()
  
  
  ## 5A. DCC(1,1): Reparameterized optimization with analytical gradient =======

  
  if (use_analytical_gradient) {
    
    ## Extract starting values
    alpha_start <- as.numeric(dcc_start_pars$alpha_1)
    beta_start <- as.numeric(dcc_start_pars$beta_1)
    
    ## Transform to unconstrained space
    unconstrained_start <- dcc11_to_unconstrained(alpha_start, beta_start)
    psi_start <- unconstrained_start["psi"]
    phi_start <- unconstrained_start["phi"]
    
    ## Build parameter vector for optimization
    if (spec$distribution == "mvt" && !is.null(dist_start_pars$shape)) {
      opt_params <- c(psi = psi_start, phi = phi_start, shape = dist_start_pars$shape)
      distribution <- "mvt"
    } else {
      opt_params <- c(psi = psi_start, phi = phi_start)
      distribution <- "mvn"
    }
    
    if (verbose) {
      cat("\nStarting values (original): alpha =", alpha_start, ", beta =", beta_start, "\n")
      cat("Starting values (reparam): psi =", psi_start, ", phi =", phi_start, "\n")
      if (distribution == "mvt") cat("Starting shape:", dist_start_pars$shape, "\n")
    }
    
    ## Define objective and gradient functions that capture diagnostics
    objective_fn <- function(params) {
      dcc11_nll(
        params = params,
        std_resid = std_residuals,
        weights = weights,
        Qbar = Qbar,
        distribution = distribution,
        use_reparam = TRUE
      )
    }
    
    gradient_fn <- function(params) {
      dcc11_gradient(
        params = params,
        std_resid = std_residuals,
        weights = weights,
        Qbar = Qbar,
        distribution = distribution,
        use_reparam = TRUE
      )
    }
    
    ## Test at starting point
    if (verbose) {
      cat("\n=== TESTING OBJECTIVE AT START PARAMS ===\n")
      test_nll <- objective_fn(opt_params)
      cat("Start NLL:", test_nll, "\n")
      test_grad <- gradient_fn(opt_params)
      cat("Start gradient:", test_grad, "\n")
    }
    
    ## Set bounds for shape parameter only (psi and phi are unconstrained)
    if (distribution == "mvt") {
      ## psi, phi unconstrained; shape > 2
      lower_bounds <- c(-Inf, -Inf, 2.01)
      upper_bounds <- c(Inf, Inf, Inf)
    } else {
      lower_bounds <- c(-Inf, -Inf)
      upper_bounds <- c(Inf, Inf)
    }
    
    ## Run optimization with analytical gradient
    opt_result <- withCallingHandlers({
      stats::optim(
        par = opt_params,
        fn = objective_fn,
        gr = gradient_fn,
        method = "L-BFGS-B",
        lower = lower_bounds,
        upper = upper_bounds,
        control = list(maxit = 500)
      )
    }, warning = function(w) {
      warnings_list <<- c(warnings_list, list(w))
      invokeRestart("muffleWarning")
    })
    
    ## Transform results back to original space
    ## Note: optim() may not preserve parameter names, so use positional indexing
    final_psi <- as.numeric(opt_result$par[1])
    final_phi <- as.numeric(opt_result$par[2])
    final_orig <- dcc11_from_unconstrained(final_psi, final_phi)
    
    dcc_pars_final <- list(
      alpha_1 = as.numeric(final_orig["alpha"]),
      beta_1 = as.numeric(final_orig["beta"])
    )
    
    if (distribution == "mvt") {
      dist_pars_final <- list(shape = as.numeric(opt_result$par[3]))
    } else {
      dist_pars_final <- dist_start_pars
    }
    
    if (verbose) {
      cat("\n=== DCC M-STEP DIAGNOSTIC (Analytical Gradient) ===\n")
      cat("Optimized (reparam): psi =", final_psi, ", phi =", final_phi, "\n")
      cat("Optimized (original): alpha =", dcc_pars_final$alpha_1, 
          ", beta =", dcc_pars_final$beta_1, "\n")
      cat("Persistence:", dcc_pars_final$alpha_1 + dcc_pars_final$beta_1, "\n")
      cat("Optimization convergence:", opt_result$convergence, "\n")
      cat("Final NLL:", opt_result$value, "\n")
      cat("Function evaluations:", opt_result$counts[1], "\n")
      cat("Gradient evaluations:", opt_result$counts[2], "\n")
    }
    
  } else {
    
    
    ## 5B. Higher-order DCC: Box-constrained with finite differences ===========
    
    ## Define objective function (existing implementation)
    weighted_dcc_loglik <- function(params, std_resid, w, spec, k, verbose = FALSE) {
      
      param_list <- as.list(params)
      names(param_list) <- names(all_stage2_pars)
      
      dcc_param_names <- names(dcc_start_pars)
      dist_param_names <- names(dist_start_pars)
      
      dcc_params_current <- param_list[dcc_param_names]
      dist_params_current <- param_list[dist_param_names]
      
      T_eff <- nrow(std_resid)
      
      ## Check stationarity using compute_dcc_persistence
      pers <- compute_dcc_persistence(dcc_params_current)
      
      if (pers$persistence >= 1 || any(pers$alphas < 0) || any(pers$betas < 0)) {
        return(1e10)
      }
      
      ## Weighted covariance (use pre-computed Qbar)
      Qbar_local <- Qbar
      
      ## DCC recursion
      recursion <- dcc_recursion(std_resid, Qbar_local, pers$alphas, pers$betas)
      
      if (!recursion$success) {
        return(1e10)
      }
      
      R <- recursion$R
      
      ## Calculate weighted log-likelihood
      ll_vec <- numeric(T_eff)
      
      if (spec$distribution == "mvn") {
        for (t in 1:T_eff) {
          R_t <- R[,,t]
          R_t_pd <- make_pd_safe(R_t)
          if (is.null(R_t_pd)) {
            ll_vec[t] <- -1e10
          } else {
            ll_vec[t] <- mvtnorm::dmvnorm(
              std_resid[t,], 
              mean = rep(0, k),
              sigma = R_t_pd, 
              log = TRUE
            )
          }
        }
      } else if (spec$distribution == "mvt") {
        shape <- dist_params_current$shape
        
        if (is.null(shape) || shape <= 2) {
          return(1e10)
        }
        
        for (t in 1:T_eff) {
          R_t <- R[,,t]
          R_t_pd <- make_pd_safe(R_t)
          if (is.null(R_t_pd)) {
            ll_vec[t] <- -1e10
          } else {
            ll_vec[t] <- mvtnorm::dmvt(
              std_resid[t,], 
              delta = rep(0, k), 
              sigma = R_t_pd, 
              df = shape, 
              log = TRUE
            )
          }
        }
      }
      
      ll_vec[!is.finite(ll_vec)] <- -1e10
      nll <- -sum(w * ll_vec, na.rm = TRUE)
      
      return(nll)
    }
    
    ## Get bounds
    bound_epsilon <- 1e-8
    lower_bounds <- numeric(length(all_stage2_pars))
    upper_bounds <- numeric(length(all_stage2_pars))
    
    for (i in seq_along(all_stage2_pars)) {
      par_name <- names(all_stage2_pars)[i]
      
      if (grepl("^(alpha|beta|gamma)_[0-9]+$", par_name)) {
        lower_bounds[i] <- bound_epsilon
        upper_bounds[i] <- 1 - bound_epsilon
      } else if (par_name == "shape") {
        lower_bounds[i] <- 2.01
        upper_bounds[i] <- 1e10
      } else {
        lower_bounds[i] <- bound_epsilon
        upper_bounds[i] <- 1 - bound_epsilon
      }
    }
    
    ## Ensure starting values are within bounds
    start_pars <- unlist(all_stage2_pars)
    for (i in seq_along(start_pars)) {
      if (start_pars[i] <= lower_bounds[i]) {
        start_pars[i] <- lower_bounds[i] + bound_epsilon
      }
      if (start_pars[i] >= upper_bounds[i]) {
        start_pars[i] <- upper_bounds[i] - bound_epsilon
      }
    }
    
    if (verbose) {
      cat("\n=== TESTING OBJECTIVE AT START PARAMS ===\n")
      test_nll <- weighted_dcc_loglik(
        params = start_pars,
        std_resid = std_residuals,
        w = weights,
        spec = spec,
        k = k,
        verbose = TRUE
      )
      cat("Start NLL:", test_nll, "\n")
    }
    
    ## Optimize with fine finite differences
    opt_result <- withCallingHandlers({
      stats::optim(
        par = start_pars,
        fn = weighted_dcc_loglik,
        lower = lower_bounds,
        upper = upper_bounds,
        method = "L-BFGS-B",
        control = list(
          ndeps = rep(1e-12, length(start_pars)),
          maxit = 500
        ),
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
    
    ## Extract results
    estimated_pars <- as.list(opt_result$par)
    names(estimated_pars) <- names(all_stage2_pars)
    
    dcc_param_names <- names(dcc_start_pars)
    dist_param_names <- names(dist_start_pars)
    
    dcc_pars_final <- estimated_pars[dcc_param_names]
    dist_pars_final <- estimated_pars[dist_param_names]
    
    if (verbose) {
      cat("\n=== DCC M-STEP DIAGNOSTIC (Finite Differences) ===\n")
      cat("Starting DCC params:\n")
      print(unlist(dcc_start_pars))
      cat("\nOptimized DCC params:\n")
      print(opt_result$par)
      cat("\nOptimization convergence:", opt_result$convergence, "\n")
      cat("Final NLL:", opt_result$value, "\n")
      cat("Function evaluations:", opt_result$counts["function"], "\n")
    }
  }
  
  
  ## 6. Return results =========================================================
  
  weighted_ll <- -opt_result$value
  
  return(list(
    dcc_pars = dcc_pars_final,
    dist_pars = dist_pars_final,
    weighted_ll = weighted_ll,
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
#' @param diagnostics diagnostics
#' @param iteration iteration
#' @param verbose verbose
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
#' @param dcc_threshold dcc_threshold
#' @param dcc_criterion dcc_criterion
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
      dcc_threshold = dcc_threshold,     ## Pass through
      dcc_criterion = dcc_criterion,     ## Pass through
      force_constant = force_constant
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
#' @param residuals residuals
#' @param weights weights
#' @param garch_pars garch_pars
#' @param dist_pars dist_pars
#' @param spec spec
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



## DCC Higher-Order Support Fix ================================================
##
## 1. Computing total persistence from all alpha and beta parameters
## 2. Implementing proper DCC(p,q) recursion
## 3. Stationarity constraint enforcement
##
## The standard DCC(p,q) model is:
##   Q_t = Qbar * (1 - sum(alpha) - sum(beta)) + 
##         sum_{j=1}^{q} alpha_j * (z_{t-j} * z_{t-j}') +
##         sum_{j=1}^{p} beta_j * Q_{t-j}
##
## Stationarity requires: sum(alpha) + sum(beta) < 1

#' Extract DCC order (p, q) from parameter names
#' 
#' @param dcc_params Named list of DCC parameters (e.g., alpha_1, alpha_2, beta_1)
#' @return Named vector c(p = ..., q = ...) where q is ARCH order (number of 
#'   alpha parameters) and p is GARCH order (number of beta parameters)
#' @keywords internal
get_dcc_order <- function(dcc_params) {
  
  if (is.null(dcc_params) || length(dcc_params) == 0) {
    return(c(p = 0, q = 0))
  }
  
  param_names <- names(dcc_params)
  
  ## Extract alpha indices (q order - ARCH-like, controls lagged z'z terms)
  alpha_names <- param_names[grepl("^alpha_", param_names)]
  if (length(alpha_names) > 0) {
    alpha_indices <- as.integer(gsub("^alpha_", "", alpha_names))
    q_order <- max(alpha_indices, na.rm = TRUE)
  } else {
    q_order <- 0
  }
  
  ## Extract beta indices (p order - GARCH-like, controls lagged Q terms)
  beta_names <- param_names[grepl("^beta_", param_names)]
  if (length(beta_names) > 0) {
    beta_indices <- as.integer(gsub("^beta_", "", beta_names))
    p_order <- max(beta_indices, na.rm = TRUE)
  } else {
    p_order <- 0
  }
  
  return(c(p = p_order, q = q_order))
}


#' Compute total persistence from DCC parameters
#' 
#' For a stationary DCC(p,q) model, we require:
#'   P = sum(alpha_1, ..., alpha_q) + sum(beta_1, ..., beta_p) < 1
#' 
#' @param dcc_params Named list of DCC parameters
#' @return List with components:
#'   \item{persistence}{Total persistence P = sum(alphas) + sum(betas)}
#'   \item{alpha_sum}{Sum of all alpha parameters}
#'   \item{beta_sum}{Sum of all beta parameters}
#'   \item{alphas}{Numeric vector of alpha values in order}
#'   \item{betas}{Numeric vector of beta values in order}
#'   \item{order}{Named vector c(p, q) giving the DCC order}
#' @keywords internal
compute_dcc_persistence <- function(dcc_params) {
  
  if (is.null(dcc_params) || length(dcc_params) == 0) {
    return(list(
      persistence = 0,
      alpha_sum = 0,
      beta_sum = 0,
      alphas = numeric(0),
      betas = numeric(0),
      order = c(p = 0, q = 0)
    ))
  }
  
  param_names <- names(dcc_params)
  
  ## Extract all alpha values in order
  alpha_names <- param_names[grepl("^alpha_", param_names)]
  if (length(alpha_names) > 0) {
    ## Sort by index to ensure correct order
    alpha_indices <- as.integer(gsub("^alpha_", "", alpha_names))
    alpha_order <- order(alpha_indices)
    alphas <- unlist(dcc_params[alpha_names[alpha_order]])
    names(alphas) <- NULL
  } else {
    alphas <- numeric(0)
  }
  
  ## Extract all beta values in order
  beta_names <- param_names[grepl("^beta_", param_names)]
  if (length(beta_names) > 0) {
    beta_indices <- as.integer(gsub("^beta_", "", beta_names))
    beta_order <- order(beta_indices)
    betas <- unlist(dcc_params[beta_names[beta_order]])
    names(betas) <- NULL
  } else {
    betas <- numeric(0)
  }
  
  alpha_sum <- sum(alphas)
  beta_sum <- sum(betas)
  
  return(list(
    persistence = alpha_sum + beta_sum,
    alpha_sum = alpha_sum,
    beta_sum = beta_sum,
    alphas = alphas,
    betas = betas,
    order = c(p = length(betas), q = length(alphas))
  ))
}


#' Check DCC stationarity constraints
#' 
#' Verifies that DCC parameters satisfy:
#' 1. All alpha_j >= 0
#' 2. All beta_j >= 0
#' 3. sum(alphas) + sum(betas) < 1
#' 
#' @param dcc_params Named list of DCC parameters
#' @param verbose Logical; if TRUE, print diagnostic messages
#' @return List with components:
#'   \item{is_stationary}{Logical indicating if constraints are satisfied}
#'   \item{persistence}{Total persistence value}
#'   \item{reason}{Character description if not stationary, NULL otherwise
#'   \item{details}{Output from compute_dcc_persistence()}
#' @keywords internal
check_dcc_stationarity <- function(dcc_params, verbose = FALSE) {
  
  pers <- compute_dcc_persistence(dcc_params)
  
  ## Check individual alpha positivity
  if (length(pers$alphas) > 0 && any(pers$alphas < 0)) {
    neg_idx <- which(pers$alphas < 0)
    reason <- sprintf("negative alpha at index %s (values: %s)", 
                      paste(neg_idx, collapse = ", "),
                      paste(round(pers$alphas[neg_idx], 6), collapse = ", "))
    if (verbose) cat("*** Stationarity violation:", reason, "***\n")
    return(list(
      is_stationary = FALSE,
      persistence = pers$persistence,
      reason = reason,
      details = pers
    ))
  }
  
  ## Check individual beta positivity
  if (length(pers$betas) > 0 && any(pers$betas < 0)) {
    neg_idx <- which(pers$betas < 0)
    reason <- sprintf("negative beta at index %s (values: %s)", 
                      paste(neg_idx, collapse = ", "),
                      paste(round(pers$betas[neg_idx], 6), collapse = ", "))
    if (verbose) cat("*** Stationarity violation:", reason, "***\n")
    return(list(
      is_stationary = FALSE,
      persistence = pers$persistence,
      reason = reason,
      details = pers
    ))
  }
  
  ## Check stationarity constraint: persistence < 1
  if (pers$persistence >= 1) {
    reason <- sprintf(
      "non-stationary (persistence = %.6f >= 1; alpha_sum = %.4f, beta_sum = %.4f)",
      pers$persistence, pers$alpha_sum, pers$beta_sum
    )
    if (verbose) cat("*** Stationarity violation:", reason, "***\n")
    return(list(
      is_stationary = FALSE,
      persistence = pers$persistence,
      reason = reason,
      details = pers
    ))
  }
  
  return(list(
    is_stationary = TRUE,
    persistence = pers$persistence,
    reason = NULL,
    details = pers
  ))
}


#' Transform DCC(1,1) parameters to reparameterized space
#' 
#' Transforms (alpha, beta) with constraint alpha + beta < 1 to 
#' (persistence, ratio) where both are in (0, 1) with no joint constraint.
#' 
#' @param alpha DCC alpha parameter
#' @param beta DCC beta parameter
#' @return Named vector c(persistence, ratio)
#' @keywords internal
dcc_to_reparam <- function(alpha, beta) {
  persistence <- alpha + beta
  
  ## Handle edge case where persistence is 0
  
  if (persistence < .Machine$double.eps) {
    ratio <- 0.5  # Arbitrary, since both are ~0
  } else {
    ratio <- alpha / persistence
  }
  
  c(persistence = persistence, ratio = ratio)
}


#' Transform reparameterized space back to DCC(1,1) parameters
#' 
#' Transforms (persistence, ratio) back to (alpha, beta).
#' Stationarity (alpha + beta < 1) is guaranteed if persistence < 1.
#' 
#' @param persistence Total persistence (alpha + beta), must be in (0, 1)
#' @param ratio Proportion allocated to alpha, must be in (0, 1)
#' @return Named vector c(alpha, beta)
#' @keywords internal
dcc_from_reparam <- function(persistence, ratio) {
  alpha <- persistence * ratio
  beta <- persistence * (1 - ratio)
  c(alpha = alpha, beta = beta)
}


#' Check if DCC specification is DCC(1,1)
#' 
#' @param dcc_params Named list of DCC parameters
#' @return Logical; TRUE if this is a DCC(1,1) specification
#' @keywords internal
is_dcc_11 <- function(dcc_params) {
  if (is.null(dcc_params) || length(dcc_params) == 0) {
    return(FALSE)
  }
  
  order <- get_dcc_order(dcc_params)
  return(order["p"] == 1 && order["q"] == 1)
}


#' Perform DCC(p,q) recursion
#' 
#' Computes the Q and R matrices for arbitrary DCC(p,q) order:
#'   Q_t = Qbar * (1 - sum(alpha) - sum(beta)) + 
#'         sum_{j=1}^{q} alpha_j * (z_{t-j} * z_{t-j}') +
#'         sum_{j=1}^{p} beta_j * Q_{t-j}
#'   R_t = diag(Q_t)^{-1/2} * Q_t * diag(Q_t)^{-1/2}
#' 
#' @param std_resid Matrix of standardized residuals (T x k)
#' @param Qbar Unconditional covariance matrix of standardized residuals (k x k)
#' @param alphas Numeric vector of alpha parameters (length q)
#' @param betas Numeric vector of beta parameters (length p)
#' @param verbose Logical; if TRUE, print diagnostic messages
#' @return List with components:
#'   \item{success}{Logical indicating if recursion completed without errors}
#'   \item{Q}{Array of Q matrices (k x k x T)}
#'   \item{R}{Array of correlation matrices (k x k x T)}
#'   \item{maxpq}{Maximum of p and q, used for burn-in period}
#'   \item{error_type}{Character describing error if success is FALSE}
#'   \item{error_time}{Time index where error occurred}
#' @keywords internal
dcc_recursion <- function(std_resid, Qbar, alphas, betas, verbose = FALSE) {
  
  T_obs <- nrow(std_resid)
  k <- ncol(std_resid)
  
  q_order <- length(alphas)  # Number of lagged z'z terms
  p_order <- length(betas)   # Number of lagged Q terms
  maxpq <- max(p_order, q_order, 1)
  
  ## Compute total persistence
  alpha_sum <- sum(alphas)
  beta_sum <- sum(betas)
  persistence <- alpha_sum + beta_sum
  
  if (verbose) {
    cat(sprintf("DCC(%d,%d) recursion: persistence = %.4f\n", 
                p_order, q_order, persistence))
  }
  
  ## Initialize arrays
  Q <- array(0, dim = c(k, k, T_obs))
  R <- array(0, dim = c(k, k, T_obs))
  
  ## Initialize first maxpq observations with Qbar
  for (t in 1:min(maxpq, T_obs)) {
    Q[,,t] <- Qbar
    Qbar_diag <- diag(Qbar)
    if (any(Qbar_diag <= 0)) {
      return(list(
        success = FALSE,
        error_type = "non_positive_Qbar_diagonal",
        error_time = 0,
        Q = Q,
        R = R,
        maxpq = maxpq
      ))
    }
    Qbar_diag_inv_sqrt <- diag(1/sqrt(Qbar_diag), k)
    R[,,t] <- Qbar_diag_inv_sqrt %*% Qbar %*% Qbar_diag_inv_sqrt
  }
  
  ## Handle case where T_obs <= maxpq (not enough data for recursion)
  if (T_obs <= maxpq) {
    return(list(
      success = TRUE,
      Q = Q,
      R = R,
      maxpq = maxpq
    ))
  }
  
  ## Main recursion starting from t = maxpq + 1
  for (t in (maxpq + 1):T_obs) {
    
    ## Intercept term: Qbar * (1 - persistence)
    Q_t <- Qbar * (1 - persistence)
    
    ## Add alpha terms: sum_{j=1}^{q} alpha_j * (z_{t-j} * z_{t-j}')
    for (j in seq_along(alphas)) {
      if (t - j >= 1) {
        z_lag <- std_resid[t - j, , drop = FALSE]
        Q_t <- Q_t + alphas[j] * (t(z_lag) %*% z_lag)
      }
    }
    
    ## Add beta terms: sum_{j=1}^{p} beta_j * Q_{t-j}
    for (j in seq_along(betas)) {
      if (t - j >= 1) {
        Q_t <- Q_t + betas[j] * Q[,,t - j]
      }
    }
    
    ## Check for numerical issues in Q
    if (any(!is.finite(Q_t))) {
      return(list(
        success = FALSE,
        error_type = "non_finite_Q",
        error_time = t,
        Q = Q,
        R = R,
        maxpq = maxpq
      ))
    }
    
    Q_diag <- diag(Q_t)
    if (any(Q_diag <= 0)) {
      return(list(
        success = FALSE,
        error_type = "non_positive_Q_diagonal",
        error_time = t,
        Q_diag_min = min(Q_diag),
        Q = Q,
        R = R,
        maxpq = maxpq
      ))
    }
    
    Q[,,t] <- Q_t
    
    ## Compute correlation matrix: R_t = diag(Q_t)^{-1/2} * Q_t * diag(Q_t)^{-1/2}
    Q_diag_inv_sqrt <- diag(1/sqrt(Q_diag), k)
    R[,,t] <- Q_diag_inv_sqrt %*% Q_t %*% Q_diag_inv_sqrt
  }
  
  return(list(
    success = TRUE,
    Q = Q,
    R = R,
    maxpq = maxpq
  ))
}


#' Reparameterized DCC(1,1) negative log-likelihood
#' 
#' This function computes the weighted negative log-likelihood for DCC(1,1)
#' using the (persistence, ratio) parameterization. This eliminates the
#' need for penalty-based stationarity enforcement since persistence < 1
#' is guaranteed by the box constraints.
#' 
#' @param reparam_pars Numeric vector c(persistence, ratio) or named list
#' @param std_resid Matrix of standardized residuals (T x k)
#' @param weights Observation weights (length T)
#' @param Qbar Unconditional covariance matrix (k x k)
#' @param distribution Character; either "mvn" or "mvt"
#' @param dist_pars Distribution parameters (e.g., shape for mvt)
#' @param verbose Logical; if TRUE, print diagnostic messages
#' @return Negative log-likelihood value (scalar)
#' @keywords internal
dcc11_nll_reparam <- function(reparam_pars, std_resid, weights, Qbar,
                              distribution = "mvn", dist_pars = NULL,
                              verbose = FALSE) {
  
  ## Extract reparameterized values
  if (is.list(reparam_pars)) {
    persistence <- reparam_pars$persistence
    ratio <- reparam_pars$ratio
  } else {
    persistence <- reparam_pars[1]
    ratio <- reparam_pars[2]
  }
  
  ## Check box constraints (should be enforced by optimizer, but be safe)
  eps <- 1e-10
  if (persistence <= eps || persistence >= 1 - eps ||
      ratio <= eps || ratio >= 1 - eps) {
    return(1e10)
  }
  
  ## Transform back to (alpha, beta)
  ## Stationarity is GUARANTEED since persistence < 1
  ab <- dcc_from_reparam(persistence, ratio)
  alpha <- ab["alpha"]
  beta <- ab["beta"]
  
  if (verbose) {
    cat(sprintf("Reparam: persistence=%.4f, ratio=%.4f -> alpha=%.4f, beta=%.4f\n",
                persistence, ratio, alpha, beta))
  }
  
  ## Perform DCC recursion
  recursion <- dcc_recursion(
    std_resid = std_resid,
    Qbar = Qbar,
    alphas = alpha,
    betas = beta,
    verbose = verbose
  )
  
  if (!recursion$success) {
    if (verbose) {
      cat("*** Recursion failed:", recursion$error_type, 
          "at t =", recursion$error_time, "***\n")
    }
    return(1e10)
  }
  
  ## Compute weighted log-likelihood
  T_obs <- nrow(std_resid)
  k <- ncol(std_resid)
  R <- recursion$R
  
  ll_vec <- numeric(T_obs)
  
  if (distribution == "mvn") {
    for (t in 1:T_obs) {
      R_t <- R[,,t]
      
      ## Check for valid correlation matrix
      if (any(!is.finite(R_t))) {
        ll_vec[t] <- -1e10
        next
      }
      
      ## Compute MVN log-likelihood
      ## For standardized residuals with correlation R:
      ## ll = -0.5 * (k*log(2*pi) + log(det(R)) + z' * R^{-1} * z)
      ## But z already has unit variance, so:
      ## ll = -0.5 * (log(det(R)) + z' * R^{-1} * z)  [up to constant]
      
      det_R <- det(R_t)
      if (det_R <= 0) {
        ll_vec[t] <- -1e10
        next
      }
      
      R_inv <- tryCatch(solve(R_t), error = function(e) NULL)
      if (is.null(R_inv)) {
        ll_vec[t] <- -1e10
        next
      }
      
      z_t <- std_resid[t, ]
      mahal <- as.numeric(t(z_t) %*% R_inv %*% z_t)
      
      ll_vec[t] <- -0.5 * (k * log(2 * pi) + log(det_R) + mahal)
    }
  } else {
    ## MVT distribution - would need shape parameter
    ## For now, fall back to MVN
    warning("\nMVT distribution not yet implemented in reparameterized version.\n")
    return(1e10)
  }
  
  ## Compute weighted NLL
  ll_vec[!is.finite(ll_vec)] <- -1e10
  nll <- -sum(weights * ll_vec)
  
  return(nll)
}


#' Estimate DCC parameters with appropriate method
#' 
#' For DCC(1,1), uses reparameterized optimization to avoid boundary
#' instabilities. For higher-order DCC(p,q), uses penalty-based optimization.
#' 
#' @param std_resid Matrix of standardized residuals (T x k)
#' @param weights Observation weights (length T)
#' @param dcc_start_pars Named list of starting DCC parameters
#' @param Qbar Unconditional covariance matrix (k x k)
#' @param distribution Character; "mvn" or "mvt"
#' @param dist_pars Distribution parameters
#' @param verbose Logical; if TRUE, print diagnostic messages
#' @return List with estimated parameters and diagnostics
#' @keywords internal
estimate_dcc_params <- function(std_resid, weights, dcc_start_pars, Qbar,
                                distribution = "mvn", dist_pars = NULL,
                                verbose = FALSE) {
  
  ## Determine DCC order
  order <- get_dcc_order(dcc_start_pars)
  
  if (verbose) {
    cat(sprintf("DCC order: p=%d, q=%d\n", order["p"], order["q"]))
  }
  
  ## Choose estimation method based on order
  if (order["p"] == 1 && order["q"] == 1) {
    ## DCC(1,1): Use reparameterized optimization
    result <- estimate_dcc11_reparam(
      std_resid = std_resid,
      weights = weights,
      dcc_start_pars = dcc_start_pars,
      Qbar = Qbar,
      distribution = distribution,
      dist_pars = dist_pars,
      verbose = verbose
    )
    result$method <- "reparameterized"
    
  } else {
    ## Higher order: Use penalty-based optimization
    result <- estimate_dcc_pq_penalty(
      std_resid = std_resid,
      weights = weights,
      dcc_start_pars = dcc_start_pars,
      Qbar = Qbar,
      distribution = distribution,
      dist_pars = dist_pars,
      verbose = verbose
    )
    result$method <- "penalty"
  }
  
  return(result)
}


#' Estimate DCC(1,1) using reparameterized optimization
#' @keywords internal
estimate_dcc11_reparam <- function(std_resid, weights, dcc_start_pars, Qbar,
                                   distribution = "mvn", dist_pars = NULL,
                                   verbose = FALSE) {
  
  ## Extract starting values
  alpha_start <- dcc_start_pars$alpha_1
  beta_start <- dcc_start_pars$beta_1
  
  ## Transform to reparameterized space
  reparam_start <- dcc_to_reparam(alpha_start, beta_start)
  
  if (verbose) {
    cat(sprintf("Starting values: alpha=%.4f, beta=%.4f\n", alpha_start, beta_start))
    cat(sprintf("Reparameterized: persistence=%.4f, ratio=%.4f\n",
                reparam_start["persistence"], reparam_start["ratio"]))
  }
  
  ## Set up bounds for reparameterized parameters
  ## Both persistence and ratio should be in (eps, 1-eps)
  eps <- 1e-6
  lower <- c(eps, eps)
  upper <- c(1 - eps, 1 - eps)
  
  ## Optimize
  opt_result <- optim(
    par = as.numeric(reparam_start),
    fn = dcc11_nll_reparam,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(ndeps = rep(1e-12, length(reparam_start))),
    std_resid = std_resid,
    weights = weights,
    Qbar = Qbar,
    distribution = distribution,
    dist_pars = dist_pars,
    verbose = FALSE
  )
  
  ## Transform back to original parameters
  final_ab <- dcc_from_reparam(opt_result$par[1], opt_result$par[2])
  
  if (verbose) {
    cat(sprintf("Optimized: persistence=%.4f, ratio=%.4f\n",
                opt_result$par[1], opt_result$par[2]))
    cat(sprintf("Final: alpha=%.4f, beta=%.4f (sum=%.4f)\n",
                final_ab["alpha"], final_ab["beta"], sum(final_ab)))
    cat(sprintf("Convergence: %d, NLL: %.4f\n", opt_result$convergence, opt_result$value))
  }
  
  return(list(
    dcc_pars = list(alpha_1 = final_ab["alpha"], beta_1 = final_ab["beta"]),
    reparam_pars = list(persistence = opt_result$par[1], ratio = opt_result$par[2]),
    nll = opt_result$value,
    convergence = opt_result$convergence,
    n_penalty_triggers = 0  # Reparameterization eliminates penalties
  ))
}


#' Estimate DCC(p,q) using penalty-based optimization
#' @keywords internal
estimate_dcc_pq_penalty <- function(std_resid, weights, dcc_start_pars, Qbar,
                                    distribution = "mvn", dist_pars = NULL,
                                    verbose = FALSE) {
  
  ## This would contain the corrected penalty-based optimization
  ## using compute_dcc_persistence() and dcc_recursion()
  ## 
  ## For now, return a placeholder indicating this needs integration
  ## with the existing weighted_dcc_loglik infrastructure
  
  warning("\nHigher-order DCC estimation requires integration with existing code.\n")
  
  return(list(
    dcc_pars = dcc_start_pars,
    nll = NA,
    convergence = -1,
    n_penalty_triggers = NA
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
