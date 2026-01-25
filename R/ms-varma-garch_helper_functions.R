## === === === === === === === === === === === === === === === === === 
## Generalized R Helper Functions for MS-ARMA-GARCH Fitting
## === === === === === === === === === === === === === === === === === 
## These functions are designed to be called from the C++ EM orchestrator
## fit_ms_varma_garch_cpp(). They contain the full logic for handling general 
## ARMA(p,q) and VAR(p) models, and interface with the tsgarch/tsmarch packages.
## See also dcc.R and cgarch.R


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
      
      # ## STEP 3: Create DCC spec from the multi_estimate object
      # final_args <- c(
      #   list(object = multi_estimate_object), 
      #   spec_args[names(spec_args) != "garch_model"]
      # )
      # final_args$distribution <- spec$distribution
      # 
      # ## Must explicitly specify dynamics = "dcc" for dynamic correlation
      # ## Without this, tsmarch defaults to constant correlation!
      # if (is.null(final_args$dynamics)) {
      #   final_args$dynamics <- "dcc"
      # }
      # 
      # garch_spec_obj <- do.call(spec_fun, final_args)
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
            
            ## Align lengths if needed
            if (length(garch_nll_vec) > length(copula_nll_vec)) {
              n_diff <- length(garch_nll_vec) - length(copula_nll_vec)
              garch_nll_vec <- garch_nll_vec[-(1:n_diff)]
            }
            
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
      ## Use the existing ICA decomposition and GARCH parameters from current_pars
      ## instead of re-estimating from scratch
      
      if (!is.null(current_pars$ica_info) && !is.null(current_pars$garch_pars)) {
        ## Use existing parameters for likelihood computation
        ll_vector <- compute_gogarch_loglik_ms(
          residuals = model_residuals,
          garch_pars = current_pars$garch_pars,
          ica_info = current_pars$ica_info,
          distribution = spec$distribution %||% "norm",
          return_vector = TRUE
        )
      } else {
        ## Fallback: estimate from scratch (first iteration or missing params)
        garch_model_fit <- suppressWarnings(estimate(garch_spec_obj, trace = FALSE))
        ll_vector <- compute_loglik_fixed(
          object = garch_model_fit,
          params = list(),
          ll_vec = TRUE
        )
      }
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

  if (model_type %in% c("dcc_modelspec")) {
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
  } else if (model_type == "cgarch_modelspec") {
    return(estimate_garch_weighted_cgarch(
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
      } else if (!is.null(variance_coeffs$correlation_type)) {
        ## Preserve correlation_type if already set (e.g., by GOGARCH or constant)
        result$correlation_type <- variance_coeffs$correlation_type
      } else {
        ## Default to constant correlation if no DCC parameters
        result$correlation_type <- "constant"
      }
      
      ## Handle copula parameters
      if (!is.null(variance_coeffs$copula_pars)) {
        for (par_name in names(variance_coeffs$copula_pars)) {
          result[[par_name]] <- variance_coeffs$copula_pars[[par_name]]
        }
      }
      
      # Handle GOGARCH ICA information
      if (!is.null(variance_coeffs$ica_info)) {
        result$ica_info <- variance_coeffs$ica_info
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
