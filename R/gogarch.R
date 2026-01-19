## ============================================================================
## tsbs_gogarch.R
## GOGARCH-specific functions for MS-VARMA-GARCH framework
## ============================================================================

#' Weighted GARCH Estimation for GOGARCH Models
#'
#' @description
#' Implements weighted maximum likelihood estimation for Generalized Orthogonal
#' GARCH (GOGARCH) models in the context of Markov-Switching frameworks.
#'
#' GOGARCH differs fundamentally from DCC and Copula GARCH in how it models
#' dynamic correlations:
#' \itemize{
#'   \item \strong{DCC/CGARCH}: Estimate explicit correlation dynamics parameters
#'     (alpha, beta) that govern how correlations evolve over time.
#'   \item \strong{GOGARCH}: Assume observations arise from independent latent
#'     factors. Time-varying correlations emerge from the time-varying volatilities
#'     of these factors, transformed through a fixed mixing matrix.
#' }
#'
#' @details
#' The estimation proceeds in two stages:
#'
#' \strong{Stage 1: ICA Decomposition}
#'
#' Independent Component Analysis extracts statistically independent components
#' from the multivariate residuals:
#' \deqn{S = Y \cdot W}
#' where \eqn{Y} is the T x k residual matrix, \eqn{W} is the unmixing matrix,
#' and \eqn{S} contains the independent components.
#'
#' The RADICAL algorithm (Learned-Miller, 2003) is used for ICA decomposition,
#' which is robust to outliers. For improved convergence with multiple restarts
#' and quality diagnostics, see \code{\link{improved_ica_decomposition}}.
#'
#' \strong{Stage 2: Component GARCH Estimation}
#'
#' Univariate GARCH models are fitted to each independent component using
#' weighted MLE. The weights come from the state probabilities in the
#' Markov-Switching framework.
#'
#' \strong{Covariance Reconstruction}
#'
#' The time-varying covariance matrix is reconstructed as:
#' \deqn{H_t = A \cdot D_t \cdot A'}
#' where \eqn{A} is the mixing matrix from ICA and \eqn{D_t = diag(\sigma^2_{1,t},
#' \ldots, \sigma^2_{k,t})} contains the component GARCH variances.
#'
#' \strong{Log-Likelihood}
#'
#' The GOGARCH log-likelihood includes a Jacobian adjustment for the ICA
#' transformation:
#' \deqn{LL = \sum_i LL_{component,i} + \log|det(K)|}
#' where \eqn{K} is the pre-whitening matrix. See \code{\link{compute_gogarch_loglik_ms}}.
#'
#' @param residuals Numeric matrix of residuals with dimensions T x k, where T is
#'   the number of observations and k is the number of series.
#' @param weights Numeric vector of state probabilities/weights with length T.
#'   These weights come from the E-step of the EM algorithm in the MS-VARMA-GARCH
#'   framework.
#' @param spec List containing the model specification with the following elements:
#'   \describe{
#'     \item{\code{garch_spec_args}}{List with:
#'       \itemize{
#'         \item \code{model}: GARCH model type (default: \code{"garch"})
#'         \item \code{order}: GARCH order as c(p, q) (default: \code{c(1, 1)})
#'         \item \code{ica}: ICA algorithm. Currently only \code{"radical"} is
#'           supported (tsmarch v1.0.0)
#'         \item \code{components}: Number of ICA components to extract
#'       }
#'     }
#'     \item{\code{distribution}}{Component distribution: \code{"norm"},
#'       \code{"nig"}, or \code{"gh"}}
#'     \item{\code{start_pars}}{List with:
#'       \itemize{
#'         \item \code{garch_pars}: List of starting GARCH parameters for each
#'           component (e.g., \code{list(list(omega=0.1, alpha1=0.1, beta1=0.8), ...)})
#'         \item \code{dist_pars}: Distribution parameters (e.g.,
#'           \code{list(shape=1, skew=0)} for NIG)
#'       }
#'     }
#'   }
#' @param diagnostics Optional diagnostics collector object for logging
#'   estimation progress and issues.
#' @param verbose Logical; if \code{TRUE}, print progress information during
#'   ICA decomposition and GARCH estimation. Default is \code{FALSE}.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{\code{coefficients}}{List containing:
#'       \itemize{
#'         \item \code{garch_pars}: List of estimated GARCH parameters for each
#'           ICA component
#'         \item \code{ica_info}: List with ICA results:
#'           \itemize{
#'             \item \code{A}: Mixing matrix (k x n_components)
#'             \item \code{W}: Unmixing matrix (n_components x k)
#'             \item \code{K}: Pre-whitening matrix (for Jacobian)
#'             \item \code{S}: Independent components matrix (T x n_components)
#'             \item \code{method}: ICA algorithm used
#'             \item \code{n_components}: Number of components extracted
#'           }
#'         \item \code{dist_pars}: Distribution parameters
#'         \item \code{correlation_type}: Always \code{"gogarch"}
#'       }
#'     }
#'     \item{\code{warnings}}{List of warning messages generated during estimation}
#'     \item{\code{diagnostics}}{Updated diagnostics object}
#'   }
#'
#' @section Comparison with DCC and CGARCH:
#' \tabular{llll}{
#'   \strong{Aspect} \tab \strong{DCC} \tab \strong{CGARCH} \tab \strong{GOGARCH} \cr
#'   Correlation source \tab alpha, beta params \tab DCC on copula \tab ICA mixing matrix \cr
#'   Dynamics \tab DCC recursion \tab DCC/ADCC \tab Component volatilities \cr
#'   Marginal treatment \tab Assumed Normal \tab Flexible (PIT) \tab ICA components \cr
#'   Key strength \tab Interpretable \tab Tail dependence \tab Non-Gaussian dependence \cr
#' }
#'
#' @section When to Use GOGARCH:
#' GOGARCH is particularly suitable when:
#' \itemize{
#'   \item The dependence structure arises from independent underlying factors
#'   \item Heavy-tailed distributions (NIG, GH) are needed for the components
#'   \item Dimension reduction is desired (fewer components than series)
#'   \item Non-linear dependence structures are suspected
#' }
#' 
#' @section ICA Quality Assessment:
#' The quality of the ICA decomposition is critical for GOGARCH reliability.
#' Key diagnostics include:
#' \itemize{
#'   \item \strong{Independence score}: Measures how uncorrelated the extracted
#'     components are (target: > 0.8)
#'   \item \strong{Reconstruction error}: How well Y = S * A' holds (target: < 1\%)
#'   \item \strong{Negentropy}: Non-Gaussianity of components (higher = better separation
#' }
#'
#' If ICA convergence is problematic, consider using \code{\link{improved_ica_decomposition}}
#' which provides multiple restarts and automatic quality assessment. For post-estimation
#' diagnostics, see \code{\link{gogarch_diagnostics}}.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{tsbs}}: Main bootstrap function
#'     (use \code{garch_spec_fun = "gogarch_modelspec"})
#'   \item \code{\link{estimate_garch_weighted_univariate_gogarch}}: Component
#'     GARCH estimation
#'   \item \code{\link{compute_gogarch_loglik_ms}}: Log-likelihood computation
#'   \item \code{\link{improved_ica_decomposition}}: Enhanced ICA with multiple
#'     restarts and quality metrics
#'   \item \code{\link{gogarch_diagnostics}}: Comprehensive model diagnostics
#'     including ICA quality, component GARCH fit, and covariance reconstruction
#'   \item \code{\link{estimate_garch_weighted_dcc}}: Alternative DCC estimator
#'   \item \code{\link{estimate_garch_weighted_cgarch}}: Alternative CGARCH estimator
#' }
#'
#' @references
#' van der Weide, R. (2002). GO-GARCH: A multivariate generalized orthogonal
#' GARCH model. \emph{Journal of Applied Econometrics}, 17(5), 549-564.
#'
#' Learned-Miller, E. G. (2003). ICA using spacings estimates of entropy.
#' \emph{Journal of Machine Learning Research}, 4, 1271-1295.
#'
#' @keywords internal
estimate_garch_weighted_gogarch <- function(
    residuals, 
    weights, 
    spec, 
    diagnostics = NULL, 
    verbose = FALSE
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
  
  ## === STEP 1: Perform ICA Decomposition ===
  ## ICA extracts independent components from the residuals
  ## The mixing matrix A and unmixing matrix W are key for GOGARCH
  
  ## Convert to xts for tsmarch compatibility
  residuals_xts <- xts::xts(residuals, order.by = seq.Date(
    from = as.Date("2000-01-01"), 
    by = "day", 
    length.out = T_obs
  ))
  colnames(residuals_xts) <- paste0("series_", 1:k)
  
  ## Get ICA parameters from spec
  ica_method <- spec$garch_spec_args$ica %||% "radical"
  n_components <- spec$garch_spec_args$components %||% k
  
  ## Perform ICA decomposition
  ## Use tsmarch's radical or fastica implementation
  # ic <- tryCatch({
  #   if (ica_method == "radical") {
  #     tsmarch::radical(residuals, components = n_components, demean = FALSE, trace = verbose)
  #   } else if (ica_method == "fastica") {
  #     tsmarch::fastica(residuals, components = n_components, demean = FALSE, trace = verbose)
  #   } else {
  #     stop("Unsupported ICA method: ", ica_method)
  #   }
  # }, error = function(e) {
  #   warning("ICA decomposition failed: ", e$message, ". Using identity transformation.")
  #   ## Fallback: identity transformation (no rotation)
  #   list(
  #     S = residuals,
  #     A = diag(k),
  #     W = diag(k),
  #     K = diag(k)
  #   )
  # })
  
  ic <- improved_ica_decomposition(
    residuals = residuals,
    method = ica_method,
    n_components = n_components,
    n_restarts = 3,  # Multiple restarts
    verbose = verbose
  )
  
  ## Extract independent components
  S <- ic$S  # T x n_components matrix of independent components
  
  if (verbose) {
    cat(sprintf("GOGARCH: ICA decomposition complete. %d components extracted.\n", ncol(S)))
  }
  
  ## === STEP 2: Estimate Univariate GARCH on Each ICA Component (Weighted) ===
  
  garch_pars_list <- list()
  warnings_list <- list()
  
  garch_model <- spec$garch_spec_args$model %||% "garch"
  garch_order <- spec$garch_spec_args$order %||% c(1, 1)
  distribution <- spec$distribution %||% "norm"
  
  for (i in 1:n_components) {
    component_data <- S[, i]
    
    ## Create spec for this component
    comp_spec <- list(
      garch_model = garch_model,
      garch_order = garch_order,
      distribution = distribution,
      start_pars = list(
        garch_pars = spec$start_pars$garch_pars[[i]] %||% 
          list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
        dist_pars = spec$start_pars$dist_pars
      )
    )
    
    ## Estimate using weighted univariate function
    uni_result <- estimate_garch_weighted_univariate_gogarch(
      residuals = component_data,
      weights = w_target,
      spec = comp_spec,
      verbose = verbose
    )
    
    garch_pars_list[[i]] <- uni_result$coefficients
    warnings_list <- c(warnings_list, uni_result$warnings)
  }
  
  ## === STEP 3: Check for Boundary Conditions ===
  omega_boundary_threshold <- 1e-8
  
  for (i in 1:n_components) {
    omega_est <- garch_pars_list[[i]]$omega
    
    if (!is.null(omega_est) && omega_est < omega_boundary_threshold) {
      warning(sprintf(
        "GOGARCH Component %d: omega near boundary (%.2e < %.2e). Volatility dynamics may be degenerate.",
        i, omega_est, omega_boundary_threshold
      ))
    }
  }
  
  ## === Return Results ===
  return(list(
    coefficients = list(
      garch_pars = garch_pars_list,
      ica_info = list(
        A = ic$A,        # Mixing matrix
        W = ic$W,        # Unmixing matrix  
        K = ic$K,        # Pre-whitening matrix
        S = ic$S,        # Independent components (for reference)
        method = ica_method,
        n_components = n_components
      ),
      dist_pars = spec$start_pars$dist_pars,
      correlation_type = "gogarch"  # Always GOGARCH structure
    ),
    warnings = warnings_list,
    diagnostics = diagnostics
  ))
}


#' Weighted Univariate GARCH Estimation for GOGARCH Components
#'
#' @description
#' Estimates univariate GARCH parameters for a single ICA component using
#' weighted maximum likelihood. This function is called by
#' \code{\link{estimate_garch_weighted_gogarch}} for each independent component
#' extracted by ICA.
#'
#' @details
#' The function implements weighted MLE for GARCH(p,q) models on ICA components.
#' The weighted log-likelihood is:
#' \deqn{LL_w = \sum_t w_t \cdot \log f(s_t | \sigma_t)}
#' where \eqn{w_t} are the state probabilities, \eqn{s_t} is the ICA component
#' value, and \eqn{\sigma_t} is the GARCH conditional volatility.
#'
#' \strong{GARCH Recursion}
#'
#' The variance recursion for GARCH(p,q) is:
#' \deqn{\sigma^2_t = \omega + \sum_{i=1}^{p} \alpha_i s^2_{t-i} +
#'   \sum_{j=1}^{q} \beta_j \sigma^2_{t-j}}
#'
#' \strong{Supported Distributions}
#'
#' \tabular{lll}{
#'   \strong{Distribution} \tab \strong{Parameters} \tab \strong{Description} \cr
#'   \code{"norm"} \tab none \tab Normal (Gaussian) \cr
#'   \code{"std"} \tab shape \tab Student-t \cr
#'   \code{"sstd"} \tab shape, skew \tab Skewed Student-t \cr
#'   \code{"nig"} \tab shape, skew \tab Normal Inverse Gaussian \cr
#'   \code{"gh"} \tab shape, skew, lambda \tab Generalized Hyperbolic \cr
#' }
#'
#' @param residuals Numeric vector of ICA component values with length T.
#' @param weights Numeric vector of weights (state probabilities) with length T.
#' @param spec List containing the component specification:
#'   \describe{
#'     \item{\code{garch_model}}{GARCH model type (default: \code{"garch"})}
#'     \item{\code{garch_order}}{GARCH order as c(p, q) (default: \code{c(1, 1)})}
#'     \item{\code{distribution}}{Component distribution}
#'     \item{\code{start_pars}}{List with \code{garch_pars} and \code{dist_pars}}
#'   }
#' @param verbose Logical; if \code{TRUE}, print progress information.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{coefficients}}{Named list of estimated parameters}
#'     \item{\code{warnings}}{List of any warnings from optimization}
#'   }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{estimate_garch_weighted_gogarch}}: Main GOGARCH estimator
#'   \item \code{\link{compute_gogarch_loglik_ms}}: Log-likelihood computation
#' }
#'
#' @keywords internal
estimate_garch_weighted_univariate_gogarch <- function(
    residuals,
    weights,
    spec,
    verbose = FALSE
) {
  
  ## Extract parameters
  start_pars <- c(spec$start_pars$garch_pars, spec$start_pars$dist_pars)
  
  if (length(start_pars) == 0) {
    return(list(coefficients = list(), warnings = list()))
  }
  
  ## Create univariate GARCH spec using tsgarch
  garch_spec <- tsgarch::garch_modelspec(
    y = xts::xts(residuals, order.by = seq.Date(
      from = as.Date("2000-01-01"),
      by = "day", 
      length.out = length(residuals)
    )),
    model = spec$garch_model %||% "garch",
    order = spec$garch_order %||% c(1, 1),
    constant = FALSE,  # ICA components are demeaned
    distribution = spec$distribution %||% "norm"
  )
  
  ## Extract bounds
  parmatrix <- garch_spec$parmatrix
  pars_to_estimate <- names(start_pars)
  
  ## Match parameters - handle naming differences
  param_mapping <- list(
    omega = "omega",
    alpha1 = "alpha1", alpha2 = "alpha2",
    beta1 = "beta1", beta2 = "beta2",
    skew = "skew", shape = "shape", lambda = "lambda"
  )
  
  lower_bounds <- numeric(length(pars_to_estimate))
  upper_bounds <- numeric(length(pars_to_estimate))
  
  for (j in seq_along(pars_to_estimate)) {
    pname <- pars_to_estimate[j]
    if (pname %in% parmatrix$parameter) {
      row_idx <- which(parmatrix$parameter == pname)
      lower_bounds[j] <- parmatrix$lower[row_idx]
      upper_bounds[j] <- parmatrix$upper[row_idx]
    } else {
      ## Default bounds
      lower_bounds[j] <- 1e-12
      upper_bounds[j] <- 1 - 1e-6
    }
  }
  
  ## Weighted log-likelihood objective
  weighted_loglik <- function(params, resid, w, dist) {
    omega <- params["omega"]
    alpha <- params[grepl("alpha", names(params))]
    beta <- params[grepl("beta", names(params))]
    
    ## Simple GARCH(p,q) variance recursion
    n <- length(resid)
    sigma2 <- rep(var(resid), n)  # Initialize with unconditional variance
    
    p <- length(alpha)
    q <- length(beta)
    maxpq <- max(p, q)
    
    for (t in (maxpq + 1):n) {
      sigma2[t] <- omega
      for (i in 1:p) {
        sigma2[t] <- sigma2[t] + alpha[i] * resid[t - i]^2
      }
      for (j in 1:q) {
        sigma2[t] <- sigma2[t] + beta[j] * sigma2[t - j]
      }
      ## Ensure positive variance
      if (sigma2[t] <= 0) sigma2[t] <- 1e-10
    }
    
    sig <- sqrt(sigma2)
    
    ## Compute log-likelihood based on distribution
    if (dist == "norm") {
      ll <- dnorm(resid, mean = 0, sd = sig, log = TRUE)
    } else if (dist == "std") {
      shape <- params["shape"]
      if (is.na(shape) || shape <= 2) shape <- 4
      ll <- tsdistributions::dstd(resid, mu = 0, sigma = sig, shape = shape, log = TRUE)
    } else if (dist == "sstd") {
      shape <- params["shape"]
      skew <- params["skew"]
      if (is.na(shape) || shape <= 2) shape <- 4
      if (is.na(skew)) skew <- 1
      ll <- tsdistributions::dsstd(resid, mu = 0, sigma = sig, shape = shape, skew = skew, log = TRUE)
    } else if (dist == "nig") {
      shape <- params["shape"]
      skew <- params["skew"]
      if (is.na(shape) || shape <= 0) shape <- 1
      if (is.na(skew)) skew <- 0
      ll <- tsdistributions::dnig(resid, mu = 0, sigma = sig, shape = shape, skew = skew, log = TRUE)
    } else if (dist == "gh") {
      shape <- params["shape"]
      skew <- params["skew"]
      lambda <- params["lambda"]
      if (is.na(shape) || shape <= 0) shape <- 1
      if (is.na(skew)) skew <- 0
      if (is.na(lambda)) lambda <- -0.5
      ll <- tsdistributions::dgh(resid, mu = 0, sigma = sig, shape = shape, skew = skew, 
                                  lambda = lambda, log = TRUE)
    } else {
      ## Fallback to normal
      ll <- dnorm(resid, mean = 0, sd = sig, log = TRUE)
    }
    
    ll[!is.finite(ll)] <- -1e10
    
    ## Return weighted negative log-likelihood
    return(-sum(w * ll, na.rm = TRUE))
  }
  
  ## Optimize
  warnings_list <- list()
  
  opt_result <- tryCatch({
    withCallingHandlers({
      stats::optim(
        par = unlist(start_pars),
        fn = weighted_loglik,
        lower = lower_bounds,
        upper = upper_bounds,
        method = "L-BFGS-B",
        control = list(ndeps = rep(1e-8, length(start_pars))),
        resid = residuals,
        w = weights,
        dist = spec$distribution %||% "norm"
      )
    }, warning = function(w) {
      warnings_list <<- c(warnings_list, list(w))
      invokeRestart("muffleWarning")
    })
  }, error = function(e) {
    warning("GARCH optimization failed: ", e$message)
    list(par = unlist(start_pars), convergence = 1)
  })
  
  estimated_coeffs <- as.list(opt_result$par)
  names(estimated_coeffs) <- names(start_pars)
  
  return(list(coefficients = estimated_coeffs, warnings = warnings_list))
}


#' Compute GOGARCH Log-Likelihood for MS Framework
#'
#' @description
#' Computes the GOGARCH log-likelihood given estimated GARCH parameters and
#' ICA transformation matrices. This function is used in the E-step of the
#' EM algorithm for Markov-Switching GOGARCH models.
#'
#' @details
#' The GOGARCH log-likelihood consists of two parts:
#'
#' \strong{1. Component Log-Likelihoods}
#'
#' For each independent component \eqn{i} and time \eqn{t}:
#' \deqn{LL_{i,t} = \log f(s_{i,t} | \sigma_{i,t})}
#'
#' \strong{2. Jacobian Adjustment}
#'
#' The ICA transformation introduces a Jacobian term:
#' \deqn{LL_{jacobian} = \log |det(K)|}
#' where \eqn{K} is the pre-whitening matrix from ICA.
#'
#' \strong{Total Log-Likelihood}
#' \deqn{LL = \sum_t \sum_i LL_{i,t} + \log |det(K)|}
#'
#' @param residuals Numeric matrix of residuals with dimensions T x k.
#' @param garch_pars List of GARCH parameters for each component.
#' @param ica_info List containing ICA transformation matrices (A, W, K).
#' @param distribution Character string specifying the component distribution.
#' @param return_vector Logical; if \code{TRUE}, return per-observation
#'   log-likelihoods. Default is \code{FALSE}.
#'
#' @return Scalar total log-likelihood, or vector of per-observation values
#'   if \code{return_vector = TRUE}.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{estimate_garch_weighted_gogarch}}: Estimation function
#'   \item \code{\link{estimate_garch_weighted_univariate_gogarch}}: Component estimation
#' }
#'
#' @keywords internal
compute_gogarch_loglik_ms <- function(
    residuals,
    garch_pars,
    ica_info,
    distribution = "norm",
    return_vector = FALSE
) {
  
  T_obs <- nrow(residuals)
  n_components <- length(garch_pars)
  
  ## Transform residuals to independent components using W
  S <- residuals %*% t(ica_info$W)
  
  ## Compute component-wise log-likelihoods
  ll_matrix <- matrix(0, nrow = T_obs, ncol = n_components)
  
  for (i in 1:n_components) {
    pars <- garch_pars[[i]]
    component <- S[, i]
    
    ## Compute GARCH variance path
    omega <- pars$omega %||% 0.1
    alpha <- unlist(pars[grepl("alpha", names(pars))])
    beta <- unlist(pars[grepl("beta", names(pars))])
    
    if (length(alpha) == 0) alpha <- 0.1
    if (length(beta) == 0) beta <- 0.8
    
    sigma2 <- rep(var(component), T_obs)
    p <- length(alpha)
    q <- length(beta)
    maxpq <- max(p, q, 1)
    
    for (t in (maxpq + 1):T_obs) {
      sigma2[t] <- omega
      for (j in 1:p) {
        sigma2[t] <- sigma2[t] + alpha[j] * component[t - j]^2
      }
      for (j in 1:q) {
        sigma2[t] <- sigma2[t] + beta[j] * sigma2[t - j]
      }
      if (sigma2[t] <= 0) sigma2[t] <- 1e-10
    }
    
    sig <- sqrt(sigma2)
    
    ## Compute log-likelihood
    if (distribution == "norm") {
      ll_matrix[, i] <- dnorm(component, mean = 0, sd = sig, log = TRUE)
    } else {
      ## For other distributions, extract distribution parameters
      shape <- pars$shape %||% 4
      skew <- pars$skew %||% 1
      
      ll_matrix[, i] <- tryCatch({
        switch(distribution,
          "std" = tsdistributions::dstd(component, mu = 0, sigma = sig, shape = shape, log = TRUE),
          "nig" = tsdistributions::dnig(component, mu = 0, sigma = sig, shape = shape, skew = skew, log = TRUE),
          "gh" = tsdistributions::dgh(component, mu = 0, sigma = sig, shape = shape, skew = skew, 
                                       lambda = pars$lambda %||% -0.5, log = TRUE),
          dnorm(component, mean = 0, sd = sig, log = TRUE)
        )
      }, error = function(e) {
        dnorm(component, mean = 0, sd = sig, log = TRUE)
      })
    }
  }
  
  ## Sum across components for each observation
  ll_vector <- rowSums(ll_matrix)
  
  ## Add Jacobian adjustment: log|det(K)|
  K <- ica_info$K
  if (nrow(K) == ncol(K)) {
    jacobian_adj <- log(abs(det(K)))
  } else {
    jacobian_adj <- log(abs(det(K %*% t(K))))
  }
  
  if (return_vector) {
    ## Distribute Jacobian across observations
    return(ll_vector + jacobian_adj / T_obs)
  } else {
    return(sum(ll_vector) + jacobian_adj)
  }
}
