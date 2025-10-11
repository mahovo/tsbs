# Helper file to create test fixtures
# Run this once to generate test data, then save the objects

#' Create Test Data for tsmarch Models
#'
#' @description Creates synthetic multivariate return data suitable for 
#' testing tsmarch models
#' @param n Number of observations
#' @param k Number of series
#' @param seed Random seed for reproducibility
#' @return An xts object with k columns and n rows
create_test_data <- function(n = 500, k = 2, seed = 123) {
  if (!requireNamespace("xts", quietly = TRUE)) {
    stop("Package 'xts' is required for test data creation")
  }
  
  set.seed(seed)
  
  # Generate correlated returns using a simple DGP
  # True parameters: alpha = 0.05, beta = 0.90
  
  # Generate standardized innovations
  z <- matrix(rnorm(n * k), ncol = k)
  
  # Generate volatilities (GARCH(1,1) with omega=0.01, alpha=0.08, beta=0.90)
  sigma <- matrix(0, nrow = n, ncol = k)
  epsilon <- matrix(0, nrow = n, ncol = k)
  
  for (j in 1:k) {
    sigma[1, j] <- sqrt(0.01 / (1 - 0.08 - 0.90))  # Unconditional variance
    epsilon[1, j] <- sigma[1, j] * z[1, j]
    
    for (i in 2:n) {
      sigma[i, j] <- sqrt(0.01 + 0.08 * epsilon[i-1, j]^2 + 0.90 * sigma[i-1, j]^2)
      epsilon[i, j] <- sigma[i, j] * z[i, j]
    }
  }
  
  # Create xts object
  dates <- seq(Sys.Date() - n + 1, Sys.Date(), by = "day")
  returns <- xts::xts(epsilon, order.by = dates)
  colnames(returns) <- paste0("series", 1:k)
  
  return(returns)
}

#' Get fixtures directory path
#' @keywords internal
.get_fixtures_dir <- function() {
  # Try different locations
  if (file.exists("tests/testthat/fixtures")) {
    return("tests/testthat/fixtures")
  } else if (file.exists("fixtures")) {
    return("fixtures")
  } else if (exists("test_path")) {
    # In test context
    return(test_path("fixtures"))
  } else {
    # Create if doesn't exist
    dir.create("tests/testthat/fixtures", recursive = TRUE, showWarnings = FALSE)
    return("tests/testthat/fixtures")
  }
}

#' Create a Simple DCC Model Estimate Object for Testing
#'
#' @description Creates a minimal DCC estimate object that can be used in tests.
#' This function estimates an actual model on synthetic data.
#' @param use_cache If TRUE, tries to load from cache; if FALSE, re-estimates
#' @return A dcc.estimate object
create_test_dcc_dynamic <- function(use_cache = TRUE) {
  
  cache_file <- file.path(.get_fixtures_dir(), "dcc_dynamic_fit.rds")
  
  # Try to load from cache
  if (use_cache && file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  
  # Check required packages
  required_packages <- c("tsmarch", "tsgarch", "xts")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      skip(paste("Package", pkg, "not available"))
    }
  }
  
  # Create test data
  returns <- create_test_data(n = 500, k = 2, seed = 123)
  
  # Estimate univariate GARCH models (keep TMB objects)
  spec1 <- garch_modelspec(returns[,1], model = "garch", order = c(1,1),
                           distribution = "norm", constant = TRUE)
  spec2 <- garch_modelspec(returns[,2], model = "garch", order = c(1,1),
                           distribution = "norm", constant = TRUE)
  
  fit1 <- estimate(spec1, keep_tmb = TRUE)
  fit2 <- estimate(spec2, keep_tmb = TRUE)
  
  # Combine into multivariate - use the proper method
  garch_fits <- list(fit1, fit2)
  names(garch_fits) <- colnames(returns)
  
  # Convert to multi_estimate using tsgarch's function
  garch_fits <- to_multi_estimate(garch_fits)
  
  # Estimate DCC model (keep hessian and TMB for completeness)
  dcc_spec <- dcc_modelspec(garch_fits, dynamics = "dcc", 
                            dcc_order = c(1,1), distribution = "mvn")
  dcc_fit <- estimate(dcc_spec, return_hessian = TRUE)
  
  # Save to cache
  if (use_cache) {
    dir.create(.get_fixtures_dir(), showWarnings = FALSE, recursive = TRUE)
    saveRDS(dcc_fit, cache_file)
  }
  
  return(dcc_fit)
}

#' Create a DCC Constant Model for Testing
create_test_dcc_constant <- function(use_cache = TRUE) {
  
  cache_file <- file.path(.get_fixtures_dir(), "dcc_constant_fit.rds")
  
  if (use_cache && file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  
  # Similar to above but with constant correlation
  required_packages <- c("tsmarch", "tsgarch", "xts")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      skip(paste("Package", pkg, "not available"))
    }
  }
  
  returns <- create_test_data(n = 500, k = 2, seed = 124)
  
  spec1 <- garch_modelspec(returns[,1], model = "garch", order = c(1,1),
                           distribution = "norm")
  spec2 <- garch_modelspec(returns[,2], model = "garch", order = c(1,1),
                           distribution = "norm")
  
  fit1 <- estimate(spec1, keep_tmb = TRUE)
  fit2 <- estimate(spec2, keep_tmb = TRUE)
  
  garch_fits <- list(fit1, fit2)
  names(garch_fits) <- colnames(returns)
  garch_fits <- to_multi_estimate(garch_fits)
  
  # Constant correlation
  dcc_spec <- dcc_modelspec(garch_fits, dynamics = "constant", 
                            distribution = "mvn")
  dcc_fit <- estimate(dcc_spec, return_hessian = TRUE)
  
  if (use_cache) {
    dir.create(.get_fixtures_dir(), showWarnings = FALSE, recursive = TRUE)
    saveRDS(dcc_fit, cache_file)
  }
  
  return(dcc_fit)
}

#' Create a Copula-GARCH Dynamic Model for Testing
create_test_copula_dynamic <- function(use_cache = TRUE) {
  
  cache_file <- file.path(.get_fixtures_dir(), "copula_dynamic_fit.rds")
  
  if (use_cache && file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  
  required_packages <- c("tsmarch", "tsgarch", "xts")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      skip(paste("Package", pkg, "not available"))
    }
  }
  
  returns <- create_test_data(n = 500, k = 2, seed = 125)
  
  spec1 <- garch_modelspec(returns[,1], model = "garch", order = c(1,1),
                           distribution = "norm")
  spec2 <- garch_modelspec(returns[,2], model = "garch", order = c(1,1),
                           distribution = "norm")
  
  fit1 <- estimate(spec1, keep_tmb = TRUE)
  fit2 <- estimate(spec2, keep_tmb = TRUE)
  
  garch_fits <- list(fit1, fit2)
  names(garch_fits) <- colnames(returns)
  garch_fits <- to_multi_estimate(garch_fits)
  
  # Copula-GARCH with parametric transformation
  cgarch_spec <- cgarch_modelspec(
    garch_fits, 
    dynamics = "dcc",
    dcc_order = c(1,1),
    copula = "mvn",
    transformation = "parametric"
  )
  cgarch_fit <- estimate(cgarch_spec, return_hessian = TRUE)
  
  if (use_cache) {
    dir.create(.get_fixtures_dir(), showWarnings = FALSE, recursive = TRUE)
    saveRDS(cgarch_fit, cache_file)
  }
  
  return(cgarch_fit)
}

#' Create a Student-t DCC Model for Testing
create_test_dcc_student <- function(use_cache = TRUE) {
  
  cache_file <- file.path(.get_fixtures_dir(), "dcc_student_fit.rds")
  
  if (use_cache && file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  
  required_packages <- c("tsmarch", "tsgarch", "xts")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      skip(paste("Package", pkg, "not available"))
    }
  }
  
  returns <- create_test_data(n = 500, k = 2, seed = 126)
  
  # Use normal for univariate when multivariate is Student-t
  spec1 <- garch_modelspec(returns[,1], model = "garch", order = c(1,1),
                           distribution = "norm")
  spec2 <- garch_modelspec(returns[,2], model = "garch", order = c(1,1),
                           distribution = "norm")
  
  fit1 <- estimate(spec1, keep_tmb = TRUE)
  fit2 <- estimate(spec2, keep_tmb = TRUE)
  
  garch_fits <- list(fit1, fit2)
  names(garch_fits) <- colnames(returns)
  garch_fits <- to_multi_estimate(garch_fits)
  
  # DCC with Student-t distribution
  dcc_spec <- dcc_modelspec(garch_fits, dynamics = "dcc", 
                            dcc_order = c(1,1), distribution = "mvt")
  dcc_fit <- estimate(dcc_spec, return_hessian = TRUE)
  
  if (use_cache) {
    dir.create(.get_fixtures_dir(), showWarnings = FALSE, recursive = TRUE)
    saveRDS(dcc_fit, cache_file)
  }
  
  return(dcc_fit)
}