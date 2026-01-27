# Compute Log-Likelihood with Fixed Parameters

Evaluates the log-likelihood for tsmarch models using fixed
(non-estimated) parameters. This function uses the existing tsmarch
estimation object and replaces the estimated parameters with
user-provided fixed values.

## Usage

``` r
compute_loglik_fixed(
  object,
  params,
  return_components = FALSE,
  ll_vec = FALSE,
  ...
)
```

## Arguments

- object:

  An estimated tsmarch object of class "dcc.estimate",
  "cgarch.estimate", or "gogarch.estimate"

- params:

  A named list of fixed parameters. The names should match the parameter
  names in the model's parmatrix:

  - For DCC models: `list(alpha_1 = 0.05, beta_1 = 0.90)` for DCC(1,1),
    or `list(alpha_1 = 0.05, alpha_2 = 0.03, beta_1 = 0.85)` for
    DCC(2,1). For ADCC models, add `gamma_1`, `gamma_2`, etc.

  - For Student-t DCC: Add `shape = 3` to specify degrees of freedom.

  - For Copula-GARCH models: Same as DCC models above. For Student-t
    copula, use `shape = 3` for degrees of freedom.

  - For GOGARCH models: Parameters for each independent component GARCH
    model, e.g.,
    `list(omega_1 = 0.01, alpha_1 = 0.05, beta_1 = 0.90, omega_2 = 0.02, ...)`.

  Parameter names use underscore notation where the number indicates the
  lag order or component index. Use `coef(object)` to see the exact
  parameter names for your model.

- return_components:

  Logical. If TRUE, returns both the total log-likelihood and its
  components (univariate GARCH + multivariate). Default is FALSE.

- ll_vec:

  Logical. If TRUE, returns per-observation log-likelihood vector
  instead of the total. Cannot be used with return_components = TRUE.
  Default is FALSE.

- ...:

  Additional arguments (currently unused)

## Value

If ll_vec = TRUE, returns a numeric vector of per-observation
log-likelihoods. The vector has length n-1 for DCC models where n is the
number of observations, because the first observation serves as
initialization for the DCC recursion (tsmarch returns a zero placeholder
for this observation which is removed). The sum of the returned vector
equals the scalar log-likelihood returned when ll_vec = FALSE, and also
equals logLik(object) when using estimated parameters.

If ll_vec = FALSE and return_components = FALSE (default), returns a
single numeric value representing the total log-likelihood.

If ll_vec = FALSE and return_components = TRUE, returns a list with
components:

- loglikTotal log-likelihood

- garch_loglikUnivariate GARCH component (DCC/Copula only)

- multivariate_loglikMultivariate component log-likelihood

## Details

This function extracts the specification from an estimated tsmarch
object, replaces the estimated parameters with the user-provided fixed
parameters, and computes the log-likelihood without re-estimating.

The function works by:

1.  Extracting the model specification from the estimated object

2.  Updating the parmatrix with the fixed parameter values

3.  Calling the internal tsmarch likelihood computation functions

## Examples

``` r
# \donttest{
# This example requires tsmarch and tsgarch packages
if (require(tsmarch) && require(tsgarch) && require(xts)) {
  # Generate small sample data
  set.seed(123)
  n <- 500
  returns <- matrix(rnorm(n * 2), ncol = 2)
  returns <- xts(returns, order.by = seq.Date(Sys.Date() - n + 1, 
                                                Sys.Date(), by = "day"))
  colnames(returns) <- c("series1", "series2")
  
  # Estimate univariate GARCH models
  spec1 <- garch_modelspec(returns[,1], model = "garch", order = c(1,1))
  spec2 <- garch_modelspec(returns[,2], model = "garch", order = c(1,1))
  fit1 <- estimate(spec1, keep_tmb = TRUE)
  fit2 <- estimate(spec2, keep_tmb = TRUE)
  
  # Combine into multivariate
  garch_fits <- to_multi_estimate(list(fit1, fit2))
  
  # Estimate DCC model
  dcc_spec <- dcc_modelspec(garch_fits, dynamics = "dcc", dcc_order = c(1,1))
  dcc_fit <- estimate(dcc_spec)
  
  # Get estimated parameters
  est_params <- coef(dcc_fit)
  
  # Compute log-likelihood at estimated parameters
  ll_at_est <- compute_loglik_fixed(dcc_fit, params = as.list(est_params))
  
  # Compute at alternative parameters
  ll_alt <- compute_loglik_fixed(dcc_fit, 
                                  params = list(alpha_1 = 0.05, beta_1 = 0.90))
  
  # The estimated parameters should give higher likelihood
  print(paste("LL at estimated:", ll_at_est))
  print(paste("LL at alternative:", ll_alt))
  print(paste("Difference:", ll_at_est - ll_alt))
}
#> Loading required package: tsmarch
#> Loading required package: tsmethods
#> Loading required package: tsgarch
#> Loading required package: xts
#> Loading required package: zoo
#> 
#> Attaching package: ‘zoo’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.Date, as.Date.numeric
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> [1] "LL at estimated: -1409.00310436415"
#> [1] "LL at alternative: -1414.19905089214"
#> [1] "Difference: 5.19594652799151"
# }

if (FALSE) { # \dontrun{
# Extended example with profile likelihood and LR tests
library(tsmarch)
library(tsgarch)
library(xts)

# Generate sample data with correlation structure
set.seed(100)
n <- 1500
# Create correlated innovations
rho <- 0.6
Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
L <- chol(Sigma)
z <- matrix(rnorm(n * 2), ncol = 2) %*% L

# Add GARCH dynamics
returns <- z
for (i in 2:n) {
  h <- 0.01 + 0.08 * returns[i-1,]^2 + 0.90 * returns[i-1,]^2
  returns[i,] <- returns[i,] * sqrt(pmax(h, 0.001))
}

returns <- xts(returns, order.by = seq(Sys.Date() - n + 1, Sys.Date(), by = "day"))
colnames(returns) <- c("asset1", "asset2")

# Estimate univariate GARCH models
spec1 <- garch_modelspec(returns[,1], model = "garch", order = c(1,1))
spec2 <- garch_modelspec(returns[,2], model = "garch", order = c(1,1))
fit1 <- estimate(spec1, keep_tmb = TRUE)
fit2 <- estimate(spec2, keep_tmb = TRUE)

# Combine into multivariate
garch_fits <- to_multi_estimate(list(fit1, fit2))

# Estimate DCC model
dcc_spec <- dcc_modelspec(garch_fits, dynamics = "dcc", dcc_order = c(1,1))
dcc_fit <- estimate(dcc_spec)

# Get estimated parameters
est_params <- coef(dcc_fit)
cat("Estimated DCC parameters:\n")
print(est_params)

# Compute log-likelihood at estimated parameters
ll_at_estimated <- compute_loglik_fixed(dcc_fit, params = as.list(est_params))
cat("\nLL at estimated parameters:", ll_at_estimated, "\n")

# Test with alternative parameter values
ll_alternative <- compute_loglik_fixed(
  dcc_fit,
  params = list(alpha_1 = 0.03, beta_1 = 0.95)
)
cat("LL at alternative parameters:", ll_alternative, "\n")

# Likelihood ratio test
lr_stat <- 2 * (ll_at_estimated - ll_alternative)
p_value <- pchisq(lr_stat, df = 2, lower.tail = FALSE)
cat("\nLikelihood Ratio Test:\n")
cat("  LR statistic:", round(lr_stat, 4), "\n")
cat("  P-value:", format.pval(p_value, digits = 4), "\n")

# Profile likelihood for alpha
# Only compute if estimated alpha is not at boundary
alpha_est <- est_params["alpha_1"]
beta_est <- est_params["beta_1"]

if (alpha_est > 0.01 && alpha_est < 0.2) {
  cat("\nComputing profile likelihood for alpha...\n")
  
  # Create grid around estimated alpha
  alpha_range <- seq(max(0.001, alpha_est - 0.03), 
                     min(0.3, alpha_est + 0.03), 
                     length.out = 20)
  
  profile_ll <- sapply(alpha_range, function(a) {
    # Skip if would violate stationarity
    if (a + beta_est >= 0.999) return(NA_real_)
    
    tryCatch({
      compute_loglik_fixed(dcc_fit, 
                           params = list(alpha_1 = a, beta_1 = beta_est))
    }, error = function(e) NA_real_)
  })
  
  # Plot profile likelihood
  valid_idx <- !is.na(profile_ll)
  if (sum(valid_idx) > 5) {
    plot(alpha_range[valid_idx], profile_ll[valid_idx], type = "b", 
         xlab = expression(alpha), ylab = "Log-Likelihood",
         main = "Profile Likelihood for Alpha Parameter",
         pch = 19, col = "blue")
    abline(v = alpha_est, col = "red", lty = 2, lwd = 2)
    
    # Add confidence interval based on chi-squared cutoff
    ll_cutoff <- max(profile_ll, na.rm = TRUE) - qchisq(0.95, 1)/2
    abline(h = ll_cutoff, col = "gray", lty = 3)
    
    legend("bottomright", 
           legend = c("Profile LL", "Estimated alpha", "95% CI cutoff"), 
           col = c("blue", "red", "gray"), 
           lty = c(1, 2, 3), pch = c(19, NA, NA), cex = 0.8)
  }
} else {
  cat("\nSkipping profile likelihood (alpha at boundary)\n")
  cat("For a better example, try a different seed or larger sample.\n")
}

# Get component-wise log-likelihoods
ll_components <- compute_loglik_fixed(
  dcc_fit,
  params = as.list(est_params),
  return_components = TRUE
)
cat("\nComponent-wise log-likelihoods:\n")
cat("  Total:", round(ll_components$loglik, 2), "\n")
cat("  GARCH:", round(ll_components$garch_loglik, 2), "\n")
cat("  DCC:  ", round(ll_components$multivariate_loglik, 2), "\n")
cat("  Sum:  ", round(ll_components$garch_loglik + 
                       ll_components$multivariate_loglik, 2), "\n")
} # }
```
