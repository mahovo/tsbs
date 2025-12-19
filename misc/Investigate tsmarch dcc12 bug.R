## =============================================================================
## Investigate tsmarch DCC(p,q) bugs
## =============================================================================
##
## tsbs primarily uses tsmarch for:
## 
## - Univariate GARCH estimation (via tsgarch)
## - Creating spec objects
## - Some internal functions
# 
## But for the DCC correlation dynamics, tsbs has its own implementation 
## (dcc_recursion, compute_dcc_persistence). This means the tsmarch bug might 
## not actually affect tsbs's higher-order DCC estimation in 
## estimate_garch_weighted_multivariate() because:
## 
## - The param names come from user's start_pars
## - The recursion uses tsbs's own dcc_recursion() function
## - The optimization is done by tsbs, not tsmarch
## 
## However, the bug DOES affect things when:
## 
## Using fit_ms_varma_garch() which creates the spec via tsmarch
## Any code path that relies on tsmarch's parmatrix
##
## Two bugs identified:
##
## BUG 1: DCC(p,q) with p≠q creates wrong number of parameters
##   - DCC(2,1) should have 2 alphas, 1 beta
##   - But tsmarch creates 2 alphas, 2 betas
##
## BUG 2: DCC(1,2) has duplicate parameter names
##   - Should have alpha_1, beta_1, beta_2
##   - But tsmarch creates alpha_1, beta_1, beta_1 (duplicate name)
##
## Standard DCC(p,q) model:
##   Q_t = (1 - sum(alpha) - sum(beta)) * Qbar 
##       + sum_{i=1}^{p} alpha_i * (eps_{t-i} * eps_{t-i}')
##       + sum_{j=1}^{q} beta_j * Q_{t-j}
##
## Where:
##   p = number of alpha parameters (ARCH-like, lagged outer products)
##   q = number of beta parameters (GARCH-like, lagged Q matrices)
##
## =============================================================================

# library(tsmarch)
# library(tsgarch)
# library(xts)

## Create test data
set.seed(123)
n <- 200
dates <- seq(as.Date("2020-01-01"), by = "day", length.out = n)
y_test <- xts(cbind(s1 = rnorm(n), s2 = rnorm(n)), order.by = dates)

## Fit univariate GARCH
uni1 <- estimate(garch_modelspec(y_test[,1], model = "garch", 
                                 order = c(1,1), distribution = "norm"),
                 keep_tmb = TRUE)
uni2 <- estimate(garch_modelspec(y_test[,2], model = "garch", 
                                 order = c(1,1), distribution = "norm"),
                 keep_tmb = TRUE)

multi_est <- to_multi_estimate(list(uni1, uni2))

## =============================================================================
## First, let's directly test .copula_parameters if accessible
## =============================================================================

cat("================================================================================\n")
cat("TESTING .copula_parameters DIRECTLY\n")
cat("================================================================================\n\n")

tryCatch({
  ## Try to access the internal function
  copula_params_fn <- tsmarch:::.copula_parameters
  
  cat("=== .copula_parameters('dcc', 'mvn', c(2,1)) ===\n")
  cat("Expected: 2 alphas, 1 beta\n")
  result <- copula_params_fn("dcc", "mvn", c(2,1))
  print(result[, c("parameter", "value", "estimate")])
  
  cat("\n=== .copula_parameters('dcc', 'mvn', c(1,2)) ===\n")
  cat("Expected: 1 alpha, 2 betas (beta_1, beta_2)\n")
  result <- copula_params_fn("dcc", "mvn", c(1,2))
  print(result[, c("parameter", "value", "estimate")])
  
  cat("\n=== .copula_parameters('dcc', 'mvn', c(3,1)) ===\n")
  cat("Expected: 3 alphas, 1 beta\n")
  result <- copula_params_fn("dcc", "mvn", c(3,1))
  print(result[, c("parameter", "value", "estimate")])
  
}, error = function(e) {
  cat("Could not access .copula_parameters directly:", e$message, "\n")
})


## =============================================================================
## Test via dcc_modelspec and check what order is stored
## =============================================================================

cat("\n================================================================================\n")
cat("CHECKING STORED DCC ORDER IN SPEC\n")
cat("================================================================================\n\n")

cat("=== DCC(2,1) - checking spec$dynamics$order ===\n")
dcc21 <- dcc_modelspec(multi_est, dynamics = "dcc", dcc_order = c(2,1), distribution = "mvn")
cat("Input dcc_order: c(2,1)\n")
cat("Stored spec$dynamics$order:", dcc21$dynamics$order, "\n")
cat("Parmatrix alpha/beta rows:\n")
print(dcc21$parmatrix[grepl("alpha|beta", dcc21$parmatrix$parameter), c("parameter", "value")])

cat("\n=== DCC(1,2) - checking spec$dynamics$order ===\n")
dcc12 <- dcc_modelspec(multi_est, dynamics = "dcc", dcc_order = c(1,2), distribution = "mvn")
cat("Input dcc_order: c(1,2)\n")
cat("Stored spec$dynamics$order:", dcc12$dynamics$order, "\n")
cat("Parmatrix alpha/beta rows:\n")
print(dcc12$parmatrix[grepl("alpha|beta", dcc12$parmatrix$parameter), c("parameter", "value")])


## =============================================================================
## Test various DCC orders and document expected vs actual
## =============================================================================

cat("\n================================================================================\n")
cat("DCC ORDER PARAMETER COUNT ANALYSIS\n")
cat("================================================================================\n\n")

cat("Standard convention: DCC(p,q) where p=#alphas, q=#betas\n\n")

## Helper to count params
count_params <- function(spec) {
  pm <- spec$parmatrix
  n_alpha <- sum(grepl("^alpha_", pm$parameter))
  n_beta <- sum(grepl("^beta_", pm$parameter))
  beta_names <- pm$parameter[grepl("^beta_", pm$parameter)]
  alpha_names <- pm$parameter[grepl("^alpha_", pm$parameter)]
  list(n_alpha = n_alpha, n_beta = n_beta, 
       alpha_names = alpha_names, beta_names = beta_names)
}

cat("=== DCC(1,1) ===\n")
cat("Expected: 1 alpha, 1 beta\n")
dcc11 <- dcc_modelspec(multi_est, dynamics = "dcc", dcc_order = c(1,1), distribution = "mvn")
p11 <- count_params(dcc11)
cat(sprintf("Actual:   %d alpha, %d beta\n", p11$n_alpha, p11$n_beta))
cat("Alpha names:", paste(p11$alpha_names, collapse = ", "), "\n")
cat("Beta names:", paste(p11$beta_names, collapse = ", "), "\n")
cat("Status: ", if(p11$n_alpha == 1 && p11$n_beta == 1) "OK" else "BUG", "\n\n")

cat("=== DCC(2,1) ===\n")
cat("Expected: 2 alphas, 1 beta\n")
dcc21 <- dcc_modelspec(multi_est, dynamics = "dcc", dcc_order = c(2,1), distribution = "mvn")
p21 <- count_params(dcc21)
cat(sprintf("Actual:   %d alpha, %d beta\n", p21$n_alpha, p21$n_beta))
cat("Alpha names:", paste(p21$alpha_names, collapse = ", "), "\n")
cat("Beta names:", paste(p21$beta_names, collapse = ", "), "\n")
cat("Status: ", if(p21$n_alpha == 2 && p21$n_beta == 1) "OK" else "BUG - wrong parameter count", "\n\n")

cat("=== DCC(1,2) ===\n")
cat("Expected: 1 alpha, 2 betas (named beta_1, beta_2)\n")
dcc12 <- dcc_modelspec(multi_est, dynamics = "dcc", dcc_order = c(1,2), distribution = "mvn")
p12 <- count_params(dcc12)
cat(sprintf("Actual:   %d alpha, %d beta\n", p12$n_alpha, p12$n_beta))
cat("Alpha names:", paste(p12$alpha_names, collapse = ", "), "\n")
cat("Beta names:", paste(p12$beta_names, collapse = ", "), "\n")
has_duplicate <- length(p12$beta_names) != length(unique(p12$beta_names))
cat("Duplicate names:", has_duplicate, "\n")
cat("Status: ", if(p12$n_alpha == 1 && p12$n_beta == 2 && !has_duplicate) "OK" else "BUG", "\n\n")

cat("=== DCC(2,2) ===\n")
cat("Expected: 2 alphas, 2 betas\n")
dcc22 <- dcc_modelspec(multi_est, dynamics = "dcc", dcc_order = c(2,2), distribution = "mvn")
p22 <- count_params(dcc22)
cat(sprintf("Actual:   %d alpha, %d beta\n", p22$n_alpha, p22$n_beta))
cat("Alpha names:", paste(p22$alpha_names, collapse = ", "), "\n")
cat("Beta names:", paste(p22$beta_names, collapse = ", "), "\n")
cat("Status: ", if(p22$n_alpha == 2 && p22$n_beta == 2) "OK" else "BUG", "\n\n")

cat("=== DCC(3,1) ===\n")
cat("Expected: 3 alphas, 1 beta\n")
dcc31 <- dcc_modelspec(multi_est, dynamics = "dcc", dcc_order = c(3,1), distribution = "mvn")
p31 <- count_params(dcc31)
cat(sprintf("Actual:   %d alpha, %d beta\n", p31$n_alpha, p31$n_beta))
cat("Alpha names:", paste(p31$alpha_names, collapse = ", "), "\n")
cat("Beta names:", paste(p31$beta_names, collapse = ", "), "\n")
cat("Status: ", if(p31$n_alpha == 3 && p31$n_beta == 1) "OK" else "BUG - wrong parameter count", "\n\n")

cat("=== DCC(1,3) ===\n")
cat("Expected: 1 alpha, 3 betas\n")
dcc13 <- dcc_modelspec(multi_est, dynamics = "dcc", dcc_order = c(1,3), distribution = "mvn")
p13 <- count_params(dcc13)
cat(sprintf("Actual:   %d alpha, %d beta\n", p13$n_alpha, p13$n_beta))
cat("Alpha names:", paste(p13$alpha_names, collapse = ", "), "\n")
cat("Beta names:", paste(p13$beta_names, collapse = ", "), "\n")
has_duplicate13 <- length(p13$beta_names) != length(unique(p13$beta_names))
cat("Duplicate names:", has_duplicate13, "\n")
cat("Status: ", if(p13$n_alpha == 1 && p13$n_beta == 3 && !has_duplicate13) "OK" else "BUG", "\n\n")


## =============================================================================
## Full parmatrix output for bug report
## =============================================================================

cat("\n================================================================================\n")
cat("FULL PARMATRIX OUTPUT\n")
cat("================================================================================\n\n")

cat("=== DCC(2,1) full parmatrix ===\n")
print(dcc21$parmatrix[, c("parameter", "lower", "upper", "value", "estimate")])

cat("\n=== DCC(1,2) full parmatrix ===\n")
print(dcc12$parmatrix[, c("parameter", "lower", "upper", "value", "estimate")])


## =============================================================================
## Minimal reproducible example for bug report
## =============================================================================

cat("\n\n================================================================================\n")
cat("MINIMAL REPRODUCIBLE EXAMPLE FOR BUG REPORT\n")
cat("================================================================================\n")
cat('
```r
library(tsmarch)
library(tsgarch)
library(xts)

set.seed(123)
n <- 200
dates <- seq(as.Date("2020-01-01"), by = "day", length.out = n)
y <- xts(cbind(s1 = rnorm(n), s2 = rnorm(n)), order.by = dates)

uni1 <- estimate(garch_modelspec(y[,1], model = "garch", order = c(1,1), 
                                  distribution = "norm"), keep_tmb = TRUE)
uni2 <- estimate(garch_modelspec(y[,2], model = "garch", order = c(1,1), 
                                  distribution = "norm"), keep_tmb = TRUE)
multi_est <- to_multi_estimate(list(uni1, uni2))

# BUG 1: DCC(2,1) creates wrong number of parameters
# Expected: 2 alphas, 1 beta
# Actual: 2 alphas, 2 betas
dcc21 <- dcc_modelspec(multi_est, dynamics = "dcc", dcc_order = c(2,1), 
                        distribution = "mvn")
print(dcc21$parmatrix[grepl("alpha|beta", dcc21$parmatrix$parameter), 
                       c("parameter", "value")])
# Shows: alpha_1, alpha_2, beta_1, beta_2 (should only have beta_1)

# BUG 2: DCC(1,2) has duplicate parameter names
# Expected: alpha_1, beta_1, beta_2
# Actual: alpha_1, beta_1, beta_1 (duplicate name)
dcc12 <- dcc_modelspec(multi_est, dynamics = "dcc", dcc_order = c(1,2), 
                        distribution = "mvn")
print(dcc12$parmatrix[grepl("alpha|beta", dcc12$parmatrix$parameter), 
                       c("parameter", "value")])
# Shows: alpha_1, beta_1, beta_1 (beta_1 appears twice)
```
')


## =============================================================================
## Session info for bug report
## =============================================================================

cat("\n\n================================================================================\n")
cat("SESSION INFO\n")
cat("================================================================================\n\n")
print(sessionInfo())


## =============================================================================
## Summary for bug report
## =============================================================================

cat("\n\n================================================================================\n")
cat("SUMMARY FOR BUG REPORT\n")
cat("================================================================================\n")
cat("
Title: DCC(p,q) parmatrix has incorrect parameter count and naming when p≠q

Description:
When creating a DCC specification with dcc_order = c(p,q) where p≠q, the parmatrix 
has incorrect parameters:

1. DCC(2,1) should create 2 alpha parameters and 1 beta parameter, but creates 
   2 alphas and 2 betas.

2. DCC(1,2) should create 1 alpha parameter and 2 beta parameters (beta_1, beta_2),
   but creates 1 alpha and 2 betas both named 'beta_1' (duplicate name).

The issue appears to be in how the parmatrix is constructed for asymmetric orders.
DCC(1,1) and DCC(2,2) work correctly.

This affects:
- Parameter estimation (wrong number of parameters)
- Coefficient extraction by name (duplicate names cause issues)

Expected behavior:
- dcc_order = c(p,q) should create exactly p alpha parameters and q beta parameters
- Beta parameters should be named beta_1, beta_2, ..., beta_q

")