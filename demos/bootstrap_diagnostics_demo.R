# =============================================================================
# Bootstrap Diagnostics Demo for tsbs Package
# =============================================================================
#
# This file demonstrates the bootstrap diagnostics system integrated into tsbs()
#

# Source the diagnostics file (would normally be part of the package)
# source("bootstrap_diagnostics.R")

library(tsbs)

# =============================================================================
# Example 1: Moving Block Bootstrap with Diagnostics
# =============================================================================

set.seed(123)

# Generate sample AR(1) data
n <- 200
x_ar1 <- arima.sim(n = n, list(ar = 0.7))

# Run moving block bootstrap with diagnostics
result_mbb <- tsbs(
  x = as.matrix(x_ar1),
  block_length = 10,
  bs_type = "moving",
  block_type = "overlapping",
  n_boot = 200,
  num_boots = 100,
  func = mean,
  return_diagnostics = TRUE
)

# View diagnostics summary
cat("\n=== MOVING BLOCK BOOTSTRAP DIAGNOSTICS ===\n")
summary(result_mbb$diagnostics)

# Plot diagnostics
cat("\n--- Plotting block length distribution ---\n")
plot(result_mbb$diagnostics, type = "block_lengths")

cat("\n--- Plotting means comparison ---\n")
plot(result_mbb$diagnostics, type = "means_comparison")

# =============================================================================
# Example 2: Stationary Bootstrap with Diagnostics
# =============================================================================

set.seed(456)

# Generate more persistent AR(1) data
x_ar_persistent <- arima.sim(n = n, list(ar = 0.9))

# Run stationary bootstrap with diagnostics
result_sb <- tsbs(
  x = as.matrix(x_ar_persistent),
  bs_type = "stationary",
  p_method = "plugin",  # Auto-select p based on autocorrelation
  n_boot = 200,
  num_boots = 100,
  func = mean,
  return_diagnostics = TRUE
)

cat("\n=== STATIONARY BOOTSTRAP DIAGNOSTICS ===\n")
summary(result_sb$diagnostics)

# The stationary bootstrap has geometrically distributed block lengths
cat("\n--- Block length distribution (should be geometric) ---\n")
plot(result_sb$diagnostics, type = "block_lengths")

# =============================================================================
# Example 3: Multivariate Bootstrap with Diagnostics
# =============================================================================

set.seed(789)

# Generate correlated multivariate data
n <- 150
Sigma <- matrix(c(1, 0.6, 0.6, 1), nrow = 2)
L <- t(chol(Sigma))
innovations <- matrix(rnorm(n * 2), ncol = 2)
x_mv <- innovations %*% t(L)

# Add some persistence
for (i in 2:n) {
  x_mv[i, ] <- 0.5 * x_mv[i-1, ] + 0.5 * x_mv[i, ]
}

# Run multivariate bootstrap with diagnostics
result_mv <- tsbs(
  x = x_mv,
  block_length = 8,
  bs_type = "moving",
  n_boot = 150,
  num_boots = 50,
  func = function(m) c(mean1 = mean(m[,1]), mean2 = mean(m[,2]), cor = cor(m[,1], m[,2])),
  apply_func_to = "all",
  return_diagnostics = TRUE
)

cat("\n=== MULTIVARIATE BOOTSTRAP DIAGNOSTICS ===\n")
summary(result_mv$diagnostics)

# Compare original vs bootstrap autocorrelations
cat("\n--- Autocorrelation comparison ---\n")
plot(result_mv$diagnostics, type = "acf_comparison")

# =============================================================================
# Example 4: Wild Bootstrap with Diagnostics
# =============================================================================

set.seed(321)

# Generate residuals with heteroskedasticity
n <- 100
sigma_t <- 1 + 0.5 * sin(2 * pi * (1:n) / 50)  # Time-varying volatility
residuals <- rnorm(n) * sigma_t

# Run wild bootstrap with diagnostics
result_wild <- tsbs(
  x = as.matrix(residuals),
  bs_type = "wild",
  num_boots = 200,
  return_diagnostics = TRUE
)

cat("\n=== WILD BOOTSTRAP DIAGNOSTICS ===\n")
summary(result_wild$diagnostics)

# Wild bootstrap doesn't have blocks, but we can still compare distributions
plot(result_wild$diagnostics, type = "means_comparison")

# =============================================================================
# Example 5: Extracting Diagnostic Data for Custom Analysis
# =============================================================================

cat("\n=== EXTRACTING DIAGNOSTIC DATA ===\n")

# Get block information as data frame
blocks_df <- extract_blocks(result_mbb$diagnostics)
cat("\nFirst few blocks:\n")
print(head(blocks_df))

# Get summary statistics
stats <- extract_summary_stats(result_mbb$diagnostics)
cat("\nOriginal mean:", stats$original$means, "\n")
cat("Bootstrap mean (avg):", stats$bootstrap$means$mean, "\n")
cat("Bootstrap mean (95% CI):", stats$bootstrap$means$quantiles[c("2.5%", "97.5%")], "\n")

# Convert to data frame for custom analysis
stats_df <- as.data.frame(result_mbb$diagnostics, what = "stats")
cat("\nReplicate statistics (first 5 rows):\n")
print(head(stats_df, 5))

# =============================================================================
# Example 6: Tapered Block Bootstrap with Diagnostics
# =============================================================================

set.seed(111)

# Run tapered block bootstrap
result_tapered <- tsbs(
  x = as.matrix(x_ar1),
  block_length = 12,
  bs_type = "moving",
  block_type = "tapered",
  taper_type = "tukey",
  tukey_alpha = 0.5,
  n_boot = 200,
  num_boots = 50,
  return_diagnostics = TRUE
)

cat("\n=== TAPERED BLOCK BOOTSTRAP DIAGNOSTICS ===\n")
summary(result_tapered$diagnostics)

# =============================================================================
# Example 7: Variable-Length Bootstrap Series (no n_boot specified)
# =============================================================================

set.seed(222)

# When n_boot is not specified, stationary bootstrap creates variable-length series
result_varlen <- tsbs(
  x = as.matrix(x_ar1),
  bs_type = "stationary",
  num_blocks = 15,  # Specify number of blocks instead of length
  num_boots = 30,
  return_diagnostics = TRUE
)

cat("\n=== VARIABLE-LENGTH BOOTSTRAP DIAGNOSTICS ===\n")
summary(result_varlen$diagnostics)

# Plot the distribution of series lengths
plot(result_varlen$diagnostics, type = "length_distribution")

cat("\n=== DEMONSTRATION COMPLETE ===\n")



