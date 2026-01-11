# tsbs: Advanced Time Series Bootstrap Tools

[![R Package](https://img.shields.io/badge/R-Package-blue.svg)](https://github.com/mahovo/tsbs)

## Overview

`tsbs` provides advanced bootstrap methods for time series data, with a focus on capturing complex temporal dependencies including regime-switching behavior and conditional heteroskedasticity. The package implements several bootstrap approaches ranging from classical block methods to sophisticated model-based resampling.

## Installation

Install from GitHub:

```r
# install.packages("devtools")
devtools::install_github("mahovo/tsbs")
```

## Bootstrap Methods

The package supports the following bootstrap types via the main `tsbs()` function:

| Method | `bs_type` | Description |
|--------|-----------|-------------|
| Moving Block | `"moving"` | Fixed-length blocks sampled with replacement |
| Stationary Block | `"stationary"` | Random block lengths (geometric distribution) |
| Hidden Markov Model | `"hmm"` | Regime-switching with discrete states |
| MS-VAR | `"msvar"` | Markov-Switching Vector Autoregression |
| MS-VARMA-GARCH | `"ms_varma_garch"` | Markov-Switching VARMA with GARCH errors |
| Wild Bootstrap | `"wild"` | Residual-based with random sign flips |

## Quick Start

```r
library(tsbs)

# Simple stationary block bootstrap
set.seed(123)
x <- arima.sim(n = 200, list(ar = 0.8))

result <- tsbs(
  x = as.matrix(x),
  bs_type = "stationary",
  block_length = 10,
  num_blocks = 20,
  num_boots = 100,
  func = mean
)

# Bootstrap distribution of the mean
hist(unlist(result$func_outs), main = "Bootstrap Distribution")
abline(v = result$func_out_means, col = "red", lwd = 2)
```

## MS-VARMA-GARCH Bootstrap

The MS-VARMA-GARCH bootstrap is designed for financial time series exhibiting both regime-switching behavior and time-varying volatility. This method:

1. Fits a Markov-Switching model with GARCH errors to the data
2. Simulates new paths from the fitted model
3. Applies user-specified functions to each bootstrap replicate

**Note:** Currently, only the **DCC (Dynamic Conditional Correlation)** specification is implemented for multivariate MS-VARMA-GARCH models.

### Example: Multivariate DCC Model

```r
library(tsbs)

# Bivariate returns data
data <- matrix(rnorm(500 * 2), ncol = 2)
colnames(data) <- c("Asset1", "Asset2")

# Define state-specific specifications
spec <- list(
 # State 1: Low volatility regime
 list(
   var_order = 1,
   garch_spec_fun = "dcc_modelspec",
   distribution = "mvn",
   garch_spec_args = list(
     dcc_order = c(1, 1),
     dynamics = "dcc",
     distribution = "mvn",
     garch_model = list(
       univariate = list(
         list(model = "garch", garch_order = c(1, 1), distribution = "norm"),
         list(model = "garch", garch_order = c(1, 1), distribution = "norm")
       )
     )
   ),
   start_pars = list(
     var_pars = rep(0.1, 6),
     garch_pars = list(
       list(omega = 0.05, alpha1 = 0.08, beta1 = 0.85),
       list(omega = 0.05, alpha1 = 0.08, beta1 = 0.85)
     ),
     dcc_pars = list(alpha_1 = 0.03, beta_1 = 0.94)
   )
 ),
 # State 2: High volatility regime
 list(
   var_order = 1,
   garch_spec_fun = "dcc_modelspec",
   distribution = "mvn",
   garch_spec_args = list(
     dcc_order = c(1, 1),
     dynamics = "dcc",
     distribution = "mvn",
     garch_model = list(
       univariate = list(
         list(model = "garch", garch_order = c(1, 1), distribution = "norm"),
         list(model = "garch", garch_order = c(1, 1), distribution = "norm")
       )
     )
   ),
   start_pars = list(
     var_pars = rep(0.1, 6),
     garch_pars = list(
       list(omega = 0.12, alpha1 = 0.15, beta1 = 0.75),
       list(omega = 0.12, alpha1 = 0.15, beta1 = 0.75)
     ),
     dcc_pars = list(alpha_1 = 0.08, beta_1 = 0.88)
   )
 )
)

# Run bootstrap
result <- tsbs(
 x = data,
 bs_type = "ms_varma_garch",
 model_type = "multivariate",
 num_states = 2,
 spec = spec,
 num_boots = 100,
 control = list(max_iter = 50, tol = 1e-4),
 return_fit = TRUE,
 collect_diagnostics = TRUE
)

# Access model fit and diagnostics
summary(result$fit$diagnostics)
```

## Documentation

### Primary Documentation

The main entry point for the package is the `tsbs()` function, which provides comprehensive documentation including detailed specification examples:

```r
?tsbs
```

### Vignettes

The package includes two vignettes for specialized topics:

| Vignette | Description |
|----------|-------------|
| `vignette("Diagnostics", package = "tsbs")` | Diagnostic system for monitoring EM convergence, parameter evolution, and numerical stability |
| `vignette("dcc_inference_guide", package = "tsbs")` | Statistical inference for DCC parameters including standard errors and hypothesis testing |

### Learning from Test Files

The package test files contain extensive examples of model specification and usage patterns. Examining these files can be a helpful way to understand how to specify models for different scenarios:

```r
# Location of test files (after installation)
system.file("tests", package = "tsbs")
```

Key test files include:
- `test-ms_varma_garch_bs.R` — Comprehensive examples of MS-VARMA-GARCH specifications
- `test-dcc_gradient.R` — DCC gradient and inference testing

### Other resources

See also [https://github.com/mahovo/tsbs/tree/dev/misc](https://github.com/mahovo/tsbs/tree/dev/misc), especially "Portfolio Optimization Demo.Rmd".

### Related Packages

`tsbs` builds on the following packages for GARCH modeling:

- [tsgarch](https://github.com/tsmodels/tsgarch) — Univariate GARCH models
- [tsmarch](https://github.com/tsmodels/tsmarch) — Multivariate GARCH (DCC, Copula)

## Key Features

- **Unified Interface**: Single `tsbs()` function for all bootstrap methods
- **Flexible Function Application**: Apply arbitrary functions to bootstrap replicates
- **Parallel Processing**: Built-in support for parallel computation
- **Comprehensive Diagnostics**: Monitor convergence and numerical stability for complex models
- **Automatic Block Length Selection**: Data-driven heuristics for block bootstrap methods

## Dependencies

Core dependencies:
- R (>= 4.0)
- Rcpp
- tsgarch
- tsmarch

## License
GPL (>= 3)

## Citation

If you use `tsbs` in your research, please cite:

```
@software{tsbs,
 title = {tsbs: Time Series Bootstrap Methods},
 author = {Vognsen, Martin Hoshi},
 url = {https://github.com/mahovo/tsbs}
}
```

## Contributing

Contributions are welcome! Please open an issue or submit a pull request on [GitHub](https://github.com/mahovo/tsbs).