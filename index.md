# tsbs: Advanced Time Series Bootstrap Tools

[![R
Package](https://img.shields.io/badge/R-Package-blue.svg)](https://github.com/mahovo/tsbs)

## Overview

`tsbs` provides advanced bootstrap methods for time series data, with a
focus on capturing complex temporal dependencies including
regime-switching behavior and conditional heteroskedasticity. The
package implements several bootstrap approaches ranging from classical
block methods to sophisticated model-based resampling.

## Installation

Install from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("mahovo/tsbs")
```

## Bootstrap Methods

The package supports the following bootstrap types via the main
[`tsbs()`](https://mahovo.github.io/tsbs/reference/tsbs.md) function:

| Method              | `bs_type`          | Description                                   |
|---------------------|--------------------|-----------------------------------------------|
| Moving Block        | `"moving"`         | Fixed-length blocks sampled with replacement  |
| Stationary Block    | `"stationary"`     | Random block lengths (geometric distribution) |
| Hidden Markov Model | `"hmm"`            | Regime-switching with discrete states         |
| MS-VAR              | `"msvar"`          | Markov-Switching Vector Autoregression        |
| MS-VARMA-GARCH      | `"ms_varma_garch"` | Markov-Switching VARMA with GARCH errors      |
| Wild Bootstrap      | `"wild"`           | Residual-based with random sign flips         |

## Quick Start

``` r
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

## Documentation

### Primary Documentation

The main entry point for the package is the
[`tsbs()`](https://mahovo.github.io/tsbs/reference/tsbs.md) function,
which provides comprehensive documentation including detailed
specification examples:

[file:///Users/mhvpbp13/git/tsbs/docs/reference/tsbs.html](file:///Users/mhvpbp13/git/tsbs/docs/reference/tsbs.md)

In R with tsbs package installed and loaded:

``` r
?tsbs
```

Full reference manual:  
<https://mahovo.github.io/tsbs/>

### Vignettes

The package includes vignettes for specialized topics:

| Vignette                                                                                                                                 | Description                                                                                   |
|------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------|
| [`vignette("portfolio-optimization-demo", package = "tsbs")`](https://mahovo.github.io/tsbs/articles/portfolio-optimization-demo.md)     | Comprehensive demo of portfolio optimization w/ tsbs                                          |
| [`vignette("bootstrap_diagnostics", package = "tsbs")`](https://mahovo.github.io/tsbs/articles/bootstrap_diagnostics.md)                 | Diagnostic system for output bootstrap series                                                 |
| `vignette("ms-varma-garch_diagnostics.Rmd", package = "tsbs")`                                                                           | Diagnostic system for monitoring EM convergence, parameter evolution, and numerical stability |
| [`vignette("ms-varma-garch_inference", package = "tsbs")`](https://mahovo.github.io/tsbs/articles/ms-varma-garch_inference.md)           | Statistical inference for DCC, CGARCH and GOGARCH parameters                                  |
| [`vignette("multivariate_garch_comparison", package = "tsbs")`](https://mahovo.github.io/tsbs/articles/multivariate_garch_comparison.md) | Diagnostic vignette comparing DCC, CGARCH and GOGARCH models                                  |

### Demos

``` r
library(tsbs)

# See what demos are available
list_package_demos()

# Get path to a demo  (choose demo name from list)
demo_path <- get_demo_path("portfolio-optimization")

# Render it (choose your output location)
demo_output_dir <- "demo-output"
rmarkdown::render(demo_path, output_dir = demo_output_dir)

# Or use in a pipeline
get_demo_path(demo_list[[1]]) |> 
  rmarkdown::render(output_dir = demo_output_dir)
```

See also <https://github.com/mahovo/tsbs/tree/main/inst/demos>

See a pre-rendered version of “portfolio-optimization.Rmd” :  
<https://github.com/mahovo/tsbs/blob/main/inst/demos/portfolio-optimization.md>

### Learning from Test Files

The package test files contain extensive examples of model specification
and usage patterns. Examining these files can be a helpful way to
understand how to specify models for different scenarios:

``` r
# Location of test files (after installation)
system.file("tests", package = "tsbs")
```

Key test files include:

- `test-ms_varma_garch_bs.R` — Comprehensive examples of MS-VARMA-GARCH
  specifications  
- `test-cgarch.R` — Testing the CGARCH integration  
- `test-gogarch.R` — Testing the GOGARCH integration  
- `test-blockBootstrap`  
- `test-hmm_bootstrap`  
- `test-wild_bootstrap`

### Related Packages

`tsbs` builds on the following packages for GARCH modeling:

- [tsgarch](https://github.com/tsmodels/tsgarch) — Univariate GARCH
  models
- [tsmarch](https://github.com/tsmodels/tsmarch) — Multivariate GARCH
  (DCC, Copula)

## Key Features

- **Unified Interface**: Single
  [`tsbs()`](https://mahovo.github.io/tsbs/reference/tsbs.md) function
  for all bootstrap methods
- **Flexible Function Application**: Apply arbitrary functions to
  bootstrap replicates
- **Parallel Processing**: Built-in support for parallel computation
- **Comprehensive Diagnostics**: Monitor convergence and numerical
  stability for complex models
- **Automatic Block Length Selection**: Data-driven heuristics for block
  bootstrap methods

## Dependencies

Core dependencies: - R (\>= 4.0) - Rcpp - tsgarch - tsmarch

## License

GPL (\>= 3)

## Citation

If you use `tsbs` in your research, please cite:

    @software{tsbs,
     title = {tsbs: Time Series Bootstrap Methods},
     author = {Vognsen, Martin Hoshi},
     url = {https://github.com/mahovo/tsbs}
    }

## Contributing

Contributions are welcome! Please open an issue or submit a pull request
on [GitHub](https://github.com/mahovo/tsbs).
