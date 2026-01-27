# Bootstrap Standard Errors for GOGARCH Parameters

Compute bootstrap standard errors for GOGARCH component GARCH parameters
using residual resampling.

## Usage

``` r
gogarch_bootstrap_se(
  residuals,
  weights,
  garch_pars,
  ica_info,
  n_boot = 200,
  method = c("residual", "parametric"),
  distribution = "norm",
  verbose = TRUE,
  seed = NULL
)
```

## Arguments

- residuals:

  T x k matrix of observed residuals

- weights:

  T-vector of observation weights

- garch_pars:

  List of GARCH parameters for each component

- ica_info:

  ICA decomposition information

- n_boot:

  Number of bootstrap replications (default 200)

- method:

  Bootstrap method: "residual" (recommended) or "parametric"

- distribution:

  GARCH distribution: "norm" or "std"

- verbose:

  Print progress

- seed:

  Random seed for reproducibility

## Value

List with bootstrap results per component
