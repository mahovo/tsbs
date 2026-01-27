# Compute Profile Likelihood Confidence Intervals

Compute Profile Likelihood Confidence Intervals

## Usage

``` r
dcc11_profile_ci(
  params,
  std_resid,
  weights,
  Qbar,
  distribution = "mvn",
  use_reparam = FALSE,
  param_idx = 1,
  level = 0.95,
  n_grid = 50
)
```

## Arguments

- params:

  MLE estimates

- std_resid:

  Standardized residuals

- weights:

  Observation weights

- Qbar:

  Unconditional covariance

- distribution:

  "mvn" or "mvt"

- use_reparam:

  Logical

- param_idx:

  Index of parameter for CI

- level:

  Confidence level (default 0.95)

- n_grid:

  Grid points (default 50)

## Value

Named vector c(lower, upper)
