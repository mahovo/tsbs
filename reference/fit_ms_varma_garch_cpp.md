# Fit a general MS-ARMA-GARCH model via the EM Algorithm in C++

Internal C++ function to orchestrate the EM estimation. Delegates all
statistical calculations (likelihood, parameter estimation) to helper
functions in R to interface with the tsgarch and tsmarch packages.

## Usage

``` r
fit_ms_varma_garch_cpp(
  y,
  M,
  spec,
  model_type,
  control,
  diagnostics = NULL,
  verbose = FALSE
)
```

## Arguments

- y:

  A (T x k) matrix of (differenced) time series data.

- M:

  The number of states.

- spec:

  A list of model specifications from R.

- model_type:

  "univariate" or "multivariate".

- control:

  A list with max_iter and tol.

## Value

A raw list with estimated parameters and results.
