# Fit GOGARCH Model

Fits a GOGARCH model to multivariate data.

## Usage

``` r
fit_gogarch_model(y, k, ica_method = "radical", distribution = "norm")
```

## Arguments

- y:

  Matrix of returns (T x k)

- k:

  Number of components

- ica_method:

  ICA algorithm ("radical" or "fastica")

- distribution:

  Component distribution

## Value

List with fitted parameters and convergence info
