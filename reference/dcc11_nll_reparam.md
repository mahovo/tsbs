# Reparameterized DCC(1,1) negative log-likelihood

This function computes the weighted negative log-likelihood for DCC(1,1)
using the (persistence, ratio) parameterization. This eliminates the
need for penalty-based stationarity enforcement since persistence \< 1
is guaranteed by the box constraints.

## Usage

``` r
dcc11_nll_reparam(
  reparam_pars,
  std_resid,
  weights,
  Qbar,
  distribution = "mvn",
  dist_pars = NULL,
  verbose = FALSE
)
```

## Arguments

- reparam_pars:

  Numeric vector c(persistence, ratio) or named list

- std_resid:

  Matrix of standardized residuals (T x k)

- weights:

  Observation weights (length T)

- Qbar:

  Unconditional covariance matrix (k x k)

- distribution:

  Character; either "mvn" or "mvt"

- dist_pars:

  Distribution parameters (e.g., shape for mvt)

- verbose:

  Logical; if TRUE, print diagnostic messages

## Value

Negative log-likelihood value (scalar)
