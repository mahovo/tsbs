# DCC(1,1) Weighted NLL Gradient (Original Parameterization)

Compute gradient of weighted NLL w.r.t. (alpha, beta).

## Usage

``` r
dcc11_gradient_original(
  alpha,
  beta,
  z,
  weights,
  Qbar,
  distribution = "mvn",
  shape = NULL
)
```

## Arguments

- alpha:

  DCC alpha parameter

- beta:

  DCC beta parameter

- z:

  T x k matrix of standardized residuals

- weights:

  T-vector of observation weights

- Qbar:

  k x k unconditional covariance

- distribution:

  "mvn" or "mvt"

- shape:

  Degrees of freedom (MVT only)

## Value

Named vector c(alpha, beta)
