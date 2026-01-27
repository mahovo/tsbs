# Copula Gradient in Original (alpha, beta) Space

Computes analytical gradient of copula NLL w.r.t. (alpha, beta)

## Usage

``` r
copula_gradient_original(
  alpha,
  beta,
  z,
  weights,
  Qbar,
  copula_dist = "mvn",
  shape = NULL
)
```

## Arguments

- alpha:

  DCC alpha parameter

- beta:

  DCC beta parameter

- z:

  T x k matrix of copula residuals

- weights:

  T-vector of observation weights

- Qbar:

  k x k unconditional covariance

- copula_dist:

  "mvn" or "mvt"

- shape:

  Degrees of freedom (MVT only)

## Value

Named vector c(alpha, beta)
