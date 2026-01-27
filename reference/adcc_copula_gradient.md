# ADCC Copula Gradient (Numerical)

Computes the gradient of ADCC copula NLL using numerical
differentiation.

## Usage

``` r
adcc_copula_gradient(
  params,
  z_matrix,
  weights,
  Qbar,
  Nbar = NULL,
  copula_dist = "mvn"
)
```

## Arguments

- params:

  Parameter vector: (alpha, gamma, beta) or with shape for MVT

- z_matrix:

  Matrix of copula residuals

- weights:

  Observation weights

- Qbar:

  Unconditional covariance matrix

- Nbar:

  Unconditional covariance of negative shocks

- copula_dist:

  Copula distribution

## Value

Gradient vector
