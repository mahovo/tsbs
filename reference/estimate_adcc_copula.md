# Estimate ADCC Copula Parameters

Estimates ADCC parameters (alpha, gamma, beta) for copula model.

## Usage

``` r
estimate_adcc_copula(
  z_matrix,
  weights,
  Qbar,
  copula_dist = "mvn",
  start_pars = NULL
)
```

## Arguments

- z_matrix:

  Matrix of copula residuals (T x k)

- weights:

  Observation weights

- Qbar:

  Unconditional covariance matrix

- copula_dist:

  Copula distribution ("mvn" or "mvt")

- start_pars:

  Optional starting values

## Value

List with alpha, gamma, beta, shape (if MVT), nll, convergence
