# ADCC Copula Negative Log-Likelihood

Computes the weighted copula NLL for ADCC model.

## Usage

``` r
adcc_copula_nll(
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

  Matrix of copula residuals (T x k)

- weights:

  Observation weights

- Qbar:

  Unconditional covariance matrix

- Nbar:

  Unconditional covariance of negative shocks (computed if NULL)

- copula_dist:

  Copula distribution ("mvn" or "mvt")

## Value

Scalar negative log-likelihood
