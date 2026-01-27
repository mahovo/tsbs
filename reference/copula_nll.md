# Copula Negative Log-Likelihood for DCC(1,1)

Computes the weighted copula negative log-likelihood

## Usage

``` r
copula_nll(
  params,
  z_matrix,
  weights,
  Qbar,
  copula_dist = "mvn",
  use_reparam = TRUE
)
```

## Arguments

- params:

  Parameter vector (psi, phi) or (psi, phi, shape)

- z_matrix:

  Matrix of copula residuals (T x k)

- weights:

  Observation weights

- Qbar:

  Unconditional covariance matrix

- copula_dist:

  Copula distribution ("mvn" or "mvt")

- use_reparam:

  Logical; if TRUE, params are in reparameterized space

## Value

Scalar negative log-likelihood
