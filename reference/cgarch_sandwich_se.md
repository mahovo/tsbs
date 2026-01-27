# CGARCH Sandwich (Robust) Standard Errors

Computes sandwich estimator for CGARCH parameters.

## Usage

``` r
cgarch_sandwich_se(
  params,
  z_matrix,
  weights,
  Qbar,
  copula_dist = "mvn",
  use_reparam = FALSE
)
```

## Arguments

- params:

  Parameter vector

- z_matrix:

  Copula residuals

- weights:

  Observation weights

- Qbar:

  Unconditional covariance

- copula_dist:

  "mvn" or "mvt"

- use_reparam:

  Use reparameterization?

## Value

List with SE results
