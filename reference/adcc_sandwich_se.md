# ADCC Sandwich (Robust) Standard Errors

Computes sandwich estimator for ADCC parameters.

## Usage

``` r
adcc_sandwich_se(params, z_matrix, weights, Qbar, Nbar, copula_dist = "mvn")
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

- Nbar:

  Average outer product of negative residuals

- copula_dist:

  "mvn" or "mvt"

## Value

List with SE results
