# Fit Standard CGARCH Copula Model

Fits a standard DCC copula model (without ADCC asymmetry).

## Usage

``` r
fit_cgarch_copula(z_matrix, weights, Qbar, copula_dist = "mvn")
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

## Value

List with fitted parameters and convergence info
