# Fit ADCC Copula Model

Fits an ADCC (Asymmetric DCC) copula model to standardized residuals.
Wraps estimate_adcc_copula if available, otherwise provides standalone
implementation.

## Usage

``` r
fit_adcc_copula(z_matrix, weights, Qbar, copula_dist = "mvn", shape = NULL)
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

- shape:

  Fixed shape parameter (if NULL and copula_dist="mvt", estimated)

## Value

List with fitted parameters and convergence info
