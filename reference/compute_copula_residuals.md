# Compute Copula Residuals from Uniform Margins

Transforms uniform margins to copula residuals (Z) based on the copula
distribution

## Usage

``` r
compute_copula_residuals(u_matrix, copula_dist = "mvn", dist_pars = NULL)
```

## Arguments

- u_matrix:

  Matrix of uniform margins (T x k)

- copula_dist:

  Copula distribution ("mvn" or "mvt")

- dist_pars:

  Distribution parameters (shape for mvt)

## Value

Matrix of copula residuals (T x k)
