# Compute NLL Surface for CGARCH(1,1)

Compute the negative log-likelihood on a grid of (alpha, beta) values
for a Copula GARCH model.

## Usage

``` r
compute_cgarch_nll_surface(
  z_matrix,
  weights,
  Qbar,
  alpha_range = c(0.01, 0.2),
  beta_range = c(0.7, 0.98),
  n_grid = 50,
  copula = "mvn",
  shape = 8
)
```

## Arguments

- z_matrix:

  T x k matrix of copula residuals

- weights:

  T-vector of observation weights

- Qbar:

  k x k unconditional correlation matrix

- alpha_range:

  Range for alpha (default c(0.01, 0.20))

- beta_range:

  Range for beta (default c(0.70, 0.98))

- n_grid:

  Number of grid points per dimension (default 50)

- copula:

  Copula distribution ("mvn" or "mvt")

- shape:

  Degrees of freedom for MVT (only used if copula="mvt")

## Value

List with alpha_grid, beta_grid, nll_surface, and mle location
