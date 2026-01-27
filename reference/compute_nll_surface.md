# Compute NLL Surface for DCC(1,1)

Compute the negative log-likelihood on a grid of (alpha, beta) values.

## Usage

``` r
compute_nll_surface(
  std_resid,
  weights,
  Qbar,
  alpha_range = c(0.01, 0.2),
  beta_range = c(0.7, 0.98),
  n_grid = 50,
  distribution = "mvn"
)
```

## Arguments

- std_resid:

  T x k matrix of standardized residuals

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

- distribution:

  "mvn" or "mvt"

## Value

List with alpha_grid, beta_grid, nll_surface, and mle location
