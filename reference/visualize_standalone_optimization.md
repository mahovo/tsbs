# Quick Visualization of DCC Optimization Path (Standalone)

Convenience function for visualizing optimization on simulated data
without needing a full MS-VARMA-GARCH fit. Useful for
testing/diagnostics.

## Usage

``` r
visualize_standalone_optimization(
  n = 500,
  true_alpha = 0.05,
  true_beta = 0.9,
  start_alpha = 0.1,
  start_beta = 0.8,
  omega = c(0.05, 0.08),
  alpha_garch = c(0.1, 0.12),
  beta_garch = c(0.85, 0.82),
  n_grid = 50,
  seed = 42
)
```

## Arguments

- n:

  Number of observations

- true_alpha:

  True DCC alpha

- true_beta:

  True DCC beta

- start_alpha:

  Starting value for alpha

- start_beta:

  Starting value for beta

- omega:

  Vector of GARCH omega parameters

- alpha_garch:

  Vector of GARCH alpha parameters

- beta_garch:

  Vector of GARCH beta parameters

- n_grid:

  Grid size for surface (default 50)

- seed:

  Random seed

## Value

List with trace data, surface, and plot
