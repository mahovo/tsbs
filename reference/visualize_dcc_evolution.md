# Visualize DCC Parameter Evolution from EM Fit

Plot the evolution of DCC parameters across EM iterations, optionally
overlaid on the likelihood surface.

## Usage

``` r
visualize_dcc_evolution(
  fit,
  state = 1,
  show_surface = TRUE,
  true_params = NULL,
  n_grid = 50
)
```

## Arguments

- fit:

  A fitted model from fit_ms_varma_garch() with diagnostics=TRUE

- state:

  Integer state index (default 1)

- show_surface:

  Logical: compute and show likelihood contours? (default TRUE)

- true_params:

  Optional true parameter values c(alpha, beta) for simulation studies

- n_grid:

  Grid size for surface computation (default 50)

## Value

List with trace data, surface (if computed), and plotly visualization
