# Visualize GOGARCH Component Evolution

Plot the evolution of GOGARCH component parameters across EM iterations.

## Usage

``` r
visualize_gogarch_evolution(diagnostics, state = 1, show_persistence = TRUE)
```

## Arguments

- diagnostics:

  An object of class `ms_diagnostics`

- state:

  Integer state index (default 1)

- show_persistence:

  Logical: show persistence rather than alpha/beta separately

## Value

plotly object with component trajectory plot
