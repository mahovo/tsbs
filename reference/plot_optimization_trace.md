# Plot Parameter Trace from Optimization

Visualize the path of parameter estimates during optimization overlaid
on the likelihood contours.

## Usage

``` r
plot_optimization_trace(
  trace_data,
  surface = NULL,
  true_params = NULL,
  title = "Optimization Trace"
)
```

## Arguments

- trace_data:

  Data frame with columns: iteration, alpha, beta, nll

- surface:

  Result from compute_nll_surface() (optional)

- true_params:

  True parameter values (optional)

- title:

  Plot title

## Value

plotly object
