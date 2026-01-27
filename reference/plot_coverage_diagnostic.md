# Plot Coverage Diagnostic

Visualize confidence interval coverage for each replication.

## Usage

``` r
plot_coverage_diagnostic(mc_result, param = "alpha", max_show = 50)
```

## Arguments

- mc_result:

  Result from run_dcc_monte_carlo()

- param:

  Which parameter ("alpha" or "beta")

- max_show:

  Maximum number of replications to show

## Value

plotly object
