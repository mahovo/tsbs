# Plot Diagnostics for MS-VARMA-GARCH Models

Creates diagnostic plots for EM algorithm convergence and parameter
evolution.

## Usage

``` r
# S3 method for class 'ms_diagnostics'
plot(
  x,
  type = c("all", "ll_evolution", "parameters", "sigma"),
  parameters = NULL,
  normalize = FALSE,
  quiet = FALSE,
  ...
)
```

## Arguments

- x:

  An object of class `ms_diagnostics` returned by the fitting procedure.

- type:

  Character string specifying which plots to produce. One of:

  `"all"`

  :   Produce all diagnostic plots (default).

  `"ll_evolution"`

  :   Plot log-likelihood evolution across EM iterations.

  `"parameters"`

  :   Plot parameter evolution across EM iterations.

  `"sigma"`

  :   Plot conditional volatility (sigma) evolution.

- parameters:

  Optional character vector of parameter names to include in the
  parameter evolution plot. If `NULL` (default), all parameters are
  plotted. Supports regex patterns when a single string containing regex
  metacharacters (`^`, `$`, `.`, `*`, `[`) is provided.

- normalize:

  Logical. If `TRUE`, normalize parameter values to \[0, 1\] within each
  parameter for easier comparison across different scales. Default is
  `FALSE`.

- quiet:

  Logical. If `TRUE`, suppress warnings from numeric coercion during
  parameter extraction. Default is `FALSE`.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns the ggplot object(s). Called primarily for side
effects (printing plots).

## Details

The function produces up to three types of diagnostic visualizations:

**Log-Likelihood Evolution** (`type = "ll_evolution"`): Two plots
showing (1) the log-likelihood value after each M-step, and (2) the
change in log-likelihood per iteration. Useful for assessing convergence
and detecting any non-monotonic behavior.

**Parameter Evolution** (`type = "parameters"`): Faceted plot showing
how each parameter evolves across EM iterations, with separate colors
for each regime state. The `parameters` argument can filter to specific
parameters of interest.

**Sigma Evolution** (`type = "sigma"`): Faceted plot showing the mean
conditional volatility (with ± 1 SD ribbon) for each series across
iterations, colored by state.

## See also

[`summary.ms_diagnostics`](https://mahovo.github.io/tsbs/reference/summary.ms_diagnostics.md)
for text summaries of diagnostics.

## Examples

``` r
if (FALSE) { # \dontrun{
# After fitting a model with diagnostics enabled
fit <- fit_ms_varma_garch(data, n_states = 2, collect_diagnostics = TRUE)
diagnostics <- attr(fit, "diagnostics")

# Plot all diagnostics
plot(diagnostics)

# Plot only log-likelihood evolution
plot(diagnostics, type = "ll_evolution")

# Plot specific parameters
plot(diagnostics, type = "parameters", parameters = c("alpha_1", "beta_1"))

# Plot parameters matching a regex pattern
plot(diagnostics, type = "parameters", parameters = "^alpha")

# Normalized parameter plot for cross-parameter comparison
plot(diagnostics, type = "parameters", normalize = TRUE)
} # }
```
