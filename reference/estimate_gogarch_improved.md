# Estimate GOGARCH with Improved ICA

Alternative to
[`estimate_garch_weighted_gogarch`](https://mahovo.github.io/tsbs/reference/estimate_garch_weighted_gogarch.md)
using
[`improved_ica_decomposition`](https://mahovo.github.io/tsbs/reference/improved_ica_decomposition.md)
for better convergence.

## Usage

``` r
estimate_gogarch_improved(
  residuals,
  weights,
  spec,
  ica_restarts = 3,
  diagnostics = NULL,
  verbose = FALSE
)
```

## Arguments

- residuals:

  T x k matrix of residuals

- weights:

  T-vector of observation weights

- spec:

  Model specification (see
  [`estimate_garch_weighted_gogarch`](https://mahovo.github.io/tsbs/reference/estimate_garch_weighted_gogarch.md))

- ica_restarts:

  Number of ICA restarts (default: 3)

- diagnostics:

  Optional diagnostics collector

- verbose:

  Print progress

## Value

Same structure as
[`estimate_garch_weighted_gogarch`](https://mahovo.github.io/tsbs/reference/estimate_garch_weighted_gogarch.md)
with ICA quality metrics in `coefficients$ica_info$quality`

## Details

Use this when ICA convergence is problematic. Provides warnings when ICA
quality is poor (independence score \< 0.8) and falls back to PCA if ICA
fails completely.

## See also

[`estimate_garch_weighted_gogarch`](https://mahovo.github.io/tsbs/reference/estimate_garch_weighted_gogarch.md),
[`improved_ica_decomposition`](https://mahovo.github.io/tsbs/reference/improved_ica_decomposition.md),
[`gogarch_diagnostics`](https://mahovo.github.io/tsbs/reference/gogarch_diagnostics.md)
