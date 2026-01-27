# GOGARCH Model Diagnostics

Comprehensive diagnostics for GOGARCH model fit.

## Usage

``` r
gogarch_diagnostics(gogarch_result, residuals, verbose = TRUE)
```

## Arguments

- gogarch_result:

  Result from
  [`estimate_garch_weighted_gogarch`](https://mahovo.github.io/tsbs/reference/estimate_garch_weighted_gogarch.md)
  or
  [`estimate_gogarch_improved`](https://mahovo.github.io/tsbs/reference/estimate_gogarch_improved.md)

- residuals:

  Original T x k residual matrix

- verbose:

  Print diagnostic report (default: TRUE)

## Value

List with:

- ica_quality:

  Independence score, negentropy, reconstruction error

- component_diagnostics:

  Per-component: persistence, Ljung-Box p-value

- mixing_matrix:

  Orthogonality error, condition number

- covariance_fit:

  RMSE, correlation vs sample covariance

## Details

Key quality thresholds:

- Independence score \> 0.8

- Reconstruction error \< 1\\

- Mixing matrix condition number \< 100

- Component persistence \< 1 (stationarity)

## See also

[`estimate_garch_weighted_gogarch`](https://mahovo.github.io/tsbs/reference/estimate_garch_weighted_gogarch.md),
[`improved_ica_decomposition`](https://mahovo.github.io/tsbs/reference/improved_ica_decomposition.md)
