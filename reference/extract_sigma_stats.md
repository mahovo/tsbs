# Extract Sigma Evolution Statistics

Extract volatility evolution summary for a specific state and series.

## Usage

``` r
extract_sigma_stats(diagnostics, state, series)
```

## Arguments

- diagnostics:

  An object of class `ms_diagnostics`

- state:

  Integer state index

- series:

  Integer series index

## Value

A data.frame with columns: iteration, mean_sigma, sd_sigma, min_sigma,
max_sigma

## Examples

``` r
if (FALSE) { # \dontrun{
sigma_stats <- extract_sigma_stats(diag, state = 1, series = 1)
plot(sigma_stats$iteration, sigma_stats$mean_sigma, type = "b")
} # }
```
