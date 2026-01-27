# Check for Sigma Update Failures

Identify cases where volatility (sigma) failed to update during
estimation.

## Usage

``` r
check_sigma_updates(diagnostics)
```

## Arguments

- diagnostics:

  An object of class `ms_diagnostics`

## Value

List with `has_failures` (logical) and details of failed updates
