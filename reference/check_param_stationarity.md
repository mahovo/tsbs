# Check Parameter Stationarity

Verify that GARCH and DCC parameters satisfy stationarity constraints.

## Usage

``` r
check_param_stationarity(diagnostics, state = NULL)
```

## Arguments

- diagnostics:

  An object of class `ms_diagnostics`

- state:

  Integer state index (if NULL, checks all states)

## Value

List with `passed` (logical) and detailed violation information

## Examples

``` r
if (FALSE) { # \dontrun{
stationarity_check <- check_param_stationarity(diag, state = 1)
if (!stationarity_check$passed) {
  print(stationarity_check$violations)
}
} # }
```
