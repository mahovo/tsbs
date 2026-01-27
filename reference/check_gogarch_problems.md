# Check GOGARCH-Specific Problems

Identify GOGARCH-specific estimation issues including ICA decomposition
problems, component GARCH instability, and mixing matrix
ill-conditioning.

## Usage

``` r
check_gogarch_problems(diagnostics, state = NULL)
```

## Arguments

- diagnostics:

  An object of class `ms_diagnostics`

- state:

  Integer state index (if NULL, checks all states)

## Value

List with `has_problems` (logical) and detailed problem information
