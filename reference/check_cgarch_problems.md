# Check CGARCH-Specific Problems

Identify CGARCH-specific estimation issues including copula parameter
problems, PIT transformation warnings, and ADCC dynamics issues.

## Usage

``` r
check_cgarch_problems(diagnostics, state = NULL)
```

## Arguments

- diagnostics:

  An object of class `ms_diagnostics`

- state:

  Integer state index (if NULL, checks all states)

## Value

List with `has_problems` (logical) and detailed problem information
