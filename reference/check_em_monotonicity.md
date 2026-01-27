# Check EM Monotonicity

Verify that the EM algorithm exhibits (near) monotonic log-likelihood
improvement.

## Usage

``` r
check_em_monotonicity(diagnostics, tolerance = 1e-04)
```

## Arguments

- diagnostics:

  An object of class `ms_diagnostics`

- tolerance:

  Numeric threshold for acceptable decreases (default: 1e-4)

## Value

List with `passed` (logical), `n_violations`, `violation_iters`

## Examples

``` r
if (FALSE) { # \dontrun{
monotonicity_check <- check_em_monotonicity(diag)
if (!monotonicity_check$passed) {
  cat("Violations at iterations:", monotonicity_check$violation_iters, "\n")
}
} # }
```
