# Check Convergence Achievement

Verify that the model achieved convergence based on tolerance criteria.

## Usage

``` r
check_convergence(diagnostics, tolerance = 1e-04)
```

## Arguments

- diagnostics:

  An object of class `ms_diagnostics`

- tolerance:

  Numeric convergence tolerance

## Value

List with `converged` (logical), `final_change`, `n_iterations`
