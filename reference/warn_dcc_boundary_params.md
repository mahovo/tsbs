# Warn About DCC Parameters Near Boundary

Issues warnings when DCC parameters are near the lower boundary, which
may indicate poorly identified correlation dynamics. Only called when
returning dynamic correlation (not when falling back to constant).

## Usage

``` r
warn_dcc_boundary_params(
  alpha_params,
  beta_params,
  threshold,
  state,
  iteration,
  diagnostics
)
```

## Arguments

- alpha_params:

  Named list of DCC alpha parameters

- beta_params:

  Named list of DCC beta parameters

- threshold:

  Numeric threshold for boundary warning (default 1e-4)

- state:

  State index (for warning message)

- iteration:

  EM iteration (for diagnostics)

- diagnostics:

  Diagnostics collector object

## Value

Updated diagnostics object
