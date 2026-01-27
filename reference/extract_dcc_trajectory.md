# Extract DCC Parameter Trajectory from MS Diagnostics

Extract the evolution of DCC parameters (alpha, beta) across EM
iterations from an ms_diagnostics object.

## Usage

``` r
extract_dcc_trajectory(diagnostics, state = 1)
```

## Arguments

- diagnostics:

  An object of class `ms_diagnostics` from fit_ms_varma_garch()

- state:

  Integer state index (default 1)

## Value

Data frame with columns: iteration, alpha, beta Returns NULL if DCC
parameters not found.
