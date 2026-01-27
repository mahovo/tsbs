# Extract CGARCH Parameter Trajectory

Extract the evolution of CGARCH-specific parameters across EM
iterations. Includes DCC/ADCC dynamics parameters and copula parameters.

Extract the evolution of CGARCH parameters (alpha, beta, gamma, shape)
across EM iterations from an ms_diagnostics object.

## Usage

``` r
extract_cgarch_trajectory(diagnostics, state = 1)

extract_cgarch_trajectory(diagnostics, state = 1)
```

## Arguments

- diagnostics:

  An object of class `ms_diagnostics` from fit_ms_varma_garch()

- state:

  Integer state index (default 1)

## Value

Data frame with columns: iteration, alpha, beta, gamma (if ADCC), shape
(if MVT copula). Returns NULL if not a CGARCH model.

Data frame with columns: iteration, alpha, beta, gamma (if ADCC), shape
(if MVT) Returns NULL if CGARCH parameters not found.

## Examples

``` r
if (FALSE) { # \dontrun{
traj <- extract_cgarch_trajectory(diag, state = 1)
plot(traj$iteration, traj$alpha, type = "b")
} # }
```
