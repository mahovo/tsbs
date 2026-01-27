# Extract GOGARCH Parameter Trajectory

Extract the evolution of GOGARCH-specific parameters across EM
iterations. Includes ICA decomposition info and component GARCH
parameters.

## Usage

``` r
extract_gogarch_trajectory(diagnostics, state = 1, component = NULL)
```

## Arguments

- diagnostics:

  An object of class `ms_diagnostics`

- state:

  Integer state index (default 1)

- component:

  Integer ICA component index (if NULL, returns summary)

## Value

Data frame with iteration and component parameter evolution. Returns
NULL if not a GOGARCH model.

## Examples

``` r
if (FALSE) { # \dontrun{
traj <- extract_gogarch_trajectory(diag, state = 1)
plot(traj$iteration, traj$component_1_persistence, type = "b")
} # }
```
