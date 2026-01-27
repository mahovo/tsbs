# Extract Parameter Trajectory

Extract the evolution of parameters across EM iterations for a given
state. Supports DCC, CGARCH, and GOGARCH model types with automatic
detection.

## Usage

``` r
extract_param_trajectory(
  diagnostics,
  state,
  param_name,
  series = NULL,
  component = NULL,
  model_type = c("auto", "dcc", "cgarch", "gogarch")
)
```

## Arguments

- diagnostics:

  An object of class `ms_diagnostics`

- state:

  Integer state index

- param_name:

  Character string specifying parameter name. For single parameter
  extraction use names like "alpha_1", "beta_1", "gamma_1", "shape",
  "omega". Use "all" to extract all relevant correlation parameters
  based on model type.

- series:

  Optional integer for series-specific GARCH parameters (omega, alpha1,
  beta1)

- component:

  Optional integer for GOGARCH component-specific parameters

- model_type:

  Character: model type override. One of "auto" (default), "dcc",
  "cgarch", or "gogarch". When "auto", attempts to detect from
  diagnostics.

## Value

For single parameter: a data.frame with columns `iteration` and `value`.
For "all": a data.frame with iteration and all relevant parameters as
columns. Returns NULL if parameter not found.

## Details

The function extracts different parameters depending on model type:

- **DCC**: alpha_1, beta_1 (correlation dynamics), shape (if MVT)

- **CGARCH**: alpha_1, beta_1, gamma_1 (if ADCC), shape (if MVT copula)

- **GOGARCH**: component GARCH parameters (omega, alpha1, beta1 per
  component)

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract single DCC parameter
alpha_traj <- extract_param_trajectory(diag, state = 1, param_name = "alpha_1")
plot(alpha_traj$iteration, alpha_traj$value, type = "b")

# Extract all correlation parameters (auto-detects model type)
all_traj <- extract_param_trajectory(diag, state = 1, param_name = "all")

# Extract GOGARCH component parameter
comp1_alpha <- extract_param_trajectory(diag, state = 1, param_name = "alpha1", 
                                         component = 1, model_type = "gogarch")

# Extract series-specific GARCH omega
omega_traj <- extract_param_trajectory(diag, state = 2, param_name = "omega", series = 1)
} # }
```
