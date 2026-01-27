# Identify Problematic States

Identify states with estimation problems based on diagnostic
information. Supports DCC, CGARCH, and GOGARCH models with
model-specific checks.

## Usage

``` r
identify_problematic_states(
  diagnostics,
  state = NULL,
  model_type = c("auto", "dcc", "cgarch", "gogarch")
)
```

## Arguments

- diagnostics:

  An object of class `ms_diagnostics`

- state:

  Integer state index. If NULL (default), checks all states.

- model_type:

  Character: model type override. One of "auto" (default), "dcc",
  "cgarch", or "gogarch". When "auto", attempts to detect from
  diagnostics.

## Value

List with:

- has_problems:

  Logical indicating if any problems were found

- n_states_affected:

  Number of states with problems

- problems:

  Named list with problem descriptions per state

## Details

The function performs different checks depending on model type:

- DCC:

  - High persistence (alpha + beta \> 0.98)

  - Constant correlation fallback

  - Parameter instability in final iterations

  - Boundary events

- CGARCH:

  - All DCC checks

  - Copula shape parameter issues (MVT: df \< 3 or \> 100)

  - ADCC gamma constraints

  - PIT transformation warnings

- GOGARCH:

  - ICA convergence failure

  - Mixing matrix ill-conditioning (condition number \> 1000)

  - Unmixing matrix near-singularity

  - Component correlation (should be \< 0.2)

  - Component GARCH high persistence

## Examples

``` r
if (FALSE) { # \dontrun{
# Check all states with auto-detection
problems <- identify_problematic_states(diag)
if (problems$has_problems) {
  print(problems$problems)
}

# Check specific state for CGARCH model
problems <- identify_problematic_states(diag, state = 1, model_type = "cgarch")
} # }
```
