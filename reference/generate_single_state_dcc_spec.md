# Generate Single-State DCC Specification

Creates a specification for a single-state (no regime switching)
DCC-GARCH model. Useful for testing parameter recovery.

## Usage

``` r
generate_single_state_dcc_spec(
  k = 2,
  var_order = 1,
  garch_order = c(1, 1),
  distribution = "mvn",
  omega = NULL,
  alpha_garch = NULL,
  beta_garch = NULL,
  alpha_dcc = NULL,
  beta_dcc = NULL,
  seed = NULL
)
```
