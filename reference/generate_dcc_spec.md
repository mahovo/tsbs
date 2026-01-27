# Generate Model Specification for DCC-GARCH

Generate a properly structured specification list for MS-VARMA-GARCH
estimation with DCC dynamics. This function is primarily intended for
testing and examples.

## Usage

``` r
generate_dcc_spec(
  M,
  k = 2,
  var_order = 1,
  garch_order = c(1, 1),
  distribution = c("mvn", "mvt"),
  seed = NULL,
  simple = FALSE
)
```

## Arguments

- M:

  Integer number of states

- k:

  Integer number of series (default: 2)

- var_order:

  Integer VAR order (default: 1)

- garch_order:

  Integer vector of length 2: GARCH(p,q) order (default: c(1,1))

- distribution:

  Character: distribution for DCC ("mvn" or "mvt")

- seed:

  Integer random seed for generating starting values

- simple:

  Logical: if TRUE, use identical starting values across states (faster
  but may not reflect regime differences). Default: FALSE

## Value

A list of length M containing properly formatted specifications

## Details

The function generates starting parameter values with mild regime
differentiation unless `simple = TRUE`. States are ordered from low to
high volatility.

For state j:

- VAR parameters: Small values near 0.1

- GARCH omega: Increasing across states (0.05 to 0.15)

- GARCH alpha: Increasing across states (0.08 to 0.15)

- GARCH beta: Decreasing across states (0.85 to 0.75)

- DCC alpha: Increasing across states (0.03 to 0.10)

- DCC beta: Decreasing across states (0.94 to 0.85)

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate 2-state specification with defaults
spec <- generate_dcc_spec(M = 2)

# 3-state specification with Student-t distribution
spec_mvt <- generate_dcc_spec(M = 3, distribution = "mvt")

# Simple specification (identical starting values)
spec_simple <- generate_dcc_spec(M = 2, simple = TRUE)

# Use with simulated data
y <- simulate_dcc_garch(n = 300, seed = 42)
fit <- fit_ms_varma_garch(
  y = y, 
  M = 2, 
  spec = spec,
  model_type = "multivariate",
  control = list(max_iter = 20, tol = 1e-4)
)
} # }
```
