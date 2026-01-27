# Detect States at Correlation Parameter Boundary

Identifies states where correlation dynamics parameters (alpha in
DCC/CGARCH) are at or near the lower boundary, indicating potential
degeneracy.

## Usage

``` r
detect_boundary_states(
  fit_result,
  threshold,
  criterion,
  model_type = "dcc_modelspec"
)
```

## Arguments

- fit_result:

  Result from fit_ms_varma_garch_cpp

- threshold:

  Numeric threshold for boundary detection

- criterion:

  Selection criterion ("threshold", "aic", or "bic")

- model_type:

  Multivariate model type ("dcc_modelspec", "cgarch_modelspec", or
  "gogarch_modelspec")

## Value

Integer vector of state indices that hit the boundary
