# Robust CGARCH Standard Errors with Edge Case Handling

High-level wrapper for CGARCH SE computation with handling for boundary
estimates and degenerate cases.

## Usage

``` r
compute_cgarch_standard_errors_robust(
  cgarch_result,
  z_matrix,
  weights,
  Qbar,
  copula_dist = "mvn",
  boundary_threshold = 1e-04,
  method = c("hessian", "sandwich")
)
```

## Arguments

- cgarch_result:

  Result from estimate_garch_weighted_cgarch() or similar

- z_matrix:

  Copula residuals (T x k)

- weights:

  Observation weights

- Qbar:

  Unconditional covariance

- copula_dist:

  "mvn" or "mvt"

- boundary_threshold:

  Threshold for boundary detection

- method:

  "hessian" or "sandwich"

## Value

List with SE results and validity flags
