# Compute ADCC Standard Errors with Robust Edge Case Handling

High-level wrapper for ADCC SE computation with validation and edge case
handling (boundary estimates, near-singular Hessians).

## Usage

``` r
compute_adcc_standard_errors_robust(
  adcc_result,
  z_matrix,
  weights,
  Qbar,
  copula_dist = "mvn",
  boundary_threshold = 1e-04,
  method = c("hessian", "sandwich")
)
```

## Arguments

- adcc_result:

  Result from estimate_adcc_copula() or similar

- z_matrix:

  T x k matrix of copula residuals

- weights:

  T-vector of observation weights

- Qbar:

  k x k unconditional correlation matrix

- copula_dist:

  "mvn" or "mvt"

- boundary_threshold:

  Threshold for boundary detection (default 1e-4)

- method:

  "hessian" or "sandwich"

## Value

List with SE results and validity information
