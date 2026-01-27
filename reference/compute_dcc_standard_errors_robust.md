# Compute DCC Standard Errors with Robust Edge Case Handling

High-level wrapper with handling for boundary estimates, constant
correlation, and non-positive-definite Hessians.

## Usage

``` r
compute_dcc_standard_errors_robust(
  dcc_result,
  std_resid,
  weights,
  Qbar,
  distribution = "mvn",
  boundary_threshold = 1e-04,
  method = c("hessian", "sandwich")
)
```

## Arguments

- dcc_result:

  Result from estimate_dcc_parameters_weighted() or similar

- std_resid:

  T x k matrix of standardized residuals

- weights:

  T-vector of observation weights

- Qbar:

  k x k unconditional correlation matrix

- distribution:

  "mvn" or "mvt"

- boundary_threshold:

  Threshold for boundary detection (default 1e-4)

- method:

  "hessian" or "sandwich"

## Value

List with:

- se:

  se

- vcov:

  vcov

- valid:

  valid

- reason:

  reason

- correlation_type:

  correlation_type

- info:

  info

- info_eigenvalues:

  info_eigenvalues
