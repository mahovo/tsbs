# Robust GOGARCH Standard Errors with Edge Case Handling

High-level wrapper for GOGARCH SE computation with validation and error
handling.

## Usage

``` r
compute_gogarch_standard_errors_robust(
  gogarch_result,
  residuals,
  weights,
  distribution = "norm",
  method = c("hessian", "sandwich")
)
```

## Arguments

- gogarch_result:

  Result from estimate_garch_weighted_gogarch()

- residuals:

  Original residuals (T x k)

- weights:

  Observation weights

- distribution:

  Component distribution

- method:

  "hessian" or "sandwich"

## Value

List with SE results and validity information
