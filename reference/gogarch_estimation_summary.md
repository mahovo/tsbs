# GOGARCH Estimation Summary

Provides a summary of GOGARCH estimation including parameter estimates
and standard errors for all components.

## Usage

``` r
gogarch_estimation_summary(
  gogarch_result,
  residuals,
  weights,
  distribution = "norm",
  level = 0.95,
  method = c("hessian", "sandwich")
)
```

## Arguments

- gogarch_result:

  Result from estimate_garch_weighted_gogarch()

- residuals:

  Original residuals

- weights:

  Observation weights

- distribution:

  Component distribution

- level:

  Confidence level

- method:

  SE method

## Value

Object of class "gogarch_summary"
