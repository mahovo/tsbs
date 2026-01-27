# Summarize DCC(1,1) Estimation Results

Summarize DCC(1,1) Estimation Results

## Usage

``` r
dcc11_estimation_summary(
  params,
  std_resid,
  weights,
  Qbar,
  distribution = "mvn",
  use_reparam = FALSE,
  level = 0.95,
  method = c("hessian", "sandwich")
)
```

## Arguments

- params:

  MLE estimates

- std_resid:

  Standardized residuals

- weights:

  Observation weights

- Qbar:

  Unconditional covariance

- distribution:

  "mvn" or "mvt"

- use_reparam:

  Logical

- level:

  Confidence level (default 0.95)

- method:

  SE method: "hessian" or "sandwich"

## Value

List of class "dcc11_summary"
