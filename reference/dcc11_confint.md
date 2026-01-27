# Compute Confidence Intervals for DCC(1,1) Parameters

Compute Confidence Intervals for DCC(1,1) Parameters

## Usage

``` r
dcc11_confint(se_result, level = 0.95)
```

## Arguments

- se_result:

  Result from dcc11_standard_errors() or dcc11_hessian()

- level:

  Confidence level (default 0.95)

## Value

Matrix with columns: estimate, se, lower, upper
