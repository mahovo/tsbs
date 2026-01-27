# Test Normality of Standardized Estimates (Generic)

Test whether (theta_hat - theta_0) / SE follows N(0,1). Works for both
CGARCH and DCC results.

## Usage

``` r
test_estimate_normality_generic(mc_result, param_names = c("alpha", "beta"))
```

## Arguments

- mc_result:

  Result from run_cgarch_monte_carlo() or run_dcc_monte_carlo()

- param_names:

  Names of parameters to test (default: c("alpha", "beta"))

## Value

List with test statistics for each parameter
