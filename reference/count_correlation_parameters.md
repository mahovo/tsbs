# Count Correlation Dynamics Parameters in Model

Counts the number of correlation dynamics parameters (alpha, beta,
gamma) across all states. Used for AIC/BIC comparison when refitting
with constant correlation.

## Usage

``` r
count_correlation_parameters(model_fits, model_type = "dcc_modelspec")
```

## Arguments

- model_fits:

  List of state-specific model fits

- model_type:

  Multivariate model type ("dcc_modelspec", "cgarch_modelspec", or
  "gogarch_modelspec")

## Value

Integer count of correlation parameters
