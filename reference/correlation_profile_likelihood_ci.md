# Profile Likelihood Confidence Intervals

Compute profile likelihood CIs for correlation parameters. Supports DCC
and CGARCH models.

## Usage

``` r
correlation_profile_likelihood_ci(
  model_type = c("dcc", "cgarch"),
  residuals,
  weights,
  Qbar,
  mle_params,
  mle_nll = NULL,
  conf_level = 0.95,
  param = c("both", "alpha", "beta"),
  n_points = 50,
  distribution = "mvn",
  verbose = TRUE
)
```

## Arguments

- model_type:

  "dcc" or "cgarch"

- residuals:

  T x k matrix of residuals

- weights:

  T-vector of observation weights

- Qbar:

  k x k unconditional correlation matrix

- mle_params:

  MLE estimates

- mle_nll:

  NLL at MLE

- conf_level:

  Confidence level (default 0.95)

- param:

  Which parameter: "alpha", "beta", or "both"

- n_points:

  Number of profile points

- distribution:

  Distribution type

- verbose:

  Print progress

## Value

List with profile likelihood results
