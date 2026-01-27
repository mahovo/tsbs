# Multivariate GARCH Parameter Estimation (Two-Stage)

Implements weighted MLE for multivariate GARCH models Stage 1: Estimate
univariate GARCH parameters for each series Stage 2: Estimate dependence
parameters (DCC/Copula) or GOGARCH rotation

## Usage

``` r
estimate_garch_weighted_multivariate(
  residuals,
  weights,
  spec,
  diagnostics = NULL,
  iteration = NULL,
  state = NULL,
  verbose = FALSE,
  dcc_threshold = 0.02,
  dcc_criterion = "bic",
  force_constant = FALSE
)
```

## Arguments

- residuals:

  residuals

- weights:

  weights

- spec:

  spec

- diagnostics:

  diagnostics

- iteration:

  iteration

- state:

  state

- verbose:

  verbose

- dcc_threshold:

  dcc_threshold

- dcc_criterion:

  dcc_criterion

- force_constant:

  force_constant
