# Estimate GARCH Parameters for Multivariate Models (R Helper)

Two-stage weighted estimation for multivariate GARCH models: Stage 1:
Univariate GARCH parameters for each series Stage 2: Multivariate
dependence parameters (DCC/Copula) or decomposition (GOGARCH)

## Usage

``` r
estimate_garch_weighted_r(
  residuals,
  weights,
  spec,
  model_type = "univariate",
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

  Matrix of residuals (T x k) for multivariate case

- weights:

  Vector of weights from E-step

- spec:

  Model specification

- model_type:

  Either "univariate" or "multivariate"

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

## Value

List with coefficients and warnings
