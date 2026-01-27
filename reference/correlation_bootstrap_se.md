# Bootstrap Standard Errors for Correlation Parameters

Unified bootstrap SE computation for DCC, CGARCH, and GOGARCH
correlation models.

## Usage

``` r
correlation_bootstrap_se(
  model_type = c("dcc", "cgarch", "gogarch"),
  residuals,
  weights,
  mle_result,
  n_boot = 200,
  method = c("residual", "parametric"),
  distribution = "mvn",
  verbose = TRUE,
  seed = NULL,
  ...
)
```

## Arguments

- model_type:

  Type of correlation model: "dcc", "cgarch", or "gogarch"

- residuals:

  T x k matrix of standardized residuals (or copula residuals)

- weights:

  T-vector of observation weights

- mle_result:

  List containing MLE estimates (model-specific structure)

- n_boot:

  Number of bootstrap replications (default 200)

- method:

  Bootstrap method: "residual" or "parametric"

- distribution:

  Distribution: "mvn", "mvt", "norm", "std"

- verbose:

  Print progress

- seed:

  Random seed for reproducibility

- ...:

  Additional model-specific arguments

## Value

List with bootstrap results
