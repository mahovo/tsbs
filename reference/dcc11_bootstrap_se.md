# Bootstrap Standard Errors for DCC(1,1) Parameters

Compute bootstrap standard errors using either parametric or residual
resampling.

## Usage

``` r
dcc11_bootstrap_se(
  std_resid,
  weights,
  Qbar,
  mle_params,
  n_boot = 200,
  method = c("residual", "parametric"),
  distribution = "mvn",
  shape = NULL,
  verbose = TRUE,
  seed = NULL
)
```

## Arguments

- std_resid:

  T x k matrix of standardized residuals

- weights:

  T-vector of observation weights

- Qbar:

  k x k unconditional correlation matrix

- mle_params:

  MLE estimates c(alpha, beta) or c(alpha, beta, shape)

- n_boot:

  Number of bootstrap replications (default 200)

- method:

  Bootstrap method: "parametric" or "residual"

- distribution:

  "mvn" or "mvt"

- shape:

  Shape parameter for MVT (required if distribution="mvt")

- verbose:

  Print progress

- seed:

  Random seed for reproducibility

## Value

List with bootstrap results
