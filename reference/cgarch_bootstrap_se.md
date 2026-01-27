# Bootstrap Standard Errors for CGARCH Parameters

Compute bootstrap standard errors for Copula GARCH correlation
parameters using residual resampling.

## Usage

``` r
cgarch_bootstrap_se(
  z_matrix,
  weights,
  Qbar,
  mle_params,
  n_boot = 200,
  method = c("residual", "parametric"),
  copula_dist = "mvn",
  verbose = TRUE,
  seed = NULL
)
```

## Arguments

- z_matrix:

  T x k matrix of copula residuals (PIT-transformed)

- weights:

  T-vector of observation weights

- Qbar:

  k x k unconditional correlation matrix

- mle_params:

  MLE estimates c(alpha, beta) or c(alpha, beta, shape)

- n_boot:

  Number of bootstrap replications (default 200)

- method:

  Bootstrap method: "residual" or "parametric"

- copula_dist:

  Copula distribution: "mvn" or "mvt"

- verbose:

  Print progress

- seed:

  Random seed for reproducibility

## Value

List with bootstrap results
