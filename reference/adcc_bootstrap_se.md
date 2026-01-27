# Bootstrap Standard Errors for ADCC Parameters

Compute bootstrap standard errors for ADCC correlation parameters
(alpha, gamma, beta) using residual or parametric resampling.

## Usage

``` r
adcc_bootstrap_se(
  z_matrix,
  weights,
  Qbar,
  mle_params,
  Nbar = NULL,
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

  MLE estimates c(alpha, gamma, beta) or with shape for MVT

- Nbar:

  k x k average outer product of negative residuals (optional)

- n_boot:

  Number of bootstrap replications (default 200)

- method:

  Bootstrap method: "residual" or "parametric"

- copula_dist:

  "mvn" or "mvt"

- verbose:

  Print progress

- seed:

  Random seed for reproducibility

## Value

List with bootstrap results including SEs for alpha, gamma, beta
