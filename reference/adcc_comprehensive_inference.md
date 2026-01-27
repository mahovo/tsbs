# Comprehensive ADCC Inference

Compute Hessian-based and bootstrap standard errors for ADCC parameters,
with comparison and diagnostics.

## Usage

``` r
adcc_comprehensive_inference(
  z_matrix,
  weights,
  adcc_result,
  Qbar = NULL,
  copula_dist = "mvn",
  n_boot = 200,
  boot_method = "residual",
  conf_level = 0.95,
  verbose = TRUE,
  seed = NULL
)
```

## Arguments

- z_matrix:

  T x k matrix of copula residuals

- weights:

  T-vector of observation weights

- adcc_result:

  Result from estimate_adcc_copula()

- Qbar:

  k x k unconditional correlation matrix

- copula_dist:

  "mvn" or "mvt"

- n_boot:

  Number of bootstrap replications

- boot_method:

  "residual" or "parametric"

- conf_level:

  Confidence level

- verbose:

  Print progress

- seed:

  Random seed

## Value

List with comprehensive inference results
