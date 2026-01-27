# Comprehensive Correlation Model Inference

Compute all three types of inference (Hessian, Bootstrap, Profile) for
DCC or CGARCH correlation parameters.

## Usage

``` r
correlation_comprehensive_inference(
  model_type = c("dcc", "cgarch", "adcc"),
  residuals,
  weights,
  Qbar = NULL,
  mle_params = NULL,
  distribution = "mvn",
  n_boot = 200,
  boot_method = "residual",
  n_profile_points = 50,
  conf_level = 0.95,
  verbose = TRUE,
  seed = NULL
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

  k x k unconditional correlation matrix (optional)

- mle_params:

  MLE estimates (optional, computed if NULL)

- distribution:

  "mvn" or "mvt"

- n_boot:

  Number of bootstrap replications

- boot_method:

  Bootstrap method

- n_profile_points:

  Points for profile likelihood

- conf_level:

  Confidence level

- verbose:

  Print progress

- seed:

  Random seed

## Value

List with all inference results
