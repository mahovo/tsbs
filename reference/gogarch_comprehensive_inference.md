# GOGARCH Comprehensive Inference

Compute bootstrap-based inference for GOGARCH model. Profile likelihood
is not available for GOGARCH due to the ICA structure.

## Usage

``` r
gogarch_comprehensive_inference(
  residuals,
  weights,
  garch_pars,
  ica_info,
  n_boot = 200,
  boot_method = "residual",
  distribution = "norm",
  conf_level = 0.95,
  verbose = TRUE,
  seed = NULL
)
```

## Arguments

- residuals:

  T x k matrix of observed residuals

- weights:

  T-vector of observation weights

- garch_pars:

  List of GARCH parameters per component

- ica_info:

  ICA decomposition information

- n_boot:

  Number of bootstrap replications

- boot_method:

  Bootstrap method

- distribution:

  "norm" or "std"

- conf_level:

  Confidence level

- verbose:

  Print progress

- seed:

  Random seed

## Value

List with inference results
