# Estimate Copula Parameters with Weighted Likelihood

Estimates copula correlation parameters using weighted maximum
likelihood. This function handles the copula-specific transformation
from standardized residuals to uniform margins before applying the
DCC-style correlation dynamics.

## Usage

``` r
estimate_copula_parameters_weighted(
  residuals,
  weights,
  garch_pars,
  dcc_start_pars,
  dist_start_pars,
  spec,
  transformation = "parametric",
  copula_dist = "mvn",
  dynamics = "dcc",
  diagnostics = NULL,
  iteration = NULL,
  state = NULL,
  verbose = FALSE
)
```

## Arguments

- residuals:

  T x k matrix of residuals

- weights:

  T-vector of observation weights

- garch_pars:

  List of GARCH parameters per series

- dcc_start_pars:

  Named list of starting DCC parameters

- dist_start_pars:

  Named list of distribution parameters (e.g., shape)

- spec:

  Model specification list

- transformation:

  Type of PIT transformation ("parametric", "empirical", "spd")

- copula_dist:

  Copula distribution ("mvn" or "mvt")

- diagnostics:

  Diagnostics object (optional)

- iteration:

  Current EM iteration (optional)

- state:

  Current regime state (optional)

- verbose:

  Logical; print diagnostic information

## Value

List with dcc_pars, dist_pars, weighted_ll, warnings, diagnostics
