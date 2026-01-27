# Estimate DCC Parameters with Weighted Likelihood

Estimates DCC parameters using weighted maximum likelihood, with
specialized handling for DCC(1,1) models (analytical gradients) and
higher-order models (finite difference with fine step size).

For DCC(1,1), the optimization is performed in reparameterized space:
psi = logit(alpha + beta) – persistence in logit space phi = log(alpha /
beta) – ratio in log space

This ensures stationarity by construction and provides smooth gradients.

## Usage

``` r
estimate_dcc_parameters_weighted(
  residuals,
  weights,
  garch_pars,
  dcc_start_pars,
  dist_start_pars,
  spec,
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

  Named list of distribution parameters

- spec:

  Model specification list

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
