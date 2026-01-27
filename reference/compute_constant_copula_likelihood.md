# Compute Constant Copula Correlation Likelihood

Computes weighted log-likelihood for constant correlation copula

## Usage

``` r
compute_constant_copula_likelihood(
  residuals,
  weights,
  garch_pars,
  dist_pars,
  spec,
  transformation = "parametric",
  copula_dist = "mvn"
)
```

## Arguments

- residuals:

  Matrix of residuals

- weights:

  Observation weights

- garch_pars:

  List of GARCH parameters

- dist_pars:

  Distribution parameters

- spec:

  Model specification

- transformation:

  PIT transformation type

- copula_dist:

  Copula distribution

## Value

List with weighted_ll, dist_pars, warnings
