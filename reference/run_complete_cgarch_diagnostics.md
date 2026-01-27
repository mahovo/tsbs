# Summary Function for Complete CGARCH Diagnostic Analysis

Run a complete diagnostic analysis including MC study, likelihood
surface, and statistical tests for CGARCH models.

## Usage

``` r
run_complete_cgarch_diagnostics(
  n_sim = 50,
  n_obs = 500,
  true_alpha = 0.05,
  true_beta = 0.9,
  omega = c(0.05, 0.08),
  alpha_garch = c(0.1, 0.12),
  beta_garch = c(0.85, 0.82),
  copula = "mvn",
  true_shape = 8,
  transformation = "parametric",
  seed = 42,
  plot = TRUE
)
```

## Arguments

- n_sim:

  Number of MC replications

- n_obs:

  Observations per replication

- true_alpha:

  True DCC alpha

- true_beta:

  True DCC beta

- omega:

  Vector of GARCH omega parameters

- alpha_garch:

  Vector of GARCH alpha parameters

- beta_garch:

  Vector of GARCH beta parameters

- copula:

  Copula distribution ("mvn" or "mvt")

- true_shape:

  True shape for MVT copula

- transformation:

  PIT transformation type

- seed:

  Random seed

- plot:

  Logical: create plots?

## Value

List with all diagnostic results
