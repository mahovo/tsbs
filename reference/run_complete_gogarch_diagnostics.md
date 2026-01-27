# Complete GOGARCH Diagnostics

Run comprehensive diagnostics for GOGARCH models.

## Usage

``` r
run_complete_gogarch_diagnostics(
  n_sim = 50,
  n_obs = 500,
  k = 3,
  omega = NULL,
  alpha_garch = NULL,
  beta_garch = NULL,
  distribution = "norm",
  ica_method = "radical",
  seed = 42,
  plot = TRUE
)
```

## Arguments

- n_sim:

  Number of MC replications

- n_obs:

  Observations per replication

- k:

  Number of series/components

- omega:

  Vector of component GARCH omega

- alpha_garch:

  Vector of component GARCH alpha

- beta_garch:

  Vector of component GARCH beta

- distribution:

  Component distribution

- ica_method:

  ICA algorithm

- seed:

  Random seed

- plot:

  Logical: create plots?

## Value

List with diagnostic results
