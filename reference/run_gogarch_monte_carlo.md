# Run Monte Carlo Study for GOGARCH Estimation

Performs a Monte Carlo simulation study to assess the accuracy of
GOGARCH parameter estimation. Since GOGARCH estimates component-wise
GARCH parameters, this evaluates estimation accuracy for each component.

## Usage

``` r
run_gogarch_monte_carlo(
  n_sim = 100,
  n_obs = 500,
  k = 3,
  true_omega = NULL,
  true_alpha = NULL,
  true_beta = NULL,
  omega = NULL,
  alpha_garch = NULL,
  beta_garch = NULL,
  distribution = "norm",
  shape = 8,
  ica_method = "radical",
  confidence_level = 0.95,
  verbose = TRUE,
  seed = 12345
)
```

## Arguments

- n_sim:

  Number of simulation replications

- n_obs:

  Number of observations per replication

- k:

  Number of series/components (default 3)

- true_omega:

  Vector of true component GARCH omega parameters (alias for omega)

- true_alpha:

  Vector of true component GARCH alpha parameters (alias for
  alpha_garch)

- true_beta:

  Vector of true component GARCH beta parameters (alias for beta_garch)

- omega:

  Vector of component GARCH omega parameters (alternative to true_omega)

- alpha_garch:

  Vector of component GARCH alpha parameters (alternative to true_alpha)

- beta_garch:

  Vector of component GARCH beta parameters (alternative to true_beta)

- distribution:

  Component distribution ("norm" or "std")

- shape:

  Degrees of freedom for "std" distribution

- ica_method:

  ICA algorithm ("radical" or "fastica")

- confidence_level:

  Confidence level for coverage (default 0.95)

- verbose:

  Print progress

- seed:

  Base seed for reproducibility

## Value

List with:

- estimates:

  List of data frames (one per component) with alpha, beta columns

- persistence:

  List of vectors (one per component)

- convergence:

  Logical vector indicating optimization convergence

- ica_converged:

  Logical vector indicating ICA convergence

- bias:

  List of named vectors (one per component)

- rmse:

  List of named vectors (one per component)

- empirical_sd:

  List of named vectors (one per component)

- coverage:

  List of named vectors (one per component)

- mixing_recovery:

  Statistics on mixing matrix recovery

- summary:

  Summary data frame
