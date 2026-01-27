# Run Monte Carlo Study for CGARCH(1,1) Estimation with ADCC Support

Performs a Monte Carlo simulation study to assess the accuracy of Copula
GARCH parameter estimation. Computes bias, RMSE, and coverage
probabilities for DCC/ADCC correlation dynamics parameters.

## Usage

``` r
run_cgarch_monte_carlo(
  n_sim = 100,
  n_obs = 500,
  k = 2,
  true_alpha = 0.05,
  true_beta = 0.9,
  true_gamma = NULL,
  omega = NULL,
  alpha_garch = NULL,
  beta_garch = NULL,
  copula = "mvn",
  true_shape = 8,
  shape = NULL,
  transformation = "parametric",
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

  Number of series (default 2)

- true_alpha:

  True DCC alpha parameter

- true_beta:

  True DCC beta parameter

- true_gamma:

  True ADCC gamma parameter (default NULL for standard DCC)

- omega:

  Vector of GARCH omega parameters (default: rep(0.05, k))

- alpha_garch:

  Vector of GARCH alpha parameters (default: rep(0.10, k))

- beta_garch:

  Vector of GARCH beta parameters (default: rep(0.85, k))

- copula:

  Copula distribution ("mvn" or "mvt")

- true_shape:

  True shape parameter for MVT copula (default 8)

- shape:

  Alias for true_shape (for backward compatibility)

- transformation:

  PIT transformation type ("parametric", "empirical", "spd")

- confidence_level:

  Confidence level for coverage (default 0.95)

- verbose:

  Print progress

- seed:

  Base seed for reproducibility

## Value

List with:

- estimates:

  Matrix of estimates (n_sim x n_params)

- std_errors:

  Matrix of standard errors (n_sim x n_params)

- bias:

  Bias for each parameter

- rmse:

  RMSE for each parameter

- coverage:

  Coverage probability for each parameter

- persistence:

  Vector of persistence values per replication

- summary:

  Summary data frame
