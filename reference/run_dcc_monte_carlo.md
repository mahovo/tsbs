# Run Monte Carlo Study for DCC(1,1) Estimation

Performs a Monte Carlo simulation study to assess the accuracy of DCC
parameter estimation. Computes bias, RMSE, and coverage probabilities.

## Usage

``` r
run_dcc_monte_carlo(
  n_sim = 100,
  n_obs = 500,
  k = 2,
  true_alpha = 0.05,
  true_beta = 0.9,
  omega = NULL,
  alpha_garch = NULL,
  beta_garch = NULL,
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

- omega:

  Vector of GARCH omega parameters (default: rep(0.05, k))

- alpha_garch:

  Vector of GARCH alpha parameters (default: rep(0.10, k))

- beta_garch:

  Vector of GARCH beta parameters (default: rep(0.85, k))

- confidence_level:

  Confidence level for coverage (default 0.95)

- verbose:

  Print progress

- seed:

  Base seed for reproducibility

## Value

List with:

- estimates:

  Matrix of estimates (n_sim x 2)

- std_errors:

  Matrix of standard errors (n_sim x 2)

- bias:

  Bias for each parameter

- rmse:

  RMSE for each parameter

- coverage:

  Coverage probability for each parameter

- summary:

  Summary data frame
