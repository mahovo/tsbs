# Simulate DCC-GARCH Data with Higher-Order Support

Simulate realistic multivariate time series data from a DCC(q,p)-GARCH
process. Supports arbitrary DCC orders through vector-valued alpha_dcc
and beta_dcc. This function is primarily intended for testing and
examples.

## Usage

``` r
simulate_dcc_garch(
  n,
  k = 2,
  omega = c(0.05, 0.08),
  alpha_garch = c(0.1, 0.12),
  beta_garch = c(0.85, 0.82),
  alpha_dcc = 0.04,
  beta_dcc = 0.93,
  Qbar = NULL,
  seed = NULL
)
```

## Arguments

- n:

  Integer number of observations to simulate

- k:

  Integer number of series (default: 2)

- omega:

  Numeric vector of length k: GARCH intercepts (default: c(0.05, 0.08))

- alpha_garch:

  Numeric vector of length k: ARCH effects (default: c(0.10, 0.12))

- beta_garch:

  Numeric vector of length k: GARCH effects (default: c(0.85, 0.82))

- alpha_dcc:

  Numeric scalar or vector: DCC alpha parameters. Length determines q
  (number of alpha lags). (default: 0.04)

- beta_dcc:

  Numeric scalar or vector: DCC beta parameters. Length determines p
  (number of beta lags). (default: 0.93)

- Qbar:

  Matrix (k x k): Unconditional correlation matrix. If NULL, uses
  moderate correlation (0.5 off-diagonal)

- seed:

  Integer random seed for reproducibility

## Value

A matrix of dimension (n x k) with simulated returns

## Details

The function simulates data from the DCC(q,p)-GARCH model: \$\$y\_{i,t}
= \sqrt{h\_{i,t}} z\_{i,t}\$\$ \$\$h\_{i,t} = \omega_i + \alpha_i
y\_{i,t-1}^2 + \beta_i h\_{i,t-1}\$\$ \$\$Q_t =
\bar{Q}(1-\sum\alpha-\sum\beta) + \sum\_{j=1}^{q} \alpha_j
z\_{t-j}z\_{t-j}' + \sum\_{j=1}^{p} \beta_j Q\_{t-j}\$\$ \$\$R_t =
\text{diag}(Q_t)^{-1/2} Q_t \text{diag}(Q_t)^{-1/2}\$\$

where \\z_t \sim N(0, R_t)\\.

For backward compatibility, scalar alpha_dcc and beta_dcc produce
DCC(1,1). To simulate DCC(q,p), provide alpha_dcc as vector of length q
and beta_dcc as vector of length p.

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate 500 observations with default DCC(1,1) parameters
y <- simulate_dcc_garch(n = 500, seed = 42)

# DCC(1,2): one alpha lag, two beta lags
y_dcc12 <- simulate_dcc_garch(
  n = 300,
  alpha_dcc = 0.05,
  beta_dcc = c(0.50, 0.40),
  seed = 123
)

# DCC(2,2): two lags each
y_dcc22 <- simulate_dcc_garch(
  n = 300,
  alpha_dcc = c(0.03, 0.02),
  beta_dcc = c(0.50, 0.40),
  seed = 456
)

# Constant correlation (set alpha_dcc = 0)
y_const <- simulate_dcc_garch(n = 200, alpha_dcc = 0, beta_dcc = 0, seed = 789)
} # }
```
