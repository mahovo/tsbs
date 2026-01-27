# Simulate CGARCH Data with ADCC Support

Simulate data from a Copula-GARCH process with DCC or ADCC dynamics.

## Usage

``` r
simulate_cgarch(
  n,
  k = 2,
  omega = NULL,
  alpha_garch = NULL,
  beta_garch = NULL,
  alpha_dcc = 0.04,
  beta_dcc = 0.93,
  gamma_dcc = NULL,
  copula = "mvn",
  shape = 8,
  Qbar = NULL,
  seed = NULL
)
```

## Arguments

- n:

  Integer number of observations

- k:

  Integer number of series (default: 2)

- omega:

  Numeric vector of GARCH omega parameters (length k)

- alpha_garch:

  Numeric vector of GARCH alpha parameters (length k)

- beta_garch:

  Numeric vector of GARCH beta parameters (length k)

- alpha_dcc:

  Numeric: DCC alpha parameter (default: 0.04)

- beta_dcc:

  Numeric: DCC beta parameter (default: 0.93)

- gamma_dcc:

  Numeric: ADCC gamma parameter for leverage (default: NULL for standard
  DCC)

- copula:

  Character: copula type ("mvn" or "mvt")

- shape:

  Numeric: degrees of freedom for MVT copula (default: 8)

- Qbar:

  Matrix: unconditional correlation matrix. If NULL, uses moderate
  correlation.

- seed:

  Integer random seed

## Value

Matrix of simulated returns (n x k)

## Details

When gamma_dcc is provided (non-NULL and non-zero), the function
simulates from an ADCC (Asymmetric DCC) process where negative shocks
have a larger impact on correlation dynamics:

Q_t = Omega + alpha \* (z_t-1 z'*t-1) + gamma \* (n*t-1 n'*t-1) + beta
\* Q*t-1

where n_t = z_t \* I(z_t \< 0) captures negative shocks only.
