# Compute Standardized Residuals from Simulated DCC-GARCH Data

Given simulated returns and known GARCH parameters, compute the
standardized residuals needed for DCC estimation.

## Usage

``` r
compute_std_residuals(y, omega, alpha_garch, beta_garch)
```

## Arguments

- y:

  T x k matrix of simulated returns

- omega:

  Vector of GARCH omega parameters (length k)

- alpha_garch:

  Vector of GARCH alpha parameters (length k)

- beta_garch:

  Vector of GARCH beta parameters (length k)

## Value

T x k matrix of standardized residuals
