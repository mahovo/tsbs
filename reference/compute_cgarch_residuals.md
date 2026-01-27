# Compute Copula Residuals from Simulated CGARCH Data

Given simulated returns and known GARCH parameters, compute the
standardized residuals, apply PIT transformation, and convert to copula
residuals needed for CGARCH estimation.

## Usage

``` r
compute_cgarch_residuals(
  y,
  omega,
  alpha_garch,
  beta_garch,
  transformation = "parametric",
  copula = "mvn",
  shape = 8
)
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

- transformation:

  PIT transformation type ("parametric", "empirical", "spd")

- copula:

  Copula distribution ("mvn" or "mvt")

- shape:

  Degrees of freedom for MVT copula (default 8)

## Value

List with:

- std_resid:

  T x k matrix of standardized residuals

- u_matrix:

  T x k matrix of PIT-transformed uniforms

- z_matrix:

  T x k matrix of copula residuals
