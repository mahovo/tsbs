# Simulate GOGARCH Data

Simulate data from a Generalized Orthogonal GARCH process.

## Usage

``` r
simulate_gogarch(
  n,
  k = 3,
  A = NULL,
  omega = NULL,
  alpha_garch = NULL,
  beta_garch = NULL,
  distribution = c("norm", "std"),
  shape = 8,
  seed = NULL
)
```

## Arguments

- n:

  Integer number of observations

- k:

  Integer number of series/components (default: 3)

- A:

  Matrix: mixing matrix (k x k). If NULL, a random orthogonal matrix is
  used.

- omega:

  Numeric vector of component GARCH omega parameters

- alpha_garch:

  Numeric vector of component GARCH alpha parameters

- beta_garch:

  Numeric vector of component GARCH beta parameters

- distribution:

  Character: component distribution ("norm" or "std")

- shape:

  Numeric: degrees of freedom for "std" distribution (default: 8)

- seed:

  Integer random seed

## Value

List with:

- y:

  Matrix of simulated returns (n x k)

- S:

  Matrix of independent components (n x k)

- A:

  Mixing matrix used

- W:

  Unmixing matrix

- true_params:

  List of true parameters

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate 500 observations from 3-series GOGARCH
result <- simulate_gogarch(n = 500, k = 3, seed = 42)
y <- result$y
} # }
```
