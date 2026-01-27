# Generate CGARCH Specification

Generate a properly structured specification list for MS-VARMA-GARCH
estimation with CGARCH (Copula-GARCH) dynamics. Supports DCC, ADCC
dynamics and MVT copula.

## Usage

``` r
generate_cgarch_spec(
  M,
  k = 2,
  var_order = 1,
  garch_order = c(1, 1),
  dynamics = c("dcc", "adcc", "constant"),
  copula = c("mvn", "mvt"),
  transformation = c("parametric", "empirical", "spd"),
  seed = NULL,
  simple = FALSE
)
```

## Arguments

- M:

  Integer number of states

- k:

  Integer number of series (default: 2)

- var_order:

  Integer VAR order (default: 1)

- garch_order:

  Integer vector of length 2: GARCH(p,q) order (default: c(1,1))

- dynamics:

  Character: correlation dynamics ("dcc", "adcc", or "constant")

- copula:

  Character: copula type ("mvn" or "mvt")

- transformation:

  Character: PIT method ("parametric", "empirical", or "spd")

- seed:

  Integer random seed for generating starting values

- simple:

  Logical: if TRUE, use identical starting values across states

## Value

A list of length M containing properly formatted specifications

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate 2-state CGARCH specification with Student-t copula
spec <- generate_cgarch_spec(M = 2, k = 3, copula = "mvt")

# ADCC dynamics with SPD transformation
spec_adcc <- generate_cgarch_spec(M = 2, dynamics = "adcc", transformation = "spd")
} # }
```
