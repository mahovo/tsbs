# Generate GOGARCH Specification

Generate a properly structured specification list for MS-VARMA-GARCH
estimation with GOGARCH (Generalized Orthogonal GARCH) dynamics.

## Usage

``` r
generate_gogarch_spec(
  M,
  k = 3,
  var_order = 1,
  garch_order = c(1, 1),
  ica_method = c("radical", "fastica"),
  n_components = k,
  distribution = c("norm", "nig", "gh"),
  seed = NULL,
  simple = FALSE
)
```

## Arguments

- M:

  Integer number of states

- k:

  Integer number of series (default: 3)

- var_order:

  Integer VAR order (default: 1)

- garch_order:

  Integer vector of length 2: GARCH(p,q) order (default: c(1,1))

- ica_method:

  Character: ICA algorithm ("radical" or "fastica")

- n_components:

  Integer: number of ICA components (default: k)

- distribution:

  Character: component distribution ("norm", "nig", or "gh")

- seed:

  Integer random seed

- simple:

  Logical: if TRUE, use identical starting values across states

## Value

A list of length M containing properly formatted specifications

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate 2-state GOGARCH specification
spec <- generate_gogarch_spec(M = 2, k = 3)

# With FastICA and NIG distribution
spec_nig <- generate_gogarch_spec(M = 2, ica_method = "fastica", distribution = "nig")
} # }
```
