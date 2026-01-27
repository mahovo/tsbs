# Compute SE for Single GOGARCH Component

Computes SE for univariate GARCH parameters on one ICA component.

## Usage

``` r
gogarch_component_se(
  pars,
  component,
  weights,
  distribution = "norm",
  method = "hessian"
)
```

## Arguments

- pars:

  GARCH parameters (omega, alpha1, beta1, possibly shape, skew)

- component:

  Vector of ICA component values

- weights:

  Observation weights

- distribution:

  Component distribution

- method:

  "hessian" or "sandwich"

## Value

List with se, vcov, valid flag
