# Jacobian of Reparameterization

Compute d(alpha, beta) / d(psi, phi) for chain rule.

## Usage

``` r
dcc11_reparam_jacobian(psi, phi)
```

## Arguments

- psi:

  Unconstrained persistence parameter

- phi:

  Unconstrained ratio parameter

## Value

2x2 matrix: rows = (alpha, beta), cols = (psi, phi)
