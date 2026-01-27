# Gradient of Copula NLL w.r.t. Shape (MVT)

Computes gradient of MVT copula NLL w.r.t. degrees of freedom. Uses
numerical differentiation for simplicity.

## Usage

``` r
copula_gradient_shape(shape, z, weights, alpha, beta, Qbar)
```

## Arguments

- shape:

  Current degrees of freedom

- z:

  T x k matrix of copula residuals

- weights:

  Observation weights

- alpha:

  DCC alpha

- beta:

  DCC beta

- Qbar:

  Unconditional covariance

## Value

Scalar gradient
