# Compute Numerical Hessian via Finite Differences

Compute Numerical Hessian via Finite Differences

## Usage

``` r
numerical_hessian(fn, params, eps = 1e-05, ...)
```

## Arguments

- fn:

  Objective function (returns scalar)

- params:

  Parameter vector at which to evaluate Hessian

- eps:

  Step size for finite differences (default 1e-5)

- ...:

  Additional arguments passed to fn

## Value

Square matrix of second derivatives
