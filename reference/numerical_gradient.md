# Numerical Gradient (for Verification)

Compute numerical gradient using central differences.

## Usage

``` r
numerical_gradient(fn, params, eps = 1e-06, ...)
```

## Arguments

- fn:

  Objective function

- params:

  Parameter vector

- eps:

  Step size (default 1e-6)

- ...:

  Additional arguments to fn

## Value

Vector of numerical gradients
