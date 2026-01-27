# Compute Numerical Hessian with Richardson Extrapolation

More accurate Hessian using Richardson extrapolation. Combines estimates
at different step sizes for higher-order accuracy.

## Usage

``` r
numerical_hessian_richardson(fn, params, eps = 1e-04, r = 2, ...)
```

## Arguments

- fn:

  Objective function

- params:

  Parameter vector

- eps:

  Base step size (default 1e-4)

- r:

  Reduction factor (default 2)

- ...:

  Additional arguments to fn

## Value

Hessian matrix with improved accuracy
