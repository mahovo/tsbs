# Copula NLL with Fixed DCC Parameters

Computes copula NLL for given shape with fixed (alpha, beta). Helper for
shape gradient computation.

## Usage

``` r
copula_nll_fixed_dcc(shape, z_matrix, weights, alpha, beta, Qbar)
```

## Arguments

- shape:

  Degrees of freedom

- z_matrix:

  Copula residuals

- weights:

  Observation weights

- alpha:

  DCC alpha

- beta:

  DCC beta

- Qbar:

  Unconditional covariance

## Value

Scalar NLL
