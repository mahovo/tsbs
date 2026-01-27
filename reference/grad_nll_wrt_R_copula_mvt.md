# Gradient of Copula NLL w.r.t. R_t (MVT)

Compute gradient of Student-t copula NLL contribution w.r.t. R_t. For
Student-t copula, the gradient has an additional scaling factor.

## Usage

``` r
grad_nll_wrt_R_copula_mvt(z_t, R_t, R_inv_t, shape)
```

## Arguments

- z_t:

  k-vector of copula residuals at time t

- R_t:

  k x k correlation matrix

- R_inv_t:

  k x k inverse correlation matrix

- shape:

  Degrees of freedom

## Value

k x k gradient matrix
