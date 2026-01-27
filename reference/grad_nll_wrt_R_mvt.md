# Gradient of NLL w.r.t. R_t (MVT)

Compute gradient of negative log-likelihood contribution w.r.t. R_t.

## Usage

``` r
grad_nll_wrt_R_mvt(z_t, R_t, R_inv_t, shape)
```

## Arguments

- z_t:

  k-vector of standardized residuals at time t

- R_t:

  k x k correlation matrix

- R_inv_t:

  k x k inverse correlation matrix

- shape:

  Degrees of freedom

## Value

k x k gradient matrix
