# Gradient of Copula NLL w.r.t. R_t (MVN)

Compute gradient of Gaussian copula NLL contribution w.r.t. R_t. For
Gaussian copula: LL_t = -0.5 \* (log\|R\| + z'R^-1z - z'z) Gradient:
d(-LL_t)/dR = 0.5 \* (R^-1 - R^-1zz'R^-1)

## Usage

``` r
grad_nll_wrt_R_copula_mvn(z_t, R_t, R_inv_t)
```

## Arguments

- z_t:

  k-vector of copula residuals at time t

- R_t:

  k x k correlation matrix

- R_inv_t:

  k x k inverse correlation matrix

## Value

k x k gradient matrix
