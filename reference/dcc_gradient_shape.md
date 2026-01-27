# Gradient of NLL w.r.t. Shape Parameter (MVT)

Compute gradient of weighted NLL w.r.t. degrees of freedom.

The MVT log-likelihood (per observation) is: ll_t = lgamma((nu+k)/2) -
lgamma(nu/2) - (k/2)*log(pi*(nu-2)) - 0.5\*log\|R_t\| -
((nu+k)/2)\*log(kappa_t)

where kappa_t = 1 + q_t/(nu-2) and q_t = z_t' R_inv z_t.

## Usage

``` r
dcc_gradient_shape(shape, z, weights, R_inv)
```

## Arguments

- shape:

  Degrees of freedom (nu)

- z:

  T x k matrix of standardized residuals

- weights:

  T-vector of observation weights

- R_inv:

  Array of inverse correlation matrices

## Value

Scalar gradient (unnamed)
