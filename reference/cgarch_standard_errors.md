# Compute Standard Errors for CGARCH Models

Computes standard errors for Copula GARCH DCC parameters using numerical
differentiation of the copula log-likelihood. Supports both Gaussian
(MVN) and Student-t (MVT) copulas.

## Usage

``` r
cgarch_standard_errors(
  params,
  z_matrix,
  weights,
  Qbar,
  copula_dist = "mvn",
  use_reparam = FALSE,
  method = c("hessian", "sandwich")
)
```

## Arguments

- params:

  Parameter vector:

  - For MVN copula: c(alpha, beta)

  - For MVT copula: c(alpha, beta, shape)

- z_matrix:

  T x k matrix of copula residuals (transformed standardized residuals).
  These are the result of the PIT transformation followed by inverse
  normal/t quantile transform.

- weights:

  T-vector of observation weights

- Qbar:

  k x k unconditional covariance matrix of z_matrix

- copula_dist:

  Copula distribution: "mvn" or "mvt"

- use_reparam:

  Logical: parameters in (psi, phi) space?

- method:

  SE method: "hessian" (default) or "sandwich"

## Value

List with:

- se:

  Standard errors

- vcov:

  Variance-covariance matrix

- params:

  Parameter values

- param_names:

  Parameter names

- method:

  Method used

- hessian:

  Hessian matrix (if method = "hessian")

- eigenvalues:

  Hessian eigenvalues (if method = "hessian")

## Details

The copula log-likelihood differs from DCC in that it subtracts the
marginal log-densities. For MVN copula: \$\$\ell_t = -0.5\[\log\|R_t\| +
z_t'R_t^{-1}z_t - z_t'z_t\]\$\$

For MVT copula: \$\$\ell_t = c(\nu) - 0.5\log\|R_t\| -
\frac{\nu+k}{2}\log(1 + q_t/(\nu-2)) - \sum_j \log f_t(z\_{jt})\$\$
where \\f_t\\ is the marginal Student-t density.

## See also

[`estimate_garch_weighted_cgarch`](https://mahovo.github.io/tsbs/reference/estimate_garch_weighted_cgarch.md),
[`copula_nll`](https://mahovo.github.io/tsbs/reference/copula_nll.md)
