# Analytical Hessian of DCC(1,1) Negative Log-Likelihood

Computes the exact Hessian matrix of the DCC(1,1) negative
log-likelihood using analytical derivatives. This is more accurate than
numerical differentiation but only implemented for the 2-parameter
(alpha, beta) case.

## Usage

``` r
dcc11_analytical_hessian(
  params,
  std_resid,
  weights,
  Qbar,
  distribution = "mvn"
)
```

## Arguments

- params:

  Parameter vector c(alpha, beta)

- std_resid:

  T x k matrix of standardized residuals from first-stage GARCH

- weights:

  T-vector of observation weights (typically all 1s)

- Qbar:

  k x k unconditional correlation matrix (sample correlation of
  std_resid)

- distribution:

  "mvn" (multivariate normal). For "mvt", use
  numerical_hessian_richardson() instead as analytical MVT Hessian is
  not yet implemented.

## Value

2 x 2 Hessian matrix of the NLL with respect to (alpha, beta)

## Details

The DCC(1,1) model specifies: \$\$Q_t = (1 - \alpha - \beta) \bar{Q} +
\alpha z\_{t-1} z\_{t-1}' + \beta Q\_{t-1}\$\$ \$\$R_t =
diag(Q_t)^{-1/2} Q_t diag(Q_t)^{-1/2}\$\$

The correlation log-likelihood contribution at time t is: \$\$\ell_t =
-\frac{1}{2}\[\log\|R_t\| + z_t' R_t^{-1} z_t - z_t' z_t\]\$\$

This function computes \\\partial^2 (-\sum_t \ell_t) / \partial \theta
\partial \theta'\\ where \\\theta = (\alpha, \beta)'\\.
