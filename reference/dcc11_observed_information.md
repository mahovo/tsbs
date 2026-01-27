# Compute Observed Information Matrix for DCC(1,1)

Compute the observed Fisher information matrix, which is the Hessian of
the negative log-likelihood evaluated at the MLE. This is a convenience
wrapper around numerical_hessian_richardson().

## Usage

``` r
dcc11_observed_information(
  params,
  std_resid,
  weights,
  Qbar,
  distribution = "mvn",
  use_reparam = FALSE,
  eps = 1e-05
)
```

## Arguments

- params:

  MLE parameter estimates (alpha, beta) or (psi, phi)

- std_resid:

  T x k matrix of standardized residuals

- weights:

  T-vector of observation weights

- Qbar:

  k x k unconditional covariance matrix

- distribution:

  "mvn" or "mvt"

- use_reparam:

  Logical: parameters in (psi, phi) space?

- eps:

  Step size for numerical differentiation (default 1e-5)

## Value

Observed information matrix (positive semi-definite if at MLE)
