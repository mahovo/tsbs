# Compute DCC(1,1) Hessian with Full Diagnostics

Core function for Hessian-based inference. Returns Hessian,
variance-covariance matrix, standard errors, and diagnostic information
(eigenvalues, condition number) useful for detecting the "flat beta
problem".

## Usage

``` r
dcc11_hessian(
  params,
  std_resid,
  weights,
  Qbar,
  distribution = "mvn",
  use_reparam = FALSE,
  hessian_method = "numerical",
  eps = 1e-05
)
```

## Arguments

- params:

  MLE parameter estimates c(alpha, beta) or c(alpha, beta, shape)

- std_resid:

  T x k matrix of standardized residuals

- weights:

  T-vector of observation weights

- Qbar:

  k x k unconditional correlation matrix

- distribution:

  "mvn" or "mvt"

- use_reparam:

  Logical: parameters in reparameterized (psi, phi) space?

- hessian_method:

  "numerical" (default) or "analytical"

- eps:

  Step size for numerical differentiation

## Value

List with:

- hessian:

  Hessian matrix of NLL

- info:

  Observed information matrix (= Hessian for NLL)

- vcov:

  Variance-covariance matrix (inverse of info)

- se:

  Standard errors

- eigenvalues:

  Eigenvalues of Hessian

- eigenvectors:

  Eigenvectors of Hessian

- condition_number:

  Condition number of Hessian

- params:

  Parameters

- param_names:

  Para,eter names

- method:

  "hessian"
