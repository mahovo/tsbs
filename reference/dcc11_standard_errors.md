# Compute Standard Errors for DCC(1,1) Parameters

Main user-facing function for SE computation. Dispatches to appropriate
method based on the `method` argument.

## Usage

``` r
dcc11_standard_errors(
  params,
  std_resid,
  weights,
  Qbar,
  distribution = "mvn",
  use_reparam = FALSE,
  method = c("hessian", "sandwich")
)
```

## Arguments

- params:

  MLE parameter estimates

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

- method:

  SE method: "hessian" (default) or "sandwich"

## Value

List with: If `method="hessian"`:

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

If `method="sandwich"`:

- se:

  Standard errors

- vcov:

  Variance-covariance matrix (inverse of info)

- bread:

  Bread

- meatinfo:

  Meat

- vcov:

  Variance-covariance matrix (inverse of info)

- params:

  Parameters

- param_names:

  Parameter names

- method:

  "sandwich"

## Details

For reliable beta inference in high-persistence DCC models, consider
using bootstrap SEs or profile likelihood CIs from `dcc_inference.R`.
