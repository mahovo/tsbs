# Compute Sandwich (Robust) Standard Errors

Heteroskedasticity-robust SEs using sandwich estimator: Var(theta) =
I^-1 J I^-1

## Usage

``` r
dcc11_sandwich_se(
  params,
  std_resid,
  weights,
  Qbar,
  distribution = "mvn",
  use_reparam = FALSE
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

## Value

List with:

- se:

  Standard errors (square root of diagonal of vcov)

- vcov:

  Sandwich variance-covariance matrix: (\\H^{-1} J H^{-1}\\

- bread:

  \\H^{âˆ’1}\\, inverse Hessian of the negative log-likelihood

- meat:

  J, outer product of score vectors: \\sum_t (grad_t)(grad_t)'\\

- params:

  Parameter values at which SEs were computed

- param_names:

  Names of parameters

- method:

  "sandwich"

## Details

The sandwich (robust) variance estimator is: \$\$Var(\hat{\theta}) =
H^{-1} J H^{-1}\$\$ where H is the Hessian and J is the outer product of
score contributions. This is consistent under heteroskedasticity when
the standard Hessian-based estimator may not be.
