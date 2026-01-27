# Full DCC(1,1) Gradient Function for Optimizer

Compute gradient of weighted NLL w.r.t. all DCC parameters. This is the
main interface for use with optim().

## Usage

``` r
dcc11_gradient(
  params,
  std_resid,
  weights,
  Qbar,
  distribution = "mvn",
  use_reparam = TRUE
)
```

## Arguments

- params:

  Named vector of parameters to optimize

- std_resid:

  T x k matrix of standardized residuals

- weights:

  T-vector of observation weights

- Qbar:

  k x k unconditional covariance

- distribution:

  "mvn" or "mvt"

- use_reparam:

  Logical: use reparameterized (psi, phi) space?

## Value

Vector of gradients (same length and order as params)
