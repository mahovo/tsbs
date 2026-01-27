# DCC(1,1) Weighted NLL Gradient (Reparameterized)

Compute gradient of weighted NLL w.r.t. (psi, phi).

## Usage

``` r
dcc11_gradient_reparam(
  psi,
  phi,
  z,
  weights,
  Qbar,
  distribution = "mvn",
  shape = NULL
)
```

## Arguments

- psi:

  Unconstrained persistence parameter

- phi:

  Unconstrained ratio parameter

- z:

  T x k matrix of standardized residuals

- weights:

  T-vector of observation weights

- Qbar:

  k x k unconditional covariance

- distribution:

  "mvn" or "mvt"

- shape:

  Degrees of freedom (MVT only)

## Value

Named vector c(psi, phi)
