# Copula Gradient for DCC(1,1) - Analytical

Computes the analytical gradient of copula NLL w.r.t. (psi, phi) or
(alpha, beta) depending on use_reparam flag.

## Usage

``` r
copula_gradient(
  params,
  z_matrix,
  weights,
  Qbar,
  copula_dist = "mvn",
  use_reparam = TRUE
)
```

## Arguments

- params:

  Parameter vector: (psi, phi) if use_reparam=TRUE, (alpha, beta)
  otherwise. For MVT copula, params\[3\] is the shape parameter.

- z_matrix:

  T x k matrix of copula residuals

- weights:

  T-vector of observation weights

- Qbar:

  k x k unconditional covariance matrix

- copula_dist:

  "mvn" or "mvt"

- use_reparam:

  Logical; if TRUE, use reparameterized (psi, phi) space

## Value

Gradient vector (same length as params)
