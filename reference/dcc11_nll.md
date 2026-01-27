# DCC(1,1) Weighted Negative Log-Likelihood

Compute weighted NLL for DCC(1,1) model.

## Usage

``` r
dcc11_nll(
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

  Named vector: c(psi, phi,
  [ggplot2::shape](https://ggplot2.tidyverse.org/reference/aes_linetype_size_shape.html))
  or c(alpha, beta,
  [ggplot2::shape](https://ggplot2.tidyverse.org/reference/aes_linetype_size_shape.html))

- std_resid:

  T x k matrix of standardized residuals

- weights:

  T-vector of observation weights

- Qbar:

  k x k unconditional covariance

- distribution:

  "mvn" or "mvt"

- use_reparam:

  Logical: use reparameterized space?

## Value

Scalar negative log-likelihood
