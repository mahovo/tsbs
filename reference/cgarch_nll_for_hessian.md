# CGARCH Negative Log-Likelihood for Hessian Computation

Wrapper around copula NLL for use with numerical_hessian.

## Usage

``` r
cgarch_nll_for_hessian(
  params,
  z_matrix,
  weights,
  Qbar,
  copula_dist = "mvn",
  use_reparam = FALSE
)
```
