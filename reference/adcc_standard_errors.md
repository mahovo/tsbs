# Compute Standard Errors for ADCC Parameters

Computes Hessian-based (or sandwich) standard errors for the 3-parameter
ADCC model: (alpha, gamma, beta), optionally with shape for MVT.

## Usage

``` r
adcc_standard_errors(
  params,
  z_matrix,
  weights,
  Qbar,
  Nbar = NULL,
  copula_dist = "mvn",
  method = c("hessian", "sandwich")
)
```

## Arguments

- params:

  Parameter vector: c(alpha, gamma, beta) or c(alpha, gamma, beta,
  shape)

- z_matrix:

  T x k matrix of copula residuals (PIT-transformed)

- weights:

  T-vector of observation weights

- Qbar:

  k x k unconditional correlation matrix

- Nbar:

  k x k average outer product of negative residuals (optional)

- copula_dist:

  "mvn" or "mvt"

- method:

  SE method: "hessian" (default) or "sandwich"

## Value

List with:

- se:

  Standard errors for each parameter

- vcov:

  Variance-covariance matrix

- hessian:

  Hessian matrix of NLL

- eigenvalues:

  Eigenvalues of Hessian (for diagnostics)

- condition_number:

  Condition number of Hessian

- param_names:

  Parameter names

- method:

  "hessian" or "sandwich"

- valid:

  Logical: are SEs valid?

## Details

For ADCC, the gamma parameter introduces asymmetric correlation
dynamics. The Hessian-based SEs may underestimate uncertainty for
high-persistence models (similar to standard DCC). Bootstrap SEs are
recommended for publication-quality inference.
