# Compute Standard Errors for GOGARCH Models

Computes standard errors for GOGARCH component GARCH parameters. GOGARCH
differs from DCC/CGARCH in that it uses ICA decomposition followed by
univariate GARCH on independent components. SE computation focuses on
the component-level GARCH parameters.

## Usage

``` r
gogarch_standard_errors(
  garch_pars,
  ica_info,
  residuals,
  weights,
  distribution = "norm",
  method = c("hessian", "sandwich")
)
```

## Arguments

- garch_pars:

  List of GARCH parameters for each component: Each element is a list
  with omega, alpha1, beta1, etc.

- ica_info:

  ICA decomposition results (A, W, K matrices, S components)

- residuals:

  Original residuals matrix (T x k)

- weights:

  Observation weights (length T)

- distribution:

  Component distribution: "norm", "std", "nig", "gh"

- method:

  SE method: "hessian" (default) or "sandwich"

## Value

List with:

- component_se:

  List of SE for each component

- vcov_blocks:

  Block-diagonal vcov matrix (component-wise)

- valid:

  Logical: all SEs computed successfully

- n_components:

  Number of components

- method:

  Method used

## Details

GOGARCH models the observation vector as: Y = A \* S, where S contains
independent components each following univariate GARCH. The
log-likelihood decomposes as: \$\$LL = \sum_i LL_i(S_i; \theta_i) +
\log\|det(K)\|\$\$

Standard errors are computed independently for each component's GARCH
parameters using the component-wise Hessian. This is justified by the
independence assumption of ICA.

Note: SEs for the ICA mixing matrix A are not provided as A is typically
treated as a fixed transformation after ICA estimation.

## See also

[`estimate_garch_weighted_gogarch`](https://mahovo.github.io/tsbs/reference/estimate_garch_weighted_gogarch.md),
[`compute_gogarch_loglik_ms`](https://mahovo.github.io/tsbs/reference/compute_gogarch_loglik_ms.md)
