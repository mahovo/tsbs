# Copula GARCH Weighted Estimation (Two-Stage)

Implements weighted MLE for Copula GARCH models. Stage 1: Estimate
univariate GARCH parameters for each series Stage 2: Estimate copula
dependence parameters (DCC dynamics + shape)

## Usage

``` r
estimate_garch_weighted_cgarch(
  residuals,
  weights,
  spec,
  diagnostics = NULL,
  iteration = NULL,
  state = NULL,
  verbose = FALSE,
  dcc_threshold = 0.02,
  dcc_criterion = "bic",
  force_constant = FALSE
)
```

## Arguments

- residuals:

  Matrix of residuals (T x k)

- weights:

  Vector of observation weights (length T)

- spec:

  Model specification list containing:

  - garch_spec_fun: "cgarch_modelspec"

  - garch_spec_args: list with dcc_order, dynamics, transformation,
    copula, garch_model (univariate specs)

  - start_pars: starting parameter values

  - distribution: copula distribution ("mvn" or "mvt")

- diagnostics:

  Optional diagnostics collector object

- iteration:

  Current EM iteration number

- state:

  Current regime state index

- verbose:

  Logical; print diagnostic information

- dcc_threshold:

  Threshold for DCC degeneracy detection (default 0.02)

- dcc_criterion:

  Selection criterion for constant vs dynamic ("bic", "aic",
  "threshold")

- force_constant:

  Logical; if TRUE, skip dynamic estimation

## Value

List with:

- coefficients: list with garch_pars, dcc_pars, dist_pars,
  correlation_type

- warnings: list of warning messages

- diagnostics: updated diagnostics object

## Details

The Copula GARCH model differs from DCC in how it handles the dependence
structure:

1.  **Marginal Transformation**: Standardized residuals are transformed
    to uniform \[0,1\] margins via probability integral transform (PIT).
    The transformation can be:

    - "parametric": Uses the estimated univariate distribution's CDF

    - "empirical": Uses the empirical CDF

    - "spd": Uses semi-parametric distribution (see
      [`fit_spd_transform`](https://mahovo.github.io/tsbs/reference/fit_spd_transform.md))

2.  **Copula Specification**: The uniform margins are then transformed
    according to the copula distribution:

    - "mvn": Multivariate Normal copula (Gaussian copula)

    - "mvt": Multivariate Student-t copula

3.  **Correlation Dynamics**: Same as DCC - can be "constant", "dcc", or
    "adcc" For ADCC, see
    [`adcc_recursion`](https://mahovo.github.io/tsbs/reference/adcc_recursion.md).

This function is called by the M-step of the EM algorithm in
[`fit_ms_varma_garch`](https://mahovo.github.io/tsbs/reference/fit_ms_varma_garch.md)
when `garch_spec_fun = "cgarch_modelspec"`.

## See also

- [`tsbs`](https://mahovo.github.io/tsbs/reference/tsbs.md): Main
  bootstrap function (use `garch_spec_fun = "cgarch_modelspec"`)

- [`estimate_copula_parameters_weighted`](https://mahovo.github.io/tsbs/reference/estimate_copula_parameters_weighted.md):
  Stage 2 estimation

- [`compute_pit_transform`](https://mahovo.github.io/tsbs/reference/compute_pit_transform.md):
  PIT transformation methods

- [`copula_nll`](https://mahovo.github.io/tsbs/reference/copula_nll.md):
  Copula log-likelihood

- [`adcc_recursion`](https://mahovo.github.io/tsbs/reference/adcc_recursion.md):
  ADCC dynamics

- [`estimate_garch_weighted_dcc`](https://mahovo.github.io/tsbs/reference/estimate_garch_weighted_dcc.md):
  Alternative DCC estimator
