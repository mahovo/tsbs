# Fit a Flexible N-State Markov-Switching Vector-ARMA-GARCH Model

This function provides a user-friendly interface to fit a general
Markov-Switching model. It handles data validation, pre-processing,
calls the C++ estimation engine, and formats the results into a
well-structured object.

## Usage

``` r
fit_ms_varma_garch(
  y,
  M,
  d = 0,
  spec,
  model_type = c("univariate", "multivariate"),
  control = list(),
  parallel = FALSE,
  num_cores = 1L,
  collect_diagnostics = FALSE,
  verbose = FALSE,
  verbose_file = NULL
)
```

## Arguments

- y:

  A numeric matrix where rows are observations and columns are the time
  series variables.

- M:

  An integer \>1 specifying the number of states in the Markov chain.

- d:

  An integer specifying the order of differencing to be applied to the
  series before fitting. This order is constant across all states.
  Defaults to 0 (no differencing), which is appropriate for returns
  data.

- spec:

  A list of model specifications, one for each of the M states. See
  Details for the required structure.

- model_type:

  A character string, either "univariate" or "multivariate". Defaults to
  "univariate".

- control:

  A list of control parameters for the EM algorithm, including:

  - max_iterAn integer specifying the maximum number of EM iterations.
    Defaults to 100.

  - tolA numeric value specifying the convergence tolerance for the
    log-likelihood. Defaults to 1e-6.

  - dcc_boundary_thresholdThreshold for detecting boundary DCC
    parameters (default 0.02).

  - dcc_boundary_criterionSelection criterion: "threshold", "aic", or
    "bic" (default "bic").

  - dcc_allow_refittingLogical; if TRUE, allow refitting with constant
    correlation when boundary parameters are detected (default TRUE).

- parallel:

  Logical.

- num_cores:

  Number of cores.

- collect_diagnostics:

  Logical. Collect diagnostics or not.

- verbose:

  Logical. If TRUE, print detailed diagnostic information during
  estimation. Default is FALSE.

- verbose_file:

  Character string specifying path to file for verbose output. If NULL
  (default), verbose output goes to console. If specified, all verbose
  output is written to this file instead. Only used if verbose = TRUE.

## Value

A list containing the full results of the estimation, including model
fits for each state, the transition matrix, smoothed probabilities, and
information criteria:

- model_fitsModel fits

- PP

- log_likelihoodLog-likelihood

- smoothed_probabilitiesSmoothed probabilities

- aicAIC

- bicBIC

- num_paramsNumber of parameters (used for AIC/BIC)

- dOrder of differencing

- yData matrix

- callReturn the exact call to fit_ms_varma_garch()

- convergenceConvergence

- dcc_boundary_criterionDCC boundary criterion used

- warningsWarnings

- diagnosticsDiagnostics

## Details

The `spec` argument is a list of length `M`, where each element
`spec[[j]]` defines the model for state `j`. This state-specific list
must contain:

- For the mean model: `arma_order` (univariate) or `var_order`
  (multivariate).

- For the variance model: `garch_model` (e.g., "garch"), `garch_order`,
  `distribution`, and for multivariate models, `garch_spec_fun`
  ("dcc_modelspec", "cgarch_modelspec", or "gogarch_modelspec") and
  `garch_spec_args`.

- Starting parameters: `start_pars`, a list containing `arma_pars` (or
  `var_pars`) and `garch_pars`.

**Multivariate Model Types:**

- `dcc_modelspec`: Dynamic Conditional Correlation model. Estimates
  time-varying correlations using the DCC recursion of Engle (2002).

- `cgarch_modelspec`: Copula GARCH model. Separates marginal
  distributions from dependence structure using copula theory. Supports
  MVN and MVT copulas with parametric, empirical, or SPD PIT
  transformations.

- `gogarch_modelspec`: Generalized Orthogonal GARCH model. Uses ICA to
  extract independent components, then estimates univariate GARCH on
  each component. Time-varying covariance arises from the fixed mixing
  matrix.

For mathematical expression of the model see
[`ms_varma_garch_bs()`](https://mahovo.github.io/tsbs/reference/ms_varma_garch_bs.md).

**Parameter Counting for Information Criteria:**

The AIC and BIC are computed using the total number of estimated
parameters, which includes: mean model parameters (VAR coefficients),
univariate GARCH parameters, correlation dynamics parameters (for
DCC/CGARCH), and transition matrix parameters.

For GOGARCH models, the ICA mixing matrix is treated as a data
transformation rather than as estimated parameters for information
criteria purposes. This is consistent with the view that ICA extracts a
fixed orthogonal rotation of the data, analogous to principal
components, rather than optimizing a likelihood contribution.
Consequently, GOGARCH models will have fewer counted parameters than
DCC/CGARCH models with the same number of series, which should be
considered when comparing models across different multivariate
specifications.

**Note** For single-regime data, since M \> 1, and assuming we set M =
2:

- both states should converge to similar parameters,

- or one state should have nearly zero probability.

Check: `fit$smoothed_probabilities` - is one state dominant?

## See also

- [`ms_varma_garch_bs`](https://mahovo.github.io/tsbs/reference/ms_varma_garch_bs.md):
  Bootstrap inference for MS-VARMA-GARCH models

- [`estimate_garch_weighted_dcc`](https://mahovo.github.io/tsbs/reference/estimate_garch_weighted_dcc.md):
  DCC estimation details

- [`estimate_garch_weighted_cgarch`](https://mahovo.github.io/tsbs/reference/estimate_garch_weighted_cgarch.md):
  Copula GARCH estimation details

- [`estimate_garch_weighted_gogarch`](https://mahovo.github.io/tsbs/reference/estimate_garch_weighted_gogarch.md):
  GOGARCH estimation details
