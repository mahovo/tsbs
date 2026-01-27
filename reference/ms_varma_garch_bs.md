# Stationary Bootstrap for a General MS-VARMA-GARCH Model

Fits a flexible \\n\\-state Markov-Switching Vector ARIMA\\(p, d, q)\\-
GARCH model and then uses the estimated state sequence to perform a
stationary block bootstrap. This generates resampled time series that
preserve the state-dependent properties of the original data.

## Usage

``` r
ms_varma_garch_bs(
  x,
  n_boot = NULL,
  num_blocks = 100,
  num_boots = 100,
  M,
  d = 0,
  spec,
  model_type = c("univariate", "multivariate"),
  control = list(),
  parallel = FALSE,
  num_cores = 1L,
  return_fit = FALSE,
  collect_diagnostics = FALSE,
  verbose = FALSE,
  verbose_file = NULL
)
```

## Arguments

- x:

  A numeric matrix where rows are observations and columns are the time
  series variables.

- n_boot:

  An integer specifying the length of each bootstrapped series. Default
  is `NULL`.

- num_blocks:

  An integer specifying the number of blocks to sample for the
  bootstrap. Defaults to 100.

- num_boots:

  An integer specifying the total number of bootstrapped series to
  generate. Defaults to 100.

- M:

  An integer specifying the number of states in the Markov chain.

- d:

  An integer specifying the order of differencing for the ARIMA model.

- spec:

  A list of model specifications, one for each of the M states.

- model_type:

  A character string, either "univariate" or "multivariate".

- control:

  A list of control parameters for the EM algorithm.

- parallel:

  A logical value indicating whether to use parallel processing.

- num_cores:

  An integer specifying the number of cores for parallel processing.

- return_fit:

  If `TRUE`, [`tsbs()`](https://mahovo.github.io/tsbs/reference/tsbs.md)
  will return model fit when `bs_type = "ms_varma_garch"`. Default is
  `return_fit = FALSE`. If `bs_type = "ms_varma_garch"`,
  `return_fit = TRUE` and `collect_diagnostics = TRUE` diagnostics can
  be extracted from `result$fit$diagnostics`. See
  `vignette("Diagnostics", package = "tsbs")`.

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

A list of bootstrapped time series matrices.

## Details

If `n_boot` is set, the last block will be trimmed when necessary. If
`n_boot` is not set, and `num_blocks` is set, the length of each
bootstrap series will be determined by the number of blocks and the
random lengths of the individual blocks for that particular series. If
neither `n_boot` nor `num_blocks` is set, `n_boot` will default to the
number of rows in `x` and the last block will be trimmed when necessary.

The fitted model is defined as:

Let \\y_t\\ be the \\k \times 1\\ vector of observations at time \\t\\.
The model assumes that the data-generating process is governed by a
latent (unobserved) state variable, \\S_t\\, which follows a first-order
Markov chain with \\M\\ states.

1.  **State Process**: The evolution of the state is described by the
    \\M \times M\\ transition probability matrix \\P\\, where the
    element \\p\_{ij}\\ is the probability of transitioning from state
    \\i\\ to state \\j\\: \$\$p\_{ij} = P(S_t = j \| S\_{t-1} = i)\$\$
    The matrix \\P\\ is structured such that \\P\_{ij} = p\_{ij}\\, and
    its rows sum to one: \\\sum\_{j=1}^{M} p\_{ij} = 1\\ for all \\i=1,
    \dots, M\\.

2.  **Observation Process**: Conditional on the system being in state
    \\S_t = j\\, each of the \\k\\ time series, \\y\_{i,t}\\ for \\i=1,
    \dots, k\\, is assumed to follow an independent ARIMA(\\p_j, d_j,
    q_j\\)-GARCH(\\q'\_j, p'\_j\\) process. The parameters for both the
    mean and variance equations are specific to the state \\j\\.

    - **Mean Equation (ARIMA)**: \$\$\phi_j(L)(1-L)^{d_j} (y\_{i,t} -
      \mu_j) = \theta_j(L) \varepsilon\_{i,t}\$\$ where \\\phi_j(L)\\
      and \\\theta_j(L)\\ are the AR and MA lag polynomials, \\\mu_j\\
      is the mean, and \\\varepsilon\_{i,t}\\ is the innovation term,
      all specific to state \\j\\.

    - **Variance Equation (GARCH)**: The innovations have a conditional
      variance \\\sigma\_{i,t}^2\\ that evolves according to:
      \$\$\varepsilon\_{i,t} = \sigma\_{i,t} z\_{i,t}, \quad z\_{i,t}
      \sim \mathcal{N}(0,1)\$\$ \$\$\sigma\_{i,t}^2 = \omega_j +
      \sum\_{l=1}^{q'\_j} \alpha\_{j,l} \varepsilon\_{i,t-l}^2 +
      \sum\_{l=1}^{p'\_j} \beta\_{j,l} \sigma\_{i,t-l}^2\$\$ where
      \\\omega_j, \alpha\_{j,l}, \beta\_{j,l}\\ are the GARCH parameters
      for state \\j\\.

Let \\\Psi_j = \\\mu_j, \phi_j, \theta_j, \omega_j, \alpha_j,
\beta_j\\\\ be the complete set of ARIMA-GARCH parameters for state
\\j\\, and let \\\Psi = \\\Psi_1, \dots, \Psi_M, P\\\\ be the full
parameter set for the entire model. The EM algorithm using a Hamilton
Filter & Kim Smoother for the E-step is used to find the Maximum
Likelihood Estimate (MLE) of \\\Psi\\.

The tsbs package uses different optimization strategies for DCC models
depending on the model order:

**DCC(1,1) - Reparameterized Optimization**

For the common DCC(1,1) case, we use a reparameterization that
transforms the constrained problem into an unconstrained one:

- Original: \\\alpha \in (0,1), \beta \in (0,1), \alpha + \beta \< 1\\

- Reparameterized: \\persistence \in (0,1), ratio \in (0,1)\\

where \\persistence = \alpha + \beta\\ and \\ratio = \alpha/(\alpha +
\beta)\\.

This eliminates the stationarity constraint since \\\alpha + \beta =
persistence \< 1\\ is automatically satisfied by the box constraint on
persistence.

Benefits:

- No penalty function discontinuities

- More stable optimization near high-persistence regions

- Eliminates "dcc_penalty" warnings from optimizer exploration

**DCC(p,q) with max(p,q) \> 1 - Penalty Method**

For higher-order DCC models, reparameterization becomes significantly
more complex (requiring softmax distributions over parameter vectors).
We therefore use the standard penalty method:

- Box constraints: \\\alpha_j, \beta_j \in (\epsilon, 1-\epsilon)\\

- Stationarity enforced via penalty when \\\sum \alpha + \sum \beta \geq
  1\\

This may result in optimizer instability warnings for models with high
persistence.

## References

Natatou Moutari, D. et al. (2021). Dependence Modeling and Risk
Assessment of a Financial Portfolio with ARMA-APARCH-EVT models based on
HACs. [arXiv:2105.09473](http://arxiv.org/abs/2105.09473)
