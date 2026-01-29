# Advanced Bootstrap for Time Series

Generates bootstrap replicates of multivariate or (multiple) univariate
time series. Supports moving and stationary block, HMM, MSVAR, MS VARMA
GARCH and wild bootstrap types.

## Usage

``` r
tsbs(
  x,
  n_boot = NULL,
  block_length = NULL,
  bs_type = c("moving", "stationary", "hmm", "msvar", "ms_varma_garch", "wild"),
  block_type = c("overlapping", "non-overlapping", "tapered"),
  taper_type = c("cosine", "bartlett", "tukey"),
  tukey_alpha = 0.5,
  num_blocks = NULL,
  num_boots = 100L,
  func = mean,
  apply_func_to = c("cols", "all"),
  p_method = c("1/n", "plugin", "cv"),
  p = NULL,
  overlap = TRUE,
  num_states = 2L,
  d = 0,
  spec = NULL,
  model_type = c("univariate", "multivariate"),
  control = list(),
  parallel = FALSE,
  num_cores = 1L,
  model_func = default_model_func,
  score_func = mse,
  stationary_max_percentile = 0.99,
  stationary_max_fraction_of_n = 0.1,
  distribution = "gaussian",
  variance_model = "sGARCH",
  micro_block_length = 1L,
  regime_basis = "market",
  seed = NULL,
  return_fit = FALSE,
  return_diagnostics = FALSE,
  fail_mode = c("predictably", "gracefully"),
  ...
)
```

## Arguments

- x:

  Numeric vector, matrix, data frame or time series observations (rows =
  time points, cols = variables).

- n_boot:

  Integer, optional desired length of each bootstrap replicate. See
  details below.

- block_length:

  Integer, length of each block. If `NULL`, an automatic heuristic is
  used: The method first calculates the average absolute first-order
  autocorrelation (\\\rho_1\\) across all time series columns. A
  candidate block length is then calculated based on this persistence
  measure using the formula \\⌊10/(1−\rho_1)⌋\\. The final block length
  is constrained to be between 5 and the square root of the series
  length. For stationary bootstrap, `block_length` is the expected block
  length, when `p_method="1/n"`.

- bs_type:

  Bootstrap type. Character string: One of `"moving"`, `"stationary"`,
  `"hmm"`, `"msvar"`, `"ms_varma_garch"`, or `"wild"`. `"msvar"` is a
  lightweight special case of `"ms_varma_garch"`. See details below.

- block_type:

  Block type. Character string: One of `"non-overlapping"`,
  `"overlapping"`, or `"tapered"`. Only affects `bs_type="moving"` and
  `bs_type="stationary"`. `block_type="tapered"` will smooth out
  transitions between blocks. It follows that `block_type="tapered"` can
  not be used when `block_length=1`.

- taper_type:

  Tapering window function. Character. One of `"cosine"`, `"bartlett"`,
  or `"tukey"`. Only affects `block_type="tapered"`. See details below.

- tukey_alpha:

  numeric, alpha parameter for `taper_type = "tukey"`. \\\alpha \in \[0,
  1\]\\ is the fraction of the window length tapered at each end, using
  indices \\0 \dots n−1\\, where \\n\\ is the block length.

- num_blocks:

  Integer number of blocks per bootstrap replicate.

- num_boots:

  Integer number of bootstrap replicates.

- func:

  A function to apply to each bootstrap replicate. The function's
  expected input depends on the `apply_func_to` parameter:

  - If `apply_func_to = "cols"`: `func` receives a single numeric vector
    (one column at a time) and should return a scalar or numeric vector.
    Examples: `mean`, `sd`, `quantile`.

  - If `apply_func_to = "all"`: `func` receives the entire bootstrap
    replicate as a numeric matrix (rows = time points, columns =
    variables) and can return:

    - A numeric scalar or vector (e.g., portfolio weights)

    - A named list of numeric values (e.g., multiple summary statistics)

    Examples: portfolio optimization functions, functions returning
    multiple statistics.

  The function outputs are collected in `func_outs` (one entry per
  replicate), and `func_out_means` contains the element-wise average
  across replicates. For list-returning functions, each named element is
  averaged separately. Default is `mean`.

- apply_func_to:

  Character string: `"cols"` to apply columnwise or `"all"` to apply on
  the full data frame or matrix.

- p_method:

  Character string to choose method for stationary bootstrap parameter:
  `"1/n"`, `"plugin"`, or `"cv"`.

- p:

  numeric \\p \in (0, 1)\\. Probability parameter for the geometric
  block length (used in Stationary Bootstrap).

- overlap:

  Logical indicating if overlapping blocks are allowed.

- num_states:

  Integer number of states for HMM, MSVAR, or MS-VARMA-GARCH models.

- d:

  Integer differencing order for the MS-VARMA-GARCH model.

- spec:

  A list defining the state-specific models for the MS-VARMA-GARCH
  bootstrap (bs_type = "ms_varma_garch"). This argument is required for
  this bootstrap type and must be a list of length `num_states`. Each
  element of the list is itself a list specifying the model for that
  state. The structure depends on whether the model is univariate or
  multivariate. Fitting relies on the tsgarch or tsmarch packages,
  respectively. See details below.

- model_type:

  Character string for MS-VARMA-GARCH: `"univariate"` or
  `"multivariate"`.

- control:

  A list of control parameters for the MS-VARMA-GARCH EM algorithm.

- parallel:

  Parallelize computation? `TRUE` or `FALSE`.

- num_cores:

  Number of cores when `parallel=TRUE`.

- model_func:

  Model-fitting function for cross-validation. See
  [k_fold_cv_ts](https://mahovo.github.io/tsbs/reference/k_fold_cv_ts.md).

- score_func:

  Scoring function for cross-validation.

- stationary_max_percentile:

  Stationary max percentile.

- stationary_max_fraction_of_n:

  Stationary max fraction of n.

- distribution:

  Character string specifying the conditional distribution for HMM
  bootstrap. One of:

  - `"gaussian"`: Gaussian HMM via depmixS4 (default, backward
    compatible)

  - `"norm"`, `"std"`, `"sstd"`, `"snorm"`, `"ged"`, `"sged"`:
    MSGARCH-based with GARCH dynamics

  - `"norm_raw"`, `"std_raw"`, `"sstd_raw"`: Raw HMM without GARCH layer

  Only used when `bs_type = "hmm"`. See Details for more information.

- variance_model:

  Character string specifying the GARCH model for MSGARCH-based HMM. One
  of `"sGARCH"`, `"eGARCH"`, `"gjrGARCH"`, or `"tGARCH"`. Default is
  `"sGARCH"`. Only used when `bs_type = "hmm"` and `distribution` is an
  MSGARCH distribution.

- micro_block_length:

  Integer specifying the block length for within-state sampling. Use 1
  for iid sampling (default), or \>1 to preserve local temporal
  dependence within states. Only used when `bs_type = "hmm"`.

- regime_basis:

  For multivariate HMM bootstrap, specifies how to identify regimes. One
  of:

  - `"market"`: Use equal-weighted average of all columns (default)

  - `"first_pc"`: Use first principal component

  - Integer: Use the specified column index

  Only used when `bs_type = "hmm"` with multivariate data.

- seed:

  Integer random seed for reproducibility of HMM bootstrap. Only used
  when `bs_type = "hmm"`.

- return_fit:

  If `TRUE`, `tsbs()` will return model fit. Default is
  `return_fit = FALSE`. If `return_fit = TRUE` and
  `bs_type = "ms_varma_garch"`, diagnostics can be extracted from
  `result$fit$diagnostics`. See
  [`?fit_ms_varma_garch`](https://mahovo.github.io/tsbs/reference/fit_ms_varma_garch.md).

- return_diagnostics:

  If `TRUE`, returns bootstrap diagnostics including block composition,
  length statistics, and quality metrics. Diagnostics are returned in
  `result$diagnostics` as a `tsbs_diagnostics` object. Use
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods on
  the diagnostics for detailed analysis. Default is `FALSE`.

- fail_mode:

  Character string. One of `"predictably"` (development/ debugging -
  fail fast on any unexpected behavior) or `"gracefully"`
  (production/robust pipelines - continue despite validation errors).

## Value

A list containing:

- bootstrap_series:

  List of bootstrap replicate matrices.

- func_outs:

  List of computed function outputs for each replicate.

- func_out_means:

  Mean of the computed outputs across replicates.

- fit:

  If `return_fit = TRUE`, a list containing model fit.

- diagnostics:

  If `return_diagnostics = TRUE`, a `tsbs_diagnostics` object containing
  block composition, statistics, and quality metrics. See
  [`compute_bootstrap_diagnostics`](https://mahovo.github.io/tsbs/reference/compute_bootstrap_diagnostics.md).

## Details

### Bootstrap types

See documentation for the individual bootstrap functions:

- For moving or stationary block bootstraps, see
  [`blockBootstrap()`](https://mahovo.github.io/tsbs/reference/blockBootstrap.md)

- [`hmm_bootstrap()`](https://mahovo.github.io/tsbs/reference/hmm_bootstrap.md)

- [`msvar_bootstrap()`](https://mahovo.github.io/tsbs/reference/msvar_bootstrap.md)

- [`ms_varma_garch_bs()`](https://mahovo.github.io/tsbs/reference/ms_varma_garch_bs.md)

- [`wild_bootstrap()`](https://mahovo.github.io/tsbs/reference/wild_bootstrap.md)

#### `bs_type="moving"`

If `n_boot` is set, the last block will be truncated when necessary to
match the length (`n_boot`) of the bootstrap series. If `n_boot` is not
set, `block_length` and `num_blocks` must be set, and `n_boot` will
automatically be set to `block_length * num_blocks`.

#### `bs_type="stationary"`, `bs_type="hmm"`, `bs_type="msvar"` If `n_boot` is

set, the last block will be truncated when necessary to match the length
(`n_boot`) of the bootstrap series. This is the only way to ensure equal
length of all bootstrap series, as the length of each block is random.
If `n_boot` is not set, `num_blocks` must be set, and the length of each
bootstrap series will be determined by the number of blocks and the
lengths of the individual blocks for that particular series. Truncation
ensures equal length of bootstrapped series for bootstrap types with
random block length. For stationary bootstrap, `block_length` is the
expected block length, when `p_method="1/n"`.

#### `bs_type="hmm"`

The HMM bootstrap uses a Hidden Markov Model to capture regime-switching
behavior in the data. The bootstrap procedure is semi-parametric:

- **Parametric**: State sequences are simulated from the fitted Markov
  chain transition matrix

- **Nonparametric**: Observations are resampled from empirical
  state-specific pools

**Distribution Options:**

The `distribution` parameter controls which backend is used:

|                  |             |           |                             |
|------------------|-------------|-----------|-----------------------------|
| **Distribution** | **Backend** | **GARCH** | **Description**             |
| `"gaussian"`     | depmixS4    | No        | Original behavior (default) |
| `"norm"`         | MSGARCH     | Yes       | Normal with GARCH           |
| `"std"`          | MSGARCH     | Yes       | Student-t with GARCH        |
| `"sstd"`         | MSGARCH     | Yes       | Skew Student-t with GARCH   |
| `"snorm"`        | MSGARCH     | Yes       | Skew normal with GARCH      |
| `"ged"`          | MSGARCH     | Yes       | GED with GARCH              |
| `"sged"`         | MSGARCH     | Yes       | Skew GED with GARCH         |
| `"norm_raw"`     | EM/fGarch   | No        | Normal, no GARCH            |
| `"std_raw"`      | EM/fGarch   | No        | Student-t, no GARCH         |
| `"sstd_raw"`     | EM/fGarch   | No        | Skew Student-t, no GARCH    |

**MSGARCH Distributions** (`"norm"`, `"std"`, `"sstd"`, etc.): These use
the MSGARCH package to fit Markov-switching GARCH models following Haas,
Mittnik & Paolella (2004). The skewed distributions use the Fernández &
Steel (1998) transformation. Requires the MSGARCH package.

**Raw Distributions** (`"norm_raw"`, `"std_raw"`, `"sstd_raw"`): These
fit a standard HMM with state-specific emission distributions but
without GARCH dynamics. Uses an EM algorithm with the fGarch package for
density evaluation.

**Multivariate Data:** For multivariate data, the `regime_basis`
parameter determines how regimes are identified:

- `"market"`: Fits the regime model to the equal-weighted average

- `"first_pc"`: Fits to the first principal component

- Integer column index: Fits to a specific column

Bootstrap samples use synchronized sampling: the same time indices are
sampled for all variables, preserving cross-sectional dependence.

**Micro-Block Sampling:** The `micro_block_length` parameter enables
block sampling within states. With `micro_block_length = 1` (default),
observations are sampled iid from state pools. With
`micro_block_length > 1`, consecutive blocks are sampled, preserving
some local autocorrelation while still allowing the Markov chain to
drive regime dynamics.

If `n_boot` is set, bootstrap series are truncated to this length.
Otherwise, the length matches the original series.

**Dependencies:**

- `"gaussian"`: Requires depmixS4

- MSGARCH distributions: Requires MSGARCH

- Raw distributions: Requires fGarch (except `"norm_raw"`)

**References:**

- Hamilton, J. D. (1989). A New Approach to the Economic Analysis of
  Nonstationary Time Series and the Business Cycle. Econometrica, 57(2),
  357-384.

- Haas, M., Mittnik, S., & Paolella, M. S. (2004). A New Approach to
  Markov-Switching GARCH Models. Journal of Financial Econometrics, 2,
  493-530.

- Fernández, C., & Steel, M. F. (1998). On Bayesian modeling of fat
  tails and skewness. Journal of the American Statistical Association,
  93(441), 359-371.

#### `bs_type="wild"`

`n_boot`, `block_length` and `num_blocks` are ignored. The length of the
bootstrap series is always the same as the original series.

#### `bs_type = "ms_varma_garch"`

`spec` list structure:

- **For univariate models (one-column x):** Each element spec\[\[j\]\]
  must be a list with the following components. For a complete list of
  all possible arguments (e.g., for different GARCH model flavors),
  please refer to the documentation for the tsgarch package
  ([`?tsgarch`](https://rdrr.io/pkg/tsgarch/man/tsgarch-package.html)
  and
  [`vignette("garch_models", package = "tsgarch")`](https://cran.rstudio.com/web/packages/tsgarch/vignettes/garch_models.pdf)).

  - `arma_order`: A numeric vector c(p,q) for the ARMA(p,q) order.

  - `garch_model`: A character string for the GARCH model type (e.g.,
    "garch", "egarch").

  - `garch_order`: A numeric vector c(g,h) for the GARCH(g,h) order.

  - `distribution`: A character string for the conditional distribution.

    - `"norm"` = Normal

    - `"snorm"` = Skew normal

    - `"std"` = Student t

    - `"sstd"` = Skew Student

    - `"ged"` = Generalized error

    - `"ged"` = Skew generalized error

    - `"ghyp"` = Generalized hyperbolic

    - `"ghst"` = Generalized hyperbolic skew Student

    - `"jsu"` = Johnson reparameterized SU

    See <https://www.nopredict.com/packages/tsdistributions>.

  - `start_pars`: A list containing the starting parameters for the
    optimization, with two named elements:

    - `arma_pars`: A named numeric vector for the ARMA parameters (e.g.,
      c(ar1 = 0.5)).

    - `garch_pars`: A named list for the GARCH parameters (e.g.,
      list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8)).

    - `dist_pars`: A named list for the starting values for univariate
      dist params (e.g. list(skew = 1.0, shape = 6.0)).

- **For multivariate models (multi-column x):** Each element
  spec\[\[j\]\] must be a list with the following components. For a
  complete list of all possible arguments, please refer to the
  documentation for the relevant tsmarch specification function. For a
  DCC model see
  [`?dcc_modelspec.tsgarch.multi_estimate`](https://rdrr.io/pkg/tsmarch/man/dcc_modelspec.html),
  for a Copula GARCH model see
  [`?cgarch_modelspec.tsgarch.multi_estimate`](https://rdrr.io/pkg/tsmarch/man/cgarch_modelspec.html).
  (See also
  [`?tsmarch`](https://rdrr.io/pkg/tsmarch/man/tsmarch-package.html),
  [`vignette("tsmarch_demo", package = "tsmarch")`](https://cran.rstudio.com/web/packages/tsmarch/vignettes/tsmarch_demo.html)
  and
  [`vignette("feasible_multivariate_garch", package = "tsmarch")`](https://cran.rstudio.com/web/packages/tsmarch/vignettes/feasible_multivariate_garch.pdf).)

  - `var_order`: An integer p for the VAR(p) order.

  - `garch_spec_fun`: A character string with the name of the tsmarch
    specification function to use (e.g., "dcc_modelspec",
    "gogarch_modelspec").

  - `distribution`: The multivariate distribution. Valid choices are
    "mvn" (for multivariate normal) and "mvt" (for multivariate
    Student's t).

  - `garch_spec_args`: A list of arguments to pass to the function
    specified in garch_spec_fun. Key arguments include:

    - For models like DCC or CGARCH, this list must also contain a
      `garch_model` definition for the underlying univariate GARCH fits.
      Note: With the exception of the Copula model, the marginal
      distributions of the univariate GARCH models should always be
      Normal, irrespective of whether a multivariate Normal or Student
      is chosen as the DCC model distribution. There are no checks
      performed for this and it is up to the user to ensure that this is
      the case.

  - `start_pars`: A list containing the starting parameters, with two
    named elements:

    - `var_pars`: A numeric vector for the VAR parameters, stacked
      column by column (intercept, lags for eq1, intercept, lags for
      eq2, ...).

    - `garch_pars`: A list of named lists for the multivariate GARCH
      parameters, e.g.
      `replicate(k, list(omega = 0.05, alpha1 = 0.08, beta1 = 0.85), simplify = FALSE)`.

    - `dcc pars`: A named list, e.g.
      `list(alpha_1 = 0.05, beta_1 = 0.90)`.

    - `dist pars`: A named list, e.g. `list(shape = 8.0)`.

#### DCC (Dynamic Conditional Correlation) Models

DCC models capture time-varying correlations between multiple series by
modeling the correlation dynamics directly on standardized residuals.
This is the standard approach for multivariate GARCH modeling in the
tsbs package.

To use DCC, set `garch_spec_fun = "dcc_modelspec"` in your spec. See
[`tsmarch::dcc_modelspec()`](https://rdrr.io/pkg/tsmarch/man/dcc_modelspec.generic.html)
for the tsmarch specification function.

**DCC-specific `garch_spec_args`:**

- `dcc_order`: Integer vector c(p, q) for DCC(p, q) order. Default is
  c(1, 1). Note: DCC(1,1) uses analytical gradients; higher orders use
  numerical differentiation.

- `dynamics`: Character. Correlation dynamics type:

  - `"constant"`: Time-invariant correlation (reduces to CCC model)

  - `"dcc"`: Standard DCC dynamics

  - `"adcc"`: Asymmetric DCC (negative shocks have larger impact)

- `distribution`: Character. Multivariate distribution:

  - `"mvn"`: Multivariate Normal

  - `"mvt"`: Multivariate Student-t (adds shape parameter)

**Important**: For DCC models, the univariate GARCH distributions should
always be Normal (`distribution = "norm"`), regardless of whether the
multivariate distribution is MVN or MVT. The multivariate distribution
applies to the joint model, not the marginals.

**Example DCC specification:**

    spec_dcc <- list(
      list(
        var_order = 1,
        garch_spec_fun = "dcc_modelspec",
        distribution = "mvn",
        garch_spec_args = list(
          dcc_order = c(1, 1),
          dynamics = "dcc",
          distribution = "mvn",
          garch_model = list(univariate = list(
            list(model = "garch", garch_order = c(1, 1), distribution = "norm"),
            list(model = "garch", garch_order = c(1, 1), distribution = "norm")
          ))
        ),
        start_pars = list(
          var_pars = rep(0, 6),
          garch_pars = list(
            list(omega = 0.05, alpha1 = 0.08, beta1 = 0.85),
            list(omega = 0.05, alpha1 = 0.08, beta1 = 0.85)
          ),
          dcc_pars = list(alpha_1 = 0.05, beta_1 = 0.90),
          dist_pars = NULL
        )
      )
    )

For ADCC (asymmetric correlation dynamics), use `dynamics = "adcc"` and
include gamma parameters in `dcc_pars`:

    dcc_pars = list(alpha_1 = 0.05, gamma_1 = 0.03, beta_1 = 0.90)

## References

Harris, F. J. (1978). "On the Use of Windows for Harmonic Analysis with
the Discrete Fourier Transform." (Proceedings of the IEEE, 66(1), pp.
66f)

## See also

- [`tsmarch::dcc_modelspec()`](https://rdrr.io/pkg/tsmarch/man/dcc_modelspec.generic.html)
  for creating DCC specifications

- [`estimate_garch_weighted_dcc()`](https://mahovo.github.io/tsbs/reference/estimate_garch_weighted_dcc.md)
  for the weighted estimation method

- [`vignette("cgarch_vs_dcc", package = "tsbs")`](https://mahovo.github.io/tsbs/articles/cgarch_vs_dcc.md)
  for model comparison

### Copula GARCH (CGARCH) Models

Copula GARCH models provide an alternative to DCC for modeling dynamic
correlations in multivariate time series. While DCC directly models
correlation dynamics, CGARCH separates marginal distributions from the
dependence structure using copulas.

To use CGARCH instead of DCC, set `garch_spec_fun = "cgarch_modelspec"`
in your spec. See
[`tsmarch::cgarch_modelspec()`](https://rdrr.io/pkg/tsmarch/man/cgarch_modelspec.generic.html)
for the tsmarch specification function.

**Key differences from DCC:**

- Uses Probability Integral Transform (PIT) to transform standardized
  residuals to uniform margins before modeling dependence

- Supports multiple transformation methods: "parametric", "empirical",
  or "spd" (semi-parametric distribution)

- Copula distribution can be "mvn" (Gaussian) or "mvt" (Student-t)

- More flexible tail dependence modeling with Student-t copula

**CGARCH-specific `garch_spec_args`:**

- `transformation`: Character. Method for PIT transformation:

  - `"parametric"`: Uses fitted univariate distribution CDFs

  - `"empirical"`: Uses empirical CDFs (non-parametric)

  - `"spd"`: Semi-parametric distribution combining kernel density in
    the center with parametric tails

- `copula`: Character. Copula distribution: `"mvn"` or `"mvt"`

- `dynamics`: Character. Correlation dynamics type:

  - `"constant"`: Time-invariant correlation

  - `"dcc"`: Dynamic Conditional Correlation

  - `"adcc"`: Asymmetric DCC (captures leverage effects where negative
    returns have larger impact on correlation)

- `dcc_order`: Integer vector c(p, q) for DCC/ADCC order

**Example CGARCH specification:**

    spec_cgarch <- list(
      list(
        var_order = 1,
        garch_spec_fun = "cgarch_modelspec",
        distribution = "mvn",
        garch_spec_args = list(
          dcc_order = c(1, 1),
          dynamics = "dcc",
          transformation = "parametric",
          copula = "mvn",
          garch_model = list(univariate = list(
            list(model = "garch", garch_order = c(1, 1), distribution = "norm"),
            list(model = "garch", garch_order = c(1, 1), distribution = "norm")
          ))
        ),
        start_pars = list(
          var_pars = rep(0, 6),
          garch_pars = list(
            list(omega = 0.05, alpha1 = 0.08, beta1 = 0.85),
            list(omega = 0.05, alpha1 = 0.08, beta1 = 0.85)
          ),
          dcc_pars = list(alpha_1 = 0.05, beta_1 = 0.90),
          dist_pars = NULL
        )
      )
    )

For ADCC (asymmetric correlation dynamics), use `dynamics = "adcc"` and
include gamma parameters:

    garch_spec_args = list(
      dynamics = "adcc",
      ...
    ),
    start_pars = list(
      ...
      dcc_pars = list(alpha_1 = 0.05, gamma_1 = 0.02, beta_1 = 0.90)
    )

- [`tsmarch::dcc_modelspec()`](https://rdrr.io/pkg/tsmarch/man/dcc_modelspec.generic.html)
  for creating DCC specifications

- [`tsmarch::cgarch_modelspec()`](https://rdrr.io/pkg/tsmarch/man/cgarch_modelspec.generic.html)
  for creating CGARCH specifications

- [`estimate_garch_weighted_dcc()`](https://mahovo.github.io/tsbs/reference/estimate_garch_weighted_dcc.md)
  for DCC weighted estimation

- [`estimate_garch_weighted_cgarch()`](https://mahovo.github.io/tsbs/reference/estimate_garch_weighted_cgarch.md)
  for CGARCH weighted estimation

- [`compute_pit_transform()`](https://mahovo.github.io/tsbs/reference/compute_pit_transform.md)
  for PIT transformation details

- [`adcc_recursion()`](https://mahovo.github.io/tsbs/reference/adcc_recursion.md)
  for asymmetric DCC dynamics

- [`vignette("cgarch_vs_dcc", package = "tsbs")`](https://mahovo.github.io/tsbs/articles/cgarch_vs_dcc.md)
  for model comparison

- `vignette("Diagnostics", package = "tsbs")` for diagnostic system

**GOGARCH-specific `garch_spec_args`:**

- `model`: Character. GARCH model type for components. Default
  `"garch"`.

- `order`: Integer vector c(p, q) for GARCH order. Default `c(1, 1)`.

- `ica`: Character. ICA algorithm:

  - `"radical"`: RADICAL algorithm (recommended, robust to outliers)

  - `"fastica"`: FastICA algorithm (currently not supported by tsbs,
    waiting for tsmarch v1.0.1 to appear on CRAN...)

- `components`: Integer. Number of ICA components to extract. Default
  equals number of series. Set lower for dimension reduction.

- `lambda_range`: Numeric vector c(min, max) for GH lambda parameter.

- `shape_range`: Numeric vector c(min, max) for GH shape parameter.

**GOGARCH Distributions:**

- `"norm"`: Normal distribution (simplest, fastest)

- `"nig"`: Normal Inverse Gaussian (captures skewness and kurtosis)

- `"gh"`: Generalized Hyperbolic (most flexible, includes NIG as special
  case)

**Example GOGARCH specification (Normal):**

    spec_gogarch <- list(
      list(
        var_order = 0,
        garch_spec_fun = "gogarch_modelspec",
        distribution = "norm",
        garch_spec_args = list(
          model = "garch",
          order = c(1, 1),
          ica = "radical",
          components = 2
        ),
        start_pars = list(
          var_pars = NULL,
          garch_pars = list(
            list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
            list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8)
          ),
          dist_pars = NULL
        )
      )
    )

**Example GOGARCH specification (NIG distribution):**

    spec_gogarch_nig <- list(
      list(
        var_order = 0,
        garch_spec_fun = "gogarch_modelspec",
        distribution = "nig",
        garch_spec_args = list(
          model = "garch",
          order = c(1, 1),
          ica = "radical",
          components = 2
        ),
        start_pars = list(
          var_pars = NULL,
          garch_pars = list(
            list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
            list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8)
          ),
          dist_pars = list(skew = 0, shape = 1)
        )
      )
    )

**Example GOGARCH specification (GH distribution):**

    spec_gogarch_gh <- list(
      list(
        var_order = 0,
        garch_spec_fun = "gogarch_modelspec",
        distribution = "gh",
        garch_spec_args = list(
          model = "garch",
          order = c(1, 1),
          ica = "radical",
          components = 2,
          lambda_range = c(-5, 5),
          shape_range = c(0.1, 25)
        ),
        start_pars = list(
          var_pars = NULL,
          garch_pars = list(
            list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
            list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8)
          ),
          dist_pars = list(skew = 0, shape = 1, lambda = -0.5)
        )
      )
    )

- [`tsmarch::gogarch_modelspec()`](https://rdrr.io/pkg/tsmarch/man/gogarch_modelspec.html)
  for creating GOGARCH specifications

- [`estimate_garch_weighted_gogarch()`](https://mahovo.github.io/tsbs/reference/estimate_garch_weighted_gogarch.md)
  for GOGARCH weighted estimation

- [`tsmarch::radical()`](https://rdrr.io/pkg/tsmarch/man/radical.html)
  for the RADICAL ICA algorithm

- [`compute_gogarch_loglik_ms()`](https://mahovo.github.io/tsbs/reference/compute_gogarch_loglik_ms.md)
  for GOGARCH log-likelihood computation

- [`vignette("cgarch_vs_dcc", package = "tsbs")`](https://mahovo.github.io/tsbs/articles/cgarch_vs_dcc.md)
  for model comparison

### Diagnostics for MS VARMA GARCH

Diagnostic System: A comprehensive diagnostic system monitors
convergence and numerical stability during estimation. Additional
parameters when `bs_type = "ms_varma_garch"`

- `collect_diagnostics` Logical. Enable diagnostics with
  `collect_diagnostics = TRUE`. Diagnostics can only be extracted if
  `return_fit = TRUE`.

- `verbose` Logical. If TRUE, print detailed diagnostic information
  during estimation. Default is FALSE.

- `verbose_file` Character string specifying path to file for verbose
  output. If NULL (default), verbose output goes to console. If
  specified, all verbose output is written to this file instead. Only
  used if verbose = TRUE.

See details in `vignette("Diagnostics", package = "tsbs")`.

### `taper_type` when `block_type="tapered"`

For block length \\n\\, and index \\i\\ within the block, and \\0 \leq i
\leq n\\,

- `"cosine"`: Hann window \$\$w(i) = \dfrac{1}{2} \left(1 - \cos \left(
  \dfrac{2 \pi i}{n - i} \right ) \right)\$\$

- `"bartlett"`: Triangular window \$\$w(i) = 1 - \left \| \dfrac{i -
  (n - 1) / 2}{(n - 1) / 2} \right \|,\\ \\ \\ \alpha \in \[0, 1\]\$\$

- `"tukey"`: Cosine tapered window. At \\\alpha = 0\\ it becomes
  rectangular, and at \\\alpha = 1\\ it becomes a Hann window.

  We implement the Tukey (tapered cosine) window using the convention
  \\i = 0, \dots, n-1\\ for a window of length \\n\\ and taper parameter
  \\\alpha \in \[0,1\]\\, consistent with common numerical libraries
  such as MATLAB and SciPy (see MATLAB tukeywin, SciPy
  scipy.signal.windows.tukey). After index shift and reflection this is
  algebraically equivalent to the classical definition in Harris (1978),
  which instead centers the index symmetrically around zero.

  Let \\N = n - 1\\ be the last index of the block, and let the running
  index be \\i = 0, 1, \dots, N\\.

  The window weights \\w\[i\]\\ are defined piecewise as follows:

  1.  Left taper region (\\0 \le i \< \alpha N / 2\\): \\ w\[i\] =
      \frac{1}{2} \left\[ 1 + \cos\left( \pi \left( \frac{2 i}{\alpha
      N} - 1 \right) \right) \right\] \\

  2.  Central (non-tapered) region (\\\alpha N / 2 \le i \le N (1 -
      \alpha/2)\\): \\ w\[i\] = 1 \\

  3.  Right taper region (\\N (1 - \alpha/2) \< i \le N\\): \\ w\[i\] =
      \frac{1}{2} \left\[ 1 + \cos\left( \pi \left( \frac{2 i}{\alpha
      N} - \frac{2}{\alpha} + 1 \right) \right) \right\] \\

(Harris, 1978)

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
x <- arima.sim(n = 100, list(ar = 0.8))
result <- tsbs(
  x = as.matrix(x),
  block_length = 10,
  bs_type = "stationary",
  num_blocks = 5,
  num_boots = 10,
  func = mean,
  apply_func_to = "cols"
)
print(result$func_out_means)

# ---- HMM Bootstrap Examples ----

# Generate regime-switching data
set.seed(42)
y <- c(rnorm(150, 0.02, 0.01), rnorm(150, -0.01, 0.03))

# Basic Gaussian HMM (original behavior)
result_gaussian <- tsbs(
  x = y,
  bs_type = "hmm",
  num_states = 2,
  num_boots = 100
)

# MSGARCH with skew Student-t distribution
result_sstd <- tsbs(
  x = y,
  bs_type = "hmm",
  num_states = 2,
  num_boots = 100,
  distribution = "sstd",
  variance_model = "sGARCH",
  seed = 123
)

# Raw HMM without GARCH
result_raw <- tsbs(
  x = y,
  bs_type = "hmm",
  num_states = 2,
  num_boots = 100,
  distribution = "norm_raw",
  seed = 123
)

# Multivariate with market-based regime identification
Y <- matrix(rnorm(600), ncol = 3)
result_multi <- tsbs(
  x = Y,
  bs_type = "hmm",
  num_states = 2,
  num_boots = 50,
  distribution = "norm",
  regime_basis = "market",
  seed = 123
)

# With micro-block sampling to preserve local dependence
result_blocks <- tsbs(
  x = y,
  bs_type = "hmm",
  num_states = 2,
  num_boots = 100,
  distribution = "norm",
  micro_block_length = 5,
  seed = 123
)
} # }
```
