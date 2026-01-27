# Stationary Bootstrap for a 2-State MS-VAR(1) Model

This function first fits a 2-state Markov-Switching Vector
Autoregressive (MS-VAR) model of order 1 to the provided multivariate
time series data. It then uses the estimated state sequence to perform a
stationary bootstrap, generating resampled time series that preserve the
state-dependent properties of the original data.

## Usage

``` r
msvar_bootstrap(
  x,
  n_boot = NULL,
  num_blocks = NULL,
  num_boots = 100,
  parallel = FALSE,
  num_cores = 1L,
  return_fit = FALSE,
  collect_diagnostics = FALSE,
  verbose = TRUE
)
```

## Arguments

- x:

  A numeric matrix or data frame where rows are observations and columns
  are the time series variables.

- n_boot:

  An integer specifying the length of each bootstrapped series. If NULL
  (the default), the length of the original series is used.

- num_blocks:

  An integer specifying the number of blocks to sample for the
  bootstrap. Defaults to 100.

- num_boots:

  An integer specifying the total number of bootstrap samples to
  generate. Defaults to 100.

- parallel:

  A logical value indicating whether to use parallel processing for
  generating bootstrap samples. Defaults to FALSE.

- num_cores:

  An integer specifying the number of cores to use for parallel
  processing. Only used if `parallel` is TRUE. Defaults to 1.

- return_fit:

  Logical. If TRUE, returns the fitted MS-VAR model along with bootstrap
  samples. Default is FALSE.

- collect_diagnostics:

  Logical. If TRUE, collects detailed diagnostic information including
  regime composition of each bootstrap replicate. Default is FALSE.

- verbose:

  Logical. If TRUE, prints fitting information. Default is TRUE.

## Value

If `return_fit = FALSE` and `collect_diagnostics = FALSE`: A list of
bootstrap replicate matrices.

Otherwise, a list containing:

- bootstrap_series:

  List of bootstrap replicate matrices

- fit:

  (if return_fit = TRUE) The fitted MS-VAR model object

- states:

  The state sequence for the original data

- smoothed_probabilities:

  (if return_fit = TRUE) Matrix of smoothed state probabilities

- diagnostics:

  (if collect_diagnostics = TRUE) A tsbs_diagnostics object

## Details

For a stationary bootstrap based on a more general \\n\\-state MS-VECTOR
ARIMA(\\p, d, q\\)-GARCH model see
[`ms_varma_garch_bs()`](https://mahovo.github.io/tsbs/reference/ms_varma_garch_bs.md).

This function:

- Fits a 2-state MS-VAR(1) model using
  [`fit_msvar()`](https://mahovo.github.io/tsbs/reference/fit_msvar.md)

- Uses smoothed probabilities to determine the most likely state
  sequence

- Samples contiguous blocks of observations within each regime

&nbsp;

- \\y_t \in \mathbb{R}^K\\ be a **\$K\$-dimensional multivariate
  response vector** at time \\t\\

- \\S_t \in {1, \dots, M}\\ be a **latent Markov chain** with \\M\\
  discrete regimes

- \\p\\ be the **lag order** of the VAR model

\\ y_t = \mu^{(S_t)} + \sum\_{i=1}^{p} A_i^{(S_t)} y\_{t-i} +
\varepsilon_t, \quad \varepsilon_t \sim \mathcal{N}(0, \Sigma^{(S_t)})
\\

Where:

- \\\mu^{(S_t)} \in \mathbb{R}^K\\ is the regime-specific intercept
  vector

- \\A_i^{(S_t)} \in \mathbb{R}^{K \times K}\\ are the **regime-specific
  autoregressive coefficient matrices**

- \\\Sigma^{(S_t)} \in \mathbb{R}^{K \times K}\\ is the regime-specific
  error covariance matrix

If `n_boot` is set, the last block will be trimmed when necessary. If
`n_boot` is not set, and `num_blocks` is set, the length of each
bootstrap series will be determined by the number of blocks and the
random lengths of the individual blocks for that particular series. If
neither `n_boot` nor `num_blocks` is set, `n_boot` will default to the
number of rows in `x` and the last block will be trimmed when necessary.

When `collect_diagnostics = TRUE`, the function records:

- Original state sequence from smoothed probabilities

- State sequence for each bootstrap replicate

- Block information (lengths, source positions)

- Source indices for probability reconstruction

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate sample data
set.seed(123)
T_obs <- 250
y1 <- arima.sim(model = list(ar = 0.7), n = T_obs)
y2 <- 0.5 * y1 + arima.sim(model = list(ar = 0.3), n = T_obs)
sample_data <- cbind(y1, y2)

# Basic bootstrap
boot_samples <- msvar_bootstrap(sample_data, num_boots = 50)

# With diagnostics for visualization
result <- msvar_bootstrap(
  sample_data, 
  num_boots = 10,
  return_fit = TRUE,
  collect_diagnostics = TRUE
)
plot_regime_composition(result, sample_data)
} # }
```
