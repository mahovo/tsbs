# Hidden Markov Model (HMM) Bootstrap for Multivariate Time Series

Fits a Hidden Markov Model (HMM) to a multivariate time series and
generates bootstrap replicates by resampling regime-specific blocks.
Supports both Gaussian emissions (via depmixS4) and non-Gaussian
emissions including skew Student-t distributions (via MSGARCH or
standalone EM).

## Usage

``` r
hmm_bootstrap(
  x,
  n_boot = NULL,
  num_states = 2,
  num_blocks = NULL,
  num_boots = 100,
  distribution = c("gaussian", "sstd", "std", "snorm", "sged", "norm", "ged", "sstd_raw",
    "std_raw", "norm_raw"),
  variance_model = c("sGARCH", "eGARCH", "gjrGARCH", "tGARCH"),
  micro_block_length = 1L,
  regime_basis = "market",
  parallel = FALSE,
  num_cores = 1L,
  return_fit = FALSE,
  collect_diagnostics = FALSE,
  verbose = FALSE,
  seed = NULL,
  ...
)
```

## Arguments

- x:

  Numeric vector or matrix representing the time series.

- n_boot:

  Integer, length of each bootstrap series. If NULL, defaults to the
  length of the original series.

- num_states:

  Integer, number of hidden states for the HMM. Default is 2.

- num_blocks:

  Integer, number of blocks to sample for each bootstrap replicate. Only
  used when `distribution = "gaussian"`.

- num_boots:

  Integer, number of bootstrap replicates to generate. Default is 100.

- distribution:

  Character, emission distribution for the HMM. One of:

  "gaussian"

  :   (Default) Gaussian emissions via depmixS4. This is the original
      behavior.

  "sstd"

  :   Skew Student-t via MSGARCH (with GARCH dynamics).

  "std"

  :   Student-t via MSGARCH (with GARCH dynamics).

  "snorm"

  :   Skew normal via MSGARCH (with GARCH dynamics).

  "sged"

  :   Skew GED via MSGARCH (with GARCH dynamics).

  "norm"

  :   Normal via MSGARCH (with GARCH dynamics).

  "ged"

  :   GED via MSGARCH (with GARCH dynamics).

  "sstd_raw"

  :   Skew Student-t without GARCH (standalone HMM-EM).

  "std_raw"

  :   Student-t without GARCH (standalone HMM-EM).

  "norm_raw"

  :   Normal without GARCH (standalone HMM-EM, equivalent to depmixS4
      but using our EM implementation).

- variance_model:

  Character, GARCH specification when using MSGARCH-based distributions.
  One of `"sGARCH"` (default), `"eGARCH"`, `"gjrGARCH"`, or `"tGARCH"`.
  Ignored for Gaussian and raw distributions.

- micro_block_length:

  Integer, block length for within-state sampling when using
  MSGARCH-based bootstrap. Use 1 (default) for iid sampling within
  states, or \>1 to preserve some local dependence. Ignored for Gaussian
  distribution which uses the traditional block sampling.

- regime_basis:

  For multivariate data with MSGARCH: how to identify regimes. One of
  `"market"` (equal-weighted average, default), `"first_pc"` (first
  principal component), or an integer column index. Ignored for Gaussian
  distribution.

- parallel:

  Logical, parallelize computation? Default is FALSE.

- num_cores:

  Integer, number of cores for parallel processing.

- return_fit:

  Logical. If TRUE, returns the fitted model along with bootstrap
  samples. Default is FALSE.

- collect_diagnostics:

  Logical. If TRUE, collects detailed diagnostic information including
  regime composition. Default is FALSE.

- verbose:

  Logical. If TRUE, prints fitting information. Default is FALSE.

- seed:

  Integer, random seed for reproducibility. Default is NULL.

- ...:

  Additional arguments passed to the underlying fitting functions.

## Value

If `return_fit = FALSE` and `collect_diagnostics = FALSE`: A list of
bootstrap replicate matrices.

If `return_fit = TRUE` or `collect_diagnostics = TRUE`: A list
containing:

- bootstrap_series:

  List of bootstrap replicate matrices

- fit:

  (if return_fit = TRUE) The fitted model object

- states:

  The Viterbi state sequence for the original data

- smoothed_probabilities:

  State probabilities matrix

- diagnostics:

  (if collect_diagnostics = TRUE) A tsbs_diagnostics object

- method:

  Character indicating which method was used

## Details

### Method Selection

When `distribution = "gaussian"`, uses depmixS4 for HMM fitting with
Gaussian emissions. This is the original behavior and uses block
resampling within identified regimes.

When `distribution` is one of the MSGARCH distributions ("sstd", "std",
etc.), uses the MSGARCH package for Markov-switching GARCH estimation.
The bootstrap is then semi-parametric: state sequences are simulated
from the fitted Markov chain, while observations are resampled from
state-specific empirical pools.

When `distribution` ends with "\_raw" (e.g., "sstd_raw"), fits an HMM
directly to the returns without the GARCH volatility layer, using a
standalone EM algorithm with the specified emission distribution.

### Literature Background

The MSGARCH-based implementation combines several established
techniques:

- **Markov-switching GARCH**: Haas, Mittnik & Paolella (2004),
  implemented via the MSGARCH package (Ardia et al., 2019)

- **Skew Student-t distribution**: Fernández & Steel (1998)
  transformation

- **Semi-parametric bootstrap**: State sequences are simulated from the
  fitted Markov chain (parametric), while innovations are resampled from
  empirical state-specific pools (nonparametric)

### Multivariate Handling

For multivariate data:

- With Gaussian: fits a joint HMM where all variables depend on the same
  hidden state sequence.

- With MSGARCH: fits the regime model to an aggregate series (see
  `regime_basis`), then samples full cross-sections synchronously to
  preserve empirical dependence structure.

## References

Ardia, D., Bluteau, K., Boudt, K., Catania, L., & Trottier, D.-A.
(2019). Markov-Switching GARCH Models in R: The MSGARCH Package. Journal
of Statistical Software, 91(4), 1-38. doi:10.18637/jss.v091.i04

Fernández, C., & Steel, M. F. (1998). On Bayesian modeling of fat tails
and skewness. Journal of the American Statistical Association, 93(441),
359-371.

Haas, M., Mittnik, S., & Paolella, M. S. (2004). A New Approach to
Markov- Switching GARCH Models. Journal of Financial Econometrics, 2,
493-530.

Hamilton, J. D. (1989). A New Approach to the Economic Analysis of
Nonstationary Time Series and the Business Cycle. Econometrica, 57(2),
357-384.

Holst, U., Lindgren, G., Holst, J. and Thuvesholmen, M. (1994),
Recursive Estimation In Switching Autoregressions With A Markov Regime.
Journal of Time Series Analysis, 15: 489-506.

## See also

[`msvar_bootstrap`](https://mahovo.github.io/tsbs/reference/msvar_bootstrap.md),
[`ms_varma_garch_bs`](https://mahovo.github.io/tsbs/reference/ms_varma_garch_bs.md)

## Examples

``` r
# \donttest{
## Example 1: Traditional Gaussian HMM bootstrap (original behavior)
set.seed(123)
x <- matrix(rnorm(500), ncol = 2)
boot_gaussian <- hmm_bootstrap(x, n_boot = 400, num_states = 2, num_boots = 50)
#> converged at iteration 184 with logLik: -687.2476 

## Example 2: Skew Student-t with MSGARCH (requires MSGARCH package)
if (requireNamespace("MSGARCH", quietly = TRUE)) {
  # Univariate
  y <- rnorm(300)
  boot_sstd <- hmm_bootstrap(y, num_states = 2, num_boots = 50,
                              distribution = "sstd")

  # With diagnostics
  result <- hmm_bootstrap(y, num_states = 2, num_boots = 20,
                           distribution = "sstd",
                           return_fit = TRUE,
                           collect_diagnostics = TRUE)
}
#> Error: Viterbi decoding assigned no observations to state(s): 2. The model effectively has fewer states than specified. Consider reducing 'num_states' or using different data.

## Example 3: Raw returns HMM without GARCH (requires fGarch)
if (requireNamespace("fGarch", quietly = TRUE)) {
  boot_raw <- hmm_bootstrap(y, num_states = 2, num_boots = 50,
                             distribution = "sstd_raw")
}
# }
```
