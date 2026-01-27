# Weighted GARCH Estimation for GOGARCH Models

Implements weighted maximum likelihood estimation for Generalized
Orthogonal GARCH (GOGARCH) models in the context of Markov-Switching
frameworks.

GOGARCH differs fundamentally from DCC and Copula GARCH in how it models
dynamic correlations:

- **DCC/CGARCH**: Estimate explicit correlation dynamics parameters
  (alpha, beta) that govern how correlations evolve over time.

- **GOGARCH**: Assume observations arise from independent latent
  factors. Time-varying correlations emerge from the time-varying
  volatilities of these factors, transformed through a fixed mixing
  matrix.

## Usage

``` r
estimate_garch_weighted_gogarch(
  residuals,
  weights,
  spec,
  diagnostics = NULL,
  verbose = FALSE
)
```

## Arguments

- residuals:

  Numeric matrix of residuals with dimensions T x k, where T is the
  number of observations and k is the number of series.

- weights:

  Numeric vector of state probabilities/weights with length T. These
  weights come from the E-step of the EM algorithm in the MS-VARMA-GARCH
  framework.

- spec:

  List containing the model specification with the following elements:

  `garch_spec_args`

  :   List with:

      - `model`: GARCH model type (default: `"garch"`). Valid choices
        are “garch” for vanilla GARCH, “gjrgarch” for asymmetric GARCH,
        “egarch” for exponential GARCH, “aparch” for asymmetric power
        ARCH, “csGARCH” for the component GARCH, “igarch” for the
        integrated GARCH.

      - `order`: GARCH order as c(p, q) (default: `c(1, 1)`)

      - `ica`: ICA algorithm. Currently only `"radical"` is supported
        (via tsmarch v1.0.0)

      - `components`: Number of ICA components to extract

  `distribution`

  :   Component distribution: `"norm"`, `"nig"`, or `"gh"`

  `start_pars`

  :   List with:

      - `garch_pars`: List of starting GARCH parameters for each
        component (e.g.,
        `list(list(omega=0.1, alpha1=0.1, beta1=0.8), ...)`)

      - `dist_pars`: Distribution parameters (e.g.,
        `list(shape=1, skew=0)` for NIG)

- diagnostics:

  Optional diagnostics collector object for logging estimation progress
  and issues.

- verbose:

  Logical; if `TRUE`, print progress information during ICA
  decomposition and GARCH estimation. Default is `FALSE`.

## Value

A list with the following components:

- `coefficients`:

  List containing:

  - `garch_pars`: List of estimated GARCH parameters for each ICA
    component

  - `ica_info`: List with ICA results:

    - `A`: Mixing matrix (k x n_components)

    - `W`: Unmixing matrix (n_components x k)

    - `K`: Pre-whitening matrix (for Jacobian)

    - `S`: Independent components matrix (T x n_components)

    - `method`: ICA algorithm used

    - `n_components`: Number of components extracted

  - `dist_pars`: Distribution parameters

  - `correlation_type`: Always `"gogarch"`

- `warnings`:

  List of warning messages generated during estimation

- `diagnostics`:

  Updated diagnostics object

## Details

The estimation proceeds in two stages:

**Stage 1: ICA Decomposition**

Independent Component Analysis extracts statistically independent
components from the multivariate residuals: \$\$S = Y \cdot W\$\$ where
\\Y\\ is the T x k residual matrix, \\W\\ is the unmixing matrix, and
\\S\\ contains the independent components.

The RADICAL algorithm (Learned-Miller, 2003) is used for ICA
decomposition, which is robust to outliers. For improved convergence
with multiple restarts and quality diagnostics, see
[`improved_ica_decomposition`](https://mahovo.github.io/tsbs/reference/improved_ica_decomposition.md).

**Stage 2: Component GARCH Estimation**

Univariate GARCH models are fitted to each independent component using
weighted MLE. The weights come from the state probabilities in the
Markov-Switching framework.

**Covariance Reconstruction**

The time-varying covariance matrix is reconstructed as: \$\$H_t = A
\cdot D_t \cdot A'\$\$ where \\A\\ is the mixing matrix from ICA and
\\D_t = \text{diag}(\sigma^2\_{1,t}, \ldots, \sigma^2\_{k,t})\\ contains
the component GARCH variances.

**Log-Likelihood**

The GOGARCH log-likelihood includes a Jacobian adjustment for the ICA
transformation: \$\$LL = \sum_i LL\_{component,i} + \log\|det(K)\|\$\$
where \\K\\ is the pre-whitening matrix. See
[`compute_gogarch_loglik_ms`](https://mahovo.github.io/tsbs/reference/compute_gogarch_loglik_ms.md).

## Comparison with DCC and CGARCH

|                    |                    |                 |                         |
|--------------------|--------------------|-----------------|-------------------------|
| **Aspect**         | **DCC**            | **CGARCH**      | **GOGARCH**             |
| Correlation source | alpha, beta params | DCC on copula   | ICA mixing matrix       |
| Dynamics           | DCC recursion      | DCC/ADCC        | Component volatilities  |
| Marginal treatment | Assumed Normal     | Flexible (PIT)  | ICA components          |
| Key strength       | Interpretable      | Tail dependence | Non-Gaussian dependence |

## When to Use GOGARCH

GOGARCH is particularly suitable when:

- The dependence structure arises from independent underlying factors

- Heavy-tailed distributions (NIG, GH) are needed for the components

- Dimension reduction is desired (fewer components than series)

- Non-linear dependence structures are suspected

## ICA Quality Assessment

The quality of the ICA decomposition is critical for GOGARCH
reliability. Key diagnostics include:

- **Independence score**: Measures how uncorrelated the extracted
  components are (target: \> 0.8)

- **Reconstruction error**: How well Y = S \* A' holds (target: \< 1\\

- **Negentropy**: Non-Gaussianity of components (higher = better
  separation

If ICA convergence is problematic, consider using
[`improved_ica_decomposition`](https://mahovo.github.io/tsbs/reference/improved_ica_decomposition.md)
which provides multiple restarts and automatic quality assessment. For
post-estimation diagnostics, see
[`gogarch_diagnostics`](https://mahovo.github.io/tsbs/reference/gogarch_diagnostics.md).

## References

van der Weide, R. (2002). GO-GARCH: A multivariate generalized
orthogonal GARCH model. *Journal of Applied Econometrics*, 17(5),
549-564.

Learned-Miller, E. G. (2003). ICA using spacings estimates of entropy.
*Journal of Machine Learning Research*, 4, 1271-1295.

## See also

- [`tsbs`](https://mahovo.github.io/tsbs/reference/tsbs.md): Main
  bootstrap function (use `garch_spec_fun = "gogarch_modelspec"`)

- [`estimate_garch_weighted_univariate_gogarch`](https://mahovo.github.io/tsbs/reference/estimate_garch_weighted_univariate_gogarch.md):
  Component GARCH estimation

- [`compute_gogarch_loglik_ms`](https://mahovo.github.io/tsbs/reference/compute_gogarch_loglik_ms.md):
  Log-likelihood computation

- [`improved_ica_decomposition`](https://mahovo.github.io/tsbs/reference/improved_ica_decomposition.md):
  Enhanced ICA with multiple restarts and quality metrics

- [`gogarch_diagnostics`](https://mahovo.github.io/tsbs/reference/gogarch_diagnostics.md):
  Comprehensive model diagnostics including ICA quality, component GARCH
  fit, and covariance reconstruction

- [`estimate_garch_weighted_dcc`](https://mahovo.github.io/tsbs/reference/estimate_garch_weighted_dcc.md):
  Alternative DCC estimator

- [`estimate_garch_weighted_cgarch`](https://mahovo.github.io/tsbs/reference/estimate_garch_weighted_cgarch.md):
  Alternative CGARCH estimator
