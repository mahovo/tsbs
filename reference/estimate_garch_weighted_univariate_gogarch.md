# Weighted Univariate GARCH Estimation for GOGARCH Components

Estimates univariate GARCH parameters for a single ICA component using
weighted maximum likelihood. This function is called by
[`estimate_garch_weighted_gogarch`](https://mahovo.github.io/tsbs/reference/estimate_garch_weighted_gogarch.md)
for each independent component extracted by ICA.

## Usage

``` r
estimate_garch_weighted_univariate_gogarch(
  residuals,
  weights,
  spec,
  verbose = FALSE
)
```

## Arguments

- residuals:

  Numeric vector of ICA component values with length T.

- weights:

  Numeric vector of weights (state probabilities) with length T.

- spec:

  List containing the component specification:

  `garch_model`

  :   GARCH model type (default: `"garch"`)

  `garch_order`

  :   GARCH order as c(p, q) (default: `c(1, 1)`)

  `distribution`

  :   Component distribution

  `start_pars`

  :   List with `garch_pars` and `dist_pars`

- verbose:

  Logical; if `TRUE`, print progress information.

## Value

A list with:

- `coefficients`:

  Named list of estimated parameters

- `warnings`:

  List of any warnings from optimization

## Details

The function implements weighted MLE for GARCH(p,q) models on ICA
components. The weighted log-likelihood is: \$\$LL_w = \sum_t w_t \cdot
\log f(s_t \| \sigma_t)\$\$ where \\w_t\\ are the state probabilities,
\\s_t\\ is the ICA component value, and \\\sigma_t\\ is the GARCH
conditional volatility.

**GARCH Recursion**

The variance recursion for GARCH(p,q) is: \$\$\sigma^2_t = \omega +
\sum\_{i=1}^{p} \alpha_i s^2\_{t-i} + \sum\_{j=1}^{q} \beta_j
\sigma^2\_{t-j}\$\$

**Supported Distributions**

|                  |                     |                         |
|------------------|---------------------|-------------------------|
| **Distribution** | **Parameters**      | **Description**         |
| `"norm"`         | none                | Normal (Gaussian)       |
| `"std"`          | shape               | Student-t               |
| `"sstd"`         | shape, skew         | Skewed Student-t        |
| `"nig"`          | shape, skew         | Normal Inverse Gaussian |
| `"gh"`           | shape, skew, lambda | Generalized Hyperbolic  |

## See also

- [`estimate_garch_weighted_gogarch`](https://mahovo.github.io/tsbs/reference/estimate_garch_weighted_gogarch.md):
  Main GOGARCH estimator

- [`compute_gogarch_loglik_ms`](https://mahovo.github.io/tsbs/reference/compute_gogarch_loglik_ms.md):
  Log-likelihood computation
