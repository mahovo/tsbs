# Compute GOGARCH Log-Likelihood for MS Framework

Computes the GOGARCH log-likelihood given estimated GARCH parameters and
ICA transformation matrices. This function is used in the E-step of the
EM algorithm for Markov-Switching GOGARCH models.

## Usage

``` r
compute_gogarch_loglik_ms(
  residuals,
  garch_pars,
  ica_info,
  distribution = "norm",
  return_vector = FALSE
)
```

## Arguments

- residuals:

  Numeric matrix of residuals with dimensions T x k.

- garch_pars:

  List of GARCH parameters for each component.

- ica_info:

  List containing ICA transformation matrices (A, W, K).

- distribution:

  Character string specifying the component distribution.

- return_vector:

  Logical; if `TRUE`, return per-observation log-likelihoods. Default is
  `FALSE`.

## Value

Scalar total log-likelihood, or vector of per-observation values if
`return_vector = TRUE`.

## Details

The GOGARCH log-likelihood consists of two parts:

**1. Component Log-Likelihoods**

For each independent component \\i\\ and time \\t\\: \$\$LL\_{i,t} =
\log f(s\_{i,t} \| \sigma\_{i,t})\$\$

**2. Jacobian Adjustment**

The ICA transformation introduces a Jacobian term: \$\$LL\_{jacobian} =
\log \|det(K)\|\$\$ where \\K\\ is the pre-whitening matrix from ICA.

**Total Log-Likelihood** \$\$LL = \sum_t \sum_i LL\_{i,t} + \log
\|det(K)\|\$\$

## See also

- [`estimate_garch_weighted_gogarch`](https://mahovo.github.io/tsbs/reference/estimate_garch_weighted_gogarch.md):
  Estimation function

- [`estimate_garch_weighted_univariate_gogarch`](https://mahovo.github.io/tsbs/reference/estimate_garch_weighted_univariate_gogarch.md):
  Component estimation
