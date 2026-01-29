# Generate Bootstrap Samples Using MSGARCH

Core bootstrap generation function. Simulates new state sequences from
the fitted Markov chain and samples observations from state-specific
empirical pools (semi-parametric bootstrap).

## Usage

``` r
generate_msgarch_bootstrap(
  fit,
  y,
  state_info,
  pools,
  n_boot = 1000L,
  micro_block_length = 1L,
  sync_sampling = TRUE,
  seed = NULL
)
```

## Arguments

- fit:

  MSGARCH fit object from `fit_msgarch_model`.

- y:

  Original data. Can be a vector (univariate) or matrix (multivariate).
  For multivariate data with synchronized sampling, the regime model is
  fit to an aggregate, and full cross-sections are sampled together.

- state_info:

  Output from `extract_msgarch_states`.

- pools:

  Output from `extract_state_pools`.

- n_boot:

  Integer, number of bootstrap replicates (default 1000).

- micro_block_length:

  Integer, block length for within-state sampling. Use 1 for iid
  sampling (default), \>1 to preserve local dependence.

- sync_sampling:

  Logical. If TRUE and y is multivariate, sample the same time indices
  across all assets to preserve cross-sectional dependence. Default is
  TRUE.

- seed:

  Integer, random seed for reproducibility. Default is NULL.

## Value

A list containing:

- samples:

  3D array of bootstrap samples with dimensions (T x N x n_boot) where T
  is series length and N is number of variables.

- states:

  Matrix of simulated state sequences (T x n_boot).

- sampled_indices:

  Matrix of sampled time indices (T x n_boot).

- state_info:

  The input state_info for reference.

- pools:

  The input pools for reference.

- params:

  List of bootstrap parameters used.

## Details

The bootstrap procedure is semi-parametric:

- **Parametric component**: State sequences are simulated from the
  fitted Markov chain with transition matrix P.

- **Nonparametric component**: Observations are resampled from empirical
  state-specific pools, preserving the actual distributional
  characteristics without parametric assumptions.

For multivariate data, the `sync_sampling = TRUE` option ensures that
the same time index is sampled for all variables, preserving whatever
cross-sectional dependence existed at that moment in the original data.
This is simpler than copula-based approaches and automatically captures
the empirical dependence structure.

## See also

[`fit_msgarch_model`](https://mahovo.github.io/tsbs/reference/fit_msgarch_model.md),
[`extract_msgarch_states`](https://mahovo.github.io/tsbs/reference/extract_msgarch_states.md),
[`extract_state_pools`](https://mahovo.github.io/tsbs/reference/extract_state_pools.md)
