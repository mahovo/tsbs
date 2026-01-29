# Extract State-Specific Innovation Pools

Builds pools of time indices for each regime based on Viterbi decoding.
These pools are used for synchronized sampling in the bootstrap
procedure.

## Usage

``` r
extract_state_pools(fit, y, state_info)
```

## Arguments

- fit:

  An MSGARCH fit object.

- y:

  Numeric vector, the original data.

- state_info:

  Output from `extract_msgarch_states`.

## Value

A list of length `n_states`, where each element contains:

- indices:

  Integer vector of time indices assigned to this state.

- n:

  Integer, count of observations in this state.

- proportion:

  Numeric, proportion of total observations.

## Details

For multivariate bootstrap with synchronized sampling, we sample entire
time indices (not individual residuals) from each state pool. This
preserves the cross-sectional dependence structure that existed at each
time point.
