# Sample from State-Specific Pools with Synchronized Sampling

For each position in a simulated state sequence, samples a time index
from the original data where that state occurred. Supports both
individual (iid) sampling and micro-block sampling within states.

## Usage

``` r
sample_from_pools_sync(sim_states, pools, micro_block_length = 1L)
```

## Arguments

- sim_states:

  Integer vector of simulated state sequence (length T).

- pools:

  State-specific pools from `extract_state_pools`.

- micro_block_length:

  Integer, length of micro-blocks for within-state sampling. Use 1 for
  iid sampling (default), \>1 to preserve local dependence.

## Value

Integer vector of length T containing sampled time indices from the
original data.

## Details

When `micro_block_length = 1`, each time point is sampled independently
from the appropriate state pool.

When `micro_block_length > 1`, the function identifies "runs" of
consecutive time points in the same state and fills them with
micro-blocks sampled from that state's pool. This preserves some local
autocorrelation in innovations while still allowing the Markov chain to
drive regime dynamics.

For multivariate data, the same time indices are used across all
variables, which preserves cross-sectional dependence.
