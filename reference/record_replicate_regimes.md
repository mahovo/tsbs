# Record Regime Composition for a Bootstrap Replicate

Records the state sequence and block-state mapping for a single
bootstrap replicate.

## Usage

``` r
record_replicate_regimes(
  diagnostics,
  replicate_idx,
  replicate_states,
  block_states = NULL,
  block_source_indices = NULL,
  source_indices = NULL
)
```

## Arguments

- diagnostics:

  A `tsbs_diagnostics` object.

- replicate_idx:

  Integer, which bootstrap replicate (1-indexed).

- replicate_states:

  Integer vector of states for this replicate.

- block_states:

  Optional integer vector of states for each block used.

- block_source_indices:

  Optional integer vector of which original block each bootstrap block
  came from.

- source_indices:

  Optional integer vector mapping each bootstrap time point to its
  source index in the original series. Used for probability lookup.

## Value

Updated diagnostics object.
