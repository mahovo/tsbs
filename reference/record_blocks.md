# Record Block Information for a Bootstrap Replicate

Record Block Information for a Bootstrap Replicate

## Usage

``` r
record_blocks(
  diagnostics,
  replicate_idx,
  block_lengths,
  start_positions,
  block_type = NA_character_
)
```

## Arguments

- diagnostics:

  A `tsbs_diagnostics` object.

- replicate_idx:

  Integer, which bootstrap replicate (1-indexed).

- block_lengths:

  Integer vector of block lengths used.

- start_positions:

  Integer vector of starting positions in original series.

- block_type:

  Character, type of block (e.g., "overlapping").

## Value

Updated diagnostics object.
