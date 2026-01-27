# Estimate Block Lengths from Bootstrap Series

For block bootstrap methods, attempts to identify block boundaries by
detecting discontinuities in the bootstrap series relative to the
original.

## Usage

``` r
estimate_block_lengths_from_series(
  diagnostics,
  bootstrap_series,
  original_data
)
```

## Arguments

- diagnostics:

  A `tsbs_diagnostics` object.

- bootstrap_series:

  List of bootstrap replicate matrices.

- original_data:

  Original data matrix.

## Value

Updated diagnostics object with estimated block information.
