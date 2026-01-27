# Enhanced Block Sampling with Diagnostic Tracking

Internal function that performs regime-based block sampling while
optionally tracking diagnostic information about each bootstrap
replicate.

## Usage

``` r
.sample_blocks_with_diagnostics(
  x,
  n_boot,
  num_blocks,
  states,
  num_boots,
  parallel = FALSE,
  num_cores = 1L,
  collect_diagnostics = FALSE
)
```

## Arguments

- x:

  Original data matrix.

- n_boot:

  Target bootstrap length.

- num_blocks:

  Number of blocks to sample.

- states:

  State sequence for original data.

- num_boots:

  Number of bootstrap replicates.

- parallel:

  Use parallel processing.

- num_cores:

  Number of cores.

- collect_diagnostics:

  Track diagnostic info.

## Value

List with samples and optional replicate_info.
