# Compute Diagnostics from Bootstrap Series

After bootstrap series have been generated, this function computes
diagnostic statistics from the actual series. Use this when block-level
tracking is not available (e.g., when using C++ backend).

## Usage

``` r
compute_bootstrap_diagnostics(
  bootstrap_series,
  original_data,
  bs_type,
  config = list()
)
```

## Arguments

- bootstrap_series:

  List of bootstrap replicate matrices.

- original_data:

  Original data matrix.

- bs_type:

  Character string specifying bootstrap type.

- config:

  Named list of configuration parameters used.

## Value

A `tsbs_diagnostics` object.
