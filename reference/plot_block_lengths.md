# Plot Block Length Distribution

Visualizes the distribution of block lengths used in bootstrap
replicates. Particularly useful for Stationary Bootstrap where lengths
are random.

## Usage

``` r
plot_block_lengths(diagnostics, show_expected = TRUE)
```

## Arguments

- diagnostics:

  A tsbs_diagnostics object from block bootstrap.

- show_expected:

  Logical. If TRUE and bs_type is "stationary", show the expected
  geometric distribution. Default TRUE.

## Value

A ggplot object (invisibly).
