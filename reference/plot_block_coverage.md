# Plot Block Coverage for Bootstrap Diagnostics

Visualizes which parts of the original series are sampled across
bootstrap replicates, showing coverage patterns and potential gaps.

## Usage

``` r
plot_block_coverage(
  diagnostics,
  type = c("heatmap", "histogram", "blocks"),
  max_replicates = 20
)
```

## Arguments

- diagnostics:

  A tsbs_diagnostics object from block bootstrap.

- type:

  Character. Type of plot:

  - "heatmap": Heatmap showing coverage intensity at each time point

  - "histogram": Histogram of block starting positions

  - "blocks": Visual representation of blocks sampled per replicate

- max_replicates:

  Integer. Maximum number of replicates to show in "blocks" plot.
  Default 20.

## Value

A ggplot object (invisibly).

## Examples

``` r
if (FALSE) { # \dontrun{
result <- blockBootstrap_with_diagnostics(x, bs_type = "moving",
                                           collect_diagnostics = TRUE)
plot_block_coverage(result$diagnostics, type = "heatmap")
plot_block_coverage(result$diagnostics, type = "histogram")
} # }
```
