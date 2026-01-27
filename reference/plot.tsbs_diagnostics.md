# Plot Method for Bootstrap Diagnostics

Creates diagnostic visualizations for bootstrap analysis.

## Usage

``` r
# S3 method for class 'tsbs_diagnostics'
plot(
  x,
  type = c("all", "block_lengths", "start_positions", "means_comparison",
    "acf_comparison", "length_distribution"),
  ...
)
```

## Arguments

- x:

  A `tsbs_diagnostics` object.

- type:

  Character string specifying plot type. One of:

  `"all"`

  :   Produce all applicable plots (default).

  `"block_lengths"`

  :   Histogram of block lengths.

  `"start_positions"`

  :   Histogram of block starting positions.

  `"means_comparison"`

  :   Compare original vs bootstrap means.

  `"acf_comparison"`

  :   Compare original vs bootstrap autocorrelation.

  `"length_distribution"`

  :   Distribution of bootstrap series lengths.

- ...:

  Additional arguments (unused).

## Value

Invisibly returns list of ggplot objects if ggplot2 available, otherwise
uses base R graphics.
