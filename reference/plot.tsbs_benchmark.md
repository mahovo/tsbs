# Plot Method for tsbs_benchmark

Creates a visualization of benchmark results showing how runtime scales
with the varying parameter across different bootstrap setups.

## Usage

``` r
# S3 method for class 'tsbs_benchmark'
plot(x, log_scale = FALSE, show_points = TRUE, show_ribbon = TRUE, ...)
```

## Arguments

- x:

  A `tsbs_benchmark` object.

- log_scale:

  Logical. If TRUE, use log scale for y-axis. Default is FALSE.

- show_points:

  Logical. If TRUE, show individual timing points. Default is TRUE.

- show_ribbon:

  Logical. If TRUE, show +/- 1 SD ribbon. Default is TRUE.

- ...:

  Additional arguments passed to ggplot.

## Value

A ggplot object (invisibly).
