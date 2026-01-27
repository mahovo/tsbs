# Block Bootstrap with Diagnostics

Performs Moving Block Bootstrap or Stationary Bootstrap while collecting
detailed diagnostic information about block composition.

## Usage

``` r
blockBootstrap_with_diagnostics(
  x,
  n_boot = NULL,
  block_length = NULL,
  bs_type = c("moving", "stationary"),
  block_type = c("overlapping", "non-overlapping"),
  num_boots = 100L,
  p = 0.1,
  collect_diagnostics = TRUE,
  taper_type = "cosine",
  tukey_alpha = 0.5
)
```

## Arguments

- x:

  Numeric matrix with rows as time points, columns as variables.

- n_boot:

  Integer. Desired length of bootstrap series. If NULL, uses length of
  original series.

- block_length:

  Integer. Block length for moving bootstrap, or expected block length
  for stationary bootstrap. If NULL, computed automatically.

- bs_type:

  Character. Either "moving" or "stationary".

- block_type:

  Character. One of "overlapping", "non-overlapping", or "tapered".

- num_boots:

  Integer. Number of bootstrap replicates.

- p:

  Numeric in (0,1). Probability parameter for geometric distribution in
  stationary bootstrap. Default 0.1.

- collect_diagnostics:

  Logical. If TRUE, collect block-level diagnostics.

- taper_type:

  Character. Taper window type if block_type = "tapered".

- tukey_alpha:

  Numeric. Alpha parameter for Tukey window.

## Value

If collect_diagnostics = FALSE, a list of bootstrap matrices. If
collect_diagnostics = TRUE, a list with:

- bootstrap_series:

  List of bootstrap matrices

- diagnostics:

  A tsbs_diagnostics object with block information

## Details

This function provides an R implementation of block bootstrap that
tracks diagnostic information. For each bootstrap replicate, it records:

- Block lengths used

- Starting positions in original series

- Source indices mapping each bootstrap observation to original

For the **Moving Block Bootstrap**, blocks have fixed length and are
sampled uniformly from all possible starting positions.

For the **Stationary Bootstrap**, block lengths are drawn from a
geometric distribution with parameter p, giving expected block length
1/p.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
x <- matrix(rnorm(500), ncol = 2)

# Moving block bootstrap with diagnostics
result <- blockBootstrap_with_diagnostics(
  x, block_length = 10, bs_type = "moving",
  num_boots = 50, collect_diagnostics = TRUE
)

# View block length distribution
summary(result$diagnostics)

# Plot coverage
plot_block_coverage(result$diagnostics)
} # }
```
