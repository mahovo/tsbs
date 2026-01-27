# Helper Function for Multivariate Stationary Block Bootstrap

Helper Function for Multivariate Stationary Block Bootstrap

## Usage

``` r
.sample_blocks(
  x,
  n_boot,
  num_blocks,
  states,
  num_boots,
  parallel = FALSE,
  num_cores = 1L
)
```

## Arguments

- x:

  The original multivariate time series as a matrix.

- n_boot:

  The target length for each bootstrap sample.

- num_blocks:

  The number of blocks to sample.

- states:

  A vector of integer states corresponding to each row of x.

- num_boots:

  The number of bootstrap series to generate.

- parallel:

  Logical, whether to use parallel processing.

- num_cores:

  The number of cores for parallel processing.

## Value

A list of bootstrapped series.
