# Find Starting Points for Consecutive Sequences

Helper function to find indices in a pool that can serve as starting
points for blocks of a specified length with consecutive values.

## Usage

``` r
.find_consecutive_starts(pool_idx, block_len)
```

## Arguments

- pool_idx:

  Integer vector of available indices.

- block_len:

  Integer, required block length.

## Value

Integer vector of valid starting points.
