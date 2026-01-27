# Extract Bootstrap Output Matrix

Converts the list of functional outputs from
[`tsbs()`](https://mahovo.github.io/tsbs/reference/tsbs.md) into a
matrix format suitable for further analysis.

## Usage

``` r
extract_func_out_matrix(func_outs)
```

## Arguments

- func_outs:

  List of functional outputs from
  [`tsbs()`](https://mahovo.github.io/tsbs/reference/tsbs.md).

## Value

A matrix where each row is a bootstrap replicate and each column is an
output dimension. Returns NULL if conversion fails.
