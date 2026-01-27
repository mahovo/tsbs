# Invalid input data

DEPRECATED

## Usage

``` r
.is_invalid_data(
  x,
  allow_null = TRUE,
  fail_mode = c("predictably", "gracefully")
)
```

## Arguments

- x:

  x, vector, matrix, data.frame or time series

- allow_null:

  NULL value may be allowed when functions calculates value of x
  automatically.

- fail_mode:

  How to handle validation errors: "predictably" (fail fast) or
  "gracefully" (return FALSE on error)

## Value

boolean

## Details

Returns TRUE if input data is not valid, FALSE otherwise.
