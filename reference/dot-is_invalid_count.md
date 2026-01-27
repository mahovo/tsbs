# Invalid input count

Returns TRUE if input count is not valid, FALSE otherwise.

## Usage

``` r
.is_invalid_count(
  n,
  allow_null = TRUE,
  fail_mode = c("predictably", "gracefully")
)
```

## Arguments

- n:

  Count value to validate

- allow_null:

  NULL value may be allowed when functions calculates value of n
  automatically.

- fail_mode:

  How to handle validation errors: "predictably" (fail fast) or
  "gracefully" (return FALSE on error)

## Value

boolean
