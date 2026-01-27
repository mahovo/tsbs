# Coerce, Validate, and Prepare Time Series Data

A robust helper function that takes a user-provided object, validates
it, and coerces it into a clean, numeric matrix suitable for downstream
analysis. It provides informative error messages or fails gracefully
based on the specified mode.

## Usage

``` r
.coerce_and_validate_data(x, fail_mode = c("predictably", "gracefully"))
```

## Arguments

- x:

  The user-provided data (e.g., vector, matrix, data.frame, ts, xts,
  zoo).

- fail_mode:

  How to handle validation errors: "predictably" (the default) will stop
  with an informative error. "gracefully" will return NULL.

## Value

A numeric matrix if the input is valid, otherwise stops or returns NULL.
