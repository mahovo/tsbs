# Value not in unit interval

Returns TRUE if `p` is not in the unit interval, FALSE otherwise.

## Usage

``` r
.not_in_unit_interval(
  p,
  allow_null = TRUE,
  interval_type = c("open", "closed"),
  fail_mode = c("predictably", "gracefully")
)
```

## Arguments

- p:

  Value to validate, e.g. a percentage (decimal fraction) .

- allow_null:

  boolean, allow NULL?

- interval_type:

  One of `"open"` or `"closed"`. Indicates if the valid interval is an
  open or closed unit interval, i.e. \\p \in (0, 1)\\ or \\p \in \[0,
  1\]\\ respectively.

- fail_mode:

  How to handle validation errors: "predictably" (fail fast) or
  "gracefully" (return FALSE on error)

## Value

boolean

## Details

If
