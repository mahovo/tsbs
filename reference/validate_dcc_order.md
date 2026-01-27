# Validate DCC order against tsmarch version

Checks if requested DCC order is supported by installed tsmarch. Issues
informative error if higher-order DCC requested with buggy tsmarch.

## Usage

``` r
validate_dcc_order(dcc_order, action = c("error", "warn"))
```

## Arguments

- dcc_order:

  Integer vector c(p, q) for DCC(p,q) order

- action:

  Character: "error" to stop, "warn" to warn and fall back to (1,1)

## Value

Validated/adjusted dcc_order
