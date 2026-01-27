# Validate DCC orders in model specification

Validates all DCC orders in a multi-state spec list. Called during model
fitting to ensure compatibility with installed tsmarch.

## Usage

``` r
validate_spec_dcc_orders(spec, action = c("error", "warn"))
```

## Arguments

- spec:

  List of state specifications

- action:

  Character: "error" to stop, "warn" to warn and fall back

## Value

Validated/adjusted spec (modified in place if needed)
