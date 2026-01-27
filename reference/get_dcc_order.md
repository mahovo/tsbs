# Extract DCC order (p, q) from parameter names

Determines (p, q) order from parameter naming convention.

## Usage

``` r
get_dcc_order(dcc_params)
```

## Arguments

- dcc_params:

  Named list of DCC parameters (e.g., alpha_1, alpha_2, beta_1)

## Value

Named vector c(p = ..., q = ...) where q is ARCH order (number of alpha
parameters) and p is GARCH order (number of beta parameters)
