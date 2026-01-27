# Check if tsmarch supports higher-order DCC

tsmarch v1.0.0 has a bug in `.copula_parameters` where
`paste0("beta_",1:order[1])` should be `paste0("beta_",1:order[2])`.
This causes DCC(p,q) with p!=q to create wrong number of parameters.
Fixed in v1.0.1.

## Usage

``` r
tsmarch_supports_higher_order_dcc()
```

## Value

Logical TRUE if tsmarch \>= 1.0.1, FALSE otherwise
