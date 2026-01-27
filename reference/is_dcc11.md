# Check if DCC Order is (1,1)

Determines if the DCC model is of order (1,1), which allows for
analytical gradient computation and reparameterized optimization.

## Usage

``` r
is_dcc11(dcc_pars)
```

## Arguments

- dcc_pars:

  Named list of DCC parameters

## Value

Logical: TRUE if DCC(1,1), FALSE otherwise
