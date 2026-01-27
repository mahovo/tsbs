# Transform reparameterized space back to DCC(1,1) parameters

Transforms (persistence, ratio) back to (alpha, beta). Stationarity
(alpha + beta \< 1) is guaranteed if persistence \< 1.

## Usage

``` r
dcc_from_reparam(persistence, ratio)
```

## Arguments

- persistence:

  Total persistence (alpha + beta), must be in (0, 1)

- ratio:

  Proportion allocated to alpha, must be in (0, 1)

## Value

Named vector c(alpha, beta)
