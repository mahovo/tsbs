# Transform DCC(1,1) parameters to reparameterized space

Transforms (alpha, beta) with constraint alpha + beta \< 1 to
(persistence, ratio) where both are in (0, 1) with no joint constraint.

## Usage

``` r
dcc_to_reparam(alpha, beta)
```

## Arguments

- alpha:

  DCC alpha parameter

- beta:

  DCC beta parameter

## Value

Named vector c(persistence, ratio)
