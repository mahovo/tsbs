# Inverse Reparameterization for DCC(1,1)

Transform unconstrained (psi, phi) back to constrained (alpha, beta).
Stationarity is guaranteed for any finite (psi, phi).

## Usage

``` r
dcc11_from_unconstrained(psi, phi)
```

## Arguments

- psi:

  Unconstrained persistence parameter (can be any real number)

- phi:

  Unconstrained ratio parameter (can be any real number)

## Value

Named vector c(alpha, beta)
