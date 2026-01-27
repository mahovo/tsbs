# Reparameterize DCC(1,1) to Unconstrained Space

Transform constrained (alpha, beta) to unconstrained (psi, phi). psi =
logit(alpha + beta) controls persistence phi = log(alpha / beta)
controls the ratio

## Usage

``` r
dcc11_to_unconstrained(alpha, beta)
```

## Arguments

- alpha:

  DCC alpha parameter (0 \< alpha)

- beta:

  DCC beta parameter (0 \< beta)

## Value

Named vector c(psi, phi)
