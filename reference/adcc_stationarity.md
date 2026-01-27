# ADCC Stationarity Constraint

Computes the ADCC stationarity constraint value. For ADCC: α + β + δ*γ
\< 1, where δ depends on the data. A simplified constraint is: α + β +
0.5*γ \< 1

## Usage

``` r
adcc_stationarity(alpha, gamma, beta, delta = 0.5)
```

## Arguments

- alpha:

  ADCC alpha

- gamma:

  ADCC gamma

- beta:

  ADCC beta

- delta:

  Asymmetry scaling factor (default 0.5)

## Value

Stationarity measure (should be \< 1 for stationarity)
