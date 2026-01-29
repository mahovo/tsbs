# Get Series for Regime Identification (Multivariate Case)

For multivariate data, determines which series (or aggregate) to use for
fitting the regime model.

## Usage

``` r
get_regime_series(y_mat, regime_basis = "market")
```

## Arguments

- y_mat:

  Numeric matrix (T x N) of returns.

- regime_basis:

  Character or integer specifying how to identify regimes:

  "market"

  :   Use equal-weighted average (rowMeans).

  "first_pc"

  :   Use first principal component.

  integer

  :   Use the specified column index.

## Value

Numeric vector of length T for regime identification.
