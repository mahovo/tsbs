# Fit MSGARCH Model

Thin wrapper around MSGARCH::FitML with sensible defaults for bootstrap
use. Provides consistent interface and error handling for the tsbs
package.

## Usage

``` r
fit_msgarch_model(
  y,
  n_states = 2L,
  variance_model = c("sGARCH", "eGARCH", "gjrGARCH", "tGARCH"),
  distribution = c("sstd", "std", "norm", "snorm", "sged", "ged"),
  ...
)
```

## Arguments

- y:

  Numeric vector of returns (univariate).

- n_states:

  Integer, number of regimes (default 2).

- variance_model:

  Character, GARCH specification: one of `"sGARCH"`, `"eGARCH"`,
  `"gjrGARCH"`, `"tGARCH"`. Default is `"sGARCH"`.

- distribution:

  Character, conditional distribution: one of `"norm"`, `"std"`,
  `"ged"`, `"snorm"`, `"sstd"`, `"sged"`. Default is `"sstd"` (skew
  Student-t).

- ...:

  Additional arguments passed to
  [`MSGARCH::FitML`](https://rdrr.io/pkg/MSGARCH/man/FitML.html).

## Value

An MSGARCH fit object (class `MSGARCH_ML_FIT`).

## Details

This function creates a Markov-switching GARCH specification using the
Haas et al. (2004) approach, which avoids the path-dependency problem by
defining K separate GARCH processes that evolve in parallel.

The skewed distributions use the Fernández & Steel (1998)
transformation, which introduces a skewness parameter \\\xi \> 0\\ where
\\\xi = 1\\ corresponds to a symmetric distribution.

## See also

[`CreateSpec`](https://rdrr.io/pkg/MSGARCH/man/CreateSpec.html),
[`FitML`](https://rdrr.io/pkg/MSGARCH/man/FitML.html)
