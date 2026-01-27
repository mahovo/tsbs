# Create Regime Bands Data Frame

Helper function to convert a state sequence into a data frame of regime
bands suitable for geom_rect().

## Usage

``` r
create_regime_bands_df(states, series_label)
```

## Arguments

- states:

  Integer vector of state assignments.

- series_label:

  Character string identifying this series.

## Value

Data frame with columns: xmin, xmax, State, Series.
