# Create Probability Data Frame for Plotting

Converts a probability matrix to long format suitable for ggplot.

## Usage

``` r
create_probability_df(probs, series_label, state_labels)
```

## Arguments

- probs:

  Matrix of probabilities (n x num_states).

- series_label:

  Character label for this series.

- state_labels:

  Character vector of state names.

## Value

Data frame with columns: Index, State, Probability, Series.
