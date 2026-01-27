# Add Probability Visualization Layer to Plot

Adds the appropriate ggplot layer for probability visualization based on
style.

## Usage

``` r
add_probability_layer(p, prob_data, style, state_colors, plot_data)
```

## Arguments

- p:

  Existing ggplot object.

- prob_data:

  Probability data frame.

- style:

  One of "ribbon", "line", "bands_alpha".

- state_colors:

  Named color vector.

- plot_data:

  Series plot data (for y-axis scaling).

## Value

Updated ggplot object.
