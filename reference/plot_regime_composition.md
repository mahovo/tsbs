# Plot Regime Composition of Bootstrap Series

Creates a visualization comparing the original series with one or more
bootstrap replicates, showing the regime structure. Can display either
discrete regime bands or smoothed state probabilities.

## Usage

``` r
plot_regime_composition(
  tsbs_result,
  original_data,
  replicate_idx = 1,
  series_idx = 1,
  show_prices = TRUE,
  initial_price = 100,
  show_probabilities = FALSE,
  probability_style = c("ribbon", "line", "bands_alpha"),
  show_bands = NULL,
  state_colors = NULL,
  title = NULL
)
```

## Arguments

- tsbs_result:

  A list returned by
  [`tsbs()`](https://mahovo.github.io/tsbs/reference/tsbs.md) with
  `collect_diagnostics = TRUE`, or a `tsbs_diagnostics` object.

- original_data:

  The original data matrix used for bootstrapping.

- replicate_idx:

  Integer or vector of integers specifying which bootstrap replicate(s)
  to plot. Default is 1 (first replicate).

- series_idx:

  Integer specifying which column/series to plot if multivariate.
  Default is 1.

- show_prices:

  Logical. If TRUE and data appears to be returns, convert to cumulative
  prices for visualization. Default is TRUE.

- initial_price:

  Numeric. Starting price for cumulative calculation. Default is 100.

- show_probabilities:

  Logical. If TRUE, overlay smoothed state probabilities on the plot
  instead of (or in addition to) discrete bands. Requires that smoothed
  probabilities are available in the diagnostics. Default is FALSE.

- probability_style:

  Character. How to display probabilities:

  - `"ribbon"`: Stacked ribbons showing probability of each state

  - `"line"`: Line plot of probability for each state

  - `"bands_alpha"`: Regime bands with alpha proportional to probability

  Default is "ribbon".

- show_bands:

  Logical. If TRUE, show discrete regime bands. If
  `show_probabilities = TRUE` and `show_bands = TRUE`, both are shown.
  Default is TRUE when `show_probabilities = FALSE`.

- state_colors:

  Optional named vector of colors for each state.

- title:

  Optional plot title.

## Value

A ggplot object (invisibly).

## Details

This function creates a faceted plot with:

- Top panel: Original series with regime information

- Bottom panel(s): Bootstrap replicate(s) with their regime information

When `show_probabilities = TRUE`, the smoothed state probabilities from
the fitted model (HMM or MS-VARMA-GARCH) are displayed. This shows the
model's uncertainty about regime membership at each time point, which is
more informative than the hard Viterbi assignments.

For bootstrap replicates, the "probabilities" shown are actually the
probabilities from the original series at the source time points of each
resampled block. This reveals how the bootstrap combines observations
from different regime-certainty periods.

## Examples

``` r
if (FALSE) { # \dontrun{
# Run HMM bootstrap with diagnostics
result <- tsbs(
  x = returns_data,
  bs_type = "hmm",
  num_states = 2,
  num_boots = 10,
  collect_diagnostics = TRUE,
  return_fit = TRUE
)

# Plot with discrete regime bands
plot_regime_composition(result, returns_data)

# Plot with smoothed probabilities as ribbons
plot_regime_composition(result, returns_data, show_probabilities = TRUE)

# Plot with probability lines
plot_regime_composition(result, returns_data, 
                        show_probabilities = TRUE, 
                        probability_style = "line")
} # }
```
