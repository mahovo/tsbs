# Plot NLL Contours with plotly

Create an interactive contour plot of the DCC likelihood surface.

## Usage

``` r
plot_nll_contours(
  surface,
  true_params = NULL,
  mle_params = NULL,
  title = "DCC(1,1) Negative Log-Likelihood Surface"
)
```

## Arguments

- surface:

  Result from compute_nll_surface()

- true_params:

  Optional vector c(alpha, beta) of true parameters

- mle_params:

  Optional vector c(alpha, beta) of MLE estimates

- title:

  Plot title

## Value

plotly object
