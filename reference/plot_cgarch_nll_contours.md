# Plot CGARCH NLL Contours

Create an interactive contour plot of the CGARCH likelihood surface.

## Usage

``` r
plot_cgarch_nll_contours(
  surface,
  true_params = NULL,
  mle_params = NULL,
  title = NULL
)
```

## Arguments

- surface:

  Result from compute_cgarch_nll_surface()

- true_params:

  Optional vector c(alpha, beta) of true parameters

- mle_params:

  Optional vector c(alpha, beta) of MLE estimates

- title:

  Plot title

## Value

plotly object
