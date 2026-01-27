# Fit SPD Transformation for a Single Series

Fits a Semi-Parametric Distribution and computes PIT values. SPD
combines:

- Parametric tails from the fitted univariate distribution

- Kernel-smoothed empirical distribution in the center

## Usage

``` r
fit_spd_transform(
  z,
  uni_fit = NULL,
  dist_pars = NULL,
  lower_threshold = 0.1,
  upper_threshold = 0.9
)
```

## Arguments

- z:

  Standardized residuals for one series

- uni_fit:

  Univariate GARCH fit object (optional, for parametric tails)

- dist_pars:

  Distribution parameters (optional)

- lower_threshold:

  Lower quantile for parametric tail (default 0.1)

- upper_threshold:

  Upper quantile for parametric tail (default 0.9)

## Value

List with u (uniform values) and spd_model (fitted SPD object)
