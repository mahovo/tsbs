# Manual SPD Computation

Computes SPD transformation manually without tsdistributions. Uses
kernel density estimation for interior and empirical tails.

## Usage

``` r
compute_spd_manual(z, lower_threshold = 0.1, upper_threshold = 0.9)
```

## Arguments

- z:

  Standardized residuals

- lower_threshold:

  Lower quantile threshold

- upper_threshold:

  Upper quantile threshold

## Value

Vector of uniform values
