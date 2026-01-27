# Compute Probability Integral Transform (PIT)

Transforms standardized residuals to uniform \[0,1\] margins

## Usage

``` r
compute_pit_transform(
  std_residuals,
  uni_fit_list,
  transformation = "parametric",
  copula_dist = "mvn",
  dist_pars = NULL
)
```

## Arguments

- std_residuals:

  Matrix of standardized residuals (T x k)

- uni_fit_list:

  List of fitted univariate GARCH models

- transformation:

  Type: "parametric", "empirical", or "spd"

- copula_dist:

  Copula distribution ("mvn" or "mvt")

- dist_pars:

  Distribution parameters

## Value

Matrix of uniform margins (T x k)
