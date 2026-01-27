# Simulate DCC-GARCH data for testing

Simulate DCC-GARCH data for testing

## Usage

``` r
simulate_dcc_garch_test_data(
  n = 200,
  k = 2,
  alpha = 0.05,
  beta = 0.9,
  seed = 123
)
```

## Arguments

- n:

  Number of observations

- k:

  Number of series

- alpha:

  DCC alpha parameter

- beta:

  DCC beta parameter

- seed:

  Random seed

## Value

List with y (data), true_params, Qbar
