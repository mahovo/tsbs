# Check DCC stationarity constraints

Verifies that DCC parameters satisfy:

1.  All alpha_j \>= 0

2.  All beta_j \>= 0

3.  sum(alphas) + sum(betas) \< 1

## Usage

``` r
check_dcc_stationarity(dcc_params, verbose = FALSE)
```

## Arguments

- dcc_params:

  Named list of DCC parameters

- verbose:

  Logical; if TRUE, print diagnostic messages

## Value

List with components:

- is_stationary:

  Logical indicating if constraints are satisfied

- persistence:

  Total persistence value

- reason:

  Character description if not stationary, NULL otherwise

- details:

  Output from compute_dcc_persistence()
