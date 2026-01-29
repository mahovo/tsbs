# Update Emission Parameters (M-step)

Updates emission distribution parameters using weighted MLE.

## Usage

``` r
.update_emission_params(y, weights, params, distribution)
```

## Arguments

- y:

  Numeric vector of observations.

- weights:

  Numeric vector of state responsibilities (gamma_k).

- params:

  Current parameter values.

- distribution:

  Distribution type.

## Value

Updated parameter list.
