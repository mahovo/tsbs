# Emission Density Evaluation

Evaluates the emission density for each observation given state
parameters.

## Usage

``` r
.emission_density(y, params, distribution)
```

## Arguments

- y:

  Numeric vector of observations.

- params:

  List with mu, sigma, nu (for t), xi (for sstd).

- distribution:

  Character, distribution type.

## Value

Numeric vector of density values.
