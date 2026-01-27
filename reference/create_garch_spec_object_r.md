# Create a GARCH Specification Object (Convenience Function)

Convenience function that handles the complex, model-specific logic for
creating both univariate and multivariate GARCH specification objects
for the purpose of log-likelihood calculation during the EM-algorithm.

## Usage

``` r
create_garch_spec_object_r(
  residuals,
  spec,
  model_type,
  current_pars,
  diagnostics = NULL,
  iteration = NULL,
  state = NULL,
  verbose = FALSE
)
```

## Arguments

- residuals:

  Numeric

- spec:

  A spec list

- model_type:

  Character string

- current_pars:

  A list of parameters

- diagnostics:

  diagnostics

- iteration:

  iteration

- state:

  state

- verbose:

  verbose

## Value

A GARCH specification object
