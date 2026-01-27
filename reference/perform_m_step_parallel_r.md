# Perform the M-Step in Parallel (R Helper)

\*\*\* THIS FUNCTION IS NOW DEPRECATED. Using perform_m_step_r() instead
\*\*\* This function is called once per EM iteration from C++. It uses
the 'future' framework to estimate the parameters for all M states in
parallel. Handles both univariate and multivariate models with proper
parameter structuring.

## Usage

``` r
perform_m_step_parallel_r(
  y,
  weights,
  spec,
  model_type,
  diagnostics = NULL,
  iteration = NULL,
  verbose = FALSE
)
```

## Arguments

- y:

  The time series data.

- weights:

  The (T x M) matrix of smoothed probabilities from the E-step.

- spec:

  The full list of model specifications.

- model_type:

  "univariate" or "multivariate".

- diagnostics:

  diagnostics

- iteration:

  iteration

- verbose:

  verbose

## Value

A list of length M containing the updated model fits for each state.
