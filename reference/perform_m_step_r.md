# Perform the M-Step (Sequential Execution)

This function is called once per EM iteration from C++. It estimates the
parameters for all M states sequentially, properly handling diagnostics
collection.

## Usage

``` r
perform_m_step_r(
  y,
  weights,
  spec,
  model_type,
  diagnostics = NULL,
  iteration = NULL,
  verbose = FALSE,
  dcc_threshold = 0.02,
  dcc_criterion = "bic"
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

  Diagnostic collector object (ms_diagnostics class).

- iteration:

  Current EM iteration number.

- verbose:

  Logical indicating whether to print progress.

- dcc_threshold:

  dcc_threshold

- dcc_criterion:

  dcc_criterion

## Value

A list of length M containing the updated model fits for each state.
