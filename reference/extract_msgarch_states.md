# Extract State Information from MSGARCH Fit

Extracts decoded state sequence, transition matrix, and state
probabilities from a fitted MSGARCH model.

## Usage

``` r
extract_msgarch_states(fit, y)
```

## Arguments

- fit:

  An MSGARCH fit object (output from `fit_msgarch_model` or
  [`MSGARCH::FitML`](https://rdrr.io/pkg/MSGARCH/man/FitML.html)).

- y:

  Numeric vector, the original data used for fitting. Required for
  computing residuals.

## Value

A list containing:

- states:

  Integer vector of Viterbi-decoded state assignments.

- transition_matrix:

  K x K matrix of transition probabilities.

- smoothed_probs:

  T x K matrix of smoothed state probabilities.

- n_states:

  Integer, number of states K.

- coefficients:

  Named vector of estimated model parameters.

- state_durations:

  Numeric vector of expected state durations (1 / (1 - p_ii) for each
  state i).
