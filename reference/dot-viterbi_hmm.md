# Viterbi Algorithm for HMM

Finds the most likely state sequence using the Viterbi algorithm.

## Usage

``` r
.viterbi_hmm(y, trans_mat, init_probs, state_params, distribution)
```

## Arguments

- y:

  Numeric vector of observations.

- trans_mat:

  Transition matrix.

- init_probs:

  Initial state probabilities.

- state_params:

  Emission parameters per state.

- distribution:

  Emission distribution type.

## Value

Integer vector of most likely states.
