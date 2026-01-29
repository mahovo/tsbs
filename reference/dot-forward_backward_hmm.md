# Forward-Backward Algorithm for HMM

Computes smoothed state probabilities using the forward-backward
algorithm.

## Usage

``` r
.forward_backward_hmm(y, trans_mat, init_probs, state_params, distribution)
```

## Arguments

- y:

  Numeric vector of observations.

- trans_mat:

  Transition matrix.

- init_probs:

  Initial state probabilities.

- state_params:

  List of emission parameters per state.

- distribution:

  Emission distribution type.

## Value

List with gamma (smoothed probs), xi (transition probs), loglik.
