# Simulate Markov Chain

Generates a realization of a discrete-state Markov chain given
transition probabilities and initial distribution.

## Usage

``` r
simulate_markov_chain(T_len, trans_mat, init_dist = NULL)
```

## Arguments

- T_len:

  Integer, length of sequence to generate.

- trans_mat:

  Square matrix of transition probabilities.

- init_dist:

  Numeric vector of initial state probabilities. If NULL, uses the
  stationary distribution.

## Value

Integer vector of length T_len containing state assignments (1 to K).
