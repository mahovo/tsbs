# Compute Stationary Distribution of Markov Chain

Computes the stationary (ergodic) distribution of a finite-state Markov
chain by solving \\\pi P = \pi\\ subject to \\\sum_k \pi_k = 1\\.

## Usage

``` r
compute_stationary_dist(trans_mat)
```

## Arguments

- trans_mat:

  Square matrix of transition probabilities. Row i, column j gives
  P(S_t+1 = j \| S_t = i). Rows must sum to 1.

## Value

Numeric vector of length K giving the stationary probabilities.

## Details

Solves the system \\(P' - I) \pi = 0\\ with the constraint \\\sum \pi =
1\\ by replacing one equation with the normalization constraint.
