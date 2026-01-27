# Reconstruct Bootstrap Probabilities from Source Blocks

For a bootstrap replicate, reconstruct the "probabilities" by looking up
the original probabilities at the source time points of each resampled
block.

## Usage

``` r
reconstruct_bootstrap_probabilities(regime_info, replicate_idx, original_probs)
```

## Arguments

- regime_info:

  Regime info from diagnostics.

- replicate_idx:

  Which replicate.

- original_probs:

  Original smoothed probability matrix.

## Value

Matrix of reconstructed probabilities or NULL.
