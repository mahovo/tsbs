# Extract State Info from HMM-SSTD Fit

Wrapper to extract state information from a `hmm_sstd_fit` object in the
same format as `extract_msgarch_states`.

## Usage

``` r
extract_hmm_sstd_states(fit, y = NULL)
```

## Arguments

- fit:

  Object of class `hmm_sstd_fit`.

- y:

  Original data (for compatibility, not used).

## Value

List compatible with `extract_msgarch_states` output.
