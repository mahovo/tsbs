# Record Regime/State Information for Bootstrap Diagnostics

Records the state sequence from the original data and the state
composition of each bootstrap replicate. This is used by HMM and
MS-VARMA-GARCH bootstrap methods to track how regimes are sampled.

## Usage

``` r
record_regime_info(
  diagnostics,
  original_states,
  num_states,
  state_labels = NULL
)
```

## Arguments

- diagnostics:

  A `tsbs_diagnostics` object.

- original_states:

  Integer vector of states for the original series.

- num_states:

  Integer, total number of possible states.

- state_labels:

  Optional character vector of state labels.

## Value

Updated diagnostics object.
