# Extract Smoothed Probabilities from Fitted Model

Attempts to extract smoothed state probabilities from a fitted HMM,
MS-VAR, or MS-VARMA-GARCH model object.

## Usage

``` r
extract_smoothed_probabilities(fit_object, regime_info)
```

## Arguments

- fit_object:

  Fitted model object (depmixS4, ms_var, or ms_varma_garch fit).

- regime_info:

  Regime info from diagnostics.

## Value

Matrix of smoothed probabilities (n x num_states) or NULL.
