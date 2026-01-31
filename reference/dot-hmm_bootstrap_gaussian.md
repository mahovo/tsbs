# Gaussian HMM Bootstran

Internal function implementing the original Gaussian HMM bootstrap via
depmixS4. This preserves the original behavior for backward
compatibility.

## Usage

``` r
.hmm_bootstrap_gaussian(
  x,
  n_boot = NULL,
  num_states = 2,
  num_blocks = NULL,
  num_boots = 100,
  parallel = FALSE,
  num_cores = 1L,
  return_fit = FALSE,
  collect_diagnostics = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  Numeric vector or matrix representing the time series.

- n_boot:

  Integer, length of each bootstrap series. If NULL, defaults to the
  length of the original series.

- num_states:

  Integer, number of hidden states for the HMM. Default is 2.

- num_blocks:

  Integer, number of blocks to sample for each bootstrap replicate. Only
  used when `distribution = "gaussian"`.

- num_boots:

  Integer, number of bootstrap replicates to generate. Default is 100.

- parallel:

  Logical, parallelize computation? Default is FALSE.

- num_cores:

  Integer, number of cores for parallel processing.

- return_fit:

  Logical. If TRUE, returns the fitted model along with bootstrap
  samples. Default is FALSE.

- collect_diagnostics:

  Logical. If TRUE, collects detailed diagnostic information including
  regime composition. Default is FALSE.

- verbose:

  Logical. If TRUE, prints fitting information. Default is FALSE.

- ...:

  Additional arguments passed to the underlying fitting functions.
