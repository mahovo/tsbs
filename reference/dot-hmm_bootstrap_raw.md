# Raw Returns HMM Bootstrap (No GARCH)

Internal function implementing HMM bootstrap using standalone EM
algorithm with non-Gaussian emissions but without GARCH dynamics.

## Usage

``` r
.hmm_bootstrap_raw(
  x,
  n_boot = NULL,
  num_states = 2,
  num_boots = 100,
  distribution = "sstd",
  micro_block_length = 1L,
  regime_basis = "market",
  return_fit = FALSE,
  collect_diagnostics = FALSE,
  verbose = FALSE,
  seed = NULL,
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

- num_boots:

  Integer, number of bootstrap replicates to generate. Default is 100.

- distribution:

  Character, emission distribution for the HMM. One of:

  "gaussian"

  :   (Default) Gaussian emissions via depmixS4. This is the original
      behavior.

  "sstd"

  :   Skew Student-t via MSGARCH (with GARCH dynamics).

  "std"

  :   Student-t via MSGARCH (with GARCH dynamics).

  "snorm"

  :   Skew normal via MSGARCH (with GARCH dynamics).

  "sged"

  :   Skew GED via MSGARCH (with GARCH dynamics).

  "norm"

  :   Normal via MSGARCH (with GARCH dynamics).

  "ged"

  :   GED via MSGARCH (with GARCH dynamics).

  "sstd_raw"

  :   Skew Student-t without GARCH (standalone HMM-EM).

  "std_raw"

  :   Student-t without GARCH (standalone HMM-EM).

  "norm_raw"

  :   Normal without GARCH (standalone HMM-EM, equivalent to depmixS4
      but using our EM implementation).

- micro_block_length:

  Integer, block length for within-state sampling when using
  MSGARCH-based bootstrap. Use 1 (default) for iid sampling within
  states, or \>1 to preserve some local dependence. Ignored for Gaussian
  distribution which uses the traditional block sampling.

- regime_basis:

  For multivariate data with MSGARCH: how to identify regimes. One of
  `"market"` (equal-weighted average, default), `"first_pc"` (first
  principal component), or an integer column index. Ignored for Gaussian
  distribution.

- return_fit:

  Logical. If TRUE, returns the fitted model along with bootstrap
  samples. Default is FALSE.

- collect_diagnostics:

  Logical. If TRUE, collects detailed diagnostic information including
  regime composition. Default is FALSE.

- verbose:

  Logical. If TRUE, prints fitting information. Default is FALSE.

- seed:

  Integer, random seed for reproducibility. Default is NULL.

- ...:

  Additional arguments passed to the underlying fitting functions.
