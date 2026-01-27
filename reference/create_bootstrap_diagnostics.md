# Create Bootstrap Diagnostic Collector

Initializes a diagnostic collector object for storing bootstrap metadata
and statistics during the bootstrap procedure.

## Usage

``` r
create_bootstrap_diagnostics(bs_type, n_original, n_vars, num_boots)
```

## Arguments

- bs_type:

  Character string specifying bootstrap type.

- n_original:

  Integer, length of original series.

- n_vars:

  Integer, number of variables/columns.

- num_boots:

  Integer, number of bootstrap replicates.

## Value

An object of class `tsbs_diagnostics`.
