# Export Diagnostics to CSV Files

Export diagnostic data to multiple CSV files for external analysis.

## Usage

``` r
export_diagnostics_csv(diagnostics, prefix = "diagnostic", path = ".")
```

## Arguments

- diagnostics:

  An object of class `ms_diagnostics`

- prefix:

  Character string prefix for output files

- path:

  Directory path for output files (default: current directory)

## Value

Character vector of created file paths
