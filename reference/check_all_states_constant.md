# Check for All-States-Constant Event

Determines if and when all states switched to constant correlation.

## Usage

``` r
check_all_states_constant(diagnostics)
```

## Arguments

- diagnostics:

  An object of class `ms_diagnostics`

## Value

List with:

- occurred:

  Logical indicating if all states became constant

- iteration:

  Iteration number when it occurred (NA if not)

- message:

  Description of the event
