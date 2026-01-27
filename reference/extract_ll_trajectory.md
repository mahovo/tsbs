# Extract Log-Likelihood Trajectory

Extract the complete log-likelihood evolution across EM iterations.

## Usage

``` r
extract_ll_trajectory(diagnostics, type = c("after", "before", "change"))
```

## Arguments

- diagnostics:

  An object of class `ms_diagnostics`

- type:

  Character: "before", "after", or "change"

## Value

Numeric vector of log-likelihood values

## Examples

``` r
if (FALSE) { # \dontrun{
ll_after <- extract_ll_trajectory(diag, type = "after")
ll_changes <- extract_ll_trajectory(diag, type = "change")
plot(ll_after, type = "b", main = "LL Evolution")
} # }
```
