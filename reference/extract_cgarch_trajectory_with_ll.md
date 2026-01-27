# Extract CGARCH Trajectory with Log-Likelihood

Extract CGARCH parameter evolution along with log-likelihood values.

## Usage

``` r
extract_cgarch_trajectory_with_ll(diagnostics, state = 1)
```

## Arguments

- diagnostics:

  An object of class `ms_diagnostics`

- state:

  Integer state index (default 1)

## Value

Data frame with columns: iteration, alpha, beta, gamma, shape, log_lik
