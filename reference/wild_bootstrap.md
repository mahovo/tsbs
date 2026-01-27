# Wild Bootstrap for Time Series Residuals

Generates wild bootstrap replicates of a vector or matrix of residuals
by multiplying each observation by a random Rademacher weight (+1 or
-1).

## Usage

``` r
wild_bootstrap(x, num_boots = 100, parallel = FALSE, num_cores = 2)
```

## Arguments

- x:

  Numeric vector or matrix of residuals.

- num_boots:

  Integer number of bootstrap replicates.

- parallel:

  Parallelize computation? `TRUE` or `FALSE`.

- num_cores:

  Number of cores.

## Value

A list of numeric matrices, each one a wild bootstrap replicate.

## Details

The wild bootstrap is often used to resample regression or model
residuals when heteroskedasticity or other non-i.i.d. errors are
present. Each replicate is constructed by multiplying every observation
by +1 or -1, where the signs are drawn randomly with equal probability.

## References

A. Colin Cameron & Jonah B. Gelbach & Douglas L. Miller, 2008.
"Bootstrap-Based Improvements for Inference with Clustered Errors", The
Review of Economics and Statistics, MIT Press, vol. 90(3), pages
414-427, August.

## Examples

``` r
set.seed(123)
resids <- rnorm(100)
boot_reps <- wild_bootstrap(resids, num_boots = 5)
length(boot_reps)           # 5 replicates
#> [1] 5
dim(boot_reps[[1]])         # 100 x 1 matrix
#> [1] 100   1
```
