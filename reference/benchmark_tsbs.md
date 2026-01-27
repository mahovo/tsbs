# Benchmark Bootstrap Methods

Compares the runtime of different
[`tsbs()`](https://mahovo.github.io/tsbs/reference/tsbs.md)
configurations by varying a single parameter (e.g., series length,
number of assets, or number of bootstrap replicates).

## Usage

``` r
benchmark_tsbs(
  x,
  setups,
  vary = c("num_boots", "n_boot", "n_assets"),
  values,
  times = 3,
  verbose = TRUE
)
```

## Arguments

- x:

  Numeric matrix of input data for bootstrapping.

- setups:

  A named list of
  [`tsbs()`](https://mahovo.github.io/tsbs/reference/tsbs.md) argument
  lists. Each element should be a list of arguments to pass to
  [`tsbs()`](https://mahovo.github.io/tsbs/reference/tsbs.md), excluding
  `x` and the varying parameter.

- vary:

  Character string specifying which parameter to vary. One of:

  `"n_boot"`

  :   Vary the length of bootstrap series

  `"num_boots"`

  :   Vary the number of bootstrap replicates

  `"n_assets"`

  :   Vary the number of assets (columns)

- values:

  Numeric vector of values for the varying parameter.

- times:

  Integer number of times to repeat each benchmark for averaging.
  Default is 3.

- verbose:

  Logical. If TRUE, print progress messages. Default is TRUE.

## Value

An object of class `tsbs_benchmark` containing:

- results:

  Data frame with columns: Setup, Parameter, Value, Time, Rep

- summary:

  Data frame with mean and sd of times by Setup and Value

- vary:

  The parameter that was varied

- setups:

  The setup configurations used

## Details

This function is useful for:

- Comparing computational cost of different bootstrap methods

- Understanding how runtime scales with data size

- Choosing appropriate settings for production use

## See also

[`plot.tsbs_benchmark`](https://mahovo.github.io/tsbs/reference/plot.tsbs_benchmark.md)
for visualization.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create test data
set.seed(123)
x <- matrix(rnorm(500 * 3), ncol = 3)

# Define setups to compare
setups <- list(
  "Plain" = list(bs_type = "moving", block_length = 1),
  "Moving" = list(bs_type = "moving", block_length = 5),
  "Stationary" = list(bs_type = "stationary")
)

# Benchmark varying number of replicates
bench <- benchmark_tsbs(
  x = x,
  setups = setups,
  vary = "num_boots",
  values = c(10, 25, 50, 100),
  times = 3
)

# View results
print(bench)
plot(bench)
} # }
```
