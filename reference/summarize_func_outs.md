# Summarize Bootstrap Functional Outputs

Computes summary statistics (mean, SD, confidence intervals) for
bootstrap functional outputs such as portfolio weights or other derived
quantities.

## Usage

``` r
summarize_func_outs(func_outs, names = NULL, probs = c(0.025, 0.975))
```

## Arguments

- func_outs:

  List of functional outputs from
  [`tsbs()`](https://mahovo.github.io/tsbs/reference/tsbs.md), or a
  matrix where each row is a bootstrap replicate and columns are output
  dimensions.

- names:

  Optional character vector of names for output dimensions.

- probs:

  Numeric vector of probabilities for quantile calculation. Defaults to
  `c(0.025, 0.975)` for 95% confidence intervals.

## Value

A data frame with columns:

- Name:

  Dimension name

- Mean:

  Bootstrap mean

- SD:

  Bootstrap standard deviation

- CI_Lower:

  Lower confidence bound

- CI_Upper:

  Upper confidence bound

- CI_Width:

  Width of confidence interval

## Examples

``` r
if (FALSE) { # \dontrun{
# After running tsbs() with a portfolio function
result <- tsbs(x, bs_type = "ms_varma_garch", func = risk_parity_portfolio, ...)

# Summarize the bootstrap weights
summary_df <- summarize_func_outs(result$func_outs, names = c("SPY", "EFA", "BND"))
print(summary_df)
} # }
```
