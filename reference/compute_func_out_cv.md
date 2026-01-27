# Compute Coefficient of Variation for Bootstrap Outputs

Calculates the coefficient of variation (CV) for each dimension of
bootstrap functional outputs, providing a measure of estimation
stability.

## Usage

``` r
compute_func_out_cv(
  func_outs,
  names = NULL,
  cv_thresholds = c(Stable = 0.3, Moderate = 0.6)
)
```

## Arguments

- func_outs:

  List of functional outputs from
  [`tsbs()`](https://mahovo.github.io/tsbs/reference/tsbs.md), or a
  matrix.

- names:

  Optional character vector of names for output dimensions.

- cv_thresholds:

  Named numeric vector with thresholds for stability classification.
  Defaults to `c(Stable = 0.3, Moderate = 0.6)`.

## Value

A data frame with columns:

- Name:

  Dimension name

- Mean:

  Bootstrap mean

- SD:

  Bootstrap standard deviation

- CV:

  Coefficient of variation (SD/Mean)

- Stability:

  Stability classification based on CV thresholds

## Details

The coefficient of variation (CV) is defined as SD/Mean. Lower CV values
indicate more stable estimates. Default thresholds classify estimates
as:

- **Stable**: CV \< 0.3

- **Moderate**: 0.3 \<= CV \< 0.6

- **Unstable**: CV \>= 0.6

## Examples

``` r
if (FALSE) { # \dontrun{
# After running tsbs() with a portfolio function
result <- tsbs(x, bs_type = "ms_varma_garch", func = risk_parity_portfolio, ...)

# Assess stability of bootstrap weights
stability_df <- compute_func_out_cv(result$func_outs, names = c("SPY", "EFA", "BND"))
print(stability_df)
} # }
```
