# Compute Robust Estimates from Bootstrap Outputs

Computes various robust estimators (mean, median, winsorized mean,
conservative quantile) from bootstrap functional outputs.

## Usage

``` r
compute_robust_estimates(
  func_outs,
  names = NULL,
  point_est = NULL,
  trim = 0.1,
  conservative_quantile = 0.25
)
```

## Arguments

- func_outs:

  List of functional outputs from
  [`tsbs()`](https://mahovo.github.io/tsbs/reference/tsbs.md), or a
  matrix.

- names:

  Optional character vector of names for output dimensions.

- point_est:

  Optional numeric vector of point estimates to include in the
  comparison.

- trim:

  Trim proportion for winsorized mean. Defaults to 0.1.

- conservative_quantile:

  Quantile for conservative estimate. Defaults to 0.25.

## Value

A data frame with columns for each robust estimator.

## Details

This function computes:

- **Point**: Original point estimate (if provided)

- **Boot_Mean**: Bootstrap mean

- **Boot_Median**: Bootstrap median

- **Winsorized**: Winsorized mean (trimmed mean)

- **Conservative**: Lower quantile estimate

For portfolio weights, the conservative estimate is renormalized to sum
to 1.
