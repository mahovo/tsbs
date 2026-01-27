# Evaluate ICA Quality

Computes quality metrics for ICA decomposition.

## Usage

``` r
evaluate_ica_quality(S)
```

## Arguments

- S:

  T x n_components matrix of independent components

## Value

List with:

- independence_score:

  1 - mean(\|cor(S)\|), higher is better

- total_negentropy:

  Sum of component negentropies

- component_kurtosis:

  Excess kurtosis per component

- correlation_matrix:

  Component correlation matrix

## See also

[`improved_ica_decomposition`](https://mahovo.github.io/tsbs/reference/improved_ica_decomposition.md),
[`gogarch_diagnostics`](https://mahovo.github.io/tsbs/reference/gogarch_diagnostics.md)
