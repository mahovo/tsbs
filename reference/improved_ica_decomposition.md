# Improved ICA Decomposition for GOGARCH

Performs ICA decomposition with multiple restarts and quality
diagnostics. Uses the RADICAL algorithm from `tsmarch`.

## Usage

``` r
improved_ica_decomposition(
  residuals,
  method = "radical",
  n_components = NULL,
  n_restarts = 3,
  demean = FALSE,
  verbose = FALSE
)
```

## Arguments

- residuals:

  T x k matrix of residuals

- method:

  Currently only `"radical"` is supported

- n_components:

  Number of components to extract (default: k)

- n_restarts:

  Number of random restarts (default: 3)

- demean:

  Whether to demean the data (default: FALSE)

- verbose:

  Print progress information

## Value

List with components:

- S:

  Independent components matrix (T x n_components)

- A:

  Mixing matrix (k x n_components)

- W:

  Unmixing matrix (n_components x k)

- K:

  Pre-whitening matrix

- method:

  "radical" or "pca_fallback"

- quality:

  Quality metrics (see Details)

- convergence_info:

  Restart diagnostics

## Details

Improvements over direct
[`radical`](https://rdrr.io/pkg/tsmarch/man/radical.html) calls:

- Multiple restarts with best-fit selection based on negentropy

- Quality metrics: independence score (target \>0.8), reconstruction
  error (target \<1\\

- Graceful PCA fallback if ICA fails

## See also

[`estimate_garch_weighted_gogarch`](https://mahovo.github.io/tsbs/reference/estimate_garch_weighted_gogarch.md),
[`gogarch_diagnostics`](https://mahovo.github.io/tsbs/reference/gogarch_diagnostics.md)
