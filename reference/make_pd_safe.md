# Make Matrix Positive Definite (Safe Version)

Ensures a matrix is positive definite by adjusting eigenvalues.

## Usage

``` r
make_pd_safe(R, min_eig = 1e-08)
```

## Arguments

- R:

  Symmetric matrix

- min_eig:

  Minimum eigenvalue (default 1e-8)

## Value

Positive definite matrix, or NULL if input contains NA/Inf
