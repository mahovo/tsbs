# Perform DCC(p,q) recursion

Computes the Q and R matrices for arbitrary DCC(p,q) order: Q_t = Qbar
\* (1 - sum(alpha) - sum(beta)) + sum_j=1^q alpha_j \* (z_t-j \*
z_t-j') + sum_j=1^p beta_j \* Q_t-j R_t = diag(Q_t)^-1/2 \* Q_t \*
diag(Q_t)^-1/2

## Usage

``` r
dcc_recursion(std_resid, Qbar, alphas, betas, verbose = FALSE)
```

## Arguments

- std_resid:

  Matrix of standardized residuals (T x k)

- Qbar:

  Unconditional covariance matrix of standardized residuals (k x k)

- alphas:

  Numeric vector of alpha parameters (length q)

- betas:

  Numeric vector of beta parameters (length p)

- verbose:

  Logical; if TRUE, print diagnostic messages

## Value

List with components:

- success:

  Logical indicating if recursion completed without errors

- Q:

  Array of Q matrices (k x k x T)

- R:

  Array of correlation matrices (k x k x T)

- maxpq:

  Maximum of p and q, used for burn-in period

- error_type:

  Character describing error if success is FALSE

- error_time:

  Time index where error occurred
