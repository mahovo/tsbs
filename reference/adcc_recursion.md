# ADCC(1,1) Recursion

Computes Q and R matrices for ADCC(1,1) model. ADCC adds asymmetric
response to negative shocks via gamma parameter.

The ADCC recursion is: Q_t = Ω + α\*(z_t-1z'*t-1) + γ\*(n*t-1n'*t-1) +
β\*Q*t-1

where n_t = z_t \* I(z_t \< 0) (element-wise negative indicator) and Ω =
(1 - α - β)*Qbar - γ*Nbar

## Usage

``` r
adcc_recursion(std_resid, Qbar, alpha, gamma, beta, Nbar = NULL)
```

## Arguments

- std_resid:

  Matrix of standardized residuals (T x k)

- Qbar:

  Unconditional covariance matrix (k x k)

- alpha:

  ADCC alpha parameter (response to shocks)

- gamma:

  ADCC gamma parameter (asymmetric response to negative shocks)

- beta:

  ADCC beta parameter (persistence)

- Nbar:

  Unconditional covariance of negative shocks (k x k), computed if NULL

## Value

List with success, Q, R, Nbar, and error info if failed
