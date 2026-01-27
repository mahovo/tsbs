# DCC(1,1) Q-Matrix Recursion with Gradient Storage

Compute Q_t matrices and store gradient information for backpropagation.

## Usage

``` r
dcc11_recursion_with_grad(z, alpha, beta, Qbar)
```

## Arguments

- z:

  T x k matrix of standardized residuals

- alpha:

  DCC alpha parameter

- beta:

  DCC beta parameter

- Qbar:

  k x k unconditional covariance matrix

## Value

List with Q, R, R_inv, log_det_R, dQ_dalpha, dQ_dbeta, success
