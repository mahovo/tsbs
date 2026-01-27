# Gradient of NLL w.r.t. Q_t (Copula)

Backpropagate gradient from R_t to Q_t through the normalization. R_t =
D^-1 Q_t D^-1 where D = diag(sqrt(diag(Q_t)))

## Usage

``` r
grad_R_to_Q_copula(grad_R, Q_t, R_t)
```

## Arguments

- grad_R:

  k x k gradient w.r.t. R_t

- Q_t:

  k x k Q matrix

- R_t:

  k x k correlation matrix

## Value

k x k gradient w.r.t. Q_t
