# Gradient of NLL w.r.t. Q_t (via R_t normalization)

Backpropagate gradient from R_t to Q_t through the normalization.

## Usage

``` r
grad_R_to_Q(grad_R, Q_t, R_t)
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
