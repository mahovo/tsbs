# Compute DCC Persistence and Extract Parameters

Extracts all alpha/beta parameters and computes total persistence. For a
stationary DCC(p,q) model, we require: P = sum(alpha_1, ..., alpha_q) +
sum(beta_1, ..., beta_p) \< 1

## Usage

``` r
compute_dcc_persistence(dcc_params)
```

## Arguments

- dcc_params:

  Named list of DCC parameters

## Value

List with components:

- persistence:

  Total persistence = sum(alphas) + sum(betas)

- alpha_sum:

  Sum of all alpha parameters

- beta_sum:

  Sum of all beta parameters

- alphas:

  Numeric vector of alpha values in order

- betas:

  Numeric vector of beta values in order

- order:

  Named vector c(p, q)
