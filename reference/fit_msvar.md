# Fit a 2-State Markov-Switching Vector Autoregressive Model (MS-VAR)

This function fits a 2-state MS-VAR(1) model using a C++ implementation
of the Expectation-Maximization (EM) algorithm. Note: Here VAR stands
for Vector AutoRegressive (not "variance", not "Value At Risk").

## Usage

``` r
fit_msvar(y, max_iter = 100, tol = 1e-06)
```

## Arguments

- y:

  A numeric matrix or data frame where rows are observations and columns
  are the time series variables.

- max_iter:

  An integer specifying the maximum number of EM iterations. Defaults to
  100.

- tol:

  A numeric value specifying the convergence tolerance for the
  log-likelihood. Defaults to 1e-6.

## Value

A list containing the estimated model parameters:

- beta1, beta2:

  VAR coefficient matrices for each state.

- sigma1, sigma2:

  Error covariance matrices for each state.

- P:

  The 2x2 transition probability matrix.

- log_likelihood:

  The final log-likelihood value.

- smoothed_probabilities:

  A matrix of smoothed state probabilities.

## Examples

``` r
# Generate sample data
set.seed(123)
T_obs <- 250
y1 <- arima.sim(model = list(ar = 0.7), n = T_obs)
y2 <- 0.5 * y1 + arima.sim(model = list(ar = 0.3), n = T_obs)
sample_data <- cbind(y1, y2)

# Fit the model (assuming the package is loaded)
# msvar_fit <- fit_msvar(sample_data)

# View results
# print(msvar_fit$P)
# plot(msvar_fit$smoothed_probabilities[, 1], type = 'l',
#      main = "Smoothed Probability of State 1", ylab = "Probability")
```
