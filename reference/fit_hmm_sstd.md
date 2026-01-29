# Fit HMM with Skew Student-t Emissions (No GARCH)

Fits a Hidden Markov Model with skew Student-t distributed emissions
directly to returns, without the GARCH volatility layer. Uses an EM
algorithm with numerical optimization for the M-step.

## Usage

``` r
fit_hmm_sstd(
  y,
  n_states = 2L,
  distribution = c("sstd", "std", "norm"),
  max_iter = 200L,
  tol = 1e-06,
  verbose = FALSE
)
```

## Arguments

- y:

  Numeric vector of returns.

- n_states:

  Integer, number of hidden states (default 2).

- distribution:

  Character, emission distribution. Currently supports `"sstd"` (skew
  Student-t, default), `"std"` (Student-t), and `"norm"` (Gaussian).

- max_iter:

  Integer, maximum EM iterations (default 200).

- tol:

  Numeric, convergence tolerance for log-likelihood (default 1e-6).

- verbose:

  Logical, print progress (default FALSE).

## Value

A list with class `"hmm_sstd_fit"` containing:

- transition_matrix:

  K x K transition probabilities.

- initial_probs:

  Length K initial state distribution.

- state_params:

  List of distribution parameters per state: mu (mean), sigma (sd), nu
  (df for t), xi (skewness for sstd).

- states:

  Viterbi-decoded state sequence.

- smoothed_probs:

  T x K smoothed state probabilities.

- loglik:

  Final log-likelihood.

- n_iter:

  Number of EM iterations.

- converged:

  Logical, whether EM converged.

## Details

This provides an alternative to MSGARCH when you want regime-switching
based purely on distributional characteristics (mean, variance,
skewness, kurtosis) without modeling GARCH dynamics.

The EM algorithm alternates between:

- **E-step**: Forward-backward algorithm to compute state probabilities
  given current parameters.

- **M-step**: Update distribution parameters by maximizing the expected
  complete-data log-likelihood.

The skew Student-t density uses the Fernández & Steel (1998)
parameterization as implemented in
[`fGarch::dsstd`](https://geobosh.github.io/fGarchDoc/reference/dist-sstd.html).

## See also

[`fit_msgarch_model`](https://mahovo.github.io/tsbs/reference/fit_msgarch_model.md)
for GARCH-based regime switching.
