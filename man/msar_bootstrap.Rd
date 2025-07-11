% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrap_functions.R
\name{msar_bootstrap}
\alias{msar_bootstrap}
\title{Markov-Switching Autoregressive (MSAR) Bootstrap for Time Series}
\usage{
msar_bootstrap(
  x,
  ar_order = 1,
  num_states = 2,
  n_boot = NULL,
  num_blocks = 100,
  num_boots = 100
)
}
\arguments{
\item{x}{Numeric vector representing the time series.}

\item{ar_order}{Integer order of the autoregressive model.}

\item{num_states}{Integer number of regimes (hidden states) in the MSAR model.}

\item{n_boot}{Length of bootstrap series.}

\item{num_blocks}{Integer number of blocks to sample per bootstrap replicate.}

\item{num_boots}{Integer number of bootstrap replicates.}
}
\value{
A list of numeric vectors, each a bootstrap replicate.
}
\description{
Fits a Markov-switching autoregressive model (MSAR) to a univariate time series
and generates bootstrap replicates by resampling regime-specific blocks.
}
\details{
This function:
\itemize{
\item Fits a Markov-switching autoregressive model using \code{MSwM::msmFit()} on an \code{lm()} fit of \code{x}.
\item Uses Viterbi decoding (\code{ms_model@Fit@smoProb}) to classify each observation into states.
\item Groups contiguous observations belonging to the same state into blocks.
\item Samples these regime-specific blocks with replacement to generate synthetic series.
}
If \code{n_boot} is set, the last block will be truncated when necessary to match
the length (\code{n_boot}) of the bootstrap series. This is the only way to
ensure equal length of all bootstrap series, as the length of each block is
random. If \code{n_boot} is not set, \code{num_blocks} must be set, and the length of
each bootstrap series will be determined by the number of blocks and the
random lengths of the individual blocks for that particular series. Note
that this almost certainly results in bootstrap series of different lengths.
Note that if \code{n_boot} and \code{num_blocks} are both set, \code{num_blocks} will be
ignored.
}
\examples{
set.seed(123)
x <- arima.sim(n = 200, list(ar = 0.7))
msar_samples <- msar_bootstrap(x, ar_order = 1, num_states = 2, num_blocks = 10, num_boots = 5)
length(msar_samples)     # 5 replicates
length(msar_samples[[1]]) # length of one bootstrap series

}
