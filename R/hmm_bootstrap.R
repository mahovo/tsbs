#' Hidden Markov Model (HMM) Bootstrap for Time Series
#'
#' Fits a hidden Markov model to the time series and generates bootstrap replicates.
#'
#' @param x A numeric matrix or data frame with time series (rows: time points, columns: series).
#' @param nstates Number of hidden states in the HMM.
#' @param num_boots Number of bootstrap samples to generate.
#'
#' @return A list of bootstrap replicates, each a matrix of the same dimension as `x`.
#'
#' @examples
#' set.seed(123)
#' returns <- arima.sim(n = 200, list(ar = 0.5))
#' boot_samples <- hmm_bootstrap(returns, nstates = 2, num_blocks = 10, num_boots = 50)
#'
#' @importFrom depmixS4 depmix fit posterior
#' @importFrom stats predict
#' @export
hmm_bootstrap <- function(returns, nstates = 2, num_blocks = 100, num_boots = 100, shuffle_blocks = TRUE) {
  # Sanity check
  if (length(returns) < nstates * 10) {
    stop("Time series too short for reliable HMM estimation.")
  }
  
  df <- data.frame(return = returns)
  
  # Fit HMM model
  model <- depmix(return ~ 1, data = df, nstates = nstates, family = gaussian())
  fit <- fit(model, verbose = FALSE)
  states <- posterior(fit, type = "viterbi")$state
  
  # Identify contiguous blocks of each regime
  r <- rle(states)
  ends <- cumsum(r$lengths)
  starts <- c(1, head(ends, -1) + 1)
  blocks <- Map(function(s, e, st) list(start = s, end = e, state = st), starts, ends, r$values)
  
  # Separate blocks by regime
  trend_blocks <- Filter(function(b) b$state == 1, blocks)
  noise_blocks <- Filter(function(b) b$state != 1, blocks)
  
  # Handle edge case: if one regime missing
  if (length(trend_blocks) == 0 || length(noise_blocks) == 0) {
    stop("HMM estimation produced only one regime — bootstrap requires at least two.")
  }
  
  get_block <- function(b) returns[b$start:b$end]
  
  # Generate bootstrap replicates
  sampled_series <- replicate(num_boots, {
    sampled_trend <- lapply(sample(trend_blocks, num_blocks, replace = TRUE), get_block)
    sampled_noise <- lapply(sample(noise_blocks, num_blocks, replace = TRUE), get_block)
    
    all_blocks <- c(sampled_trend, sampled_noise)
    if (shuffle_blocks) {
      all_blocks <- sample(all_blocks)
    }
    
    unlist(all_blocks)
  }, simplify = FALSE)
  
  return(sampled_series)
}
