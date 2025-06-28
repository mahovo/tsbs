#' Hidden Markov Model (HMM) Bootstrap for Time Series
#'
#' Fits a Gaussian Hidden Markov Model (HMM) to a univariate time series
#' and generates bootstrap replicates by resampling regime-specific blocks.
#'
#' @param x Numeric vector representing the time series.
#' @param num_states Integer number of hidden states for the HMM.
#' @param num_blocks Integer number of blocks to sample for each bootstrap replicate.
#' @param num_boots Integer number of bootstrap replicates to generate.
#'
#' @return A list of numeric vectors, each one a bootstrap replicate.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Fits a Gaussian HMM to `x` using `depmixS4::depmix()` and `depmixS4::fit()`.
#'   \item Uses Viterbi decoding (`posterior(fit, type = "viterbi")$state`)
#'     to assign each observation to a state.
#'   \item Samples contiguous blocks of observations belonging to each state.
#' }
#'
#' @examples
#' set.seed(123)
#' x <- arima.sim(n = 200, list(ar = 0.5))
#' hmm_samples <- hmm_bootstrap(x, num_states = 2, num_blocks = 10, num_boots = 5)
#' length(hmm_samples)     # 5 bootstrap replicates
#' length(hmm_samples[[1]]) # length of one bootstrap series
#'
#' @importFrom depmixS4 depmix fit posterior
#' @importFrom stats gaussian
#' @export
hmm_bootstrap <- function(x, num_states = 2, num_blocks = 100, num_boots = 100) {
  if (!is.numeric(x)) stop("`x` must be a numeric vector.")
  if (!is.numeric(num_states) || num_states < 1) stop("`num_states` must be a positive integer.")
  if (!is.numeric(num_blocks) || num_blocks < 1) stop("`num_blocks` must be a positive integer.")
  if (!is.numeric(num_boots) || num_boots < 1) stop("`num_boots` must be a positive integer.")
  
  df <- data.frame(return = x)
  
  # Fit HMM
  model <- depmixS4::depmix(return ~ 1, data = df, num_states = num_states, family = gaussian())
  fit <- depmixS4::fit(model, verbose = FALSE)
  states <- depmixS4::posterior(fit, type = "viterbi")$state
  
  # Extract regime-specific blocks
  r <- rle(states)
  ends <- cumsum(r$lengths)
  starts <- c(1, head(ends, -1) + 1)
  blocks <- Map(
    function(s, e) list(start = s, end = e),
    starts, ends
  )
  
  # Sample `num_blocks` blocks with replacement
  get_block <- function(b) x[b$start:b$end]
  sampled_series <- replicate(num_boots, {
    sampled_blocks <- lapply(sample(blocks, num_blocks, replace = TRUE), get_block)
    unlist(sampled_blocks)
  }, simplify = FALSE)
  
  sampled_series
}


#' Markov-Switching Autoregressive (MSAR) Bootstrap for Time Series
#'
#' Fits a Markov-switching autoregressive model (MSAR) to a univariate time series
#' and generates bootstrap replicates by resampling regime-specific blocks.
#'
#' @param x Numeric vector representing the time series.
#' @param ar_order Integer order of the autoregressive model.
#' @param num_states Integer number of regimes (hidden states) in the MSAR model.
#' @param num_blocks Integer number of blocks to sample per bootstrap replicate.
#' @param num_boots Integer number of bootstrap replicates.
#'
#' @return A list of numeric vectors, each a bootstrap replicate.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Fits a Markov-switching autoregressive model using `MSwM::msmFit()` on an `lm()` fit of `x`.
#'   \item Uses Viterbi decoding (`ms_model@Fit@smoProb`) to classify each observation into states.
#'   \item Groups contiguous observations belonging to the same state into blocks.
#'   \item Samples these regime-specific blocks with replacement to generate synthetic series.
#' }
#'
#' @examples
#' set.seed(123)
#' x <- arima.sim(n = 200, list(ar = 0.7))
#' msar_samples <- msar_bootstrap(x, ar_order = 1, num_states = 2, num_blocks = 10, num_boots = 5)
#' length(msar_samples)     # 5 replicates
#' length(msar_samples[[1]]) # length of one bootstrap series
#'
#' @importFrom MSwM msmFit
#' @importFrom stats lm
#' @export
msar_bootstrap <- function(x, ar_order = 1, num_states = 2, num_blocks = 100, num_boots = 100) {
  if (!is.numeric(x)) stop("`x` must be a numeric vector.")
  #if (!is.numeric(ar_order) || ar_order < 0) stop("`ar_order` must be a non-negative integer.")
  invalid_lengths(ar_order, allow_null = FALSE) stop("`ar_order` must be a non-negative integer.")
  #if (!is.numeric(num_states) || num_states < 1) stop("`num_states` must be a positive integer.")
  invalid_lengths(num_states, allow_null = FALSE) stop("`num_states` must be a positive integer.")
  #if (!is.numeric(num_blocks) || num_blocks < 1) stop("`num_blocks` must be a positive integer.")
  invalid_lengths(num_blocks, allow_null = FALSE) stop("`num_blocks` must be a positive integer.")
  #if (!is.numeric(num_boots) || num_boots < 1) stop("`num_boots` must be a positive integer.")
  invalid_lengths(num_boots, allow_null = FALSE) stop("`num_boots` must be a positive integer.")
  
  # Prepare lagged design matrix
  df <- data.frame(y = x)
  for (i in seq_len(ar_order)) df[[paste0("lag", i)]] <- dplyr::lag(x, i)
  df <- na.omit(df)
  
  # Fit base AR model
  base_model <- lm(y ~ ., data = df)
  
  # Fit MSAR model
  ms_model <- MSwM::msmFit(base_model, k = num_states, sw = rep(TRUE, length(coef(base_model)) + 1))
  
  # State sequence
  states <- ms_model@Fit@smoProb
  state_seq <- apply(states, 1, which.max)
  
  # Get contiguous blocks per state
  r <- rle(state_seq)
  ends <- cumsum(r$lengths)
  starts <- c(1, head(ends, -1) + 1)
  blocks <- Map(
    function(s, e) list(start = s, end = e),
    starts, ends
  )
  
  get_block <- function(b) x[b$start:b$end]
  
  # Sample num_blocks blocks per bootstrap
  sampled_series <- replicate(num_boots, {
    sampled_blocks <- lapply(sample(blocks, num_blocks, replace = TRUE), get_block)
    unlist(sampled_blocks)
  }, simplify = FALSE)
  
  sampled_series
}


#' Wild Bootstrap for Time Series Residuals
#'
#' Generates wild bootstrap replicates of a vector or matrix of residuals
#' by multiplying each observation by a random Rademacher weight (+1 or -1).
#'
#' @param x Numeric vector or matrix of residuals.
#' @param num_boots Integer number of bootstrap replicates.
#'
#' @return A list of numeric matrices, each one a wild bootstrap replicate.
#'
#' @details
#' The wild bootstrap is often used to resample regression or model residuals
#' when heteroskedasticity or other non-i.i.d. errors are present.
#' Each replicate is constructed by multiplying every observation by +1 or -1,
#' where the signs are drawn randomly with equal probability.
#'
#' @examples
#' set.seed(123)
#' resids <- rnorm(100)
#' boot_reps <- wild_bootstrap(resids, num_boots = 5)
#' length(boot_reps)           # 5 replicates
#' dim(boot_reps[[1]])         # 100 x 1 matrix
#'
#' @export
wild_bootstrap <- function(x, num_boots = 100) {
  if (!is.numeric(x)) stop("`x` must be numeric.")
  if (!is.numeric(num_boots) || num_boots < 1) stop("`num_boots` must be a positive integer.")
  
  # Coerce vector to a column matrix
  if (is.vector(x)) {
    x <- matrix(x, ncol = 1)
  }
  
  n <- nrow(x)
  d <- ncol(x)
  
  bootstrap_samples <- vector("list", num_boots)
  
  for (b in seq_len(num_boots)) {
    # Rademacher weights (+1 or -1)
    v <- sample(c(-1, 1), size = n, replace = TRUE)
    # Multiply each column by v
    bootstrap_samples[[b]] <- x * v
  }
  
  bootstrap_samples
}
