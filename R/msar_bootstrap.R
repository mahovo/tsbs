#' Markov Switching Autoregressive (MSAR) Bootstrap for Time Series
#'
#' Fits an MSAR model and generates bootstrap replicates.
#'
#' @param x A numeric matrix with one column (univariate time series).
#' @param k Number of states in the MSAR model.
#' @param order AR order of each regime.
#' @param num_boots Number of bootstrap samples to generate.
#'
#' @return A list of bootstrap replicates (numeric vectors).
#'
#' @examples
#' set.seed(456)
#' returns <- arima.sim(n = 200, list(ar = 0.7))
#' boot_samples <- msar_bootstrap(returns, ar_order = 1, nstates = 2, num_blocks = 10, num_boots = 50)
#'
#' @importFrom MSwM msmFit
#' @importFrom stats predict
#' @export
msar_bootstrap <- function(returns, ar_order = 1, nstates = 2, num_blocks = 100, num_boots = 100, shuffle_blocks = TRUE) {
  if (length(returns) < (ar_order + 1) * nstates * 10) {
    stop("Time series too short for reliable MSAR estimation.")
  }
  
  df <- data.frame(y = returns)
  for (i in 1:ar_order) df[[paste0("lag", i)]] <- dplyr::lag(returns, i)
  df <- na.omit(df)
  
  base_model <- lm(y ~ ., data = df)
  ms_model <- msmFit(base_model, k = nstates, sw = rep(TRUE, length(coef(base_model)) + 1))
  
  states <- ms_model@Fit@smoProb
  state_seq <- apply(states, 1, which.max)
  
  r <- rle(state_seq)
  ends <- cumsum(r$lengths)
  starts <- c(1, head(ends, -1) + 1)
  blocks <- Map(function(s, e, st) list(start = s, end = e, state = st), starts, ends, r$values)
  
  trend_blocks <- Filter(function(b) b$state == 1, blocks)
  noise_blocks <- Filter(function(b) b$state != 1, blocks)
  
  if (length(trend_blocks) == 0 || length(noise_blocks) == 0) {
    stop("MSAR estimation produced only one regime — bootstrap requires at least two.")
  }
  
  get_block <- function(b) returns[b$start:b$end]
  
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
