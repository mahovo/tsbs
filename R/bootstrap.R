#' General Bootstrap Function for Time Series with Flexible Block Structures
#'
#' Performs time-series-aware block bootstrap using moving, stationary, hidden Markov model (HMM),
#' Markov-switching autoregressive (MSAR), or wild bootstrap. Includes automatic estimation of block parameters.
#'
#' @param x A numeric vector, matrix, or data frame of time series data.
#' @param n_length Desired length of the bootstrap series. If \code{num_blocks} is given, \code{n_length} is overridden.
#' @param block_length Length of individual blocks (ignored for MSAR and HMM types).
#' @param block_type Type of block sampling: "moving", "stationary", "hmm", "msar", or "wild".
#' @param func A function to apply to each bootstrap sample (e.g. mean, standard deviation).
#' @param apply_func_to Whether to apply function to each column ("cols") or full sample ("df").
#' @param num_blocks Number of blocks per bootstrap sample.
#' @param num_boots Number of bootstrap replicates.
#' @param p_method Method to estimate stationary bootstrap p: "1/n", "plugin", or "cross validation".
#' @param p Optional value for p. If \code{NULL}, uses \code{p_method}.
#' @param overlap Allow overlapping blocks (for moving/stationary block types).
#' @param ar_order AR order for MSAR model.
#' @param nstates Number of regimes for HMM/MSAR models.
#' @param model_func Model function used when \code{p_method="cross validation"}.
#' @param score_func Scoring function used when \code{p_method="cross validation"}.
#'
#' @return A list containing:
#' \describe{
#'   \item{bootstrap_series}{List of bootstrapped time series.}
#'   \item{func_outs}{List of outputs of \code{func} applied to each sample.}
#'   \item{func_out_means}{Mean output across bootstrap samples.}
#' }
#'
#' @examples
#' set.seed(42)
#' x <- arima.sim(n = 100, list(ar = 0.8))
#' result <- bootstrap(
#'   x = x,
#'   block_length = 10,
#'   block_type = "stationary",
#'   num_blocks = 5,
#'   num_boots = 100,
#'   func = mean,
#'   apply_func_to = "cols",
#'   p_method = "plugin"
#' )
#' print(result$func_out_means)
#'
#' @export
bootstrap <- function(
    x, 
    n_length = 100,
    block_length = 10, 
    block_type = "stationary",
    func = mean, 
    apply_func_to = "cols",
    num_blocks = 10, 
    num_boots = 100, 
    p_method = "1/n",
    p = NULL, 
    overlap = TRUE, 
    ar_order = 1, 
    nstates = 2,
    model_func = NULL, 
    score_func = NULL) {
  
  block_type <- match.arg(block_type, choices = c("moving", "stationary", "hmm", "msar", "wild"))
  apply_func_to <- match.arg(apply_func_to, choices = c("cols", "df"))
  p_method <- match.arg(p_method, choices = c("1/n", "cross validation", "plugin"))
  
  # Ensure matrix
  x <- as.matrix(x)
  
  if (is.vector(x)) x <- as.matrix(x)
  if (is.data.frame(x)) x <- as.matrix(x)
  
  n <- nrow(x)
  d <- ncol(x)
  
  if (block_type == "stationary" && is.null(p)) {
    if (p_method == "1/n") {
      p <- 1 / block_length  
    } else if (p_method == "plugin") {
      ac1 <- acf(x, lag.max = 1, plot = FALSE)$acf[2,,1]  
      p <- 1 - abs(mean(ac1, na.rm = TRUE))  
    } else if (p_method == "cross validation") {
      if (is.null(model_func)) model_func <- default_model_func
      if (is.null(score_func)) score_func <- mse
      y <- x[, 1]  
      p <- k_fold_cv_ts(x, y, model_func = model_func, score_func = score_func)
    }
  }
  
  bootstrap_series <-
    if (block_type %in% c("moving", "stationary")) {
      blockBootstrap(x, n_length, block_length, num_blocks, num_boots, block_type, p, overlap)
    } else if (block_type == "hmm") {
      lapply(hmm_bootstrap(x[, 1], nstates = nstates,
                           num_blocks = num_blocks, num_boots = num_boots),
             function(series) matrix(series, ncol = 1))
    } else if (block_type == "msar") {
      lapply(msar_bootstrap(x[, 1], ar_order = ar_order, nstates = nstates,
                            num_blocks = num_blocks, num_boots = num_boots),
             function(series) matrix(series, ncol = 1))
    } else if (block_type == "wild") {
      wild_bootstrap(x, num_boots)
    } else {
      stop("Unsupported block_type.")
    }
  
  func_outs <- lapply(bootstrap_series, function(sampled) {
    if (apply_func_to == "df") {
      func(as.data.frame(sampled))
    } else {
      apply(sampled, 2, func)
    }
  })
  
  func_out_means <- Reduce("+", func_outs) / length(func_outs)
  
  return(list(
    bootstrap_series = bootstrap_series,
    func_outs = func_outs,
    func_out_means = func_out_means
  ))
}
