#' Flexible Block Bootstrap for Time Series
#'
#' Generates block bootstrap replicates of a numeric time series or multivariate time series.
#' Supports moving, stationary, HMM, MSAR, and wild block types.
#'
#' @param x Numeric vector, matrix, or data frame of time series observations (rows = time points, cols = variables).
#' @param n_length Integer, optional desired length of each bootstrap replicate.
#' @param block_length Integer length of each block; if `-1`, an automatic heuristic is used.
#' @param block_type Character string: one of `"moving"`, `"stationary"`, `"hmm"`, `"msar"`, or `"wild"`.
#' @param func A summary function to apply to each bootstrap replicate or column.
#' @param apply_func_to Character string: `"cols"` to apply columnwise or `"df"` to apply on the full data frame.
#' @param num_blocks Integer number of blocks per bootstrap replicate.
#' @param num_boots Integer number of bootstrap replicates.
#' @param p_method Character string to choose method for stationary bootstrap parameter: `"1/n"`, `"plugin"`, or `"cross validation"`.
#' @param p Optional numeric value for stationary bootstrap `p`.
#' @param overlap Logical indicating if overlapping blocks are allowed.
#' @param ar_order Integer AR order for MSAR (`block_type="msar"`).
#' @param nstates Integer number of states for HMM or MSAR.
#' @param model_func Model-fitting function for cross-validation.
#' @param score_func Scoring function for cross-validation.
#'
#' @return A list containing:
#' \describe{
#'   \item{bootstrap_series}{List of bootstrap replicate matrices.}
#'   \item{func_outs}{List of computed function outputs for each replicate.}
#'   \item{func_out_means}{Mean of the computed outputs across replicates.}
#' }
#' @examples
#' set.seed(123)
#' x <- arima.sim(n = 100, list(ar = 0.8))
#' result <- tsbs(
#'   x = x,
#'   block_length = 10,
#'   block_type = "stationary",
#'   num_blocks = 5,
#'   num_boots = 10,
#'   func = mean,
#'   apply_func_to = "cols"
#' )
#' print(result$func_out_means)
#'
#' @importFrom stats acf ar
#' @export
tsbs <- function(
    x,
    n_length = NULL,
    block_length = 10,
    block_type = c("moving", "stationary", "hmm", "msar", "wild"),
    func = mean,
    apply_func_to = c("cols", "df"),
    num_blocks = 10,
    num_boots = 100,
    p_method = c("1/n", "plugin", "cross validation"),
    p = NULL,
    overlap = TRUE,
    ar_order = 1,
    nstates = 2,
    model_func = default_model_func,
    score_func = mse) {
  
  # --- Input validation & coercion ---
  if (!is.function(func)) stop("`func` must be a valid function.")
  
  block_type <- match.arg(block_type)
  apply_func_to <- match.arg(apply_func_to)
  p_method <- match.arg(p_method)
  
  if (is.vector(x)) x <- as.matrix(x)
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.numeric(x)) stop("`x` must be numeric or coercible to numeric matrix.")
  
  n <- nrow(x); d <- ncol(x)
  if (is.null(n_length)) n_length <- n
  
  if (!is.numeric(num_blocks) || num_blocks < 1) stop("`num_blocks` must be a positive integer.")
  if (!is.numeric(num_boots) || num_boots < 1) stop("`num_boots` must be a positive integer.")
  
  # Estimate p if needed
  if (block_type == "stationary" && is.null(p)) {
    if (p_method == "1/n") {
      p <- 1 / block_length
    } else if (p_method == "plugin") {
      ac1 <- acf(x, lag.max = 1, plot = FALSE)$acf[2,,1]
      p <- 1 - abs(mean(ac1, na.rm = TRUE))
    } else if (p_method == "cross validation") {
      if (is.null(model_func) || is.null(score_func)) {
        stop("For p_method = 'cross validation', please provide `model_func` and `score_func`.")
      }
      y <- x[, 1]
      p <- k_fold_cv_ts(x, y, model_func, score_func)
    }
  }
  
  # Dispatch to appropriate bootstrap mechanism
  bootstrap_series <- switch(
    block_type,
    moving = blockBootstrap(x, n_length, block_length, num_blocks, num_boots, "moving", p, overlap),
    stationary = blockBootstrap(x, n_length, block_length, num_blocks, num_boots, "stationary", p, overlap),
    hmm = lapply(hmm_bootstrap(x[,1], nstates = nstates, num_blocks = num_blocks, num_boots = num_boots), function(s) matrix(s, ncol = 1)),
    msar = lapply(msar_bootstrap(x[,1], ar_order = ar_order, nstates = nstates, num_blocks = num_blocks, num_boots = num_boots), function(s) matrix(s, ncol = 1)),
    wild = wild_bootstrap(x, num_boots),
    stop("Unsupported block_type")
  )
  
  func_outs <- lapply(bootstrap_series, function(sampled) {
    if (apply_func_to == "df") {
      func(as.data.frame(sampled))
    } else {
      apply(sampled, 2, func)
    }
  })
  func_out_means <- Reduce(`+`, func_outs) / length(func_outs)
  
  list(
    bootstrap_series = bootstrap_series,
    func_outs = func_outs,
    func_out_means = func_out_means
  )
}
