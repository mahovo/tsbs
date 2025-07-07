#' Flexible Block Bootstrap for Time Series
#'
#' Generates block bootstrap replicates of a numeric time series or multivariate 
#' time series. Supports moving, stationary, HMM, MSAR, and wild block types.
#'
#' @param x Numeric vector, matrix, or data frame of time series observations 
#'   (rows = time points, cols = variables).
#' @param n_boot Integer, optional desired length of each bootstrap replicate.
#' @param block_length Integer length of each block; if `-1`, an automatic h
#'   euristic is used.
#' @param block_type Character string: one of `"moving"`, `"stationary"`, `
#'   "hmm"`, `"msar"`, or `"wild"`. See details below.
#' @param num_blocks Integer number of blocks per bootstrap replicate.
#' @param num_boots Integer number of bootstrap replicates.
#' @param func A summary function to apply to each bootstrap replicate or column.
#' @param apply_func_to Character string: `"cols"` to apply columnwise or `"df"` 
#'   to apply on the full data frame.
#' @param p_method Character string to choose method for stationary bootstrap 
#'   parameter: `"1/n"`, `"plugin"`, or `"cross validation"`.
#' @param p Optional numeric value for stationary bootstrap `p`.
#' @param overlap Logical indicating if overlapping blocks are allowed.
#' @param ar_order Integer AR order for MSAR (`block_type="msar"`).
#' @param num_states Integer number of states for HMM or MSAR.
#' @param model_func Model-fitting function for cross-validation.
#' @param score_func Scoring function for cross-validation.
#' 
#' @datails
#' `block_type="moving"`: If `n_boot` is set, the last block will be truncated
#'   when necessary to match the length (`n_boot`) of the bootstrap series. If 
#'   `n_boot` is not set, `block_length` and `num_blocks` must be set, and  
#'   `n_boot` will automatically be set to `block_length * num_blocks`.
#' `block_type="stationary"`, `block_type="hmm"`, `block_type="msar"`: If 
#'   `n_boot` is set, the last block will be truncated when necessary to match 
#'   the length (`n_boot`) of the bootstrap series. This is the only way to 
#'   ensure equal length of all bootstrap series, as the length of each block is 
#'   random. If `n_boot` is not set, `num_blocks` must be set, and the length of 
#'   each bootstrap series will be determined by the number of blocks and the 
#'   random lengths of the individual blocks for that particular series. Note 
#'   that this typically results in bootstrap series of different lengths.
#' : 
#' `block_type="wild"`: 
#'   
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
#' 
#' @useDynLib tsbs, .registration = TRUE
#' @export
tsbs <- function(
    x,
    n_boot = NULL,
    block_length = NULL,
    block_type = c("moving", "stationary", "hmm", "msar", "wild"),
    num_blocks = NULL,
    num_boots = 100L,
    func = mean,
    apply_func_to = c("cols", "df"),
    p_method = c("1/n", "plugin", "cross validation"),
    p = NULL,
    overlap = TRUE,
    ar_order = 1,
    num_states = 2,
    model_func = default_model_func,
    score_func = mse
) {
  
 
  
  ## --- Validation ---
  ## match.arg() picks the first element in the vector as default
  block_type <- match.arg(block_type)
  apply_func_to <- match.arg(apply_func_to)
  p_method <- match.arg(p_method)
  
  if(.invalid_length(x))
    stop("No valued x value provided.")
  if (is.vector(x) || is.data.frame(x) || is.ts(x)) x <- as.matrix(x)

  n <- nrow(x)
  
  if (!is.function(func))
    stop("`func` must be a valid function.")

  ## Note: If NULL, value is calculated automatically below
  if (.invalid_length(n_boot))
    stop("`n_boot` must be a positive integer or NULL.")
  if (.invalid_length(block_length))
    stop("`block_length` must be a positive integer or NULL.")
  if (.invalid_length(num_blocks))
    stop("`num_blocks` must be a positive integer or NULL.")

  ## Need to provide either n_boot or num_blocks.
  ## If missing, block_length will be calculated automatically.
  ## n_boot = block_length * num_blocks.
  if(.invalid_length(n_boot) && .invalid_length(num_blocks)) {
    stop("Must provide either n_boot or num_blocks.")
  } 
  if(.invalid_length(n_boot)) {
    if (is.null(block_length)) {block_length <- compute_default_block_length(x)}
    n_boot <- num_blocks * block_length
  }
  if(.invalid_length(num_blocks)) {
    if (is.null(block_length)) {block_length <- compute_default_block_length(x)}
    num_blocks <- n_boot / block_length
  }
  
    
  ## Fails if NULL. Value is not calculated automatically.
  if (.invalid_length(ar_order, allow_null = FALSE))
    stop("`ar_order` must be a positive integer.")
  if (.invalid_length(num_states, allow_null = FALSE))
    stop("`num_states` must be a positive integer.")
  if (.invalid_length(num_boots, allow_null = FALSE))
    stop("`num_boots` must be a positive integer.")
  
  ## Validate p
  if (.invalid_length(p) ) {
    stop("`p` must be a single number in (0,1) or NULL.")
  }
  if (!is.null(p) && ( length(p) > 1 || p <= 0 || p >= 1)) {
    stop("`p` must be a single number in (0,1) or NULL.")
  }

  bootstrap_series <- switch(
    block_type,
    moving = {
      blockBootstrap(x, n_boot, block_length, num_blocks, num_boots, "moving", p, overlap)
    },
    stationary = {
      if(is.null(p)) {p <- .estimate_p(x, p_method, block_length, model_func, score_func)}
      blockBootstrap(x, n_boot, block_length, num_blocks, num_boots, "stationary", p, overlap)
    },
    hmm = lapply(hmm_bootstrap(x[,1], num_states = num_states, num_blocks = num_blocks, num_boots = num_boots),
                 function(s) matrix(s, ncol = 1)),
    msar = lapply(msar_bootstrap(x[,1], ar_order = ar_order, num_states = num_states,
                                 num_blocks = num_blocks, num_boots = num_boots),
                  function(s) matrix(s, ncol = 1)),
    wild = wild_bootstrap(x, num_boots),
    stop("Unsupported block_type.")
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



#' Invalid Length
#'
#' @param x x, vector, matrix or data.frame
#' @param allow_null NULL value may be allowed when functions calculates value
#'   of x automatically.
#' 
#' @details
#' As long as exists() is before the first "or", this will return TRUE,
#' when x doesn't exist.
#'
#' @returns
#' @export
#'
#' @examples
.invalid_length <- function(x, allow_null = TRUE) {
  ## deparse(substitute(x)) converts variable name to string.
  
  !exists(deparse(substitute(x)), where = parent.frame()) || 
  if(allow_null) {
    !is.null(x) && (!is.numeric(x) || length(x) < 1) 
  } else {
    !is.numeric(x) || length(x) < 1
  } || 
  any(is.na(x)) || 
  any(!is.finite(x))
}




#' Estimate p
#'
#' @param x x
#' @param p_method p estimation method
#' @param block_length Block length
#' @param model_func Model function for k fold cross validation
#' @param score_func Score function for k fold cross validation
#'
#' @returns p
#'
#' @examples
.estimate_p <- function(x, p_method, block_length, model_func, score_func) {
  if (p_method == "1/n") {
    p <- 1 / if (is.null(block_length)) compute_default_block_length(x) else block_length
  } else if (p_method == "plugin") {
    ## acf[2,,1] is the autocorrelation for the 1st order (index 2) the first 
    ## variable (col 1).
    #ac1 <- acf(x, lag.max = 1, plot = FALSE)$acf[2,,1]
    
    ## acf[2,,] is the vector of autocorrelations for the 1st order (index 2) 
    ## of each variable (if multivariate).
    ac1 <- acf(x, lag.max = 1, plot = FALSE)$acf[2,,]
    
    ## The larger the average 1st order autocorrelation, the smaller the p
    p <- 1 - abs(mean(ac1, na.rm = TRUE))
  } else if (p_method == "cross validation") {
    if (is.null(model_func) || is.null(score_func))
      stop("For p_method = 'cross validation', provide `model_func` and `score_func`.")
    p <- k_fold_cv_ts(x, x[,1], model_func, score_func)
  }
}