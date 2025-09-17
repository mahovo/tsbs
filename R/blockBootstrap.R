#' Block Bootstrap for Time Series (Rcpp Implementation)
#'
#' Generates bootstrap replicates of a multivariate time series using either the
#' Moving Block Bootstrap (MBB) or Stationary Bootstrap (SB) methods.
#' Not intended to be called directly by user, only to be called internally by 
#' tsbs().
#'
#' @param x A numeric matrix with rows representing time points and columns 
#'   representing variables.
#' @param n_length Integer or `NA`. Desired number of time points in the 
#'   bootstrap sample. If `NA`, the original series length is used.
#' @param block_length Integer or `NA`. Length of each block. If `NA`, a 
#'   heuristic is used based on the average lag-1 autocorrelation across variables.
#' @param num_blocks Integer or `NA`. Number of blocks to use. If `NA`, 
#'   this is inferred from `n_length`.
#' @param num_boots Number of bootstrap replicates to generate.
#' @param bs_type Bootstrap type. Character string: One of `"moving"`, 
#'   `"stationary"`, `"hmm"`, `"msar"`, `"ms_varma_garch"`, or `"wild"`. 
#' @param block_type Character. Block type. Either `"non-overlapping"`, 
#'   `"overlapping"` or `"tapered"`. Default `"overlapping"`.
#' @param taper_type Tapering window function. Character. One of `"cosine"`, 
#'   `"bartlett"`, or `"tukey"`.
#' @param p Probability parameter for the geometric block length (used in 
#'   Stationary Bootstrap).
#' @param stationary_max_percentile Stationary max percentile.
#' @param stationary_max_fraction_of_n Stationary max fraction of n.
#'
#' @return A list of matrices, each one a bootstrap replicate with dimensions 
#'   approximately `n_length` × `ncol(x)`.
#'
#' @details
#' The Moving Block Bootstrap (MBB) resamples fixed-length overlapping blocks
#' from the original time series. The Stationary Bootstrap (SB) samples
#' blocks with geometrically distributed lengths and can optionally enforce
#' non-overlap constraints.
#'
#' If `block_length_spec` is `NA`, an automatic heuristic is used that depends
#' on the average lag-1 autocorrelation across variables. The function
#' supports flexible specifications of output length (`n_length`) and
#' block structure (`num_blocks`).
#'
#' @examples
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol = 1)
#' boots <- blockBootstrap(x, type = "stationary", num_boots = 5)
#' dim(boots[[1]])
#'
#' @importFrom Rcpp sourceCpp
#' @export
blockBootstrap <- function(
  x,
  n_boot = NULL,
  block_length = NULL,
  bs_type = "moving",
  block_type = "overlapping",
  taper_type = "cosine",
  tukey_alpha = 0.5, 
  num_blocks = NULL,
  num_boots = 20L,
  p = 0.1,
  stationary_max_percentile = 0.99,
  stationary_max_fraction_of_n = 0.10
) {

  .Call(`_tsbs_blockBootstrap_cpp`,
    as.matrix(x),
    n_boot,
    block_length,
    as.character(bs_type),
    as.character(block_type),
    as.character(taper_type),
    tukey_alpha,
    num_blocks,
    as.integer(num_boots),
    p,
    stationary_max_percentile,
    stationary_max_fraction_of_n
  )
}


