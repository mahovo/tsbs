#' Block Bootstrap for Time Series (Rcpp Implementation)
#'
#' Generates bootstrap replicates of a multivariate time series using either the
#' Moving Block Bootstrap (MBB) or Stationary Bootstrap (SB) methods.
#' Not intended to be called directly by user, only to be called internally by 
#' tsbs().
#'
#' @param x A numeric matrix with rows representing time points and columns
#'   representing variables.
#' @param n_boot Integer or `NULL`. Desired number of time points in the
#'   bootstrap sample. If `NULL`, the original series length is used.
#' @param block_length Integer or `NULL`. Length of each block. If `NULL`, a
#'   heuristic is used based on the average lag-1 autocorrelation across
#'   variables.
#' @param bs_type Bootstrap type. Character string: One of `"moving"` or
#'   `"stationary"`. Default `"moving"`.
#' @param block_type Character. Block type. Either `"non-overlapping"`,
#'   `"overlapping"` or `"tapered"`. Default `"overlapping"`.
#' @param taper_type Tapering window function. Character. One of `"cosine"`,
#'   `"bartlett"`, or `"tukey"`. Default `"cosine"`.
#' @param tukey_alpha `alpha` when `block_type = "tapered"` and  `taper_type = "tukey"`.
#'   Default 0.5.
#' @param num_blocks Integer or `NULL`. Number of blocks to use. If `NULL`, this 
#'   is inferred from `n_boot`.
#' @param num_boots Number of bootstrap replicates to generate.
#' @param p numeric \eqn{p \in (0, 1)}. Probability parameter for the geometric
#'   block length (used in Stationary Bootstrap). Default is 0.1.
#' @param stationary_max_percentile Stationary max percentile. Default is 0.99.
#' @param stationary_max_fraction_of_n Stationary max fraction of n. Default is
#'   0.10.
#'
#' @return A list of matrices, each one a bootstrap replicate with dimensions
#'   approximately `n_length` × `ncol(x)`.
#'
#' @details The Moving Block Bootstrap (MBB) resamples fixed-length overlapping
#' blocks from the original time series. The Stationary Bootstrap (SB) samples
#' blocks with geometrically distributed lengths and can optionally enforce
#' non-overlap constraints.
#' 
#' If `n_boot` is set, the last block will be trimmed when necessary.
#' If `n_boot` is not set, and `num_blocks` is set, the length of each 
#' bootstrap series will be determined by the number of blocks and the 
#' lengths of the individual blocks for that particular series.
#' If neither `n_boot` nor `num_blocks` is set, `n_boot` will default to the
#' number of rows in `x` and the last block will be trimmed when necessary.
#' 
#' @references 
#' Politis, D., & Romano, J. (1994). The Stationary Bootstrap. Journal of the 
#'   American Statistical Association, 89, 1303-1313.
#'   [http://dx.doi.org/10.1080/01621459.1994.10476870](http://dx.doi.org/10.1080/01621459.1994.10476870)
#'
#' @examples
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol = 1)
#' boots <- blockBootstrap(x, bs_type = "stationary", num_boots = 5L)
#' dim(boots[[1]])
#'
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


