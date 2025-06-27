#' Block Bootstrap for Time Series (Rcpp Implementation)
#'
#' Generates bootstrap replicates of a multivariate time series using either the
#' Moving Block Bootstrap (MBB) or Stationary Bootstrap (SB) methods.
#'
#'
#' @param x A numeric matrix with rows representing time points and columns representing variables.
#' @param n_length_spec Integer or `NA`. Desired number of time points in the bootstrap sample.
#'        If `NA`, the original series length is used.
#' @param block_length_spec Integer or `NA`. Length of each block.
#'        If `NA`, a heuristic is used based on the average lag-1 autocorrelation across variables.
#' @param num_blocks_spec Integer or `NA`. Number of blocks to use.
#'        If `NA`, this is inferred from `n_length_spec`.
#' @param num_boots Number of bootstrap replicates to generate.
#' @param block_type Character. Either `"moving"` (Moving Block Bootstrap) or
#'        `"stationary"` (Stationary Bootstrap).
#' @param p Probability parameter for the geometric block length (used in Stationary Bootstrap).
#' @param overlap Logical. If `FALSE`, stationary blocks are chosen from
#'        non-overlapping segments.
#'
#' @return A list of matrices, each one a bootstrap replicate with dimensions approximately
#' `n_length_spec` × `ncol(x)`.
#'
#' @details
#' The Moving Block Bootstrap (MBB) resamples fixed-length overlapping blocks
#' from the original time series. The Stationary Bootstrap (SB) samples
#' blocks with geometrically distributed lengths and can optionally enforce
#' non-overlap constraints.
#'
#' If `block_length_spec` is `NA`, an automatic heuristic is used that depends
#' on the average lag-1 autocorrelation across variables. The function
#' supports flexible specifications of output length (`n_length_spec`) and
#' block structure (`num_blocks_spec`).
#'
#' @examples
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol = 1)
#' boots <- blockBootstrap(x, block_type = "stationary", num_boots = 5)
#' dim(boots[[1]])
#'
#' @importFrom Rcpp sourceCpp
#' @export
blockBootstrap <- function(x,
                           n_length = NA_integer_,
                           block_length = NA_integer_,
                           num_blocks = NA_integer_,
                           num_boots = 20L,
                           block_type = "stationary",
                           p = 0.1,
                           overlap = TRUE) {

  if (length(p) != 1 || is.na(p) || !is.finite(p) || p <= 0 || p >= 1) {
    stop("'p' must be a single number in (0,1). Got: ", p)
  }
  
  .Call(`_tsbs_blockBootstrap_cpp`,
        as.matrix(x),
        as.numeric(n_length),
        as.numeric(block_length),
        as.numeric(num_blocks),
        as.integer(num_boots),
        as.character(block_type),
        as.numeric(p),
        as.logical(overlap))
}


