#' Block Bootstrap for Time Series (Rcpp Implementation)
#'
#' Generates bootstrap replicates of a multivariate time series using either the
#' Moving Block Bootstrap (MBB) or Stationary Bootstrap (SB) methods.
#'
#' @name blockBootstrap
#' @rdname blockBootstrap
#'
#' @param x A numeric matrix with rows representing time points and columns representing variables.
#' @param n_length Integer. Desired number of time points in the bootstrap sample. If `-1`, uses the original series length.
#' @param block_length Integer. Length of each block. If `-1`, a default is computed based on the lag-1 autocorrelation of the series.
#' @param num_blocks Optional number of blocks to use instead of setting `n_length`.
#' @param num_boots Number of bootstrap replicates to generate.
#' @param block_type Character. Type of block bootstrap: either `"moving"` (Moving Block Bootstrap) or `"stationary"` (Stationary Bootstrap).
#' @param p Probability parameter for the geometric distribution (used in Stationary Bootstrap).
#' @param overlap Logical. If `FALSE`, stationary blocks are chosen from non-overlapping segments.
#'
#' @return A list of matrices, each one a bootstrap replicate with dimensions approximately `n_length × ncol(x)`.
#'
#' @details
#' The Moving Block Bootstrap (MBB) resamples fixed-length overlapping blocks from the original time series.
#' The Stationary Bootstrap (SB) samples blocks with geometrically distributed lengths and can optionally avoid overlap.
#'
#' If `block_length` is set to `-1`, an automatic heuristic is used that depends on the average lag-1 autocorrelation across variables.
#'
#' @examples
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol = 1)
#' boots <- blockBootstrap(x, block_type = "stationary", num_boots = 5)
#' dim(boots[[1]])
#'
#' @useDynLib tsbs, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
NULL


#' Compute Default Block Length Based on Autocorrelation (Internal Utility)
#'
#' Computes a heuristic default block length for a multivariate time series
#' based on the average lag-1 autocorrelation.
#'
#' @name compute_default_block_length
#' @rdname compute_default_block_length
#'
#' @param x A numeric matrix with rows as time points and columns as variables.
#'
#' @return An integer block length.
#'
#' @examples
#' x <- matrix(rnorm(100), ncol = 2)
#' compute_default_block_length(x)
#'
#' @useDynLib tsbs, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
NULL