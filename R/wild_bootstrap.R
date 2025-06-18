#' Wild Bootstrap for Residual Resampling
#'
#' Performs wild bootstrap by multiplying residuals by random weights (+1 or -1).
#'
#' @param x Numeric vector or matrix.
#' @param num_boots Number of bootstrap replicates.
#'
#' @return A list of bootstrapped samples.
#'
#' @examples
#' set.seed(789)
#' resids <- rnorm(100)
#' boot_samples <- wild_bootstrap(resids, num_boots = 50)
#'
#' @export
wild_bootstrap <- function(x, num_boots) {
  if (is.vector(x)) x <- matrix(x, ncol = 1)
  
  n <- nrow(x)
  d <- ncol(x)
  bootstrap_samples <- vector("list", num_boots)
  
  for (b in seq_len(num_boots)) {
    v <- sample(c(-1, 1), size = n, replace = TRUE)
    bootstrap_samples[[b]] <- x * v
  }
  bootstrap_samples
}
