#' tsbs: Advanced Bootstrap for Time Series
#'
#' Generates bootstrap replicates of multivariate or (multiple) univariate time
#' series. Supports moving and stationary block, HMM, MSVAR, MS VARMA GARCH and 
#' wild bootstrap types.
#'
#' @section Main Function:
#' \itemize{
#'   \item \code{\link{tsbs}}: Main user-facing function
#' }
#'
#'
#' @seealso
#' Vignettes:
#' \itemize{
#'   \item \code{vignette("Diagnostics", package = "tsbs")} - User guide for 
#'     MS-VARMA-GARCH model fitting diagnostic system 
#' }
#'
#' @import tsgarch
#' @import tsmarch
#' @import tsmethods
#' @import ggplot2
#' @import data.table
#' @importFrom stats acf ar
#' @importFrom tsmethods tsfilter
#' @importFrom Rcpp sourceCpp
#' @keywords internal
"_PACKAGE"

