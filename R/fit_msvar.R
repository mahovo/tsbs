#' Fit a 2-State Markov-Switching Vector Autoregressive Model
#'
#' This function fits a 2-state MS-VAR(1) model using a C++ implementation of 
#' the Expectation-Maximization (EM) algorithm. 
#' Note: Here VAR stands for Vector AutoRegressive (not "variance", not
#' "Value At Risk").
#'
#' @param y A numeric matrix or data frame where rows are observations and
#'   columns are the time series variables.
#' @param max_iter An integer specifying the maximum number of EM iterations.
#'   Defaults to 100.
#' @param tol A numeric value specifying the convergence tolerance for the
#'   log-likelihood. Defaults to 1e-6.
#'
#' @return A list containing the estimated model parameters:
#'   \item{beta1, beta2}{VAR coefficient matrices for each state.}
#'   \item{sigma1, sigma2}{Error covariance matrices for each state.}
#'   \item{P}{The 2x2 transition probability matrix.}
#'   \item{log_likelihood}{The final log-likelihood value.}
#'   \item{smoothed_probabilities}{A matrix of smoothed state probabilities.}
#'
#' @export
#' @examples
#' # Generate sample data
#' set.seed(123)
#' T_obs <- 250
#' y1 <- arima.sim(model = list(ar = 0.7), n = T_obs)
#' y2 <- 0.5 * y1 + arima.sim(model = list(ar = 0.3), n = T_obs)
#' sample_data <- cbind(y1, y2)
#'
#' # Fit the model (assuming the package is loaded)
#' # msvar_fit <- fit_msvar(sample_data)
#'
#' # View results
#' # print(msvar_fit$P)
#' # plot(msvar_fit$smoothed_probabilities[, 1], type = 'l',
#' #      main = "Smoothed Probability of State 1", ylab = "Probability")
fit_msvar <- function(y, max_iter = 100, tol = 1e-6) {
  
  # Ensure y is a matrix
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }
  
  # Call the C++ function directly
  # Rcpp creates the R bindings automatically during package compilation.
  results <- fit_msvar_em_cpp(y, max_iter, tol)
  
  # Add variable names to coefficients and matrices for better output
  var_names <- colnames(y)
  if (is.null(var_names)) {
    var_names <- paste0("y", 1:ncol(y))
  }
  
  # Assumes a VAR(1) model as in the C++ implementation
  lag_names <- paste0(rep(var_names, each = 1), "_lag1")
  param_names <- c("const", lag_names)
  
  rownames(results$beta1) <- rownames(results$beta2) <- param_names
  colnames(results$beta1) <- colnames(results$beta2) <- var_names
  rownames(results$sigma1) <- colnames(results$sigma1) <- var_names
  rownames(results$sigma2) <- colnames(results$sigma2) <- var_names
  colnames(results$smoothed_probabilities) <- c("State1_Prob", "State2_Prob")
  
  return(results)
}