#' Mean Squared Error (MSE)
#'
#' @param actual Actual values.
#' @param predicted Predicted values.
#' 
#' @return MSE
#' 
#' @export
mse <- function(pred, actual) {
  mean((pred - actual)^2)
}

#' Mean Absolute Error (MAE)
#'
#' @param actual Actual values.
#' @param predicted Predicted values.
#' 
#' @return MAE
#' 
#' @export
mae <- function(pred, actual) {
  mean(abs(pred - actual))
}
