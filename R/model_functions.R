#' Default Model Fitting Function
#'
#' Fits a simple mean model that returns the mean of `y_train`. The `predict()` method returns this mean for all new observations.
#'
#' @param x_train Matrix or data frame of training features.
#' @param y_train Numeric vector of training targets.
#'
#' @return A list with class `"default_model"` containing the mean of the training target.
#' @examples
#' default_model_func(matrix(rnorm(100), ncol = 1), rnorm(100))
#' @export
default_model_func <- function(x_train, y_train) {
  structure(
    list(mean_y = mean(y_train)),
    class = "default_model"
  )
}

#' @export
predict.default_model <- function(object, newdata, ...) {
  rep(object$mean_y, nrow(newdata))
}
