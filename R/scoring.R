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

#' K-Fold Cross-Validation for Time Series Models
#'
#' Performs k-fold cross-validation while preserving temporal structure.
#'
#' @param x Numeric matrix or time series.
#' @param model_func A function taking a training set and returning a prediction function.
#' @param k Number of folds.
#' @param score_func A function to score prediction accuracy (e.g., MSE, MAE).
#'
#' @return Cross-validation score.
#' 
#' @importFrom dplyr bind_rows
#' @export
k_fold_cv_ts <- function(x, y, k = 5, p_grid = seq(0.01, 0.5, by = 0.01), model_func, score_func) {
  n <- nrow(x)
  fold_size <- floor(n / k)
  ends <- seq(fold_size, by = fold_size, length.out = k)
  starts <- c(1, head(ends, -1) + 1)
  
  results <- sapply(p_grid, function(p_val) {
    scores <- mapply(function(start, end) {
      val_idx <- start:end
      train_idx <- setdiff(1:n, val_idx)
      
      x_train <- x[train_idx, , drop = FALSE]
      y_train <- y[train_idx]
      x_val   <- x[val_idx, , drop = FALSE]
      y_val   <- y[val_idx]
      
      model <- model_func(x_train, y_train)
      preds <- predict_ar_model(model, newdata = as.data.frame(x_val))
      score_func(preds, y_val)
    }, starts, ends)
    mean(unlist(scores))
  })
  
  best_p <- p_grid[which.min(results)]
  return(best_p)
}
