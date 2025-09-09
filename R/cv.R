#' K-Fold Cross-Validation for Time Series
#'
#' Time-series-aware k-fold cross-validation over a grid of candidate `p` 
#'   values, using a user-specified model and scoring function.
#'
#' @param x Numeric matrix or vector of time series data (rows = time points, 
#'   columns = variables).
#' @param y Numeric vector of target values to predict (must match number of 
#'   rows of `x`).
#' @param model_func A function that accepts training `x_train`, `y_train` and 
#'   returns a fitted model object.
#' @param score_func A scoring function of the form `function(preds, y_val)` 
#'   returning a numeric score.
#' @param k Integer number of folds.
#' @param p_grid Numeric vector of candidate `p` values to evaluate.
#'
#' @return Optimal `p` value that minimizes the average validation score.
#'
#' @details
#' This function partitions the time series into contiguous validation folds,
#'   fits the model on the training set for each fold, and evaluates the 
#'   predictive performance on the validation set using `score_func`.
#'
#' @examples
#' # Generate toy time-series data
#' set.seed(123)
#' x <- matrix(rnorm(200), ncol = 2)
#' y <- 0.5 * x[,1] - 0.3 * x[,2] + rnorm(100)
#'
#' # Perform time-series k-fold CV using default model and scoring functions
#' best_p <- k_fold_cv_ts(
#'   x = x,
#'   y = y,
#'   p_list = c(1, 5, 10),
#'   k = 5,
#'   model_func = default_model_func,
#'   score_func = default_score_func
#' )
#' print(best_p)
#'
#' @export
k_fold_cv_ts <- function(
    x, 
    y, 
    k = 5, 
    p_grid = seq(0.01, 0.5, by = 0.01), 
    model_func = default_model_func, 
    score_func = mse) {
  # Input validation
  if (!is.numeric(x)) stop("`x` must be numeric vector, matrix, or data frame.")
  if (is.vector(x)) x <- as.matrix(x)
  if (nrow(x) != length(y)) stop("Number of rows in `x` must match length of `y`.")
  if (!is.function(model_func)) stop("`model_func` must be a valid function.")
  if (!is.function(score_func)) stop("`score_func` must be a valid function.")
  
  n <- nrow(x)
  fold_size <- floor(n / k)
  ends <- seq(fold_size, by = fold_size, length.out = k)
  starts <- c(1, head(ends, -1) + 1)
  
  # Evaluate each p value
  mean_scores <- sapply(p_grid, function(p_val) {
    fold_scores <- mapply(
      function(start, end) {
        val_idx <- seq(from = start, to = end)
        train_idx <- setdiff(seq_len(n), val_idx)
        
        x_train <- x[train_idx, , drop = FALSE]
        y_train <- y[train_idx]
        x_val <- x[val_idx, , drop = FALSE]
        y_val <- y[val_idx]
        
        model <- model_func(x_train, y_train)
        preds <- predict(model, newdata = as.data.frame(x_val))  # assumes model supports predict
        score_func(preds, y_val)
      },
      starts, ends
    )
    mean(unlist(fold_scores))
  })
  
  # Choose the p that minimizes average validation error
  best_p <- p_grid[which.min(mean_scores)]
  best_p
}

