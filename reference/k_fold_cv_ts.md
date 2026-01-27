# K-Fold Cross-Validation for Time Series

Time-series-aware k-fold cross-validation over a grid of candidate `p`
values, using a user-specified model and scoring function.

## Usage

``` r
k_fold_cv_ts(
  x,
  y,
  k = 5,
  p_grid = seq(0.01, 0.5, by = 0.01),
  model_func = default_model_func,
  score_func = mse
)
```

## Arguments

- x:

  Numeric matrix or vector of time series data (rows = time points,
  columns = variables).

- y:

  Numeric vector of target values to predict (must match number of rows
  of `x`).

- k:

  Integer number of folds.

- p_grid:

  Numeric vector of candidate `p` values to evaluate.

- model_func:

  A function that accepts training `x_train`, `y_train` and returns a
  fitted model object.

- score_func:

  A scoring function of the form `function(preds, y_val)` returning a
  numeric score.

## Value

Optimal `p` value that minimizes the average validation score.

## Details

This function partitions the time series into contiguous validation
folds, fits the model on the training set for each fold, and evaluates
the predictive performance on the validation set using `score_func`.

## Examples

``` r
# Generate toy time-series data
set.seed(123)
x <- matrix(rnorm(200), ncol = 2)
y <- 0.5 * x[,1] - 0.3 * x[,2] + rnorm(100)

# Perform time-series k-fold CV using default model and scoring functions
best_p <- k_fold_cv_ts(
  x = x,
  y = y,
  k = 5,
  p_grid = seq(0.01, 0.5, by = 0.01),
  model_func = default_model_func,
  score_func = mse
)
print(best_p)
#> [1] 0.01
```
