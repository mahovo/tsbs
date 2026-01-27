# Default Model Fitting Function

Fits a simple mean model that returns the mean of `y_train`. The
[`predict()`](https://rdrr.io/r/stats/predict.html) method returns this
mean for all new observations.

## Usage

``` r
default_model_func(x_train, y_train)
```

## Arguments

- x_train:

  Matrix or data frame of training features.

- y_train:

  Numeric vector of training targets.

## Value

A list with class `"default_model"` containing the mean of the training
target.

## Examples

``` r
default_model_func(matrix(rnorm(100), ncol = 1), rnorm(100))
#> $mean_y
#> [1] 0.1482784
#> 
#> attr(,"class")
#> [1] "default_model"
```
