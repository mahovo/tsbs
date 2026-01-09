test_that("default_model_func returns model object", {
  x <- matrix(rnorm(50), ncol = 1)
  y <- rnorm(50)
  model <- default_model_func(x, y)
  expect_s3_class(model, "default_model")
  expect_true(is.numeric(model$mean_y))
  expect_length(model$mean_y, 1)
  
  # Test predict method
  preds <- predict(model, newdata = as.data.frame(x))
  expect_true(is.numeric(preds))
  expect_equal(length(preds), nrow(x))
})

test_that("mse returns numeric scalar", {
  actual <- rnorm(50)
  pred <- actual + rnorm(50, sd = 0.1)
  result <- mse(pred, actual)
  expect_true(is.numeric(result))
  expect_length(result, 1)
  expect_true(result >= 0)  # MSE is always non-negative
})
