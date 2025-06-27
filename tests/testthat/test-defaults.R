test_that("default_model_func returns numeric vector of correct length", {
  x <- rnorm(50)
  fit <- default_model_func(x, param = list(p = 2))
  expect_true(is.numeric(fit))
  expect_equal(length(fit), length(x))
})

test_that("default_score_func returns numeric scalar", {
  x <- rnorm(50); fit <- x + rnorm(50)
  expect_true(is.numeric(default_score_func(x, fit)))
  expect_length(default_score_func(x, fit), 1)
})
