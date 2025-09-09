test_that("tsbs returns list with model and score results", {
  skip_on_cran()
  x <- rnorm(50)
  result <- tsbs(x)
  expect_true(is.list(result))
  expect_named(result, c("model", "score"))
})

test_that("tsbs respects user-supplied model_func and score_func", {
  skip_on_cran()
  x <- rnorm(50)
  result <- tsbs(
    x,
    model_func = default_model_func,
    score_func = default_score_func
  )
  expect_true(is.list(result))
})

test_that("bootstrap returns list with num_boots elements", {
  skip_on_cran()
  x <- matrix(rnorm(100), ncol = 2)
  out <- tsbs(
    x, num_boots = 2, block_length = 5, func = mean,
    apply_func_to = "cols"
  )
  expect_length(out, 2)
  expect_true(all(sapply(out, is.numeric)))
})

test_that("bootstrap supports p_method options", {
  skip_on_cran()
  x <- matrix(rnorm(100), ncol = 2)
  expect_silent(
    tsbs(x, num_boots = 3, p_method = "percentile")
  )
  expect_error(
    bootstrap(x, num_boots = 3, p_method = "unknown"),
    regexp = "Unsupported"
  )
})

test_that("bootstrap stops on non-matrix x", {
  skip_on_cran()
  x <- data.frame(matrix(rnorm(100), ncol = 2))
  expect_error(
    tsbs(x, num_boots = 2), "Not a matrix"
  )
})
