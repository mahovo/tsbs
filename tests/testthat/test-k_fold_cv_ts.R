test_that("k_fold_cv_ts returns numeric optimal param", {
  x <- rnorm(50)
  p_list <- c(1, 2, 3)
  opt <- k_fold_cv_ts(
    x, p_list, k = 2, model_func = default_model_func, score_func = mse
  )
  expect_true(opt %in% p_list)
})

test_that("k_fold_cv_ts errors if k is larger than length(x)", {
  x <- rnorm(10)
  expect_error(
    k_fold_cv_ts(x, p_list = c(1,2), k = 20),
    "cannot split into"
  )
})
