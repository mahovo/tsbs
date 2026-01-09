test_that("k_fold_cv_ts returns numeric optimal param", {
  x <- matrix(rnorm(100), ncol = 2)  # x should be matrix
  y <- rnorm(50)                      # y is required
  p_grid <- c(1, 2, 3)                # named p_grid, not p_list
  opt <- k_fold_cv_ts(
    x = x, 
    y = y, 
    k = 2, 
    p_grid = p_grid,
    model_func = default_model_func, 
    score_func = mse
  )
  expect_true(opt %in% p_grid)
})

test_that("k_fold_cv_ts errors if k is larger than nrow(x)", {
  x <- matrix(rnorm(20), ncol = 2)  # 10 rows
  y <- rnorm(10)
  expect_error(
    k_fold_cv_ts(x = x, y = y, p_grid = c(1, 2), k = 20),
    "Cannot split"
  )
})
