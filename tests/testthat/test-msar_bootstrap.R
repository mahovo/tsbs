test_that("msar_bootstrap returns list of correct length and shapes", {
  skip_if_not_installed("MSwM")
  df <- data.frame(return = rnorm(50))
  boots <- msar_bootstrap(df, num_boots = 2)
  expect_length(boots, 2)
  expect_true(all(sapply(boots, is.data.frame)))
})

test_that("msar_bootstrap fails gracefully if MSwM not installed", {
  if (!requireNamespace("MSwM", quietly = TRUE)) {
    df <- data.frame(return = rnorm(50))
    expect_error(
      msar_bootstrap(df), "MSwM"
    )
  }
})
