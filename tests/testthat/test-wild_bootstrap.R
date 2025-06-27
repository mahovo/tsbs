test_that("wild_bootstrap returns list of length num_boots", {
  x <- rnorm(100)
  out <- wild_bootstrap(x, num_boots = 2, wild_type = "mammen")
  expect_length(out, 2)
  expect_true(all(sapply(out, is.numeric)))
})

test_that("wild_bootstrap throws error on invalid wild_type", {
  expect_error(
    wild_bootstrap(rnorm(50), num_boots = 2, wild_type = "unknown")
  )
})
