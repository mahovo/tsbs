test_that("blockBootstrap returns list of correct length and shapes", {
  set.seed(123)
  x <- matrix(rnorm(100), ncol = 2)
  boots <- blockBootstrap(x, num_boots = 5, block_length = 10, block_type = "moving")
  
  expect_type(boots, "list")
  expect_length(boots, 5)
  for (b in boots) {
    expect_true(is.matrix(b))
    expect_equal(ncol(b), 2)
  }
})


test_that("blockBootstrap supports stationary type and optional shuffling", {
  set.seed(123)
  x <- matrix(rnorm(100), ncol = 2)
  boots <- blockBootstrap(
    x, num_boots = 3, block_length = 5,
    block_type = "stationary", overlap = TRUE
  )
  expect_length(boots, 3)
})


test_that("blockBootstrap errors on bad block_type", {
  x <- matrix(rnorm(100), ncol = 2)
  expect_error(blockBootstrap(x, block_type = "invalid"))
})


test_that("blockBootstrap computes default block_length when -1 passed", {
  x <- matrix(rnorm(100), ncol = 2)
  boots <- blockBootstrap(x, block_length = NA, num_boots = 2)
  expect_length(boots, 2)
  expect_true(nrow(boots[[1]]) == nrow(x))
})


test_that("blockBootstrap works with defaults", {
  set.seed(123)
  x <- matrix(rnorm(100), ncol = 1)
  result <- blockBootstrap(x,
                           n_length = NA,
                           block_length = NA,
                           num_blocks = NA,
                           num_boots = 2,
                           block_type = "moving",
                           p = 0.1,
                           overlap = TRUE)
  
  expect_type(result, "list")
  expect_length(result, 2)
  expect_equal(nrow(result[[1]]), nrow(x))
})


test_that("blockBootstrap works with explicit block_length", {
  x <- matrix(rnorm(100), ncol = 1)
  result <- blockBootstrap(x,
                           n_length = NA,
                           block_length = 10,
                           num_blocks = NA,
                           num_boots = 1,
                           block_type = "moving",
                           p = 0.1,
                           overlap = TRUE)
  
  expect_equal(nrow(result[[1]]), nrow(x))
})


test_that("blockBootstrap works with explicit num_blocks", {
  x <- matrix(rnorm(100), ncol = 1)
  result <- blockBootstrap(x,
                           n_length = NA,
                           block_length = 10,
                           num_blocks = 3,
                           num_boots = 1,
                           block_type = "moving",
                           p = 0.1,
                           overlap = TRUE)
  
  expect_equal(nrow(result[[1]]), 30) # 3 blocks * block_length (10) = 30
})


test_that("blockBootstrap works with stationary block_type", {
  x <- matrix(rnorm(100), ncol = 1)
  result <- blockBootstrap(x,
                           n_length = NA,
                           block_length = 10,
                           num_blocks = NA,
                           num_boots = 1,
                           block_type = "stationary",
                           p = 0.1,
                           overlap = FALSE)
  
  expect_equal(nrow(result[[1]]), 100) # Default length
})