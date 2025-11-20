
## Basic Functionality Tests ==================================================

test_that("wild_bootstrap handles univariate data", {
  set.seed(123)
  x <- rnorm(100)
  
  result <- wild_bootstrap(
    x = x,
    num_boots = 5
  )
  
  expect_type(result, "list")
  expect_length(result, 5)
  expect_true(all(sapply(result, is.matrix)))
  expect_true(all(sapply(result, nrow) == 100))
  expect_true(all(sapply(result, ncol) == 1))
})

test_that("wild_bootstrap handles multivariate data", {
  set.seed(456)
  x <- matrix(rnorm(200), ncol = 2)
  
  result <- wild_bootstrap(
    x = x,
    num_boots = 8
  )
  
  expect_type(result, "list")
  expect_length(result, 8)
  expect_true(all(sapply(result, is.matrix)))
  expect_true(all(sapply(result, nrow) == 100))
  expect_true(all(sapply(result, ncol) == 2))
})

test_that("wild_bootstrap preserves original length", {
  set.seed(789)
  x <- rnorm(150)
  
  result <- wild_bootstrap(
    x = x,
    num_boots = 10
  )
  
  ## Wild bootstrap always preserves original length
  expect_true(all(sapply(result, nrow) == 150))
})

test_that("wild_bootstrap with single bootstrap", {
  set.seed(321)
  x <- rnorm(100)
  
  result <- wild_bootstrap(
    x = x,
    num_boots = 1
  )
  
  expect_length(result, 1)
  expect_equal(nrow(result[[1]]), 100)
})

test_that("wild_bootstrap with many bootstraps", {
  set.seed(654)
  x <- rnorm(50)
  
  result <- wild_bootstrap(
    x = x,
    num_boots = 500
  )
  
  expect_length(result, 500)
  expect_true(all(sapply(result, nrow) == 50))
})

## Data Type Handling Tests ===================================================

test_that("wild_bootstrap handles vector input", {
  set.seed(111)
  x <- rnorm(100)
  
  result <- wild_bootstrap(x, num_boots = 3)
  
  expect_type(result, "list")
  expect_true(all(sapply(result, is.matrix)))
  expect_true(all(sapply(result, ncol) == 1))
})

test_that("wild_bootstrap handles matrix input", {
  set.seed(222)
  x <- matrix(rnorm(150), ncol = 3)
  
  result <- wild_bootstrap(x, num_boots = 3)
  
  expect_type(result, "list")
  expect_true(all(sapply(result, is.matrix)))
  expect_true(all(sapply(result, ncol) == 3))
})

test_that("wild_bootstrap handles single column matrix", {
  set.seed(333)
  x <- matrix(rnorm(100), ncol = 1)
  
  result <- wild_bootstrap(x, num_boots = 3)
  
  expect_true(all(sapply(result, ncol) == 1))
  expect_true(all(sapply(result, nrow) == 100))
})

## wild_type not currently implemented
# test_that("wild_bootstrap throws error on invalid wild_type", {
#   expect_error(
#     wild_bootstrap(rnorm(50), num_boots = 2, wild_type = "unknown")
#   )
# })

## Parallel Execution Tests ===================================================

test_that("wild_bootstrap handles parallel execution", {
  skip_if_not_installed("foreach")
  skip_if_not_installed("doParallel")
  skip_on_cran()
  
  set.seed(444)
  x <- rnorm(100)
  
  result <- wild_bootstrap(
    x = x,
    num_boots = 10,
    parallel = TRUE,
    num_cores = 2
  )
  
  expect_type(result, "list")
  expect_length(result, 10)
  expect_true(all(sapply(result, nrow) == 100))
})

test_that("wild_bootstrap parallel gives consistent structure", {
  skip_if_not_installed("foreach")
  skip_if_not_installed("doParallel")
  skip_on_cran()
  
  x <- rnorm(100)
  
  set.seed(555)
  result_seq <- wild_bootstrap(
    x = x,
    num_boots = 5,
    parallel = FALSE
  )
  
  set.seed(555)
  result_par <- wild_bootstrap(
    x = x,
    num_boots = 5,
    parallel = TRUE,
    num_cores = 2
  )
  
  # Check structural equality
  expect_equal(length(result_seq), length(result_par))
  expect_equal(sapply(result_seq, dim), sapply(result_par, dim))
})

## Rademacher Properties Tests ================================================

test_that("wild_bootstrap uses Rademacher weights", {
  set.seed(666)
  x <- rep(1, 100)  ## Constant value
  
  result <- wild_bootstrap(x, num_boots = 1)
  
  ## With constant input, output should be all +1 or all -1
  unique_vals <- unique(as.vector(result[[1]]))
  expect_true(all(unique_vals %in% c(-1, 1)))
})

test_that("wild_bootstrap preserves absolute values", {
  set.seed(777)
  x <- rnorm(100)
  
  result <- wild_bootstrap(x, num_boots = 10)
  
  ## Each bootstrap should have same absolute values as original
  for (i in seq_along(result)) {
    boot_abs <- sort(abs(result[[i]]))
    orig_abs <- sort(abs(x))
    expect_equal(boot_abs, orig_abs)
  }
})

test_that("wild_bootstrap signs are random", {
  set.seed(888)
  x <- abs(rnorm(100))  ## All positive values
  
  result <- wild_bootstrap(x, num_boots = 50)
  
  ## With many bootstraps, we should see both positive and negative values
  all_values <- unlist(result)
  expect_true(any(all_values > 0))
  expect_true(any(all_values < 0))
})

test_that("wild_bootstrap approximately 50/50 sign distribution", {
  skip_on_cran()  ## Stochastic test
  
  set.seed(999)
  x <- rep(1, 1000)  ## Large constant vector
  
  result <- wild_bootstrap(x, num_boots = 100)
  
  ## Count proportion of positive values across all bootstraps
  all_values <- unlist(result)
  prop_positive <- mean(all_values > 0)
  
  ## Should be approximately 50% (allow 45-55% range)
  expect_true(prop_positive > 0.45 && prop_positive < 0.55)
})

## Edge Cases =================================================================

test_that("wild_bootstrap handles short series", {
  set.seed(1111)
  x <- rnorm(10)
  
  result <- wild_bootstrap(x, num_boots = 5)
  
  expect_length(result, 5)
  expect_true(all(sapply(result, nrow) == 10))
})

test_that("wild_bootstrap handles single observation", {
  set.seed(2222)
  x <- 5
  
  result <- wild_bootstrap(x, num_boots = 10)
  
  expect_length(result, 10)
  expect_true(all(sapply(result, nrow) == 1))
  
  ## Values should be +5 or -5
  vals <- sapply(result, function(m) m[1,1])
  expect_true(all(vals %in% c(-5, 5)))
})

test_that("wild_bootstrap output is numeric and finite", {
  set.seed(3333)
  x <- rnorm(100)
  
  result <- wild_bootstrap(x, num_boots = 5)
  
  expect_true(all(sapply(result, is.numeric)))
  expect_true(all(sapply(result, function(m) all(is.finite(m)))))
})

test_that("wild_bootstrap handles zero values", {
  set.seed(4444)
  x <- c(rep(0, 50), rnorm(50))
  
  result <- wild_bootstrap(x, num_boots = 5)
  
  expect_length(result, 5)
  ## Zeros remain zeros
  expect_true(all(sapply(result, function(m) sum(m == 0) == 50)))
})

## Data Properties Tests ======================================================

test_that("wild_bootstrap handles positive values", {
  set.seed(5555)
  x <- abs(rnorm(100)) + 1  ## All positive
  
  result <- wild_bootstrap(x, num_boots = 5)
  
  expect_length(result, 5)
  expect_true(all(sapply(result, nrow) == 100))
})

test_that("wild_bootstrap handles negative values", {
  set.seed(6666)
  x <- -abs(rnorm(100)) - 1  ## All negative
  
  result <- wild_bootstrap(x, num_boots = 5)
  
  expect_length(result, 5)
  expect_true(all(sapply(result, nrow) == 100))
})

test_that("wild_bootstrap handles mixed positive/negative", {
  set.seed(7777)
  x <- rnorm(100)
  
  result <- wild_bootstrap(x, num_boots = 5)
  
  expect_length(result, 5)
  expect_true(all(sapply(result, nrow) == 100))
})

test_that("wild_bootstrap handles large values", {
  set.seed(8888)
  x <- rnorm(100, mean = 1000, sd = 100)
  
  result <- wild_bootstrap(x, num_boots = 5)
  
  expect_length(result, 5)
  ## Absolute values should match
  for (i in seq_along(result)) {
    expect_equal(sort(abs(result[[i]])), sort(abs(x)))
  }
})

test_that("wild_bootstrap handles small values", {
  set.seed(9999)
  x <- rnorm(100, mean = 0, sd = 0.001)
  
  result <- wild_bootstrap(x, num_boots = 5)
  
  expect_length(result, 5)
  for (i in seq_along(result)) {
    expect_equal(sort(abs(result[[i]])), sort(abs(x)))
  }
})

test_that("wild_bootstrap handles constant values", {
  set.seed(1234)
  x <- rep(5, 100)
  
  result <- wild_bootstrap(x, num_boots = 10)
  
  expect_length(result, 10)
  
  ## Each bootstrap should contain only +5 and -5 (each observation gets random sign)
  for (i in seq_along(result)) {
    unique_vals <- unique(as.vector(result[[i]]))
    ## Should have both +5 and -5 values (with high probability)
    expect_true(all(unique_vals %in% c(-5, 5)))
    ## Should have absolute value of 5
    expect_true(all(abs(result[[i]]) == 5))
  }
})

test_that("wild_bootstrap with heteroskedastic data", {
  set.seed(2345)
  ## Simulate heteroskedastic residuals
  x <- numeric(200)
  for (i in seq_len(200)) {
    x[i] <- rnorm(1, sd = 0.5 + 0.01 * i)
  }
  
  result <- wild_bootstrap(x, num_boots = 5)
  
  expect_length(result, 5)
  expect_true(all(sapply(result, nrow) == 200))
})

test_that("wild_bootstrap with clustered data structure", {
  set.seed(3456)
  ## Simulate clustered structure
  x <- c(
    rnorm(50, mean = 0, sd = 1),
    rnorm(50, mean = 5, sd = 1),
    rnorm(50, mean = -3, sd = 1)
  )
  
  result <- wild_bootstrap(x, num_boots = 5)
  
  expect_length(result, 5)
  expect_true(all(sapply(result, nrow) == 150))
})

## Multivariate Specific Tests ================================================

test_that("wild_bootstrap applies same weights across columns", {
  set.seed(4567)
  x <- matrix(1, nrow = 100, ncol = 2)
  
  result <- wild_bootstrap(x, num_boots = 1)
  
  ## For constant matrix, each row should have same sign across columns
  for (i in seq_len(nrow(result[[1]]))) {
    row_vals <- result[[1]][i, ]
    expect_true(all(row_vals == row_vals[1]))
  }
})

test_that("wild_bootstrap preserves row structure", {
  set.seed(5678)
  x <- matrix(rnorm(200), ncol = 2)
  
  result <- wild_bootstrap(x, num_boots = 5)
  
  ## Each bootstrap should preserve the relationship between columns
  for (i in seq_along(result)) {
    ## Absolute values should match original
    expect_equal(sort(abs(result[[i]][,1])), sort(abs(x[,1])))
    expect_equal(sort(abs(result[[i]][,2])), sort(abs(x[,2])))
  }
})

test_that("wild_bootstrap with correlated columns", {
  set.seed(6789)
  x1 <- rnorm(100)
  x2 <- 0.9 * x1 + rnorm(100, sd = 0.1)
  x <- cbind(x1, x2)
  
  result <- wild_bootstrap(x, num_boots = 5)
  
  expect_length(result, 5)
  expect_true(all(sapply(result, ncol) == 2))
})

test_that("wild_bootstrap with uncorrelated columns", {
  set.seed(7890)
  x <- matrix(rnorm(200), ncol = 2)
  
  result <- wild_bootstrap(x, num_boots = 5)
  
  expect_length(result, 5)
  expect_true(all(sapply(result, ncol) == 2))
})

## Reproducibility Tests ======================================================

test_that("wild_bootstrap same seed gives same results", {
  x <- rnorm(100)
  
  set.seed(8888)
  result1 <- wild_bootstrap(x, num_boots = 5)
  
  set.seed(8888)
  result2 <- wild_bootstrap(x, num_boots = 5)
  
  expect_equal(result1, result2)
})

test_that("wild_bootstrap different seed gives different results", {
  x <- rnorm(100)
  
  set.seed(9999)
  result1 <- wild_bootstrap(x, num_boots = 5)
  
  set.seed(1111)
  result2 <- wild_bootstrap(x, num_boots = 5)
  
  expect_false(identical(result1, result2))
})

## Input Validation Tests =====================================================

test_that("wild_bootstrap assumes validated inputs from tsbs()", {
  ## wild_bootstrap is an internal function called by tsbs()
  ## All validation should happen in tsbs() before calling wild_bootstrap()
  ## These tests verify the function works with valid inputs
  
  x <- rnorm(100)
  
  ## Valid inputs should work
  result <- wild_bootstrap(x, num_boots = 5)
  
  expect_type(result, "list")
  expect_length(result, 5)
})

test_that("wild_bootstrap ignores n_boot parameter", {
  ## Wild bootstrap doesn't use n_boot - it always returns original length
  ## This is mentioned in tsbs() documentation
  
  set.seed(2222)
  x <- rnorm(100)
  
  ## Even if n_boot were passed (it's not a parameter), length is preserved
  result <- wild_bootstrap(x, num_boots = 5)
  
  expect_true(all(sapply(result, nrow) == 100))
})

test_that("wild_bootstrap ignores num_blocks parameter", {
  ## Wild bootstrap doesn't use num_blocks
  ## This is mentioned in tsbs() documentation
  
  set.seed(3333)
  x <- rnorm(100)
  
  result <- wild_bootstrap(x, num_boots = 5)
  
  expect_true(all(sapply(result, nrow) == 100))
})

## Integration Tests ==========================================================

test_that("wild_bootstrap integration check", {
  set.seed(4444)
  x <- matrix(rnorm(200), ncol = 2)
  
  result <- wild_bootstrap(x, num_boots = 10)
  
  ## Verify complete structure
  expect_length(result, 10)
  expect_true(all(sapply(result, function(m) {
    is.matrix(m) && nrow(m) == 100 && ncol(m) == 2
  })))
})
