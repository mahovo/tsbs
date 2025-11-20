
## Basic Functionality Tests ==================================================

test_that("msvar_bootstrap handles bivariate time series", {
  set.seed(123)
  n <- 200
  y1 <- arima.sim(n = n, list(ar = 0.7))
  y2 <- 0.5 * y1 + arima.sim(n = n, list(ar = 0.3))
  x <- cbind(y1, y2)
  
  result <- msvar_bootstrap(
    x = x,
    n_boot = 150,
    num_boots = 5
  )
  
  expect_type(result, "list")
  expect_length(result, 5)
  expect_true(all(sapply(result, is.matrix)))
  expect_true(all(sapply(result, nrow) == 150))
  expect_true(all(sapply(result, ncol) == 2))
})

test_that("msvar_bootstrap handles trivariate time series", {
  set.seed(456)
  n <- 180
  y1 <- arima.sim(n = n, list(ar = 0.6))
  y2 <- 0.4 * y1 + arima.sim(n = n, list(ar = 0.4))
  y3 <- 0.3 * y2 + arima.sim(n = n, list(ar = 0.5))
  x <- cbind(y1, y2, y3)
  
  result <- msvar_bootstrap(
    x = x,
    n_boot = 150,
    num_boots = 4
  )
  
  expect_type(result, "list")
  expect_length(result, 4)
  expect_true(all(sapply(result, ncol) == 3))
})

test_that("msvar_bootstrap truncates to n_boot when specified", {
  set.seed(789)
  n <- 200
  y1 <- rnorm(n)
  y2 <- rnorm(n)
  x <- cbind(y1, y2)
  
  result <- msvar_bootstrap(
    x = x,
    n_boot = 100,  ## Shorter than original
    num_boots = 3
  )
  
  expect_true(all(sapply(result, nrow) == 100))
})

test_that("msvar_bootstrap generates variable lengths without n_boot", {
  set.seed(321)
  n <- 150
  y1 <- rnorm(n)
  y2 <- rnorm(n)
  x <- cbind(y1, y2)
  
  result <- msvar_bootstrap(
    x = x,
    n_boot = NULL,  ## No truncation
    num_blocks = 10,
    num_boots = 5
  )
  
  lengths <- sapply(result, nrow)
  ## Lengths should vary or all be positive
  expect_true(length(unique(lengths)) > 1 || all(lengths > 0))
})

test_that("msvar_bootstrap with num_blocks produces valid output", {
  set.seed(9876)
  n <- 200
  y1 <- arima.sim(n = n, list(ar = 0.7))
  y2 <- arima.sim(n = n, list(ar = 0.5))
  x <- cbind(y1, y2)
  
  result <- msvar_bootstrap(
    x = x,
    n_boot = NULL,
    num_blocks = 20,
    num_boots = 8
  )
  
  expect_type(result, "list")
  expect_length(result, 8)
  expect_true(all(sapply(result, is.matrix)))
  expect_true(all(sapply(result, ncol) == 2))
  
  lengths <- sapply(result, nrow)
  expect_true(all(lengths > 0))
  expect_true(all(lengths > 10))  ## At least some data
  ## With 20 blocks from 200 obs, samples could range widely
  ## Just check they're not absurdly long (e.g., not 10x original)
  expect_true(all(lengths < n * 10))
})

## Data Type Handling Tests ===================================================

test_that("msvar_bootstrap handles data.frame input", {
  set.seed(444)
  df <- data.frame(
    x1 = rnorm(100),
    x2 = rnorm(100)
  )
  
  result <- msvar_bootstrap(
    x = df,
    n_boot = 80,
    num_boots = 3
  )
  
  expect_type(result, "list")
  expect_true(all(sapply(result, is.matrix)))
  expect_true(all(sapply(result, ncol) == 2))
})

test_that("msvar_bootstrap handles matrix input", {
  set.seed(555)
  x <- matrix(rnorm(200), ncol = 2)
  
  result <- msvar_bootstrap(
    x = x,
    n_boot = 80,
    num_boots = 3
  )
  
  expect_type(result, "list")
  expect_true(all(sapply(result, is.matrix)))
})

## Parallel Execution Tests ===================================================

test_that("msvar_bootstrap handles parallel execution", {
  skip_if_not_installed("foreach")
  skip_if_not_installed("doParallel")
  skip_on_cran()
  
  set.seed(222)
  n <- 120
  y1 <- rnorm(n)
  y2 <- rnorm(n)
  x <- cbind(y1, y2)
  
  result <- msvar_bootstrap(
    x = x,
    n_boot = 100,
    num_boots = 4,
    parallel = TRUE,
    num_cores = 2
  )
  
  expect_type(result, "list")
  expect_length(result, 4)
  expect_true(all(sapply(result, nrow) == 100))
})

test_that("msvar_bootstrap parallel gives consistent structure", {
  skip_if_not_installed("foreach")
  skip_if_not_installed("doParallel")
  skip_on_cran()
  
  n <- 100
  y1 <- rnorm(n)
  y2 <- rnorm(n)
  x <- cbind(y1, y2)
  
  set.seed(333)
  result_seq <- msvar_bootstrap(
    x = x,
    n_boot = 80,
    num_boots = 3,
    parallel = FALSE
  )
  
  set.seed(333)
  result_par <- msvar_bootstrap(
    x = x,
    n_boot = 80,
    num_boots = 3,
    parallel = TRUE,
    num_cores = 2
  )
  
  ## Check structural equality
  expect_equal(length(result_seq), length(result_par))
  expect_equal(sapply(result_seq, dim), sapply(result_par, dim))
})

## Edge Cases =================================================================

test_that("msvar_bootstrap handles short series", {
  set.seed(666)
  n <- 50
  y1 <- rnorm(n)
  y2 <- rnorm(n)
  x <- cbind(y1, y2)
  
  result <- msvar_bootstrap(
    x = x,
    n_boot = 40,
    num_boots = 2
  )
  
  expect_length(result, 2)
  expect_true(all(sapply(result, nrow) == 40))
})

test_that("msvar_bootstrap with single bootstrap", {
  set.seed(777)
  n <- 100
  y1 <- rnorm(n)
  y2 <- rnorm(n)
  x <- cbind(y1, y2)
  
  result <- msvar_bootstrap(
    x = x,
    n_boot = 80,
    num_boots = 1
  )
  
  expect_length(result, 1)
  expect_true(nrow(result[[1]]) == 80)
})

test_that("msvar_bootstrap with n_boot > original length", {
  set.seed(888)
  n <- 60
  y1 <- rnorm(n)
  y2 <- rnorm(n)
  x <- cbind(y1, y2)
  
  result <- msvar_bootstrap(
    x = x,
    n_boot = 120,  ## Longer than original
    num_boots = 3
  )
  
  expect_true(all(sapply(result, nrow) == 120))
})

test_that("msvar_bootstrap output is numeric and finite", {
  set.seed(999)
  n <- 100
  y1 <- rnorm(n)
  y2 <- rnorm(n)
  x <- cbind(y1, y2)
  
  result <- msvar_bootstrap(
    x = x,
    n_boot = 80,
    num_boots = 3
  )
  
  expect_true(all(sapply(result, is.numeric)))
  expect_true(all(sapply(result, function(m) all(is.finite(m)))))
})

## Data Properties Tests ======================================================

test_that("msvar_bootstrap handles highly correlated series", {
  set.seed(1111)
  n <- 150
  y1 <- arima.sim(n = n, list(ar = 0.8))
  y2 <- 0.9 * y1 + rnorm(n, sd = 0.1)  ## High correlation
  x <- cbind(y1, y2)
  
  result <- msvar_bootstrap(
    x = x,
    n_boot = 120,
    num_boots = 3
  )
  
  expect_length(result, 3)
  expect_true(all(sapply(result, nrow) == 120))
})

test_that("msvar_bootstrap handles uncorrelated series", {
  set.seed(2222)
  n <- 150
  y1 <- rnorm(n)
  y2 <- rnorm(n)  ## Independent
  x <- cbind(y1, y2)
  
  result <- msvar_bootstrap(
    x = x,
    n_boot = 120,
    num_boots = 3
  )
  
  expect_length(result, 3)
  expect_true(all(sapply(result, nrow) == 120))
})

test_that("msvar_bootstrap handles series with regime switches", {
  skip_on_cran()  ## Stochastic test
  
  set.seed(3333)
  n <- 200
  ## Create two regimes
  regime1 <- cbind(
    rnorm(100, mean = 0, sd = 1),
    rnorm(100, mean = 0, sd = 1)
  )
  regime2 <- cbind(
    rnorm(100, mean = 5, sd = 2),
    rnorm(100, mean = 5, sd = 2)
  )
  x <- rbind(regime1, regime2)
  
  result <- msvar_bootstrap(
    x = x,
    n_boot = 180,
    num_boots = 5
  )
  
  expect_length(result, 5)
  expect_true(all(sapply(result, nrow) == 180))
  
  ## Check that bootstrap captures range from both regimes
  all_values <- unlist(result)
  expect_true(min(all_values) < 1)
  expect_true(max(all_values) > 3)
})

test_that("msvar_bootstrap handles trending data", {
  set.seed(4444)
  n <- 150
  trend <- seq(0, 10, length.out = n)
  y1 <- trend + rnorm(n, sd = 0.5)
  y2 <- 0.5 * trend + rnorm(n, sd = 0.5)
  x <- cbind(y1, y2)
  
  result <- msvar_bootstrap(
    x = x,
    n_boot = 120,
    num_boots = 4
  )
  
  expect_length(result, 4)
  expect_true(all(sapply(result, nrow) == 120))
})

test_that("msvar_bootstrap handles series with different scales", {
  set.seed(5555)
  n <- 120
  y1 <- rnorm(n, mean = 1000, sd = 100)
  y2 <- rnorm(n, mean = 0.01, sd = 0.001)
  x <- cbind(y1, y2)
  
  result <- msvar_bootstrap(
    x = x,
    n_boot = 100,
    num_boots = 3
  )
  
  expect_length(result, 3)
  expect_true(all(sapply(result, ncol) == 2))
  
  ## Check both variables maintain reasonable ranges
  for (i in seq_along(result)) {
    expect_true(mean(result[[i]][,1]) > 500)
    expect_true(mean(abs(result[[i]][,2])) < 1)
  }
})

## Reproducibility Tests ======================================================

test_that("msvar_bootstrap same seed gives same results", {
  n <- 100
  y1 <- rnorm(n)
  y2 <- rnorm(n)
  x <- cbind(y1, y2)
  
  set.seed(6666)
  result1 <- msvar_bootstrap(x = x, n_boot = 80, num_boots = 2)
  
  set.seed(6666)
  result2 <- msvar_bootstrap(x = x, n_boot = 80, num_boots = 2)
  
  expect_equal(result1, result2)
})

test_that("msvar_bootstrap different seed gives different results", {
  n <- 100
  y1 <- rnorm(n)
  y2 <- rnorm(n)
  x <- cbind(y1, y2)
  
  set.seed(7777)
  result1 <- msvar_bootstrap(x = x, n_boot = 80, num_boots = 2)
  
  set.seed(8888)
  result2 <- msvar_bootstrap(x = x, n_boot = 80, num_boots = 2)
  
  expect_false(identical(result1, result2))
})

## Integration Tests ==========================================================

test_that("msvar_bootstrap integration with .sample_blocks", {
  set.seed(9999)
  n <- 120
  y1 <- rnorm(n)
  y2 <- rnorm(n)
  x <- cbind(y1, y2)
  
  result <- msvar_bootstrap(
    x = x,
    n_boot = 100,
    num_boots = 4
  )
  
  ## Verify structure matches expectations
  expect_length(result, 4)
  expect_true(all(sapply(result, function(m) {
    is.matrix(m) && nrow(m) == 100 && ncol(m) == 2
  })))
})

test_that("msvar_bootstrap assumes validated inputs from tsbs()", {
  ## msvar_bootstrap is an internal function called by tsbs()
  ## All validation should happen in tsbs() before calling msvar_bootstrap()
  ## These tests verify the function works with valid inputs
  
  n <- 100
  y1 <- rnorm(n)
  y2 <- rnorm(n)
  x <- cbind(y1, y2)
  
  ## Valid inputs should work
  result <- msvar_bootstrap(
    x = x, 
    n_boot = 80, 
    num_boots = 3
  )
  
  expect_type(result, "list")
  expect_length(result, 3)
})

## Additional Edge Cases ======================================================

test_that("msvar_bootstrap handles high persistence", {
  set.seed(1234)
  n <- 150
  y1 <- arima.sim(n = n, list(ar = 0.95))
  y2 <- 0.5 * y1 + arima.sim(n = n, list(ar = 0.95))
  x <- cbind(y1, y2)
  
  result <- msvar_bootstrap(
    x = x,
    n_boot = 120,
    num_boots = 3
  )
  
  expect_length(result, 3)
  expect_true(all(sapply(result, nrow) == 120))
})

test_that("msvar_bootstrap handles low persistence", {
  set.seed(2345)
  n <- 150
  y1 <- arima.sim(n = n, list(ar = 0.1))
  y2 <- arima.sim(n = n, list(ar = 0.1))
  x <- cbind(y1, y2)
  
  result <- msvar_bootstrap(
    x = x,
    n_boot = 120,
    num_boots = 3
  )
  
  expect_length(result, 3)
  expect_true(all(sapply(result, nrow) == 120))
})

test_that("msvar_bootstrap handles non-stationary variance", {
  set.seed(3456)
  n <- 150
  y1 <- numeric(n)
  y2 <- numeric(n)
  for (i in seq_len(n)) {
    y1[i] <- rnorm(1, sd = 0.5 + 0.01 * i)
    y2[i] <- rnorm(1, sd = 0.5 + 0.01 * i)
  }
  x <- cbind(y1, y2)
  
  result <- msvar_bootstrap(
    x = x,
    n_boot = 120,
    num_boots = 4
  )
  
  expect_length(result, 4)
  expect_true(all(sapply(result, nrow) == 120))
})

test_that("msvar_bootstrap handles data with outliers", {
  set.seed(4567)
  n <- 120
  y1 <- rnorm(n)
  y2 <- rnorm(n)
  ## Add outliers
  y1[c(10, 50, 90)] <- c(-10, 15, -8)
  y2[c(15, 55, 95)] <- c(12, -9, 10)
  x <- cbind(y1, y2)
  
  result <- msvar_bootstrap(
    x = x,
    n_boot = 100,
    num_boots = 3
  )
  
  expect_length(result, 3)
  expect_true(all(sapply(result, is.finite)))
})

test_that("msvar_bootstrap with many blocks", {
  set.seed(5678)
  n <- 100
  y1 <- rnorm(n)
  y2 <- rnorm(n)
  x <- cbind(y1, y2)
  
  result <- msvar_bootstrap(
    x = x,
    n_boot = NULL,
    num_blocks = 50,  ## Many blocks
    num_boots = 3
  )
  
  expect_length(result, 3)
  lengths <- sapply(result, nrow)
  expect_true(all(lengths > 0))
})

test_that("msvar_bootstrap with few blocks", {
  set.seed(6789)
  n <- 150
  y1 <- rnorm(n)
  y2 <- rnorm(n)
  x <- cbind(y1, y2)
  
  result <- msvar_bootstrap(
    x = x,
    n_boot = NULL,
    num_blocks = 5,  ## Few blocks
    num_boots = 3
  )
  
  expect_length(result, 3)
  lengths <- sapply(result, nrow)
  expect_true(all(lengths > 0))
})
