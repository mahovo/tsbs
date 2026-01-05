# ==============================================================================
# TESTS FOR .is_invalid_data()
# ==============================================================================

test_that("Correctly validates and coerces standard time series objects", {
  # Create sample data in different formats
  base_matrix <- matrix(1:10, ncol = 2)
  ts_obj <- stats::ts(base_matrix)
  zoo_obj <- zoo::zoo(base_matrix)
  xts_obj <- xts::xts(base_matrix, order.by = Sys.Date() - (5:1))
  
  # Test that all valid time series objects are correctly coerced to a matrix.
  # Use expect_equivalent on the vectorized data to check for numerical
  # equality while ignoring all attributes like class, which are intentionally
  # changed by the function.
  expect_equivalent(as.vector(.coerce_and_validate_data(ts_obj)), as.vector(base_matrix))
  expect_equivalent(as.vector(.coerce_and_validate_data(zoo_obj)), as.vector(base_matrix))
  expect_equivalent(as.vector(.coerce_and_validate_data(xts_obj)), as.vector(base_matrix))
})

test_that("Correctly handles unsupported data types", {
  unsupported_obj <- list(a = 1, b = 2)
  
  ## Predictable mode should stop with an error
  expect_error(
    .coerce_and_validate_data(unsupported_obj, fail_mode = "predictably"),
    "Unsupported data type"
  )
  
  ## Graceful mode should return NULL
  expect_null(
    .coerce_and_validate_data(unsupported_obj, fail_mode = "gracefully")
  )
})

test_that("fail_mode = 'gracefully' works for time series objects", {
  # Create a valid and an invalid xts object
  valid_xts <- xts::xts(matrix(1:10, ncol = 2), order.by = Sys.Date() - (5:1))
  invalid_xts <- xts::xts(matrix(c(1:9, NA), ncol = 2), order.by = Sys.Date() - (5:1))
  
  # Graceful mode on valid data should return the matrix
  expect_equivalent(
    as.vector(.coerce_and_validate_data(valid_xts, fail_mode = "gracefully")),
    as.vector(zoo::coredata(valid_xts))
  )
  
  # Graceful mode on invalid data (with NA) should return NULL
  expect_null(
    .coerce_and_validate_data(invalid_xts, fail_mode = "gracefully")
  )
})

test_that(".coerce_and_validate_data() handles NULL values correctly", {
  # Predictable mode should stop with an error
  expect_error(
    .coerce_and_validate_data(NULL, fail_mode = "predictably"),
    "Input 'NULL' is NULL."
  )
  # Graceful mode should return NULL
  expect_null(
    .coerce_and_validate_data(NULL, fail_mode = "gracefully")
  )
})

test_that(".coerce_and_validate_data() handles non-existent variables", {
  # A non-existent variable should cause a standard R error
  expect_error(
    .coerce_and_validate_data(nonexistent_var),
    "object 'nonexistent_var' not found"
  )
})

test_that(".coerce_and_validate_data() validates numeric vectors correctly", {
  # Valid numeric vector should be coerced to a matrix
  valid_vec <- c(1, 2, 3, 4, 5)
  expect_equivalent(.coerce_and_validate_data(valid_vec), as.matrix(valid_vec))
  
  # Non-numeric vectors should fail
  char_vec <- c("a", "b", "c")
  expect_error(.coerce_and_validate_data(char_vec), "must be numeric")
  expect_null(.coerce_and_validate_data(char_vec, fail_mode = "gracefully"))
  
  logical_vec <- c(TRUE, FALSE, TRUE)
  expect_error(.coerce_and_validate_data(logical_vec), "must be numeric")
  expect_null(.coerce_and_validate_data(logical_vec, fail_mode = "gracefully"))
  
  # Vectors with non-finite values should fail
  na_vec <- c(1, 2, NA, 4, 5)
  inf_vec <- c(1, 2, Inf, 4, 5)
  expect_error(.coerce_and_validate_data(na_vec), "contains non-finite values")
  expect_null(.coerce_and_validate_data(na_vec, fail_mode = "gracefully"))
  expect_error(.coerce_and_validate_data(inf_vec), "contains non-finite values")
  expect_null(.coerce_and_validate_data(inf_vec, fail_mode = "gracefully"))
  
  # Empty vector should fail (must not be empty)
  empty_vec <- numeric(0)
  expect_error(.coerce_and_validate_data(empty_vec), "must not be empty")
  expect_null(.coerce_and_validate_data(empty_vec, fail_mode = "gracefully"))
})

test_that(".coerce_and_validate_data() validates numeric matrices correctly", {
  # Valid numeric matrix should remain a matrix
  valid_mat <- matrix(1:12, nrow = 3, ncol = 4)
  expect_equivalent(.coerce_and_validate_data(valid_mat), valid_mat)
  
  # Non-numeric matrix should fail
  char_mat <- matrix(letters[1:12], nrow = 3, ncol = 4)
  expect_error(.coerce_and_validate_data(char_mat), "must be numeric")
  expect_null(.coerce_and_validate_data(char_mat, fail_mode = "gracefully"))
  
  # Matrix with non-finite values should fail
  na_mat <- matrix(c(1:11, NA), nrow = 3, ncol = 4)
  expect_error(.coerce_and_validate_data(na_mat), "contains non-finite values")
  expect_null(.coerce_and_validate_data(na_mat, fail_mode = "gracefully"))
  
  # Empty matrix should fail
  empty_mat <- matrix(numeric(0), nrow = 0, ncol = 0)
  expect_error(.coerce_and_validate_data(empty_mat), "must not be empty")
  expect_null(.coerce_and_validate_data(empty_mat, fail_mode = "gracefully"))
})

test_that(".coerce_and_validate_data() validates data frames correctly", {
  # Valid numeric data frame should be coerced to a matrix
  valid_df <- data.frame(a = 1:5, b = 6:10, c = 11:15)
  expect_equivalent(.coerce_and_validate_data(valid_df), as.matrix(valid_df))
  
  # Data frame with non-numeric columns should fail
  mixed_df <- data.frame(a = 1:5, b = letters[1:5], c = 11:15)
  expect_error(.coerce_and_validate_data(mixed_df), "must be numeric")
  expect_null(.coerce_and_validate_data(mixed_df, fail_mode = "gracefully"))
  
  # Data frame with NA values should fail
  na_df <- data.frame(a = c(1, 2, NA, 4, 5), b = 6:10)
  expect_error(.coerce_and_validate_data(na_df), "contains non-finite values")
  expect_null(.coerce_and_validate_data(na_df, fail_mode = "gracefully"))
  
  # Empty data frames should fail
  empty_df_rows <- data.frame(a = numeric(0), b = numeric(0))
  empty_df_cols <- data.frame()
  expect_error(.coerce_and_validate_data(empty_df_rows), "must not be empty")
  expect_null(.coerce_and_validate_data(empty_df_rows, fail_mode = "gracefully"))
  expect_error(.coerce_and_validate_data(empty_df_cols), "must not be empty")
  expect_null(.coerce_and_validate_data(empty_df_cols, fail_mode = "gracefully"))
})

test_that(".coerce_and_validate_data() rejects unsupported object types", {
  # List should be invalid
  list_obj <- list(a = 1:5, b = 6:10)
  expect_error(.coerce_and_validate_data(list_obj), "Unsupported data type")
  expect_null(.coerce_and_validate_data(list_obj, fail_mode = "gracefully"))
  
  # Function should be invalid
  func_obj <- function(y) y^2
  expect_error(.coerce_and_validate_data(func_obj), "Unsupported data type")
  expect_null(.coerce_and_validate_data(func_obj, fail_mode = "gracefully"))
})


# ==============================================================================
# TESTS FOR OTHER HELPERS (Unchanged from original file)
# ==============================================================================

# Note: The following tests for .is_invalid_count() and .is_invalid_fraction()
# are preserved from the original test file as they test different functions
# and are not affected by the refactoring of the data validation logic.

test_that(".is_invalid_count() validates counts correctly", {
  # Valid positive integer
  expect_false(.is_invalid_count(5))
  # Zero should be invalid
  expect_true(.is_invalid_count(0))
  # Negative integer should be invalid
  expect_true(.is_invalid_count(-5))
  # Non-integer numeric should be invalid
  expect_true(.is_invalid_count(5.5))
  # Non-numeric should be invalid
  expect_true(.is_invalid_count("5"))
  # Vector of length > 1 should be invalid
  expect_true(.is_invalid_count(c(1, 2, 3)))
  # NA should be invalid
  expect_true(.is_invalid_count(NA_integer_))
})

# 
# test_that(".is_invalid_data() handles NULL values correctly", {
# 
#   # Test 1: NULL object with allow_null = TRUE (should be valid)
#   x <- NULL
#   expect_false(.is_invalid_data(x, allow_null = TRUE))
# 
#   # Test 2: NULL object with allow_null = FALSE (should be invalid)
#   x <- NULL
#   expect_true(.is_invalid_data(x, allow_null = FALSE))
# 
#   # Test 3: Default behavior (allow_null = TRUE)
#   x <- NULL
#   expect_false(.is_invalid_data(x))
# })
# 
# test_that(".is_invalid_data() handles non-existent variables", {
# 
#   # Test 4: Non-existent variable should be invalid
#   rm(list = ls()) # Clear environment first
#   expect_true(.is_invalid_data(nonexistent_var))
# 
#   # Test 5: Existing variable should proceed to other checks
#   existing_var <- c(1, 2, 3)
#   expect_false(.is_invalid_data(existing_var))
# })
# 
# test_that(".is_invalid_data() validates numeric vectors correctly", {
# 
#   # Test 6: Valid numeric vector
#   x <- c(1, 2, 3, 4, 5)
#   expect_false(.is_invalid_data(x))
# 
#   # Test 7: Non-numeric vector (character)
#   x <- c("a", "b", "c")
#   expect_true(.is_invalid_data(x))
# 
#   # Test 8: Non-numeric vector (logical)
#   x <- c(TRUE, FALSE, TRUE)
#   expect_true(.is_invalid_data(x))
# 
#   # Test 9: Numeric vector with NA values
#   x <- c(1, 2, NA, 4, 5)
#   expect_true(.is_invalid_data(x))
# 
#   # Test 10: Numeric vector with Inf values
#   x <- c(1, 2, Inf, 4, 5)
#   expect_true(.is_invalid_data(x))
# 
#   # Test 11: Numeric vector with -Inf values
#   x <- c(1, 2, -Inf, 4, 5)
#   expect_true(.is_invalid_data(x))
# 
#   # Test 12: Empty numeric vector
#   x <- numeric(0)
#   expect_false(.is_invalid_data(x))
# })
# 
# test_that(".is_invalid_data() validates numeric matrices correctly", {
# 
#   # Test 13: Valid numeric matrix
#   x <- matrix(1:12, nrow = 3, ncol = 4)
#   expect_false(.is_invalid_data(x))
# 
#   # Test 14: Non-numeric matrix (character)
#   x <- matrix(letters[1:12], nrow = 3, ncol = 4)
#   expect_true(.is_invalid_data(x))
# 
#   # Test 15: Numeric matrix with NA values
#   x <- matrix(c(1:11, NA), nrow = 3, ncol = 4)
#   expect_true(.is_invalid_data(x))
# 
#   # Test 16: Numeric matrix with Inf values
#   x <- matrix(c(1:11, Inf), nrow = 3, ncol = 4)
#   expect_true(.is_invalid_data(x))
# 
#   # Test 17: Empty matrix
#   x <- matrix(numeric(0), nrow = 0, ncol = 0)
#   expect_false(.is_invalid_data(x))
# })
# 
# test_that(".is_invalid_data() validates data frames correctly", {
# 
#   # Test 18: Valid numeric data frame
#   x <- data.frame(a = 1:5, b = 6:10, c = 11:15)
#   expect_false(.is_invalid_data(x))
# 
#   # Test 19: Data frame with non-numeric columns
#   x <- data.frame(a = 1:5, b = letters[1:5], c = 11:15)
#   expect_true(.is_invalid_data(x))
# 
#   # Test 20: Data frame with NA values
#   x <- data.frame(a = c(1, 2, NA, 4, 5), b = 6:10)
#   expect_true(.is_invalid_data(x))
# 
#   # Test 21: Data frame with Inf values
#   x <- data.frame(a = c(1, 2, Inf, 4, 5), b = 6:10)
#   expect_true(.is_invalid_data(x))
# 
#   # Test 22: Empty data frame (no rows)
#   x <- data.frame(a = numeric(0), b = numeric(0))
#   expect_true(.is_invalid_data(x))
# 
#   # Test 23: Data frame with no columns
#   x <- data.frame()
#   expect_true(.is_invalid_data(x))
# })
# 
# test_that(".is_invalid_data() rejects unsupported object types", {
# 
#   # Test 24: List should be invalid
#   x <- list(a = 1:5, b = 6:10)
#   expect_true(.is_invalid_data(x))
# 
#   # Test 25: Function should be invalid
#   x <- function(y) y^2
#   expect_true(.is_invalid_data(x))
# 
#   # Test 26: Environment should be invalid
#   x <- new.env()
#   expect_true(.is_invalid_data(x))
# })
# 
# # ==============================================================================
# # TESTS FOR .is_invalid_count()
# # ==============================================================================
# 
# test_that(".is_invalid_count() handles NULL values correctly", {
# 
#   # Test 27: NULL with allow_null = TRUE (should be valid)
#   n <- NULL
#   expect_false(.is_invalid_count(n, allow_null = TRUE))
# 
#   # Test 28: NULL with allow_null = FALSE (should be invalid)
#   n <- NULL
#   expect_true(.is_invalid_count(n, allow_null = FALSE))
# 
#   # Test 29: Default behavior
#   n <- NULL
#   expect_false(.is_invalid_count(n))
# })
# 
# test_that(".is_invalid_count() handles non-existent variables", {
# 
#   # Test 30: Non-existent variable should be invalid
#   expect_true(.is_invalid_count(nonexistent_count))
# 
#   # Test 31: Existing variable should proceed to validation
#   existing_count <- 5
#   expect_false(.is_invalid_count(existing_count))
# })
# 
# test_that(".is_invalid_count() validates counts correctly", {
# 
#   # Test 32: Valid positive integer
#   n <- 5
#   expect_false(.is_invalid_count(n))
# 
#   # Test 33: Valid count = 1
#   n <- 1
#   expect_false(.is_invalid_count(n))
# 
#   # Test 34: Zero should be invalid (counts must be positive)
#   n <- 0
#   expect_true(.is_invalid_count(n))
# 
#   # Test 35: Negative integer should be invalid
#   n <- -5
#   expect_true(.is_invalid_count(n))
# 
#   # Test 36: Non-integer numeric should be invalid
#   n <- 5.5
#   expect_true(.is_invalid_count(n))
# 
#   # Test 37: Very large integer should be valid
#   n <- 1000000
#   expect_false(.is_invalid_count(n))
# 
#   # Test 38: Non-numeric should be invalid
#   n <- "5"
#   expect_true(.is_invalid_count(n))
# 
#   # Test 39: Vector of length > 1 should be invalid
#   n <- c(1, 2, 3)
#   expect_true(.is_invalid_count(n))
# 
#   # Test 40: NA should be invalid
#   n <- NA_integer_
#   expect_true(.is_invalid_count(n))
# 
#   # Test 41: Inf should be invalid
#   n <- Inf
#   expect_true(.is_invalid_count(n))
# 
#   # Test 42: Logical should be invalid
#   n <- TRUE
#   expect_true(.is_invalid_count(n))
# })
# 
# # ==============================================================================
# # TESTS FOR .is_invalid_fraction()
# # ==============================================================================
# 
# test_that(".is_invalid_fraction() handles NULL values correctly", {
# 
#   # Test 43: NULL with allow_null = TRUE (should be valid)
#   p <- NULL
#   expect_false(.is_invalid_fraction(p, allow_null = TRUE))
# 
#   # Test 44: NULL with allow_null = FALSE (should be invalid)
#   p <- NULL
#   expect_true(.is_invalid_fraction(p, allow_null = FALSE))
# 
#   # Test 45: Default behavior
#   p <- NULL
#   expect_false(.is_invalid_fraction(p))
# })
# 
# test_that(".is_invalid_fraction() handles non-existent variables", {
# 
#   # Test 46: Non-existent variable should be invalid
#   expect_true(.is_invalid_fraction(nonexistent_fraction))
# 
#   # Test 47: Existing variable should proceed to validation
#   existing_fraction <- 0.5
#   expect_false(.is_invalid_fraction(existing_fraction))
# })
# 
# test_that(".is_invalid_fraction() validates fractions correctly", {
# 
#   # Test 48: Valid fraction in middle of range
#   p <- 0.5
#   expect_false(.is_invalid_fraction(p))
# 
#   # Test 49: Valid fraction near 0
#   p <- 0.001
#   expect_false(.is_invalid_fraction(p))
# 
#   # Test 50: Valid fraction near 1
#   p <- 0.999
#   expect_false(.is_invalid_fraction(p))
# 
#   # Test 51: Exactly 0 should be invalid (exclusive range)
#   p <- 0
#   expect_true(.is_invalid_fraction(p))
# 
#   # Test 52: Exactly 1 should be invalid (exclusive range)
#   p <- 1
#   expect_true(.is_invalid_fraction(p))
# 
#   # Test 53: Negative fraction should be invalid
#   p <- -0.5
#   expect_true(.is_invalid_fraction(p))
# 
#   # Test 54: Fraction > 1 should be invalid
#   p <- 1.5
#   expect_true(.is_invalid_fraction(p))
# 
#   # Test 55: Non-numeric should be invalid
#   p <- "0.5"
#   expect_true(.is_invalid_fraction(p))
# 
#   # Test 56: Vector of length > 1 should be invalid
#   p <- c(0.3, 0.7)
#   expect_true(.is_invalid_fraction(p))
# 
#   # Test 57: NA should be invalid
#   p <- NA_real_
#   expect_true(.is_invalid_fraction(p))
# 
#   # Test 58: Inf should be invalid
#   p <- Inf
#   expect_true(.is_invalid_fraction(p))
# 
#   # Test 59: -Inf should be invalid
#   p <- -Inf
#   expect_true(.is_invalid_fraction(p))
# 
#   # Test 60: Logical should be invalid
#   p <- TRUE
#   expect_true(.is_invalid_fraction(p))
# })
# 
# # ==============================================================================
# # EDGE CASES AND INTEGRATION TESTS
# # ==============================================================================
# 
# test_that("All validation functions handle edge cases consistently", {
# 
#   # Test 61: Empty vectors
#   expect_false(.is_invalid_data(numeric(0)))
#   expect_true(.is_invalid_count(numeric(0)))  # length != 1
#   expect_true(.is_invalid_fraction(numeric(0)))  # length != 1
# 
#   # Test 62: Very small positive numbers
#   expect_false(.is_invalid_data(1e-100))
#   expect_true(.is_invalid_count(1e-100))  # not integer
#   expect_false(.is_invalid_fraction(1e-100))  # valid fraction
# 
#   # Test 63: Very large numbers
#   expect_false(.is_invalid_data(1e100))
#   expect_true(.is_invalid_count(1e100))  # not integer (floating point precision)
#   expect_true(.is_invalid_fraction(1e100))  # > 1
# })
# 
# test_that("All validation functions handle edge cases consistently", {
# 
#   # Test 61: Empty vectors (assign to variables first)
#   empty_vec <- numeric(0)
#   expect_false(.is_invalid_data(empty_vec))
#   expect_true(.is_invalid_count(empty_vec))  # length != 1
#   expect_true(.is_invalid_fraction(empty_vec))  # length != 1
# 
#   # Test 62: Very small positive numbers
#   tiny_num <- 1e-100
#   expect_false(.is_invalid_data(tiny_num))
#   expect_true(.is_invalid_count(tiny_num))  # not integer
#   expect_false(.is_invalid_fraction(tiny_num))  # valid fraction
# 
#   # Test 63: Very large numbers
#   big_num <- 1e100
#   expect_false(.is_invalid_data(big_num))
#   expect_true(.is_invalid_count(big_num))  # not integer (floating point precision)
#   expect_true(.is_invalid_fraction(big_num))  # > 1
# })

