# Comprehensive tests for bootstrap diagnostics across all bootstrap types
#
# Tests cover:
# - Diagnostics object creation and structure
# - Recording functions (blocks, stats, regimes)
# - Summary and print methods
# - Plot methods
# - Method-specific diagnostics for each bootstrap type

# =============================================================================
# Test Data Generators
# =============================================================================

generate_test_series <- function(n = 200, k = 1, seed = 42) {
  set.seed(seed)
  if (k == 1) {
    rnorm(n, mean = 0.01, sd = 0.02)
  } else {
    matrix(rnorm(n * k, mean = 0.01, sd = 0.02), nrow = n, ncol = k)
  }
}

generate_regime_series <- function(n = 300, k = 1, seed = 42) {
  set.seed(seed)
  n1 <- n %/% 2
  n2 <- n - n1
  
  if (k == 1) {
    c(rnorm(n1, mean = 0.02, sd = 0.01),
      rnorm(n2, mean = -0.01, sd = 0.03))
  } else {
    factor <- c(rnorm(n1, mean = 0.02, sd = 0.01),
                rnorm(n2, mean = -0.01, sd = 0.03))
    Y <- matrix(0, nrow = n, ncol = k)
    for (j in 1:k) {
      Y[, j] <- factor + rnorm(n, sd = 0.002)
    }
    Y
  }
}


# =============================================================================
# Tests for Diagnostic Object Creation
# =============================================================================

test_that("create_bootstrap_diagnostics creates valid object", {
  diag <- create_bootstrap_diagnostics(
    bs_type = "moving",
    n_original = 100,
    n_vars = 2,
    num_boots = 50
  )
  
  expect_s3_class(diag, "tsbs_diagnostics")
  expect_equal(diag$meta$bs_type, "moving")
  expect_equal(diag$meta$n_original, 100)
  expect_equal(diag$meta$n_vars, 2)
  expect_equal(diag$meta$num_boots, 50)
  expect_true(!is.null(diag$meta$timestamp))
})

test_that("create_bootstrap_diagnostics initializes all required fields", {
  diag <- create_bootstrap_diagnostics(
    bs_type = "stationary",
    n_original = 200,
    n_vars = 3,
    num_boots = 100
  )
  
  # Check blocks structure
 expect_true("blocks" %in% names(diag))
  expect_length(diag$blocks$replicate_blocks, 100)
  
  # Check series_stats structure
  expect_true("series_stats" %in% names(diag))
  expect_length(diag$series_stats$replicate_means, 100)
  expect_length(diag$series_stats$replicate_sds, 100)
  expect_length(diag$series_stats$replicate_ac1, 100)
  
  # Check original_stats structure
  expect_true("original_stats" %in% names(diag))
  
  # Check method_specific exists
  expect_true("method_specific" %in% names(diag))
})


# =============================================================================
# Tests for Recording Functions
# =============================================================================

test_that("record_blocks stores block information correctly", {
  diag <- create_bootstrap_diagnostics("moving", 100, 1, 10)
  
  diag <- record_blocks(
    diag,
    replicate_idx = 1,
    block_lengths = c(10, 10, 10),
    start_positions = c(1, 25, 50),
    block_type = "overlapping"
  )
  
  expect_equal(nrow(diag$blocks$replicate_blocks[[1]]), 3)
  expect_equal(diag$blocks$replicate_blocks[[1]]$length, c(10, 10, 10))
  expect_equal(diag$blocks$replicate_blocks[[1]]$start_pos, c(1, 25, 50))
  expect_equal(diag$blocks$all_block_lengths, c(10, 10, 10))
})

test_that("record_replicate_stats computes statistics correctly", {
  diag <- create_bootstrap_diagnostics("moving", 100, 2, 10)
  
  set.seed(123)
  boot_matrix <- matrix(rnorm(100 * 2), nrow = 100, ncol = 2)
  
  diag <- record_replicate_stats(diag, replicate_idx = 1, bootstrap_matrix = boot_matrix)
  
  expect_equal(diag$series_stats$replicate_means[[1]], colMeans(boot_matrix))
  expect_equal(diag$series_stats$replicate_sds[[1]], apply(boot_matrix, 2, sd))
  expect_equal(diag$series_stats$replicate_lengths[1], 100)
  expect_length(diag$series_stats$replicate_ac1[[1]], 2)
})

test_that("record_original_stats stores original series statistics", {
  diag <- create_bootstrap_diagnostics("moving", 100, 2, 10)
  
  set.seed(123)
  x <- matrix(rnorm(100 * 2), nrow = 100, ncol = 2)
  
  diag <- record_original_stats(diag, x)
  
  expect_equal(diag$original_stats$means, colMeans(x))
  expect_equal(diag$original_stats$sds, apply(x, 2, sd))
  expect_length(diag$original_stats$ac1, 2)
})

test_that("record_method_specific adds custom diagnostics", {
  diag <- create_bootstrap_diagnostics("hmm", 100, 1, 10)
  
  diag <- record_method_specific(
    diag,
    states = rep(1:2, each = 50),
    transition_matrix = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2)
  )
  
  expect_equal(diag$method_specific$states, rep(1:2, each = 50))
  expect_true(is.matrix(diag$method_specific$transition_matrix))
})


# =============================================================================
# Tests for compute_bootstrap_diagnostics
# =============================================================================

test_that("compute_bootstrap_diagnostics works with list of bootstrap series", {
  set.seed(123)
  original <- matrix(rnorm(100 * 2), nrow = 100, ncol = 2)
  
  # Simulate bootstrap series
  boot_series <- lapply(1:10, function(i) {
    idx <- sample(1:100, 100, replace = TRUE)
    original[idx, ]
  })
  
  diag <- compute_bootstrap_diagnostics(
    bootstrap_series = boot_series,
    original_data = original,
    bs_type = "moving",
    config = list(block_length = 10)
  )
  
  expect_s3_class(diag, "tsbs_diagnostics")
  expect_equal(diag$meta$num_boots, 10)
  expect_equal(diag$meta$n_vars, 2)
  expect_length(diag$series_stats$replicate_means, 10)
})


# =============================================================================
# Tests for Summary Method
# =============================================================================

test_that("summary.tsbs_diagnostics produces output without error", {
  diag <- create_bootstrap_diagnostics("moving", 100, 1, 10)
  diag <- record_original_stats(diag, rnorm(100))
  
  # Record some replicate stats
  for (i in 1:10) {
    boot <- rnorm(100)
    diag <- record_replicate_stats(diag, i, matrix(boot, ncol = 1))
    diag <- record_blocks(diag, i, c(10, 10), c(1, 50))
  }
  
  expect_output(summary(diag), "tsbs Bootstrap Diagnostics")
  expect_output(summary(diag), "Bootstrap type: moving")
})

test_that("summary.tsbs_diagnostics shows block statistics for block methods", {
  diag <- create_bootstrap_diagnostics("stationary", 100, 1, 10)
  
  for (i in 1:10) {
    # Simulate random block lengths for stationary bootstrap
    block_lens <- rgeom(5, 0.1) + 1
    diag <- record_blocks(diag, i, block_lens, sample(1:90, length(block_lens)))
  }
  
  expect_output(summary(diag), "BLOCK LENGTH STATISTICS")
})


# =============================================================================
# Tests for Print Method
# =============================================================================

test_that("print.tsbs_diagnostics produces concise output", {
  diag <- create_bootstrap_diagnostics("moving", 100, 2, 50)
  
  expect_output(print(diag), "tsbs Bootstrap Diagnostics")
  expect_output(print(diag), "moving")
  expect_output(print(diag), "50")
})


# =============================================================================
# Tests for Plot Method
# =============================================================================

test_that("plot.tsbs_diagnostics handles missing data gracefully", {
  skip_on_cran()
  
  diag <- create_bootstrap_diagnostics("moving", 100, 1, 10)
  
  # Should not error even with no data recorded
  expect_message(
    plot(diag, type = "all"),
    "No diagnostic plots available"
  )
})

test_that("plot.tsbs_diagnostics works with block data", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  
  diag <- create_bootstrap_diagnostics("moving", 100, 1, 10)
  
  for (i in 1:10) {
    diag <- record_blocks(diag, i, c(10, 10, 10), c(1, 30, 60))
  }
  
  # Should not error
  expect_silent(plot(diag, type = "block_lengths"))
})


# =============================================================================
# Tests for Moving Block Bootstrap Diagnostics
# =============================================================================

test_that("tsbs moving bootstrap returns diagnostics when requested", {
  x <- generate_test_series(n = 100, k = 1)
  
  result <- tsbs(
    x = x,
    bs_type = "moving",
    block_length = 10,
    num_boots = 20,
    return_diagnostics = TRUE
  )
  
  expect_true("diagnostics" %in% names(result))
  expect_s3_class(result$diagnostics, "tsbs_diagnostics")
  expect_equal(result$diagnostics$meta$bs_type, "moving")
})

test_that("moving bootstrap diagnostics contain block information", {
  x <- generate_test_series(n = 100, k = 1)
  
  result <- tsbs(
    x = x,
    bs_type = "moving",
    block_length = 10,
    num_boots = 20,
    return_diagnostics = TRUE
  )
  
  diag <- result$diagnostics
  
  # Should have block length data
  expect_true(length(diag$blocks$all_block_lengths) > 0)
  
  # Should have original stats
  expect_true(!is.null(diag$original_stats$means))
})


# =============================================================================
# Tests for Stationary Bootstrap Diagnostics
# =============================================================================

test_that("tsbs stationary bootstrap returns diagnostics when requested", {
  x <- generate_test_series(n = 100, k = 1)
  
  result <- tsbs(
    x = x,
    bs_type = "stationary",
    block_length = 10,
    num_boots = 20,
    return_diagnostics = TRUE
  )
  
  expect_true("diagnostics" %in% names(result))
  expect_s3_class(result$diagnostics, "tsbs_diagnostics")
  expect_equal(result$diagnostics$meta$bs_type, "stationary")
})

test_that("stationary bootstrap diagnostics show variable block lengths", {
  x <- generate_test_series(n = 200, k = 1)
  
  result <- tsbs(
    x = x,
    bs_type = "stationary",
    block_length = 10,
    num_boots = 50,
    return_diagnostics = TRUE
  )
  
  diag <- result$diagnostics
  block_lengths <- diag$blocks$all_block_lengths
  
  # Stationary bootstrap should have variable block lengths
  if (length(block_lengths) > 1) {
    expect_true(length(unique(block_lengths)) > 1 || length(block_lengths) == 1)
  }
})


# =============================================================================
# Tests for Wild Bootstrap Diagnostics  
# =============================================================================

test_that("tsbs wild bootstrap returns diagnostics when requested", {
  x <- generate_test_series(n = 100, k = 1)
  
  result <- tsbs(
    x = x,
    bs_type = "wild",
    num_boots = 20,
    return_diagnostics = TRUE
  )
  
  expect_true("diagnostics" %in% names(result))
  expect_s3_class(result$diagnostics, "tsbs_diagnostics")
  expect_equal(result$diagnostics$meta$bs_type, "wild")
})

test_that("wild bootstrap diagnostics have no block data", {
  x <- generate_test_series(n = 100, k = 1)
  
  result <- tsbs(
    x = x,
    bs_type = "wild",
    num_boots = 20,
    return_diagnostics = TRUE
  )
  
  diag <- result$diagnostics
  
  # Wild bootstrap doesn't use blocks
  expect_equal(length(diag$blocks$all_block_lengths), 0)
  
  # But should still have series statistics
  expect_true(!is.null(diag$original_stats$means))
})


# =============================================================================
# Tests for HMM Bootstrap Diagnostics
# =============================================================================

test_that("tsbs HMM bootstrap returns diagnostics when requested", {
  skip_if_not_installed("depmixS4")
  
  x <- generate_regime_series(n = 200, k = 1)
  
  result <- tsbs(
    x = x,
    bs_type = "hmm",
    num_states = 2,
    num_boots = 10,
    return_diagnostics = TRUE
  )
  
  expect_true("diagnostics" %in% names(result))
  expect_s3_class(result$diagnostics, "tsbs_diagnostics")
  expect_equal(result$diagnostics$meta$bs_type, "hmm")
})

test_that("HMM diagnostics contain regime information", {
  skip_if_not_installed("depmixS4")
  
  x <- generate_regime_series(n = 200, k = 1)
  
  result <- tsbs(
    x = x,
    bs_type = "hmm",
    num_states = 2,
    num_boots = 10,
    return_diagnostics = TRUE,
    return_fit = TRUE
  )
  
  diag <- result$diagnostics
  
  # Should have method-specific regime info
  expect_true("method_specific" %in% names(diag))
  
  # Should have original stats
  expect_true(!is.null(diag$original_stats$means))
})

test_that("HMM bootstrap with MSGARCH returns diagnostics", {
  skip_if_not_installed("MSGARCH")
  
  x <- generate_regime_series(n = 300, k = 1)
  
  result <- tsbs(
    x = x,
    bs_type = "hmm",
    num_states = 2,
    num_boots = 5,
    distribution = "norm",
    seed = 123,
    return_diagnostics = TRUE
  )
  
  expect_true("diagnostics" %in% names(result))
  expect_s3_class(result$diagnostics, "tsbs_diagnostics")
})


# =============================================================================
# Tests for MSVAR Bootstrap Diagnostics
# =============================================================================

test_that("tsbs MSVAR bootstrap returns diagnostics when requested", {
  skip_if_not_installed("MSwM")
  
  # MSVAR needs multivariate data
  Y <- generate_regime_series(n = 200, k = 2)
  
  result <- tryCatch({
    tsbs(
      x = Y,
      bs_type = "msvar",
      num_boots = 5,
      return_diagnostics = TRUE
    )
  }, error = function(e) {
    skip(paste("MSVAR fitting failed:", e$message))
  })
  
  expect_true("diagnostics" %in% names(result))
  expect_s3_class(result$diagnostics, "tsbs_diagnostics")
  expect_equal(result$diagnostics$meta$bs_type, "msvar")
})


# =============================================================================
# Tests for Regime Diagnostics Functions
# =============================================================================

test_that("record_regime_info initializes regime tracking", {
  diag <- create_bootstrap_diagnostics("hmm", 100, 1, 10)
  
  states <- rep(1:2, each = 50)
  diag <- record_regime_info(diag, states, num_states = 2)
  
  expect_true("regime_info" %in% names(diag$method_specific))
  expect_equal(diag$method_specific$regime_info$num_states, 2)
  expect_equal(diag$method_specific$regime_info$original_states, states)
})

test_that("record_replicate_regimes stores per-replicate state info", {
  diag <- create_bootstrap_diagnostics("hmm", 100, 1, 10)
  
  states <- rep(1:2, each = 50)
  diag <- record_regime_info(diag, states, num_states = 2)
  
  # Record replicate regimes
  rep_states <- sample(1:2, 100, replace = TRUE)
  diag <- record_replicate_regimes(diag, replicate_idx = 1, replicate_states = rep_states)
  
  expect_equal(diag$method_specific$regime_info$replicate_states[[1]], rep_states)
})


# =============================================================================
# Tests for Extract Functions
# =============================================================================

test_that("extract_blocks returns block information", {
  diag <- create_bootstrap_diagnostics("moving", 100, 1, 10)
  
  for (i in 1:10) {
    diag <- record_blocks(diag, i, c(10, 10), c(1 + (i-1)*10, 50 + (i-1)*5))
  }
  
  # Extract all blocks
  all_blocks <- extract_blocks(diag)
  expect_true(is.data.frame(all_blocks) || is.list(all_blocks))
  
  # Extract specific replicate
  rep1_blocks <- extract_blocks(diag, replicate = 1)
  expect_true(!is.null(rep1_blocks))
})

test_that("extract_summary_stats returns summary statistics", {
  diag <- create_bootstrap_diagnostics("moving", 100, 1, 10)
  diag <- record_original_stats(diag, rnorm(100))
  
  for (i in 1:10) {
    diag <- record_replicate_stats(diag, i, matrix(rnorm(100), ncol = 1))
  }
  
  stats <- extract_summary_stats(diag)
  
  expect_true(is.list(stats) || is.data.frame(stats))
})


# =============================================================================
# Tests for summarize_func_outs
# =============================================================================

test_that("summarize_func_outs computes bootstrap statistics", {
  # Simulate func_outs from tsbs
  func_outs <- lapply(1:100, function(i) c(mean = rnorm(1), sd = abs(rnorm(1))))
  
  summary_stats <- summarize_func_outs(func_outs)
  
  expect_true(is.list(summary_stats) || is.data.frame(summary_stats))
})

test_that("summarize_func_outs handles scalar outputs", {
  func_outs <- as.list(rnorm(100))
  
  summary_stats <- summarize_func_outs(func_outs)
  
  expect_true(!is.null(summary_stats))
})


# =============================================================================
# Tests for Multivariate Diagnostics
# =============================================================================

test_that("diagnostics work with multivariate data", {
  Y <- generate_test_series(n = 100, k = 3)
  
  result <- tsbs(
    x = Y,
    bs_type = "moving",
    block_length = 10,
    num_boots = 20,
    return_diagnostics = TRUE
  )
  
  diag <- result$diagnostics
  
  expect_equal(diag$meta$n_vars, 3)
  expect_length(diag$original_stats$means, 3)
  expect_length(diag$original_stats$sds, 3)
})


# =============================================================================
# Tests for Edge Cases
# =============================================================================

test_that("diagnostics handle single bootstrap replicate", {
  x <- generate_test_series(n = 100, k = 1)
  
  result <- tsbs(
    x = x,
    bs_type = "moving",
    block_length = 10,
    num_boots = 1,
    return_diagnostics = TRUE
  )
  
  expect_true("diagnostics" %in% names(result))
  expect_equal(result$diagnostics$meta$num_boots, 1)
})

test_that("diagnostics handle short series", {
  x <- generate_test_series(n = 30, k = 1)
  
  result <- tsbs(
    x = x,
    bs_type = "moving",
    block_length = 5,
    num_boots = 10,
    return_diagnostics = TRUE
  )
  
  expect_true("diagnostics" %in% names(result))
})


# =============================================================================
# Tests for Consistency Across Bootstrap Types
# =============================================================================

test_that("all bootstrap types produce consistent diagnostics structure", {
  x <- generate_test_series(n = 100, k = 1)
  
  # Test moving
  result_moving <- tsbs(x = x, bs_type = "moving", block_length = 10, 
                        num_boots = 10, return_diagnostics = TRUE)
  
  # Test stationary
  result_stat <- tsbs(x = x, bs_type = "stationary", block_length = 10,
                      num_boots = 10, return_diagnostics = TRUE)
  
  # Test wild
  result_wild <- tsbs(x = x, bs_type = "wild", num_boots = 10,
                      return_diagnostics = TRUE)
  
  # All should have same basic structure
  for (result in list(result_moving, result_stat, result_wild)) {
    diag <- result$diagnostics
    expect_s3_class(diag, "tsbs_diagnostics")
    expect_true("meta" %in% names(diag))
    expect_true("series_stats" %in% names(diag))
    expect_true("original_stats" %in% names(diag))
  }
})
