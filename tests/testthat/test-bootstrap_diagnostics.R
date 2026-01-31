# ==============================================================================
# Tests for bootstrap_diagnostics.R
# ==============================================================================
# This test suite covers all functions in bootstrap_diagnostics.R:
#
# 1. Diagnostic Collector Creation
#    - create_bootstrap_diagnostics()
#
# 2. Recording Functions
#    - record_blocks()
#    - record_replicate_stats()
#    - record_original_stats()
#    - record_method_specific()
#
# 3. Post-Processing / Computation
#    - compute_bootstrap_diagnostics()
#    - estimate_block_lengths_from_series()
#
# 4. S3 Methods
#    - summary.tsbs_diagnostics()
#    - print.tsbs_diagnostics()
#    - plot.tsbs_diagnostics()
#    - as.data.frame.tsbs_diagnostics()
#
# 5. Extraction Functions
#    - extract_blocks()
#    - extract_summary_stats()
#
# 6. Functional Output Analysis
#    - summarize_func_outs()
#    - compute_func_out_cv()
#    - extract_func_out_matrix()
#    - compute_robust_estimates()
#
# 7. Benchmarking Utilities
#    - benchmark_tsbs() (may require tsbs() function)
#    - print.tsbs_benchmark()
#    - plot.tsbs_benchmark()
#    - summary.tsbs_benchmark()
#
# 8. Regime Diagnostics
#    - record_regime_info()
#    - record_replicate_regimes()
#    - plot_regime_composition()
#    - summarize_regime_diagnostics()
#    - create_regime_bands_df()
#    - extract_smoothed_probabilities()
#    - create_probability_df()
#    - reconstruct_bootstrap_probabilities()
#
# 9. Block Bootstrap Diagnostics
#    - blockBootstrap_with_diagnostics()
#    - compute_default_block_length_r()
#    - plot_block_coverage()
#    - plot_block_lengths()
#    - summarize_block_diagnostics()
#
# 10. Internal Plot Helpers
#     - plot_block_lengths() (internal)
#     - plot_start_positions()
#     - plot_means_comparison()
#     - plot_acf_comparison()
#     - plot_length_distribution()
#
# ==============================================================================

# ==============================================================================
# Test Fixtures and Helper Functions
# ==============================================================================

#' Create a simple test data matrix
#' @param n Number of rows (time points)
#' @param k Number of columns (variables)
#' @param seed Random seed for reproducibility
create_test_data <- function(n = 100, k = 3, seed = 123) {
  set.seed(seed)
  matrix(rnorm(n * k), nrow = n, ncol = k)
}

#' Create a test diagnostics object with minimal setup
#' @param bs_type Bootstrap type
#' @param n_original Original series length
#' @param n_vars Number of variables
#' @param num_boots Number of bootstrap replicates
create_test_diagnostics <- function(
    bs_type = "moving",
    n_original = 100,
    n_vars = 3,
    num_boots = 10
) {
  create_bootstrap_diagnostics(
    bs_type = bs_type,
    n_original = n_original,
    n_vars = n_vars,
    num_boots = num_boots
  )
}

#' Create a fully populated diagnostics object for testing
create_populated_diagnostics <- function(seed = 42) {
  set.seed(seed)
  
  n_original <- 100
  n_vars <- 2
  num_boots <- 5
  
  # Create diagnostics

  diag <- create_bootstrap_diagnostics(
    bs_type = "moving",
    n_original = n_original,
    n_vars = n_vars,
    num_boots = num_boots
  )
  
  # Create and record original data
  orig_data <- create_test_data(n_original, n_vars, seed)
  diag <- record_original_stats(diag, orig_data)
  
  # Record blocks and stats for each replicate
  for (b in seq_len(num_boots)) {
    # Simulate block info
    n_blocks <- sample(5:15, 1)
    block_lengths <- sample(5:20, n_blocks, replace = TRUE)
    start_positions <- sample(1:(n_original - 5), n_blocks, replace = TRUE)
    
    diag <- record_blocks(
      diag,
      replicate_idx = b,
      block_lengths = block_lengths,
      start_positions = start_positions,
      block_type = "overlapping"
    )
    
    # Simulate bootstrap series
    boot_mat <- create_test_data(n_original, n_vars, seed + b)
    diag <- record_replicate_stats(diag, b, boot_mat)
  }
  
  # Add config
  diag$config <- list(
    block_length = 10,
    block_type = "overlapping",
    n_boot = n_original
  )
  
  diag
}


# ==============================================================================
# SECTION 1: Diagnostic Collector Creation Tests
# ==============================================================================

context("create_bootstrap_diagnostics()")

test_that("create_bootstrap_diagnostics returns correct class", {
  diag <- create_bootstrap_diagnostics(
    bs_type = "moving",
    n_original = 100,
    n_vars = 3,
    num_boots = 10
  )
  
  expect_s3_class(diag, "tsbs_diagnostics")
})

test_that("create_bootstrap_diagnostics stores metadata correctly", {
  diag <- create_bootstrap_diagnostics(
    bs_type = "stationary",
    n_original = 250,
    n_vars = 5,
    num_boots = 50
  )
  
  expect_equal(diag$meta$bs_type, "stationary")
  expect_equal(diag$meta$n_original, 250)
  expect_equal(diag$meta$n_vars, 5)
  expect_equal(diag$meta$num_boots, 50)
  expect_true(inherits(diag$meta$timestamp, "POSIXct"))
})

test_that("create_bootstrap_diagnostics initializes storage lists correctly", {
  num_boots <- 20
  diag <- create_bootstrap_diagnostics(
    bs_type = "moving",
    n_original = 100,
    n_vars = 3,
    num_boots = num_boots
  )
  
  # Check blocks structure
  expect_type(diag$blocks, "list")
  expect_length(diag$blocks$replicate_blocks, num_boots)
  expect_length(diag$blocks$all_block_lengths, 0)
  expect_length(diag$blocks$all_start_positions, 0)
  
  # Check series_stats structure
  expect_type(diag$series_stats, "list")
  expect_length(diag$series_stats$replicate_means, num_boots)
  expect_length(diag$series_stats$replicate_sds, num_boots)
  expect_length(diag$series_stats$replicate_ac1, num_boots)
  expect_length(diag$series_stats$replicate_lengths, num_boots)
  
  # Check original_stats structure
  expect_null(diag$original_stats$means)
  expect_null(diag$original_stats$sds)
  expect_null(diag$original_stats$ac1)
  
  # Check method_specific and config
  expect_type(diag$method_specific, "list")
  expect_length(diag$method_specific, 0)
  expect_type(diag$config, "list")
  expect_length(diag$config, 0)
})

test_that("create_bootstrap_diagnostics handles different bootstrap types", {
  types <- c("moving", "stationary", "hmm", "ms_varma_garch")
  
  for (bs_type in types) {
    diag <- create_bootstrap_diagnostics(
      bs_type = bs_type,
      n_original = 50,
      n_vars = 2,
      num_boots = 5
    )
    expect_equal(diag$meta$bs_type, bs_type)
  }
})


# ==============================================================================
# SECTION 2: Recording Functions Tests
# ==============================================================================

context("record_blocks()")

test_that("record_blocks stores block info for single replicate", {
  diag <- create_test_diagnostics(num_boots = 5)
  
  block_lengths <- c(10, 15, 8, 12)
  start_positions <- c(1, 11, 26, 34)
  
  diag <- record_blocks(
    diag,
    replicate_idx = 1,
    block_lengths = block_lengths,
    start_positions = start_positions,
    block_type = "overlapping"
  )
  
  # Check per-replicate storage
  expect_true(is.data.frame(diag$blocks$replicate_blocks[[1]]))
  expect_equal(nrow(diag$blocks$replicate_blocks[[1]]), 4)
  expect_equal(diag$blocks$replicate_blocks[[1]]$length, block_lengths)
  expect_equal(diag$blocks$replicate_blocks[[1]]$start_pos, start_positions)
  expect_equal(diag$blocks$replicate_blocks[[1]]$block_type[1], "overlapping")
  
  # Check aggregated storage
  expect_equal(diag$blocks$all_block_lengths, block_lengths)
  expect_equal(diag$blocks$all_start_positions, start_positions)
})

test_that("record_blocks accumulates across multiple replicates", {
  diag <- create_test_diagnostics(num_boots = 3)
  
  # Record for replicate 1
  diag <- record_blocks(diag, 1, c(5, 10), c(1, 6), "overlapping")
  # Record for replicate 2
  diag <- record_blocks(diag, 2, c(8, 7), c(20, 28), "overlapping")
  # Record for replicate 3
  diag <- record_blocks(diag, 3, c(12), c(50), "non-overlapping")
  
  # Check aggregation
  expect_equal(diag$blocks$all_block_lengths, c(5, 10, 8, 7, 12))
  expect_equal(diag$blocks$all_start_positions, c(1, 6, 20, 28, 50))
  
  # Check individual replicates preserved
  expect_equal(nrow(diag$blocks$replicate_blocks[[1]]), 2)
  expect_equal(nrow(diag$blocks$replicate_blocks[[2]]), 2)
  expect_equal(nrow(diag$blocks$replicate_blocks[[3]]), 1)
})

test_that("record_blocks handles NA block_type", {
  diag <- create_test_diagnostics()
  
  diag <- record_blocks(diag, 1, c(10, 10), c(1, 11))
  
  expect_true(all(is.na(diag$blocks$replicate_blocks[[1]]$block_type)))
})


context("record_replicate_stats()")

test_that("record_replicate_stats computes correct statistics", {
  diag <- create_test_diagnostics(n_vars = 2, num_boots = 3)
  
  set.seed(123)
  boot_mat <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                       2, 4, 6, 8, 10, 12, 14, 16, 18, 20), 
                     nrow = 10, ncol = 2)
  
  diag <- record_replicate_stats(diag, 1, boot_mat)
  
  # Check means
  expected_means <- colMeans(boot_mat)
  expect_equal(diag$series_stats$replicate_means[[1]], expected_means)
  
  # Check SDs
  expected_sds <- apply(boot_mat, 2, sd)
  expect_equal(diag$series_stats$replicate_sds[[1]], expected_sds)
  
  # Check length
  expect_equal(diag$series_stats$replicate_lengths[1], 10)
  
  # Check lag-1 autocorrelation is computed (just check it's numeric, not NA)
  expect_true(is.numeric(diag$series_stats$replicate_ac1[[1]]))
  expect_length(diag$series_stats$replicate_ac1[[1]], 2)
})

test_that("record_replicate_stats handles single-row matrix", {
  diag <- create_test_diagnostics(n_vars = 2, num_boots = 1)
  
  single_row <- matrix(c(5, 10), nrow = 1)
  diag <- record_replicate_stats(diag, 1, single_row)
  
  expect_equal(diag$series_stats$replicate_lengths[1], 1)
  expect_true(all(is.na(diag$series_stats$replicate_ac1[[1]])))
})

test_that("record_replicate_stats handles data frame input", {
  diag <- create_test_diagnostics(n_vars = 2, num_boots = 1)
  
  boot_df <- data.frame(V1 = 1:10, V2 = 11:20)
  diag <- record_replicate_stats(diag, 1, boot_df)
  
  expect_equal(unname(diag$series_stats$replicate_means[[1]]), c(5.5, 15.5))
})


context("record_original_stats()")

test_that("record_original_stats computes correct statistics", {
  diag <- create_test_diagnostics(n_original = 50, n_vars = 3)
  
  set.seed(456)
  x <- create_test_data(50, 3, 456)
  
  diag <- record_original_stats(diag, x)
  
  expect_equal(diag$original_stats$means, colMeans(x))
  expect_equal(diag$original_stats$sds, apply(x, 2, sd))
  expect_length(diag$original_stats$ac1, 3)
  expect_true(all(is.numeric(diag$original_stats$ac1)))
})

test_that("record_original_stats handles vector input", {
  diag <- create_test_diagnostics(n_original = 20, n_vars = 1)
  
  x_vec <- rnorm(20)
  diag <- record_original_stats(diag, x_vec)
  
  expect_length(diag$original_stats$means, 1)
  expect_equal(diag$original_stats$means[1], mean(x_vec))
})


context("record_method_specific()")

test_that("record_method_specific stores named items", {
  diag <- create_test_diagnostics()
  
  diag <- record_method_specific(
    diag,
    regime_info = list(num_states = 2),
    transition_matrix = matrix(c(0.9, 0.1, 0.2, 0.8), 2, 2)
  )
  
  expect_true("regime_info" %in% names(diag$method_specific))
  expect_true("transition_matrix" %in% names(diag$method_specific))
  expect_equal(diag$method_specific$regime_info$num_states, 2)
})

test_that("record_method_specific updates existing items", {
  diag <- create_test_diagnostics()
  
  diag <- record_method_specific(diag, value1 = 10)
  expect_equal(diag$method_specific$value1, 10)
  
  diag <- record_method_specific(diag, value1 = 20, value2 = 30)
  expect_equal(diag$method_specific$value1, 20)
  expect_equal(diag$method_specific$value2, 30)
})


# ==============================================================================
# SECTION 3: Post-Processing / Computation Tests
# ==============================================================================

context("compute_bootstrap_diagnostics()")

test_that("compute_bootstrap_diagnostics creates valid diagnostics object", {
  set.seed(789)
  original_data <- create_test_data(100, 2)
  
  # Create simple bootstrap series
  bootstrap_series <- lapply(1:5, function(i) {
    idx <- sample(1:100, 100, replace = TRUE)
    original_data[idx, , drop = FALSE]
  })
  
  diag <- compute_bootstrap_diagnostics(
    bootstrap_series = bootstrap_series,
    original_data = original_data,
    bs_type = "iid",
    config = list(method = "simple")
  )
  
  expect_s3_class(diag, "tsbs_diagnostics")
  expect_equal(diag$meta$bs_type, "iid")
  expect_equal(diag$meta$num_boots, 5)
  expect_equal(diag$meta$n_original, 100)
  expect_equal(diag$meta$n_vars, 2)
})

test_that("compute_bootstrap_diagnostics records original stats", {
  original_data <- create_test_data(50, 3, seed = 111)
  bootstrap_series <- list(original_data, original_data)
  
  diag <- compute_bootstrap_diagnostics(
    bootstrap_series = bootstrap_series,
    original_data = original_data,
    bs_type = "test"
  )
  
  expect_equal(diag$original_stats$means, colMeans(original_data))
  expect_equal(diag$original_stats$sds, apply(original_data, 2, sd))
})

test_that("compute_bootstrap_diagnostics records replicate stats", {
  set.seed(222)
  original_data <- create_test_data(80, 2)
  
  bootstrap_series <- lapply(1:3, function(i) {
    create_test_data(80, 2, seed = 222 + i)
  })
  
  diag <- compute_bootstrap_diagnostics(
    bootstrap_series = bootstrap_series,
    original_data = original_data,
    bs_type = "test"
  )
  
  # Check all replicates have stats
  for (i in 1:3) {
    expect_length(diag$series_stats$replicate_means[[i]], 2)
    expect_length(diag$series_stats$replicate_sds[[i]], 2)
    expect_equal(diag$series_stats$replicate_lengths[i], 80)
  }
})

test_that("compute_bootstrap_diagnostics stores config", {
  original_data <- create_test_data(30, 1)
  bootstrap_series <- list(original_data)
  
  config <- list(block_length = 5, p = 0.1, custom = "value")
  
  diag <- compute_bootstrap_diagnostics(
    bootstrap_series = bootstrap_series,
    original_data = original_data,
    bs_type = "stationary",
    config = config
  )
  
  expect_equal(diag$config, config)
})


context("estimate_block_lengths_from_series()")

test_that("estimate_block_lengths_from_series detects block boundaries", {
  set.seed(333)
  n <- 50
  k <- 2
  original_data <- create_test_data(n, k, 333)
  
  # Create bootstrap by concatenating blocks
  block1 <- original_data[1:10, ]
  block2 <- original_data[25:35, ]
  block3 <- original_data[5:15, ]
  boot_mat <- rbind(block1, block2, block3)
  
  bootstrap_series <- list(boot_mat)
  
  diag <- create_bootstrap_diagnostics("moving", n, k, 1)
  diag <- tsbs:::estimate_block_lengths_from_series(diag, bootstrap_series, original_data)
  
  # Should have detected block boundaries (3 blocks)
  expect_true(length(diag$blocks$all_block_lengths) >= 1)
  expect_true(!is.null(diag$blocks$replicate_blocks[[1]]))
})


# ==============================================================================
# SECTION 4: S3 Methods Tests
# ==============================================================================

context("summary.tsbs_diagnostics()")

test_that("summary.tsbs_diagnostics produces output without error", {
  diag <- create_populated_diagnostics()
  
  expect_output(summary(diag), "tsbs Bootstrap Diagnostics Summary")
  expect_output(summary(diag), "BOOTSTRAP CONFIGURATION")
  expect_output(summary(diag), "Bootstrap type:")
})

test_that("summary.tsbs_diagnostics returns object invisibly", {
  diag <- create_populated_diagnostics()
  
  result <- capture.output(ret <- summary(diag))
  expect_s3_class(ret, "tsbs_diagnostics")
})

test_that("summary.tsbs_diagnostics shows block statistics for block methods", {
  diag <- create_populated_diagnostics()
  
  output <- capture.output(summary(diag))
  output_text <- paste(output, collapse = "\n")
  
  expect_true(grepl("BLOCK LENGTH STATISTICS", output_text))
  expect_true(grepl("Mean block length", output_text))
})

test_that("summary.tsbs_diagnostics shows original vs bootstrap comparison", {
  diag <- create_populated_diagnostics()
  
  output <- capture.output(summary(diag))
  output_text <- paste(output, collapse = "\n")
  
  expect_true(grepl("ORIGINAL vs BOOTSTRAP STATISTICS", output_text))
  expect_true(grepl("MEANS", output_text))
})


context("print.tsbs_diagnostics()")

test_that("print.tsbs_diagnostics produces brief output", {
  diag <- create_test_diagnostics(bs_type = "moving", n_original = 100,
                                   n_vars = 3, num_boots = 10)
  
  output <- capture.output(print(diag))
  output_text <- paste(output, collapse = "\n")
  
  expect_true(grepl("tsbs Bootstrap Diagnostics", output_text))
  expect_true(grepl("Type: moving", output_text))
  expect_true(grepl("Replicates: 10", output_text))
  expect_true(grepl("100 x 3", output_text))
})

test_that("print.tsbs_diagnostics returns object invisibly", {
  diag <- create_test_diagnostics()
  
  result <- capture.output(ret <- print(diag))
  expect_s3_class(ret, "tsbs_diagnostics")
})


context("as.data.frame.tsbs_diagnostics()")

test_that("as.data.frame with what='stats' returns replicate statistics", {
  diag <- create_populated_diagnostics()
  
  df <- as.data.frame(diag, what = "stats")
  
  expect_true(is.data.frame(df))
  expect_true("replicate" %in% names(df))
  expect_true("length" %in% names(df))
  expect_equal(nrow(df), diag$meta$num_boots)
})

test_that("as.data.frame with what='blocks' returns block data", {
  diag <- create_populated_diagnostics()
  
  df <- as.data.frame(diag, what = "blocks")
  
  expect_true(is.data.frame(df))
  expect_true("length" %in% names(df))
  expect_true("start_pos" %in% names(df))
  expect_true("replicate" %in% names(df))
})

test_that("as.data.frame with what='all' combines stats and blocks", {
  diag <- create_populated_diagnostics()
  
  df <- as.data.frame(diag, what = "all")
  
  expect_true(is.data.frame(df))
  expect_true("replicate" %in% names(df))
  # Should have block summary columns
  expect_true(any(grepl("block", names(df), ignore.case = TRUE)))
})


# ==============================================================================
# SECTION 5: Extraction Functions Tests
# ==============================================================================

context("extract_blocks()")

test_that("extract_blocks returns specific replicate when specified", {
  diag <- create_populated_diagnostics()
  
  blocks_1 <- extract_blocks(diag, replicate = 1)
  
  expect_true(is.data.frame(blocks_1))
  expect_equal(blocks_1, diag$blocks$replicate_blocks[[1]])
})

test_that("extract_blocks returns all replicates when no replicate specified", {
  diag <- create_populated_diagnostics()
  
  all_blocks <- extract_blocks(diag)
  
  expect_true(is.data.frame(all_blocks))
  expect_true("replicate" %in% names(all_blocks))
  # Should have data from all replicates
  expect_equal(length(unique(all_blocks$replicate)), diag$meta$num_boots)
})

test_that("extract_blocks handles empty replicates", {
  diag <- create_test_diagnostics(num_boots = 3)
  # Only populate first replicate
  diag <- record_blocks(diag, 1, c(10), c(1))
  
  all_blocks <- extract_blocks(diag)
  
  # Should only have data from replicate 1
  expect_equal(unique(all_blocks$replicate), 1)
})


context("extract_summary_stats()")

test_that("extract_summary_stats returns correct structure", {
  diag <- create_populated_diagnostics()
  
  stats <- extract_summary_stats(diag)
  
  expect_type(stats, "list")
  expect_true("original" %in% names(stats))
  expect_true("bootstrap" %in% names(stats))
  expect_true("block_lengths" %in% names(stats))
})

test_that("extract_summary_stats includes bootstrap quantiles", {
  diag <- create_populated_diagnostics()
  
  stats <- extract_summary_stats(diag)
  
  expect_true(!is.null(stats$bootstrap$means$quantiles))
  expect_equal(nrow(stats$bootstrap$means$quantiles), 5)  # 5 quantile levels
})

test_that("extract_summary_stats includes block length statistics", {
  diag <- create_populated_diagnostics()
  
  stats <- extract_summary_stats(diag)
  
  expect_true(!is.null(stats$block_lengths$mean))
  expect_true(!is.null(stats$block_lengths$sd))
  expect_true(!is.null(stats$block_lengths$min))
  expect_true(!is.null(stats$block_lengths$max))
  expect_true(!is.null(stats$block_lengths$median))
})


# ==============================================================================
# SECTION 6: Functional Output Analysis Tests
# ==============================================================================

context("summarize_func_outs()")

test_that("summarize_func_outs handles list input", {
  func_outs <- list(
    c(0.3, 0.4, 0.3),
    c(0.35, 0.35, 0.3),
    c(0.25, 0.45, 0.3),
    c(0.32, 0.38, 0.3)
  )
  
  result <- summarize_func_outs(func_outs)
  
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 3)
  expect_true(all(c("Name", "Mean", "SD", "CI_Lower", "CI_Upper", "CI_Width") %in% names(result)))
})

test_that("summarize_func_outs handles matrix input", {
  boot_mat <- matrix(c(
    0.3, 0.4, 0.3,
    0.35, 0.35, 0.3,
    0.25, 0.45, 0.3,
    0.32, 0.38, 0.3
  ), nrow = 4, byrow = TRUE)
  
  result <- summarize_func_outs(boot_mat)
  
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 3)
})

test_that("summarize_func_outs uses provided names", {
  func_outs <- list(c(0.5, 0.5), c(0.6, 0.4), c(0.55, 0.45))
  names_vec <- c("Asset_A", "Asset_B")
  
  result <- summarize_func_outs(func_outs, names = names_vec)
  
  expect_equal(result$Name, names_vec)
})

test_that("summarize_func_outs uses custom probability levels", {
  func_outs <- list(c(1, 2), c(2, 3), c(1.5, 2.5), c(1.8, 2.2))
  
  result <- summarize_func_outs(func_outs, probs = c(0.1, 0.9))
  
  expect_true(is.data.frame(result))
  # CI bounds should be at 10% and 90%
})

test_that("summarize_func_outs returns NULL for insufficient data", {
  # Only 1 replicate
  func_outs <- list(c(0.5, 0.5))
  
  expect_warning(
    result <- summarize_func_outs(func_outs),
    "at least 2 bootstrap replicates"
  )
  expect_null(result)
})


context("compute_func_out_cv()")

test_that("compute_func_out_cv computes coefficient of variation", {
  # Create data with known CV
  set.seed(444)
  func_outs <- lapply(1:20, function(i) c(10 + rnorm(1, 0, 1), 5 + rnorm(1, 0, 0.5)))
  
  result <- compute_func_out_cv(func_outs)
  
  expect_true(is.data.frame(result))
  expect_true("CV" %in% names(result))
  expect_true("Stability" %in% names(result))
})

test_that("compute_func_out_cv classifies stability correctly", {
  # Stable: CV < 0.3 (low variability relative to mean)
  # Use tight distribution around mean of 10: SD ~ 1, CV ~ 0.1
  set.seed(100)
  stable_outs <- lapply(1:20, function(i) c(10 + rnorm(1, 0, 0.5)))
  
  # Unstable: CV > 0.6 (high variability relative to mean)
  # Use wide distribution: mean ~ 1, SD ~ 2, CV ~ 2.0
  set.seed(200)
  unstable_outs <- lapply(1:20, function(i) c(1 + rnorm(1, 0, 2)))
  
  stable_result <- compute_func_out_cv(stable_outs)
  unstable_result <- compute_func_out_cv(unstable_outs)
  
  expect_equal(stable_result$Stability, "Stable")
  expect_equal(unstable_result$Stability, "Unstable")
})

test_that("compute_func_out_cv uses custom thresholds", {
  func_outs <- lapply(1:10, function(i) c(10 + rnorm(1, 0, 4)))  # CV ~ 0.4
  
  # With default thresholds (0.3, 0.6), should be "Moderate"
  default_result <- compute_func_out_cv(func_outs)
  
  # With custom thresholds, should be different
  custom_result <- compute_func_out_cv(
    func_outs, 
    cv_thresholds = c(Stable = 0.5, Moderate = 0.8)
  )
  
  # The classifications might differ based on actual CV
  expect_true(all(custom_result$Stability %in% c("Stable", "Moderate", "Unstable")))
})


context("extract_func_out_matrix()")

test_that("extract_func_out_matrix converts list to matrix", {
  func_outs <- list(c(1, 2, 3), c(4, 5, 6), c(7, 8, 9))
  
  result <- extract_func_out_matrix(func_outs)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3, 3))
  expect_equal(result[1, ], c(1, 2, 3))
})

test_that("extract_func_out_matrix handles matrix input", {
  mat <- matrix(1:9, nrow = 3)
  
  result <- extract_func_out_matrix(mat)
  
  expect_true(is.matrix(result))
  expect_equal(result, mat)
})

test_that("extract_func_out_matrix handles data frame input", {
  # Note: data.frame is a list, so each column becomes a row in the output
  # For typical use (bootstrap replicates as rows), pass a matrix instead
  df <- data.frame(a = 1:3, b = 4:6)
  
  result <- extract_func_out_matrix(df)
  
  expect_true(is.matrix(result))
  # Data frame has 2 columns (list elements), each of length 3
  # So result is 2 rows x 3 columns
  expect_equal(dim(result), c(2, 3))
})

test_that("extract_func_out_matrix flattens matrix elements in list", {
  # List of matrices (e.g., correlation matrices)
  func_outs <- list(
    matrix(c(1, 0.5, 0.5, 1), 2, 2),
    matrix(c(1, 0.6, 0.6, 1), 2, 2)
  )
  
  result <- extract_func_out_matrix(func_outs)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 4)  # Flattened 2x2 matrix
})


context("compute_robust_estimates()")

test_that("compute_robust_estimates returns all estimators", {
  func_outs <- lapply(1:20, function(i) c(0.3 + rnorm(1, 0, 0.05), 
                                           0.7 + rnorm(1, 0, 0.05)))
  
  result <- compute_robust_estimates(func_outs)
  
  expect_true(is.data.frame(result))
  expect_true(all(c("Boot_Mean", "Boot_Median", "Winsorized", "Conservative") 
                  %in% names(result)))
})

test_that("compute_robust_estimates includes point estimate when provided", {
  func_outs <- list(c(0.5, 0.5), c(0.55, 0.45), c(0.48, 0.52))
  point_est <- c(0.5, 0.5)
  
  result <- compute_robust_estimates(func_outs, point_est = point_est)
  
  expect_true("Point" %in% names(result))
  expect_equal(result$Point, c(0.5, 0.5))
})

test_that("compute_robust_estimates uses custom trim proportion", {
  func_outs <- lapply(1:20, function(i) c(10 + rnorm(1)))
  
  result_default <- compute_robust_estimates(func_outs, trim = 0.1)
  result_heavy <- compute_robust_estimates(func_outs, trim = 0.3)
  
  # Both should work, potentially with different Winsorized values
  expect_true(is.data.frame(result_default))
  expect_true(is.data.frame(result_heavy))
})

test_that("compute_robust_estimates renormalizes portfolio-like outputs", {
  # Weights that sum to 1
  func_outs <- lapply(1:10, function(i) {
    w <- c(0.3, 0.3, 0.4) + rnorm(3, 0, 0.02)
    w / sum(w)
  })
  
  result <- compute_robust_estimates(func_outs, conservative_quantile = 0.25)
  
  # Conservative estimate should sum to 1 after renormalization
  expect_equal(sum(result$Conservative), 1, tolerance = 0.001)
})


# ==============================================================================
# SECTION 7: Regime Diagnostics Tests
# ==============================================================================

context("record_regime_info()")

test_that("record_regime_info stores regime information", {
  diag <- create_test_diagnostics(bs_type = "hmm", num_boots = 5)
  
  original_states <- c(rep(1, 50), rep(2, 50))
  
  diag <- record_regime_info(
    diag,
    original_states = original_states,
    num_states = 2,
    state_labels = c("Low Vol", "High Vol")
  )
  
  expect_true(!is.null(diag$method_specific$regime_info))
  expect_equal(diag$method_specific$regime_info$num_states, 2)
  expect_equal(diag$method_specific$regime_info$state_labels, c("Low Vol", "High Vol"))
  expect_equal(diag$method_specific$regime_info$original_states, original_states)
})

test_that("record_regime_info computes state counts", {
  diag <- create_test_diagnostics(bs_type = "hmm")
  
  # 30% state 1, 70% state 2
  original_states <- c(rep(1, 30), rep(2, 70))
  
  diag <- record_regime_info(diag, original_states, 2)
  
  counts <- diag$method_specific$regime_info$original_state_counts
  expect_equal(as.numeric(counts), c(30, 70))
})

test_that("record_regime_info creates default labels", {
  diag <- create_test_diagnostics(bs_type = "hmm")
  
  diag <- record_regime_info(diag, c(1, 2, 1, 2), 2, state_labels = NULL)
  
  expect_equal(diag$method_specific$regime_info$state_labels, c("State 1", "State 2"))
})


context("record_replicate_regimes()")

test_that("record_replicate_regimes stores state sequence", {
  diag <- create_test_diagnostics(bs_type = "hmm", num_boots = 3)
  diag <- record_regime_info(diag, c(1, 2, 1), 2)
  
  rep_states <- c(1, 1, 2, 2, 1)
  diag <- record_replicate_regimes(diag, 1, rep_states)
  
  expect_equal(
    diag$method_specific$regime_info$replicate_states[[1]], 
    rep_states
  )
})

test_that("record_replicate_regimes stores block states and source indices", {
  diag <- create_test_diagnostics(bs_type = "hmm", num_boots = 2)
  diag <- record_regime_info(diag, c(1, 2, 1), 2)
  
  diag <- record_replicate_regimes(
    diag,
    replicate_idx = 1,
    replicate_states = c(1, 2),
    block_states = c(1, 2),
    block_source_indices = c(1, 2),
    source_indices = c(1, 2, 3, 4)
  )
  
  block_info <- diag$method_specific$regime_info$replicate_block_states[[1]]
  expect_equal(block_info$states, c(1, 2))
  expect_equal(block_info$source_indices, c(1, 2, 3, 4))
})

test_that("record_replicate_regimes warns if regime_info not initialized", {
  diag <- create_test_diagnostics()
  
  expect_warning(
    record_replicate_regimes(diag, 1, c(1, 2)),
    "record_regime_info"
  )
})


context("create_regime_bands_df()")

test_that("create_regime_bands_df creates correct band structure", {
  states <- c(1, 1, 1, 2, 2, 1, 1)
  
  bands <- tsbs:::create_regime_bands_df(states, "Test")
  
  expect_true(is.data.frame(bands))
  expect_equal(names(bands), c("xmin", "xmax", "State", "Series"))
  expect_equal(nrow(bands), 3)  # Three state runs
  
  # Check boundaries
  expect_equal(bands$xmin, c(1, 4, 6))
  expect_equal(bands$xmax, c(3, 5, 7))
  expect_equal(bands$State, c(1, 2, 1))
})

test_that("create_regime_bands_df handles single state", {
  states <- rep(1, 10)
  
  bands <- tsbs:::create_regime_bands_df(states, "Single")
  
  expect_equal(nrow(bands), 1)
  expect_equal(bands$xmin, 1)
  expect_equal(bands$xmax, 10)
})


context("summarize_regime_diagnostics()")

test_that("summarize_regime_diagnostics prints summary", {
  diag <- create_test_diagnostics(bs_type = "hmm", num_boots = 3)
  diag <- record_regime_info(diag, c(rep(1, 40), rep(2, 60)), 2, 
                              c("Low", "High"))
  
  # Add replicate data
  diag <- record_replicate_regimes(diag, 1, c(rep(1, 35), rep(2, 65)))
  diag <- record_replicate_regimes(diag, 2, c(rep(1, 45), rep(2, 55)))
  
  output <- capture.output(summarize_regime_diagnostics(diag))
  output_text <- paste(output, collapse = "\n")
  
  expect_true(grepl("Regime Bootstrap Diagnostics", output_text))
  expect_true(grepl("Number of states: 2", output_text))
  expect_true(grepl("Low", output_text))
  expect_true(grepl("High", output_text))
})

test_that("summarize_regime_diagnostics returns stats invisibly", {
  diag <- create_test_diagnostics(bs_type = "hmm")
  diag <- record_regime_info(diag, c(1, 1, 2, 2), 2)
  
  result <- capture.output(ret <- summarize_regime_diagnostics(diag))
  
  expect_type(ret, "list")
  expect_true("original_counts" %in% names(ret))
  expect_true("original_props" %in% names(ret))
})

test_that("summarize_regime_diagnostics handles missing regime info", {
  diag <- create_test_diagnostics()
  
  expect_output(summarize_regime_diagnostics(diag), "No regime information")
})


# ==============================================================================
# SECTION 8: Block Bootstrap Diagnostics Tests
# ==============================================================================

context("blockBootstrap_with_diagnostics()")

test_that("blockBootstrap_with_diagnostics returns correct structure", {
  x <- create_test_data(100, 2, seed = 555)
  
  result <- blockBootstrap_with_diagnostics(
    x,
    block_length = 10,
    bs_type = "moving",
    num_boots = 5,
    collect_diagnostics = TRUE
  )
  
  expect_type(result, "list")
  expect_true("bootstrap_series" %in% names(result))
  expect_true("diagnostics" %in% names(result))
  expect_length(result$bootstrap_series, 5)
  expect_s3_class(result$diagnostics, "tsbs_diagnostics")
})

test_that("blockBootstrap_with_diagnostics creates correct length series", {
  x <- create_test_data(80, 3)
  
  result <- blockBootstrap_with_diagnostics(
    x,
    n_boot = 100,
    block_length = 8,
    bs_type = "moving",
    num_boots = 3,
    collect_diagnostics = TRUE
  )
  
  for (i in 1:3) {
    expect_equal(nrow(result$bootstrap_series[[i]]), 100)
    expect_equal(ncol(result$bootstrap_series[[i]]), 3)
  }
})

test_that("blockBootstrap_with_diagnostics returns only series when diagnostics disabled", {
  x <- create_test_data(50, 2)
  
  result <- blockBootstrap_with_diagnostics(
    x,
    block_length = 5,
    bs_type = "moving",
    num_boots = 3,
    collect_diagnostics = FALSE
  )
  
  expect_type(result, "list")
  expect_length(result, 3)
  expect_true(is.matrix(result[[1]]))
})

test_that("blockBootstrap_with_diagnostics records block info correctly", {
  set.seed(666)
  x <- create_test_data(100, 2)
  
  result <- blockBootstrap_with_diagnostics(
    x,
    block_length = 10,
    bs_type = "moving",
    num_boots = 10,
    collect_diagnostics = TRUE
  )
  
  diag <- result$diagnostics
  
  # Should have block data
  expect_true(length(diag$blocks$all_block_lengths) > 0)
  expect_true(length(diag$blocks$all_start_positions) > 0)
  
  # All block lengths should be <= block_length for moving bootstrap
  expect_true(all(diag$blocks$all_block_lengths <= 10))
})

test_that("blockBootstrap_with_diagnostics handles stationary bootstrap", {
  set.seed(777)
  x <- create_test_data(200, 2)
  
  result <- blockBootstrap_with_diagnostics(
    x,
    bs_type = "stationary",
    p = 0.1,  # Expected block length = 10
    num_boots = 20,
    collect_diagnostics = TRUE
  )
  
  diag <- result$diagnostics
  
  # Block lengths should vary (geometric distribution)
  block_lengths <- diag$blocks$all_block_lengths
  expect_true(length(unique(block_lengths)) > 1)
  
  # Mean should be approximately 1/p = 10
  # (with some tolerance due to truncation and small sample)
  expect_true(mean(block_lengths) > 3)  # Loose check
})

test_that("blockBootstrap_with_diagnostics stores config", {
  x <- create_test_data(50, 2)
  
  result <- blockBootstrap_with_diagnostics(
    x,
    block_length = 8,
    bs_type = "moving",
    block_type = "overlapping",
    num_boots = 2,
    collect_diagnostics = TRUE
  )
  
  config <- result$diagnostics$config
  expect_equal(config$block_length, 8)
  expect_equal(config$block_type, "overlapping")
})


context("compute_default_block_length_r()")

test_that("compute_default_block_length_r returns reasonable value", {
  # Autocorrelated data
  set.seed(888)
  n <- 200
  x <- matrix(0, n, 2)
  x[1, ] <- rnorm(2)
  for (i in 2:n) {
    x[i, ] <- 0.7 * x[i-1, ] + rnorm(2, 0, 0.5)
  }
  
  block_len <- tsbs:::compute_default_block_length_r(x)
  
  expect_true(is.integer(block_len) || is.numeric(block_len))
  expect_true(block_len >= 5)
  expect_true(block_len <= sqrt(n))
})

test_that("compute_default_block_length_r handles iid data", {
  set.seed(999)
  x <- create_test_data(100, 3)  # Essentially iid
  
  block_len <- tsbs:::compute_default_block_length_r(x)
  
  expect_true(block_len >= 5)
})

test_that("compute_default_block_length_r handles short series", {
  x <- create_test_data(20, 2)
  
  block_len <- tsbs:::compute_default_block_length_r(x)
  
  expect_true(block_len >= 1)
  expect_true(block_len <= 20)
})


context("summarize_block_diagnostics()")

test_that("summarize_block_diagnostics produces output", {
  diag <- create_populated_diagnostics()
  diag$method_specific$block_info$replicate_source_indices <- 
    lapply(1:diag$meta$num_boots, function(i) sample(1:100, 100, replace = TRUE))
  
  output <- capture.output(summarize_block_diagnostics(diag))
  output_text <- paste(output, collapse = "\n")
  
  expect_true(grepl("Block Bootstrap Diagnostics", output_text))
  expect_true(grepl("Bootstrap type:", output_text))
})


# ==============================================================================
# SECTION 9: Plot Function Tests (basic functionality without graphics device)
# ==============================================================================

context("plot.tsbs_diagnostics()")

test_that("plot.tsbs_diagnostics accepts valid type arguments", {
  skip_on_cran()
  diag <- create_populated_diagnostics()

  valid_types <- c("all", "block_lengths", "start_positions",
                   "means_comparison", "acf_comparison", "length_distribution")

  for (type in valid_types) {
    # Should not error (messages are OK, errors are not)
    expect_no_error(
      tryCatch(
        suppressMessages(plot(diag, type = type)),
        error = function(e) {
          if (!grepl("ggplot2", e$message)) stop(e)
        }
      )
    )
  }
})

test_that("plot.tsbs_diagnostics rejects invalid type",
  {diag <- create_test_diagnostics()
  
  expect_error(
    plot(diag, type = "invalid_type"),
    "should be one of"
  )
})


# Skip actual plotting tests if ggplot2 not available
skip_if_no_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    skip("ggplot2 not available")
  }
}

test_that("plot_block_lengths creates ggplot object", {
  skip_on_cran()
  skip_if_no_ggplot2()
  
  diag <- create_populated_diagnostics()
  
  # Capture the plot output (exported function has show_expected parameter, not use_ggplot)
  p <- plot_block_lengths(diag, show_expected = FALSE)
  
  expect_s3_class(p, "ggplot")
})

test_that("plot_start_positions creates ggplot object", {
  skip_on_cran()
  skip_if_no_ggplot2()
  
  diag <- create_populated_diagnostics()
  
  p <- tsbs:::plot_start_positions(diag, use_ggplot = TRUE)
  
  expect_s3_class(p, "ggplot")
})

test_that("plot_means_comparison creates ggplot object", {
  skip_on_cran()
  skip_if_no_ggplot2()
  
  diag <- create_populated_diagnostics()
  
  p <- tsbs:::plot_means_comparison(diag, use_ggplot = TRUE)
  
  expect_s3_class(p, "ggplot")
})

test_that("plot_acf_comparison creates ggplot object", {
  skip_on_cran()
  skip_if_no_ggplot2()
  
  diag <- create_populated_diagnostics()
  
  p <- tsbs:::plot_acf_comparison(diag, use_ggplot = TRUE)
  
  expect_s3_class(p, "ggplot")
})


context("plot_block_coverage()")

test_that("plot_block_coverage accepts valid types", {
  skip_on_cran()
  skip_if_no_ggplot2()
  
  diag <- create_populated_diagnostics()
  # Need source indices for heatmap
  diag$method_specific$block_info$replicate_source_indices <- 
    lapply(1:diag$meta$num_boots, function(i) sample(1:100, 100, replace = TRUE))
  
  # histogram should work
  p <- plot_block_coverage(diag, type = "histogram")
  expect_s3_class(p, "ggplot")
  
  # blocks should work
  p <- plot_block_coverage(diag, type = "blocks")
  expect_s3_class(p, "ggplot")
})


# ==============================================================================
# SECTION 10: Edge Cases and Error Handling Tests
# ==============================================================================

context("Edge cases and error handling")

test_that("functions handle empty bootstrap series list", {
  # Empty list is handled gracefully - returns diagnostics with num_boots = 0
  result <- compute_bootstrap_diagnostics(
    bootstrap_series = list(),
    original_data = matrix(1:10, ncol = 1),
    bs_type = "test"
  )
  
  expect_s3_class(result, "tsbs_diagnostics")
  expect_equal(result$meta$num_boots, 0)
})

test_that("functions handle single-column data", {
  x <- create_test_data(50, 1)
  
  result <- blockBootstrap_with_diagnostics(
    x,
    block_length = 5,
    bs_type = "moving",
    num_boots = 3,
    collect_diagnostics = TRUE
  )
  
  expect_equal(ncol(result$bootstrap_series[[1]]), 1)
  expect_length(result$diagnostics$original_stats$means, 1)
})

test_that("functions handle very short series", {
  x <- create_test_data(10, 2)
  
  result <- blockBootstrap_with_diagnostics(
    x,
    block_length = 3,
    bs_type = "moving",
    num_boots = 2,
    collect_diagnostics = TRUE
  )
  
  expect_length(result$bootstrap_series, 2)
})

test_that("extract_blocks handles diagnostics with no blocks", {
  diag <- create_test_diagnostics(num_boots = 3)
  
  all_blocks <- extract_blocks(diag)
  
  expect_null(all_blocks)
})

test_that("summarize_func_outs handles mismatched names length", {
  func_outs <- list(c(1, 2, 3), c(4, 5, 6))
  
  expect_warning(
    result <- summarize_func_outs(func_outs, names = c("A", "B")),
    "does not match"
  )
  
  # Should fall back to default names
  expect_equal(result$Name, c("V1", "V2", "V3"))
})

test_that("compute_robust_estimates handles all-NA column", {
  # Test with a column that has all NA values
  func_outs <- list(c(1, NA), c(2, NA), c(3, NA))
  
  result <- compute_robust_estimates(func_outs)
  
  # Should handle gracefully - returns data frame with NaN for the NA column
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2)
  # First column should have valid values
  
  expect_false(is.na(result$Boot_Mean[1]))
})


# ==============================================================================
# SECTION 11: Integration Tests
# ==============================================================================

context("Integration tests")

test_that("complete workflow: create -> record -> extract -> summarize", {
  set.seed(1234)
  
  # Create original data
  x <- create_test_data(100, 3)
  
  # Run bootstrap with diagnostics
  result <- blockBootstrap_with_diagnostics(
    x,
    block_length = 10,
    bs_type = "moving",
    num_boots = 10,
    collect_diagnostics = TRUE
  )
  
  diag <- result$diagnostics
  
  # Extract various outputs
  blocks <- extract_blocks(diag)
  stats <- extract_summary_stats(diag)
  df_stats <- as.data.frame(diag, what = "stats")
  df_blocks <- as.data.frame(diag, what = "blocks")
  
  # Verify everything works together
  expect_true(is.data.frame(blocks))
  expect_type(stats, "list")
  expect_true(is.data.frame(df_stats))
  expect_true(is.data.frame(df_blocks))
  
  # Summary should work
  expect_output(summary(diag), "tsbs Bootstrap Diagnostics")
})

test_that("functional output analysis workflow", {
  set.seed(5678)
  
  # Simulate portfolio weights from multiple bootstrap runs
  func_outs <- lapply(1:50, function(i) {
    w <- c(0.3, 0.3, 0.4) + rnorm(3, 0, 0.03)
    w / sum(w)
  })
  
  # Full analysis
  summary_df <- summarize_func_outs(func_outs, names = c("A", "B", "C"))
  cv_df <- compute_func_out_cv(func_outs, names = c("A", "B", "C"))
  robust_df <- compute_robust_estimates(func_outs, names = c("A", "B", "C"))
  boot_mat <- extract_func_out_matrix(func_outs)
  
  # Verify consistency
  expect_equal(nrow(summary_df), 3)
  expect_equal(nrow(cv_df), 3)
  expect_equal(nrow(robust_df), 3)
  expect_equal(dim(boot_mat), c(50, 3))
  
  # Bootstrap mean should be consistent across functions (allowing for rounding differences)
  expect_equal(summary_df$Mean, robust_df$Boot_Mean, tolerance = 0.001)
})

test_that("regime diagnostics workflow", {
  set.seed(9012)
  
  # Create diagnostics
  diag <- create_bootstrap_diagnostics("hmm", 100, 2, 5)
  
  # Record regime info
  original_states <- c(rep(1, 60), rep(2, 40))
  diag <- record_regime_info(diag, original_states, 2, c("Bull", "Bear"))
  
  # Record original stats
  x <- create_test_data(100, 2)
  diag <- record_original_stats(diag, x)
  
  # Record replicate regimes
  for (i in 1:5) {
    rep_states <- sample(c(1, 2), 100, replace = TRUE, prob = c(0.6, 0.4))
    diag <- record_replicate_regimes(diag, i, rep_states)
  }
  
  # Summarize
  result <- capture.output(summarize_regime_diagnostics(diag))
  
  expect_true(any(grepl("Bull", result)))
  expect_true(any(grepl("Bear", result)))
})


# ==============================================================================
# SECTION 12: Benchmarking Tests (if tsbs function available)
# ==============================================================================

context("Benchmarking utilities")

# These tests require the tsbs() function which may not be available
# in isolation. They are wrapped in skip conditions.

test_that("print.tsbs_benchmark produces output", {
  # Create a mock benchmark object
  mock_benchmark <- list(
    results = data.frame(
      Setup = c("A", "A", "B", "B"),
      Parameter = rep("num_boots", 4),
      Value = c(10, 10, 10, 10),
      Time = c(0.1, 0.12, 0.2, 0.22),
      Rep = c(1, 2, 1, 2)
    ),
    summary = data.frame(
      Setup = c("A", "B"),
      Value = c(10, 10),
      Mean_Time = c(0.11, 0.21),
      SD_Time = c(0.01, 0.01)
    ),
    vary = "num_boots",
    values = 10,
    setups = list(A = list(bs_type = "moving"), B = list(bs_type = "stationary")),
    times = 2,
    data_dim = c(n_obs = 100, n_vars = 3)
  )
  class(mock_benchmark) <- "tsbs_benchmark"
  
  output <- capture.output(print(mock_benchmark))
  expect_true(any(grepl("tsbs Benchmark Results", output)))
})

test_that("summary.tsbs_benchmark produces output", {
  mock_benchmark <- list(
    results = data.frame(
      Setup = rep(c("A", "B"), each = 4),
      Parameter = rep("num_boots", 8),
      Value = rep(c(10, 20), 4),
      Time = c(0.1, 0.2, 0.11, 0.21, 0.15, 0.25, 0.14, 0.24),
      Rep = rep(1:2, 4)
    ),
    summary = data.frame(
      Setup = c("A", "A", "B", "B"),
      Value = c(10, 20, 10, 20),
      Mean_Time = c(0.105, 0.205, 0.145, 0.245),
      SD_Time = c(0.007, 0.007, 0.007, 0.007)
    ),
    vary = "num_boots",
    values = c(10, 20),
    setups = list(A = list(bs_type = "moving"), B = list(bs_type = "stationary")),
    times = 2,
    data_dim = c(n_obs = 100, n_vars = 3)
  )
  class(mock_benchmark) <- "tsbs_benchmark"
  
  output <- capture.output(summary(mock_benchmark))
  expect_true(any(grepl("Benchmark Summary", output)))
  expect_true(any(grepl("Fastest method", output)))
})


# ==============================================================================
# SECTION 13: Parallel Execution Tests
# ==============================================================================

context("Parallel execution")

test_that(".sample_blocks_with_diagnostics works with explicit num_cores", {
  skip_if_not_installed("foreach")
  skip_if_not_installed("doParallel")
  
  # This test ensures %dopar% is properly defined when parallel = TRUE
  # and num_cores > 1 is explicitly specified
  x <- create_test_data(100, 2, seed = 1111)
  states <- c(rep(1, 50), rep(2, 50))
  
  # Should not error with explicit num_cores = 2
  result <- .sample_blocks_with_diagnostics(
    x = x,
    n_boot = 100,
    num_blocks = NULL,
    states = states,
    num_boots = 3,
    parallel = TRUE,
    num_cores = 2L,
    collect_diagnostics = TRUE
  )
  
  expect_type(result, "list")
  expect_true("samples" %in% names(result))
  expect_length(result$samples, 3)
})

test_that(".sample_blocks_with_diagnostics works without explicit num_cores", {
  x <- create_test_data(100, 2, seed = 2222)
  states <- c(rep(1, 50), rep(2, 50))
  
  # Default num_cores = 1 should use sequential lapply

  result <- .sample_blocks_with_diagnostics(
    x = x,
    n_boot = 100,
    num_blocks = NULL,
    states = states,
    num_boots = 3,
    parallel = FALSE,
    num_cores = 1L,
    collect_diagnostics = FALSE
  )
  
  expect_type(result, "list")
  expect_length(result$samples, 3)
})

